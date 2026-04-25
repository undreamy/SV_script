import pandas as pd
import numpy as np
import re
from scipy import stats
from statsmodels.stats.multitest import multipletests
import seaborn as sns
import matplotlib.pyplot as plt
import os
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import psutil
import sys

warnings.filterwarnings('ignore')


# =========================
# 1. 优化数据读取
# =========================

def load_data():
    sv_file = "target_SV_with_sample.txt"
    expr_file = "tpmss.tsv"

    print("加载SV数据...")
    sv_df = pd.read_csv(sv_file, sep="\t")

    print("加载表达量数据...")
    # 先读取前100行获取列信息
    with open(expr_file, 'r') as f:
        lines = []
        for i, line in enumerate(f):
            if line.startswith('#'):
                continue
            lines.append(line)
            if i >= 100:  # 读取前100行非注释行
                break

    # 从样本行推断数据类型
    sample_df = pd.read_csv(pd.io.common.StringIO(''.join(lines)), sep='\t')

    # 构建dtype字典：Gene ID和Gene Name为object，其他为float32
    dtype_dict = {}
    for col in sample_df.columns:
        if col in ["Gene ID", "Gene Name"]:
            dtype_dict[col] = 'object'
        else:
            dtype_dict[col] = 'float32'

    # 分块读取大文件
    chunks = []
    chunk_size = 10000

    for chunk in pd.read_csv(expr_file, sep='\t', comment='#', dtype=dtype_dict,
                             chunksize=chunk_size, low_memory=False):
        chunks.append(chunk)

    expr_df = pd.concat(chunks, ignore_index=True)

    print(f"数据加载完成: SV形状={sv_df.shape}, 表达量形状={expr_df.shape}")
    return sv_df, expr_df


# =========================
# 2. SV数据处理优化
# =========================

def build_gene_sv_mapping(sv_df):
    sv_df["Samples"] = sv_df["Samples"].fillna("")

    gene_sv_dict = {}
    for _, row in sv_df.iterrows():
        gene = row["GeneID"]
        samples = row["Samples"]

        if gene not in gene_sv_dict:
            gene_sv_dict[gene] = set()

        # 一次性处理所有样本
        lines = {str(s).upper().strip() for s in samples.split(",") if str(s).strip()}
        gene_sv_dict[gene].update(lines)

    return gene_sv_dict


# =========================
# 3. 表达矩阵转长表优化
# =========================

def melt_expression_matrix(expr_df, gene_sv_dict):
    id_cols = ["Gene ID", "Gene Name"]
    value_cols = [c for c in expr_df.columns if c not in id_cols]

    print("转换表达矩阵为长格式...")

    # 分批处理melt，避免内存溢出
    chunk_size = 500  # 每次处理500个样本列
    long_chunks = []

    for i in range(0, len(value_cols), chunk_size):
        chunk_cols = value_cols[i:i + chunk_size]
        chunk_df = expr_df[id_cols + chunk_cols]

        chunk_long = chunk_df.melt(
            id_vars=id_cols,
            value_vars=chunk_cols,
            var_name="Sample",
            value_name="TPM"
        )

        long_chunks.append(chunk_long)

    # 合并所有chunk
    long_df = pd.concat(long_chunks, ignore_index=True)

    # 拆分Sample列
    print("处理样本信息...")
    split_data = long_df["Sample"].str.split(",", n=1, expand=True)
    long_df["Line"] = split_data[0].str.upper().str.strip()
    long_df["Tissue"] = split_data[1].str.strip().fillna("unknown")

    long_df.rename(columns={"Gene ID": "GeneID"}, inplace=True)

    # 添加Group列
    print("添加SV/Ref分组...")

    def assign_group(row):
        gene = row["GeneID"]
        line = row["Line"]

        if gene in gene_sv_dict and line in gene_sv_dict[gene]:
            return "SV"
        return "Ref"

    long_df["Group"] = long_df.apply(assign_group, axis=1)

    # 删除不必要的列节省内存
    if "Gene Name" in long_df.columns:
        long_df.drop(columns=["Gene Name", "Sample"], inplace=True)

    return long_df


# =========================
# 4. 单组织分析
# =========================

def analyze_single_tissue_optimized(df, tissue="leaf", outdir="leaf_plots",
                                    max_plots=50, plot_top_only=True):
    os.makedirs(outdir, exist_ok=True)

    # 筛选组织
    sub = df[df["Tissue"] == tissue].copy()

    if sub.empty:
        print(f"警告: 组织 '{tissue}' 无数据")
        return pd.DataFrame()

    results = []
    genes_to_plot = []

    print(f"分析{tissue}组织，共{len(sub['GeneID'].unique())}个基因...")

    # 使用groupby一次计算
    for gene, gdf in tqdm(sub.groupby("GeneID"), desc=f"分析{tissue}"):
        gdf = gdf.dropna(subset=["TPM"])

        sv_vals = gdf[gdf["Group"] == "SV"]["TPM"].values
        ref_vals = gdf[gdf["Group"] == "Ref"]["TPM"].values

        if len(sv_vals) < 2 or len(ref_vals) < 2:
            continue

        # 使用Mann-Whitney U检验（等价于ranksums但更快）
        try:
            stat, p = stats.mannwhitneyu(sv_vals, ref_vals, alternative='two-sided')
        except:
            continue

        result = {
            "GeneID": gene,
            "mean_SV": np.mean(sv_vals),
            "mean_Ref": np.mean(ref_vals),
            "P_value": p
        }
        results.append(result)

        # 只记录需要绘图的基因
        if plot_top_only and len(genes_to_plot) < max_plots:
            genes_to_plot.append((gene, p, gdf))

    # 批量绘图
    if genes_to_plot:
        print(f"绘制前{len(genes_to_plot)}个基因的箱线图...")

        # 按P值排序
        genes_to_plot.sort(key=lambda x: x[1])

        for gene, p, gdf in tqdm(genes_to_plot, desc="绘图"):
            plt.figure(figsize=(4, 5))

            # 使用兼容性更好的绘图方式
            try:
                # 尝试多种绘图方法
                sns.boxplot(data=gdf, x="Group", y="TPM")
                # 添加散点
                for i, group in enumerate(gdf["Group"].unique()):
                    group_data = gdf[gdf["Group"] == group]["TPM"]
                    x_pos = i + np.random.normal(0, 0.05, len(group_data))  # 轻微抖动
                    plt.scatter(x_pos, group_data, color='black', alpha=0.5, s=20)
            except Exception as e:
                # 如果seaborn失败，使用matplotlib直接绘制
                groups = gdf["Group"].unique()
                data = [gdf[gdf["Group"] == group]["TPM"] for group in groups]

                plt.boxplot(data, labels=groups)
                # 添加散点
                for i, group_data in enumerate(data):
                    x_pos = i + 1 + np.random.normal(0, 0.05, len(group_data))
                    plt.scatter(x_pos, group_data, color='black', alpha=0.5, s=20)

            plt.title(f"{gene} ({tissue})\nP={p:.3e}")
            plt.tight_layout()
            plt.savefig(f"{outdir}/{gene}.png", dpi=150, bbox_inches='tight')
            plt.close()

    res_df = pd.DataFrame(results)
    if not res_df.empty:
        output_file = f"{tissue}_results.csv"
        res_df.to_csv(output_file, index=False)
        print(f"结果已保存到: {output_file}")

    return res_df


# =========================
# 5. Fisher合并优化
# =========================

def fisher_method(pvals):
    """向量化Fisher方法"""
    pvals = np.array(pvals)
    pvals = pvals[~np.isnan(pvals) & (pvals > 0) & (pvals <= 1)]

    if len(pvals) == 0:
        return np.nan

    chi2 = -2 * np.sum(np.log(pvals))
    df = 2 * len(pvals)

    return stats.chi2.sf(chi2, df)


# =========================
# 6. 并行化全组织分析
# =========================

def analyze_gene_parallel(args):
    gene, tissue_groups, group_col = args

    pvals = []
    tissue_means = {}

    for tissue, tdf in tissue_groups:
        tdf = tdf.dropna(subset=["TPM"])

        sv_vals = tdf[tdf[group_col] == "SV"]["TPM"].values
        ref_vals = tdf[tdf[group_col] == "Ref"]["TPM"].values

        if len(sv_vals) < 2 or len(ref_vals) < 2:
            continue

        try:
            stat, p = stats.mannwhitneyu(sv_vals, ref_vals, alternative='two-sided')
            pvals.append(p)

            tissue_means[f"{tissue}_SV"] = np.mean(sv_vals)
            tissue_means[f"{tissue}_Ref"] = np.mean(ref_vals)
        except:
            continue

    combined_p = fisher_method(pvals)

    if np.isnan(combined_p):
        return None

    row = {
        "GeneID": gene,
        "Combined_P": combined_p,
        "Num_Tissues": len(pvals)
    }
    row.update(tissue_means)

    return row


def analyze_all_tissues_parallel(df, max_workers=None):
    """并行化全组织分析"""
    if max_workers is None:
        max_workers = max(1, psutil.cpu_count(logical=False) - 1)  # 物理核心数-1

    print(f"使用{max_workers}个进程进行并行分析...")

    # 预处理：按基因和组织分组
    print("数据预处理...")
    gene_groups = []
    for gene, gdf in df.groupby("GeneID"):
        tissue_groups = list(gdf.groupby("Tissue"))
        if len(tissue_groups) >= 2:  # 至少两个组织才分析
            gene_groups.append((gene, tissue_groups, "Group"))

    print(f"共{len(gene_groups)}个基因需要分析")

    results = []

    # 使用进程池并行处理
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(analyze_gene_parallel, args): args[0]
                   for args in gene_groups}

        for future in tqdm(as_completed(futures), total=len(futures), desc="并行分析"):
            try:
                result = future.result(timeout=30)  # 30秒超时
                if result is not None:
                    results.append(result)
            except Exception as e:
                gene_id = futures[future]
                print(f"基因{gene_id}分析失败: {e}")
                continue

    if not results:
        print("无有效结果")
        return pd.DataFrame()

    res_df = pd.DataFrame(results)

    # 多重检验校正
    print("进行多重检验校正...")
    res_df = res_df.dropna(subset=["Combined_P"])

    if not res_df.empty:
        res_df["FDR"] = multipletests(res_df["Combined_P"], method="fdr_bh")[1]

        # 按FDR排序
        res_df = res_df.sort_values("FDR")

        output_file = "gene_level_sv_expression_results.csv"
        res_df.to_csv(output_file, index=False)
        print(f"结果已保存到: {output_file}")
        print(f"显著基因(FDR<0.05): {len(res_df[res_df['FDR'] < 0.05])}个")

    return res_df


# =========================
# 7. 修复多组织可视化（最终版）
# =========================

def plot_multi_tissue_fixed(df, genes, outdir="multi_tissue_plots", max_genes=5):
    os.makedirs(outdir, exist_ok=True)

    if len(genes) > max_genes:
        print(f"只绘制前{max_genes}个基因")
        genes = genes[:max_genes]

    for idx, gene in enumerate(tqdm(genes, desc="绘制多组织图")):
        gdf = df[df["GeneID"] == gene].copy()

        if gdf.empty:
            print(f"警告: 基因 {gene} 无数据")
            continue

        # 深度清洗数据
        gdf = gdf.dropna(subset=["TPM", "Tissue", "Group"])
        gdf = gdf[gdf["TPM"].notna() & (gdf["TPM"] >= 0)]

        if gdf.empty:
            print(f"警告: 基因 {gene} 清洗后无有效数据")
            continue

        # 检查Group列是否有有效值
        valid_groups = gdf["Group"].unique()
        if len(valid_groups) < 2:
            print(f"警告: 基因 {gene} 只有{len(valid_groups)}个有效组别")
            continue

        # 获取组织和组别
        tissues = sorted(gdf["Tissue"].unique())
        n_tissues = len(tissues)

        if n_tissues == 0:
            continue

        # 创建图形
        fig, ax = plt.subplots(figsize=(max(8, n_tissues * 1.5), 6))

        # 准备数据
        sv_data_by_tissue = []
        ref_data_by_tissue = []

        for tissue in tissues:
            tissue_data = gdf[gdf["Tissue"] == tissue]
            sv_data = tissue_data[tissue_data["Group"] == "SV"]["TPM"].values
            ref_data = tissue_data[tissue_data["Group"] == "Ref"]["TPM"].values

            # 清洗NaN和无穷大值
            sv_data = sv_data[~np.isnan(sv_data)]
            sv_data = sv_data[~np.isinf(sv_data)]
            ref_data = ref_data[~np.isnan(ref_data)]
            ref_data = ref_data[~np.isinf(ref_data)]

            sv_data_by_tissue.append(sv_data)
            ref_data_by_tissue.append(ref_data)

        # 设置位置
        positions_sv = np.arange(n_tissues) - 0.2
        positions_ref = np.arange(n_tissues) + 0.2

        # 绘制SV组
        for i, sv_data in enumerate(sv_data_by_tissue):
            if len(sv_data) > 0:
                # 绘制箱线图
                bp_sv = ax.boxplot([sv_data], positions=[positions_sv[i]], widths=0.3,
                                   patch_artist=True, boxprops=dict(facecolor='lightcoral'))
                # 添加散点
                sv_x = np.full(len(sv_data), positions_sv[i]) + np.random.normal(0, 0.03, len(sv_data))
                ax.scatter(sv_x, sv_data, color='red', alpha=0.6, s=30, edgecolor='black', linewidth=0.5, zorder=5)

        # 绘制Ref组
        for i, ref_data in enumerate(ref_data_by_tissue):
            if len(ref_data) > 0:
                # 绘制箱线图
                bp_ref = ax.boxplot([ref_data], positions=[positions_ref[i]], widths=0.3,
                                    patch_artist=True, boxprops=dict(facecolor='lightblue'))
                # 添加散点
                ref_x = np.full(len(ref_data), positions_ref[i]) + np.random.normal(0, 0.03, len(ref_data))
                ax.scatter(ref_x, ref_data, color='blue', alpha=0.6, s=30, edgecolor='black', linewidth=0.5, zorder=5)

        # 设置x轴标签
        ax.set_xticks(np.arange(n_tissues))
        ax.set_xticklabels(tissues, rotation=45, ha='right')

        # 手动创建图例
        from matplotlib.patches import Patch

        legend_elements = [
            Patch(facecolor='lightcoral', edgecolor='red', label='SV'),
            Patch(facecolor='lightblue', edgecolor='blue', label='Ref')
        ]
        ax.legend(handles=legend_elements, loc='best')

        ax.set_title(f"Gene: {gene}")
        ax.set_ylabel("TPM")
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        output_path = f"{outdir}/{gene}.png"
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()

        print(f"已保存: {output_path}")


# =========================
# 8. 极简多组织可视化
# =========================

def plot_simple_multi_tissue_fixed(df, genes, outdir="simple_multi_tissue_plots_fixed", max_genes=5):
    """极简但最稳定的多组织可视化"""
    os.makedirs(outdir, exist_ok=True)

    if len(genes) > max_genes:
        print(f"只绘制前{max_genes}个基因")
        genes = genes[:max_genes]

    for gene in tqdm(genes, desc="绘制极简多组织图"):
        gdf = df[df["GeneID"] == gene].copy()

        if gdf.empty:
            continue

        # 彻底清洗数据
        gdf = gdf.dropna(subset=["TPM", "Tissue", "Group"])
        gdf = gdf[gdf["TPM"].notna() & (gdf["TPM"] >= 0)]

        if gdf.empty:
            continue

        # 获取组织和组别
        tissues = sorted(gdf["Tissue"].unique())
        n_tissues = len(tissues)

        fig, ax = plt.subplots(figsize=(max(6, n_tissues * 1.2), 5))

        for i, tissue in enumerate(tissues):
            tissue_data = gdf[gdf["Tissue"] == tissue]

            # SV组
            sv_data = tissue_data[tissue_data["Group"] == "SV"]["TPM"].values
            sv_data = sv_data[~np.isnan(sv_data)]
            sv_data = sv_data[~np.isinf(sv_data)]

            if len(sv_data) > 0:
                x_pos = i - 0.15
                # 绘制均值
                sv_mean = np.mean(sv_data)
                sv_std = np.std(sv_data)

                plt.plot([x_pos - 0.1, x_pos + 0.1], [sv_mean, sv_mean], 'r-', linewidth=2,
                         label='SV Mean' if i == 0 else "")
                # 绘制标准差
                plt.plot([x_pos, x_pos], [sv_mean - sv_std, sv_mean + sv_std], 'r-', linewidth=1, alpha=0.7)
                plt.plot([x_pos - 0.05, x_pos + 0.05], [sv_mean - sv_std, sv_mean - sv_std], 'r-', linewidth=1,
                         alpha=0.7)
                plt.plot([x_pos - 0.05, x_pos + 0.05], [sv_mean + sv_std, sv_mean + sv_std], 'r-', linewidth=1,
                         alpha=0.7)

                # 绘制散点
                sv_x = np.full(len(sv_data), x_pos) + np.random.normal(0, 0.02, len(sv_data))
                plt.scatter(sv_x, sv_data, color='red', alpha=0.5, s=20, label='SV Data' if i == 0 else "")

            # Ref组
            ref_data = tissue_data[tissue_data["Group"] == "Ref"]["TPM"].values
            ref_data = ref_data[~np.isnan(ref_data)]
            ref_data = ref_data[~np.isinf(ref_data)]

            if len(ref_data) > 0:
                x_pos = i + 0.15
                # 绘制均值
                ref_mean = np.mean(ref_data)
                ref_std = np.std(ref_data)

                plt.plot([x_pos - 0.1, x_pos + 0.1], [ref_mean, ref_mean], 'b-', linewidth=2,
                         label='Ref Mean' if i == 0 else "")
                # 绘制标准差
                plt.plot([x_pos, x_pos], [ref_mean - ref_std, ref_mean + ref_std], 'b-', linewidth=1, alpha=0.7)
                plt.plot([x_pos - 0.05, x_pos + 0.05], [ref_mean - ref_std, ref_mean - ref_std], 'b-', linewidth=1,
                         alpha=0.7)
                plt.plot([x_pos - 0.05, x_pos + 0.05], [ref_mean + ref_std, ref_mean + ref_std], 'b-', linewidth=1,
                         alpha=0.7)

                # 绘制散点
                ref_x = np.full(len(ref_data), x_pos) + np.random.normal(0, 0.02, len(ref_data))
                plt.scatter(ref_x, ref_data, color='blue', alpha=0.5, s=20, label='Ref Data' if i == 0 else "")

        plt.xticks(range(n_tissues), tissues, rotation=45, ha='right')
        plt.title(f"Gene: {gene}")
        plt.ylabel("TPM")
        plt.grid(True, alpha=0.3)

        # 简化图例
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        if by_label:
            plt.legend(by_label.values(), by_label.keys())

        plt.tight_layout()

        output_path = f"{outdir}/{gene}_simple.png"
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()

        print(f"已保存: {output_path}")


# =========================
# 9. 主流程
# =========================

def main():
    """主函数"""
    print("=" * 50)
    print("SV表达差异分析")
    print("=" * 50)

    # 1. 加载数据
    sv_df, expr_df = load_data()

    # 2. 构建SV映射
    print("\n构建基因-SV映射...")
    gene_sv_dict = build_gene_sv_mapping(sv_df)

    # 3. 转换表达矩阵
    print("\n处理表达矩阵...")
    long_df = melt_expression_matrix(expr_df, gene_sv_dict)

    # 释放内存
    del expr_df
    del sv_df

    print(f"长格式数据形状: {long_df.shape}")
    print(f"唯一基因数: {long_df['GeneID'].nunique()}")
    print(f"组织类型: {long_df['Tissue'].unique()}")

    # 检查数据质量
    print(f"\n数据质量检查:")
    print(f"TPM列NaN值数量: {long_df['TPM'].isna().sum()}")
    print(f"TPM列无穷大值数量: {np.isinf(long_df['TPM']).sum()}")
    print(f"TPM列负值数量: {(long_df['TPM'] < 0).sum()}")

    # 清理数据
    long_df_clean = long_df.dropna(subset=["TPM", "Tissue", "Group"]).copy()
    long_df_clean = long_df_clean[~np.isinf(long_df_clean["TPM"])]
    long_df_clean = long_df_clean[long_df_clean["TPM"] >= 0]

    print(f"清洗后数据形状: {long_df_clean.shape}")
    print(f"清洗后有效组织数: {long_df_clean['Tissue'].nunique()}")

    # 4. 单组织分析（可选）
    print("\n" + "=" * 30)
    print("单组织分析 (leaf)")
    print("=" * 30)

    # 检查是否有leaf组织
    if "leaf" in long_df_clean["Tissue"].unique():
        leaf_res = analyze_single_tissue_optimized(
            long_df_clean,
            tissue="leaf",
            outdir="leaf_plots",
            max_plots=20,
            plot_top_only=True
        )

        if not leaf_res.empty:
            print(f"leaf组织分析完成，共{len(leaf_res)}个基因")
            print(f"最显著基因: {leaf_res.iloc[0]['GeneID']} (P={leaf_res.iloc[0]['P_value']:.2e})")
    else:
        print("无leaf组织数据，跳过单组织分析")

    # 5. 全组织分析（并行）
    print("\n" + "=" * 30)
    print("全组织分析")
    print("=" * 30)

    gene_res = analyze_all_tissues_parallel(long_df_clean)

    if not gene_res.empty:
        # 6. 可视化显著基因
        sig_genes = gene_res[gene_res["FDR"] < 0.05]["GeneID"].tolist()

        print(f"\n发现{len(sig_genes)}个显著差异基因(FDR<0.05)")

        if sig_genes:
            # 取FDR最小的前几个
            top_genes = gene_res.head(min(5, len(gene_res)))["GeneID"].tolist()

            print(f"Top {len(top_genes)}显著基因: {top_genes}")

            print("\n选项1: 绘制多组织图...")
            plot_multi_tissue_fixed(long_df_clean, top_genes, max_genes=5)

            print("\n选项2: 绘制极简版多组织图...")
            plot_simple_multi_tissue_fixed(long_df_clean, top_genes[:3], max_genes=3)

        # 7. 保存结果摘要
        summary = {
            "总分析基因数": len(gene_res),
            "显著基因数(FDR<0.05)": len(sig_genes),
            "最小FDR": gene_res["FDR"].min() if not gene_res.empty else np.nan,
            "涉及组织数": long_df_clean["Tissue"].nunique(),
            "SV样本数": len(long_df_clean[long_df_clean["Group"] == "SV"]["Line"].unique()),
            "Ref样本数": len(long_df_clean[long_df_clean["Group"] == "Ref"]["Line"].unique()),
            "清洗后总样本数": len(long_df_clean),
            "原始TPM NaN数": long_df["TPM"].isna().sum(),
            "原始TPM 负值数": (long_df["TPM"] < 0).sum()
        }

        summary_df = pd.DataFrame([summary])
        summary_df.to_csv("analysis_summary.csv", index=False)
        print(f"\n分析摘要已保存到: analysis_summary.csv")

    print("\n" + "=" * 30)
    print("分析完成!")
    print("=" * 30)

    # 返回重要结果路径
    result_files = []
    if os.path.exists("gene_level_sv_expression_results.csv"):
        result_files.append("gene_level_sv_expression_results.csv")
    if os.path.exists("analysis_summary.csv"):
        result_files.append("analysis_summary.csv")

    if result_files:
        print(f"\n生成的结果文件:")
        for f in result_files:
            print(f"  - {f}")


if __name__ == "__main__":
    main()
