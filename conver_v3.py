import sys

# -----------------------------
# 1 输入输出文件
# -----------------------------

input_file = "annotation_table.txt"
output_file = "maize_long_table.txt"

print("Step 1: 读取并解析数据...")

all_rows = []

with open(input_file, 'r', encoding='utf-8') as fin:

    header = next(fin)

    for line_num, line in enumerate(fin, 2):

        parts = line.rstrip('\n').split('\t', 8)

        if len(parts) < 9:
            continue

        chrom = parts[0]
        pos = parts[1]
        end = parts[2]
        sv_id = parts[3]
        ref = parts[4]
        alt = parts[5]

        gene_str = parts[6]
        effect_str = parts[7]
        impact_str = parts[8]

        # -----------------------------
        # 数值转换保护
        # -----------------------------

        try:
            pos = int(pos)
        except:
            pos = 0

        try:
            end = int(end)
        except:
            end = pos

        # -----------------------------
        # 拆分 annotation
        # -----------------------------

        genes = [] if gene_str == '.' else gene_str.split(',')
        effects = [] if effect_str == '.' else effect_str.split(',')
        impacts = [] if impact_str == '.' else impact_str.split(',')

        max_len = max(len(genes), len(effects), len(impacts))

        for i in range(max_len):

            gene = genes[i] if i < len(genes) else '.'
            effect = effects[i] if i < len(effects) else '.'
            impact = impacts[i] if i < len(impacts) else '.'

            if gene == '.' and effect == '.' and impact == '.':
                continue

            all_rows.append({
                'chrom': chrom,
                'start': pos,
                'end': end,
                'id': sv_id,
                'gene': gene,
                'effect': effect,
                'impact': impact
            })

print(f"解析完成，共 {len(all_rows)} 行")

# -----------------------------
# 2 染色体自然排序
# -----------------------------

def get_sort_key(row):

    c = row['chrom']
    p = row['start']

    clean = c.lower().replace('chr', '')

    try:
        chrom_num = int(clean)
        is_number = 0
    except ValueError:
        chrom_num = clean
        is_number = 1

    return (is_number, chrom_num, p)

print("Step 2: 排序中...")
all_rows.sort(key=get_sort_key)

# -----------------------------
# 3 写出结果
# -----------------------------

print("Step 3: 写入文件...")

with open(output_file, 'w', encoding='utf-8') as fout:

    fout.write("CHROM\tSTART\tEND\tID\tGeneID\tEffect\tImpact\n")

    for r in all_rows:

        fout.write(
            f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['id']}\t{r['gene']}\t{r['effect']}\t{r['impact']}\n"
        )

print("完成！")