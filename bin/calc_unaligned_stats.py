#!/usr/bin/env python3
import pandas as pd
import re
import sys

if len(sys.argv) != 4:
    print("Usage: calc_unaligned_stats.py <sample_id> <unaligned_file> <alignment_tsv>", file=sys.stderr)
    sys.exit(1)

sample_id, unaligned_file, alignment_tsv = sys.argv[1], sys.argv[2], sys.argv[3]

# --- Read unaligned contigs ---
with open(unaligned_file) as f:
    unaligned_contigs = {line.strip() for line in f if line.strip()}

# --- Parse all_alignments TSV ---
data = []
with open(alignment_tsv) as f:
    for line in f:
        m = re.search(r'(NODE_\d+_length_\d+_cov_[\d.]+)', line)
        if m:
            contig = m.group(1)
            mm = re.search(r'length_(\d+)_cov_([\d.]+)', contig)
            if mm:
                length = float(mm.group(1))
                cov = float(mm.group(2))
                total = length * cov
                data.append((contig, length, cov, total))

if not data:
    print(f"No contigs parsed for {sample_id}")
    sys.exit(0)

df = pd.DataFrame(data, columns=['contig', 'length', 'cov', 'bases'])

total_bases = df['bases'].sum()
unaligned_bases = df[df['contig'].isin(unaligned_contigs)]['bases'].sum()
pct = unaligned_bases / total_bases * 100 if total_bases > 0 else 0

# --- Write CSV ---
df_out = pd.DataFrame([{
    'sample_id': sample_id,
    'total_bases': f"{int(round(total_bases))}",        
    'unaligned_bases': f"{int(round(unaligned_bases))}", 
    'unaligned_pct': f"{pct:.8f}"                        
}])

df_out.to_csv(f"{sample_id}_unaligned_summary.csv", index=False)
