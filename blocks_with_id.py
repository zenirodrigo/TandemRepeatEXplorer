#!/usr/bin/env python3
import sys

bed_file=sys.argv[1]
anchors_file=sys.argv[2]
output_file=sys.argv[3]

gene_pos={}
with open(bed_file) as f:
    for line in f:
        if not line.strip():
            continue
        cols=line.strip().split()
        if len(cols)<4:
            continue
        gene_pos[cols[3]]=(cols[0],int(cols[1]),int(cols[2]))

print(f"[INFO] Loaded {len(gene_pos)} genes")

total=matched=written=0
with open(anchors_file) as f, open(output_file,"w") as out:
    for i,line in enumerate(f):
        if not line.strip():
            continue
        cols=line.strip().split()
        if len(cols)<4:
            continue
        total+=1
        g_start,g_end=cols[0],cols[1]
        if g_start not in gene_pos or g_end not in gene_pos:
            continue
        matched+=1
        chrom1,s1,e1=gene_pos[g_start]
        chrom2,s2,e2=gene_pos[g_end]
        out.write(f"{chrom1}\t{min(s1,s2)}\t{max(e1,e2)}\tblock_{i}\n")
        written+=1

print(f"[INFO] Total lines: {total}")
print(f"[INFO] Matched: {matched}")
print(f"[INFO] Written: {written}")
