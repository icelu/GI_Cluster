Input files for WH8102 (downloaded from NCBI):
GCF_000195975.1_ASM19597v1.fna  -- whole genome sequence
GCF_000195975.1_ASM19597v1_protein.faa  -- protein sequence
GCF_000195975.1_ASM19597v1_cds_from_genomic.fna -- DNA sequence of CDS regions

Input files for VcRC9 (downloaded from NCBI):
NZ_ACHX00000000.fna  -- DMA sequences of genomic contigs


Sample commands for running GI-Cluster on genome WH8102:
```
gnome=WH8102
organism=GCF_000195975.1_ASM19597v1
pred_prog=ncbi
seg_prog=window
output_dir=/home/b/bingxin/cyanobacteria/genome/$gnome
prog_dir=/home/b/bingxin/GI-Cluster
nohup /usr/bin/time sh $prog_dir/GI-Cluster.sh -s $prog_dir -o $output_dir -n $organism -m $seg_prog -p $pred_prog -d 16 > std_"$seg_prog"_"$pred_prog" 2>&1 &
```

Sample commands for running GI-Cluster on genome VcRC9:
```
gnome=VcRC9
organism=NZ_ACHX00000000
pred_prog=prodigal
seg_prog=window
prog_dir=/home/b/bingxin/GI-Cluster
nohup /usr/bin/time sh $prog_dir/GI-Cluster.sh -s /home/b/bingxin/GI-Cluster -o $output_dir -n $organism -m $seg_prog -p $pred_prog -d 16 -b 1 > std_"$seg_prog"_"$pred_prog" 2>&1 &
```


To save space, the following output files are removed:
1. The output folders which contatin the output of BLAST for COGs (Clusters of Orthologous Groups) query, including BLASTcogn, BLASTff, BLASTno, and BLASTss;
2. The output folder which contatins the output of seperate clusterings:  $pred_prog/$seg_prog/comp_content/average/1/1/sepcluster;
3. The output folder which contatins files mainly regarding the input genome: $pred_prog/genome;
4. The output folder which contatins files for gene features: $pred_prog/feature.
5. The output folder which contatins files for boundary features: boundary.
6. The output folder which contatins files for window features and visualization: $seg_prog.
