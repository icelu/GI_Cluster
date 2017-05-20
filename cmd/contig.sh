# find the distribution of contig size

faCount Roseobacter_species_M05.fasta
# https://github.com/hyattpd/Prodigal/issues/11
/home/b/bingxin/GI-Cluster/bin/prodigal -p meta -m -i Roseobacter_species_M05.fasta -o Roseobacter_species_M05.gbk -a Roseobacter_species_M05.faa -d Roseobacter_species_M05.fna



source ~/.bashrc

gnome=Roseobacter
organism=Roseobacter_species_M05
pred_prog=prodigal
seg_prog=window
prog_dir=/home/b/bingxin/GI-Cluster
output_dir=/home/b/bingxin/genome/draft/$gnome
mode=2
sh $prog_dir/GI_Segmentation.sh -s $prog_dir -o $output_dir -n $organism -m $seg_prog -p $pred_prog -b $mode

sh $prog_dir/GI_Feature.sh -s $prog_dir -o $output_dir -n $organism -m $seg_prog -p $pred_prog

nohup /usr/bin/time sh /home/b/bingxin/GI-Cluster/GI-Cluster.sh -s /home/b/bingxin/GI-Cluster -o /home/b/bingxin/genome/draft/$gnome -n $organism -m $seg_prog -p $pred_prog -d 16 > std_"$seg_prog"_"$pred_prog" 2>&1 &


repseek -l 15 -O 0 -r Roseobacter_species_M05.repseek Roseobacter_species_M05.fna


organism=Roseobacter_species_M05
tRNAscan-SE -B --frag $organism.trna_frag -o $organism.trna_pred -m $organism.trna_stat --brief $organism.fna



prog_dir=/home/b/bingxin/GI-Cluster
output_dir=/home/b/bingxin/genome/draft/$gnome
$prog_dir/GI_Feature.sh -s $prog_dir -o $output_dir -n $organism -m $seg_prog -p $pred_prog



python $prog_dir/evaluation/eval4contigs.py -q gicluster.gilist -r NZ_ACHX00000000.refgi


python $prog_dir/evaluation/eval4contigs.py -q gicluster_unannotated.gilist -r NZ_ACHX00000000.refgi


python $prog_dir/evaluation/eval4contigs.py -q islandviewer4.gilist -r NZ_ACHX00000000.refgi



# Find genes on each contig
for i in {1..9..1}
do
  less NZ_ACHX0100000"$i".gbk | grep locus_tag | sort | uniq > gene_c"$i"
done

for i in {10..10..1}
do
  less NZ_ACHX010000"$i".gbk | grep locus_tag | sort | uniq > gene_c"$i"
done
