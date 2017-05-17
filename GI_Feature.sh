# Script for extracting features related to genomic islands in the unit of genes
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg


########## usage ##########
# ./GI_Feature.sh -s  $prog_dir -o $output_dir -n $organism -m $seg_prog -b $mode
# e.g. ./GI_Feature.sh -s ./GIFilter -o ./research/data/species/cft73 -n NC_004431  -m ./research/software/HGT/mjsd -b 0


software=$(basename $0)

function usage() {
  echo -e "GI_Feature: extracting features related to genomic islands in a genomic region"
  echo "Version 1.0
Usage: $software [options] -s [the directory containing all the scripts] -o [the output directory]
-n [the name of the organism (e.g. NC_003198)] -m [programs for genome segmation (e.g. mjsd, gcprofile, gisvm, alienhunter)] -p [programs for gene prediction (e.g. prodigal, ncbi)]

OPTIONS	Default	DESCIPTION
-b	0	: mode of running: 0 for complete genome, 1 for contigs without gene predictions, 2 for contigs with gene predictions.

-e	1e-5	: e-value used during identification of phage-related genes, i.e., blastp against PHAST.
-r	1e-5	: e-value used during identification of virulence factors, i.e., blastp against VFDB.
-a	1e-5	: e-value used during identification of antibiotic resistance genes, i.e., blastp against CARD.

-d	4	: number of threads used by blast.
-u	16 : number of CPUs used by cmsearch.

-h 	----	: print this help
  "
  exit -1
}


phage_evalue=1e-5
virdb_evalue=1e-5
arg_evalue=1e-5

num_threads=4
num_cpus=16

mode=0


while getopts "b:s:o:n:m:p:c:e:r:a:d:u:f:h" OPT; do
  case $OPT in
    b) mode=$OPTARG || exit 1;;

    s) prog_dir=$OPTARG || exit 1;;
    o) output_dir=$OPTARG || exit 1;;
    n) organism=$OPTARG || exit 1;;
    m) seg_prog=$OPTARG || exit 1;;
    p) pred_prog=$OPTARG || exit 1;;

    #t) mobgene_evalue=$OPTARG || exit 1;;
    e) phage_evalue=$OPTARG || exit 1;;
    r) virdb_evalue=$OPTARG || exit 1;;
    a) arg_evalue=$OPTARG || exit 1;;

    d) num_threads=$OPTARG || exit 1;;
    u) num_cpus=$OPTARG || exit 1;;

    h) usage && exit;;
  esac
done


if [ $mode != 1 ]
then
############## compute known features related to genomic islands #########################
if [ ! -d $output_dir/$pred_prog/feature/ ]
then
  mkdir $output_dir/$pred_prog/feature/
fi

############## compute compositional bias #########################
# GC bias
echo "##########################################"
echo "Analyzing GC Content for each gene/ORF"
if [ ! -f $output_dir/$pred_prog/feature/$organism.feature.gc ]
then
  python $prog_dir/feature/analyze_GC.py -g $output_dir/$pred_prog/genome/$organism.ffn -o $output_dir/$pred_prog/feature/$organism.feature.gc
fi

# codon usage
echo "##########################################"
echo "Analyzing codon usage for each gene/ORF"
if [ ! -f $output_dir/$pred_prog/feature/$organism.feature.codon ]
then
if [ ! -f $output_dir/$pred_prog/feature/codonw.out ]
then
  $prog_dir/bin/codonW/codonw -all_indices -nomenu $output_dir/$pred_prog/genome/$organism.ffn $output_dir/$pred_prog/feature/codonw.out $output_dir/$pred_prog/feature/codonw.blk
fi
  less $output_dir/$pred_prog/feature/codonw.out | cut -f6-8 | sed '1d' > $output_dir/$pred_prog/feature/codonw.out.f3
  python $prog_dir/feature/analyze_codon.py -g $output_dir/$pred_prog/genome/$organism.ffn -o $output_dir/$pred_prog/feature/codon.out
  paste $output_dir/$pred_prog/feature/codon.out $output_dir/$pred_prog/feature/codonw.out.f3 > $output_dir/$pred_prog/feature/$organism.feature.codon
fi

# k-mer Frequency
echo "##########################################"
echo "Analyzing k-mer frequency for each gene/ORF"
if [ ! -f $output_dir/$pred_prog/feature/$organism.feature.kmer ]
then
  python $prog_dir/feature/analyze_kmer.py -i $output_dir/$organism.fna -g $output_dir/$pred_prog/genome/$organism.ffn -o $output_dir/$pred_prog/feature/$organism.kmer
  # Add addtional extraction step in case that distance other than covariance is used
  less $output_dir/$pred_prog/feature/$organism.kmer.covariance | cut -f2- > $output_dir/$pred_prog/feature/$organism.feature.kmer
fi

######################## compute gene content ##############################
echo "##########################################"
echo "Predicting ncRNA"
if [ ! -f $output_dir/$pred_prog/feature/$organism.cmscan ]
then
  #cmsearch --cpu $num_cpus --nohmmonly --rfam --cut_ga --tblout $output_dir/$pred_prog/feature/$organism.tbl $prog_dir/db/RNA/CMs/Rfam.cm  $output_dir/$organism.fna > $output_dir/$pred_prog/feature/$organism.cmsearch
  cmscan --cpu $num_cpus --rfam --cut_ga --nohmmonly --tblout $output_dir/$pred_prog/feature/$organism.cmscan.tbl --oskip --fmt 2 $prog_dir/db/RNA/CMs/Rfam.cm  $output_dir/$organism.fna > $output_dir/$pred_prog/feature/$organism.cmscan
fi
#perl $prog_dir/parse_infernal_output.py -i $output_dir/$pred_prog/feature/$organism.cmscan.tbl -o $output_dir/$pred_prog/feature/$organism.rna

# mobility-related gene
echo "##########################################"
echo "Identifying mobility-related genes"
if [ ! -f $prog_dir/db/MOB/Pfam_mobgene.hmm ]
then
  python $prog_dir/feature/extract_mobgene.py -i $prog_dir/db/MOB/Pfam-A.hmm.dat -o $prog_dir/db/MOB/mobgene.list
  hmmfetch -f $prog_dir/db/MOB/Pfam-A.hmm $prog_dir/db/MOB/mobgene.list > $prog_dir/db/MOB/Pfam_mobgene.hmm
fi
if [ ! -f $output_dir/$pred_prog/feature/$organism.feature.mobgene ]
then
  hmmsearch $prog_dir/db/MOB/Pfam_mobgene.hmm $output_dir/$pred_prog/genome/$organism.faa > $output_dir/$pred_prog/feature/$organism.hit.mobgene
  python $prog_dir/feature/parse_hmmer_output.py -d $output_dir/$pred_prog/feature/$organism.hit.mobgene -g $output_dir/$pred_prog/genome/$organism.gene_id -o $output_dir/$pred_prog/feature/$organism.feature.mobgene
fi

# phage-related gene
echo "##########################################"
echo "Identifying phage-related genes"
if [ ! -f $prog_dir/db/PHAST/prophage_virus.pin ]
then
  makeblastdb -in $prog_dir/db/PHAST/prophage_virus.db -title prophage_virus -dbtype prot -out $prog_dir/db/PHAST/prophage_virus
fi
if [ ! -f $output_dir/$pred_prog/feature/$organism.feature.phage ]
then
  blastp -query $output_dir/$pred_prog/genome/$organism.faa -db $prog_dir/db/PHAST/prophage_virus -out $output_dir/$pred_prog/feature/$organism.hit.phage -evalue $phage_evalue -outfmt 6  -num_alignments 1 -num_threads $num_threads
  python $prog_dir/feature/parse_blast_output.py -b $output_dir/$pred_prog/feature/$organism.hit.phage -g $output_dir/$pred_prog/genome/$organism.gene_id -o $output_dir/$pred_prog/feature/$organism.feature.phage
fi

# virulence factor
echo "##########################################"
echo "Identifying virulence factors"
# setB refers to the full set, setA contains only experimentally validated VFs
if [ ! -f $prog_dir/db/VFDB/VFDB_setB_pro.pin ]
then
  makeblastdb -in $prog_dir/db/VFDB/VFDB_setB_pro.fas -title VFDB_setB_pro -dbtype prot -out $prog_dir/db/VFDB/VFDB_setB_pro
fi
if [ ! -f $output_dir/$pred_prog/feature/$organism.feature.vfdb ]
then
  blastp -query $output_dir/$pred_prog/genome/$organism.faa -db $prog_dir/db/VFDB/VFDB_setB_pro -out $output_dir/$pred_prog/feature/$organism.hit.vfdb -evalue $virdb_evalue -outfmt 6  -num_alignments 1 -num_threads $num_threads
  python $prog_dir/feature/parse_blast_output.py -b $output_dir/$pred_prog/feature/$organism.hit.vfdb -g $output_dir/$pred_prog/genome/$organism.gene_id -o $output_dir/$pred_prog/feature/$organism.feature.vfdb
fi


# antibiotic resistance gene
echo "##########################################"
echo "Identifying antibiotic resistance genes"
if [ ! -f $prog_dir/db/CARD/CARD_prot_homolog.pin ]
then
  makeblastdb -in $prog_dir/db/CARD/protein_fasta_protein_homolog_model.fasta  -title CARD_prot_homolog -dbtype prot -out $prog_dir/db/CARD/CARD_prot_homolog
fi
if [ ! -f $output_dir/$pred_prog/feature/$organism.feature.AR ]
then
  blastp -query $output_dir/$pred_prog/genome/$organism.faa -db $prog_dir/db/CARD/CARD_prot_homolog  -out $output_dir/$pred_prog/feature/$organism.hit.AR -evalue $arg_evalue -outfmt 6  -num_alignments 1 -num_threads $num_threads
  python $prog_dir/feature/parse_blast_output.py -b $output_dir/$pred_prog/feature/$organism.hit.AR -g $output_dir/$pred_prog/genome/$organism.gene_id -o $output_dir/$pred_prog/feature/$organism.feature.AR
fi


# Novel gene
echo "##########################################"
echo "Identifying novel genes"
if [ ! -f $prog_dir/db/COG/COGs.pin ]
then
  makeblastdb -in $prog_dir/db/COG/prot2003-2014.fa -dbtype prot -out $prog_dir/db/COG/COGs
fi

if [ ! -f $prog_dir/db/COG/COG.p2o.csv ]
then
  less $prog_dir/db/COG/cog2003-2014.csv  | cut -d',' -f3 > $prog_dir/db/COG/cog_protid
  less $prog_dir/db/COG/cog2003-2014.csv  | cut -d',' -fe > $prog_dir/db/COG/cog_gnomeid
  paste -d',' $prog_dir/db/COG/cog_protid $prog_dir/db/COG/cog_gnomeid > $prog_dir/db/COG/COG.p2o.csv
fi

if [ ! -f $output_dir/$pred_prog/genome/"$organism"_db.pin ]
then
  makeblastdb -in $output_dir/$pred_prog/genome/$organism.faa -dbtype prot -out $output_dir/$pred_prog/genome/"$organism"_db
fi

if [ ! -f $output_dir/$pred_prog/feature/"$organism".feature.ngene ]
then
  if [ ! -d $output_dir/$pred_prog/BLASTss ]
  then
    mkdir $output_dir/$pred_prog/BLASTss
    mkdir $output_dir/$pred_prog/BLASTno
    mkdir $output_dir/$pred_prog/BLASTff
    mkdir $output_dir/$pred_prog/BLASTcogn
  fi
  # it may take a long time
  if [ ! -f $output_dir/$pred_prog/BLASTss/QuerySelf.tab ]
  then
    psiblast -query $output_dir/$pred_prog/genome/$organism.faa -db $output_dir/$pred_prog/genome/"$organism"_db -show_gis -outfmt 7 -num_alignments 10 -dbsize 100000000 -comp_based_stats F -seg no -out $output_dir/$pred_prog/BLASTss/QuerySelf.tab -num_threads $num_threads
  fi

  if [ ! -f $output_dir/$pred_prog/BLASTno/QueryCOGs.tab ]
  then
    psiblast -query $output_dir/$pred_prog/genome/$organism.faa -db $prog_dir/db/COG/COGs -show_gis -outfmt 7 -num_alignments 1000 -dbsize 100000000 -comp_based_stats F -seg no -out $output_dir/$pred_prog/BLASTno/QueryCOGs.tab -num_threads $num_threads
  fi

  if [ ! -f $output_dir/$pred_prog/BLASTff/QueryCOGs.tab ]
  then
    psiblast -query $output_dir/$pred_prog/genome/$organism.faa -db $prog_dir/db/COG/COGs -show_gis -outfmt 7 -num_alignments 1000 -dbsize 100000000 -comp_based_stats T -seg yes -out $output_dir/$pred_prog/BLASTff/QueryCOGs.tab -num_threads $num_threads
  fi

  echo "$organism" > "$output_dir/$pred_prog/genome/$organism".genomeid
  protnum=`less  "$output_dir/$pred_prog/genome/$organism".gene_id | wc -l`
  awk -v size="$protnum" '{for(i=0;i<size;i++) print}' "$output_dir/$pred_prog/genome/$organism".genomeid > "$output_dir/$pred_prog/genome/$organism".gid
  rm "$output_dir/$pred_prog/genome/$organism".genomeid

  if [ "$pred_prog" == "ncbi" ]
  then
    less $output_dir/$pred_prog/genome/$organism.gene_id | cut -d'|' -f2 > $output_dir/$pred_prog/genome/$organism.protid
    paste -d',' $output_dir/$pred_prog/genome/$organism.protid $output_dir/$pred_prog/genome/$organism.gid > $output_dir/$pred_prog/genome/$organism.p2o.csv
  else # For genes predicted by Prodigal
    paste -d',' $output_dir/$pred_prog/genome/$organism.gene_id $output_dir/$pred_prog/genome/$organism.gid > $output_dir/$pred_prog/genome/$organism.p2o.csv
    sed -i 's/>//g' $output_dir/$pred_prog/genome/$organism.p2o.csv
  fi

  cat $output_dir/$pred_prog/genome/"$organism".p2o.csv $prog_dir/db/COG/COG.p2o.csv > $output_dir/$pred_prog/genome/tmp.p2o.csv

  $prog_dir/bin/COGSoft/COGmakehash/COGmakehash -i=$output_dir/$pred_prog/genome/tmp.p2o.csv -o=$output_dir/$pred_prog/BLASTcogn -s="," -n=1

  $prog_dir/bin/COGSoft/COGreadblast/COGreadblast -d=$output_dir/$pred_prog/BLASTcogn -u=$output_dir/$pred_prog/BLASTno -f=$output_dir/$pred_prog/BLASTff -s=$output_dir/$pred_prog/BLASTss -e=0.1 -q=2 -t=2

  $prog_dir/bin/COGSoft/COGcognitor/COGcognitor -i=$output_dir/$pred_prog/BLASTcogn -t=$prog_dir/db/COG/cog2003-2014.csv -q=$output_dir/$pred_prog/genome/"$organism".p2o.csv -o=$output_dir/$pred_prog/genome/"$organism".COG.csv

  if [ "$pred_prog" == "ncbi" ]
  then
    python $prog_dir/feature/parse_COG.py -g $output_dir/$pred_prog/genome/"$organism".protid  -c $output_dir/$pred_prog/genome/"$organism".COG.csv  -o  $output_dir/$pred_prog/feature/"$organism".feature.ngene
  else # For genes predicted by Prodigal
    python $prog_dir/feature/parse_COG.py -g $output_dir/$pred_prog/genome/"$organism".gene_id  -c $output_dir/$pred_prog/genome/"$organism".COG.csv  -o  $output_dir/$pred_prog/feature/"$organism".feature.ngene
  fi
fi

fi # for a complete genome

############################ merge files ################################################
echo "##########################################"
echo "Merging features for each gene/ORF"
paste $output_dir/$pred_prog/genome/$organism.glist $output_dir/$pred_prog/feature/$organism.feature.gc $output_dir/$pred_prog/feature/$organism.feature.codon $output_dir/$pred_prog/feature/$organism.feature.kmer $output_dir/$pred_prog/feature/$organism.feature.mobgene $output_dir/$pred_prog/feature/$organism.feature.phage $output_dir/$pred_prog/feature/$organism.feature.vfdb $output_dir/$pred_prog/feature/$organism.feature.AR  $output_dir/$pred_prog/feature/$organism.feature.ngene > $output_dir/$pred_prog/feature/$organism.feature.multi


############################ Analyzing boundary features ################################################
echo "##########################################"
echo "Predicting tRNA"
if [ ! -d $output_dir/boundary ]
then
  mkdir $output_dir/boundary
fi

if [ ! -f $output_dir/boundary/$organism.trna_pred ]
then
  # tRNA is only dependant on the original genome sequence
  tRNAscan-SE -B --frag $output_dir/boundary/$organism.trna_frag -o $output_dir/boundary/$organism.trna_pred -m $output_dir/boundary/$organism.trna_stat --brief $output_dir/$organism.fna
fi
less $output_dir/boundary/$organism.trna_pred | cut -f1-4 > $output_dir/boundary/$organism.pred_trna


if [ $mode == 0 ]  # Predict repeats only for finished complete genomes
then
echo "##########################################"
echo "Finding repeats"
if [ ! -f $output_dir/boundary/$organism.repseek ]
then
# repeat is only dependant on the original genome sequence
repseek -l 15 -O 0 -r $output_dir/boundary/$organism.repseek $output_dir/$organism.fna
fi
fi
