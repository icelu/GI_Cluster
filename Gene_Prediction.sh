# Script for predicting genes from (complete/draft) genome sequence
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg

########## usage ##########
# ./Gene_Prediction.sh -s  $prog_dir -o $output_dir -n $organism -m $seg_prog -b $mode
# e.g. ./Gene_Prediction.sh -s ./GIFilter -o ./research/data/species/cft73 -n NC_004431  -m ./research/software/HGT/mjsd -b 0


software=$(basename $0)

function usage() {
  echo -e "Gene_Prediction: predicting genes from (complete/draft) genome sequence"
  echo "Version 1.0
Usage: $software [options] -s [the directory containing all the scripts] -o [the output directory]
-n [the name of the organism (e.g. NC_003198)] -m [programs for genome segmation (e.g. mjsd, gcprofile, gisvm, alienhunter)] -p [programs for gene prediction (e.g. prodigal, ncbi)]

OPTIONS	Default	DESCIPTION
-b	0	: mode of running: 0 for complete genome, 1 for contigs without gene predictions, 2 for contigs with gene predictions.

-h 	----	: print this help
  "
  exit -1
}


while getopts "b:s:o:n:m:p:h" OPT; do
  case $OPT in
    b) mode=$OPTARG || exit 1;;

    s) prog_dir=$OPTARG || exit 1;;
    o) output_dir=$OPTARG || exit 1;;
    n) organism=$OPTARG || exit 1;;
    m) seg_prog=$OPTARG || exit 1;;
    p) pred_prog=$OPTARG || exit 1;;

    h) usage && exit;;
  esac
done


if [ ! -d $output_dir/$pred_prog/genome/ ]
then
  mkdir -p $output_dir/$pred_prog/genome/
fi

if [ "$pred_prog" == "prodigal" ]
then
  echo "##########################################"
  echo "Predicting genes with prodigal"
  if [ ! -f $output_dir/$pred_prog/genome/"$organism".ffn ]
  then
    if [ $mode == 0 ]
    then
      $prog_dir/bin/prodigal -m -i $output_dir/$organism.fna -o $output_dir/$pred_prog/genome/"$organism".gbk -a $output_dir/$pred_prog/genome/"$organism".faa -d $output_dir/$pred_prog/genome/"$organism".ffn
    fi
    if [ $mode == 2 ]
    then
      $prog_dir/bin/prodigal -p meta -m -i $output_dir/$organism.fna -o $output_dir/$pred_prog/genome/"$organism".gbk -a $output_dir/$pred_prog/genome/"$organism".faa -d $output_dir/$pred_prog/genome/"$organism".ffn
    fi
  fi
fi

if [ "$pred_prog" == "ncbi" ]
then
  echo "##########################################"
  echo "Getting gene predictions from NCBI"
  if [ ! -f $output_dir/$pred_prog/genome/"$organism".ffn ]
  then
    cp $output_dir/"$organism".ffn $output_dir/$pred_prog/genome/"$organism".ffn
    cp $output_dir/"$organism".faa $output_dir/$pred_prog/genome/"$organism".faa
  fi
fi

# Get gene/orf list: name, loc, strand
if [ ! -f $output_dir/$pred_prog/genome/$organism.gene_id ]
then
  # File $organism.gene_id is used later when extracting features
  less $output_dir/$pred_prog/genome/"$organism".faa | grep '^>'  | cut -d' ' -f1 > $output_dir/$pred_prog/genome/$organism.gene_id
  # for each gene, a record in a row
  if [ "$pred_prog" == "prodigal" ]
  then
    less $output_dir/$pred_prog/genome/"$organism".ffn | grep '^>'  | cut -d'#' -f2-4  > $output_dir/$pred_prog/genome/$organism.gene_locus
    python $prog_dir/feature/parse_geneloc.py -l $output_dir/$pred_prog/genome/$organism.gene_locus  -o $output_dir/$pred_prog/genome/$organism.glist
  else
    # For   ncbi annotation files
    less $output_dir/$pred_prog/genome/"$organism".ffn | grep '^>'  | cut -d':' -f2 | cut -d' ' -f1 > $output_dir/$pred_prog/genome/$organism.gene_locus
    python $prog_dir/feature/parse_geneloc.py -n -l $output_dir/$pred_prog/genome/$organism.gene_locus  -o $output_dir/$pred_prog/genome/$organism.glist
  fi
fi
