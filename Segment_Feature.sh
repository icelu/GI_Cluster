# Script for extracting features related to genomic islands in a genomic region
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg

########## Usage ##########
# ./Segment_Feature.sh -s  $prog_dir -o $output_dir -n $organism -m $seg_prog -b $mode
# e.g. ./Segment_Feature.sh -s ./GIFilter -o ./research/data/species/cft73 -n NC_004431  -m ./research/software/HGT/mjsd -b 0


software=$(basename $0)

function usage() {
  echo -e "Segment_Feature: Extracting features related to genomic islands in a genomic region"
  echo "Version 1.0
Usage: $software [options] -s [the directory containing all the scripts] -o [the output directory]
-n [the name of the organism (e.g. NC_003198)] -m [programs for genome segmation (e.g. mjsd, gcprofile, gisvm, alienhunter)] -p [programs for gene prediction (e.g. prodigal, ncbi)]

OPTIONS	Default	DESCIPTION
-b	0	: mode of running: 0 for complete genome, 1 for draft genome (contigs).
-t  1 : availablity of gene predictions: 1: with gene predictions, 0: without gene predictions.
-h 	----	: print this help
  "
  exit -1
}

mode=0
gene_prediction=1

while getopts "b:t:s:o:n:m:p:h" OPT; do
  case $OPT in
    b) mode=$OPTARG || exit 1;;
    t) gene_prediction=$OPTARG || exit 1;;

    s) prog_dir=$OPTARG || exit 1;;
    o) output_dir=$OPTARG || exit 1;;
    n) organism=$OPTARG || exit 1;;
    m) seg_prog=$OPTARG || exit 1;;
    p) pred_prog=$OPTARG || exit 1;;

    h) usage && exit;;
  esac
done

# For features not related to gene predictions
if [ ! -d $output_dir/$seg_prog/feature ]
then
  mkdir -p $output_dir/$seg_prog/feature
fi

# For features related to gene predictions
if [ "$pred_prog" != "none" ]
then
  if [ ! -d $output_dir/$pred_prog/$seg_prog ]
  then
    mkdir -p $output_dir/$pred_prog/$seg_prog
  fi
fi

echo "##########################################"
echo "Finding tRNA genes around each segment"
if [ ! -f $output_dir/$seg_prog/feature/$organism.boundary.trna ]
then
  if [ $mode == 0 ]  # For finished complete genomes
  then
    python $prog_dir/feature/find_trnas.py -i $output_dir/$seg_prog/$organism.segment -t $output_dir/boundary/$organism.pred_trna -o $output_dir/$seg_prog/feature/$organism.boundary.trna
  else  # For contigs
    python $prog_dir/feature/find_trnas.py -c -g $output_dir/$organism.fna -i $output_dir/$seg_prog/$organism.segment -t $output_dir/boundary/$organism.pred_trna -o $output_dir/$seg_prog/feature/$organism.boundary.trna
  fi
fi

if [ $mode == 0 ]  # Only for finished complete genomes, not dependant on gene predictions
then
  echo "##########################################"
  echo "Finding repeats around each segment"
  if [ ! -f $output_dir/$seg_prog/feature/$organism.boundary.repeat ]
  then
  # The output is based on the genomic segments
  python $prog_dir/feature/find_repeats.py -i $output_dir/$seg_prog/$organism.segment -r $output_dir/boundary/$organism.repseek -o $output_dir/$seg_prog/feature/$organism.boundary.repeat
  fi
fi


if [ $gene_prediction == 0 ]  # Only for sequences without gene predictions
then
  echo "##########################################"
  echo "Analyzing GC Content for each segment"
  if [ ! -f $output_dir/$seg_prog/feature/$organism.gc ]
  then
    if [ $mode == 0 ] # For complete genome
    then
      python $prog_dir/feature/analyze_GC.py -r -i $output_dir/$organism.fna -s $output_dir/$seg_prog/$organism.segment -o $output_dir/$seg_prog/feature/$organism.gc
    fi
    if [ $mode == 1 ] # For draft genome
    then
      python $prog_dir/feature/analyze_GC.py -c -i $output_dir/$organism.fna -s $output_dir/$seg_prog/$organism.segment -o $output_dir/$seg_prog/feature/$organism.gc
    fi
  fi
fi


echo "##########################################"
echo "Analyzing k-mer frequency for each segment"
if [ ! -f $output_dir/$seg_prog/feature/$organism.kmer.covariance ]
then
  if [ $mode == 0 ] # For complete genome
  then
    python $prog_dir/feature/analyze_kmer.py -r -i $output_dir/$organism.fna -s $output_dir/$seg_prog/$organism.segment -o $output_dir/$seg_prog/feature/$organism.kmer
  else # For contigs, not dependant on gene predictions
    python $prog_dir/feature/analyze_kmer.py -c -i $output_dir/$organism.fna -s $output_dir/$seg_prog/$organism.segment -o $output_dir/$seg_prog/feature/$organism.kmer
  fi
fi


echo "##########################################"
echo "Merging features for each segment"
segs=$output_dir/$seg_prog/$organism.segment
seg_feature1=$output_dir/$pred_prog/$seg_prog/$organism.seg_feature.multi
seg_feature1_percentage=$output_dir/$pred_prog/$seg_prog/$organism.feature.multi.percentage
seg_feature2=$output_dir/$seg_prog/feature/$organism.feature.multi.percentage
if [ $gene_prediction == 1 ] # With gene predictions
then
  if [ $mode == 1 ] # For contigs with gene predictions
  then
    python $prog_dir/feature/mergeFeature.py -g $output_dir/$pred_prog/feature/$organism.feature.multi -o $seg_feature1 -s $segs -r $output_dir/$pred_prog/feature/$organism.cmscan.tbl -c -m $output_dir/$organism.fna -d $output_dir/$pred_prog/genome/$organism.gene_id
  fi
  if [ $mode == 0 ] # For compute genome with gene predictions
  then
      python $prog_dir/feature/mergeFeature.py -g $output_dir/$pred_prog/feature/$organism.feature.multi -o $seg_feature1 -s $segs -r $output_dir/$pred_prog/feature/$organism.cmscan.tbl
  fi
  # Form the feature matrix for clustering. For kmer,  use the value for each segment directly. No need to distinguish mode
  python $prog_dir/feature/countFeature.py -a -i $seg_feature1 -o $seg_feature1_percentage -r $output_dir/$seg_prog/feature/$organism.boundary.repeat -t $output_dir/$seg_prog/feature/$organism.boundary.trna -k $output_dir/$seg_prog/feature/$organism.kmer.covariance
else  # Without gene predictions
  python $prog_dir/feature/countFeature.py -i $segs -o $seg_feature2 -r $output_dir/$seg_prog/feature/$organism.boundary.repeat -t $output_dir/$seg_prog/feature/$organism.boundary.trna -k  $output_dir/$seg_prog/feature/$organism.kmer.covariance -g $output_dir/$seg_prog/feature/$organism.gc
fi
