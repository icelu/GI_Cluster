#!/bin/bash

# Script for predicting genes from (complete/incomplete) genome sequence
#
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
#

########## Usage ##########
# ./Gene_Prediction.sh -s  $prog_dir -o $output_dir -n $organism -m $seg_prog  -p prodigal -b $mode
# e.g. ./Gene_Prediction.sh -s ./GIFilter -o ./research/data/species/cft73 -n NC_004431  -m ./research/software/HGT/mjsd -p prodigal -b 0


software=$(basename $0)

function usage() {
  echo -e "Gene_Prediction: predicting genes from (complete/incomplete) genome sequence"
  echo "Version 1.0
Usage: $software [options] -s [the directory containing all the scripts] -o [the output directory]
-n [the name of the organism (e.g. NC_003198)] -m [programs for genome segmation (e.g. mjsd, gcprofile, gisvm, alienhunter)] -p [programs for gene prediction (e.g. prodigal, ncbi, ncbi_old)]

OPTIONS	Default	DESCIPTION
-b	0	: mode of running: 0 for complete genome, 1 for incomplete genome (contigs).

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
      prodigal -m -i $output_dir/$organism.fna -o $output_dir/$pred_prog/genome/"$organism".gbk -a $output_dir/$pred_prog/genome/"$organism".faa -d $output_dir/$pred_prog/genome/"$organism".ffn
    fi
    if [ $mode == 1 ]
    then
      prodigal -p meta -m -i $output_dir/$organism.fna -o $output_dir/$pred_prog/genome/"$organism".gbk -a $output_dir/$pred_prog/genome/"$organism".faa -d $output_dir/$pred_prog/genome/"$organism".ffn
    fi
  fi
fi

if [ "$pred_prog" == "ncbi_old" ]
then
  echo "##########################################"
  echo "Getting gene predictions from NCBI old files"
  if [ ! -f $output_dir/$pred_prog/genome/"$organism".ffn ]
  then
    cp $output_dir/"$organism".ffn $output_dir/$pred_prog/genome/"$organism".ffn
    cp $output_dir/"$organism".faa $output_dir/$pred_prog/genome/"$organism".faa
  fi
fi

if [ "$pred_prog" == "ncbi" ]
then
  echo "##########################################"
  echo "Getting gene predictions from NCBI"
  if [ ! -f $output_dir/$pred_prog/genome/"$organism".ffn ]
  then
    cp $output_dir/"$organism"_cds_from_genomic.fna $output_dir/$pred_prog/genome/"$organism"_cds_from_genomic.fna
    cp $output_dir/"$organism"_protein.faa $output_dir/$pred_prog/genome/"$organism".faa
  fi
fi

if [ "$pred_prog" == "custom" ]
then
  echo "##########################################"
  echo "Getting gene predictions from custom predictions"
  if [ ! -f $output_dir/$pred_prog/genome/"$organism".ffn ]
  then
    cp $output_dir/"$organism".ffn $output_dir/$pred_prog/genome/"$organism".ffn
    cp $output_dir/"$organism".faa $output_dir/$pred_prog/genome/"$organism".faa
    cp $output_dir/"$organism".glist $output_dir/$pred_prog/genome/$organism.glist 
  fi
fi

# Get gene/orf list: name, loc, strand
if [ ! -f $output_dir/$pred_prog/genome/$organism.gene_id ]
then
  # File $organism.gene_id is used later when extracting features
  less $output_dir/$pred_prog/genome/"$organism".faa | grep '^>'  | cut -d' ' -f1 > $output_dir/$pred_prog/genome/$organism.gene_id
  # $output_dir/$pred_prog/genome/$organism.glist -- ID, start, end, strand. This file is mainly used for the feature table of genes
  if [ "$pred_prog" == "prodigal" ]
  then
    less $output_dir/$pred_prog/genome/"$organism".ffn | grep '^>'  | cut -d'#' -f2-4  > $output_dir/$pred_prog/genome/$organism.gene_locus
    python $prog_dir/feature/parse_geneloc.py -l $output_dir/$pred_prog/genome/$organism.gene_locus  -o $output_dir/$pred_prog/genome/$organism.glist
  fi
  if [ "$pred_prog" == "ncbi_old" ]
  then
    # For NCBI annotation files (old version)
    # e.g. >gi|33864539|ref|NC_005070.1|:174-1331 Synechococcus sp. WH 8102, complete genome
    less $output_dir/$pred_prog/genome/$organism.gene_id | cut -d'|' -f2 > $output_dir/$pred_prog/genome/$organism.protid
    less $output_dir/$pred_prog/genome/"$organism".ffn | grep '^>'  | cut -d':' -f2 | cut -d' ' -f1 > $output_dir/$pred_prog/genome/$organism.gene_locus
    python $prog_dir/feature/parse_geneloc.py -n -l $output_dir/$pred_prog/genome/$organism.gene_locus  -o $output_dir/$pred_prog/genome/$organism.glist
  fi
  if [ "$pred_prog" == "ncbi" ]
  then
    # For NCBI annotation files (new version)
    less $output_dir/$pred_prog/genome/"$organism".faa | grep '^>'  | cut -d'>' -f2 | cut -d' ' -f1 > $output_dir/$pred_prog/genome/$organism.protid
    python $prog_dir/feature/parse_ncbi.py -i $output_dir/$pred_prog/genome/"$organism"_cds_from_genomic.fna -p $output_dir/$pred_prog/genome/$organism.protid -d $output_dir/$pred_prog/genome/$organism.gene_id -l $output_dir/$pred_prog/genome/$organism.glist -g $output_dir/$pred_prog/genome/$organism.ffn
  fi 
fi
