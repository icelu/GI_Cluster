#!/bin/bash

# Script for segmenting a complete genome or a set of contigs
#
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
#

########## Usage ##########
# ./GI_Segmentation.sh -s  $prog_dir -o $output_dir -n $organism -m $seg_prog -p prodigal
# e.g. ./GI_Segmentation.sh -s ./GIFilter -o ./research/data/species/cft73 -n NC_004431  -m ./research/software/HGT/mjsd -p prodigal


software=$(basename $0)

function usage() {
  echo -e "GI_Segmentation: segmenting a complete genome or a set of contigs"
  echo "Version 1.0
Usage: $software [options] -s [the directory containing all the scripts] -o [the output directory]
-n [the name of the organism (NCBI Accession, e.g. NC_003198)] -m [programs for genome segmation (e.g. mjsd, gcprofile, gisvm, alienhunter)] -p [programs for gene prediction (e.g. prodigal, ncbi)]

OPTIONS	Default	DESCIPTION
-b	0	: mode of running: 0 for complete genome, 1 for incomplete genome (contigs).
-w 5000: the window size.
-j 5000: the step size.
-h 	----	: print this help.
  "
  exit -1
}

width=5000
step=5000

while getopts "b:s:o:n:m:p:w:j:h" OPT; do
  case $OPT in
    b) mode=$OPTARG || exit 1;;

    s) prog_dir=$OPTARG || exit 1;;
    o) output_dir=$OPTARG || exit 1;;
    n) organism=$OPTARG || exit 1;;
    m) seg_prog=$OPTARG || exit 1;;
    p) pred_prog=$OPTARG || exit 1;;

    w) width=$OPTARG || exit 1;;
    j) step=$OPTARG || exit 1;;

    h) usage && exit;;
  esac
done



if [ ! -d $output_dir/$seg_prog ]
then
  mkdir -p $output_dir/$seg_prog
fi
p0="window"
p1="mjsd"
p2="gcprofile"
p3="gisvm"
p4="alienhunter"

if [ ! -f $output_dir/$seg_prog/$organism."$seg_prog" ]
then
  if [ "$seg_prog" == "$p0" ]
  then
    if [ $mode == 0 ]  # For finished complete genomes
    then
      python $prog_dir/segmentation/generate_segs.py -w $width -s $step -i $output_dir/$organism.fna -o $output_dir/$seg_prog/$organism."$seg_prog"."$width"."$step"
    else  # For contigs
      python $prog_dir/segmentation/generate_segs.py -w $width -s $step -c -i $output_dir/$organism.fna -o $output_dir/$seg_prog/$organism."$seg_prog"."$width"."$step"
    fi
  fi

  if [ $mode == 0 ]  # The following programs are only suitable to finished complete genomes
  then
    if [ "$seg_prog" == "$p1" ]
    then
      $prog_dir/bin/mjsd/so_jensen -f $output_dir/$organism.fna -s 0.1 -o 2 > $output_dir/$seg_prog/$organism.$seg_prog.orig
      python $prog_dir/segmentation/parse_MJSD_segs.py -f $output_dir/$organism.fna -m $width -s $output_dir/$seg_prog/$organism.$seg_prog.orig   > $output_dir/$seg_prog/$organism."$seg_prog"
    fi

    if [ "$seg_prog" == "$p2" ]
    then
      python $prog_dir/segmentation/parse_gcprofile_segs.py -t 10 -m $width -f $output_dir/$organism.fna -p $prog_dir/bin/gcprofile/linux/ -o $output_dir/$seg_prog/$organism."$seg_prog"
    fi

    if [ "$seg_prog" == "$p3" ]
    then
      python $prog_dir/bin/gisvm/GI_SVM.py -N 0.9 -t 2 -k 1 -K 8 -d $output_dir/$seg_prog -c 5 -a -q $output_dir/$organism.fna
      mv $output_dir/$seg_prog/mergedRes_auto_*_5000_2500 $output_dir/$seg_prog/$organism."$seg_prog"
    fi

    if [ "$seg_prog" == "$p4" ]
    then
      if [ ! -f $output_dir/$seg_prog/$organism."$seg_prog".orig ]
      then
        $prog_dir/bin/alienhunter/alien_hunter $output_dir/$organism.fna  $output_dir/$seg_prog/$organism."$seg_prog".orig
      fi
      less $output_dir/$seg_prog/$organism."$seg_prog".orig | grep 'misc_feature'| sed 's/\s\+/ /g' | cut -d' ' -f3 > $output_dir/$seg_prog/$organism."$seg_prog"
      sed -i 's/\.\./\t/g' $output_dir/$seg_prog/$organism."$seg_prog"
    fi

  fi
fi

cp $output_dir/$seg_prog/"$organism"."$seg_prog"."$width"."$step" $output_dir/$seg_prog/$organism.segment
