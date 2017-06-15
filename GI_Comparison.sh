# Script for GI_Comparison: Visualizing predicted GIs by different methods for a microbial genome in a Circos plot
#
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
#
########## Usage ##########
# ./GI_Comparison.sh -s $prog_dir -o $output_dir -n $organism
# e.g. sh /home/b/bingxin/GI-Cluster/GI_Comparison.sh  -s /home/b/bingxin/GI-Cluster -o /home/b/bingxin/genome/StypiCT18 -n NC_003198


software=$(basename $0)

function usage() {
  echo -e "GI_Comparison: Visualizing the results of different GI prediction methods"
  echo "Version 1.0
Usage: $software -s [the directory containing all the scripts] -o [the output directory] -n [the name of the organism (NCBI Accession, e.g. NC_003198)]

OPTIONS	Default	DESCIPTION
-h 	----	: print this help
  "
  exit -1
}


while getopts "s:o:n:m:p:h" OPT; do
  case $OPT in
    s) prog_dir=$OPTARG || exit 1;;
    o) output_dir=$OPTARG || exit 1;;
    n) organism=$OPTARG || exit 1;;

    h) usage && exit;;
  esac
done


###################### visualization  ##############################
# required input file:
# $prog.gilist -- predicted intervals from different programs
echo "##########################################"
echo "Visualizing the results of different GI prediction methods"

if [ ! -d $output_dir/visualization/data ]
then
  mkdir -p $output_dir/visualization/data
  mkdir -p $output_dir/visualization/etc
  mkdir -p $output_dir/visualization/img
fi
# copy template files
cp $prog_dir/visualization/etc/* $output_dir/visualization/etc

# Add the name of the programs here
prog1=refgi
prog2=gicluster
prog3=gihunter
prog4=merged_islandviewer4
prog5=gicluster_gisvm
prog6=gisvm

# declare an array variable
declare -a arr=("$prog1" "$prog2" "$prog3" "$prog4" "$prog5" "$prog6")

# now loop through the above array
for prog in "${arr[@]}"
do
  echo "$prog"
  # To make sure all the input files are available
  if [ $prog == "refgi" ]
  then
    echo "Creating highlight file"
    python $prog_dir/visualization/prepare_intervals.py -g $output_dir/$organism.fna -i $output_dir/comparison/$prog.gilist -o $output_dir/visualization/data/$prog.gilist -f $output_dir/visualization/data/$organism.refgi.highlight
  else
    python $prog_dir/visualization/prepare_intervals.py -g $output_dir/$organism.fna -i $output_dir/comparison/$prog.gilist -o $output_dir/visualization/data/$prog.gilist
  fi
done

python $prog_dir/visualization/prepare_interval_config.py -i $output_dir/visualization/etc/circos.gi.conf.template -o $output_dir/visualization/etc/circos.gi.conf -p "$prog1","$prog2","$prog3","$prog4","$prog5","$prog6" -f "$organism"_gi_cmp -n $organism

# run circos at etc folder since the output directory is relative
cd $output_dir/visualization/etc
circos -conf $output_dir/visualization/etc/circos.gi.conf
cd $output_dir
