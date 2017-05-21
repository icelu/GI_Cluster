########## usage ##########
# ./GI_Visualization.sh $prog_dir $output_dir $organism
# e.g. ./GI_Visualization.sh $output_dir/GIFilter /home/ice/vmshare/research/data/species/cft73 NC_004431


software=$(basename $0)

function usage() {
  echo -e "GI_Comparison: Visualizing the results of different GI prediction methods"
  echo "Version 0.1
Usage: $software -s [the directory containing all the scripts] -o [the output directory]
-n [the name of the organism (NCBI Accession, e.g. NC_003198)]
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
cmp=1
if [ "$show_figure" -eq "$cmp" ]
then
  if [ ! -d $output_dir/visualization/data ]
  then
    mkdir -p $output_dir/visualization/data
    mkdir -p $output_dir/visualization/etc
    mkdir -p $output_dir/visualization/img
  fi
  # copy template files
  cp $prog_dir/visualization/etc/* $output_dir/visualization/etc

  # Add the name of the programs here
  prog1=gihunter
  prog2=merged_islandviewer4
  prog3=gisvm
  prog4=gicluster
  prog4=gicluster_gisvm

  ## declare an array variable
  declare -a arr=("prog1" "prog2" "prog3" "prog4")

  ## now loop through the above array
  for prog in "${arr[@]}"
  do
     echo "$prog"
     python $prog_dir/visualization/prepare_intervals.py -g $output_dir/$organism.fna -i $output_dir/$prog.gilist -o $output_dir/visualization/data/$prog.gilist
  done

  python $prog_dir/visualization/prepare_interval_config.py -i $output_dir/visualization/etc/circos.gi.conf.template -n $organism -o $output_dir/visualization/etc/circos.gi.conf


  # run circos at etc folder since the output directory is relative
  cd $output_dir/visualization/etc
  circos -conf $output_dir/visualization/etc/circos.gifeature.conf
  cd $output_dir
fi
