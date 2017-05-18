########## usage ##########
# ./GI_Visualization.sh $prog_dir $output_dir $organism
# e.g. ./GI_Visualization.sh $output_dir/GIFilter /home/ice/vmshare/research/data/species/cft73 NC_004431

# Contact: bingxin@comp.nus.edu.sg

software=$(basename $0)

function usage() {
  echo -e "GI_Visualization: Visualizing the results of different GI prediction methods"
  echo "Version 0.1
Usage: $software -s [the directory containing all the scripts] -o [the output directory]
-n [the name of the organism (NCBI Accession, e.g. NC_003198)] -m [programs for genome segmation (e.g. mjsd, gcprofile, gisvm, alienhunter)] -p [programs for gene prediction (e.g. prodigal, ncbi)]
-h 	----	: print this help
  "
  exit -1
}


while getopts "s:o:n:m:p:h" OPT; do
  case $OPT in
    s) prog_dir=$OPTARG || exit 1;;
    o) output_dir=$OPTARG || exit 1;;
    n) organism=$OPTARG || exit 1;;
    m) seg_prog=$OPTARG || exit 1;;
    p) pred_prog=$OPTARG || exit 1;;

    h) usage && exit;;
  esac
done



###################### visualization  ##############################
# required input file:
# $output_dir/$pred_prog/$seg_prog/$organism.feature.multi.percentage -- feature table for each segments
# $output_dir/$pred_prog/$seg_prog/merged_GI -- final GI candidates
echo "##########################################"
echo "Visualizing the results of different GI prediction methods"
cmp=1
if [ "$show_figure" -eq "$cmp" ]
then
  if [ ! -d $output_dir/$pred_prog/$seg_prog/visualization/data ]
  then
    mkdir -p $output_dir/$pred_prog/$seg_prog/visualization/data
    mkdir -p $output_dir/$pred_prog/$seg_prog/visualization/etc
    mkdir -p $output_dir/$pred_prog/$seg_prog/visualization/img
  fi
  # copy template files
  cp $prog_dir/visualization/etc/* $output_dir/$pred_prog/$seg_prog/visualization/etc

  gifile=$output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/merged_"$organism"_refined_GI


  prog1=gihunter
  prog2=mergedRes_islandviewer3
  prog3=gisvm
  prog4=gifile

  ## declare an array variable
  declare -a arr=("prog1" "prog2" "prog3" "prog4")

  ## now loop through the above array
  for prog in "${arr[@]}"
  do
     echo "$prog"
     python $prog_dir/visualization/prepare_intervals.py -g $output_dir/$organism.fna -i $output_dir/$prog.gilist -o $output_dir/$pred_prog/$seg_prog/visualization/data/$prog.gilist
  done

  python $prog_dir/visualization/prepare_interval_config.py -i $output_dir/$seg_prog/feature/visualization/etc/circos.gi.conf.template -n $organism -o $output_dir/$seg_prog/feature/visualization/etc/circos.gi.conf


  # run circos at etc folder since the output directory is relative
  cd $output_dir/$pred_prog/$seg_prog/visualization/etc
  circos -conf $output_dir/$pred_prog/$seg_prog/visualization/etc/circos.gifeature.conf
  cd $output_dir
fi
