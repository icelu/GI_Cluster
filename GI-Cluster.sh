# Script for GI-Cluster: detecting genomic islands in newly sequenced microbial genomes by consensus clustering on multiple features
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg

########## Input ##########
# The genome sequence of a newly sequenced bacterium: fasta file format (NCBI file: $organism.fna)
# OPTIONAL INPUT: $organism.ffn, $organism.faa

########## usage ##########
# ./GI-Cluster.sh $prog_dir $output_dir $organism -b 0
# e.g. ./GI-Cluster.sh ./GIFilter ./research/data/species/cft73 NC_004431

####################################################
software=$(basename $0)

function usage() {
  echo -e "GI-Cluster: detecting genomic islands in newly sequenced microbial genomes by consensus clustering on multiple features"
  echo "Version 1.0
Usage: $software [options] -s [the directory containing all the scripts] -o [the output directory]
-n [the name of the organism (NCBI Accession, e.g. NC_003198)] -m [programs for genome segmation (e.g. window, mjsd, gcprofile, gisvm, alienhunter)] -p [programs for gene prediction (e.g. prodigal, ncbi, none)]

OPTIONS	Default	DESCIPTION
-b	0	: mode of running: 0 for complete genome, 1 for contigs without gene predictions, 2 for contigs with gene predictions.

-e	1e-5	: e-value used during identification of phage-related genes, i.e., blastp against PHAST.
-r	1e-5	: e-value used during identification of virulence factors, i.e., blastp against VFDB.
-a	1e-5	: e-value used during identification of antibiotic resistance genes, i.e., blastp against CARD.

-d	4	: number of threads used by blast.
-u	16	: number of CPUs used by cmsearch.

-g	2500	: gaps allowed when merging adjacent segments.

-f 	1	: show the visualization of features related to genomic island
-h ---- : print this help
-v ---- : Version 1.0
  "
  exit -1
}

# -t	1e-5	: e-value used during identification of mobgenes, i.e., hmmsearch against pfam database of mobgenes.
pred_prog=none
# Set default values for optional parameters
phage_evalue=1e-5
virdb_evalue=1e-5
arg_evalue=1e-5

num_threads=4
num_cpus=16

gap=5000

mode=0

show_figure=1


while getopts "b:s:o:n:m:p:g:c:e:r:a:d:u:f:h" OPT; do
  case $OPT in
    b) mode=$OPTARG || exit 1;;

    s) prog_dir=$OPTARG || exit 1;;
    o) output_dir=$OPTARG || exit 1;;
    n) organism=$OPTARG || exit 1;;
    m) seg_prog=$OPTARG || exit 1;;
    p) pred_prog=$OPTARG || exit 1;;

    g) gap=$OPTARG || exit 1;;

    #t) mobgene_evalue=$OPTARG || exit 1;;
    e) phage_evalue=$OPTARG || exit 1;;
    r) virdb_evalue=$OPTARG || exit 1;;
    a) arg_evalue=$OPTARG || exit 1;;

    d) num_threads=$OPTARG || exit 1;;
    u) num_cpus=$OPTARG || exit 1;;

    f) show_figure=$OPTARG || exit 1;;

    h) usage && exit;;
  esac
done


echo "##########################################"
if [ $mode == 0 ]
then
echo "Running GI-Cluster for a complete genome"
fi
if [ $mode == 1 ]
then
echo "Running GI-Cluster for contigs (without gene predictions)"
fi
if [ $mode == 2 ]
then
echo "Running GI-Cluster for contigs (with gene predictions)"
fi


if [ $mode != 1 ]
then
############## Predict genes from raw genome sequence ##################
# 3 types of required output files: fna (DNA sequence), ffn (gene sequence), faa (protein sequence)
  sh $prog_dir/Gene_Prediction.sh -s $prog_dir -o $output_dir -n $organism -m $seg_prog -p $pred_prog -b $mode

############## Compute known features related to genomic islands in the unit of genes #########################
  if [ "$pred_prog" != "none" ]
  then
  sh $prog_dir/GI_Feature.sh -s $prog_dir -o $output_dir -n $organism -m $seg_prog -p $pred_prog -b $mode -e $phage_evalue -r $virdb_evalue -a $arg_evalue -d $num_threads -u $num_cpus
  wait
  fi
fi # for complete genomes


############################## genome segmentation ####################################
echo "##########################################"
echo "Getting genome segements"
if [ ! -f $output_dir/$seg_prog/$organism.segment ]
then
  sh $prog_dir/GI_Segmentation.sh -s $prog_dir -o $output_dir -n $organism -m $seg_prog -p $pred_prog -b $mode
  wait
fi

#################################### summarize features of each genomic region ##############################
echo "##########################################"
echo "Getting features of each genome segement"
if [ ! -f $output_dir/$pred_prog/$seg_prog/$organism.feature.multi.percentage ]
then
  sh $prog_dir/Segment_Feature.sh -s $prog_dir -o $output_dir -n $organism -m $seg_prog -p $pred_prog -b $mode
  wait
fi

# Get the labels for each segment when there are references
# have to run this after changing references
# python $prog_dir/IntervalIntersection.py --interval1 $output_dir/$seg_prog/$organism.segment --interval2 $output_dir/$organism.refgi -o $output_dir/$seg_prog/"$seg_prog"_refgi -f 0
#
# python $prog_dir/assign_cluster.py --interval1 $output_dir/$seg_prog/feature/$organism.feature.multi.percentage --interval2 $output_dir/$seg_prog/"$seg_prog"_refgi_interval1 -o $output_dir/$seg_prog/feature/$organism.feature.multi.percentage.labeled
# python $prog_dir/assign_cluster.py -c --interval1 $output_dir/$seg_prog/feature/$organism.feature.multi.percentage.labeled --interval2 $output_dir/$seg_prog/"$seg_prog"_pos_overlap_interval1 --interval3 $output_dir/$seg_prog/"$seg_prog"_neg_overlap_interval1 -o $output_dir/$seg_prog/feature/$organism.feature.multi.percentage.labeled2


###################### Run Rscript for consensus clustering ##############################
echo "##########################################"
echo "Running consensus clustering on the feature matrix"
# TODO: add parameters to input
# remember to put the separator at the end of output_dir
pFeature=1
feature=comp_content
method=average
rep=1
postprocess=true
if [ $mode == 1 ] # Not dependant on gene predictions
then
  # TODO: make R code more flexible
  echo
else  # Dependant on gene predictions
  if [ ! -f $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/"$organism"_GI ]
  then
    mkdir -p $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/
    ffile=$output_dir/$pred_prog/$seg_prog/$organism.feature.multi.percentage
    nohup Rscript $prog_dir/GI_Clustering.R -f $ffile -o $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/ -a "$organism" -l $prog_dir/clustering -k 2 -K 3 -r $rep -e $pFeature -s $seg_prog -v true -C hclust -P method=$method -m true -d $feature -S $postprocess 2>&1 > $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/rscript_std
  fi
fi

###################### Boundary refinement ##############################
echo "##########################################"
echo "Refining the boundary of predicted GIs"
# TODO: add parameters to input
dist=1000
# No need to distinguish the mode.
python $prog_dir/boundary/refine_boundary.py -i $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/"$organism"_GI -r $output_dir/$seg_prog/feature/$organism.boundary.repeat -t $output_dir/$seg_prog/feature/$organism.boundary.trna -o  $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/"$organism"_refined_GI -c $dist


# gap=2500
# post process GI candidates
if [ $mode == 0 ] # For complete genomes
then
  python $prog_dir/boundary/postprocess_segments.py -l $gap -i $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/"$organism"_refined_GI -a -g $output_dir/$pred_prog/genome/$organism.glist
else  # TODO: utilize gene predictions
  python $prog_dir/boundary/postprocess_segments.py -c -l $gap -i $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/"$organism"_refined_GI -g $output_dir/$pred_prog/genome/$organism.glist
fi

# python $prog_dir/boundary/convert_contig_GIs.py -i $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/merged_"$organism"_refined_GI -o $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/reformated_merged_"$organism"_refined_GI

##################### Get the features of predicted GIs #####################
echo "##########################################"
echo "Getting the features of predicted GIs"
gifile=$output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/merged_"$organism"_refined_GI

if [ $mode == 2 ] # For contigs with gene predictions
then
  python $prog_dir/feature/mergeFeature.py -g $output_dir/$pred_prog/feature/$organism.feature.multi -o $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/$organism.gi.feature -s $gifile -r $output_dir/$pred_prog/feature/$organism.cmscan.tbl -c -m $output_dir/$organism.fna -d $output_dir/$pred_prog/genome/$organism.gene_id

  python $prog_dir/feature/countFeature.py -a -i $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/$organism.gi.feature -o $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/$organism.gi.feature.percentage -r $output_dir/$seg_prog/feature/$organism.boundary.repeat -t $output_dir/$seg_prog/feature/$organism.boundary.trna -k $output_dir/$seg_prog/feature/$organism.kmer.covariance
fi

if [ $mode == 1 ] # TODO: For contigs without gene predictions
then
  python $prog_dir/feature/mergeFeature.py -g $output_dir/$pred_prog/feature/$organism.feature.multi -o $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/$organism.gi.feature  -s $gifile -r $output_dir/$pred_prog/feature/$organism.cmscan.tbl
fi

if [ $mode == 0 ] # For compute genome with gene predictions
then
  python $prog_dir/feature/mergeFeature.py -g $output_dir/$pred_prog/feature/$organism.feature.multi -o $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/$organism.gi.feature  -s $gifile -r $output_dir/$pred_prog/feature/$organism.cmscan.tbl

  # Form the feature matrix for clustering. For kmer,  use the value for each segment directly
  python $prog_dir/feature/countFeature.py -a -i $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/$organism.gi.feature -o $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/$organism.gi.feature.percentage -r $output_dir/$seg_prog/feature/$organism.boundary.repeat -t $output_dir/$seg_prog/feature/$organism.boundary.trna -k $output_dir/$seg_prog/feature/$organism.kmer.covariance
fi

##################### Visualization  ##############################
if [ $mode == 0 ] # Only support complete genomes
  then
  echo "##########################################"
  echo "Create figures to visualize the predicted genomic islands and related features"
  # sh $prog_dir/GI_Visualization.sh -s $prog_dir -o $output_dir -n $organism -m $seg_prog -p $pred_prog
  show=1
  if [ "$show_figure" -eq "$show" ]
  then
    if [ ! -d $output_dir/$seg_prog/feature/visualization/data ]
    then
      mkdir -p $output_dir/$seg_prog/feature/visualization/data
      mkdir -p $output_dir/$seg_prog/feature/visualization/etc
      mkdir -p $output_dir/$seg_prog/feature/visualization/img
    fi
    # Copy template files
    cp $prog_dir/visualization/etc/* $output_dir/$seg_prog/feature/visualization/etc

    python $prog_dir/visualization/prepare_intervals.py -g $output_dir/$organism.fna -i $gifile -o $output_dir/$seg_prog/feature/visualization/data/$organism.gi -c $output_dir/$seg_prog/feature/visualization/data/$organism.chr -f $output_dir/$seg_prog/feature/visualization/data/$organism.highlight
    # Use features of all the input segments
    python $prog_dir/visualization/prepare_features.py -a -i $output_dir/$pred_prog/$seg_prog/$organism.feature.multi.percentage -o $output_dir/$seg_prog/feature/visualization/data

    prefix="feature.multi.percentage"
    python $prog_dir/visualization/prepare_feature_config.py -i $output_dir/$seg_prog/feature/visualization/etc/circos.gifeature.conf.template -n $organism -o $output_dir/$seg_prog/feature/visualization/etc/circos.gifeature.conf -p $prefix -f "$organism"_gifeature

    # Run circos at etc folder since the output directory is relative
    cd $output_dir/$seg_prog/feature/visualization/etc
    circos -conf $output_dir/$seg_prog/feature/visualization/etc/circos.gifeature.conf
    cd $output_dir
  fi
fi

if [ $comparison == 1 ] # Compare with the predictions from other programs
  then
  echo "##########################################"
  echo "Create figures to visualize the predicted genomic islands of different GI prediction methods"
  # sh $prog_dir/GI_Visualization.sh -s $prog_dir -o $output_dir -n $organism -m $seg_prog -p $pred_prog
  sh $prog_dir/GI_Comparison.sh -s $prog_dir -o $output_dir -n $organism -m $seg_prog -p $pred_prog -b $mode
fi
