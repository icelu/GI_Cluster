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
-b	0	: mode of running: 0 for complete genome, 1 for draft genome (contigs).
-t  1 : availablity of gene predictions: 1: with gene predictions, 0: without gene predictions.

-e	1e-5	: e-value used during identification of phage-related genes, i.e., blastp against PHAST.
-r	1e-5	: e-value used during identification of virulence factors, i.e., blastp against VFDB.
-a	1e-5	: e-value used during identification of antibiotic resistance genes, i.e., blastp against CARD.

-d	4	: number of threads used by blast.
-u	16	: number of CPUs used by cmsearch.

-g	5000	: gaps allowed when merging adjacent segments.

-f 	1	: show the visualization of features related to genomic islands
-q 	0	: show the visualizations of genomic islands from different methods
-h ---- : print this help
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
gene_prediction=1

show_figure=1
comparison=0

while getopts "b:t:s:o:n:m:p:g:c:e:r:a:d:u:f:q:h" OPT; do
  case $OPT in
    b) mode=$OPTARG || exit 1;;
    t) gene_prediction=$OPTARG || exit 1;;

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
    q) comparison=$OPTARG || exit 1;;

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
echo "Running GI-Cluster for a draft genome"
fi

if [ $gene_prediction == 1 ] # With gene predictions
then
  echo "Running GI-Cluster with gene predictions"
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
if [ $gene_prediction == 0 ] # Not dependant on gene predictions
then
  if [ ! -f $output_dir/$seg_prog/feature/$organism.feature.multi.percentage ]
  then
    sh $prog_dir/Segment_Feature.sh -s $prog_dir -o $output_dir -n $organism -m $seg_prog -p $pred_prog -b $mode -t 0
    wait
  fi
else
  if [ ! -f $output_dir/$pred_prog/$seg_prog/$organism.feature.multi.percentage ]
  then
    sh $prog_dir/Segment_Feature.sh -s $prog_dir -o $output_dir -n $organism -m $seg_prog -p $pred_prog -b $mode
    wait
  fi
fi

###################### Run Rscript for consensus clustering ##############################
echo "##########################################"
echo "Running consensus clustering on the feature matrix"
# TODO: add parameters to input
# remember to put the separator at the end of output_dir
pFeature=1
method=average
rep=1
if [ $gene_prediction == 0 ] # Not dependant on gene predictions
then
  echo "Using features related to GC and k-mer"
  feature=gc_kmer
  postprocess=false
  if [ ! -d $output_dir/unannotated/$seg_prog/$feature/$method/$pFeature/$rep/ ]
  then
    mkdir -p $output_dir/unannotated/$seg_prog/$feature/$method/$pFeature/$rep/
  fi
  if [ ! -f $output_dir/unannotated/$seg_prog/$feature/$method/$pFeature/$rep/"$organism"_GI ]
  then
    ffile=$output_dir/$seg_prog/feature/$organism.feature.multi.percentage
    # Not use preprocessing
    nohup Rscript $prog_dir/GI_Clustering.R -f $ffile -o $output_dir/unannotated/$seg_prog/$feature/$method/$pFeature/$rep/ -a "$organism" -l $prog_dir/clustering -k 2 -K 3 -r $rep -e $pFeature -s $seg_prog -v true -C hclust -P method=$method -m true -d $feature -S $postprocess -E false -g 0 2>&1 > $output_dir/unannotated/$seg_prog/$feature/$method/$pFeature/$rep/rscript_std
  fi
else  # Dependant on gene predictions
  echo "Using features related to sequence compostion and gene functions"
  feature=comp_content
  postprocess=true
  if [ ! -d $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/ ]
  then
    mkdir -p $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/
  fi
  if [ ! -f $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/"$organism"_GI ]
  then
    ffile=$output_dir/$pred_prog/$seg_prog/$organism.feature.multi.percentage
    # Use preprocessing by default
    nohup Rscript $prog_dir/GI_Clustering.R -f $ffile -o $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/ -a "$organism" -l $prog_dir/clustering -k 2 -K 3 -r $rep -e $pFeature -s $seg_prog -v true -C hclust -P method=$method -m true -d $feature -S $postprocess 2>&1 > $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/rscript_std
  fi
fi

###################### Boundary refinement ##############################
echo "##########################################"
echo "Refining the boundary of predicted GIs"
# TODO: add parameters to input
dist=1000
# No need to distinguish the mode.
if [ $gene_prediction == 1 ] # Dependant on gene predictions
then
  orig_gifile=$output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/"$organism"_GI
  refined_gifile=$output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/"$organism"_refined_GI
else
  orig_gifile=$output_dir/unannotated/$seg_prog/$feature/$method/$pFeature/$rep/"$organism"_GI
  refined_gifile=$output_dir/unannotated/$seg_prog/$feature/$method/$pFeature/$rep/"$organism"_refined_GI
fi
python $prog_dir/postprocess/refine_boundary.py -i $orig_gifile -r $output_dir/$seg_prog/feature/$organism.boundary.repeat -t $output_dir/$seg_prog/feature/$organism.boundary.trna -o $refined_gifile -c $dist


# post process GI candidates
if [ $mode == 0 ] # For complete genomes
then
  if [ $gene_prediction == 1 ] # With gene predictions
  then
    python $prog_dir/postprocess/postprocess_segments.py -l $gap -i $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/"$organism"_refined_GI -a -g $output_dir/$pred_prog/genome/$organism.glist
  else
    python $prog_dir/postprocess/postprocess_segments.py -l $gap -i $output_dir/unannotated/$seg_prog/$feature/$method/$pFeature/$rep/"$organism"_refined_GI
  fi
else
  if [ $gene_prediction == 1 ] # With gene predictions
  then
    final_dir=$output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep
    python $prog_dir/postprocess/postprocess_segments.py -c -l $gap -i $final_dir/"$organism"_refined_GI -g $output_dir/$pred_prog/genome/$organism.glist -a  -m $output_dir/$organism.fna -d $output_dir/$pred_prog/genome/$organism.gene_id
    python $prog_dir/postprocess/convert_contig_GIs.py -i $final_dir/merged_"$organism"_refined_GI -o $final_dir/reformated_merged_"$organism"_refined_GI -g $output_dir/$organism.fna
  else
    final_dir=$output_dir/unannotated/$seg_prog/$feature/$method/$pFeature/$rep
    python $prog_dir/postprocess/postprocess_segments.py -c -l $gap -i $final_dir/"$organism"_refined_GI -g $output_dir/$pred_prog/genome/$organism.glist
    python $prog_dir/postprocess/convert_contig_GIs.py -i $final_dir/merged_"$organism"_refined_GI -o $final_dir/reformated_merged_"$organism"_refined_GI -g $output_dir/$organism.fna
  fi
fi


##################### Get the features of predicted GIs #####################
echo "##########################################"
echo "Getting the features of predicted GIs"
gifile1=$output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/merged_"$organism"_refined_GI
gifile2=$output_dir/unannotated/$seg_prog/$feature/$method/$pFeature/$rep/merged_"$organism"_refined_GI
seg_feature1=$output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/$organism.gi.feature
seg_feature1_percentage=$output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/$organism.gi.feature.percentage
seg_feature2=$output_dir/unannotated/$seg_prog/$feature/$method/$pFeature/$rep/$organism.gi.feature.percentage
if [ $gene_prediction == 1 ] # With gene predictions
then
  segs=$gifile1
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
  segs=$gifile2
  python $prog_dir/feature/countFeature.py -i $segs -o $seg_feature2 -r $output_dir/$seg_prog/feature/$organism.boundary.repeat -t $output_dir/$seg_prog/feature/$organism.boundary.trna -k  $output_dir/$seg_prog/feature/$organism.kmer.covariance -g $output_dir/$seg_prog/feature/$organism.gc
fi

##################### Visualization  ##############################
if [ $mode == 0 ] # Only support complete genomes and gene predictions
  then
    if [ $gene_prediction == 1 ]
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

      python $prog_dir/visualization/prepare_intervals.py -g $output_dir/$organism.fna -i $gifile1 -o $output_dir/$seg_prog/feature/visualization/data/$organism.gi -c $output_dir/$seg_prog/feature/visualization/data/$organism.chr -f $output_dir/$seg_prog/feature/visualization/data/$organism.highlight
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
    sh $prog_dir/GI_Comparison.sh -s $prog_dir -o $output_dir -n $organism -m $seg_prog -p $pred_prog
  fi

fi
