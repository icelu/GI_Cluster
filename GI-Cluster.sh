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

############## predict genes from raw genome sequence ##################
# 3 types of required output files: fna (DNA sequence), ffn (gene sequence), faa (protein sequence)
if [ $mode != 1 ]
then
if [ ! -d $output_dir/$pred_prog/genome/ ]
then
  mkdir -p $output_dir/$pred_prog/genome/
fi

if [ "$pred_prog" == "prodigal" ]
then
  echo "##########################################"
  echo "Predict genes with prodigal"
  if [ ! -f $output_dir/$pred_prog/genome/"$organism".ffn ]
  then
    # $prog_dir/bin/prodigalrunner $output_dir/$organism.fna
    # mv $output_dir/"$organism"_prodigal.orf.fna $output_dir/$pred_prog/genome/"$organism".ffn
    # mv $output_dir/"$organism"_prodigal.orf.fsa $output_dir/$pred_prog/genome/"$organism".faa
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
  echo "Gene predictions from NCBI"
  if [ ! -f $output_dir/$pred_prog/genome/"$organism".ffn ]
  then
    cp $output_dir/"$organism".ffn $output_dir/$pred_prog/genome/"$organism".ffn
    cp $output_dir/"$organism".faa $output_dir/$pred_prog/genome/"$organism".faa
  fi
fi

# get gene/orf list: name, loc, strand
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

############## compute known features related to genomic islands #########################
if [ "$pred_prog" != "none" ]
then
sh $prog_dir/GI_Feature.sh -s $prog_dir -o $output_dir -n $organism -m $seg_prog -p $pred_prog -b $mode -e $phage_evalue -r $virdb_evalue -a $arg_evalue -d $num_threads -u $num_cpus
wait
fi
fi # for complete genomes


############################## genome segmentation ####################################
echo "##########################################"
echo "Get genome segements"
if [ ! -f $output_dir/$seg_prog/$organism.segment ]
then
sh $prog_dir/GI_Segmentation.sh -s $prog_dir -o $output_dir -n $organism -m $seg_prog -p $pred_prog -b $mode
wait
fi

#################################### summarize features of each genomic region ##############################
# the features are different for different gene predictions
echo "##########################################"
echo "Get features of each genome segement"
if [ ! -d $output_dir/$seg_prog/feature ]
then
  mkdir -p $output_dir/$seg_prog/feature
fi
if [ ! -d $output_dir/$pred_prog/$seg_prog ]
then
if [ "$pred_prog" != "none" ]
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
  else
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

if [ $mode == 1 ] # Only for contigs without gene predictions, not dependant on gene predictions
then
echo "##########################################"
echo "Analyzing GC Content for each segment"
if [ ! -f $output_dir/$seg_prog/feature/$organism.gc ]
then
  python $prog_dir/feature/analyze_GC.py -c -i $output_dir/$organism.fna -s $output_dir/$seg_prog/$organism.segment -o $output_dir/$seg_prog/feature/$organism.gc
fi
fi

echo "##########################################"
echo "Analyzing k-mer frequency for each segment"
if [ ! -f $output_dir/$seg_prog/feature/$organism.kmer.covariance ]
then
  if [ $mode == 0 ]
  then
    python $prog_dir/feature/analyze_kmer.py -r -i $output_dir/$organism.fna -s $output_dir/$seg_prog/$organism.segment -o $output_dir/$seg_prog/feature/$organism.kmer
  else # For contigs, not dependant on gene predictions
    python $prog_dir/feature/analyze_kmer.py -c -i $output_dir/$organism.fna -s $output_dir/$seg_prog/$organism.segment -o $output_dir/$seg_prog/feature/$organism.kmer
  fi
fi


echo "##########################################"
echo "Merging features for each segment"
if [ $mode == 2 ] # For contigs with gene predictions
then
python $prog_dir/feature/mergeFeature.py -g $output_dir/$pred_prog/feature/$organism.feature.multi -o $output_dir/$pred_prog/$seg_prog/$organism.seg_feature.multi  -s $output_dir/$seg_prog/$organism.segment -r $output_dir/$pred_prog/feature/$organism.cmscan.tbl -c -m $output_dir/$organism.fna -d $output_dir/$pred_prog/genome/$organism.gene_id
fi

if [ $mode == 1 ] # For contigs without gene predictions
then
python $prog_dir/feature/mergeFeature.py -g $output_dir/$pred_prog/feature/$organism.feature.multi -o $output_dir/$seg_prog/feature/$organism.seg_feature.multi  -s $output_dir/$seg_prog/$organism.segment -r $output_dir/$pred_prog/feature/$organism.cmscan.tbl
fi

if [ $mode == 0 ] # Dependant on gene predictions
then
python $prog_dir/feature/mergeFeature.py -g $output_dir/$pred_prog/feature/$organism.feature.multi -o $output_dir/$pred_prog/$seg_prog/$organism.seg_feature.multi  -s $output_dir/$seg_prog/$organism.segment -r $output_dir/$pred_prog/feature/$organism.cmscan.tbl
fi


# Form the feature matrix for clustering, for gc and kmer use the value for each segment directly
python $prog_dir/feature/countFeature.py -i $output_dir/$seg_prog/feature/$organism.seg_feature.multi -o $output_dir/$seg_prog/feature/$organism.feature.multi.percentage -r $output_dir/$seg_prog/feature/$organism.boundary.repeat -t $output_dir/$seg_prog/feature/$organism.boundary.trna -k $output_dir/$seg_prog/$organism.kmer.covariance

# Get the labels for each segment when there are references
# have to run this after changing references
# python $prog_dir/IntervalIntersection.py --interval1 $output_dir/$seg_prog/$organism.segment --interval2 $output_dir/$organism.refgi -o $output_dir/$seg_prog/"$seg_prog"_refgi -f 0
#
# python $prog_dir/assign_cluster.py --interval1 $output_dir/$seg_prog/feature/$organism.feature.multi.percentage --interval2 $output_dir/$seg_prog/"$seg_prog"_refgi_interval1 -o $output_dir/$seg_prog/feature/$organism.feature.multi.percentage.labeled
# python $prog_dir/assign_cluster.py -c --interval1 $output_dir/$seg_prog/feature/$organism.feature.multi.percentage.labeled --interval2 $output_dir/$seg_prog/"$seg_prog"_pos_overlap_interval1 --interval3 $output_dir/$seg_prog/"$seg_prog"_neg_overlap_interval1 -o $output_dir/$seg_prog/feature/$organism.feature.multi.percentage.labeled2


###################### run Rscript for consensus clustering ##############################
echo "##########################################"
echo "Run consensus clustering on the feature matrix"
# remember to put the separator at the end of output_dir
pFeature=1
feature=comp_content
method=average
rep=1
postprocess=true

if [ ! -f $output_dir/$seg_prog/feature/$feature/$method/$pFeature/$rep/"$organism"_GI ]
then
  mkdir -p $output_dir/$seg_prog/feature/$feature/$method/$pFeature/$rep/
  ffile=$output_dir/$seg_prog/feature/$organism.feature.multi.percentage
  nohup Rscript $prog_dir/GI_Clustering.R -f $ffile -o $output_dir/$seg_prog/feature/$feature/$method/$pFeature/$rep/ -a "$organism" -l $prog_dir/clustering -k 2 -K 3 -r $rep -e $pFeature -s $seg_prog -v true -C hclust -P method=$method -m true -d $feature -S $postprocess 2>&1 > $output_dir/$seg_prog/feature/$feature/$method/$pFeature/$rep/rscript_std
fi

###################### boundary refinement ##############################
echo "##########################################"
echo "Refine the boundary of predicted GIs"
dist=1000
python $prog_dir/boundary/refine_boundary.py -i $output_dir/$seg_prog/feature/$feature/$method/$pFeature/$rep/"$organism"_GI -r $output_dir/$seg_prog/feature/$organism.boundary.repeat -t $output_dir/$seg_prog/feature/$organism.boundary.trna -o  $output_dir/$seg_prog/feature/$feature/$method/$pFeature/$rep/"$organism"_refined_GI -c $dist

gifile=$output_dir/$seg_prog/feature/$feature/$method/$pFeature/$rep/merged_"$organism"_refined_GI

# post process GI candidates
python $prog_dir/boundary/postprocess_segments.py  -l $gap -a -i $gifile -g $output_dir/$pred_prog/genome/$organism.glist

##################### visualization  ##############################
# get the features of predicted GIs
python $prog_dir/feature/mergeFeature.py -g $output_dir/$pred_prog/feature/$organism.feature.multi -s $gifile -r $output_dir/$pred_prog/feature/$organism.cmscan.tbl -o $output_dir/$seg_prog/feature/$feature/$method/$pFeature/$rep/$organism.gi.feature

# # get the features of FPs
# python $prog_dir/mergeFeature.py -g $output_dir/$pred_prog/feature/$organism.feature.multi -o $output_dir/$seg_prog/feature/$feature/$method/$pFeature/$rep/eval_FPs_GI_feature  -s $output_dir/$seg_prog/feature/$feature/$method/$pFeature/$rep/eval_FPs_GI -r $output_dir/$pred_prog/feature/$organism.cmscan.tbl

echo "##########################################"
echo "Create figures to visualize the predicted genomic island and related features"
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
  # copy template files
  cp $prog_dir/visualization/etc/* $output_dir/$seg_prog/feature/visualization/etc

  python $prog_dir/visualization/prepareGIForCircos.py -g $output_dir/$organism.fna -i $gifile -o $output_dir/$seg_prog/feature/visualization/data/$organism.gi -c $output_dir/$seg_prog/feature/visualization/data/$organism.chr -f $output_dir/$seg_prog/feature/visualization/data/$organism.highlight
  # use features of all the input segments
  python $prog_dir/visualization/prepareGIFeatureForCircos.py -a -i $output_dir/$seg_prog/feature/$organism.feature.multi.percentage -o $output_dir/$seg_prog/feature/visualization/data

  prefix="feature.multi.percentage"
  python $prog_dir/visualization/prepareConfigForCircos.py -i $output_dir/$seg_prog/feature/visualization/etc/circos.gifeature.conf.template -n $organism -o $output_dir/$seg_prog/feature/visualization/etc/circos.gifeature.conf -p $prefix -f "$organism"_gifeature

  # run circos at etc folder since the output directory is relative
  cd $output_dir/$seg_prog/feature/visualization/etc
  circos -conf $output_dir/$seg_prog/feature/visualization/etc/circos.gifeature.conf
  cd $output_dir
fi
