########## Predict genomic islands (GIs) via integrating multiple GI-related features from a newly sequenced (microbial) genome ##########

########## Input ##########
# The genome sequence of a newly sequenced bacterium: fasta file format
# OPTIONAL INPUT: $organism.fna, $organism.ffn, $organism.faa

########## prebuilt database ##########
# PFAM (ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz)
# RFAM (ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.tar.gz)
# PHAST (http://phast.wishartlab.com/phage_finder/DB/prophage_virus.db)
# VFDB (http://www.mgc.ac.cn/VFs/download.htm)
# CARD (https://card.mcmaster.ca/download)
# COG (https://www.ncbi.nlm.nih.gov/COG/)

########## external tools ##########
# prodigal  -- gene prediction
# sigihmm -- codon bias
# blast -- database search
# hmmmer -- database search
# circos -- visualization
# tRNAscan-SE
# Infernal cmscan  -- ncRNA

########## usage ##########
# ./GI-Cluster.sh $prog_dir $output_dir $organism
# e.g. ./GI-Cluster.sh $output_dir/GIFilter /home/ice/vmshare/research/data/species/cft73 NC_004431

# Contact: bingxin@comp.nus.edu.sg

software=$(basename $0)

function usage() {
  echo -e "GI-Cluster: detecting genomic islands in newly sequenced microbial genomes by consensus clustering on multiple features"
  echo "Version 0.1
Usage: $software [options] -s [the directory containing all the scripts] -o [the output directory]
-n [the name of the organism (NCBI Accession, e.g. NC_003198)] -m [programs for genome segmation (e.g. mjsd, gcprofile, gisvm, alienhunter)] -p [programs for gene prediction (e.g. prodigal, ncbi)]

OPTIONS	Default	DESCIPTION
-g	2500	: gaps allowed when merging adjacent segments.

-t	1e-5	: e-value used during identification of mobgenes, i.e., hmmsearch against pfam database of mobgenes.
-e	1e-5	: e-value used during identification of phage-related genes, i.e., blastp against PHAST.

-r	1e-5	: e-value used during identification of virulence factors, i.e., blastp against VFDB.
-a	1e-5	: e-value used during identification of antibiotic resistance genes, i.e., blastp against CARD.

-d	4	: number of threads used by blast.
-u	16	: number of cpus used by cmsearch.

-f 	1	: show the visualization of features related to genomic islands
-h 	----	: print this help
  "
  exit -1
}


phage_evalue=1e-5
virdb_evalue=1e-5
arg_evalue=1e-5
num_threads=4
show_figure=1
num_cpus=16
len=2500

while getopts "s:o:n:m:p:g:c:e:r:a:d:f:u:h" OPT; do
  case $OPT in
    s) prog_dir=$OPTARG || exit 1;;
    o) output_dir=$OPTARG || exit 1;;
    n) organism=$OPTARG || exit 1;;
    m) seg_prog=$OPTARG || exit 1;;
    p) pred_prog=$OPTARG || exit 1;;

    g) len=$OPTARG || exit 1;;

    #t) mobgene_evalue=$OPTARG || exit 1;;
    e) phage_evalue=$OPTARG || exit 1;;

    r) virdb_evalue=$OPTARG || exit 1;;
    a) arg_evalue=$OPTARG || exit 1;;

    d) num_threads=$OPTARG || exit 1;;
    f) show_figure=$OPTARG || exit 1;;
    u) num_cpus=$OPTARG || exit 1;;
    h) usage && exit;;
  esac
done


#echo $prog_dir
#echo $output_dir
#echo $organism

############## predict genes from raw genome sequence ##################
# 3 types of required output files: fna, ffn, faa
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
    $prog_dir/bin/prodigalrunner $output_dir/$organism.fna
    mv $output_dir/"$organism"_prodigal.orf.fna $output_dir/$pred_prog/genome/"$organism".ffn
    mv $output_dir/"$organism"_prodigal.orf.fsa $output_dir/$pred_prog/genome/"$organism".faa
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
  less $output_dir/$pred_prog/genome/"$organism".faa | grep '^>'  | cut -d' ' -f1 > $output_dir/$pred_prog/genome/$organism.gene_id
  # for each gene, a record in a row
  if [ "$pred_prog" == "prodigal" ]
  then
    python $prog_dir/parse_geneloc.py -l $output_dir/$pred_prog/genome/$organism.gene_id  -o $output_dir/$pred_prog/genome/$organism.glist
  else
    # ncbi annotation files
    less $output_dir/$pred_prog/genome/"$organism".ffn | grep '^>'  | cut -d':' -f2 | cut -d' ' -f1 > $output_dir/$pred_prog/genome/$organism.gene_locus
    python $prog_dir/parse_geneloc.py -n -l $output_dir/$pred_prog/genome/$organism.gene_locus  -o $output_dir/$pred_prog/genome/$organism.glist
  fi
fi

############## compute known features related to genomic islands #########################
$prog_dir/GI_Feature.sh -s $prog_dir -o $output_dir -n $organism -m $seg_prog -p $pred_prog
wait
############################## genome segmentation ####################################
echo "##########################################"
echo "Get genome segements"
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
    python $prog_dir/generate_segs.py -i $output_dir/$organism.fna -o $output_dir/$seg_prog/$organism."$seg_prog"
  fi
  if [ "$seg_prog" == "$p1" ]
  then
    # run mjsd, -s 0.99
    $prog_dir/bin/mjsd/so_jensen -f $output_dir/$organism.fna -s 0.1 -o 2 > $output_dir/$seg_prog/$organism.$seg_prog.orig
    python $prog_dir/parse_MJSD_segs.py -f $output_dir/$organism.fna -m 5000 -s $output_dir/$seg_prog/$organism.$seg_prog.orig   > $output_dir/$seg_prog/$organism."$seg_prog"
  fi

  if [ "$seg_prog" == "$p2" ]
  then
    # -t 25
    python $prog_dir/bin/gcprofile/parse_gcprofile_segs.py -t 10 -m 5000 -f $output_dir/$organism.fna -p $prog_dir/bin/gcprofile/linux/ -o $output_dir/$seg_prog/$organism."$seg_prog"
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

cp $output_dir/$seg_prog/"$organism"."$seg_prog" $output_dir/$seg_prog/$organism.segment

#################################### summarize features of each genomic region ##############################
# the features are different for different gene predictions
echo "##########################################"
echo "Get features of each genome segement"
if [ ! -d $output_dir/$pred_prog/$seg_prog ]
then
  mkdir -p $output_dir/$pred_prog/$seg_prog
fi

if [ ! -f $output_dir/$pred_prog/$seg_prog/$organism.boundary.trna ]
then
python $prog_dir/find_trnas.py -i $output_dir/$seg_prog/$organism.segment -t $output_dir/boundary/$organism.pred_trna -o $output_dir/$pred_prog/$seg_prog/$organism.boundary.trna
fi

if [ ! -f $output_dir/$pred_prog/$seg_prog/$organism.boundary.repeat ]
then
# The output is based on the genomic segments
python $prog_dir/find_repeats.py -i $output_dir/$seg_prog/$organism.segment -r $output_dir/boundary/$organism.repseek -o $output_dir/$pred_prog/$seg_prog/$organism.boundary.repeat
fi

if [ ! -f $output_dir/$seg_prog/$organism.kmer.covariance ]
then
python $prog_dir/analyze_kmer.py -r -i $output_dir/$organism.fna -s $output_dir/$seg_prog/$organism.segment -o $output_dir/$seg_prog/$organism.kmer
fi

python $prog_dir/mergeFeature.py -g $output_dir/$pred_prog/feature/$organism.feature.multi -o $output_dir/$pred_prog/$seg_prog/$organism.seg_feature.multi  -s $output_dir/$seg_prog/$organism.segment -r $output_dir/$pred_prog/feature/$organism.cmscan.tbl


# Form the feature matrix for clustering, for gc and kmer use the value for each segment directly
python $prog_dir/countFeature.py -i $output_dir/$pred_prog/$seg_prog/$organism.seg_feature.multi -o $output_dir/$pred_prog/$seg_prog/$organism.feature.multi.percentage -r $output_dir/$pred_prog/$seg_prog/$organism.boundary.repeat -t $output_dir/$pred_prog/$seg_prog/$organism.boundary.trna -k $output_dir/$seg_prog/$organism.kmer.covariance

# Get the labels for each segment when there are references
# have to run this after changing references
# python $prog_dir/IntervalIntersection.py --interval1 $output_dir/$seg_prog/$organism.segment --interval2 $output_dir/$organism.refgi -o $output_dir/$seg_prog/"$seg_prog"_refgi -f 0
#
# python $prog_dir/assign_cluster.py --interval1 $output_dir/$pred_prog/$seg_prog/$organism.feature.multi.percentage --interval2 $output_dir/$seg_prog/"$seg_prog"_refgi_interval1 -o $output_dir/$pred_prog/$seg_prog/$organism.feature.multi.percentage.labeled
# python $prog_dir/assign_cluster.py -c --interval1 $output_dir/$pred_prog/$seg_prog/$organism.feature.multi.percentage.labeled --interval2 $output_dir/$seg_prog/"$seg_prog"_pos_overlap_interval1 --interval3 $output_dir/$seg_prog/"$seg_prog"_neg_overlap_interval1 -o $output_dir/$pred_prog/$seg_prog/$organism.feature.multi.percentage.labeled2


###################### run Rscript for consensus clustering ##############################
echo "##########################################"
echo "Run consensus clustering on the feature matrix"
# remember to put the separator at the end of output_dir
pFeature=1
feature=comp_content
method=average
rep=1
postprocess=true

if [ ! -f $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/"$organism"_GI ]
then
  mkdir -p $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/
  ffile=$output_dir/$pred_prog/$seg_prog/$organism.feature.multi.percentage
  nohup Rscript $prog_dir/GI_Clustering.R -f $ffile -o $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/ -a "$organism" -l $prog_dir/clustering -k 2 -K 3 -r $rep -e $pFeature -s $seg_prog -v true -C hclust -P method=$method -m true -d $feature -S $postprocess 2>&1 > $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/rscript_std
fi

###################### boundary refinement ##############################
echo "##########################################"
echo "Refine the boundary of predicted GIs"
dist=1000
python $prog_dir/refine_boundary.py -i $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/"$organism"_GI -r $output_dir/$pred_prog/$seg_prog/$organism.boundary.repeat -t $output_dir/$pred_prog/$seg_prog/$organism.boundary.trna -o  $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/"$organism"_refined_GI -c $dist

# post process GI candidates
# python $prog_dir/postprocess_segments.py  -l $gap -a -i $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/"$organism"_GI -g $output_dir/$pred_prog/genome/$organism.glist
gap=5000
python $prog_dir/postprocess_segments.py  -l $gap -a -i $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/"$organism"_refined_GI -g $output_dir/$pred_prog/genome/$organism.glist

##################### visualization  ##############################
gifile=$output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/merged_"$organism"_refined_GI
# get the features of predicted GIs
python $prog_dir/mergeFeature.py -g $output_dir/$pred_prog/feature/$organism.feature.multi -s $gifile -r $output_dir/$pred_prog/feature/$organism.cmscan.tbl -o $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/$organism.gi.feature

# get the features of FPs
python $prog_dir/mergeFeature.py -g $output_dir/$pred_prog/feature/$organism.feature.multi -o $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/eval_FPs_GI_feature  -s $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep/eval_FPs_GI -r $output_dir/$pred_prog/feature/$organism.cmscan.tbl


echo "##########################################"
echo "Create figures to visualize the predicted genomic island and related features"
# sh $prog_dir/GI_Visualization.sh -s $prog_dir -o $output_dir -n $organism -m $seg_prog -p $pred_prog
show=1
if [ "$show_figure" -eq "$show" ]
then
  if [ ! -d $output_dir/$pred_prog/$seg_prog/visualization/data ]
  then
    mkdir -p $output_dir/$pred_prog/$seg_prog/visualization/data
    mkdir -p $output_dir/$pred_prog/$seg_prog/visualization/etc
    mkdir -p $output_dir/$pred_prog/$seg_prog/visualization/img
  fi
  # copy template files
  cp $prog_dir/visualization/etc/* $output_dir/$pred_prog/$seg_prog/visualization/etc

  python $prog_dir/visualization/prepareGIForCircos.py -g $output_dir/$organism.fna -i $gifile -o $output_dir/$pred_prog/$seg_prog/visualization/data/$organism.gi -c $output_dir/$pred_prog/$seg_prog/visualization/data/$organism.chr -f $output_dir/$pred_prog/$seg_prog/visualization/data/$organism.highlight
  # use features of all the input segments
  python $prog_dir/visualization/prepareGIFeatureForCircos.py -a -i $output_dir/$pred_prog/$seg_prog/$organism.feature.multi.percentage -o $output_dir/$pred_prog/$seg_prog/visualization/data

  # prefix="intervals.feature.percentage"
  prefix="feature.multi.percentage"
  python $prog_dir/visualization/prepareConfigForCircos.py -i $output_dir/$pred_prog/$seg_prog/visualization/etc/circos.gifeature.conf.template -n $organism -o $output_dir/$pred_prog/$seg_prog/visualization/etc/circos.gifeature.conf -p $prefix -f "$organism"_gifeature

  # run circos at etc folder since the output directory is relative
  cd $output_dir/$pred_prog/$seg_prog/visualization/etc
  circos -conf $output_dir/$pred_prog/$seg_prog/visualization/etc/circos.gifeature.conf
  cd $output_dir

  # compare with reference GIs
fi
