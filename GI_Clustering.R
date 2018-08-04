#!/usr/bin/env Rscript

##=========================================================
# Script for running Consensus Clustering on genomic regions with extracted GI-related features
#
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
#

##=============================================================

# Utility functions
if(!require(R.utils)){
    install.packages("R.utils")
    library(R.utils)
}

#   for parsing parameters
if(!require(optparse)){
    install.packages("optparse")
    library(optparse)
}
# for existsFunction
if(!require(methods)){
    install.packages("methods")
    library(methods)
}
# for hclust
if(!require(stats)){
    install.packages("stats")
    library(stats)
}

# for 1D kmeans
if(!require(Ckmeans.1d.dp)){
    install.packages("Ckmeans.1d.dp")
    library(Ckmeans.1d.dp)
}

# for density estimation
if(!require(mclust)){
    install.packages("mclust")
    library(mclust)
}

# for cluster statistics
if(!require(fpc)){
    install.packages("fpc")
    library(fpc)
}

##=============================================
# Read the arguments
##=============================================
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="input dataset file name"),
  make_option(c("-o", "--outdir"), type="character", default="out.txt",
              help="output file name [default= %default]"),
  make_option(c("-l", "--libdir"), type="character", default="/home/b/bingxin/GI-Cluster/clustering",
              help="output file name [default= %default]"),
  make_option(c("-s", "--seg_prog"), type="character", default="window",
              help="source for input segments [default= %default]"),
  make_option(c("-k", "--clnum"), type="integer", default="2",
              help="maximum cluster number to evaluate [default= %default]"),
  make_option(c("-K", "--internalK"), type="integer", default="3",
              help="number of clusters for clustering on each feature[default= %default]"),
  make_option(c("-r", "--reps"), type="integer", default="50",
              help="number of subsamples [default= %default]"),
  make_option(c("-t", "--pItem"), type="double", default="1",
              help="proportion of items to sample [default= %default]"),
  make_option(c("-e", "--pFeature"), type="double", default="0.6",
              help="proportion of features to sample [default= %default]"),
  make_option(c("-d", "--feature"), type="character", default="comp_content",
              help="group of feature to use [default= %default]"),
  make_option(c("-c", "--clusterAlg"), type="character", default="Ckmeans.1d.dp,densityMclust",
              help="cluster algorithm used [default= %default]"),
  make_option(c("-p", "--alParams"), type="character", default="",
              help="parameters for the cluster algorithms used [default= %default]"),
  make_option(c("-C", "--mergeAlg"), type="character", default="hclust",
              help="cluster algorithm used on consensus matrix [default= %default]"),
  make_option(c("-P", "--mergeParams"), type="character", default="method=average",
              help="cluster algorithm used [default= %default]"),
  make_option(c("-a", "--organism"), type="character", default="NC_010170",
              help="accession number of the genome [default= %default]"),
  make_option(c("-E", "--pre_process"), type="character", default="true",
              help="output more information for reference [default= %default]"),
  make_option(c("-S", "--post_process"), type="character", default="false",
              help="output more information for reference [default= %default]"),
  make_option(c("-v", "--verbose"), type="character", default="false",
              help="output more information for reference [default= %default]"),
  make_option(c("-m", "--heatmap"), type="character", default="false",
              help="output consensus matrix as a heatmap [default= %default]"),
  make_option(c("-g", "--gene_prediction"), type="integer", default="1",
              help="availablity of gene predictions: 1: with gene predictions, 0: without gene predictions")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

file = opt$file
outdir = opt$outdir
libdir = opt$libdir
gene_prediction = opt$gene_prediction
clnum = opt$clnum
internalK = opt$internalK
seg_prog = opt$seg_prog

reps = opt$reps
pItem = opt$pItem
pFeature = opt$pFeature
feature = opt$feature

clusterAlg = as.list(strsplit(opt$clusterAlg, ",")[[1]])
if(opt$alParams!=""){
  aParams = as.list(strsplit(opt$alParams, ",")[[1]])
  alParams=list()
  for(i in 1:length(alParams)){
    param = strsplit(alParams[[i]],"=")[[1]]
    plist=list()
    key = param[1]
    plist[[key]]=param[2]
    alParams[[i]]=plist
  }
}else{
  alParams=list()
}

mergeAlg = as.list(strsplit(opt$mergeAlg, ",")[[1]])

# mergeP = "method='average',method='single'"
mParams = as.list(strsplit(opt$mergeParams, ",")[[1]])
mergeParams=list()
for(i in 1:length(mParams)){
  param = strsplit(mParams[[i]],"=")[[1]]
  plist=list()
  key = param[1]
  plist[[key]]=param[2]
  mergeParams[[i]]=plist
}

organism = opt$organism

pre_process = opt$pre_process
post_process = opt$post_process
verbose = opt$verbose
heatmap = opt$heatmap

if (verbose=="true"){
  cat("Parameters used here:\n")
  cat("clnum:", clnum, "\n")
  cat("seg_prog:", seg_prog, "\n")
  cat("reps:", reps, "\n")
  cat("pItem:", pItem, "\n")
  cat("pFeature:", pFeature, "\n")
  cat("clusterAlg", str(clusterAlg), "\n\n")
  cat("alParams:", str(alParams), "\n")
  cat("mergeAlg", str(mergeAlg), "\n\n")
  cat("mergeParams", str(mergeParams), "\n\n")
}

R.utils::sourceDirectory(libdir, modifiedOnly=FALSE)

# The format of the feature table:
# For genome with gene annotations: (31 columns)
# (ID, start, end, size, has_trna, has_repeat, cdst, rnat, gc, gc1, gc3, cub, aab, chi, cai, cbi, fop, kmer2, kmer3, kmer4, kmer5, kmer6, kmer7, kmer8, mbp, phagep, vfp, arp, ngp, cdsp, inter_gene_dist)
# For genome without gene annotations: (14 columns)
# (ID, start, end, size, has_trna, has_repeat, gc, kmer2, kmer3, kmer4, kmer5, kmer6, kmer7, kmer8)
ft0 <- read.table(file)

# Remove regions with 0 genes or with too many ncRNAs
{
if(pre_process=="true"){
    # V7 -- Number of genes, V8 -- Number of ncRNAs
    fte = ft0[ft0$V7==0|ft0$V7!=0&(ft0$V8/(ft0$V7+ft0$V8))>0.5,]
    ft = ft0[ft0$V7!=0&(ft0$V8/(ft0$V7+ft0$V8))<=0.5,]
}
else{
  ft = ft0;
}
}

{
# gc, gc1, gc3, cub, aab, chi, cai, cbi, fop, kmer2, kmer3, kmer4, kmer5, kmer6, kmer7, kmer8, mbp, phagep, vfp, arp, ngp, cdsp, inter_gene_dist
if(gene_prediction==1){
  features=ft[,9:24]

  # Transform features related to GC content into distance values
  # For each value, compute its chisquare value (x-mean)^2/mean
  start=9;
  for(i in 1:4){
    if(i==3) next;
    fv = ft[,start+i]*100
    mu = mean(fv)
    # cat("mu ", paste0("V",start+i), ": ", mu, "\n")
    features[, paste0("V",start+i)] = ((fv-mu)^2)/mu
    # cat("feature ", paste0("V",start+i), ": ", str(features[, paste0("V",start+i)]), "\n")
  }

  # Features related to gene content (phage, vf, ar, ng)
  for(i in 1:4){
    features[, paste0("V",25+i)] = ft[,25+i]
    # cat("feature ", paste0("V",25+i), ": ", str(features[, paste0("V",25+i)]), "\n")
  }

  {
  if (feature=="gc") {data=features[1:3]}
  else if (feature=="codon")
  {
    # data=features[4:9]
    data=features[4:7]

  }
  else if (feature=="kmer"){
    data=features[10:16]
  }
  else if (feature=="content")
  {
    data=features[17:20]
  }
  else if (feature=="gc_kmer")
  {
    data1=features[1:3]
    data2=features[10:16]
    data=cbind(data1, data2)
  }
  # Only use compositional features for clustering
  else if (feature=="comp")
  {
    data1=features[1:3]
    data2=features[4:7]
    data3=features[10:16]
    data=cbind(data1, data2, data3)
  }
  else if (feature=="comp_content")
  {
    data1=features[1:3]
    data2=features[4:7]
    data3=features[10:16]
    data4=features[17:20]
    # data4=features[17]
    data=cbind(data1, data2, data3, data4)
  }
  else
  {
    data1=features[1:3]
    data2=features[4:7]
    data3=features[10:16]
    data4=features[17:20]
    data=cbind(data1, data2, data3, data4)
  }

  }
} # With gene predictions
else{
# gc, kmer2, kmer3, kmer4, kmer5, kmer6, kmer7, kmer8
features=ft[,7:14]

fv = ft[,7]*100
mu = mean(fv)
# cat("mu ", paste0("V",start+i), ": ", mu, "\n")
features[, "V7"] = ((fv-mu)^2)/mu

if (feature=="gc") {data=features[1]}
else if (feature=="kmer"){
  data=features[2:8]
}
else  # gc_kmer
{
  data=features
}
}
}

cat("Number of selected features:", floor(pFeature*ncol(data)), "\n")
# Set 3 clusters as default to divide each feature into 3 groups: native, alien, ambiguous
outdir_sep=paste0(outdir, "sepcluster/")
ifelse(!dir.exists(outdir_sep), dir.create(outdir_sep), FALSE)
res_list=clustering.1D(data, clnum=internalK, pCol=pFeature, algorithms=clusterAlg,alparams=alParams, rep=reps, verbose=verbose, outdir=outdir_sep, ft=ft)

# Compute the average concensus matrix
avg_cm=matrix(0,nrow=dim(data)[1],ncol=dim(data)[1],dimnames=list(row.names(data),row.names(data)));
for(i in 1:length(res_list)){
  avg_cm = avg_cm+res_list[[i]]
}
avg_cm = avg_cm/length(res_list)
str(avg_cm)

thisPal <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928",
             "#bd18ea", #magenta
             "#2ef4ca", #aqua
             "#f4cced", #pink,
             "#f4cc03", #lightorange
             "#05188a", #navy,
             "#e5a25a", #light brown
             "#06f106", #bright green
             "#85848f", #med gray
             "#000000", #black
             "#076f25", #dark green
             "#93cd7f",#lime green
             "#4d0776", #dark purple
             "#ffffff" #white
)

myPal = function(n=10){
	#returns n colors
	seq = rev(seq(0,255,by=255/(n)))
	palRGB = cbind(seq,seq,255)
	rgb(palRGB,maxColorValue=255)
}

########################################### run clustering on consensus matrix ############################################
for(j in 1:length(mergeAlg)){
  algo = mergeAlg[[j]]
  message('clustering on consensus matrix by algorithm ',algo)
  current_params = c(mergeParams[[j]],list(diss=TRUE));
  algo_clmem <- paste(algo,'_clmem',sep='');
  current_algo <- get(algo_clmem);

  x = as.dist(1- avg_cm);

  cls_res <- current_algo(x,clnum,params=current_params);
  clmem <- cls_res$clmem;

  if( heatmap == "true")
  {
    figfile=paste(outdir, organism, "_Ccutoff_", clnum, ".pdf",  sep="")
    pdf(figfile, width=11, height=11*0.618);
    hc=cls_res$clust_obj;
    ct = cutree(hc,clnum);
    names(ct) = row.names(clmem);
    pc = avg_cm;
    pc=pc[hc$order,] # pc is matrix for plotting. It is similar to c but it is row-ordered and has names and extra row of zeros.

    # Set colors for heatmap
    colBreaks=10;
    tmyPal = myPal(colBreaks);
    newColors = thisPal[ct];
    colori=2;
    colorList = list(newColors,colori,unique(newColors));

    pc = rbind(pc,0)
    #former with tree:
    heatmap(pc, Colv=as.dendrogram(hc), Rowv=NA, symm=FALSE, scale='none', col=tmyPal, na.rm=TRUE,labRow=F,labCol=F,mar=c(5,5),main=paste("consensus matrix k=",clnum,sep="") , ColSideCol=colorList[[1]]);
    dev.off()
  }

  Clab2= as.vector(t(clmem))
  c1=which(grepl(1,Clab2))
  c2=which(grepl(2,Clab2))

{
  if(seg_prog == "window"){
    if(length(c1)>length(c2)){ # suppose less regions are GI for output from genome segmentation
      gi = c2
      ngi = c1
      # Assign gi to be 1, ngi to 0 for comparison
      Clab2 = replace(Clab2, Clab2==1, 0)
      Clab2 = replace(Clab2, Clab2==2, 1)
    } else
    {
      gi = c1
      ngi = c2
      # Assign gi to be 1, ngi to 0 for comparison
      Clab2 = replace(Clab2, Clab2==2, 0)
    }
  }
}

{
  if(seg_prog != "window"){
    if(length(c1)<length(c2)){ # suppose more regions are GIs for output from GI prediction tools
      gi = c2
      ngi = c1
      # Assign gi to be 1, ngi to 0 for comparison
      Clab2 = replace(Clab2, Clab2==1, 0)
      Clab2 = replace(Clab2, Clab2==2, 1)
    } else
    {
      gi = c1
      ngi = c2
      # Assign gi to be 1, ngi to 0 for comparison
      Clab2 = replace(Clab2, Clab2==2, 0)
    }
  }
}
  # If there are known labeled for each item
  # if("V32" %in% colnames(ft))
  # {
  #   message('comparison with segments labeled as GIs in L-data:')
  #   print(table(ft$V32,Clab2))
  # }
  #
  # if("V33" %in% colnames(ft))
  # {
  #   message('comparison with segments labeled as GIs in C-data:')
  #   print(table(ft$V33,Clab2))
  # }
  # Compute cluster stats
  message('cluster stats:')
  clust_stats <- cluster.stats(d = x, as.vector(t(Clab2)))
  # Corrected Rand index
  # print(clust_stats)
  cat("silwidth: ", clust_stats$avg.silwidth, "\n")
  cat("dunn: ", clust_stats$dunn, "\n")
  cat("dunn2: ", clust_stats$dunn2, "\n")

  # message('internal clustering validation:')
  # x=as.matrix(1-avg_cm)
  # part=as.integer(Clab2)
  # print(intCriteria(x, part,"all"))
  cat('#########################\n')
}


############################### postprocess to remove potential FPs and FNs ########################
gi0 = ft[gi,]
ngi0 = ft[ngi,]

{
if (post_process=="false")  # For mode 2, post_process should be "false"
{
  gi = gi0
  ngi = ngi0
}
# The rules used in postprocessing can be modified to adapt different requirements
else {
  if(seg_prog=="window"){
    # find GIs from ngi by other evidence
    # regions with a tRNA at the boundary, mobility evidence
    gi1 = subset(ngi0,(V5==1)&(V25>0))
    ngi1 = subset(ngi0,! V1 %in% as.vector(gi1$V1))
    # Check gene content information
    # Regions with mostly phage-related genes or novel genes (prophage clusters)
    gi2 = subset(ngi1, ((V26>0&V26+V29>=0.8)&V7>=5))
    # gi2 = subset(ngi1, ((V26>0&V26+V29>=0.8))|(V27==1)|(V28==1)|(V29==1)|(V25>0&V26+V29>=0.6))  # Relaxed rule
    # Exclude FP regions: no evidence of mobility and other genes enriched in islands
    cutoff_density=quantile(ft$V30,c(0.5))[1]
    cutoff_distance=quantile(ft$V31,c(0.5))[1]
    cutoff = 0.1
    cutoff_phage=max(quantile(ft$V26,c(0.5))[1], cutoff)
    cutoff_vf=max(quantile(ft$V27,c(0.5))[1], cutoff)
    cutoff_ar=max(quantile(ft$V28,c(0.5))[1], cutoff)
    cutoff_ng=max(quantile(ft$V29,c(0.5))[1], cutoff)
    cat("Cutoffs to determine normal regions: ", cutoff_phage, cutoff_vf, cutoff_ar, cutoff_ng, cutoff_density, cutoff_distance, fill = TRUE)
    # Exlude ncRNA clusters, or regions with almost all features are 0
    content_normal = subset(gi0, (V5==0&V8>5&V25==0)|(V5==0&V25==0&V26<=cutoff_phage&V27<=cutoff_vf&V28<=cutoff_ar&V29<=cutoff_ng&V30>cutoff_density&V31<=cutoff_distance))
    # content_normal = subset(gi0, (V5==0&V8>5&V25==0)|(V5==0&V25==0&V26<=cutoff_phage&V27<=cutoff_vf&V28<=cutoff_ar&V29<=cutoff_ng))
    gir=subset(gi0, ! V1 %in% as.vector(content_normal$V1))
    ngic=rbind(content_normal)
    if(verbose=="true"){
      cat("new non-GIs: ", dim(content_normal), "\n")
      if("V32" %in% colnames(ft))
      {
        ngic_fn = subset(ngic, V32==1)
        cat("FNs in L-dataset: ",dim(ngic_fn), fill = TRUE)
        print(ngic_fn[,c(1,2,3)])
      }
      if("V33" %in% colnames(ft))
      {
        ngic_fn = subset(ngic, V32==1|V33==1)
        cat("FNs in L-dataset and C-dataset: ",dim(ngic_fn), fill = TRUE)
        print((ngic_fn[,c(1,2,3)]))
      }
    }

    gic = rbind(gi1,gi2)
    if(verbose=="true"){
      cat("new GI: ", dim(gic), "\n")
      if("V32" %in% colnames(ft))
      {
        gic_tp = subset(gic, V32==1)
        cat("new TPs in L-dataset: ",dim(gic_tp), fill = TRUE)
      }
      if("V33" %in% colnames(ft))
      {
        gic_tp = subset(gic, V32==1|V33==1)
        cat("new TPs in L-dataset and C-dataset: ",dim(gic_tp), fill = TRUE)
      }
    }
    gi = rbind(gir, gic)
    ngi2 = subset(ngi1, ! V1 %in% as.vector(gi2$V1))
    ngi = rbind(ngi2, ngic, fte)
  } else {
    gi1 = subset(ngi0,(V5==1|V6==1)&(V25>0))
    ngi1 = subset(ngi0,  ! V1 %in% as.vector(gi1$V1))
    # Regions with mostly phage-related genes, or virulence factors, or antibiotic resistance genes, or novel genes
    gi2 = subset(ngi1, (V26>=0.6|V27>=0.6|V28>=0.6|V29>=0.6|(V26>0&V26+V29>=0.6)))
    gic=rbind(gi1,gi2)
    if(verbose=="true"){
      cat("new GI: ", dim(gic), "\n")
      if("V32" %in% colnames(ft))
      {
        gic_tp = subset(gic, V32==1)
        cat("new TPs in L-dataset: ",dim(gic_tp), fill = TRUE)
      }
      if("V33" %in% colnames(ft))
      {
        gic_tp = subset(gic, V32==1|V33==1)
        cat("new TPs in L-dataset and C-dataset: ",dim(gic_tp), fill = TRUE)
      }
    }
    # Since not all genomic segments are included, it is not meaningful to use median as cutoff here
    content_normal = subset(gi0, (V5==0&V8>5&V25==0)|(V5==0&V25==0&V26<=0.2&V27<=0.2&V28<=0.2&V29<=0.2))
    gir = subset(gi0, ! V1 %in% as.vector(content_normal$V1))
    ngic = rbind(content_normal)
    if(verbose=="true"){
      cat("FNs: ", dim(content_normal), "\n")
      if("V32" %in% colnames(ft))
      {
        ngic_fn = subset(ngic, V32==1)
        cat("FNs in L-dataset: ",dim(ngic_fn), fill = TRUE)
        print(ngic_fn[,c(1,2,3)])
      }
      if("V33" %in% colnames(ft))
      {
        ngic_fn = subset(ngic, V32==1|V33==1)
        cat("FNs in L-dataset and C-dataset: ",dim(ngic_fn), fill = TRUE)
        print(ngic_fn[,c(1,2,3)])
      }
    }

    gi = rbind(gir, gic)
    ngi2 = subset(ngi1, ! V1 %in% as.vector(gi2$V1))
    ngi = rbind(ngi2, ngic, fte)

    cat("Number of Intervals with no genes or many tRNAs(rRNAs): ",dim(fte), fill = TRUE)
  } # for else
}
cat("Number of Intervals: ",dim(ft0), fill = TRUE)
cat("Number of GIs: ",dim(gi), fill = TRUE)
cat("Number of non-GIs: ",dim(ngi), fill = TRUE)
}


############################### write the coordinates to file ########################
fout1=paste(outdir, organism, "_GI",  sep="")
fout2=paste(outdir, organism, "_NonGI",  sep="")

i1=format(gi[c(-1)], digits=3)
i2=format(ngi[c(-1)], digits=3)

write.table(i1, fout1, sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(i2, fout2, sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
