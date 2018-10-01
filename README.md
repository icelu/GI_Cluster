# INTRODUCTION
GI-Cluster is a program for detecting genomic islands (GIs) in a genome by consensus clustering on multiple features.
It includes a sets of scripts for extracting GI-related features from a genome sequence,
performing consensus clustering on the obtained feature matrix to get potential GIs,
and visualizing the predictions.

Author: Bingxin Lu
Affiliation : National University of Singapore
E-mail : bingxin@comp.nus.edu.sg


# CITATION
Bingxin Lu and HonWai Leong (2018). GI-Cluster: detecting genomic islands via consensus clustering on multiple features. Journal of bioinformatics and computational biology, 16(03), 1840010.

# Content of GI-Cluster program
* clustering -- scripts to perform clustering 
* evalution -- data and script for evaluating different GI prediction tools 
* feature -- scripts to extract GI-related features 
* postprocess -- scripts for postprocessing 
* segmentation -- scripts to segment large genome sequences into short intervals
* util -- scripts of general usage 
* visualization -- scripts for visualizing predictions
* Gene_Prediction.sh -- scripts for predicting genes from genome sequences
* GI_Clustering.R -- scripts for running Consensus Clustering on genomic regions with extracted GI-related features
* GI_Comparison.sh -- scripts for visualizing predictions from different methods 
* GI_Feature.sh -- scripts for extracting features related to genomic islands in the unit of genes
* GI_Segmentation.sh -- scripts for splitting a genome sequence into a set of segments
* GI-Cluster.sh  -- the main program
* Segment_Feature.sh -- scripts for extracting features related to genomic islands in a genomic region
* install.sh -- sample commands to insall external tools and packages
* README.md -- this file


# SYSTEM REQUIREMENTS
Unix-based systems.

Since GI-Cluster uses bash scripts to connect each step, it is not convenient to run the program on Windows system.


# SOFTWARE REQUIREMENTS
GI-Cluster is written with Python 2.7, R and bash.
There are several external tools and packages that GI-Cluster depends on.
See install.sh for sample commands.

## Prebuilt database 
* [PFAM](https://pfam.xfam.org/)  
* [RFAM](http://rfam.xfam.org/)  
* [PHAST](http://phast.wishartlab.com)  
* [VFDB](http://www.mgc.ac.cn/VFs/)  
* [CARD](https://card.mcmaster.ca)  
  Download CARD Data on the webpage https://card.mcmaster.ca/download.
* [COG](https://www.ncbi.nlm.nih.gov/COG/)  


## External tools 
### Required
* [CodonW](http://codonw.sourceforge.net/) -- codon analysis
* [Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi) -- database search
* [Hmmmer](http://hmmer.org/) -- database search
* [Circos](http://circos.ca/) -- visualization
* [Infernal cmscan](http://eddylab.org/infernal/) -- ncRNA prediction
* [tRNAscan-SE](http://lowelab.ucsc.edu/tRNAscan-SE/) -- tRNA prediction
* [Repseek](http://wwwabi.snv.jussieu.fr/public/RepSeek/) -- repeat detection
* [COG software](https://www.ncbi.nlm.nih.gov/COG/)

### Optional
* [Prodigal](https://github.com/hyattpd/Prodigal)  -- gene prediction
* [MJSD](https://sourceforge.net/projects/mjsd/) -- genome segmentation
* [GI-SVM](https://github.com/icelu/GI_Prediction/tree/master/GI_SVM) -- genome segmentation
* [Alienhunter](https://www.sanger.ac.uk/science/tools/alien-hunter) -- genome segmentation

Please follow related documenation if you want to use one of these programs.

## Python packages to install 
* [scipy](https://www.scipy.org/) -- used for computing chisquare

## R packages to install 
Required R packages should be installed automatically.
If problems occur, please install manually.


# INSTALLATION
1. Download the source code of GI-Cluster and decompress.
2. Install required tools and packages which are described above.

# USAGE
## INPUT FILE
A FASTA file of genome sequence(s) (ended by '.fna', e.g. NC_010161.fna)

Input when using annotation files from NCBI ("$organism" refers to the name of the organism of interest):
1. The required old NCBI files include:
  - "$organism".fna, genomic DNA sequence,
  - "$organism".ffn, gene sequence -- used for analyzing sequence composition,
  - "$organism".faa, protein sequence -- used for analyzing gene function
2. The required new NCBI files include:
  - "$organism"_genomic.fna, genomic DNA sequence,
  - "$organism"_cds_from_genomic.fna, gene sequence -- used for analyzing sequence composition,
  - "$organism"_protein.faa, protein sequence -- used for analyzing gene function

Input when using custom annotation files:
1. The required files include:
  - "$organism".fna, genomic DNA sequence,
  - "$organism".ffn, gene sequence -- used for analyzing sequence composition,
  - "$organism".faa, protein sequence -- used for analyzing gene function
  - "$organism".glist, gene locations -- a tsv file with 4 columns: ID, start, end, strand (e.g. 1       25      500     F)

  Suppose the header of "$organism".ffn has a format as:
  ">C_RS25945 | C_RS25945 | tRNA/rRNA methyltransferase | 25:711 Reverse", 
  then one can use the following command to get "$organism".glist:   
  `less "$organism".ffn | grep '^>' | cut -d'|' -f4 | sed 's/:/ /g' | awk '{if($3=="Forward") lcol="F"; else lcol="R"; print NR,$1,$2,lcol;}' | sed 's/ /\t/g' > "$organism".glist`


## OUTPUT FILE
There are multiple folders and files generated by the program.
1. Folder $pred_prog -- folders containing gene predictions, gene features and predictions of genomic islands
2. Folder $seg_prog -- folders containing segmentation results and segment features
3. Folder boundary -- files containing predicted tRNAs and repeats in a genome

The predictions of genomic islands are obtained by consensus clustering on multiple features. They depend on several parameters, mainly including feature (the type of features to use for clustering, including "gc", "codon", "kmer", "content", "gc_kmer", "comp", "comp_content"), method (the clustering method), pFeature (the percentage of features to use for clustering), and rep (the number of replicated clusterings). See GI_Clustering.R for all the parameters used in consensus clustering. The corresponding GI predictions are in folders named after the four most import parameters.
Suppose feature=comp_content, method=average, pFeature=1, rep=1, then:
  * When running gene prediction,
  the directory including the final GI candidates is $output_dir/$pred_prog/$seg_prog/$feature/$method/$pFeature/$rep.
  * When not running gene prediction,
  the directory including the final GI candidates is $output_dir/unannotated/$seg_prog/$feature/$method/$pFeature/$rep.
  * The file name including the final GI candidates is merged_"$organism"_refined_GI.

## Running
The programs can be called by following commands.

### Run on complete genome with gene predictions
```
gnome=BtribCIP105476   
organism=NC_010161  
pred_prog=prodigal  
seg_prog=window 
output_dir=/home/b/bingxin/genome/$gnome  
prog_dir=/home/b/bingxin/GI-Cluster 
nohup /usr/bin/time sh $prog_dir/GI-Cluster.sh -s $prog_dir -o $output_dir -n "$organism" -m $seg_prog -p $pred_prog -d 16 > std_"$seg_prog"_"$pred_prog" 2>&1 &
```

### Run on complete genome without gene predictions
```
gnome=BtribCIP105476  
organism=NC_010161  
pred_prog=none  
seg_prog=window 
output_dir=/home/b/bingxin/genome/$gnome  
prog_dir=/home/b/bingxin/GI-Cluster 
nohup /usr/bin/time sh $prog_dir/GI-Cluster.sh -s $prog_dir -o $output_dir -n "$organism" -m $seg_prog -p $pred_prog -d 16 -t 0 > std_"$seg_prog"_"$pred_prog"_unannotated 2>&1 &
```

### Run on incomplete genome with gene predictions
```
gnome=Vibrio_cholerae_RC9_uid55789  
organism=NZ_ACHX00000000  
pred_prog=prodigal  
seg_prog=window 
output_dir=/home/b/bingxin/genome/incomplete/$gnome 
prog_dir=/home/b/bingxin/GI-Cluster 
nohup /usr/bin/time sh $prog_dir/GI-Cluster.sh -s $prog_dir -o $output_dir -n "$organism" -m $seg_prog -p $pred_prog -d 16 -b 1 > std_"$seg_prog"_"$pred_prog" 2>&1 &
```

### Run on incomplete genome without gene predictions
```
gnome=Vibrio_cholerae_RC9_uid55789  
organism=NZ_ACHX00000000  
pred_prog=none  
seg_prog=window 
output_dir=/home/b/bingxin/genome/incomplete/$gnome 
prog_dir=/home/b/bingxin/GI-Cluster 
nohup /usr/bin/time sh $prog_dir/GI-Cluster.sh -s $prog_dir -o $output_dir -n "$organism" -m $seg_prog -p $pred_prog -d 16 -b 1 -t 0 > std_"$seg_prog"_"$pred_prog"_unannotated 2>&1 &
```
When using the old NCBI annotation files, please use option ""-p ncbi_old", namely pred_prog=ncbi_old.

### Run on complete genome with custom gene predictions
```
gnome=BtribCIP105476   
organism=NC_010161  
pred_prog=custom  
seg_prog=window 
output_dir=/home/b/bingxin/genome/$gnome  
prog_dir=/home/b/bingxin/GI-Cluster 
nohup /usr/bin/time sh $prog_dir/GI-Cluster.sh -s $prog_dir -o $output_dir -n "$organism" -m $seg_prog -p $pred_prog -d 16 > std_"$seg_prog"_"$pred_prog" 2>&1 &
```


### Note: 
Suppose the input files are in $output_dir, the names of the input genome file should be "$organism".fna.
* When "pred_prog=ncbi_old", the names of the input file should be "$organism".fna,  "$organism".faa, "$organism".ffn.
* When "pred_prog=ncbi", the names of the input file should be "$organism".fna, "$organism"_protein.faa, "$organism"_cds_from_genomic.fna.

In term of running time, GI-Cluster is fast in most steps, except database searching and consensus clustering. It may take a long time to find novel genes when searching against COG databases. 

There are multiple intermediate files generated. If you suspect some files are not correct, please delete them and rerun GI_Cluster.sh. Or else, the program will assume the current files are correct and go on to the next steps.
