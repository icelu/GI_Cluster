#!/usr/bin/env Rscript

# This script provides interfaces to call different clustering algorithms in R
# This is based on implementations of clusterCons (https://github.com/cran/clusterCons)

# To include kmeans, pam, clara
if(!require(cluster)){
    install.packages("cluster")
    library(cluster)
}
# Affinity Propagation Clustering
if(!require(apcluster)){
    install.packages("apcluster")
    library(apcluster)
}
# For model-based clustering
if(!require(mclust)){
    install.packages("mclust")
    library(mclust)
}
# To include kmeanspp, or kmeans++
if(!require(LICORS)){
    install.packages("LICORS")
    library(LICORS)
}
# To include Kernel Clustering
if(!require(kernlab)){
    install.packages("kernlab")
    library(kernlab)
}

# For CLUster Ensembles
if(!require(clues)){
    install.packages("clues")
    library(clues)
}

#agnes
agnes_clmem <- function(x,clnum,params=list()){
	#assign the parameters to the data.frame
	params$x <- x
	#perform the clustering
	hc <- do.call(agnes,params)
	#grab the membership object
	clmem <- data.frame(cutree(hc,clnum))
	#move the row.names across
	row.names(clmem) <- row.names(x)
	#tidy the column header
	names(clmem) <-c('cm')
	#return the membership object
	res_list=list('clmem'=clmem,'clust_obj'=hc);
	#return(clmem);
	return(res_list);
}
#examples
#1. with nothing but the data.frame and the clnum
#agnes_clmem(cluster_example_data2,clnum=2)
#2. with both specified (k will be ignored as this will be set to the current clnum in the loop)
#agnes_clmem(cluster_example_data2,clnum=2,params=list(metric='euclidean',method='average'))

#pam
#be careful to use a valid metric name for pam as the function does not check and will default to euclidean
pam_clmem <- function(x,clnum,params=list()){
	#the data is not a dissimilarity matrix so :-
	#check if params have been passed if not set up the default
	#the k value should not be specfied as this comes from the loop
	if(length(params)==0){
		#pam will automatically generate a distance matrix from a data frame, but here we explicitly state the distance metric for clarity
		params <- list(metric='euclidean',k=clnum)
	}
	else{
		#If k is set then replace with clnum. If diss has not been set default to 'euclidean'
		#it checks the list with regexpr checks if any are equal to 1 (T/F) return then sums the logicals
		if(sum(regexpr('^metric$',names(params))==1)<1){
			params$metric='euclidean'
		}
		params$k <- clnum
	}
	params$x <- x;
	pam_obj <- do.call(pam,params);
	clmem <- data.frame(pam_obj$clustering);
	names(clmem) <- c('cm');
	res_list=list('clmem'=clmem,'clust_obj'=pam_obj);
	#return(clmem);
	return(res_list);
}
#examples
#1. with no dissimilarity matrix or cluster number specified
#pam_clmem(cluster_example_data2,clnum=2)
#2. with both specified (k will be ignored as this will be set to the current clnum in the loop)
#pam_clmem(cluster_example_data2,clnum=2,params=list(metric='euclidean',k=2))

#kmeans
kmeans_clmem <- function(x,clnum,params=list()){
	#remove the diss argument as kmeans doesn't take this argument
	params$diss=NULL;
	#set the MacQueen algorithm as default we have seen random failures with Hartigan-Wong
	params$algorithm='Mac';
	params$nstart=25
	#if object is dist then convert back to matrix
	if(class(x)=='dist'){
		x <- as.matrix(x);
	}
	if(length(params)==0){
		params <- list(centers=clnum);
	}
	else{
		#add or override the centers number in existing params to the current clnum (resolves potential conflict)
		params$centers=clnum;
	}
	params$x <- x;
	km <- do.call(kmeans,params);
	clmem <- data.frame(km$cluster);
	names(clmem) <- c('cm');

	res_list=list('clmem'=clmem,'clust_obj'=km);
	#return(clmem);
	return(res_list);
}
#examples
#1. simple
#kmeans_clmem(cluster_example_data2,clnum=4)
#2. with change of kmeans algorithm, note clnum overides centers
#kmeans_clmem(cluster_example_data2,clnum=6,params=list(algorithm="Lloyd",centers=23))


#kmeans++
kmeanspp_clmem <- function(x,clnum,params=list()){
	params$nstart=25;
	#add or override the centers number in existing params to the current clnum (resolves potential conflict)
	if(length(params)==0){
	  params <- list(k=clnum);
	}
	else{
	  #add or override the centers number in existing params to the current clnum (resolves potential conflict)
	  params$k=clnum;
	}
	params$data <- x;
	km <- do.call(kmeanspp, params);
	clmem <- data.frame(km$cluster);
	names(clmem) <- c('cm');

	res_list=list('clmem'=clmem,'clust_obj'=km);
	#return(clmem);
	return(res_list);
}



# Kernel k-means
kkmeans_clmem <- function(x,clnum,params=list()){
	params$kernel = "vanilladot";

	if(length(params)==0){
		params <- list(centers=clnum);
	}
	else{
		#add or override the centers number in existing params to the current clnum (resolves potential conflict)
		params$centers=clnum;
	}
	params$x <- x;
	km <- do.call(kkmeans, params);
	clmem <- data.frame(km@.Data);
	names(clmem) <- c('cm');

	res_list=list('clmem'=clmem,'clust_obj'=km);
	#return(clmem);
	return(res_list);
}

#hclust
hclust_clmem <- function(x,clnum,params=list()){
	#remove any erroneous diss parameters, hclust only takes 'dist' objects
	params$diss=NULL;
	#if the data is not a dissimilarity already make one. The user should not use hclust without providing the
	#distance matrix as it only takes distance matrices.
	if(class(x)!='dist'){
		dm <- dist(x);
		params$d=dm;
	}
	else{
		dm <- x;
		params$d=dm;
	}
	#perform the clustering
	hc <- do.call(hclust,params);
	#grab the membership object
	clmem <- data.frame(cutree(hc,clnum));
	#tidy the column header
	names(clmem) <-c('cm');
	# print(str(clmem));
	res_list=list('clmem'=clmem,'clust_obj'=hc);
	#return(clmem);
	return(res_list);
}
#examples
#1. simple
#hclust_clmem(cluster_example_data2,clnum=4)
#2. with change of distance algorithm
#hclust_clmem(cluster_example_data2,clnum=6,params=list(diss='manhattan'))

#diana
#be careful to use a valid metric name for pam as the function does not check and will default to euclidean
diana_clmem <- function(x,clnum,params=list()){
	#the data is not a dissimilarity matrix so :-
	#check if params have been passed if not set up the default
	if(length(params)!=0){
		#check whether diss has been set, if not then set to default 'euclidean'
		if(sum(regexpr('^metric$',names(params))==1)<1){
			params$metric='euclidean'
		}
	}
	#pam will automatically generate a distance matrix from a data frame, but here we explicitly state the distance metric for clarity
	else{
		params <- list(metric='euclidean')
	}
	params$x <- x;
	hc <- do.call(diana,params);
	#grab the membership object
	clmem <- data.frame(cutree(hc,clnum))
	#move the row.names across
	row.names(clmem) <- row.names(x)
	#tidy the column header
	names(clmem) <-c('cm')
	#return the membership object
	res_list=list('clmem'=clmem,'clust_obj'=hc);
	#return(clmem);
	return(res_list);
}
#examples
#1. simple
#diana_clmem(cluster_example_data2,clnum=4)
#2. with change of distance algorithm
#diana_clmem(cluster_example_data2,clnum=6,params=list(metric='manhattan'))

#apcluster
#clustering using the affinity propagation method of Frey and Dueck Science (2007)
#using the apclusterK method from the apcluster package that allows you to specify the desired cluster number
apcluster_clmem <- function(x,clnum,params=list()){
	#convert the data.frame to a matrix
	x <- as.matrix(x);
	#the diss flag is unecessary so remove
	params$diss <- NULL;
	#note that apcluster requires a square matrix entry, in reality this means a similarity matrix
	#have to assume the user knows how to use the method so just check that the input is a square matrix
	if(class(x)=='matrix' ){
	   # && dim(x)[1] == dim(x)[2]){
		#ensure that the process returns the desired number of clusters
		params$prc=0;
		#add the data to the params
		params$s <- negDistMat(x);
		params$x <- x;
		#specify the cluster number
		params$K <- clnum;
		#perform the clustering
		apc <- (do.call(apclusterK,params));
		#grab the membership objects from the cluster list
		cl_list <- apc@clusters;
		#reconfigure the output to a standard cluster membership list
		for(i in 1:length(cl_list)){
			current <- names(cl_list[[i]]);
			for(j in 1:length(current)){
				if(i==1 && j==1){
					clmem <- data.frame(i);
					names(clmem) <- c('cm');
					rownames(clmem)<-current[j];
				}
				else{
					#add the next row
					add<- data.frame(i);
					names(add) <- c('cm');
					rownames(add) <- current[j];
					clmem <- rbind(clmem,add);
				}
			}
		}
		res_list=list('clmem'=clmem,'clust_obj'=apc);
		#return(clmem);
		return(res_list);
	}
	else{
		stop('apcluster requires a square numeric matrix for clustering, this is supposed to be a similarity matrix');
	}
}
#examples
#1. simple
#apcluster_clmem(data,clnum=4);
#2. increasing the number of bisection steps to perfrom
#apcluster_clmem(data,clnum=4,params=list(bimaxit=100));

#OPTIONAL ADDITIONAL FUNCTIONS (uncomment for use)

##ADD YOUR OWN ALGORITHM CALLING FUNCTIONS HERE

##generic new function
#new_clmem(x,clnum,params=list()){
#	#perform parameter checks to ensure that the call is consitent with the algorithm function
#	#assign the paramters to the data frame
#	params$x <- x;
#	#perform the clustering
#	cl_obj <- data.frame(do.call(new,params));
#	#extract or convert to the cluster membership object and convert to data.frame for example :-
#	#for direct assignment algorithms
#	clmem <- data.frame(cl_obj$clustering) # or whatever the field name for the assuigned clusters is in the returned object
#	#or for heirarchical cluster results
#	#clmem <- data.frame(cutree(cl_obj,clnum));
#	names(clmem) <- c('cm');
#	return(clmem);
#}


#dbscan
dbscan_clmem <- function(x,clnum,params=list()){
#perform parameter checks to ensure that the call is consitent with the algorithm function
#assign the paramters to the data frame
params$data <- x;
params$eps <- 5;
#perform the clustering
dbscan_obj <- (do.call(dbscan,params));
#extract or convert to the cluster membership object and convert to data.frame for example :-
#for direct assignment algorithms
clmem <- data.frame(dbscan_obj$cluster); # or whatever the field name for the assuigned clusters is in the returned object
#or for heirarchical cluster results
#clmem <- data.frame(cutree(cl_obj,clnum));
row.names(clmem) <- row.names(x);
names(clmem) <- c('cm');
res_list=list('clmem'=clmem,'clust_obj'=dbscan_obj);
#return(clmem);
return(res_list);
}


# Mclust
Mclust_clmem <- function(x,clnum,params=list()){
  #perform parameter checks to ensure that the call is consitent with the algorithm function
  #assign the paramters to the data frame
  params$data <- x;
  #params$eps <- 5;
  #perform the clustering
  mclust_obj <- (do.call(Mclust,params));
  #extract or convert to the cluster membership object and convert to data.frame for example :-
  #for direct assignment algorithms
  clmem <- data.frame(mclust_obj$classification); # or whatever the field name for the assuigned clusters is in the returned object
  row.names(clmem) <- row.names(x);
  names(clmem) <- c('cm');

  res_list=list('clmem'=clmem,'clust_obj'=mclust_obj);
  #return(clmem);
  return(res_list);
}


#clara
clara_clmem <- function(x,clnum,params=list()){
  #remove the diss argument as kmeans doesn't take this argument
  params$diss=NULL;
  params$stand=TRUE;
  #if object is dist then convert back to matrix
  if(class(x)=='dist'){
    x <- as.matrix(x);
  }
  if(length(params)==0){
    params <- list(centers=clnum);
  }
  else{
    #add or override the centers number in existing params to the current clnum (resolves potential conflict)
    params$k=clnum;
  }
  params$x <- x;
  clara_obj <- do.call(clara,params);
  clmem <- data.frame(clara_obj$clustering);
  names(clmem) <- c('cm');

  res_list=list('clmem'=clmem,'clust_obj'=clara_obj);
  #return(clmem);
  return(res_list);
}


clues_clmem <- function(x,clnum,params=list()){
  #perform parameter checks to ensure that the call is consitent with the algorithm function
  #assign the paramters to the data frame
  params$y <- x;
  params$n0 <- clnum;
  #perform the clustering
  clues_obj <- (do.call(clues,params));
  #extract or convert to the cluster membership object and convert to data.frame for example :-
  #for direct assignment algorithms
  clmem <- data.frame(clues_obj$mem); # or whatever the field name for the assuigned clusters is in the returned object
  row.names(clmem) <- row.names(x);
  names(clmem) <- c('cm');

  res_list=list('clmem'=clmem,'clust_obj'=clues_obj);
  #return(clmem);
  return(res_list);
}


#kmeans
Ckmeans.1d.dp_clmem <- function(x,clnum,params=list()){
  #remove the diss argument as kmeans doesn't take this argument
  params$diss=NULL;

  if(length(params)==0){
    params <- list(k=clnum);
  }
  else{
    #add or override the centers number in existing params to the current clnum (resolves potential conflict)
    params$k=clnum;
  }
  params$x <- x;
  km <- do.call(Ckmeans.1d.dp,params);
  clmem <- data.frame(km$cluster);
  names(clmem) <- c('cm');

  res_list=list('clmem'=clmem,'clust_obj'=km);
  #return(clmem);
  return(res_list);
}

# Mclust
densityMclust_clmem <- function(x,clnum,params=list()){
  #perform parameter checks to ensure that the call is consitent with the algorithm function
  #assign the paramters to the data frame
  params$data <- x;
  #params$eps <- 5;
  #perform the clustering
  mclust_obj <- (do.call(densityMclust, params));
  #extract or convert to the cluster membership object and convert to data.frame for example :-
  #for direct assignment algorithms
  clmem <- data.frame(mclust_obj$classification); # or whatever the field name for the assuigned clusters is in the returned object
  row.names(clmem) <- row.names(x);
  names(clmem) <- c('cm');

  res_list=list('clmem'=clmem,'clust_obj'=mclust_obj);
  #return(clmem);
  return(res_list);
}
