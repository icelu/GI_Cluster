# Perform 1D clustering on each feature;
# Generate a consensus matrix (CM) on all these clusterings;
# Run clustering on CM to get final groups
#
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg

library(Ckmeans.1d.dp)
library(mclust)

# x -- input data matrix, clmem -- cluster assignment
getConn<-function(x, clmem){
  # create a current instance of the zeroes indexed array (without pseudocount)
  conn <- matrix(0,nrow=dim(x)[1]+1,ncol=dim(x)[1]+1,dimnames=list(c('x',row.names(x)),c('x',row.names(x))))
  # put in the cluster memberships in the x row,col
  conn[row.names(clmem),1] <- as.matrix(clmem)
  conn[1,row.names(clmem)] <- as.matrix(clmem)

  for(i in 2:dim(conn)[1]){
    for(j in 2:dim(conn)[1]){
      if((conn[i,1] == conn[1,j]) & conn[i,1] != 0){
        conn[i,j] <- 1
        conn[j,i] <- 1
      }
      else{
        conn[i,j] <-0
        conn[j,i] <-0
      }
    }
  }
  # cut to the connectivity matrix
  final <- conn[2:dim(conn)[1],2:dim(conn)[1]]

  return(final)
}

# Sampling both items and features
# Function call:
# sample_x = sampleCols( d, pCol, pRow, weightsCol, weightsRow )
# d -- a data frame
# pCol -- feature
# pRow -- sample
# weightsCol -- probability vector for each column
# weightsRow -- probability vector for each row
# Return a list with the sampled sub-matrix, sampled rows and sampled columns
sampleCols <- function( d,
                        pCol=NULL,
                        pRow=NULL,
                        weightsCol=NULL,
                        weightsRow=NULL ){
  space <- ifelse( inherits( d, "dist" ), ncol( as.matrix(d) ), ncol(d) )
  sampleN <- floor(space*pCol)
  sampCols <- sort( sample(space, sampleN, replace = FALSE, prob = weightsCol) )
  this_sample <- sampRows <- NA

  if ( (!is.null(pRow)) &&
       ((pRow <= 1) || (! is.null( weightsRow ) ) ) ) {
    space = nrow(d)
    sampleN = floor(space*pRow)
    sampRows = sort( sample(space, sampleN, replace = FALSE, prob = weightsRow) )
    this_sample <- d[sampRows,sampCols]
  }
  else{
    this_sample <- d[,sampCols]
  }

  return( list(submat=this_sample,
               subrows=sampRows,
               subcols=sampCols ) )
}



# data: row -- item, column -- feature
clustering.1D <- function(data, clnum=3, algorithms=list('Ckmeans.1d.dp','densityMclust'),alparams=list(), pCol=1, pRow=1, weightsCol=NULL, weightsRow=NULL, reps=2, verbose="false", outdir="", ft=data){
  if(class(algorithms)!='list'){stop('You must pass the algorithm names in a list')};
  # check that all of the algorithms specified are available
  for(i in 1:length(algorithms)){
    if(existsFunction(algorithms[[i]])!=1){stop('One or more of the specified algorithms does not exist')} else {}
  }

  # list to hold all of the consmatrix objects
  clmem_list <- list();
  # for each algorithm
  for(a in 1:length(algorithms)){
    algo = algorithms[[a]];
    # message("Algorithm ", algo)
    # find the clustering for each data, store them into a list
    for(d in 1:length(data))
    {
      # message("Feature ", d)
      x=data[,d];

      if(length(alparams)!=0){
        current_params <- alparams[[a]];
      }
      else{
        current_params = list();
      }
      algo_clmem <- paste(algo,'_clmem',sep='');
      current_algo <- get(algo_clmem);
      cls_res <- current_algo(x, clnum, params=current_params);
      clmem <- data.frame(cls_res$clmem);
      # move the row.names across
      row.names(clmem) <- row.names(data)
      # tidy the column header
      names(clmem) <-c('cm');

      if (verbose=="true"){
        message("cluster assignment:")
        print(table(clmem));
        # print the clustering results for check
        # find the number of real clusters
        if(outdir!=""){
          num_clust <- unique(clmem);
          for(c in 1:nrow(num_clust)){
            cid <- num_clust[c,];
            Clab= as.vector(t(clmem));
            cls=which(grepl(cid,Clab));
            segs = ft[cls,];
            fsegs=format(segs, digits=3)
            fout=paste(outdir, 'f',d,'_',algo,'_k',cid,  sep="")
            if(!file.exists(fout)){
              write.table(fsegs, fout, sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
            }
          }
        }
      }

      final = getConn(data, clmem);

      id = final;
      id[] <- 0;
      id[row.names(clmem),row.names(clmem)] <- 1;
      key = paste('f',d,'_',algo,'_k',clnum,sep='');
      clmem_list[[key]]=list('conn'=final, 'id'=id);
    }
  }

  cm_list = list()
  for(i in 1:reps){
    if (verbose=="true"){
      message("Iteration ", i)
    }
    # here features are in columns
    if(pCol<1 || pRow<1){
      samples = sampleCols(data, pCol=pCol, pRow=pRow, weightsCol=weightsCol, weightsRow=weightsRow)$subcols;
    }
    else{
      samples = seq(1, ncol(data))
    }
    if (verbose=="true"){
      message("Selected features ", str(samples))
    }
    final_conn <- matrix(0,nrow=dim(data)[1],ncol=dim(data)[1],dimnames=list(row.names(data),row.names(data)));
    final_id <- final_conn;
    for(a in 1:length(algorithms)){
      algo = algorithms[[a]];
      for(d in 1:length(samples))
      {
        feature = samples[[d]]
        # retrive the clustering and id
        key = paste('f',feature,'_',algo,'_k',clnum,sep='');

        conn_id = clmem_list[[key]];
        final = conn_id[['conn']];
        id = conn_id[['id']];

        final_conn = final_conn+final;
        final_id = final_id+id;
      }
    }
    # this is the consensus matrix for one clustering algorithm
    consensus = final_conn/final_id;
    # print(str(consensus))
    # check for NAs (this is safer than pseudo counting, replace NAs with 0 i.e. they are never drawn together and/or never connected)
    ind <- is.na(consensus);
    consensus[ind] <- 0;
    key = paste('rep',i,'_k',clnum,sep='');
    cm_list[[key]] <- consensus
  }
  return(cm_list);
}
