# Data preparation -------------------------------------------------------------

.onLoad <- function(libname, pkgname) {
    data(tcga, package = pkgname, envir = parent.env(environment()))
    data(ccle, package = pkgname, envir = parent.env(environment()))
    data(gene_info, package = pkgname, envir = parent.env(environment()))
    data(map.ensembl_symbol_hg19, package = pkgname, envir = parent.env(environment()))
    data(map.ensembl_symbol_hg38, package = pkgname, envir = parent.env(environment()))
}

.translate <- function(eset, refset, genome = "hg38"){

  if (genome == "hg19"){
    map <- get("map.ensembl_symbol_hg19")
  } else if (genome == "hg38"){
    map <- get("map.ensembl_symbol_hg38")
  } else {
    stop("Unsupported genome")
  }

  map <- map[map$external_gene_name!="",]
  expr <- exprs(eset)
  expr <- merge(expr,unique(map[,c("external_gene_name","ensembl_gene_id")]),by.x="row.names",by.y="external_gene_name")
  expr <- expr[expr$ensembl_gene_id %in% rownames(refset),]

  duplicate_genes <- names(which(table(expr$Row.names)>1))
  if (length(duplicate_genes)>0){
    idx_duplicate <- expr$Row.names %in% duplicate_genes
    expr.non_duplicate <- expr[!idx_duplicate,]
    expr.duplicate <- list()
    for (dup_gene in duplicate_genes){
      expr.duplicate[[dup_gene]] <- expr[expr$Row.names==dup_gene,][which.max(expr[expr$Row.names==dup_gene,]$count),]
    }
    expr.duplicate <- do.call("rbind",expr.duplicate)
    expr <- rbind(expr.non_duplicate,expr.duplicate)
  }
  expr$Row.names <- NULL
  rownames(expr) <- expr$ensembl_gene_id
  expr$ensembl_gene_id <- NULL

  eset <- .matrix_to_eset(expr)

  return(eset)
}

.deduplicate_map <- function(map){
  duplicate_genes <- names(which(table(map$gene_id)>1))
  idx_duplicate <- map$gene_id %in% duplicate_genes
  map.non_duplicate <- map[!idx_duplicate,]
  map.duplicate <- list()
  for (dup_gene in duplicate_genes){
    map.duplicate[[dup_gene]] <- max(map[map$gene_id==dup_gene,]$transcript_length)
  }
  map.duplicate <- do.call("rbind",map.duplicate)
  map.duplicate <- data.frame(gene_id=rownames(map.duplicate),
                              transcript_length=map.duplicate)
  map <- rbind(map.non_duplicate,map.duplicate)
  return(map)
}

.get_gene_info <- function(genes, genome, symbols = F){

  if (any(table(genes)>1)){
    stop("Duplicate gene names in input are not supported.")
  }

  if (genome == "hg19"){
    map <- get("map.ensembl_symbol_hg19")
  } else if (genome == "hg38"){
    map <- get("map.ensembl_symbol_hg38")
  } else {
    stop("Unsupported genome")
  }

  if (substr(head(genes)[1],1,4) != "ENSG"){
    stop("If symbols = F, provided gene names need to be human Ensembl identifiers (ENSG...)")
  }
  map <- map[which(map$ensembl_gene_id %in% genes),]
  map$external_gene_name <- NULL
  colnames(map)[colnames(map)=="ensembl_gene_id"] <- "gene_id"
  map <- unique(map)

  map <- .deduplicate_map(map)

  genes_isect <- intersect(map$gene_id , genes)

  lengths <- sapply(genes_isect, function(x) {
    max(map$transcript_length[map$gene_id == x])
  })
  #stopifnot(identical(names(lengths), genes_ens))
  stopifnot(identical(names(lengths), genes_isect))
  g_info <- as.data.frame(lengths)
  colnames(g_info) <- "length"

  return(g_info)
}

.rpkm <- function(x_col,expr.len){
  #% RPKM normalization - non-rubost version (as opposed to the TCGA data - but does not affect correlations)
  depth <- sum(x_col)
  x_col.rpm <- (x_col/depth)*10^6
  x_col.rpkm <- (x_col.rpm/expr.len)*10^3
  return(x_col.rpkm)
}

.rpkm_eset <- function(eset, gene_info){
  stopifnot(identical(rownames(eset),rownames(gene_info)))
  exprs(eset) <- apply(exprs(eset),2,function(x){.rpkm(x,gene_info$length)})
  return(eset)
}

.matrix_to_eset <- function(mat){
  sample_info <- data.frame(sample_barcode=colnames(mat),stringsAsFactors = F)
  rownames(sample_info) <- sample_info$sample_barcode
  gene_info <- data.frame(ens_id=rownames(mat),stringsAsFactors = F)
  rownames(gene_info) <- gene_info$ens_id
  pdat <- Biobase::AnnotatedDataFrame(data=sample_info)
  fdat <- Biobase::AnnotatedDataFrame(data=gene_info)
  stopifnot(identical(rownames(mat),rownames(fdat)))
  stopifnot(identical(colnames(mat),rownames(pdat)))
  eset <- Biobase::ExpressionSet(assayData=as.matrix(mat),phenoData=pdat,
                                 featureData=fdat)
  return(eset)
}

.deduplicate_htseq <- function(y){
  duplicate_genes <- names(which(table(y$gene_id)>1))
  idx_duplicate <- y$gene_id %in% duplicate_genes
  y.non_duplicate <- y[!idx_duplicate,]
  y.duplicate <- list()
  for (dup_gene in duplicate_genes){
    y.duplicate[[dup_gene]] <- max(y[y$gene_id==dup_gene,]$count)
  }
  y.duplicate <- do.call("rbind",y.duplicate)
  y.duplicate <- data.frame(gene_id=rownames(y.duplicate),
                            count=y.duplicate)
  y <- rbind(y.non_duplicate,y.duplicate)
  rownames(y) <- y$gene_id
  y$gene_id <- NULL
  return(y)
}

.read_htseq <- function(indir=NULL,file=NULL,extension = ".gene_counts", recursive = T){
  if (!is.null(indir) && is.null(file)){
    fnames <- list.files(indir,pattern=extension,full.names = T,recursive = recursive)
  } else if (is.null(indir) && !is.null(file)){
    fnames <- file
  } else {
    stop("Either indir or file must be specified, but not both.")
  }
  expr <- list()
  for (fname in fnames){
    sname <- stringr::str_split_fixed(basename(fname),"\\.",2)[,1]
    expr[[sname]] <- as.data.frame(data.table::fread(fname),stringsAsFactors=F)
    expr[[sname]] <- expr[[sname]][1:(nrow(expr[[sname]])-5),]
    expr[[sname]]$V1 <- stringr::str_split_fixed(expr[[sname]]$V1,"\\.",2)[,1]
  }
  stopifnot(length(unlist(lapply(expr,function(x){grep(x$V1,pattern="^_")})))==0)

  expr <- lapply(expr,function(x){
    y <- data.frame(count=x$V2,stringsAsFactors = F);
    y$gene_id <- x$V1;
    y <- .deduplicate_htseq(y);
    return(y)
  })
  snames <- names(expr)
  expr <- do.call("cbind",expr)
  colnames(expr) <- snames

  eset <- .matrix_to_eset(expr)
  return(eset)
}

.read_matrix <- function(pth,sep="\t",header = T, row.names = 1){
  expr <- read.table(file = pth, sep = sep, header = header, row.names = row.names)
  eset <- .matrix_to_eset(expr)
  return(eset)
}

.setup_expression <- function(input, extension = ".gene_counts", is_matrix = F){
  if (class(input)=="character"){
    if (file.exists(input) && !dir.exists(input) & is_matrix){ # matrix file
      print("Assuming path to matrix with unnormalized read counts")
      eset <- .read_matrix(input)
    } else if (file.exists(input) && dir.exists(input)){ # htseq-count files
      print(paste0("Assuming htseq-count files with the extension ",extension))
      eset <- .read_htseq(indir=input, extension = extension)
    } else if (file.exists(input) && !dir.exists(input)  & !is_matrix){ # single htseq-count file
      print("Assuming single htseq-count file")
      eset <- .read_htseq(file=input, extension = extension)
    } else if (input %in% c("tcga","ccle")){ # if neither path to directory nor matrix file
      eset <- get(input)
    } else {
      stop("Input string is neither path to directory nor file and does not represent an included dataset.")
    }
  } else if (class(input)=="ExpressionSet"){
    eset <- input
  } else if ("matrix" %in% class(input)) {
    print("Assuming matrix with unnormalized read counts")
    eset <- .matrix_to_eset(input)
  } else {
    stop("Unsupported input class")
  }
  return(eset)
}

.subset_samples <- function(refset,subset_to){
  if (length(subset_to)>0){
    refset <- refset[,which(Biobase::pData(refset)$cancer_type %in% subset_to)]
  } else {
    stop("Category to subset not specified")
  }
  return(refset)
}

.write_correlation_table <- function(r, inhouse, refset, outdir){
  expr_inhouse.data_norm <- Biobase::exprs(inhouse)
  tab <- c()
  for (j in 1:ncol(r)){
    rw <- colnames(expr_inhouse.data_norm)[j]
    tab <- c(tab,rw)
    print(rw)
    idx <- order(r[,j],decreasing=T)
    for (i in 1:10){
      rw <- paste0(i,": ","cor. = ",r[idx[i], j], ": ",Biobase::pData(refset)$cancer_type[idx[i]], ": ",
                   colnames(refset)[idx[i]])
      tab <- c(tab,rw)
      print(rw)
    }
  }
  write.table(tab,file=paste0(outdir,'top10_similar.txt'),quote=F,row.names = F,col.names = F)
}

# Calculations -----------------------------------------------------------------

.calculate_correlations <- function(inhouse,refset){
  stopifnot(identical(rownames(inhouse),rownames(refset)))

  expr_inhouse.data_norm.coding <- Biobase::exprs(inhouse)
  expr_pancan.coding <- Biobase::exprs(refset)
  stopifnot(identical(rownames(expr_inhouse.data_norm.coding),rownames(expr_pancan.coding)))

  r <- t(cor(expr_inhouse.data_norm.coding ,expr_pancan.coding, method = "spearman"))

  return(r)
}

.knn <- function(r, refset, k = 6){
  options(warn=2)
  stopifnot(nrow(Biobase::pData(refset))==nrow(r))
  predicted <- c()
  for (j in 1:ncol(r)){
    for (kk in seq(from=k,to=1,by=-1)){
      idx <- order(r[,j],decreasing = T)
      x <- Biobase::pData(refset)$cancer_type[idx[1:kk]]
      y <- unique(x)
      n <- rep(0,length(y))
      for (iy in 1:length(y)){
        n[iy] <- length(which(y[iy]==x))
      }
      mx <- max(n)
      itemp <- which(n==mx)
      if (length(itemp)==1){
        predicted[j] <- y[itemp]
        break
      }
    }
  }
  options(warn=1)
  tab <- data.frame(Sample=colnames(r),predicted)
  return(tab)
}

.cv.eval.k <- function(r.tcga,refset.custom,params,k){
  res.class.tcga <- list()
  snames <- colnames(r.tcga)
  for (cn in params){
    pData(refset.custom)$cancer_type <- pData(refset.custom)[[cn]]
    res.class.tcga[[cn]] <- .knn.cv(r = r.tcga,refset = refset.custom, k = k)$tab
    rownames(res.class.tcga[[cn]]) <- res.class.tcga[[cn]]$Sample
    res.class.tcga[[cn]] <- res.class.tcga[[cn]][snames,]
    res.class.tcga[[cn]]$Sample <- NULL
    colnames(res.class.tcga[[cn]]) <- cn
  }
  res.class.tcga <- do.call("cbind",res.class.tcga)

  cv.tab <- list()
  for (cn in colnames(res.class.tcga)){
    tmp_1 <- res.class.tcga[,cn]
    names(tmp_1) <- rownames(res.class.tcga)
    tmp_2 <- pData(refset.custom)[[cn]]
    names(tmp_2) <- rownames(pData(refset.custom))
    stopifnot(identical(names(tmp_1),names(tmp_2)))
    cv.tab[[cn]] <- caret::confusionMatrix(factor(tmp_1),factor(tmp_2))
  }

  cv.df <- do.call("rbind",lapply(cv.tab,function(x){x$overall}))
  cv.df.melt <- reshape2::melt(cv.df)

  return(list(res.class.tcga=res.class.tcga,cv.df=cv.df,cv.df.melt=cv.df.melt))
}

.knn.cv <- function (r, refset, k = 6)
{
  options(warn = 2)
  r.base <- r
  refset.base <- refset
  stopifnot(ncol(r.base)==nrow(r.base))
  stopifnot(identical(colnames(r.base),rownames(r.base)))
  stopifnot(nrow(Biobase::pData(refset)) == nrow(r))
  predicted <- c()
  for (j in 1:ncol(r)) {
    r <- r.base[-j,]
    refset <- refset.base[,-j]
    for (kk in seq(from = k, to = 1, by = -1)) {
      idx <- order(r[, j], decreasing = T)
      x <- Biobase::pData(refset)$cancer_type[idx[1:kk]]
      y <- unique(x)
      n <- rep(0, length(y))
      for (iy in 1:length(y)) {
        n[iy] <- length(which(y[iy] == x))
      }
      mx <- max(n)
      itemp <- which(n == mx)
      if (length(itemp) == 1) {
        predicted[j] <- y[itemp]
        break
      }
    }
  }
  options(warn = 1)
  tab <- data.frame(Sample = colnames(r), predicted)
  return(list(predicted = predicted, tab = tab))
}

.tsne_core <- function(y_combined, outdir, plot_title, niter = 1000, perplexity = 30){

  tsne_combined <- Rtsne::Rtsne(t(Biobase::exprs(y_combined)), dims = 2, perplexity=perplexity,
                                verbose=TRUE, max_iter = niter, num_threads = 50, partial_pca = TRUE)

  # Assign parameters for plotting
  df_combined <- as.data.frame(tsne_combined$Y,stringsAsFactors=F)
  rm(tsne_combined)
  colnames(df_combined) <- c("Dimension 1","Dimension 2")
  df_combined$color <- Biobase::pData(y_combined)$color

  df_combined$sample_of_interest <- Biobase::pData(y_combined)$sample_barcode
  df_combined$meanX <- NA
  df_combined$meanY <- NA
  df_combined$size <- 2
  df_combined$segment.size <- 0.25
  df_combined$fontface <- 'plain'

  # Calculate center coordinates for each cancer type and assign labels
  unique_cancer_types <- setdiff(unique(df_combined$color),"unknown")
  for (j in 1:length(unique_cancer_types)){
    meanX <- median(df_combined[which(df_combined$color==unique_cancer_types[j]),][["Dimension 1"]])
    meanY <- median(df_combined[which(df_combined$color==unique_cancer_types[j]),][["Dimension 2"]])
    idx <- which.min(abs(df_combined[which(df_combined$color==unique_cancer_types[j]),][["Dimension 1"]]-meanX) +
                       abs(df_combined[which(df_combined$color==unique_cancer_types[j]),][["Dimension 2"]]-meanY))
    df_combined[which(df_combined$color==unique_cancer_types[j]),]$meanX[idx] <- meanX
    df_combined[which(df_combined$color==unique_cancer_types[j]),]$meanY[idx] <- meanY

    df_combined[which(df_combined$color==unique_cancer_types[j]),]$sample_of_interest[idx] <- unique_cancer_types[j]
    df_combined[which(df_combined$color==unique_cancer_types[j]),]$size[idx] <- 3
    df_combined[which(df_combined$color==unique_cancer_types[j]),]$segment.size[idx] <- NA
    df_combined[which(df_combined$color==unique_cancer_types[j]),]$fontface[idx] <- 'bold'
  }

  saveRDS(df_combined,paste0(outdir,"/",plot_title,"_","tsne_df.rda"))

  pdf(file=paste0(outdir,"/",plot_title,"_","tsne.pdf"),width=6,height=6)
  g <- ggplot(as.data.frame(df_combined),aes(x=`Dimension 1`,y=`Dimension 2`,group=color)) +
    geom_point(aes(colour=color),size=0.25) + theme_bw() + theme(legend.position="none")
  g <- g + ggrepel::geom_text_repel(aes(label=sample_of_interest),
                                    min.segment.length = unit(0.15, 'lines'),
                                    size=df_combined$size,segment.size=df_combined$segment.size,
                                    fontface=df_combined$fontface,
                                    box.padding = 0,
                                    point.padding=NA,
                                    force=10)
  print(g)
  dev.off()
}

# Main functions ---------------------------------------------------------------

#' @export
cv <- function(annot=NULL,params=NULL,refset,max_k=20){

  if (is.null(annot) && is.null(params)){
    stop("params and annot cannot both be empty")
  }
  if (is.null(params)){
    params <- setdiff(colnames(annot),"sample_barcode")
  }

  refset.custom <- refset
  if (!is.null(annot)){
    tmp <- merge(pData(refset.custom),annot,by.x="sample_barcode",by.y="row.names",all.x=T)
    rownames(tmp) <- tmp$sample_barcode
    #stopifnot(identical(rownames(tmp[rownames(pData(refset.custom)),]),rownames(pData(refset.custom))))
    pData(refset.custom) <- tmp[rownames(pData(refset.custom)),]
  }

  r.tcga <- .calculate_correlations(inhouse = refset, refset = refset)

  cv.res.k <- list()
  for (k in 1:max_k){
    cv.res.k[[k]] <- cv.eval.k(r.tcga,refset.custom,params,k=k)
  }

  cv.res.k.acc <- data.frame(1:length(cv.res.k),
                             unlist(lapply(cv.res.k,function(x){
                               mean(as.data.frame(x$cv.df)[["Accuracy"]])})))
  colnames(cv.res.k.acc) <- c("k","Accuracy")

  cv.res.k.acc <- cv.res.k.acc[order(cv.res.k.acc$Accuracy,decreasing = T),]
  k.best <- min(cv.res.k.acc[cv.res.k.acc$Accuracy==max(cv.res.k.acc$Accuracy),]$k)
  g.k_vs_accuracy <- ggplot(cv.res.k.acc,aes(x=k,y=Accuracy)) + geom_point(stat="identity") +
    theme_classic() + geom_vline(xintercept = k.best,linetype = "longdash",colour="gray")

  cv.res.k.best.acc <- cv.res.k[[k.best]]$cv.df.melt[cv.res.k[[k.best]]$cv.df.melt$Var2 %in%
                                                       c("Accuracy","AccuracyPValue"),]
  colnames(cv.res.k.best.acc) <- c("category","measure","value")
  g.accuracy_best_k <- ggplot(cv.res.k.best.acc,aes(x=category,y=value)) + geom_bar(stat="identity") +
    theme_classic() + facet_wrap(~ measure) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  return(list(k.best = k.best,
              cv.res.k.acc = cv.res.k.acc,
              cv.res.k.best.acc = cv.res.k.best.acc,
              g.k_vs_accuracy = g.k_vs_accuracy,
              g.accuracy_best_k = g.accuracy_best_k,
              refset.custom = refset.custom))
}

#' @export
run_tsne <- function(inhouse, refset, outdir, include_all=F, niter=1000, perplexity=30, color_by = NULL){

  dir.create(outdir,showWarnings = F,recursive = T)
  if (is.null(color_by)){
    colr <- Biobase::pData(refset)$cancer_type
  } else if (any(colnames(Biobase::pData(refset)) == color_by)){
    colr <- Biobase::pData(refset)[[color_by]]
  } else {
    stop("color_by ",color_by, " is not a column in the provided reference dataset.")
  }

  sample_info <- data.frame(sample_barcode=c(rep(NA,ncol(refset)),colnames(inhouse)),
             cancer_type=c(Biobase::pData(refset)$cancer_type,rep("unknown",ncol(inhouse))),
             color=c(colr,rep("unknown",ncol(inhouse))),stringsAsFactors = F)
  rownames(sample_info) <- c(colnames(refset),colnames(inhouse))
  stopifnot(identical(rownames(tail(sample_info,n=ncol(inhouse))),
                      tail(sample_info,n=ncol(inhouse))$sample_barcode))

  stopifnot(identical(rownames(Biobase::fData(refset)),rownames(Biobase::fData(inhouse))))

  pdat <- Biobase::AnnotatedDataFrame(data=sample_info)
  fdat <- Biobase::AnnotatedDataFrame(data=Biobase::fData(refset))
  y_combined <- cbind(Biobase::exprs(refset),Biobase::exprs(inhouse))
  y_combined_transformed <- t(log2(t(y_combined)+1))
  y_combined <- Biobase::ExpressionSet(assayData = y_combined_transformed, phenoData = pdat,
                                       featureData = fdat)

  if (include_all){
    plot_title <- "all"
    .tsne_core(y_combined = y_combined, outdir = outdir, plot_title = "all", niter = niter,
              perplexity = perplexity)
  } else {
    idxnr_unknown <- which(Biobase::pData(y_combined)$cancer_type=="unknown")
    stopifnot(length(idxnr_unknown)==ncol(inhouse))

    idxnr_refset <- 1:(ncol(y_combined)-length(idxnr_unknown))
    stopifnot(max(idxnr_refset)==ncol(refset))

    for (idxnr in idxnr_unknown){
      y_combined.subset <- y_combined[,c(idxnr_refset,idxnr)]
      .tsne_core(y_combined = y_combined.subset, outdir = outdir,
                plot_title = rownames(Biobase::pData(y_combined.subset))[ncol(y_combined.subset)],
                niter = niter, perplexity = perplexity)
    }
  }
}

#' @export
classify <- function(input_test, input_train = "tcga", outdir = NULL, k = 6,
                     subset_to=c(), genome = NULL, labels_train = NULL,
                     extension_test = ".gene_counts", extension_train = ".gene_counts",
                     override_labels = F, is_matrix = F){

  # Change: labels_train neeeds to be a data.frame, with rownames as sample barcodes and
  # a column titled cancer_type

  inhouse <- .setup_expression(input_test, extension = extension_test, is_matrix)
  refset <- .setup_expression(input_train, extension = extension_train, is_matrix)

  if (!is.null(labels_train) && ((! "cancer_type" %in% colnames(pData(refset))) || (override_labels==T))){
    #pData(refset)$cancer_type <- labels_train
    tmp <- pData(refset)
    tmp <- merge(tmp, labels_train, by="row.names",all.x=T,all.y=F)
    rownames(tmp) <- tmp$Row.names
    tmp$Row.names <- NULL

    stopifnot(all(rownames(tmp) %in% rownames(pData(refset))))
    tmp <- tmp[rownames(pData(refset)),]
    pData(refset) <- tmp
    refset <- refset[,!is.na(pData(refset)$cancer_type)]
  }

  if (all(grepl("^ENSG",rownames(refset))) && ! all(grepl("^ENSG",rownames(inhouse)))){
    inhouse <- .translate(inhouse, refset)
  }

  res <- sync_expression_sets(eset_1=refset, eset_2=inhouse)
  refset <- res$eset_1
  inhouse <- res$eset_2

  if (!(identical(rownames(inhouse),rownames(gene_info)))){
    if (genome %in% c("hg19","hg38")){
      g_info <- .get_gene_info(genes=rownames(inhouse), genome=genome)
    } else {
      stop("Input gene names do not fully overlap provded default reference dataset and no value for the genome parameter was provided")
    }
  } else {
    g_info <- gene_info
  }

  inhouse <- inhouse[rownames(g_info),]
  refset <- refset[rownames(g_info),]
  stopifnot(identical(rownames(inhouse),rownames(refset)))

  inhouse <- .rpkm_eset(inhouse, gene_info = g_info)
  refset <- .rpkm_eset(refset, gene_info = g_info)

  if (length(subset_to)>0){
    refset <- .subset_samples(refset = refset, subset_to = subset_to)
  }

  r <- .calculate_correlations(inhouse = inhouse, refset = refset)
  knn.res <- .knn(r=r, refset = refset, k = k)

  if (!is.null(outdir)){
    dir.create(outdir,showWarnings = F,recursive = T)
    .write_correlation_table(r = r, inhouse = inhouse, refset = refset, outdir = outdir)
    write.table(knn.res,file=paste0(outdir,'kNN_prediction.txt'),quote=F,row.names=F,col.names = F)
    saveRDS(knn.res,file=paste0(outdir,'kNN_res.rda'))
    saveRDS(r,file=paste0(outdir,'corr_mat.rda'))
  }

  return(list(knn.res = knn.res,
              r = r,
              inhouse = inhouse,
              refset = refset,
              outdir = outdir))
}
