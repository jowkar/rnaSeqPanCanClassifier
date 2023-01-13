.load_tcga <- function(){
  gencode <- data.table::fread("/home/joakim/proj/pdac/Pipelines/rna/classification/gencode19/R/gencode.txt",data.table = F,stringsAsFactors = F)
  gencode$geneId <- stringr::str_split_fixed(gencode$geneId,"\\.",2)[,1]
  classes.code <- data.table::fread("/home/joakim/proj/pdac/Pipelines/rna/classification/gencode19/R/classes.code.txt",data.table = F,stringsAsFactors = F)
  classes.code$code <- stringr::str_split_fixed(classes.code$code,"\\.",2)[,1]

  gene_info <- data.table::fread("/home/joakim/proj/pdac/Pipelines/rna/classification/nilsson_olofsson_nneighbor/mat/R/expr__info_1.txt",data.table = F,stringsAsFactors = F)
  colnames(gene_info) <- c("gene_id","gene_symbol","length")
  gene_info$gene_id <- stringr::str_split_fixed(gene_info$gene_id,"\\.",2)[,1]
  rownames(gene_info) <- gene_info$gene_id
  sample_info <- data.table::fread("/home/joakim/proj/pdac/Pipelines/rna/classification/nilsson_olofsson_nneighbor/mat/R/expr__info_2.txt",data.table = F,stringsAsFactors = F)
  sample_info <- sample_info[,c(1,3)]
  colnames(sample_info) <- c("cancer_type","sample_barcode")
  rownames(sample_info) <- sample_info$sample_barcode

  expr_pancan <- data.table::fread("/home/joakim/proj/pdac/Pipelines/rna/classification/nilsson_olofsson_nneighbor/mat/R/expr__data_norm.csv",data.table = F,stringsAsFactors = F)
  rownames(expr_pancan) <- gene_info$gene_id
  colnames(expr_pancan) <- sample_info$sample_barcode

  # create eSet
  stopifnot(identical(colnames(expr_pancan),rownames(sample_info)))
  stopifnot(identical(rownames(expr_pancan),rownames(gene_info)))

  gene_info$is_coding <- F
  gene_info$is_coding[gene_info$gene_id %in% classes.code$code] <- T

  ens_sex <- gencode$geneId[which(gencode$chrom %in% c("chrX","chrY"))]
  gene_info$sex_chromosome <- F
  gene_info$sex_chromosome[gene_info$gene_id %in% ens_sex] <- T

  pdat <- Biobase::AnnotatedDataFrame(data=sample_info)
  fdat <- Biobase::AnnotatedDataFrame(data=gene_info)
  tcga <- Biobase::ExpressionSet(assayData=as.matrix(expr_pancan),phenoData=pdat,featureData=fdat)

  # Subset to autosomes
  tcga <- tcga[which(! Biobase::fData(tcga)$sex_chromosome),]

  # Subset to coding
  tcga <- tcga[which(Biobase::fData(tcga)$is_coding),]

  # gencode <- gencode[which(gencode$geneId %in% rownames(tcga)),]
  # gencode$is_coding <- F
  # gencode$is_coding[gencode$geneId %in% classes.code$code] <- T

  return(tcga)
}

.load_ccle <- function(){
  library(stringr)
  tab <- data.table::fread("~/proj/ganesh/Investigations/data/CCLE/CCLE_DepMap_18Q2_RNAseq_RPKM_20180502.gct",sep="\t",
                           data.table = F,stringsAsFactors = F)
  tab$Name <- stringr::str_split_fixed(tab$Name,"\\.",2)[,1]
  rownames(tab) <- tab$Name
  tab$Name <- NULL

  gene_info <- data.frame(ens_id=rownames(tab),genes_symb=tab$Description,stringsAsFactors = F)
  tab$Description <- NULL
  rownames(gene_info) <- gene_info$ens_id

  sample_info <- data.frame(cancer_type=stringr::str_split_fixed(colnames(tab),"_",2)[,2],
             sample_barcode=colnames(tab),
             stringsAsFactors = F)
  rownames(sample_info) <- sample_info$sample_barcode

  pdat <- Biobase::AnnotatedDataFrame(data=sample_info)
  fdat <- Biobase::AnnotatedDataFrame(data=gene_info)
  stopifnot(identical(rownames(tab),rownames(fdat)))
  stopifnot(identical(colnames(tab),rownames(pdat)))

  stopifnot(min(tab,na.rm=T)==0)
  tab[is.na(tab)] <- 0

  ccle <- Biobase::ExpressionSet(assayData=as.matrix(tab),phenoData=pdat,featureData=fdat)

  return(ccle)
}

#' @export
sync_expression_sets <- function(eset_1,eset_2,eset_3=NULL){
  if (!is.null(eset_3)){
    isect <- intersect(intersect(rownames(eset_1),rownames(eset_2)),rownames(eset_3))
  } else {
    isect <- intersect(rownames(eset_1),rownames(eset_2))
  }
  if (length(isect)==0){
    stop("No overlapping features between test and reference set")
  }
  eset_1 <- eset_1[isect,]
  eset_2 <- eset_2[isect,]

  if (!is.null(eset_3)){
    eset_3 <- eset_3[isect,]
    stopifnot(identical(rownames(eset_1),rownames(eset_3)))
    stopifnot(identical(rownames(eset_2),rownames(eset_3)))
    out <- list(eset_1=eset_1,eset_2=eset_2,eset_3=eset_3)
  } else {
    stopifnot(identical(rownames(eset_1),rownames(eset_2)))
    out <- list(eset_1=eset_1,eset_2=eset_2)
  }

  return(out)
}

.load_gene_info <- function(){
  library(biomaRt)
  ensembl_hg19 <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                     host = "grch37.ensembl.org", path = "/biomart/martservice",
                     dataset = "hsapiens_gene_ensembl")
  map.ensembl_symbol_hg19 <- getBM(attributes = c("ensembl_gene_id",
                                                  "external_gene_name",
                                                  "transcript_length"),
                                   mart = ensembl_hg19)

  ensembl_hg38 <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  map.ensembl_symbol_hg38 <- getBM(attributes = c("ensembl_gene_id",
                                                  "external_gene_name",
                                                  "transcript_length"),
                                   mart = ensembl_hg38)

  return(list(map.ensembl_symbol_hg19 = map.ensembl_symbol_hg19,
              map.ensembl_symbol_hg38 = map.ensembl_symbol_hg38))
}

.savedata <- function(){
  tcga <- .load_tcga()
  ccle <- .load_ccle()

  res <- sync_expression_sets(eset_1 = tcga, eset_2 = ccle)
  tcga <- res$eset_1
  ccle <- res$eset_2

  gene_info <- data.frame(length=Biobase::fData(tcga)$length)
  rownames(gene_info) <- rownames(Biobase::fData(tcga))

  tmp <- .load_gene_info()
  map.ensembl_symbol_hg19 <- tmp$map.ensembl_symbol_hg19
  map.ensembl_symbol_hg38 <- tmp$map.ensembl_symbol_hg38

  usethis::use_data(tcga, ccle, gene_info,
                    map.ensembl_symbol_hg19, map.ensembl_symbol_hg38,
                    internal = F, overwrite = T,version = 3)
}
