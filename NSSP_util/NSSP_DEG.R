

# DEG for NSSP

library(TCC)
fn <- "NSSP_Data/smRNAs.xls"
data <- read.table(fn, header = TRUE, row.names = 1)
data_rc <- data[,startsWith(colnames(data), "ReadCnt")]
colnames(data_rc) <- gsub('ReadCnt.','',colnames(data_rc))

deg_name <- 'DEG003'
sel_col <- c(2,3)
group <- c(1,2)
indata_rc <- data_rc[,sel_col]
run(indata_rc, group)

run <- function(indata_rc, group) {
  
  tcc <- new("TCC", indata_rc, group)
  tcc <- filterLowCountGenes(tcc = tcc, low.count = 100)
  tcc <- calcNormFactors(tcc, norm.method="tmm", test.method="edger", iteration=3, FDR=0.1, floorPDEG=0.05)
  tcc <- estimateDE(tcc, test.method="edger", FDR=0.1)
  
  indata_rc_norm <- getNormalizedData(tcc)
  gene_id <- row.names(indata_rc_norm)
  indata_rc_norm2 <- cbind(gene_id, indata_rc_norm)
  row.names(indata_rc_norm2) <- NULL
  indata_rc <- indata_rc[c(row.names(indata_rc_norm)),]
  gene_id <- row.names(indata_rc)
  indata_rc2 <- cbind(gene_id, indata_rc)
  row.names(indata_rc2) <- NULL
  result <- getResult(tcc, sort=TRUE)
  pre_data_result <- merge(as.data.frame(result), as.data.frame(indata_rc2), by=c("gene_id"))
  data_result <- merge(pre_data_result, as.data.frame(indata_rc_norm2), by=c("gene_id"))
  colnames(data_result) <- gsub(".x", ".RawReadCnt", colnames(data_result))
  colnames(data_result) <- gsub(".y", ".NormReadCnt", colnames(data_result))
  
  write.table(data_result, paste("NSSP_Data/smRNAs",deg_name,"xls",sep = '.'), sep="\t", quote=F, col.names=T, row.names=T)
  png(paste("NSSP_Data/smRNAs",deg_name,"png",sep = '.'), width=640, height = 640)
  plot(tcc, median.lines = TRUE, cex=0.4)
  dev.off()
  
}

tcc <- new("TCC", indata_rc, group)
dim(tcc$count)
tcc <- filterLowCountGenes(tcc = tcc, low.count = 100)
dim(tcc$count)

tcc <- calcNormFactors(tcc, norm.method="tmm", test.method="edger", iteration=3, FDR=0.1, floorPDEG=0.05)
tcc <- estimateDE(tcc, test.method="edger", FDR=0.1)

indata_rc_norm <- getNormalizedData(tcc)
gene_id <- row.names(indata_rc_norm)
indata_rc_norm2 <- cbind(gene_id, indata_rc_norm)
row.names(indata_rc_norm2) <- NULL

indata_rc <- indata_rc[c(row.names(indata_rc_norm)),]
gene_id <- row.names(indata_rc)
indata_rc2 <- cbind(gene_id, indata_rc)
row.names(indata_rc2) <- NULL

result <- getResult(tcc, sort=TRUE)


pre_data_result <- merge(as.data.frame(result), as.data.frame(indata_rc2), by=c("gene_id"))
data_result <- merge(pre_data_result, as.data.frame(indata_rc_norm2), by=c("gene_id"))
colnames(data_result) <- gsub(".x", ".RawReadCnt", colnames(data_result))
colnames(data_result) <- gsub(".y", ".NormReadCnt", colnames(data_result))


write.table(data_result, paste("NSSP_Data/smRNAs",deg_name,"xls",sep = '.'), sep="\t", quote=F, col.names=T, row.names=T)
png(paste("NSSP_Data/smRNAs",deg_name,"png",sep = '.'), width=640, height = 640)
plot(tcc, median.lines = TRUE, cex=0.4)
dev.off()
