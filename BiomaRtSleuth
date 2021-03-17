library("sleuth")

#Esta función permite mapear, a partir de la base de datos de
tx2gene <- function(){
  
  #     Dataset you want to use. To see the different datasets available within a biomaRt yo$
  #     host
  
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = "celegans_gene_ensembl")
  
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                       "external_gene_name"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = "ensembl_transcript_id",
                       ens_gene = "ensembl_gene_id", ext_gene = "external_gene_name")
  return(t2g)
}

t2g <- tx2gene()



base_dir<-"~/Dropbox/"

#All samples
#samples <- paste0("sample_", c("14BE01A_R","14BE02A_R","14BE08A_R","17BE01A_R",
#                                "17BE02A_R","17BE03A_R"))

#Selected samples

samples <- paste0("sample_LG00", c("1_S13_R","2_S14_R","3_S15_R",
                               "4_S16_R","6_S18_R"))
                               

kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))

#All samples
#s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("2014","2014","2014", "2017",
#
#"2017","2017"), stringsAsFactors=FALSE)

#Selected samples

s2c <- data.frame(path=kal_dirs, sample=samples, muestras = c("N2","gcy-25","N2", "gcy-25",
                                                              "gcy-25"), stringsAsFactors=FALSE)

so <- sleuth_prep(s2c, ~muestras, target_mapping = t2g,extra_bootstrap_summary = TRUE)

so <- sleuth_fit(so)
so <- sleuth_wt(so, which_beta="muestrasN2")                   
sleuth_live(so)

setwd("~/Dropbox/")
resultados<-read.table("test_table.csv",sep=",",
                       header=TRUE)
significativos<-which(resultados$qval<0.1)
significativos<-resultados[significativos,]
upregulated<-which(significativos$b>0)
upregulated<-significativos[upregulated,]
downregulated<-which(significativos$b<0)
downregulated<-significativos[downregulated,]

write.table(upregulated,file="~/Dropbox/Upregulated_N2vsCyg-25.txt",sep="\t")
write.table(downregulated,file="~/Dropbox/Downregulated_N2vsCyg-25.txt",sep="\t")

upregulated$ext_gene


