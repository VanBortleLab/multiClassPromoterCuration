rm(list = ls())
require(dplyr)
require(data.table)


  #RNACENTRAL
  rnacentral=fread("./rnacentral_masterfiles_final/TOTAL_rnacentral.bed",sep="\t")
  rnacentral=rnacentral[!grepl("_",rnacentral$Chr),] ; rnacentral=rnacentral[!grepl("\\.",rnacentral$Chr),] ## removing bad chromosome with "_" and "."
  
  #REPEATS
  repets=fread("REPEATMASKER_processed.bed")

  ##PCG
  pcg=fread("./ENSEMBL_processed.bed",sep="\t")
  
  
  #COMBINING ALL
  
  FINAL=rbind(rnacentral,repets,pcg)  
  FINAL = as.data.frame(FINAL) %>% rowwise() %>%
    mutate(V1=Chr,V2=ifelse(Strand=="+",TSS_median-50,TSS_median-150),
           V3=ifelse(Strand=="+",TSS_median+150,TSS_median+50)) %>% 
    select(V1,V2,V3,everything())
  FINAL= FINAL %>% rowwise() %>% mutate(V2=max(0,V2))  ##Making negative V2 as 0
  FINAL$V2=as.integer(FINAL$V2);FINAL$V3=as.integer(FINAL$V3) ##Making V2 and V3 as integer
  FINAL=FINAL %>% dplyr::rename(TSS_Start=Start,TSS_Stop=Stop) %>% dplyr::arrange(V1,V2)
  
  fwrite(FINAL,"MultiClassPromoter.tsv",row.names = F,col.names = T,sep="\t",quote=F,append=F)
  remove(rnacentral,repets,pcg)

  
  

