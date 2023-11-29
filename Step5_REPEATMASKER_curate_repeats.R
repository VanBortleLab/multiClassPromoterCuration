rm(list = ls())
require(data.table)
require(dplyr)


a=fread("REPEATMASKER.tsv",sep="\t",header = T)
# b=a[sort(sample(c(1:nrow(a)),as.integer(nrow(a)/25))),]
# fwrite(b,"REPEATMASKERR.tsv",sep="\t",col.names = T,row.names = F,quote = F,append=F)
###!!!!!!! FILTER ON THE BASIS OF SW SCORE??????
###!!!!!!! GET RID of DNA?, RC? and SINE? ,,tho they are very few

#converting into the format of masterfiles_final
a=a %>% transmute(Chr=genoName,Start=genoStart,Stop=genoEnd,TSS="R",TSS_median=ifelse(strand=="+",genoStart,genoEnd),URSID="R",Strand=strand,RNA_type=paste0(repClass,"_r"),TOTAL_INFO=repName)
a= a %>% arrange(Chr,Start) %>% filter(!grepl("_",Chr))

table(a$RNA_type)
table(a$Chr)

write.table(a,"REPEATMASKER_processed.bed",sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE,append=F)

