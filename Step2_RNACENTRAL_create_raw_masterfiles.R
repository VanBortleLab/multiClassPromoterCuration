require(here)
setwd(here::here())
#Caution!!!! 
##Make_sure, you replace all "Ensembl/Gencode" with "EnsemblGencode" in Grch38.bed file since the "/" sign had issues naming variables later
##Make sure you put important columns of each database in IMP_COLUMNS_for_each_DB2.R 
##Make sure you have "homo_sapiens.GRCh38.bed" in same directory ./
##Make sure you have all databases in ./databases/


#################################################################################
###############################PART A###############################################
############ANALYZING GRCH38 and IMP_COLUMNS.R#####################
#################################################################################
#################################################################################
require(data.table)
require(dplyr)
####Analyazing GRCH bed file
{
G38= read.table(file="./homo_sapiens.GRCh38.bed", sep="\t",quote="",comment.char="")
G38=G38 %>% arrange(V1,V2)
# set.seed(12)
# G38=G38[runif(nrow(G38)*0.18,1,nrow(G38)),]
# write.table(G38,"homo_sapiens.GRCh38_trimmed.bed",sep="\t",row.names = F,col.names = F,append = F,quote = F)
G38=subset(G38,select=c(V1,V2,V3,V4,V6,V14,V15))
colnames(G38)=c("Chr","Start","Stop","URSID","Strand","RNA_type","Databases")
###########Finding types of RNA mentioned in the GRCh38 file.
TypesOfRNA=unique(G38$RNA_type)
#########Dividing G38 into multiple files for each RNA type
for(i in TypesOfRNA)
{
  do.call('<-',list(paste0("G38_",i),subset(G38,RNA_type==i)))
  # print(paste0("G38_",i," is created "))
}
########################Finding databases referred by each of these newly created G38_RNAtype files
for(i in unique(G38$RNA_type))
{
  do.call('<-',list(paste0(i,"_db"),unique(unlist(strsplit(unique(subset(G38,RNA_type==i)$Databases),',')))))
  # print(paste0(i,"_db is created"))
}
#list of total databases used in G38
all_db=unique(unlist(strsplit(unique(G38$Databases),',')))

##Checking if we have all the required databases
print("Please have the following databases in folder ./database and have information column information of these databases in IMP_COLUMNS_for_each_DB2.R");print(paste0(all_db,".tsv"))
}

##OPening IMP_COLUMNS.R
eval(parse("./Step1_RNACENTRAL_column_information.R"))

###Checking for necessary input files
{missing=all_db[!file.exists(paste0("./database/",all_db,".tsv"))]
if(length(missing)>0){print(paste0("ERROR !! ./database/",missing,".tsv is missing"))}
for(i in all_db)
{ if(exists(paste0("imp_",i))==F) print(paste0("Variable imp_",i," is missing"))}
}##luka##

########################Referring to the "Step1_RNACENTRAL_column_information.R" and finding those databases that has useful details in it.
{filter=c()  ###list that stores index of unimportant databases
for (i in seq_along(all_db))
{if(length(eval(parse(text=paste0("imp_",all_db[i]))))==0) {filter=c(filter,i)}}
imp_db=all_db[-c(filter)]  #removing unimportant databases
print(paste0("These are important databases that have useful details about RNA->",paste0(imp_db,collapse="       ")))
unimp_db=setdiff(all_db,imp_db) 
print(paste0("These are unimportant databases we wont use here ->",paste0(unimp_db,collapse="    ")))
for(i in unique(G38$RNA_type))
{do.call("<-",list(paste0(i,"_db_imp"),setdiff(eval(parse(text=paste0(i,"_db"))),unimp_db)))}
}##luka##

#################################################################################
#################################PART B##########################################
############LOADING AND PREPROCESSING THE IMPORTANT DATABASES####################
#################################################################################
#################################################################################
##Loading those important/useful databases and removing unimportant columns from them
for(i in imp_db)
{ 
  print(paste0("Currently loading database is ",i))
  do.call('<-',list(paste0("impdb_",i),fread(file=paste0("./database/",tolower(i),".tsv"),sep='\t',quote="")))
  print(" Now, deleting unimportant columns")
  do.call('<-',list(paste0("impdb_trim_",i),eval(parse(text=paste0("impdb_",i,"[,c(imp_",i,"),with=FALSE]")))))
  eee=eval(parse(text=paste0("impdb_trim_",i)))
  colnames(eee)=c("URSID","RNA_type","Info")
  do.call("<-",list(paste0("impdb_trim_",i),eee))
}##luka##

#######################Removing rows from databases(eg. impdb_trim_ENA) that doesn't have valuable information and collapsing rows with same URSID using semicolon ";"
for(i in imp_db)
{ 
  print(paste0("Trimming ",i))
  kkk=eval(parse(text=paste0("impdb_trim_",i)))
  ttt=subset(kkk,kkk$Info!="")
  tttt=ttt %>% group_by(URSID) %>% summarize(RNA_type=RNA_type[1],Info=paste(unique(Info),collapse=";"))
  print(paste0("no of rows reduced = ",nrow(kkk)-nrow(tttt)))
 do.call("<-",list(paste0("impdb_trim_",i),tttt))
}##luka##

############Storing_processed_databases(i.e. FIRST: rows with empty info deleted, SECOND: rows with same URSID collapsed)
system("rm -r  databases_processed")
system("mkdir databases_processed")
print("Storing processed databases in ./databases_processed")
for(i in imp_db)
{write.table(eval(parse(text=paste0("impdb_trim_",i))),file=paste0("./databases_processed/",i),row.names = FALSE,col.names = TRUE,quote=FALSE,append=FALSE,sep="\t")}


#################################################################################
#################################PART C##########################################
#####################MAKING raw MASTER FILES ######################
#################################################################################
#################################################################################

#######MAKING MAKEMASTER function
MAKEMASTER= function(type,master)
{   
  total=nrow(master)
  for(i in eval(parse(text=paste0(type,"_db_imp"))))
  { ddb=as.data.frame(eval(parse(text=paste0("impdb_trim_",i,"_",type))))
  a1a=ddb[,c("URSID","Info")]
  colnames(a1a)=c("URSID",i)
  if(nrow(a1a)!=0){a1a[,1]=paste0(a1a[,1],"_9606")}
  master=left_join(master,a1a,by=c("URSID"))
  }
  return(master)
}##luka##

#####################Calling MAKEMASTER function for each RNA type
system("rm -r rnacentral_masterfiles_raw")
system("mkdir rnacentral_masterfiles_raw")   ####making a new folder to keep our master file
for(i in TypesOfRNA)
{ 
  print(paste("making masterfile for",i))
  #trimming database to contain current RNA type only
  for (j in eval(parse(text=paste0(i,"_db_imp"))))
  {do.call("<-",list(paste0("impdb_trim_",j,"_",i),subset(eval(parse(text=paste0("impdb_trim_",j))),RNA_type==i)))}
  #Calling MAKEMASTER function
  finally=MAKEMASTER(i,eval(parse(text=paste0("G38_",i))))
  finally=finally[,sapply(c(1:ncol(finally)),function(x) sum(!is.na(finally[,x])))!=0] #removing columns that have NA value in all rows
  write.table(finally,file=paste0("./rnacentral_masterfiles_raw/",i,"_masterfile.bed"),quote=F,row.names=F,sep="\t")
  #now removing RNA_type specific databases
  for (j in eval(parse(text=paste0(i,"_db_imp"))))
  {remove(list=c(paste0("impdb_trim_",j,"_",i)))}
}##luka##


