setwd("/Users/rkc/Downloads/FLOW/MAKEMASTER_new")

###############################PART A###########################################
############PREPROCESSING GRCH38 and creating many variables#####################
#################################################################################
#Caution!!!!
##Make_sure, you replace all "Ensembl/Gencode" with "EnsemblGencode" in Grch38.bed file since the "/" sign had issues naming variables later
##Make sure you put important columns of each database in IMP_COLUMNS_for_each_DB2.R 
##Make sure you have "homo_sapiens.GRCh38.bed" in same directory ./
##Make sure you have all databases in ./databases/
require(data.table)
require(dplyr)

eval(parse("./IMP_COLUMNS_for_each_DB2.R"))

###########Loading and keeping only the important columns of_G3rch8_FILE
G38= read.table(file="./homo_sapiens.GRCh38.bed", sep="\t",quote="",comment.char="")
G38=subset(G38,select=c(V1,V2,V3,V4,V6,V14,V15))
colnames(G38)=c("Chr","Start","Stop","URSID","Strand","RNA_type","Databases")
###########Finding types of RNA mentioned in the GRCh38 file.
TypesOfRNA=unique(G38$RNA_type)
#########Dividing G38 into multiple files for each RNA type
for(i in TypesOfRNA)
{
  do.call('<-',list(paste0("G38_",i),subset(G38,RNA_type==i)))
  print(paste0("G38_",i," is created "))
}
########################Finding databases referred by each of these newly created G38_RNAtype files
for(i in unique(G38$RNA_type))
{
  do.call('<-',list(paste0(i,"_db"),unique(unlist(strsplit(unique(subset(G38,RNA_type==i)$Databases),',')))))
  print(paste0(i,"_db is created"))
}
#list of total databases used in G38
all_db=unique(unlist(strsplit(unique(G38$Databases),',')))     
print("all_db is created. NOW. make sure you have all the databases in all_db list with you in database folder and in IMP_COLUMNS_for_each_DB2.R file")




########################Referring to the "IMP_COLUMNS_for_each_DB2.R" and storing only those databases that has useful details in it.
filter=c()  ###list that stores index of unimportant databases
for (i in seq_along(all_db))
{
  if(length(eval(parse(text=paste0("imp_",all_db[i]))))==0)
  {filter=c(filter,i)}
}
imp_db=all_db[-c(filter)]  #removing unimportant databases
print(paste0("These are important databases that have useful extra details->",paste0(imp_db,collapse="       ")))
unimp_db=setdiff(all_db,imp_db) 
print(paste0("These are unimportant databases we wont use here ->",paste0(unimp_db,collapse="    ")))
for(i in unique(G38$RNA_type))
{
  do.call("<-",list(paste0(i,"_db_imp"),setdiff(eval(parse(text=paste0(i,"_db"))),unimp_db)))
  print(paste0(i,"_db_imp is created"))
}


#################################PART B##########################################
############LOADING AND PREPROCESSING THE IMPORTANT DATABASES####################
#################################################################################
######################Loading those important/useful databases and removing unimportant columns from them
##### I had tough time removing unimportant columns because of technical difficulties, finally I figured out,I was missing "with=FALSE"
for(i in imp_db)
{ 
  print(paste0("Currently loading database is ",i))
  do.call('<-',list(paste0("impdb_",i),fread(file=paste0("./database/",tolower(i),".tsv"),sep='\t',quote="")))
  print(" Now, deleting unimportant columns")
  do.call('<-',list(paste0("impdb_trim_",i),eval(parse(text=paste0("impdb_",i,"[,c(imp_",i,"),with=FALSE]")))))
  eee=eval(parse(text=paste0("impdb_trim_",i)))
  colnames(eee)=c("URSID","RNA_type","Info")
  do.call("<-",list(paste0("impdb_trim_",i),eee))
}  
print("Important databases with important columns are saved as -> impdb_trim_tRNA impdb_trim_miRNA .........")

#######################Removing rows from databases(eg. impdb_trim_ENA) that doesn't have valuable information and collapsing using dplyr
for(i in imp_db)
{ 
  print(paste0("Trimming ",i))
  kkk=eval(parse(text=paste0("impdb_trim_",i)))
  ttt=subset(kkk,kkk$Info!="")
  tttt=ttt %>% group_by(URSID) %>% summarize(RNA_type=RNA_type[1],Info=paste(unique(Info),collapse=";"))
  print(paste0("no of rows reduced = ",nrow(kkk)-nrow(tttt)))
 do.call("<-",list(paste0("impdb_trim_",i),tttt))
}


############Storing_processed_databases(first rows with empty info deleted, second rows with same URSID collapsed)
system("mkdir databases_processed")
print("Storing processed databases in ./databases_processed")
for(i in imp_db)
{print(i)
  write.table(eval(parse(text=paste0("impdb_trim_",i))),file=paste0("./databases_processed/",i),row.names = FALSE,col.names = TRUE,quote=FALSE,append=FALSE,sep="\t")
}


#################################################################################
#################################PART C##########################################
#####################MAKING AND CALLING MAKEMASTER FUNCTION######################
#################################################################################
#################################################################################

#######MAKING MAKEMASTER function
MAKEMASTER= function(type,master)
{   
  #master[eval(parse(text=paste0(type,"_db_imp")))]=NA  ###Adding columns for databases that will be used in making this master file
  total=nrow(master)
  for(i in eval(parse(text=paste0(type,"_db_imp"))))
  { ddb=as.data.frame(eval(parse(text=paste0("impdb_trim_",i,"_",type))))
  a1a=ddb[,c("URSID","Info")]
  colnames(a1a)=c("URSID",i)
  if(nrow(a1a)!=0){a1a[,1]=paste0(a1a[,1],"_9606")}
  master=left_join(master,a1a,by=c("URSID"))
  }
  return(master)
}

#####################Calling MAKEMASTER functin for each RNA type
system("rm -r masterfiles_new")
system("mkdir masterfiles_new")   ####making a new folder to keep our master file
for(i in TypesOfRNA)
{ 
  #if(i=="piRNA" | i=="lncRNA") next
  print(paste("making masterfile for",i))
  #trimming database to contain current RNA type only
  for (j in eval(parse(text=paste0(i,"_db_imp"))))
  {do.call("<-",list(paste0("impdb_trim_",j,"_",i),subset(eval(parse(text=paste0("impdb_trim_",j))),RNA_type==i)))}
  #Calling Make_master function
  do.call("<-",list(paste0("masterlist_",i),MAKEMASTER(i,eval(parse(text=paste0("G38_",i))))))
  write.table(eval(parse(text=paste0("masterlist_",i))),file=paste0("./masterfiles_new/",i,"_masterfile.bed"),quote=FALSE,row.names=FALSE,sep="\t")
}


