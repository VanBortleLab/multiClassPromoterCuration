rm(list = ls())
require(dplyr)
require(data.table)
#########################################OUTPUT FOLDER !!Caution::It deletes existing folders too
system("rm -r mfinal.bed pfinal.bed cat.bed merged.bed sorted.bed intersected.bed")
system("rm -r rnacentral_masterfiles_final")
system("mkdir rnacentral_masterfiles_final")  

######################################### COLLAPSING TSS within 120 nt range
range=120
####+ve and -ve strand collapsed separately
for(i in list.files("./rnacentral_masterfiles_raw",full.names = TRUE,pattern = "*masterfile.bed"))
{
  #i="./masterfiles_new//rRNA_masterfile.bed"
  print(i)
  x=read.table(file=i,sep='\t',header=TRUE,quote="")
  if(ncol(x)<10) next
  x=x %>% mutate(start_median=Start,stop_median=Stop) %>% select(Chr,Start,Stop,URSID,Strand,RNA_type,Databases,start_median,stop_median,everything())
  columns=colnames(x)
  db_used=colnames(x)[10:ncol(x)]
  ###########################################SEPARATING STRANDS AND ELONGATE/SHORTEN THE INTERVAL TO 120
  plus=subset(x,Strand=="+")
  plus$Stop=as.integer(plus$Start)+range
  minus=subset(x,Strand=="-")
  minus$Start=as.integer(minus$Stop)-range
  for( t1 in list("plus","minus"))
  {
       # t1="plus"
     t2=eval((parse(text=paste(t1))))
     if(nrow(t2!=0))
     {
    ###########################################SORTING+MERGING+INTERSECTING
    sorted=t2[order(t2$Chr,t2$Start),,drop=FALSE] #sorting in R because bedtools showed up error while sorting
    write.table(sorted,file="sorted.bed",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE,append=FALSE)
    system(" bedtools merge -i sorted.bed > merged.bed")
    system(" bedtools intersect -a merged.bed -b sorted.bed -wb -wa > intersected.bed")
    ##########################################COLLAPSING ALL INFORMATION IN INTERSECTED.BED FILE
    y=read.table(file="intersected.bed",sep='\t',quote="")
    
    
    ###########################MAKING CODE TO COLLAPSE
    if(t1=="plus")
    {
      code='z=y %>% group_by(V1,V2,V3) %>% summarize(,TSS=paste(V5,collapse="$"),start_median=median(V11),stop_median=median(V12),URSID=paste(na.omit(V7),collapse="$"),Strand=V8[1],RNA_type=V9[1]'
    }else
      {code='z=y %>% group_by(V1,V2,V3) %>% summarize(,TSS=paste(V6,collapse="$"),start_median=median(V11),stop_median=median(V12),URSID=paste(na.omit(V7),collapse="$"),Strand=V8[1],RNA_type=V9[1]'}
    
    for(i2 in (1:length(db_used)))
    { 
      code=paste0(code,",",db_used[i2],'=paste(na.omit(V',12+i2,'),collapse="$")')
      if(i2==length(db_used)) code=paste0(code,")")
    }
    ##########################EXECUTING CODE TO COLLAPSE
    eval(parse(text=code))
    
    
    ###################Shortening 3' end by 120nts
    if(t1=="plus")
    {z$V3=z$V3-range}else
    {z$V2=z$V2+range}	
    
    
    ##collapsing start_median and stop_median, to TSS_median... coz start and stop are related by +-120. TSS is either start or stop
    z= z %>% 
      mutate(TSS_median=as.integer(ifelse(t1=="plus",start_median,stop_median))) %>%
      mutate(STOP_median=as.integer(ifelse(t1=="minus",start_median,stop_median))) %>%
      select(V1,V2,V3,TSS,TSS_median,STOP_median,everything()) %>%
      select(-c(start_median,stop_median))
     
    write.table(z,file=ifelse(t1=="plus","pfinal.bed","mfinal.bed"),sep="\t",row.names = FALSE,col.names = FALSE,quote=FALSE)
     }else
     {
       write.table("",file=ifelse(t1=="plus","pfinal.bed","mfinal.bed"),sep="\t",row.names = FALSE,col.names = FALSE,quote=FALSE)
       }
    
    }
  ###########################################CONCATENATING STRANDS AND NAMING THE COLUMN AND SORTING LEXOGRAPHICALLY
  system( "cat pfinal.bed mfinal.bed > cat.bed")
  ff2=read.table("cat.bed",sep="\t",header=FALSE,quote="")
  colnames(ff2)=c("Chr","Start","Stop","TSS","TSS_median","STOP_median","URSID","Strand","RNA_type",db_used)
  sorted=ff2 %>% select(-c(STOP_median)) %>% arrange(Chr,Start)  #here we remove STOP_median
  write.table(sorted,paste0("./rnacentral_masterfiles_final/",substring(i,30,nchar(i)-4),"_final.bed"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE) ##!!MANUAL
  system("rm -r mfinal.bed pfinal.bed cat.bed merged.bed sorted.bed intersected.bed") #######removing intermediate files
  
}##luka##


#########COLLAPSING all database info into single column i.e. TOTAL_INFO
####!!!If you have STOP_median included in the table above, then you may have to make some changes below while selecting columns to collapse.
FLAG=0
for (i in list.files("./rnacentral_masterfiles_final/",full.names = T))
{ 
  #i="./rnacentral_masterfiles_final/snRNA_masterfile_final.bed" 
  print(i)
  t=fread(i,sep="\t")
  t= t  %>% rowwise() %>% transmute(across(Chr:RNA_type),TOTAL_INFO=paste(c_across(9:ncol(t)),collapse="%")) ##manually specify what are the columns to collapse together to form TOTAL_INFO
  
  if(FLAG==0) {TOTAL=t} else
  {
    TOTAL=rbind(TOTAL,t)
  }
  FLAG=FLAG+1
}##luka##


##########WRITING files
TOTAL$TSS_median=as.integer(TOTAL$TSS_median) #### some medians are in points thus, making them integer
# write.table(TOTAL,"./rnacentral_masterfiles_final/TOTAL_rnacentral_raw.bed",sep="\t",
#             append=F,quote=F,row.names = F,col.names = T)
TOTAL2=TOTAL[grepl("[a-zA-Z]",TOTAL$TOTAL_INFO),]   #####removing those entries that doesnt have any information in any database, thus their total info is simply %%%%%%% .
write.table(TOTAL2,"./rnacentral_masterfiles_final/TOTAL_rnacentral.bed",sep="\t",
            append=F,quote=F,row.names = F,col.names = T) 

