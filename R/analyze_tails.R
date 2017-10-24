#load libraries
library(data.table)
library(ggplot2)
library(plyr)
library(dplyr)

melt_data_localization<-'/home/smaegol/storage/analyses/tail_seq_3/new/LINE_1_RACE_seq_analysis/flowcell2/fastq/processing_out_sabre/test.tsv'

#read tailing information to data.frame using data.table
tails_data_melt<-fread(melt_data_localization,sep="\t",header=T,stringsAsFactors=T,data.table=F,showProgress=TRUE) # read data

tails_data_melt$tailed<-tails_data_melt$tail_length>0 #mark tailed reads
tails_data_melt$mapped<-tails_data_melt$mapping_position!=-1 #mark mapped reads
tails_data_mapped<-tails_data_melt[tails_data_melt$mapped!=0,] #discard unmapped reads

#head(tails_data_mapped)
tails_data_mapped$ref_name_R5<-as.character(tails_data_mapped$ref_name_R5)
tails_data_mapped$ref_name_R3<-as.character(tails_data_mapped$ref_name_R3)
#tails_data_mapped_same_ref<-tails_data_mapped[tails_data_mapped$ref_name_R5==tails_data_mapped$ref_name_R3,]

#mark uridylated reads
tails_data_mapped$uridylated<-FALSE
tails_data_mapped[tails_data_mapped$Utail_length>0,]$uridylated=TRUE

#in fither analyses use only those read which got CTGAC delimiter identified in the clipped fragment
tails_data_mapped_true<-tails_data_mapped[tails_data_mapped$CTGAC_R5>0,]
tails_data_mapped_true$ref_name=tails_data_mapped_true$ref_name_R5 #use ref_name_R5 as ref_name
#remove heterogenous tails from analysis
tails_data_mapped_true_no_hetero = tails_data_mapped_true[-grep("hetero",tails_data_mapped_true$tail_type),]
#remove other type tails from the analysis
tails_data_mapped_true_no_hetero_no_other = tails_data_mapped_true_no_hetero[-grep("other",tails_data_mapped_true_no_hetero$tail_type),]


#treat all AG,UG or UA tails as other_no_tail
tails_data_mapped_true_no_hetero_no_other$tail_type = as.character(tails_data_mapped_true_no_hetero_no_other$tail_type)
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$tail_type=='AG',]$tail_type<-"other_no_tail"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$tail_type=='UG',]$tail_type<-"other_no_tail"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$tail_type=='UA',]$tail_type<-"other_no_tail"
tails_data_mapped_true_no_hetero_no_other<-tails_data_mapped_true_no_hetero_no_other[-grep("other",tails_data_mapped_true_no_hetero_no_other$tail_type),]


#create classes for A-tail lengths (0,1,2-5,6-10,11-20,21-30,30+)
tails_data_mapped_true_no_hetero_no_other$A_length=''
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Atail_length==0,]$A_length="0"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Atail_length==1,]$A_length="1"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Atail_length %in% seq(2,5,1),]$A_length="2-5"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Atail_length %in% seq(6,10,1),]$A_length="6-10"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Atail_length %in% seq(11,20,1),]$A_length="11-20"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Atail_length %in% seq(21,30,1),]$A_length="21-30"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Atail_length>30,]$A_length="30+"
tails_data_mapped_true_no_hetero_no_other$A_length<-factor(tails_data_mapped_true_no_hetero_no_other$A_length,levels=c("0","1","2-5","6-10","11-20","21-30","30+"))

#create classes for U-tail lengths (0,1,2,3-5,6-10,10+)
tails_data_mapped_true_no_hetero_no_other$U_length=''
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Utail_length==0,]$U_length="0"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Utail_length==1,]$U_length="1"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Utail_length==2,]$U_length="2"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Utail_length %in% seq(3,5,1),]$U_length="3-5"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Utail_length %in% seq(6,10,1),]$U_length="6-10"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Utail_length>10,]$U_length="10+"
tails_data_mapped_true_no_hetero_no_other$U_length<-factor(tails_data_mapped_true_no_hetero_no_other$U_length,levels=c("10+","6-10","3-5","2","1","0"))


#modify levels of tail_types to have U_only,A-only,AU or no_tail
tails_data_mapped_true_no_hetero_no_other_tails <- tails_data_mapped_true_no_hetero_no_other
tails_data_mapped_true_no_hetero_no_other_tails$tail_type<-as.character(tails_data_mapped_true_no_hetero_no_other_tails$tail_type)
tails_data_mapped_true_no_hetero_no_other_tails$tail_type<-factor(tails_data_mapped_true_no_hetero_no_other_tails$tail_type,levels=c("U_only","AU","no_tail","A_only"))


#create dataframe with PA1 cell_line data
tails_data_mapped_true_no_hetero_no_other_PA1<-tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$cell_line=='PA1',]
#create dataframe with PA1 cell_line data for knockdown conditions
tails_data_mapped_true_no_hetero_no_other_PA1_tails_KD<-tails_data_mapped_true_no_hetero_no_other_PA1[tails_data_mapped_true_no_hetero_no_other_PA1$condition!='NT',]
#filter out tails longer than 64nt - due to a limit of spike-ins used
tails_data_mapped_true_no_hetero_no_other_PA1_tails_KD<-tails_data_mapped_true_no_hetero_no_other_PA1_tails_KD[tails_data_mapped_true_no_hetero_no_other_PA1_tails_KD$tail_length<=64,]

#create dataframe with PA1 cell_line data for untreated
tails_data_mapped_true_no_hetero_no_other_PA1_tails_NT<-tails_data_mapped_true_no_hetero_no_other_PA1[tails_data_mapped_true_no_hetero_no_other_PA1$condition=='NT',]
#filter out tails longer than 64nt - due to a limit of spike-ins used
tails_data_mapped_true_no_hetero_no_other_PA1_tails_NT<-tails_data_mapped_true_no_hetero_no_other_PA1_tails_NT[tails_data_mapped_true_no_hetero_no_other_PA1_tails_NT$tail_length<=64,]

#create dataframe with data for untreated (irrespective of cell line used)
tails_data_mapped_true_no_hetero_no_other_tails_NT<-tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$condition=='NT' & tails_data_mapped_true_no_hetero_no_other_tails$localization=='TOTAL',]
tails_data_mapped_true_no_hetero_no_other_tails_NT<-tails_data_mapped_true_no_hetero_no_other_tails_NT[tails_data_mapped_true_no_hetero_no_other_tails_NT$tail_length<=64,]


tails_data_mapped_true_no_hetero_no_other_tails_NT1<-tails_data_mapped_true_no_hetero_no_other_tails[(tails_data_mapped_true_no_hetero_no_other_tails$condition=='NT' & tails_data_mapped_true_no_hetero_no_other_tails$localization=='TOTAL') | grepl("CYTO_NT",tails_data_mapped_true_no_hetero_no_other_tails$sample_name) ,]
#filter out tails longer than 64nt - due to a limit of spike-ins used
tails_data_mapped_true_no_hetero_no_other_tails_NT1<-tails_data_mapped_true_no_hetero_no_other_tails_NT1[tails_data_mapped_true_no_hetero_no_other_tails_NT1$tail_length<=64,]
