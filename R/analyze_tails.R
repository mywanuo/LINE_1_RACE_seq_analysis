#######################################################################################
###                                                                                 ###
###     Copyright (C) 2017  Pawel Krawczyk (p.krawczyk@ibb.waw.pl)                  ###
###                                                                                 ###
###     This program is free software: you can redistribute it and/or modify        ###
###     it under the terms of the GNU General Public License as published by        ###
###     the Free Software Foundation, either version 3 of the License, or           ###
###     (at your option) any later version.                                         ###
###                                                                                 ###
###     This program is distributed in the hope that it will be useful,             ###
###     but WITHOUT ANY WARRANTY; without even the implied warranty of              ###
###     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               ###
###     GNU General Public License for more details.                                ###
###                                                                                 ###
###     You should have received a copy of the GNU General Public License           ###
###     along with this program. If not, see <http://www.gnu.org/licenses/>.        ###
###                                                                                 ###
#######################################################################################


# load libraries
library(data.table)
library(ggplot2)
library(plyr)
library(dplyr)
library(ggseqlogo)

melt_data_localization <- "/home/smaegol/storage/analyses/tail_seq_3/new/LINE_1_RACE_seq_analysis/flowcells_all/processing_out/test_fin.tsv"


# read tailing information to data.frame using data.table
tails_data_melt <- fread(melt_data_localization, sep = "\t", header = T, stringsAsFactors = T,
    data.table = F, showProgress = TRUE)  # read data


#tails_data_melt[tails_data_melt$localization=='JADRO',]$localization<-"TOTAL"
tails_data_melt$tailed <- tails_data_melt$tail_length > 0  #mark tailed reads
tails_data_melt$mapped <- tails_data_melt$mapping_position != -1  #mark mapped reads
tails_data_mapped <- tails_data_melt[tails_data_melt$mapped != 0, ]  #discard unmapped reads

#head(tails_data_mapped$transcript)

# fix of flase_no_tail - should be assigned no tail cause CTGAC was found in the
# clipped fragment - to be fixed in python script
tails_data_mapped$tail_type <- as.character(tails_data_mapped$tail_type)
tails_data_mapped[tails_data_mapped$tail_type == "false_no_tail_no_CTGAC", ]$tail_type <- "no_tail"

tails_data_mapped$ref_name_R5 <- as.character(tails_data_mapped$ref_name_R5)
tails_data_mapped$ref_name_R3 <- as.character(tails_data_mapped$ref_name_R3)
# tails_data_mapped_same_ref<-tails_data_mapped[tails_data_mapped$ref_name_R5==tails_data_mapped$ref_name_R3,]

# mark uridylated reads
tails_data_mapped$uridylated <- FALSE
tails_data_mapped[tails_data_mapped$Utail_length > 0, ]$uridylated = TRUE

#convert T to U in terminal nucleotides (for seqlogo)
tails_data_mapped$terminal_nucleotides<-as.character(tails_data_mapped$terminal_nucleotides)
tails_data_mapped$terminal_nucleotides<-gsub("T","U",tails_data_mapped$terminal_nucleotides)

# in further analyses use only those read which got CTGAC delimiter identified in
# the clipped fragment
# for reporter analyses all mapped reads got CTGAC_R5 variable = 1 (because of short reads) so they will be included in the analysis
tails_data_mapped_true <- tails_data_mapped[tails_data_mapped$CTGAC_R5 > 0, ]
tails_data_mapped_true$ref_name = tails_data_mapped_true$ref_name_R5  #use ref_name_R5 as ref_name
# remove heterogenous tails from analysis
tails_data_mapped_true_no_hetero = tails_data_mapped_true[-grep("hetero", tails_data_mapped_true$tail_type),
    ]
# remove other type tails (for which we can suspect they are not tails but rather origin from improper mapping/repeatmasker) from the analysis
tails_data_mapped_true_no_hetero_no_other = tails_data_mapped_true_no_hetero[-grep("other",
    tails_data_mapped_true_no_hetero$tail_type), ]


# treat all AG,UG or UA tails as other_no_tail
tails_data_mapped_true_no_hetero_no_other$tail_type = as.character(tails_data_mapped_true_no_hetero_no_other$tail_type)
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$tail_type ==
    "AG", ]$tail_type <- "other_no_tail"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$tail_type ==
    "UG", ]$tail_type <- "other_no_tail"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$tail_type ==
    "UA", ]$tail_type <- "other_no_tail"
tails_data_mapped_true_no_hetero_no_other <- tails_data_mapped_true_no_hetero_no_other[-grep("other",
    tails_data_mapped_true_no_hetero_no_other$tail_type), ] #remove all other from analysis


# create classes for A-tail lengths (0,1,2-5,6-10,11-20,21-30,30+)
tails_data_mapped_true_no_hetero_no_other$A_length = ""
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Atail_length ==
    0, ]$A_length = "0"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Atail_length ==
    1, ]$A_length = "1"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Atail_length %in%
    seq(2, 5, 1), ]$A_length = "2-5"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Atail_length %in%
    seq(6, 10, 1), ]$A_length = "6-10"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Atail_length %in%
    seq(11, 20, 1), ]$A_length = "11-20"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Atail_length %in%
    seq(21, 30, 1), ]$A_length = "21-30"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Atail_length >
    30, ]$A_length = "30+"
tails_data_mapped_true_no_hetero_no_other$A_length <- factor(tails_data_mapped_true_no_hetero_no_other$A_length,
    levels = c("0", "1", "2-5", "6-10", "11-20", "21-30", "30+"))

# create classes for U-tail lengths (0,1,2,3-5,6-10,10+)
tails_data_mapped_true_no_hetero_no_other$U_length = ""
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Utail_length ==
    0, ]$U_length = "0"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Utail_length ==
    1, ]$U_length = "1"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Utail_length ==
    2, ]$U_length = "2"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Utail_length %in%
    seq(3, 5, 1), ]$U_length = "3-5"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Utail_length %in%
    seq(6, 10, 1), ]$U_length = "6-10"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Utail_length >
    10, ]$U_length = "10+"
tails_data_mapped_true_no_hetero_no_other$U_length <- factor(tails_data_mapped_true_no_hetero_no_other$U_length,
    levels = c("10+", "6-10", "3-5", "2", "1", "0"))


# modify levels of tail_types to have U_only,A-only,AU or no_tail
tails_data_mapped_true_no_hetero_no_other_tails <- tails_data_mapped_true_no_hetero_no_other
tails_data_mapped_true_no_hetero_no_other_tails$tail_type <- as.character(tails_data_mapped_true_no_hetero_no_other_tails$tail_type)
tails_data_mapped_true_no_hetero_no_other_tails$tail_type <- factor(tails_data_mapped_true_no_hetero_no_other_tails$tail_type,
    levels = c("U_only", "AU", "no_tail", "A_only"))

# create dataframe with PA1 cell_line data
tails_data_mapped_true_no_hetero_no_other_PA1 <- tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$cell_line ==
    "PA1", ]
# create dataframe with PA1 cell_line data for knockdown conditions
tails_data_mapped_true_no_hetero_no_other_PA1_tails_KD <- tails_data_mapped_true_no_hetero_no_other_PA1[tails_data_mapped_true_no_hetero_no_other_PA1$condition !=
    "NT", ]
# filter out tails longer than 64nt - due to a limit of spike-ins used
tails_data_mapped_true_no_hetero_no_other_PA1_tails_KD <- tails_data_mapped_true_no_hetero_no_other_PA1_tails_KD[tails_data_mapped_true_no_hetero_no_other_PA1_tails_KD$tail_length <=
    64, ]

# create dataframe with PA1 cell_line data for untreated
tails_data_mapped_true_no_hetero_no_other_PA1_tails_NT <- tails_data_mapped_true_no_hetero_no_other_PA1[tails_data_mapped_true_no_hetero_no_other_PA1$condition ==
    "NT", ]
# filter out tails longer than 64nt - due to a limit of spike-ins used
tails_data_mapped_true_no_hetero_no_other_PA1_tails_NT <- tails_data_mapped_true_no_hetero_no_other_PA1_tails_NT[tails_data_mapped_true_no_hetero_no_other_PA1_tails_NT$tail_length <=
    64, ]

# create dataframe with data for untreated (irrespective of cell line used)
tails_data_mapped_true_no_hetero_no_other_tails_NT <- tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$condition ==
    "NT" & tails_data_mapped_true_no_hetero_no_other_tails$localization == "TOTAL",
    ]
tails_data_mapped_true_no_hetero_no_other_tails_NT <- tails_data_mapped_true_no_hetero_no_other_tails_NT[tails_data_mapped_true_no_hetero_no_other_tails_NT$tail_length <=
    64, ]


tails_data_mapped_true_no_hetero_no_other_tails_NT1 <- tails_data_mapped_true_no_hetero_no_other_tails[(tails_data_mapped_true_no_hetero_no_other_tails$condition ==
    "NT" & tails_data_mapped_true_no_hetero_no_other_tails$localization %in% c("TOTAL","JADRO")),
    ]
# filter out tails longer than 64nt - due to a limit of spike-ins used
tails_data_mapped_true_no_hetero_no_other_tails_NT1 <- tails_data_mapped_true_no_hetero_no_other_tails_NT1[tails_data_mapped_true_no_hetero_no_other_tails_NT1$tail_length <=
    64, ]


tails_data_for_analysis_reporter_overexp <- tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$transcript=="REPORTERL1_overexp",]
tails_data_for_analysis_reporter_overexp <- tails_data_for_analysis_reporter_overexp[tails_data_for_analysis_reporter_overexp$tail_length <=
    64, ]

tails_data_for_analysis_reporter_KD <- tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$transcript=="REPORTERL1",]
tails_data_for_analysis_reporter_KD <- tails_data_for_analysis_reporter_KD[tails_data_for_analysis_reporter_KD$tail_length <=
    64, ]


### functions for summarizing data ###

summarize_tails_by_tail_type <- function(tail_data2) {
    # function summarizing tails by tail type
    tail_data <- tail_data2
    tail_data_summarized <- ddply(tail_data, .(transcript, condition, cell_line,
        localization, primer_name, tail_type), summarise, N = length(primer_name))
    tail_data_summarized <- ddply(tail_data_summarized, .(transcript, condition,
        cell_line, localization, primer_name), transform, freq = N/sum(N))
    return(tail_data_summarized)
}



tails_data_mapped_reporter_overexp_for_terminal <- tails_data_mapped[tails_data_mapped$transcript=='REPORTERL1_overexp',]
tails_data_mapped_reporter_overexp_for_terminal$condition<-as.character(tails_data_mapped_reporter_overexp_for_terminal$condition)
tails_data_mapped_reporter_overexp_for_terminal$condition<-as.factor(tails_data_mapped_reporter_overexp_for_terminal$condition)


tails_data_mapped_gapdh_overexp_for_terminal <- tails_data_mapped[tails_data_mapped$transcript=='GAPDH' & tails_data_mapped$condition %in% c("CNTRL","TUT7WT","MOV10"),]
tails_data_mapped_gapdh_overexp_for_terminal$condition<-as.character(tails_data_mapped_gapdh_overexp_for_terminal$condition)
tails_data_mapped_gapdh_overexp_for_terminal$condition<-as.factor(tails_data_mapped_gapdh_overexp_for_terminal$condition)


for (cond in levels(tails_data_mapped_reporter_overexp_for_terminal$condition)) {
  print(cond)
  terminal_nucleotides_cond_name=paste("terminal_nucleotides_reporter_overexp",cond,sep="_")
  #temp2=eval(as.symbol(summary_tail_lengths_table_name))
  assign(terminal_nucleotides_cond_name,tails_data_mapped_reporter_overexp_for_terminal[tails_data_mapped_reporter_overexp_for_terminal$condition==cond,]$terminal_nucleotides)
}

terminal_nucleotides_reporter_overexp=list(CTRL=terminal_nucleotides_reporter_overexp_CNTRL,TUT7WT=terminal_nucleotides_reporter_overexp_TUT7WT,MOV10=terminal_nucleotides_reporter_overexp_MOV10,TUT4WT=terminal_nucleotides_reporter_overexp_TUT4WT,TUT7MUT=terminal_nucleotides_reporter_overexp_TUT7MUT,TUT4MUT=terminal_nucleotides_reporter_overexp_TUT4MUT)

for (cond in levels(tails_data_mapped_gapdh_overexp_for_terminal$condition)) {
  print(cond)
  terminal_nucleotides_cond_name=paste("terminal_nucleotides_gapdh_overexp",cond,sep="_")
  #temp2=eval(as.symbol(summary_tail_lengths_table_name))
  assign(terminal_nucleotides_cond_name,tails_data_mapped_gapdh_overexp_for_terminal[tails_data_mapped_gapdh_overexp_for_terminal$condition==cond,]$terminal_nucleotides)
}

terminal_nucleotides_gapdh_overexp=list(CTRL=terminal_nucleotides_gapdh_overexp_CNTRL,TUT7WT=terminal_nucleotides_gapdh_overexp_TUT7WT,MOV10=terminal_nucleotides_gapdh_overexp_MOV10)

setwd("/home/smaegol/")
pdf("fig_4a_reporter.pdf")
print(ggseqlogo(terminal_nucleotides_reporter_overexp,ncol=2,method='prob'))
dev.off()

pdf("fig_4a_gapdh.pdf")
print(ggseqlogo(terminal_nucleotides_gapdh_overexp,ncol=2,method='prob'))
dev.off()


pdf("fig_4B.pdf")
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_tails_by_tail_type(tails_data_for_analysis_reporter_overexp[tails_data_for_analysis_reporter_overexp$condition %in%
    c("CNTRL","TUT7WT","TUT4WT","MOV10"), ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,levels=c("CNTRL","TUT7WT","TUT4WT","MOV10"))
plot_tail_lengths <- ggplot(summary_tail_types_table, aes(x = as.factor(condition),
    fill = tail_type, colours = tail_type)) + geom_bar(aes(y = freq), position = position_stack(),
    stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
    scale_fill_grey() + ggtitle(paste("OVEREXPRESSION"))
print(plot_tail_lengths)
dev.off()

pdf("fig_4C.pdf")
summary_tail_types_table_name <- paste("repo_KD_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_tails_by_tail_type(tails_data_for_analysis_reporter_KD[tails_data_for_analysis_reporter_KD$condition %in%
    c("CNTRLKD","TUT7KD","TUT4KD","TUT4TUT7KD"), ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,levels=c("CNTRLKD","TUT7KD","TUT4KD","TUT4TUT7KD"))
plot_tail_lengths <- ggplot(summary_tail_types_table, aes(x = as.factor(condition),
    fill = tail_type, colours = tail_type)) + geom_bar(aes(y = freq), position = position_stack(),
    stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
    scale_fill_grey() + ggtitle(paste("DEPLETION"))
print(plot_tail_lengths)
dev.off()


pdf("fig_4D.pdf")
summary_tail_types_table_name <- paste("genomic_L1_summarized_tails_types_by_cell_lines")
assign(summary_tail_types_table_name, summarize_tails_by_tail_type(tails_data_mapped_true_no_hetero_no_other_tails_NT1[tails_data_mapped_true_no_hetero_no_other_tails_NT1$transcript=="ENDOL1" & tails_data_mapped_true_no_hetero_no_other_tails_NT1$cell_line %in% c("293T","H9","HELAHA","NPC","PA1","MYSZ"), ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$cell_line <- as.character(summary_tail_types_table$cell_line)
summary_tail_types_table[summary_tail_types_table$cell_line=='MYSZ',]$cell_line="MOUSE_TESTIS"
summary_tail_types_table$cell_line <- factor(summary_tail_types_table$cell_line,levels=c("293T","HELAHA","NPC","PA1","H9","MOUSE_TESTIS"))
plot_tail_lengths <- ggplot(summary_tail_types_table, aes(x = as.factor(cell_line),
    fill = tail_type, colours = tail_type)) + geom_bar(aes(y = freq), position = position_stack(),
    stat = "identity") + scale_y_continuous() + xlab("cell_line") + ylab("fraction of transcripts") +
    scale_fill_grey() + ggtitle(paste("DEPLETION"))
print(plot_tail_lengths)
dev.off()pdf("fig_4F.pdf")
summary_tail_types_table_name <- paste("genomic_L1_summarized_tails_types_by_condition_KD")
assign(summary_tail_types_table_name, summarize_tails_by_tail_type(tails_data_mapped_true_no_hetero_no_other_PA1_tails_KD[tails_data_mapped_true_no_hetero_no_other_PA1_tails_KD$condition %in%
    c("CNTRLKD","MOV10KD","TUT4TUT7KD") & tails_data_mapped_true_no_hetero_no_other_PA1_tails_KD$transcript=="ENDOL1", ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,levels=c("CNTRLKD","TUT4TUT7KD","MOV10KD"))
plot_tail_lengths <- ggplot(summary_tail_types_table, aes(x = as.factor(condition),
    fill = tail_type, colours = tail_type)) + geom_bar(aes(y = freq), position = position_stack(),
    stat = "identity") + scale_y_continuous() + xlab("condition") + ylab("fraction of transcripts") +
    scale_fill_grey() + ggtitle(paste("DEPLETION (PA-1)"))
print(plot_tail_lengths)
dev.off()
