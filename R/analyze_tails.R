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
library(RColorBrewer)

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
    "NT" & tails_data_mapped_true_no_hetero_no_other_PA1$localization %in% c("CYTO","NUC") & tails_data_mapped_true_no_hetero_no_other_PA1$primer_name %in% c("L1NGS0","GAPDH"), ]
# filter out tails longer than 64nt - due to a limit of spike-ins used
tails_data_mapped_true_no_hetero_no_other_PA1_tails_NT <- tails_data_mapped_true_no_hetero_no_other_PA1_tails_NT[tails_data_mapped_true_no_hetero_no_other_PA1_tails_NT$tail_length <=
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


summarize_tails_by_tail_type_loc <- function(tail_data2) {
    # function summarizing tails by tail type
    tail_data <- tail_data2
    tail_data_summarized <- ddply(tail_data, .(transcript, condition, cell_line,
        localization, primer_name, tail_type), summarise, N = length(primer_name))
    tail_data_summarized <- ddply(tail_data_summarized, .(transcript, condition,
        cell_line, primer_name), transform, freq = N/sum(N))
    return(tail_data_summarized)
}

summarize_Utails_calculate_means <- function(tail_data2) {
    # function summarizing Utail length data for given dataset
    tail_data <- tail_data2
    tail_data_summarized <- ddply(tail_data, .(transcript, condition, cell_line,
        primer_name), summarise, N = length(Utail_length), mean_Utail = mean(Utail_length,
        na.rm = T), sd_utail = sd(Utail_length))
    tail_data_summarized <- ddply(tail_data_summarized, .(transcript, condition,
        cell_line, primer_name), transform, se_utail = sd_utail/sqrt(N))
    return(tail_data_summarized)
}


summarize_Utails_lengths <- function(tail_data2) {
    # function summarizing tail length data for given dataset
    tail_data <- tail_data2
    tail_data_summarized <- ddply(tail_data, .(transcript, condition, cell_line,
        localization, primer_name, U_length), summarise, N = length(U_length))
    tail_data_summarized <- ddply(tail_data_summarized, .(transcript, condition,
        cell_line, localization, primer_name), transform, freq = N/sum(N))
    return(tail_data_summarized)
}



summarize_tails_by_tail_type_collapse_short <- function(tail_data2) {
    # function summarizing tail length data for given dataset
    tail_data <- tail_data2
    tail_data$tail_length <- as.numeric(tail_data$tail_length)
    max_tail_length <- summary(tail_data2$tail_length)[6]
    print(max_tail_length)
    if (max_tail_length >= 60) {
        tail_data[tail_data$tail_length >= 60, ]$tail_length <- "60+"
    }
    if (max_tail_length >= 50) {
        for (temp_len in seq(50, 59, 1)) {
            if (any(tail_data$tail_length == temp_len)) {
                tail_data[tail_data$tail_length == temp_len, ]$tail_length <- "50-59"
            }
        }
    }
    if (max_tail_length >= 40) {
        for (temp_len in seq(40, 49, 1)) {
            if (any(tail_data$tail_length == temp_len)) {
                tail_data[tail_data$tail_length == temp_len, ]$tail_length <- "40-49"
            }
        }
    }
    if (max_tail_length >= 30) {
        for (temp_len in seq(30, 39, 1)) {
            if (any(tail_data$tail_length == temp_len)) {
                tail_data[tail_data$tail_length == temp_len, ]$tail_length <- "30-39"
            }
        }
    }
    if (max_tail_length >= 20) {
        for (temp_len in seq(20, 29, 1)) {
            if (any(tail_data$tail_length == temp_len)) {
                tail_data[tail_data$tail_length == temp_len, ]$tail_length <- "20-29"
            }
        }
    }
    for (temp_len in seq(10, 19, 1)) {
        if (any(tail_data$tail_length == temp_len)) {
            tail_data[tail_data$tail_length == temp_len, ]$tail_length <- "10-19"
        }
    }

    for (temp_len in seq(1, 9, 1)) {
        if (any(tail_data$tail_length == temp_len)) {
            tail_data[tail_data$tail_length == temp_len, ]$tail_length <- "1-9"
        }
    }
    tail_data$tail_length <- as.factor(tail_data$tail_length)
    tail_data$tail_length <- factor(tail_data$tail_length, levels = c("0", "1-9",
        "10-19", "20-29", "30-39", "40-49", "50-59", "60+"))
    tail_data_summarized <- ddply(tail_data, .(transcript, condition, tail_length,
        cell_line, localization, primer_name, tail_type), summarise, N = length(tail_length))
    tail_data_summarized <- ddply(tail_data_summarized, .(transcript, condition,
        cell_line, localization, primer_name), transform, freq = N/sum(N))
    return(tail_data_summarized)
}



## FIGURES ##


tails_data_mapped_reporter_overexp_for_terminal <- tails_data_mapped[tails_data_mapped$transcript ==
    "REPORTERL1_overexp", ]
tails_data_mapped_reporter_overexp_for_terminal$condition <- as.character(tails_data_mapped_reporter_overexp_for_terminal$condition)
tails_data_mapped_reporter_overexp_for_terminal$condition <- as.factor(tails_data_mapped_reporter_overexp_for_terminal$condition)


tails_data_mapped_gapdh_overexp_for_terminal <- tails_data_mapped[tails_data_mapped$transcript ==
    "GAPDH" & tails_data_mapped$condition %in% c("CNTRL", "TUT7WT", "MOV10"), ]
tails_data_mapped_gapdh_overexp_for_terminal$condition <- as.character(tails_data_mapped_gapdh_overexp_for_terminal$condition)
tails_data_mapped_gapdh_overexp_for_terminal$condition <- as.factor(tails_data_mapped_gapdh_overexp_for_terminal$condition)


for (cond in levels(tails_data_mapped_reporter_overexp_for_terminal$condition)) {
    print(cond)
    terminal_nucleotides_cond_name = paste("terminal_nucleotides_reporter_overexp",
        cond, sep = "_")
    # temp2=eval(as.symbol(summary_tail_lengths_table_name))
    assign(terminal_nucleotides_cond_name, tails_data_mapped_reporter_overexp_for_terminal[tails_data_mapped_reporter_overexp_for_terminal$condition ==
        cond, ]$terminal_nucleotides)
}

terminal_nucleotides_reporter_overexp = list(CTRL = terminal_nucleotides_reporter_overexp_CNTRL,
    TUT7WT = terminal_nucleotides_reporter_overexp_TUT7WT, MOV10 = terminal_nucleotides_reporter_overexp_MOV10,
    TUT4WT = terminal_nucleotides_reporter_overexp_TUT4WT, TUT7MUT = terminal_nucleotides_reporter_overexp_TUT7MUT,
    TUT4MUT = terminal_nucleotides_reporter_overexp_TUT4MUT)

for (cond in levels(tails_data_mapped_gapdh_overexp_for_terminal$condition)) {
    print(cond)
    terminal_nucleotides_cond_name = paste("terminal_nucleotides_gapdh_overexp",
        cond, sep = "_")
    # temp2=eval(as.symbol(summary_tail_lengths_table_name))
    assign(terminal_nucleotides_cond_name, tails_data_mapped_gapdh_overexp_for_terminal[tails_data_mapped_gapdh_overexp_for_terminal$condition ==
        cond, ]$terminal_nucleotides)
}

terminal_nucleotides_gapdh_overexp = list(CTRL = terminal_nucleotides_gapdh_overexp_CNTRL,
    TUT7WT = terminal_nucleotides_gapdh_overexp_TUT7WT, MOV10 = terminal_nucleotides_gapdh_overexp_MOV10)

setwd("/home/smaegol/")
pdf("fig_4a_reporter.pdf")
print(ggseqlogo(terminal_nucleotides_reporter_overexp, ncol = 2, method = "prob"))
dev.off()

pdf("fig_4a_gapdh.pdf")
print(ggseqlogo(terminal_nucleotides_gapdh_overexp, ncol = 2, method = "prob"))
dev.off()


pdf("fig_4B.pdf")
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_tails_by_tail_type(tails_data_for_analysis_reporter_overexp[tails_data_for_analysis_reporter_overexp$condition %in%
    c("CNTRL", "TUT7WT", "TUT4WT", "MOV10"), ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
    levels = c("CNTRL", "TUT7WT", "TUT4WT", "MOV10"))
plot_tail_lengths <- ggplot(summary_tail_types_table, aes(x = as.factor(condition),
    fill = tail_type, colours = tail_type)) + geom_bar(aes(y = freq), position = position_stack(),
    stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
    scale_fill_grey() + ggtitle(paste("OVEREXPRESSION"))
print(plot_tail_lengths)
dev.off()

pdf("fig_4C.pdf")
summary_tail_types_table_name <- paste("repo_KD_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_tails_by_tail_type(tails_data_for_analysis_reporter_KD[tails_data_for_analysis_reporter_KD$condition %in%
    c("CNTRLKD", "TUT7KD", "TUT4KD", "TUT4TUT7KD"), ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
    levels = c("CNTRLKD", "TUT7KD", "TUT4KD", "TUT4TUT7KD"))
plot_tail_lengths <- ggplot(summary_tail_types_table, aes(x = as.factor(condition),
    fill = tail_type, colours = tail_type)) + geom_bar(aes(y = freq), position = position_stack(),
    stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
    scale_fill_grey() + ggtitle(paste("DEPLETION"))
print(plot_tail_lengths)
dev.off()


pdf("fig_4D.pdf")
summary_tail_types_table_name <- paste("genomic_L1_summarized_tails_types_by_cell_lines")
assign(summary_tail_types_table_name, summarize_tails_by_tail_type(tails_data_mapped_true_no_hetero_no_other_tails_NT1[tails_data_mapped_true_no_hetero_no_other_tails_NT1$transcript ==
    "ENDOL1" & tails_data_mapped_true_no_hetero_no_other_tails_NT1$cell_line %in%
    c("293T", "H9", "HELAHA", "NPC", "PA1", "MYSZ"), ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$cell_line <- as.character(summary_tail_types_table$cell_line)
summary_tail_types_table[summary_tail_types_table$cell_line == "MYSZ", ]$cell_line = "MOUSE_TESTIS"
summary_tail_types_table$cell_line <- factor(summary_tail_types_table$cell_line,
    levels = c("293T", "HELAHA", "NPC", "PA1", "H9", "MOUSE_TESTIS"))
plot_tail_lengths <- ggplot(summary_tail_types_table, aes(x = as.factor(cell_line),
    fill = tail_type, colours = tail_type)) + geom_bar(aes(y = freq), position = position_stack(),
    stat = "identity") + scale_y_continuous() + xlab("cell_line") + ylab("fraction of transcripts") +
    scale_fill_grey() + ggtitle(paste("DEPLETION"))
print(plot_tail_lengths)
dev.off()

pdf("fig_4F.pdf")
summary_tail_types_table_name <- paste("genomic_L1_summarized_tails_types_by_condition_KD")
assign(summary_tail_types_table_name, summarize_tails_by_tail_type(tails_data_mapped_true_no_hetero_no_other_PA1_tails_KD[tails_data_mapped_true_no_hetero_no_other_PA1_tails_KD$condition %in%
    c("CNTRLKD", "MOV10KD", "TUT4TUT7KD") & tails_data_mapped_true_no_hetero_no_other_PA1_tails_KD$transcript ==
    "ENDOL1", ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
    levels = c("CNTRLKD", "TUT4TUT7KD", "MOV10KD"))
plot_tail_lengths <- ggplot(summary_tail_types_table, aes(x = as.factor(condition),
    fill = tail_type, colours = tail_type)) + geom_bar(aes(y = freq), position = position_stack(),
    stat = "identity") + scale_y_continuous() + xlab("condition") + ylab("fraction of transcripts") +
    scale_fill_grey() + ggtitle(paste("DEPLETION (PA-1)"))
print(plot_tail_lengths)
dev.off()


pdf("fig_4G.pdf")
summary_tail_types_table_name <- paste("genomic_L1_and_GAPDH_summarized_tails_types_by_subcellular_localization")
assign(summary_tail_types_table_name, summarize_tails_by_tail_type_loc(tails_data_mapped_true_no_hetero_no_other_PA1_tails_NT))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
plot_tail_lengths <- ggplot(summary_tail_types_table, aes(x = as.factor(localization),
    fill = tail_type, colours = tail_type)) + geom_bar(aes(y = freq), position = position_stack(),
    stat = "identity") + scale_y_continuous() + xlab("localization") + ylab("fraction of transcripts") +
    scale_fill_grey() + facet_grid(. ~ transcript)
print(plot_tail_lengths)
dev.off()


pdf("fig_S5A.pdf")
summary_mean_Utail_table_name <- paste("reporter_mean_Utails")
assign(summary_mean_Utail_table_name, summarize_Utails_calculate_means(tails_data_for_analysis_reporter_overexp[tails_data_for_analysis_reporter_overexp$condition %in%
    c("CNTRL", "TUT7WT", "TUT7MUT", "MOV10") & tails_data_for_analysis_reporter_overexp$uridylated ==
    TRUE, ]))
summary_mean_Utail_table <- eval(as.symbol(summary_mean_Utail_table_name))
summary_mean_Utail_table$condition <- as.character(summary_mean_Utail_table$condition)
summary_mean_Utail_table$condition <- factor(summary_mean_Utail_table$condition,
    levels = c("CNTRL", "TUT7WT", "TUT7MUT", "MOV10"))
plot_tail_lengths <- ggplot(summary_mean_Utail_table, aes(x = as.factor(condition))) +
    geom_bar(aes(y = mean_Utail), position = position_dodge(), stat = "identity") +
    xlab("condition") + ylab("Utail mean length") + scale_fill_grey() + geom_errorbar(aes(ymin = mean_Utail -
    se_utail, ymax = mean_Utail + se_utail), colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()


pdf("fig_S5B.pdf")
summary_Utail_lengths_table_name <- paste("reporter_overexp_summarized_Utails_lengths_by_condition")
assign(summary_Utail_lengths_table_name, summarize_Utails_lengths(tails_data_for_analysis_reporter_overexp))
summary_Utail_lengths_table <- eval(as.symbol(summary_Utail_lengths_table_name))
summary_Utail_lengths_table$condition <- as.character(summary_Utail_lengths_table$condition)
summary_Utail_lengths_table$condition <- factor(summary_Utail_lengths_table$condition,
    levels = c("CNTRL", "MOV10", "TUT7WT", "TUT7MUT", "TUT4WT", "TUT4MUT"))
plot_tail_lengths <- ggplot(summary_Utail_lengths_table, aes(x = as.factor(condition),
    fill = U_length, colours = U_length)) + geom_bar(aes(y = freq), position = position_stack(),
    stat = "identity") + xlab("cell line") + ylab("fraction of transcripts") + scale_fill_grey()
print(plot_tail_lengths)
dev.off()

pdf("fig_S5C-G")
for (transcript in c("ENDOL1")) {
    for (cell_line in c("H9", "HELAHA", "293T", "NPC", "PA1")) {
        summary_tail_lengths_table_name <- paste("_summarized_tails_lengths_by_type")
        assign(summary_tail_lengths_table_name, summarize_tails_by_tail_type_collapse_short(tails_data_mapped_true_no_hetero_no_other_tails_NT1[tails_data_mapped_true_no_hetero_no_other_tails_NT1$transcript ==
            transcript & tails_data_mapped_true_no_hetero_no_other_tails_NT1$cell_line ==
            cell_line & tails_data_mapped_true_no_hetero_no_other_tails_NT1$primer_name ==
            "L1NGS0", ]))
        summary_tail_lengths_table <- eval(as.symbol(summary_tail_lengths_table_name))
        plot_tail_lengths <- ggplot(summary_tail_lengths_table, aes(x = as.factor(tail_length),
            fill = tail_type, colours = tail_type)) + geom_bar(aes(y = freq), position = position_stack(),
            stat = "identity") + scale_y_continuous(labels = percent) + xlab("Atail_length") +
            ylab("% of 3' ends") + facet_grid(. ~ cell_line) + scale_fill_brewer(palette = "Spectral")
        print(plot_tail_lengths)
    }
}
dev.off()
