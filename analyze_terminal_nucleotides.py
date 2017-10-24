#!/usr/bin/env python


#######################################################################################
###																					###
###    Copyright (C) 2017  Pawel Krawczyk (p.krawczyk@ibb.waw.pl)					###
###																					###
###    This program is free software: you can redistribute it and/or modify			###
###    it under the terms of the GNU General Public License as published by			###
###    the Free Software Foundation, either version 3 of the License, or			###
###    (at your option) any later version.											###
###																					###
###    This program is distributed in the hope that it will be useful,				###
###    but WITHOUT ANY WARRANTY; without even the implied warranty of				###
###    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the				###
###    GNU General Public License for more details.									###
###																					###
###    You should have received a copy of the GNU General Public License			###
###    along with this program. If not, see <http://www.gnu.org/licenses/>.			###
###																					###
#######################################################################################


import os, sys
import argparse

#srcipt path is required to find the location of files required for analysis (indexes and other scripts)
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))

#parse command line arguments
parser = argparse.ArgumentParser(description='analyze terminal nucleotides of reads')
parser.add_argument('--inputdir', dest='inputdir', action='store', help='Input dir(required)',required=True)
parser.add_argument('--output', dest='output', action='store', help='Output tsv file (required)',required=True)
parser.add_argument('--window', dest='window', action='store', help='Window size [nucleotides] for 3prime end nucleotides analysis (required)',required=True)
parser.add_argument('--glob', dest='glob', action='store', help='Custom glob for files to analyze (optional)',required=False)
parser.add_argument('--samplesheet', dest='samplesheet', action='store', help='Custom samplesheet location (optional)',required=False)

args = parser.parse_args()


from Bio import SeqIO
import re
import glob
import numpy as np
import pandas as pd
from Bio.Seq import Seq
import subprocess
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#get samplesheet from command-line
samplesheet_location=script_path+'/flowcell2/flowcell2_analysis_samplesheet.csv'
if(args.samplesheet):
	samplesheet_location=args.samplesheet
#read samplesheet into pandas dataframe
data = pd.DataFrame.from_csv(samplesheet_location, sep='\t')

#function returns n last nucleotides of sequence (where n is specified with --windows option)
def get_3end_nucleotides(sequence,window_size):

	terminal_nucleotides = sequence[-int(window_size):]
	return(terminal_nucleotides)


#main processing function
def analyze_tails(R1,R2,transcript,sample_name,localization,replicate,condition,cell_line,primer_name,person):
	#index R2(3'-end) reads
	R2_reads = SeqIO.index(R2, "fastq")

	#dict storing results which will be saved in tsv file
	tails_results = {}
	final_results = {}

	#read all fastq records from R1 (R5) file
	for record in SeqIO.parse(R1, "fastq"):

		seq_id = record.id #get id of read
		tails_results[seq_id]={} #create dict for storing results of pair
		final_results[seq_id]={} #create dict for storing results of pair
		R5_seq = record.seq #get seq of R5 read (after clipping)

		#check if mate is present in the R2 reads file (can be absent in case of rmasker
		if(str(seq_id) in R2_reads):
			record2 = R2_reads[str(seq_id)] # get R3 read from rmasker output
			R3_seq=record2.seq # get seq of R3 read (after clipping)
			R3_seq=R3_seq.reverse_complement()
			seq_R3_length = len(record2.seq)
		else:
			record2=SeqRecord(Seq(''),description=">a0000:00000000:0000:0:0:\tclip5: \tclip3: \tpos: -1\tref: -1")
			R3_seq=''
			seq_R3_length=0

		#get terminal nucleotides information:
		terminal_nucleotides=get_3end_nucleotides(R3_seq,args.window)
		number_U_in_terminal=terminal_nucleotides.count("T")
		number_A_in_terminal=terminal_nucleotides.count("A")
		number_C_in_terminal=terminal_nucleotides.count("C")
		number_G_in_terminal=terminal_nucleotides.count("G")
		ratio_A_in_terminal=number_A_in_terminal/int(args.window)
		ratio_U_in_terminal=number_U_in_terminal/int(args.window)
		ratio_C_in_terminal=number_C_in_terminal/int(args.window)
		ratio_G_in_terminal=number_G_in_terminal/int(args.window)

		#store final results
		final_results[seq_id]['transcript']=transcript
		final_results[seq_id]['cell_line']=cell_line
		final_results[seq_id]['person']=person
		final_results[seq_id]['localization']=localization
		final_results[seq_id]['condition']=condition
		final_results[seq_id]['replicate']=replicate
		final_results[seq_id]['sample_name']=sample_name
		final_results[seq_id]['primer_name']=primer_name
		final_results[seq_id]['number_U']=number_U_in_terminal
		final_results[seq_id]['number_A']=number_A_in_terminal
		final_results[seq_id]['number_C']=number_C_in_terminal
		final_results[seq_id]['number_G']=number_G_in_terminal
		final_results[seq_id]['ratio_U']=ratio_U_in_terminal
		final_results[seq_id]['ratio_A']=ratio_A_in_terminal
		final_results[seq_id]['ratio_C']=ratio_C_in_terminal
		final_results[seq_id]['ratio_G']=ratio_G_in_terminal
		final_results[seq_id]['window_size']=args.window
		final_results[seq_id]['terminal_nucleotides']=terminal_nucleotides


	return final_results

analyzed = 0
os.chdir(args.inputdir)

#define default files to search
files_to_search = "*R5.fastq"

#get files to search from command-line (if present)
if (args.glob):
	files_to_search = args.glob

for R5_file in glob.glob(files_to_search):
#iterate through R5 files
	analyzed = analyzed + 1 #increment number of analyzed files
	print(R5_file)
	file_parts = re.search("(?P<path>.*//|)(?P<basename>.*fastq)",R5_file)
	file_basename=file_parts.group("basename") #get file basename
	file_path = file_parts.group("path")
	file_parts = re.search("(?P<prefix>.*)_R5(?P<suffix>.*)",file_basename)
	file_prefix = file_parts.group("prefix")
	file_suffix = file_parts.group("suffix")
	R5_fasta_file = file_prefix + "_.fasta"
	R3_file = file_prefix + "_R3" + file_suffix #create R3 file name
	#get information about sample pair from samplesheet:
	temp = data[data["R5_file"]==R5_file]
	transcript=temp['transcript'][0] #transcript name
	cell_line=temp['cell_line'][0] #cell line
	condition=temp['condition'][0] #condition
	person=temp['person'][0] #person - who prepared the library
	localization=temp['localization'][0] #subcellular localization (fractionation experiment)
	replicate=temp['replicate'][0] #replicate number
	sample_name=temp['Sample_Name'][0] #sample name
	primer_name=temp['primer_name'][0] #name of the primer used in RACE
	#if bowtie was not run before - run bowtie2 on R5 and R3 files:
	#check if there is any sequence in R5 file, if not - skip analysis for this pair
	R5_records = list(SeqIO.parse(R5_file, "fastq"))

	if(len(R5_records)==0):
		print("skipping analysis because no sequences found in the R5 file")
	else:
		#Run the analysis of tails:
		paired_results = analyze_tails(R5_file,R3_file,transcript,sample_name,localization,replicate,condition,cell_line,primer_name,person)
		#Create pandas data frame with result
		tails_df = pd.DataFrame.from_dict(paired_results,orient='index')

		#make sure data are saved after each library processed:
		if (analyzed > 1):
			tails_df.to_csv(args.output, mode='a', sep='\t',header=False)
		else:
			tails_df.to_csv(args.output, mode='w', sep='\t',header=True)


print("all " + str(analyzed) + " samples analyzed succesfully\n")

## END ##
