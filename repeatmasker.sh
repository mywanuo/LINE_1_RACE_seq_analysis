#!/bin/bash

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

#make sure fastx_toolkit is in the path
module load fastx_toolkit

#before start - set paths for software and libraries:
repeat_masker=`which RepeatMasker`
repeat_masker_parsing="/home/smaegol/storage/soft/repeatmasker/Parsing-RepeatMasker-Outputs/parseRM_simple.pl -fast"
clip_rmasker=`pwd`"/identify_LINE_repeatmasker_softclip.py"
clip_rmasker_R3=`pwd`"/identify_LINE_repeatmasker_softclip_R3.py"
lib_location="/home/smaegol/storage/soft/repeatmasker/RepeatMasker/Libraries/Dfam_2.0/homo_sapiens/LINEs/masklib.hmm"

#options for repeatmasker
threads=5
max_divergence=10

#find sequence files in the currect folder based on name (for LINE RACE seqs - begin with L1)
#we start processing with R5 files
for f in `find . -name "L1*R5.fastq"`
do
 echo "Processing file: $f "
 FILENAME_PREFIX=`expr match "$f" '\(.*\)R5'`
 FILENAME_SUFFIX=`expr match "$f" '.*R5\(.*\)'`
 fasta=$FILENAME_PREFIX"R5.fasta"
 rmasker_out=$fasta".out"
 rmasker_parsed=$rmasker_out".parsed1/"$rmasker_out".length.tab"
 clipped_output=$f".sam.clipped.fasta"
 touch $clipped_output
 clipped_rmasker_output=$f".sam.clipped.fasta.rmasker.fasta"
 #convert fastq to fasta using fastx_toolkit
 echo "Converting $f to fasta ($fasta)"
 fastq_to_fasta -i $f -o $fasta -Q33
 #run repeatmasker using specified library
 echo "Running repeatmasker on input fasta file"
 $repeat_masker -lib $lib_location -div $max_divergence -u -source -xsmall -norna -nolow -qq -pa $threads $fasta
 #parse repeatmasker output using ParseRM_simple.pl from Parsing-RepeatMasker-Outputs
 echo "Parsing repeatmasker output"
 $repeat_masker_parsing -RMout $rmasker_out -genfile $fasta
 #clip fragments outside of identified LINE
 echo "Getting soft-clipped sequences"
 $clip_rmasker --input $rmasker_parsed --fasta $fasta --output $clipped_rmasker_output

 #process R3 file - the same way like for R5

 R3_FILE=$FILENAME_PREFIX"R3"$FILENAME_SUFFIX
 echo $R3_FILE
 fasta_R3=$FILENAME_PREFIX"R3.fasta"
 rmasker_out_R3=$fasta_R3".out"
 rmasker_parsed_R3=$rmasker_out_R3".parsed1/"$rmasker_out_R3".length.tab"
 clipped_rmasker_output_R3=$R3_FILE".sam.clipped.fasta.rmasker.fasta"
 clipped_output_R3=$R3_FILE".sam.clipped.fasta"
 touch $clipped_output_R3
 #convert fastq to fasta using fastx_toolkit
 echo "Converting $R3_FILE to fasta ($fasta_R3)"
 fastq_to_fasta -i $R3_FILE -o $fasta_R3 -Q33
 #run repeatmasker using specified library
 echo "Running repeatmasker on input fasta file"
 $repeat_masker -lib $lib_location -div $max_divergence -u -source -xsmall -norna -nolow -qq -pa $threads $fasta_R3
 #parse repeatmasker output using ParseRM_simple.pl from Parsing-RepeatMasker-Outputs
 echo "Parsing repeatmasker output"
 $repeat_masker_parsing -RMout $rmasker_out_R3 -genfile $fasta_R3
 #clip fragments outside of identified LINE
 echo "Getting soft-clipped sequences"
 $clip_rmasker_R3 --input $rmasker_parsed_R3 --fasta $fasta_R3 --output $clipped_rmasker_output_R3
done

echo "DONE ALL"
