#!/bin/bash


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




#specify location of sabre:
sabre_script=`which sabre`
sabre_max_mismatches='2' # max number of mismatches allowed in primer/barcode

#get input dir from command line
input_dir=$1

#output dir
out_dir='processing_out_sabre'
#nmae of barcode files

mkdir $out_dir


for R5_file in `find $input_dir -name '*R5.fastq'`
do
	echo "processing $R5_file"
	FILENAME_PREFIX=`expr match "$R5_file" '\(.*\)R5'`
	PREFIX_BASENAME=`expr match "$FILENAME_PREFIX" '.*\/\(.*\)'`
	R3_file=$FILENAME_PREFIX"R3.fastq"
	sabre_barcodes_file=$PREFIX_BASENAME"barcodes.txt"
	sabre_output=$PREFIX_BASENAME"sabre_out.txt" # this file will contain all demultiplex statistics
	output_R5_untrimmed=$out_dir'/'$PREFIX_BASENAME"R5_untrimmed.fastq" #location of untrimmed R5 reads
	output_R3_untrimmed=$out_dir'/'$PREFIX_BASENAME"R3_untrimmed.fastq"	#location of untrimmed R3 reads
	
	#perform demultiplexing
	$sabre_script pe -f $R5_file -r $R3_file -b $sabre_barcodes_file -u $output_R5_untrimmed -w $output_R3_untrimmed -m $sabre_max_mismatches > $sabre_output
	
done
