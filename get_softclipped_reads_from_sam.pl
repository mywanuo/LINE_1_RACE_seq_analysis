#!/usr/bin/perl

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



#load required libraries
use strict;
#use warnings;
use Getopt::Long;
use File::Copy;
use File::Basename;
use Bio::SeqIO;
use Bio::Seq::Quality;
use Term::ProgressBar 2.00;

my $current_dir = `pwd`;

print $current_dir; #print current dir

my $input_SAM = '';
my $output_file = '';

#Get options from command line
GetOptions ('input=s' => \$input_SAM, 'output=s' => \$output_file);


if ($input_SAM eq '') {
	die "no input file specified\n"
}
if ($output_file eq '') {
	die "no output file specified\n"
}

#process input files:
get_softclipped_reads_with_ref($input_SAM,$output_file);




### SUBROUTINES ###

sub get_softclipped_reads_with_ref {
#get mapping positions of mapped reads - find where transcript ends in RACE sequencing
#sam - mapping of reads with tails clipped

	my $input_sam_file = shift; #sam file with mappings
	my $output_fastq_file = shift; #sam file with mappings
	my $input_fastq_file = shift;

	print "sam: $input_sam_file\n\n";

	#parameters read from SAM file
	my $query_name;
	my $bit_flag;
	my $ref_name;
	my $pos;
	my $mapq;
	my $cigar;
	my $rnext;
	my $pnext;
	my $tlen;
	my $seq;
	my $qual;


	my %reads_descriptions = ();

	my %positions = ();
	my $positions_file = $input_sam_file."_positions.txt";
	my %positions_tails = ();
	my $positions_file_tails = $input_sam_file."_positions_tails.txt";
	my %positions_notails = ();
	my $positions_file_notails = $input_sam_file."_positions_notails.txt";

	if ($input_fastq_file) {
		#specify input file:
		my $sequences_in = Bio::SeqIO->new(
								-file   => "<$input_fastq_file",
								-format => 'fastq',
								);
		while (my $data = $sequences_in->next_dataset) {
				my $read =  Bio::Seq::Quality->new(%$data); #read sequence with quality info
				my $read_seq = $read->seq; #get actual sequence
				my $read_qual = $read->qual;
				my $read_id = $read->id;
				my $readdesc = $read->desc;
				$reads_descriptions{$read_id}=$readdesc;
		}
	}

	#create fasta for output files
	my $sequences_out = Bio::SeqIO->new(
							-file   => ">$output_fastq_file",
							-format => 'fasta',
							);

	open (INPUT,"<$input_sam_file"); #open SAM file for reading
	my $i=0;

	my $no_sam_lines = `wc -l $input_sam_file`;
	$no_sam_lines =~ /(\d+)\s+/;
	$no_sam_lines = $1;

	my $progress = Term::ProgressBar->new({name => 'Reading SAM file', count => $no_sam_lines, ETA   => 'linear', });
	my $next_update=0;
	while (my $line = <INPUT>) {
	#read all lines of SAM file
		$i++;
		if (! ($line =~ /^\@/)) { #ignore header lines
			#read all fileds in SAM:
			($query_name,$bit_flag,$ref_name,$pos,$mapq,($cigar),$rnext,$pnext,$tlen,$seq,$qual) = ($line =~ /^([!-?A-~]{1,255})\s+(\d+)\s+(\*|[!-()+-<>-~][!-~]*)\s+(\d+)\s+(\d+)\s+(\*|[0-9MIDNSHPX]+)\s+(\*|=|.+)\s+(\d+)\s+(-{0,1}\d+)\s+(\*|[A-Za-z=.]+)\s+([!-~]+)\s+/);

			if ($ref_name ne '*') { #if read was mapped - get position of mapping
				my $seq_length = length($seq); #get length of mapped read
				if ($cigar=="*") {
				#if no clipping occured:
					my $readdesc = $reads_descriptions{$query_name};
					$readdesc.="\tclip5: \tclip3: \tpos: -1\tref: -1";
					my $seq_obj_out = Bio::Seq->new(
									-seq => $seq,
									-id => $query_name,
									-desc => $readdesc);
					$sequences_out->write_seq($seq_obj_out);
				}
				else {
				#get sequence parts from CIGAR:
					my @cigar = ($cigar =~ m/(\d+\D)/g);
					my $nocigars = scalar @cigar; # get number of cigar operations
					my $cum_length = 0;
					my $clip5;
					my $clip3;
					my $matched_read = $seq;
					my $cig;
					my $cig_type;
					my $cig_length;
					my $match_start=0;
					my $match_end=$seq_length;
					if ($nocigars>1) {
						#if there is more than 1 cigar operation - look for soft clipping
						for (my $z=0;$z<$nocigars;$z++){
							$cig = $cigar[$z];
							my ($cig_length,$cig_type) = ($cig =~ m/(\d+)(\D)/);

							if ($z==0) {
								#get first CIGAR operation - if S (softclipping) - output as 5'-clipping
								if ($cig_type eq "S") {
									$clip5 = substr($seq,0,$cig_length);
									$match_start=$cig_length;
									$match_end-=$cig_length;
								}
							}
							if ($z==$nocigars-1) {
								#get last CIGAR operation - if S (softclipping) - output as 3'-clipping
								if ($cig_type eq "S") {
									$clip3 = substr($seq,$cum_length,$cig_length);
									$match_end-=$cig_length;
								}
							}
							$cum_length+=$cig_length;
						}
						$matched_read=substr($seq,$match_start,$match_end); #get matched read (without clipped parts)

					}
					my $mapping_pos = $pos + $match_end - 1; #get mapping position
					($positions{$mapping_pos}) and ($positions{$mapping_pos}++) or ($positions{$mapping_pos}=1);
					if ($clip3) {
						($positions_tails{$mapping_pos}) and ($positions_tails{$mapping_pos}++) or ($positions_tails{$mapping_pos}=1);
					}
					else {
						($positions_notails{$mapping_pos}) and ($positions_notails{$mapping_pos}++) or ($positions_notails{$mapping_pos}=1);
					}
					my $readdesc = $reads_descriptions{$query_name}; #create read description with clipping info
					$readdesc.="\tclip5: $clip5\tclip3: $clip3\tpos: $mapping_pos\tref: $ref_name";
					#my $seq_obj_out = Bio::Seq::Quality->new(
					my $seq_obj_out = Bio::Seq->new(
										-seq => $matched_read,
										-id => $query_name,
										-desc => $readdesc);
					$sequences_out->write_seq($seq_obj_out);
				}
			} #end if
			else {
				#read was unmapped
				my $readdesc = $reads_descriptions{$query_name};
				$readdesc.="\tclip5: \tclip3: \tpos: -1\tref: -1";
				my $seq_obj_out = Bio::Seq->new(
									-seq => $seq,
									-id => $query_name,
									-desc => $readdesc);
				$sequences_out->write_seq($seq_obj_out);

			}
		} #end if


		#manage progress bar
		$next_update = $progress->update($i)
		if $i >= $next_update;
	} #end while

	open (POS,">$positions_file");
	foreach my $pos (sort keys %positions) {
		print POS "$pos\t$positions{$pos}\n";
	}
	close(POS);

	open (POS,">$positions_file_tails");
	foreach my $pos (sort keys %positions_tails) {
		print POS "$pos\t$positions_tails{$pos}\n";
	}
	close(POS);

	open (POS,">$positions_file_notails");
	foreach my $pos (sort keys %positions_notails) {
		print POS "$pos\t$positions_notails{$pos}\n";
	}
	close(POS);

} # end sub


exit;
