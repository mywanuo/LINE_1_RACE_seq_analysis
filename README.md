# RACE-seq

## Overview of the procedure

* Run tailseeker3 (https://github.com/hyeshik/tailseeker) on your data
* if required - perform additional demultiplexing based on primer sequences (we used modified sabre (https://github.com/najoshi/sabre) for this purpose, available as a for at https://github.com/smaegol/sabre)
* create samplesheet describing files to be analyzed, and they important features
* for LINE sequences (or other repetitive elements) identification run repeatmasker using repeatmasker.sh script (may take a long time)
* run analyze_tails.py script to get tails analysis done
* perform additional analyzes in R using attached scripts

## Requirements

- tailseeker (https://github.com/hyeshik/tailseeker)
- repeatmasker (http://www.repeatmasker.org/) with the repbase libraries downloaded and installed (http://www.girinst.org/server/RepBase/index.php)
- Parsing-RepeatMasker-Outputs scripts downloaded (https://github.com/4ureliek/Parsing-RepeatMasker-Outputs)
- fastx toolkit (http://hannonlab.cshl.edu/fastx_toolkit/)
- bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- Python libraries installed: 
  - biopython
  - pandas
  - numpy
  - subprocess
- Perl libraries installed:
  - bioperl
  - Getopt
  - Term::ProgressBar
- R together with:
  - ggplot2
  - plyr
  - dplyr
  
## Detailed procedure

1. Run tailseeker using tailseeker.yaml file provided in the flowcell1/2 folders. Modify the tailseeker.yaml file (dir) to provide the proper path to flowcell data.
	Tailseeker will generate R5 and R3 files, which are deduplicated (based on UMI sequences), and provide information about lengths of A-tails (and possible additions at the end of a tail)

		tseek -j 
	 
2. Perform demultiplexing using `demultiplex_sabre.sh` (applies to flowcell2, which is multiplexed using primer sequences)
	* Script should be run inside the `fastq` folder created by tailseeker3
	* barcode files should be located inside the `fastq` folder
	
3. Run repeatmasker	to identify LINE1 sequences in RACE libraries prepeared using LINE1-specific primers. 
Fastq sequences are first converted to fasta. Then, RepeatMasker is run over the LINE1-specific database. Obtained hits are parsed using Parsing-RepeatMasker-Outputs and analyzed using `identify_LINE_repeatmasker_softclip.py` to get information about location of LINE1 in the saeuencing reads and about non-templated nucleotides (possible tails)
Scripts `identify_LINE_repeatmasker_softclip.py` and `identify_LINE_repeatmasker_softclip_R3.py` must be copied to the `processing_out_sabre` folder. `repeatmasker.sh` should be run in the same folder.

		./repeatmasker.sh

4. Run `analyze_race_seq_flowcell2.py`



## How to get LINE-specific repeatmasker library:

* in the `util` subfolder in the RepeatMasker base directory run

      ./queryRepeatDatabase.pl -class "LINE" -species "Homo sapiens" > LINE_sequences.fasta
  
* extract LINEs names:

      cat LINE_sequences.fasta | grep "^>"  | cut -d' ' -f1 | sed -r 's/>(.*)/\1/' > LINE_names.txt

* extract hmms for LINEs from homo sapiens LINEs library (located in `Libraries/Dfam_2.0/homo_sapiens/`) using hmmfetch

      hmmfetch -f masklib.hmm lines_names.txt > LINEs.hmm

* press obtained hmm database

      hmmpress LINEs.hmm
      
* use this database in repeatmasker by specifying the exact path with the -d option  
