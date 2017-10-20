# RACE-seq

## Procedure

* Run tailseeker3 (https://github.com/hyeshik/tailseeker) on your data
* if required - perform additional demultiplexing based on primer sequences (we used modified sabre (https://github.com/najoshi/sabre) for this purpose)
* create samplesheet describing files to be analyzed, and they important features
* for LINE sequences (or other repetitive elements) identification run repeatmasker using repeatmasker.sh script
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
