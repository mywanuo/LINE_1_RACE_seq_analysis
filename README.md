# RACE-seq

## Overview of the procedure

- Run [tailseeker3](https://github.com/hyeshik/tailseeker) on your data
- if required - perform additional demultiplexing based on primer sequences (we used modified [sabre](https://github.com/najoshi/sabre) for this purpose, available as a fork at <https://github.com/smaegol/sabre>)
- create samplesheet describing files to be analyzed, and they important features
- for LINE sequences (or other repetitive elements) identification run repeatmasker using `repeatmasker.sh` script (may take a long time)
- run `analyze_race_seq_flowcell2` script to get tails analysis done
- run `analyze_terminal_nucleotides.py` to analyze 3'-terminome
- perform additional analyzes in R using attached scripts

## Requirements

- tailseeker (<https://github.com/hyeshik/tailseeker>)
- repeatmasker (<http://www.repeatmasker.org/>) with the repbase libraries downloaded and installed (<http://www.girinst.org/server/RepBase/index.php>)
- Parsing-RepeatMasker-Outputs scripts downloaded (<https://github.com/4ureliek/Parsing-RepeatMasker-Outputs>)
- fastx toolkit (<http://hannonlab.cshl.edu/fastx_toolkit/>)
- bowtie2 (<http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>)
- Python 3
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

1. **Basecalling, deduplication, identification of A-tails**

    This step is done using [tailseeker software](https://github.com/hyeshik/tailseeker). Configuration of tailseeker run is provided in the tailseeker.yaml file which can be found in the flowcell1/flowcell2 folders. Before running the tailseeker modify the file to provide the proper path to flowcell data (dir variable). Tailseeker will generate R5 and R3 files, which are deduplicated (based on UMI sequences), and provide information about lengths of A-tails (and possible additions at the end of a tail). Tailseeker is run using command:

  ```
   tseek -j
  ```

2. **Demultiplexing (only for flowcell2)**

    Demultiplexing is done using [sabre](https://github.com/najoshi/sabre). The code of sabre was modified to include primer sequences in the output (apropriate [pull request](https://github.com/najoshi/sabre/pull/8) sent to the sabre developer).

    Script prepared for this purpose require the barcode files used for demultiplex (which can be found in `flowcell2/sabre_barcodes`) are located in the same folder as fastq files which will be demultiplexed (folder `fastq` in the output of tailseeker).

    To perform demultiplexing copy `demultiplex_sabre.sh` to the `fastq` folder and run:

  ```
   ./demultiplex_sabre.sh
  ```

3. **LINE1 sequences identification**

    To identify LINE1 sequences in demultiplexed reads [RepeatMasker](http://www.repeatmasker.org/) is used.

    Fastq sequences are first converted to fasta using `fastq_to_fasta` from fastx_toolkit. Then, RepeatMasker is run over the LINE1-specific database. Obtained hits are parsed using `parseRM_simple.pl` from [Parsing-RepeatMasker-Outputs](https://github.com/4ureliek/Parsing-RepeatMasker-Outputs) and analyzed using `identify_LINE_repeatmasker_softclip.py` to get information about location of LINE1 in the sequencing reads and about non-templated nucleotides (possible tails).

    Scripts `identify_LINE_repeatmasker_softclip.py` and `identify_LINE_repeatmasker_softclip_R3.py` must be copied to the `processing_out_sabre` folder. `repeatmasker.sh` should be run in the same folder.

  ```
   ./repeatmasker.sh
  ```

    This part of analysis can be time-consuming.

4. **Tails analysis**

    In the next step the actual analysis is done. For the LINE1 sequences the information about non-templated nucleotides (possible tails) is already obtained. For other sequences (GAPDH, reporter LINE1) it is retrieved by mapping using `bowtie2` with `--very-sensitive-local` option to get soft-clipping. Soft-clipped fragments are then retrieved using `get_softclipped_reads_from_sam.pl` script.

    Analysis is run on the all files with names ending with `_R5.fastq` in the folder specified with the `--inputdir` option. More specific selection can be done with the `--glob` option of the analysis script.

    For the all analyzed files a samplesheet is required, which contains all information regarding the samples, including experimental conditions, primer used, transcript, etc. This file must be prepared before the analysis. Example is located in the `flowcell2` folder (`samplesheet.csv`). Default path to the samplesheet is provided in the script, but the alternative one can be provided with the `--samplesheet` option.

    As the output (specified with the `--output` option) a tsv file is generated, containing tailing information for each sequence analysed, as well as additional data regarding the procedure.

    Before running the script it is required to customize settings (at the beginning of the script), like the path to `bowtie2`, number of threads it can use, names and locations of bowtie2 indexes.

    Analysis for the flowcell2 can be run using:

  ```
   ./analyze_race_seq_flowcell2.py --inputdir processing_out_sabre/ --output flowcell2_output.tsv
  ```

6. **Statistical analysis, plots**

  Further analysis is done using R scripts.

## How to get LINE-specific repeatmasker library

For the identification of LINE1 sequences we used RepeatMasker run with the library containing only human LINE1 sequences. It can be easily prepared using scripts available in the RepeatMasker distribution and HMMER. The detailed procedure is shown below:

- in the `util` subfolder in the RepeatMasker base directory run

  ```
    ./queryRepeatDatabase.pl -class "LINE" -species "Homo sapiens" > LINE_sequences.fasta
  ```

- extract LINEs names:

  ```
  cat LINE_sequences.fasta | grep "^>"  | cut -d' ' -f1 | sed -r 's/>(.*)/\1/' > LINE_names.txt
  ```

- extract hmms for LINEs from homo sapiens LINEs library (located in `Libraries/Dfam_2.0/homo_sapiens/`) using hmmfetch

  ```
    hmmfetch -f masklib.hmm lines_names.txt > LINEs.hmm
  ```

- press obtained hmm database

  ```
    hmmpress LINEs.hmm
  ```

- use this database in repeatmasker by specifying the exact path with the -d option
