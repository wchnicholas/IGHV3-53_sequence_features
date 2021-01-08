### Characterization of sequence features in IGHV3-53/66 antibodies that bind to SARS-CoV-2 RBD
# Dependencies #
* [Python](https://www.python.org/) version 2.7
* [PEAR](https://github.com/tseemann/PEAR)

## Installation ##

Easiest way to install PEAR is through [Anaconda](https://anaconda.org/bioconda/pear)

And the following to install the needed package:
```
conda install -c bioconda pear

```
## INPUT FILES

## Steps ##
1. merge overlapping paired-end reads
    - pear -f Sample10_TAGCTTAT_L001_R1_001.fastq.gz -r Sample10_TAGCTTAT_L001_R2_001.fastq.gz -o fastq/S10
    
2. run analysis for the assembled files
    - go into ```/script/``` folder
    -run python script ```/script/fastq2count.py``` to count the frequence of each variant of CDRH3 loop
    -run python script ```/script/merge&calculate.py``` to calculate the relative enrichment of each variant
