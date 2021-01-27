### Characterization of sequence features in IGHV3-53/66 antibodies that bind to SARS-CoV-2 RBD
# Dependencies #
* [Python](https://www.python.org/) version 2.7
* [R](https://www.r-project.org/) version 3.6
* [PEAR](https://github.com/tseemann/PEAR)

## Installation ##

Easiest way to install PEAR is through [Anaconda](https://anaconda.org/bioconda/pear)

And the following to install the needed package:
```
conda install -c bioconda pear

```
## INPUT FILES
* [./data/VH3-53_RBD_antibodies.tsv](./data/VH3-53_RBD_antibodies.tsv): List of IGHV3-53/3-66 RBD antibodies from literatures
* [./data/IGKV1-9_CDRH3.txt](./data/IGKV1-9_CDRH3.txt): List of CDR H3 from IGHV3-53/3-66 RBD antibodies that pair with IGKV1-9

## Analyzing IGHV3-53/3-66 RBD antibodies from literature ##
1. categorize antibodies based on CDR H3 length and light chain germline
    - run python script ```script/summarize_LC_vs_CDRH3_length.py```
    - run R script ```script/plot_LC_vs_CDRH3_length.R``` for plotting

2. analyzing SHM for antibodies with different CDR H3 length
    - run python script ```script/format_SHM_table.py```
    - run R script ```script/plot_SHM_heatmap.R``` for plotting

## Analyzing yeast display screening NGS data ##
1. merge overlapping paired-end reads
    - pear -f [READ1 FASTQ FILE] -r [READ2 FASTQ FILE] -o fastq/S10
    
2. run analysis for the assembled files
    - run python script ```script/fastq2count.py``` to count the frequence of each variant of CDRH3 loop
    - run python script ```script/merge&calculate.py``` to calculate the relative enrichment of each variant

3. plotting results
    - run R script ```script/plot_yeast_display.R```
