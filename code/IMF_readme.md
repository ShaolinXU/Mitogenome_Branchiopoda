# Installation
1. blatq
    1. [GitHub - calacademy-research/BLATq: blatq is a modified version of BLAT—The BLAST-Like Alignment Tool to allow direct use of fastq files](https://github.com/calacademy-research/BLATq)
2. spades
    1. install with conda: `conda install spades`
3. seqkit
    1. install with conda: `conda install seqkit`

# Usage
## Quality control
Remove adapter and low quality reads in your NGS datasets, using tool like [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
## Seed sequence
### Source
1. closely related species
2. contigs from pre-assembly
3. Sanger sequence from PCR, such as the barcoding region in COX1


## iterative assembly

```bash
iterative_price_fishing.sh \
-1 forward.fq \ # forward reads
-2 reverse.fq \ #reverse reads
 -r seed_sequence.fa \ #seed sequence to be extended
 -i 20 \ # number of iteration to try
 -l 20000 \ # the length threshold
 -t 1 \ # number of cpu used for seqkit extraction step
 -T 2 \ # number of cpu used by metaspades.py
 1>logfile.txt 2>&1 # save stdout and possible error message to logfile.txt
```

## Parameters details
- values within parentheses are the default value
- `-1` and `-2` must used at the same time

```bash
-1       The forward pair-end sequence file
-2       The reverse pair-end sequence file
-r       The reference fasta file, only one sequence allowed
-t (4)   The number of threads which seqkit use to pick out seqeunce
-i (4)   The number of iteration 
-T (4)   The number of threads for SPAdes
-I (95)  The minium identity (in percent) between reads and reference for blatq
-G (1)	The max-gap allowed in alignment
-l    	The predefined length to stop that specific contigs' extension
```

## Tips
1. `-G` and `-I` can be played with for better extension results
