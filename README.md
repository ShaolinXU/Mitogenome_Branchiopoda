# Mitogenome_Branchiopoda
scripts used for the mitogenome paper of Cladoceran by Shao-Lin Xu

## Script
1. `code/iterative_metaspades_fishing.sh`: extend seed sequence by  iterative target assembly
    1. `code/IMF_readme.md`: readme file for `iterative_metaspades_fishing.sh`
2. `code/gcskwew.py`: calculate compositional statistic of dna sequence (input as fasta format)
    1. this script was modified from [SKewIT](https://github.com/jenniferlu717/SkewIT), all credit goes to the original author <sup>[1](#myfootnote1)</sup>
3. `code/zorro.py`: filter ambiguious align region
    1. This script was modified from Schwentner 2018 <sup>[2](#myfootnote2)</sup>

## Data
1. `data/Appdendix.zip`: all appdendix of this paper
2. new 'matrix-7' data:
    1. `data/OG_168.fa`: amino acid sequence of 168 orthologous genes in the new 'matrix-7' in the paper
    2. `data/OG_168.partition`: the partition file of `OG_168.fasta`
3. `data/MG15.gb`: genebank file of 15 newly assembled mitogenomes of Cladocera


# Reference
<a name="myfootnote1">1</a>: Lu, Jennifer, and Steven L. Salzberg. "SkewIT: The Skew Index Test for large-scale GC Skew analysis of bacterial genomes." PLoS computational biology 16.12 (2020): e1008439.

<a name="myfootnote2">2</a>: Schwentner, M., S. Richter, D. C. Rogers, and G. Giribet. 2018. Tetraconatan phylogeny with special focus on Malacostraca and Branchiopoda: highlighting the strength of taxon-specific matrices in phylogenomics. Proceedings of the Royal Society B: Biological Sciences 285:20181524.
