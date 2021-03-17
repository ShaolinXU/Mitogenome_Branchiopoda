#!/bin/bash

seqkit_threads=4
spades_threads=4
number_iteration=4
minimum_id=95
max_gap=1
mode=dna
how_long=15000

function usage () {
   cat <<EOF
Author: Shaolin-Xu
Date: 2018-04-12
This script need metaspades blatq seqkit to do the analysis
 have to be in the <<<<<  PATH  >>>>>
-------- blatq to seeds -> metaspades -> new_ref -> compared to the seeds -> iteration.  --------
Usage: $scriptname [-1] [-2] [-r] [-t] [-i] [-T] [-I] [-G] [-l]
   -1       The forward pair-end sequence file
   -2       The reverse pair-end sequence file
   -r       The reference fasta file, only one sequence allowed
   -t (4)   The number of threads which seqkit use to pick out seqeunce
   -i (4)   The number of iteration 
   -T (4)   The number of threads for SPAdes
   -I (95)  The minium identity (in percent) between reads and reference for blatq
   -G (1)	The max-gap allowed in alignment
   -l    	The predefined length to stop that specific contigs' assembly
EOF
   exit 0
}
# As I use virtual environment, so I have to source the env first, if you don't need, comment the line below
#source activate bioinfo

# check if seqkit was installed
if [ -z `which seqkit` ] ; then
	echo "\nInstall seqkit first"
	exit 1
fi
# check if blatq was installed
if [ -z `which blatq` ] ; then
	echo "\nInstall blatq first"
	exit 1
fi
# check if spades was installed
if [ -z `which metaspades.py` ] ; then
	echo "\nInstall metaspades.py first"
	exit 1
fi
# check if usearch was installed
#if [ -z `which usearch` ] ; then
#	echo "\nInstall usearch first"
#	exit 1
#fi

while getopts :1:2:r:t:i:T:I:G:l: opt;do
    case $opt in
        1)  pair1=$OPTARG  ;;
        2)  pair2=$OPTARG  ;;
        r)  reference=$OPTARG  ;;
        t)  seqkit_threads=$OPTARG  ;;
        i)  number_iteration=$OPTARG  ;;
        T)  spades_threads=$OPTARG ;;
		I)	minimum_id=$OPTARG ;;
		G)	max_gap=$OPTARG ;;
		l)	how_long=$OPTARG ;;
        \?) usage
            exit 1  ;;
        h)	usage
      			exit 0	;;
    esac
done

if [ -z ${pair1+x} ]; then echo "pair1 (-1) is unset" ; exit 1; fi
if [ -z ${pair2+x} ]; then echo "pair2 (-2) is unset" ; exit 1; fi
if [ -z ${reference+x} ]; then echo "reference (-r) is unset" ; exit 1; fi
if [ -z ${seqkit_threads+x} ]; then echo "seqkit_threads (-t) is unset" ; exit 1; fi
if [ -z ${number_iteration+x} ]; then echo "number_iteration (-i) is unset" ; exit 1; fi
if [ -z ${spades_threads+x} ]; then echo "spades_threads (-T) is unset" ; exit 1; fi
if [ -z ${minimum_id+x} ]; then echo "minimum_id (-I) is unset" ; exit 1; fi
if [ -z ${max_gap+x} ]; then echo "max_gap (-G) is unset" ; exit 1; fi
# if [ ${number_iteration} -gt 4 ]; then echo "number_iteration can not greater than 4" ; exit 1; fi

index_log=index.log

[[ -e $index_log ]] && rm -f index.log

# put all ref in new_reference.fa 
seqkit --quiet seq $reference > new_reference.fa

touch pair12_backup.id

## create a file to preserve the total assembled contigs length within a metaspades execution
touch sum_assembled_contigs.txt

## reference history files
mkdir reference_files

## read id history found by blat
mkdir blat_files

## keep change of fq files
mkdir fq_files

START0=$(date +%s)

for (( i = 0; i < $number_iteration; i++ )); do # use only 2 file to assembly
	printf "\n====================  <<  The Analysis run_${i} started  >>  ======================\n==> $(date) ---- \n"
	START=$(date +%s) ## record the starting time of blatq
	#### starting blatq analysis
	pair1_basename=`basename $pair1`
	pair2_basename=`basename $pair2`

	blatq new_reference.fa $pair1 -mask=lower -q=dna -t=dna -minIdentity=$minimum_id -maxGap=$max_gap -fastMap -out=blast8 ${pair1_basename%.*}.blatq &
	blatq new_reference.fa $pair2 -mask=lower -q=dna -t=dna -minIdentity=$minimum_id -maxGap=$max_gap -fastMap -out=blast8 ${pair2_basename%.*}.blatq &

	#### check if the blatq analysis is finished or not

	wait 

	END=$(date +%s)
	DIFF=$(( $END - $START ))
	printf "\n==> $(date) ---- blatq finished, using $DIFF seconds\n"

	#### get the id of hited sequence

	START=$(date +%s) ## record the starting time of seqkit

	cat ${pair1_basename%.*}.blatq ${pair2_basename%.*}.blatq | awk -v minIdentity="$minimum_id" '$3>=minIdentity{print $1}' | sort | uniq > pair12.id
	cat pair12.id pair12_backup.id | sort | uniq > pair12_2.id
	mv pair12_2.id pair12.id
	cp -f pair12.id pair12_backup.id

	## keep the change of id file in the blat_files dir
	cp pair12_backup.id blat_files/pair12_run${i}_id.txt

	seqkit grep -i -f pair12.id $pair1 -o forward_pair.fq --threads $seqkit_threads &
	seqkit grep -i -f pair12.id $pair2 -o reverse_pair.fq --threads $seqkit_threads &

	## keep the changes of fq file
	if [[ !${i} -eq 0 ]]; then
		cp forward_pair.fq fq_files/forward_run${i}.fq
		cp reverse_pair.fq fq_files/reverse_run${i}.fq
	fi

	#### check if seqkit is finished or not

	wait

	END=$(date +%s)
	DIFF=$(( $END - $START ))
	printf "\n==> $(date) ---- seqkit finished, using $DIFF seconds\n"


	#### use spades to assembly the result
	START=$(date +%s) ## record the starting time of spades
	
	metaspades.py -1 forward_pair.fq -2 reverse_pair.fq -o metaSPAdes_run_${i} -t ${spades_threads} -k 47,67,87,107,127 --phred-offset 33 >> spades.log

	blatq metaSPAdes_run_${i}/scaffolds.fasta new_reference.fa -q=dna -t=dna -minIdentity=97 -maxGap=$max_gap -out=blast8 new_reference.blatq 


	##########################################
	## check if the run should be stoped-step01
	if [[ ! -s new_reference.blatq ]]; then
		printf "\nAnalysis filed at contig blating process of run_${i} \n\n"
		exit 1
	fi

	cat new_reference.blatq | awk '$4>=100&&$3>=97{print $2}' > better_contig.id # only get new contigs overlap with old reference more than 100 bp, at the same time, those contigs also need to have a identity score no less than 97% 
	

	if [[ ! -s better_contig.id ]]; then
		printf "\nAnalysis filed at contig blating filtering process of run_${i} \n\n"
		exit 1
	fi
	##########################################

	cp new_reference.fa reference_files/run_${i}_reference.fa

	seqkit grep -i -f better_contig.id metaSPAdes_run_${i}/scaffolds.fasta | seqkit --quiet seq -m 1000 -g > new_reference.fa1
	seqkit --quiet seq -m ${how_long} new_reference.fa1 -g >> Long_enough.fa
	seqkit --quiet seq -M ${how_long} new_reference.fa1 -g > new_reference.fa


	##########################################
	## check if the run should be stoped-step02

	if [[ !${i} -eq "0" ]]; then
		current_run_sum_contigs=`seqkit stat -T new_reference.fa | csvtk -t cut -f 5 | tail -n 1`
		if grep -Fxq "$current_run_sum_contigs" sum_assembled_contigs.txt; then
			printf "\nAnalysis filed because of assembled total length was exist at run_${i} \n which means this assembled were very likely from the same set of reads\n\n"
			exit 1
		fi
	fi

	seqkit stat -T new_reference.fa | csvtk -t cut -f 5 | tail -n 1 >> sum_assembled_contigs.txt


	printf "==> $(date) ---- This is the ${i}_th iteration \n " >> index.log
	read_num_pair=`wc -l pair12.id | awk '{print $1}'`
	printf "==> $(date) ---- The reads in pair12.id is ${read_num_pair} ! \n " >> index.log

	END=$(date +%s)
	DIFF=$(( $END - $START ))
	printf "\n==> $(date) ---- SPAdes finished, using $DIFF seconds\n\n"


	done

DIFF=$(( $END - $START0 ))
printf "\n==> $(date) ---- The whole analysis took $DIFF seconds\n"
