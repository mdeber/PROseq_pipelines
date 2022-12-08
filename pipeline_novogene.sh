#!/usr/bin/env bash

# PRO-seq alignment (for Novogene data releases)
# Mike DeBerardine
# Updated October 20 2022
#
# This script is for aligning PRO-seq data, going from unprocessed fastq files to 
# UMI-deduplicated bam files.
#
# Code is written to find fastq files and user-submitted sample names as they would be
# organized in a data release from Novogene: in a folder named "raw_data" with subfolders 
# for each sample that contain files ending in _1.fq/_2.fq or _1.fq.gz/_2.fq.gz
#
# requirements:
#   pigz/unpigz if fastq are compressed (or can be changed to gzip/gunzip)
#   fastqc (can be safely skipped; mostly redundant with fastp's qc)
#   fastp
#   bowtie2
#   umi_tools
#   samtools
#
# This script to be run within the parent directory of the pipeline directory
# that contains this file (so will be run as 'bash pipeline/pipeline.sh')

NTHREADS=64;

# These arguments are passed to fastp. UMI defaults are based on the REV3 and REV5 adapters, 
# as used in the qPRO-seq protocol (UMI_LEN=6 and UMI_SKIP=1). The adapter sequences themselves 
# are Tru-seq small RNA (although reverse complemented in PRO-seq)
UMI_LEN=6;
UMI_SKIP=1;
ADAPTER1="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC";
ADAPTER2="GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT";

# bowtie2 genomes
GENOME="/home/mdd238/genomes/bowtie2/GRCh38_noalt_decoy_as/GRCh38_noalt_decoy_as";
RDNA="/home/mdd238/genomes/rDNA_repeats/rDNAhuman_bt2/rDNAhuman";
MAPQ=30;


# ----
# Notes on alternative UMI structures:
# 
# If a UMI is present in only 1 read, the correct arguments for fastp are straightforward.
# See the fastp manual for that (and change the --umi_loc argument).
# 
# If two UMIs are used, but they are different lengths, or if the length of non-UMI bases
# between the UMI and the insert are different for each adapter, then it is very difficult
# (or maybe impossible) for fastp to correctly extract the two UMIs from each read. A number
# of arguments may seem possible, but the issue is whether fastp will correctly perform overlap
# analysis and remove all adapter bases from the opposing read. You cannot, for instance, extract
# a base from the Illumina adapter, because it's already been trimmed in pre-processing from one
# of the reads. The --trim_front1 or --trim_front2 arguments do not (from what I can tell in tests) 
# remove that base from the opposing read. You also cannot do adapter trimming in 2 steps, because
# this also prevents overlap analysis necessary for identifying those bases in the opposing read.
# As far as I know, you will either have to leave an adapter base on the opposing read (for reads 
# in which the two pairs overlap entirely, which is common in PRO-seq), or remove one of the genomic 
# bases from every read.
# 
# A simple exception is if the UMIs are different lengths, but if this is compensated by invariant
# bases that can be extracted as part of the UMI.
# ----

echo "creating output directories if they don't already exist..."
[ -e logs ] || mkdir logs;
[ -e logs/fastqc ] || mkdir logs/fastqc;
[ -e logs/fastp ] || mkdir logs/fastp;
[ -e logs/bowtie2_rRNA ] || mkdir logs/bowtie2_rRNA;
[ -e logs/bowtie2 ] || mkdir logs/bowtie2;
[ -e logs/bam_dedup ] || mkdir logs/bam_dedup;
[ -e trimmedFastq ] || mkdir trimmedFastq;
[ -e bam ] || mkdir bam;
[ -e bam_dedup ] || mkdir bam_dedup;

echo "unzipping any compressed fastq files...";
if [ $(find raw_data -type f -name "*fq.gz" | wc -l) -gt 0 ];
then
    unpigz -p $NTHREADS $(find raw_data -type f -name "*fq.gz");
fi;

echo "running fastqc...";
fastqc \
    $(find raw_data -type f -name "*fq") \
    -o logs/fastqc \
    --threads $NTHREADS \
    --quiet;


echo "trimming adapters & depleting rRNA reads for sample:";
for sample_name in $(ls raw_data | grep -v Undetermined | grep -v Readme);
do
    fq_basename=$(ls raw_data/${sample_name}/*fq | sed "s/_..fq//" | uniq);

    echo "  ${sample_name}...";
	
	(fastp \
		--in1 ${fq_basename}_1.fq \
		--in2 ${fq_basename}_2.fq \
		--stdout \
		-c \
		--overlap_len_require=15 \
		--adapter_sequence $ADAPTER1 \
		--adapter_sequence_r2 $ADAPTER2 \
		--umi \
		--umi_loc=per_read \
		--umi_len=${UMI_LEN} \
		--umi_skip=${UMI_SKIP} \
		--umi_prefix="UMI" \
		--html logs/fastp/fastp_${sample_name}.html \
		--json logs/fastp/fastp_${sample_name}.json \
		-w $(echo ${NTHREADS}/3 | bc) \
		2> logs/fastp/fastp_${sample_name}.log) |
	(bowtie2 \
		--fast-local \
		--un-conc trimmedFastq/${sample_name}.fq \
		--interleaved - \
		-x $RDNA \
		--threads $(echo ${NTHREADS}/5*3 | bc) \
		2> logs/bowtie2_rRNA/bt2_rRNA_${sample_name}.log) \
		> /dev/null;
    
    mv trimmedFastq/${sample_name}.1.fq trimmedFastq/${sample_name}_1.fq;
    mv trimmedFastq/${sample_name}.2.fq trimmedFastq/${sample_name}_2.fq;
done

echo "finished trimming & rRNA depletion. compressing original/untrimmed fastq files...";
pigz -p $NTHREADS $(find raw_data -type f -name "*fq");

echo "aligning to genome for sample:"
for sample_name in $(ls -1 trimmedFastq | sed "s/_..fq//" | uniq)
do
    echo "  ${sample_name}...";

	# ----
	# Note that local alignment mode is being used. While "end-to-end" might seem to be
	# more consistent with basepair resolution, it's safer in principle to drop non-genomic
	# bases from edges of aligned sequences, as they could be from a mistake in adapter
	# trimming (like failed identification of the reverse complement of a UMI due to an
	# insufficient overlap) or another adapter/trimming error. The aligned (genome-matching)
	# bases are the best bet for determining the insert with high precision.
	# ----
    
    (bowtie2 \
		--very-sensitive-local \
		--threads $(echo ${NTHREADS}/3*2 | bc) \
		-x $GENOME \
		-1 trimmedFastq/${sample_name}_1.fq \
		-2 trimmedFastq/${sample_name}_2.fq \
		--no-unal \
		--no-mixed \
		--no-discordant \
		2> logs/bowtie2/bt2_${sample_name}.log) |
	samtools view -b -q ${MAPQ} |
	samtools sort -@ $(echo ${NTHREADS}/3 | bc) -o bam/${sample_name}.bam;
    samtools index -@ $NTHREADS bam/${sample_name}.bam;
done

echo "finished aligning to genome. compressing trimmed fastq files...";
pigz -p $NTHREADS $(find trimmedFastq -type f -name "*fq");

echo "removing UMI duplicate reads for sample:";
for sample_name in $(ls bam | grep "bam$" | sed "s/.bam//");
do
    echo "  ${sample_name}...";
    umi_tools dedup \
		-I bam/${sample_name}.bam \
		--paired \
		--umi-separator=":UMI_" \
		--log logs/bam_dedup/dedup_${sample_name}.log \
		-S bam_dedup/${sample_name}_dedup.bam;
    samtools index -@ $NTHREADS bam_dedup/${sample_name}_dedup.bam;
done

echo "finished removing UMI duplicate reads";

# ----
# Below is mostly a copy of Julius Judd's code to generate a summary statistics table,
# with various modifications to suit this pipeline. Namely, this pipeline aligns to rRNA 
# separately and exports a log for that, and so the number of rRNA reads is directly measured. 
# Additionally, the percent rRNA reads is calculated based on the total trimmed/filtered reads 
# (after fastp), while Julius's code uses the original (unfiltered) number of read pairs. 
# 
# This pipeline also aligns to a single genome and does not address the spike-in. This is because
# I align to a combined genome for competitive alignment of spike-in vs. non-spike-in reads, 
# and I separate out the chromosomes later on in R.
# ----

echo "creating the alignment summary infoTable.tsv document";

if [ ! -s logs/infoTable.tsv ];
then
    touch logs/infoTable.tsv;

    echo -e Name'\t'\
			RawReads'\t'\
			TrimmedReads'\t'\
			%PassedTrimming'\t'\
			rRNAreads'\t'\
			%rRNA'\t'\
			nonRRNA'\t'\
			bowtieUniqConcordant'\t'\
			bowtieMulti'\t'\
			bowtieUnal'\t'\
			bowtieOverallMap%'\t'\
			bowtieUniqConcordant%'\t'\
			bowtieMulti%'\t'\
			bowtieUnal%'\t'\
			UniqueMapped'\t'\
			UniqueMappedNondup'\t'\
			%PCRdups >> logs/infoTable.tsv;

    for SAMPLE in $(ls bam | grep "bam$" | sed "s/.bam//");
    do
		NAME=${SAMPLE};
		RAW_READS=$(cat logs/fastp/fastp_${SAMPLE}.log |
			grep "total reads:" | head -n 1 |
			awk '{print $3}');
		TRIMMED_READS=$(cat logs/fastp/fastp_${SAMPLE}.log |
			grep "total reads:" | tail -n 1 |
			awk '{print $3}');	
		PER_PASSED_TRIM=$(echo "("${TRIMMED_READS}"/"${RAW_READS}")*100" | bc -l)%;
		
		NON_RRNA=$(cat logs/bowtie2/bt2_${SAMPLE}.log |
			grep "reads; of these:$" |
			awk '{print $1}');
		RRNA=$(echo ${TRIMMED_READS}"-"${NON_RRNA} | bc );
		PER_RRNA=$(echo ${RRNA}"/"${TRIMMED_READS}"*100" | bc -l)%;
		
		B_CONC=$(cat logs/bowtie2/bt2_${SAMPLE}.log |
			grep "aligned concordantly exactly 1 time$" |
			awk '{print $1}');
		B_MULTI=$(cat logs/bowtie2/bt2_${SAMPLE}.log |
			grep "aligned concordantly >1 times$" |
			awk '{print $1}');
		B_UNAL=$(cat logs/bowtie2/bt2_${SAMPLE}.log |
			grep "aligned concordantly 0 times$" |
			awk '{print $1}');
		B_OAP=$(cat logs/bowtie2/bt2_${SAMPLE}.log |
			grep "overall alignment rate$" |
			awk '{print $1}');
		B_CONC_PER=$(echo ${B_CONC}"/"${NON_RRNA}"*100" | bc -l)%;
		B_MULTI_PER=$(echo ${B_MULTI}"/"${NON_RRNA}"*100" | bc -l)%;
		B_UNAL_PER=$(echo ${B_UNAL}"/"${NON_RRNA}"*100" | bc -l)%;

		UNIQ_MAPPED=$(cat logs/bam_dedup/dedup_${SAMPLE}.log |
			grep "Input Reads:" | awk '{print $10}');
		UNIQ_MAPPED_DEDUP=$(cat logs/bam_dedup/dedup_${SAMPLE}.log |
			grep "Number of reads out:" | awk '{print $8}');
		PER_DUPS=$(echo "(1-"${UNIQ_MAPPED_DEDUP}"/"${UNIQ_MAPPED}")*100" | bc -l)%;

		echo -e $NAME'\t'\
				$RAW_READS'\t'\
				$TRIMMED_READS'\t'\
				$PER_PASSED_TRIM'\t'\
				$RRNA'\t'\
				$PER_RRNA'\t'\
				$NON_RRNA'\t'\
				$B_CONC'\t'\
				$B_MULTI'\t'\
				$B_UNAL'\t'\
				$B_OAP'\t'\
				$B_CONC_PER'\t'\
				$B_MULTI_PER'\t'\
				$B_UNAL_PER'\t'\
				$UNIQ_MAPPED'\t'\
				$UNIQ_MAPPED_DEDUP'\t'\
				$PER_DUPS >> logs/infoTable.tsv;
    done;
fi;

echo "script finished.";
exit 0;
