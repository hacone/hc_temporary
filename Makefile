
### Makefile records only the types (type-theoretic) of procedures.

## TODO: check sanity at each step as possible
## TODO; under construction...

ACC=SRR3189743
#ACC=SRR1997411

all: aligned_sam aligned_fastq

$(ACC)_1.fastq:
	echo "pass"
	# fastq-dump for ACC in sratoolkit, too slow?

$(ACC)_2.fastq:
	echo "pass"
	# fastq-dump for ACC in sratoolkit, too slow?

# how this works?
#join: $(ACC).join.fq $(ACC).un1.fq $(ACC).un2.fq
join: $(ACC).join.fq

$(ACC).join.fq: $(ACC)_1.fastq $(ACC)_2.fastq
	fastq-join -v ' ' $(ACC)_1.fastq $(ACC)_2.fastq -r $(ACC).join.log -o $(ACC).%.fq

## Alignment
align: $(ACC).join.single.sam $(ACC).un1.single.sam $(ACC).un2.single.sam

$(ACC).join.single.sam: $(ACC).join.fq
	bwa mem data/monomers/single_mon.fa ./$(ACC).join.fq > $(ACC).join.single.sam 2> $(ACC).join.single.log

$(ACC).un1.single.sam: $(ACC).un1.fq
	bwa mem data/monomers/single_mon.fa ./$(ACC).un1.fq > $(ACC).un1.single.sam 2> $(ACC).un1.single.log

$(ACC).un2.single.sam: $(ACC).un2.fq
	bwa mem data/monomers/single_mon.fa ./$(ACC).un2.fq > $(ACC).un2.single.sam 2> $(ACC).un2.single.log

## Get Sorted SAM/BAM for Aligned Reads # TODO; need alignment; should it be clear?
aligned_sam: $(ACC).join.aligned.sort.sam $(ACC).un.aligned.sort.sam

$(ACC).join.aligned.sort.sam: $(ACC).join.single.sam
	samtools view -@ 12 -b -F 4 $(ACC).join.single.sam > $(ACC).join.aligned.bam
	samtools sort -@ 12 -o $(ACC).join.aligned.sort.bam $(ACC).join.aligned.bam && \
	rm $(ACC).join.aligned.bam
	samtools view $(ACC).join.aligned.sort.bam > $(ACC).join.aligned.sort.sam

$(ACC).un.aligned.sort.sam: $(ACC).un1.single.sam $(ACC).un2.single.sam
	samtools view -@ 12 -b -F 4 $(ACC).un1.single.sam > $(ACC).un1.aligned.bam
	samtools view -@ 12 -b -F 4 $(ACC).un2.single.sam > $(ACC).un2.aligned.bam
	samtools cat -o $(ACC).un.aligned.bam $(ACC).un1.aligned.bam $(ACC).un2.aligned.bam && \
	rm $(ACC).un{1,2}.aligned.bam
	samtools sort -@ 12 -o $(ACC).un.aligned.sort.bam $(ACC).un.aligned.bam && \
	rm $(ACC).un.aligned.bam
	samtools view $(ACC).un.aligned.sort.bam > $(ACC).un.aligned.sort.sam

## Get Fastq for Aligned Reads
aligned_fastq: $(ACC).join.aligned.fq $(ACC).un1.aligned.fq $(ACC).un2.aligned.fq

$(ACC).join.aligned.fq: $(ACC).join.single.sam
	samtools fastq -@ 12 -F 4 -0 $(ACC).join.aligned.fq $(ACC).join.single.sam

# TODO: why I cant set -1 ??
#$(ACC).un1.aligned.fq: $(ACC).un1.single.sam
$(ACC).un1.aligned.fq:
	samtools fastq -@ 12 -F 4 -0 $(ACC).un1.aligned.fq $(ACC).un1.single.sam

#$(ACC).un2.aligned.fq: $(ACC).un2.single.sam
$(ACC).un2.aligned.fq:
	samtools fastq -@ 12 -F 4 -0 $(ACC).un2.aligned.fq $(ACC).un2.single.sam

.PHONY: all join aligned_sam aligned_fastq
