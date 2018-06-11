# ACC=SRR3189743
# ACC=SRR1997411

all: aligned_sam aligned_fastq

## delegate tasks for each dataset to Makefiles in their dirs 
makefile: data/$(ACC)/Makefile

data/$(ACC)/Makefile: Makefile.tmpl
	echo "ACC=$(ACC)" > data/$(ACC)/Makefile
	cat Makefile.tmpl >> data/$(ACC)/Makefile

## Check reads are present
reads: makefile
	cd data/$(ACC); make reads

## Join paired-end reads
join: makefile
	cd data/$(ACC); make join

## Alignment using whole set of monomers & reads
align: makefile
	cd data/$(ACC); make align

## Get Sorted SAM/BAM for Aligned Reads
aligned_sam: makefile
	cd data/$(ACC); make aligned_sam

## Get Fastq for Aligned Reads
aligned_fastq: makefile
	cd data/$(ACC); make aligned_fastq

subalign: makefile
	cd data/$(ACC); make subalign

sub_bam: makefile
	cd data/$(ACC); make sub_bam

s14align: makefile
	cd data/$(ACC); make s14align

s14bam: makefile
	cd data/$(ACC); make s14bam

.PHONY: all makefile reads join aligned_sam aligned_fastq subalign sub_bam s14align s14bam
