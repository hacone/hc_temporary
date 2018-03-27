#!/bin/bash

#URL here

fastq-dump.2.9.0 SRR1997411 --split-files --skip-technical

# maybe this one is corrupted...
# fastq-dump.2.9.0 SRR3189743 --split-files --skip-technical

# fastq-join -v ' ' SRR3189743_1.fastq SRR3189743_2.fastq -r join.log -o SRR3189743.%.fq
