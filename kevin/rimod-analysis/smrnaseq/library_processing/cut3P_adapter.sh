## Lane 1
# cut P7 adapter
cutadapt -a ATCTCGTATGCCGTCTTCTGCTTGX -o all_lanes_trimmed3p.fastq --untrimmed-output=no_adapter.fastq ./basecalled_fastq/smRNAseq_all_lanes.fastq.gz > cutadapt_info.txt
