# cut P7 adapter
cutadapt  -a ATCTCGTATGCCGTCTTCTGCTTG -o trimmed_3prime.fastq --untrimmed-output=no_adapter.fastq Undetermined_S0_L003_R1_001.fastq

# perform the demultiplexing
cutadapt -a file:barcodes.fasta --untrimmed-o no_barcode.fastq -o sample_{name}.fastq trimmed_3prime.fastq

# unzip fastq.gz files for concatenating
for f in *.fastq.gz; do gunzip $f; done;

# concatenate 

# cut the 3'-adapter
cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC -o trimmed_sample.fastq --untrimed-output=no_adapter.fastq 

# cut the 4 random 3'-bases
cutadapt -u -4 -o trimmed2_sample10316.fastq trimmed_sample_10316.fastq

# Shorten 5'-end by removing N's and clipping read to maximal 25 bp length
cutadapt --trim-n -l -25 
