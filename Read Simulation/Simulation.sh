# I used chr1 of hg38 as the sorce sequence. 
#SNP information is from UCSC Gene Browser and filtered that it only contains info about chr1 and every line has 26 columns. 
#FASTQ info is from SRR1534766.

# Get chr1 for simulation
sed -n '1,4979130p' hg38.fa > hg38_chr1.fa

# Simulate methylation rates for every cytosine in chr1
fasta-methyl-sim hg38_chr1.fa > meth.fa

# Get 1 million length-85 chunks, simulate bisulfite conversion, and simulate sequencing errors
fasta-random-chunks -n1000000 -s85 poly.fa |
fasta-bisulf-sim |
fastq-sim - sample.fastq > reads.fastq
