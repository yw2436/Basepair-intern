#prepare the genome, the genome and gtf filename might change

STAR --runMode genomeGenerate --runThreadN 16 --genomeDir genome/GRCz10/ --genomeFastaFiles genome/GRCz10/danRer10.fa --sjdbGTFfile gtfFile/GRCz10/Danio_rerio.GRCz10.91.gtf
STAR --runMode genomeGenerate --runThreadN 16 --genomeDir genome/GRCz11/ --genomeFastaFiles genome/GRCz11/Danio_rerio.GRCz11.dna.primary_assembly.fa --sjdbGTFfile gtfFile/GRCz11/Danio_rerio.GRCz11.92.chr.gtf

#reads QC

cd reads/
mkdir fastqc
for i in $(ls); do echo $i; fastqc -t 8 -o fastqc $i; done

#Align
cd ..
STAR --runMode alignReads --runThreadN 16 --genomeDir genome/GRCz10/ --readFilesIn reads/wt_fasta/Miller_wt80.fastq --outFileNamePrefix reads/STAR_wt/wt80
STAR --runMode alignReads --runThreadN 16 --genomeDir genome/GRCz10/ --readFilesIn reads/wt_fasta/Miller_nhslWT.fastq --outFileNamePrefix reads/STAR_wt/nhslWT
STAR --runMode alignReads --runThreadN 16 --genomeDir genome/GRCz10/ --readFilesIn reads/wt_fasta/Miller_hox80.fastq  --outFileNamePrefix reads/STAR_wt/hox80
STAR --runMode alignReads --runThreadN 16 --genomeDir genome/GRCz10/ --readFilesIn reads/wt_fasta/Miller_nhslMUT.fastq --outFileNamePrefix reads/STAR_wt/nhslMUT

samtools mpileup -ABuf /mnt/bioinfo/RNAmapper_test2/genome/GRCz10/Danio_rerio.GRCz10.dna.toplevel.fa /mnt/bioinfo/RNAmapper_test2/reads/STAR_wt/nhslwt_sorted.bam | bcftools view -O b --threads 8 > /mnt/bioinfo/RNAmapper_test2/reads/STAR_wt/nhslwt_raw.bcf 
samtools mpileup -ABuf /mnt/bioinfo/RNAmapper_test2/genome/GRCz10/Danio_rerio.GRCz10.dna.toplevel.fa /mnt/bioinfo/RNAmapper_test2/reads/STAR_wt/wt80_sorted.bam | bcftools view -O b --threads 8 > /mnt/bioinfo/RNAmapper_test2/reads/STAR_wt/wt80_raw.bcf 
samtools mpileup -ABuf /mnt/bioinfo/RNAmapper_test2/genome/GRCz10/Danio_rerio.GRCz10.dna.toplevel.fa /mnt/bioinfo/RNAmapper_test2/read/STAR_mut/hox80_sorted.bam | bcftools view -O b --threads 8 > /mnt/bioinfo/RNAmapper_test2/reads/STAR_mut/hox80_raw.bcf 
samtools mpileup -ABuf /mnt/bioinfo/RNAmapper_test2/genome/GRCz10/Danio_rerio.GRCz10.dna.toplevel.fa /mnt/bioinfo/RNAmapper_test2/read/STAR_mut/nhslmut_sorted.bam | bcftools view -O b --threads 8 > /mnt/bioinfo/RNAmapper_test2/reads/STAR_mut/nhslmut_raw.bcf   ##--outSAMstrandField intronMotif 

for i in {1..25}; do bcftools view nhslwt_raw.bcf | grep -w ^$i > nhslwt_chr$i.vcf; done &
for i in {1..25}; do bcftools view wt80_raw.bcf | grep -w ^$i > wt80_chr$i.vcf; done &
for i in {1..25}; do bcftools view hox80_raw.bcf  | grep -w ^$i > hox80_chr$i.vcf; done &
for i in {1..25}; do bcftools view nhslmut_raw.bcf  | grep -w ^$i > nhslmut_chr$i.vcf; done &
