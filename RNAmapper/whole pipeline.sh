mkdir wt_fasta mut_fasta  #to store all the reads
mkdir genome gtfFiles #store the genome and gtf files
mkdir alignment
mkdir alignment/STAR_wt alignment/STAR_mut  #to store the bam files
mkdir vcf_files #to store the vcf files for every chr
mkdir cufflinks



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
STAR --runMode alignReads --runThreadN 16 --genomeDir genome/GRCz10/ --readFilesIn reads/mut_fasta/Miller_hox80.fastq  --outFileNamePrefix reads/STAR_wt/hox80
STAR --runMode alignReads --runThreadN 16 --genomeDir genome/GRCz10/ --readFilesIn reads/mut_fasta/Miller_nhslMUT.fastq --outFileNamePrefix reads/STAR_wt/nhslMUT

samtools mpileup -ABuf /mnt/bioinfo/RNAmapper_test2/genome/GRCz10/Danio_rerio.GRCz10.dna.toplevel.fa /mnt/bioinfo/RNAmapper_test2/reads/STAR_wt/nhslwt_sorted.bam | bcftools view -O b --threads 8 > /mnt/bioinfo/RNAmapper_test2/reads/STAR_wt/nhslwt_raw.bcf 
samtools mpileup -ABuf /mnt/bioinfo/RNAmapper_test2/genome/GRCz10/Danio_rerio.GRCz10.dna.toplevel.fa /mnt/bioinfo/RNAmapper_test2/reads/STAR_wt/wt80_sorted.bam | bcftools view -O b --threads 8 > /mnt/bioinfo/RNAmapper_test2/reads/STAR_wt/wt80_raw.bcf 
samtools mpileup -ABuf /mnt/bioinfo/RNAmapper_test2/genome/GRCz10/Danio_rerio.GRCz10.dna.toplevel.fa /mnt/bioinfo/RNAmapper_test2/read/STAR_mut/hox80_sorted.bam | bcftools view -O b --threads 8 > /mnt/bioinfo/RNAmapper_test2/reads/STAR_mut/hox80_raw.bcf 
samtools mpileup -ABuf /mnt/bioinfo/RNAmapper_test2/genome/GRCz10/Danio_rerio.GRCz10.dna.toplevel.fa /mnt/bioinfo/RNAmapper_test2/read/STAR_mut/nhslmut_sorted.bam | bcftools view -O b --threads 8 > /mnt/bioinfo/RNAmapper_test2/reads/STAR_mut/nhslmut_raw.bcf   ##--outSAMstrandField intronMotif 

for i in {1..25}; do bcftools view nhslwt_raw.bcf | grep -w ^$i > nhslwt_chr$i.vcf; done &
for i in {1..25}; do bcftools view wt80_raw.bcf | grep -w ^$i > wt80_chr$i.vcf; done &
for i in {1..25}; do bcftools view hox80_raw.bcf  | grep -w ^$i > hox80_chr$i.vcf; done &
for i in {1..25}; do bcftools view nhslmut_raw.bcf  | grep -w ^$i > nhslmut_chr$i.vcf; done &

#need to move all of the *chr*files to a directory in order to run the RNAmappe.R script
#also need to make a directory called _settings, and a file called mappeRsettings.txt containing a string '“coverage=25” “zygosity=25”' within that directory
for i in {1..25}; do Rscript --vanilla /mnt/bioinfo/RNAmapper_test2/scripts/RNAmapperScripts/RNAmappe.R \
wt80_chr$i.vcf,hox80_chr$i.vcf; done
for i in {1..25}; do Rscript --vanilla /mnt/bioinfo/RNAmapper_test2/scripts/RNAmapperScripts/RNAmappe.R \
nhslwt_chr$i.vcf,nhslmut_chr$i.vcf; done



#run cufflink
#first set up directories sepreately for different samples
cufflinks -o _cufflinks_wt/wt80 -p 8 -g /mnt/bioinfo/RNAmapper_test2/gtfFile/GRCz10/Danio_rerio.GRCz10.91.gtf -u /mnt/bioinfo/RNAmapper_test2/reads/STAR_wt_2/wt80/Aligned.sortedByCoord.out.bam
cufflinks -o _cufflinks_wt/nhslwt -p 8 -g /mnt/bioinfo/RNAmapper_test2/gtfFile/GRCz10/Danio_rerio.GRCz10.91.gtf -u /mnt/bioinfo/RNAmapper_test2/reads/STAR_wt_2/nhslWT/Aligned.sortedByCoord.out.bam
cufflinks -o _cufflinks_mut/nhslmut -p 8 -g /mnt/bioinfo/RNAmapper_test2/gtfFile/GRCz10/Danio_rerio.GRCz10.91.gtf -u /mnt/bioinfo/RNAmapper_test2/reads/STAR_mut_2/nhsl_mut/Aligned_sorted.bam
cufflinks -o _cufflinks_mut/hox80 -p 8 -g /mnt/bioinfo/RNAmapper_test2/gtfFile/GRCz10/Danio_rerio.GRCz10.91.gtf -u /mnt/bioinfo/RNAmapper_test2/reads/STAR_mut_2/hox80/Aligned.sortedByCoord.out.bam



