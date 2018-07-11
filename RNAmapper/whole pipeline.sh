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

