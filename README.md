# CompBio_Pipeline
##Step 1: 
The objective of this project is to compare HCMV transcriptomes 2- and 6-days post-infection (dpi). 
All of the following transcriptomes below were from two patient donors from SRA and converted to paired-end fastq files. Go to the SRA and under runs click the lit up SRA link. This link will take you to a sequence archive, where the website address can be found under DATA access. 

Donor 1 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896360

Donor 1 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896363 

Donor 3 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896374 

Donor 3 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896375

Example for one of the samples: 
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030
fastq-dump -I --split-files SRR5660030



##Packages used 

Os

Glob 

Subprocess 


From biopython: 

Entrez 

SeqIO 

**Spades, Bowtie and NCBI are all tools used in the code - Spades and Bowtie are software tools that need to be installed but NCBI needs an email and username.** 

##Copy respository 
If you would like to clone this repository the command is

git clone https://github.com/smjablonska/CompBio_Pipeline.git 

##Sample input data 
Command used: head -n 40000 data.fastq > sampledata.fastq
File attached in sampledata.fastq file

