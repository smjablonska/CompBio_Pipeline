import os
from Bio import Entrez
from Bio import SeqIO
from zipfile import ZipFile

#Making paths as I go so everything is nice and cleanly kept
base_path = os.getcwd() 
base_path = base_path+'/PipelineProject_Run'
os.makedirs(base_path)
log_file = open(base_path+'/PipelineProject.log', 'w+') #Getting log file ready so it can be open until the end of the code 


#First start out with creating index for reference genome
def index(paths):    
    fileloco= paths+'/bowtie_fasta' 
    #Want a new directory for HCMV
    os.makedirs(paths+'/bowtie_index')
    indexloco = paths+'/bowtie_index'    
    #Now we are searching for the wanted HCMV accession id and putting records into our reference directory 
    handle = Entrez.efetch(db="nucleotide", id= 'NC_006273.2', rettype='fasta')
    records = list(SeqIO.parse(handle, 'fasta')) 
    with open("NCref.fasta", "w") as output_handle: 
        SeqIO.write(records, fileloco, "fasta")
    #Want to build the bowtie index where all the HCMV files are 
    os.system('bowtie2-build '+fileloco+' '+indexloco+'/HCMV')
    

def mapping(paths): #Now we need to map HCMV to files in directory   
    path = os.getcwd()+'/inputdata'
    os.makedirs(paths+'/mapping')
    #for each file in the sample data fasta folder
     #open the file and append each read to a list
    with open("sampledata.fasta", 'r') as f:
        reads = []
        for line in f:
            if line.startswith("@"): #Go through sample and append file name because we have a mix of sample data 
                reads.append(line)
    
    reads = sorted(reads) 
    #this is the path to the HCMV indexes
    index_path = paths+'/bowtie_index/NCref' #Have the index pathway 
    file_name = paths+'/mapping'     #And file output out so only the mapped reads will go out 
    #for each name of a read - run the bowtie2 command and output  
    for i in reads:
        command = 'bowtie2 --quiet -x '+index_path+' -1 '+path+'/'+reads[i]+' -2 '+path+'/'+reads[i]+' -S '+file_name+' --al-conc-gz '+path+'/aligned_reads/'+reads[i]+'_mapped_%.fq.gz'
        return command 
    
    os.system(command)
  
# #this counts the reads before and after the bowtie filtering
# Function to count reads in a FASTQ file
def count_reads_fastq_transcripts(paths):
    #list of files before  
    prev_path = os.getcwd()+'/inputdata' #Where all fasta files are kept
    prev_files = [] 
    unfiltered_numbers = [] 
    #Next I want to associate each file name with length of reads before/after 
    #We will start with before first and sorting was recommended to me since its easy to align file names with outputs  
    for filename in os.listdir(prev_path):
         f = os.path.join(prev_path, filename)
         prev_files.append(filename)      
    #Want to sort names in order first 
    prev_files = sorted(prev_files)
    #for each file name in the before cateogory, we want to know how many reads there are
    for line in prev_files:
        summary = sum(1 for line in f) // 4 #Easier for before to just split reads into 4 
        unfiltered_numbers.append([len(summary)])       

###After filtering section! 
    #Getting into the after filtering section of reads 
    after_path = paths+'/alignedreads'
    after_files = []
    filt_numbers = []
    #gets the path for where the unzipped filed are
    unzipped_path = paths+'/alignedreadsunzipped'
    #directory for unzipping 
    os.makedirs(paths+'/alignedreads_unzip')
    with ZipFile('.zip', 'r') as f:
        #extract in different directory
        f.extractall('alignedreads_unzip')
    #Now we need to unzip all the files that were filtered, otherwise number will not be accurate 
    for filename in os.listdir("alignedreads.unzip"): #Again use sorting because filename can correspond to numbers 
        f =  os.path.join(after_path, filename)
        after_files.append(filename)
        after_files = sorted(after_files)

    #gets the number of reads of the filtered file
    for i in after_files:
        file = unzipped_path+'/'+i
        with open(file, 'r') as f:
            reads = []
            for line in f:
                if line.startswith("@"):
                    reads.append(line)
        
        filt_numbers.append([len(reads)])
        
    #this writes out the number of reads in the instructed way 
    log_file.write('Donor 1 (2dpi) had '+str(unfiltered_numbers[0])+' read pairs before Bowtie2 filtering and '+str(filt_numbers[1])+' read pairs after.'+'\n')
    log_file.write('Donor 1 (6dpi) had '+str(unfiltered_numbers[1])+' read pairs before Bowtie2 filtering and '+str(filt_numbers[2])+' read pairs after.'+'\n')
    log_file.write('Donor 3 (2dpi) had '+str(unfiltered_numbers[2])+' read pairs before Bowtie2 filtering and '+str(filt_numbers[3])+' read pairs after.'+'\n')
    log_file.write('Donor 3 (6dpi) had '+str(unfiltered_numbers[3])+' read pairs before Bowtie2 filtering and '+str(filt_numbers[4])+' read pairs after.'+'\n')

    #Using rnaspades.py because we are looking at a virus and it is faster to use
    files = after_files 
    commandspades = '/home/sjablonska/software/spades/rnaspades.py' 
    output = paths+'/spades_assembly'
    spadescommand = 'python3 '+commandspades+' --pe-1 1 '+files[0]+' --pe-2 1 '+files[1]+' --pe-1 2 '+files[2]+' --pe-2 2 '+files[3]+' --pe-1 3 '+files[4]+' --pe-2 3 '+files[5]+' --pe-1 4 '+files[6]+' --pe-2 4 '+files[7]+' -o '+output
    os.system(spadescommand)
    #Writing out command to log file 
    log_file.write("Spades command: "+spadescommand+'\n')
  
def numberofcontigs(paths):
    #Make the path to all seperated fasta files 
    filename = paths+'/spades_raw/transcriptomes.fasta'  
    records = SeqIO.parse(filename, "fasta")
    base_counter = 0 #Initiate 2 counters to know the contig count and base count 
    contig_counter = 0 
    #Use seq record to make it easier for contigs to be counted 
    for record in records:       
        if len(record.seq) > 1000:
            contig_counter+=1 #Contig presence is measured through occurance,so 1 is added 
            base_counter+=len(record.seq) #Full number of base pairs added up from each iteration in seq      
    #Writes out informaton to log file 
    log_file.write('There are '+str(contig_counter)+' contigs > 1000 bp in the assembly.'+ '\n')
    log_file.write('There are '+str(base_counter)+' bp in the assembly.'+'\n')
    
  
def blasteverything(paths):    
    #making a different folder for only blast information 
    os.makedirs(paths+'/blastdirectory')
    #Path to raw spades output 
    filename = paths+'/spades_raw/transcriptomes.fasta'
    #this parses the fasta file
    records = SeqIO.parse(filename, "fasta")
    max_len = 0
    max_contig = ''
    #We want the longest contig in our sequence to blast against 
    for record in records:
        if len(record.seq) > max_len:
            max_len = len(record.seq)
            #Keeps appending new number until longest contig is found 
            max_contig = str(record.seq) 
    
    #Need HCMV file in the same directory as blast database to use it to blast against 
    HCMVdb = open(paths + "/blast/HCMV.fasta", "w")

    #Write longest contig into separate fasta file
    longest_contig = open(paths+'/blastdirectory/longest_contig.fasta' ,'w')
    longest_contig.write(max_contig) #Sequence is written to outfile to show what the contig is 
    longest_contig.close()
    
    #Getting only nucleotide records from Betaherpesvirinae family
    handle = Entrez.esearch(db = "nucleotide", term = ("Betaherpesvirinae[All Fields]"), retmax = 2000) 
    record = Entrez.read(handle)
    
    Beta_ids = record["IdList"]
    handle = Entrez.efetch(db = "nucleotide", id = Beta_ids, rettype = "fasta")
    
    records = list(SeqIO.parse(handle, format = "fasta"))

    for record in records: #for each record, write to fasta in HCMV database
        HCMVdb.write(str(record.seq) + "\n")
        
    HCMVdb.close()

    #All commands need to be at the end because we just made a pipeline for all the information. 
    #Making blastDB 
    db_path = paths+'/blast/HCMV.fasta'
    output_db = paths+'/blast/blastdb/blastresults'
    db_command = 'makeblastdb -in '+db_path+' -out '+output_db+' -title HCMV'+' -dbtype nucl'
    os.system(db_command)

    #Blast longest contig file and write out results
    query_seqfile = paths+'/blast/longcontig.fasta'
    output_path = paths+'/blast/blast_outfile'
    blast_command = 'blastn -query '+query_seqfile+' -db '+output_db+" -num_threads 2 -max_target_seqs 10 -out "+output_path+" -outfmt='6 sacc pident length qstart qend sstart send bitscore evalue stitle'"
    log_file.write(blast_command)
    os.system(blast_command)


index(base_path)
mapping(base_path)
count_reads_fastq_transcripts(base_path)
numberofcontigs(base_path)
blasteverything(base_path)

log_file.close()
