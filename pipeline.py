import os
import csv
import argparse
from Bio import Entrez
from Bio import SeqIO

'''
Lists to be used throughout code
'''
#List of accession #'s (D1-2dpi, D1-6dpi, D3-2dpi, D3-6dpi)
SRRs = ['SRR5660030.1', 'SRR5660033.1', 'SRR5660044.1', 'SRR5660045.1']
#List of conditions, corresponding to SRRs order
conditions = ['2dpi', '6dpi', '2dpi', '6dpi']
#List of what the paired-end fastq files will be named
fastqs = [['SRR5660030.1_1.fastq', 'SRR5660030.1_2.fastq'],['SRR5660033.1_1.fastq','SRR5660033.1_2.fastq'],['SRR5660044.1_1.fastq', 'SRR5660044.1_2.fastq'],['SRR5660045.1_1.fastq', 'SRR5660045.1_2.fastq']]
bt_fastqs = [['SRR5660030.1_mapped_1.fq', 'SRR5660030.1_mapped_2.fq'],['SRR5660033.1_mapped_1.fq','SRR5660033.1_mapped_2.fq'],['SRR5660044.1_mapped_1.fq', 'SRR5660044.1_mapped_2.fq'],['SRR5660045.1_mapped_1.fq', 'SRR5660045.1_mapped_2.fq']]

'''
Functions to run through steps
'''
#Step 1
def get_fastq(main_wd):
    #read in transcriptome urls from txt file into list, which corresponds to SRRs list
    file = main_wd + '/sra_links.txt'
    SRR_infile = open(file)
    urls = SRR_infile.read().splitlines()
    SRR_infile.close

    #get SRA files using wget cmd for each url
    for url in urls:
        cmd = 'wget ' + url
        os.system(cmd)

    #fastq-dump to convert to paired-end fastq files
    #establish generic beginning string to cmd
    fq_dump = 'fastq-dump -I --split-files '
    #use for loop to produce fq dump cmds for each SRR file
    for srr in SRRs:
        fqd_cmd = fq_dump + srr
        os.system(fqd_cmd)

#Step 2
def get_cds():
    Entrez.email = 'anovak9@luc.edu'

    #get fasta file of all cds in HCMV's genbank entry
    handle = Entrez.efetch(db = 'nucleotide', id = 'EF999921', rettype = 'fasta_cds_na', retmode = 'text')
    #store cds in list
    cds = list(SeqIO.parse(handle, 'fasta'))
    #write # of cds to log file
    outfile = open('miniProject.log', 'a')
    outfile.write('The HCMV genome (EF999921) has ' + str(len(cds)) + ' CDS.\n\n')
    outfile.close()
    #write list of cds to file (for running kallisto)
    SeqIO.write(cds, 'hcmv_cds.fasta', 'fasta')

#Step 3
def run_kall_sleuth(curr_dir):
    #get lists
    global SRRs
    global conditions
    global fastqs

    #create reference index for HCMV
    #don't need to change kmer defualt size b/c reads ~150 bp in length
    os.system('time kallisto index -i hcmv.idx hcmv_cds.fasta')

    #quantify each sample
    for i in range(len(SRRs)):
        cmd = 'time kallisto quant -i hcmv.idx -o ' + SRRs[i] + '_results -b 30 -t 2 ' + fastqs[i][0] + ' ' + fastqs[i][1]
        os.system(cmd)


    #make table for sleuth
    header = ['sample', 'condition', 'path']
    rows = []
    #for each sample
    for i in range(len(SRRs)):
        #path to sample's result file
        path = curr_dir + '/' + SRRs[i] + '_results'
        #make list of sample, its condition, & its path
        row = [SRRs[i], conditions[i], path]
        rows.append(row)

    #write table to csv file
    with open('input_table.csv', 'w') as csvfile:
        #create csv writer object
        csvwriter = csv.writer(csvfile)
        #write headers
        csvwriter.writerow(header)
        #write sample rows
        csvwriter.writerows(rows)

    #run sleuth script
    cmd = 'Rscript ' + curr_dir + '/sleuth_code.R'
    os.system(cmd)

#Step 4
def run_bowtie():
    global SRRs
    global conditions
    global fastqs
    global bt_fastqs

    #get # of reads in each transcriptome before mapping
    before = []
    for fastq in fastqs:
        file = fastq[0]
        reads = len(list(SeqIO.parse(file, 'fastq')))
        before.append(reads)

    #build hcmv index for bowtie2
    os.system('bowtie2-build hcmv_cds.fasta hcmv_bowtie_idx')

    #map each sample
    for i in range(len(SRRs)):
        #--al-conc option will only save reads that did map to hcmv index into .fq files
        cmd = 'bowtie2 --quiet -x hcmv_bowtie_idx -1 ' + fastqs[i][0] + ' -2 ' + fastqs[i][1] + ' -S ' + SRRs[i] + '_map.sam --al-conc ' + SRRs[i] + '_mapped_%.fq'
        os.system(cmd)

    #get # of reads in each transcriptome after mapping
    after = []
    for fastq in bt_fastqs:
        file = fastq[0]
        reads = len(list(SeqIO.parse(file, 'fastq')))
        after.append(reads)

    #make before & after lists str datatype
    before = list(map(str, before))
    after = list(map(str, after))

    #get strings to write # reads to file
    outputs = []
    for i in range(len(SRRs)):
        #Donor 1
        if i == 0 or i == 1:
            output = 'Donor 1 (' + conditions[i] + ') had ' + before[i] +\
                 ' read pairs before Bowtie2 filtering and ' + after[i] + ' read pairs after.'
        #Donor 3
        else:
            output = 'Donor 3 (' + conditions[i] + ') had ' + before[i] +\
                 ' read pairs before Bowtie2 filtering and ' + after[i] + ' read pairs after.'
        outputs.append(output)
    #write to log file
    outfile = open('miniProject.log', 'a')
    outfile.write('\n\n'+'\n'.join(outputs)+'\n\n')
    outfile.close()

#Step 5
def run_spades(kmers):
    global bt_fastqs

    #write spades assembly command to make read-pair libraries
    cmd = 'spades -k ' + kmers + ' -t 2 --only-assembler '
    #specify paired-end libraries
    for i in range(4):
        string = '--pe' + str(i+1) + '-1 ' + bt_fastqs[i][0] + ' --pe' + str(i+1) + '-2 ' + bt_fastqs[i][1] + ' '
        cmd += string
    cmd += '-o hcmv_assembly/'

    #write spades cmd to log file
    outfile = open('miniProject.log', 'a')
    outfile.write(cmd + '\n\n')
    outfile.close()

    #run assembly
    os.system(cmd)

#Step 6/7
def get_contigs():

    #get contigs fasta file from assembly
    contigs = list(SeqIO.parse('hcmv_assembly/contigs.fasta', 'fasta'))
    #filter contigs >1000 bp
    filt_contigs = []
    #for each contig in file
    for contig in contigs:
        #if length of seq greater than 1000 bp
        if len(contig.seq) > 1000:
            #add to filtered list
            filt_contigs.append(contig)
    #print num of contiges to log file
    output = 'There are ' + str(len(filt_contigs)) + ' contigs > 1000 bp in the assembly.'

    #calculate length of assembly
    num_bp = 0
    for contig in filt_contigs:
        length = len(contig.seq)
        num_bp += length

    #add to output
    output += '\nThere are ' + str(num_bp) + ' bp in the assembly.\n\n'

    #write output to log file
    outfile = open('miniProject.log', 'a')
    outfile.write(output)
    outfile.close()

    #return list of filtered contigs
    return filt_contigs

#To use in step 8
#Returns list of top 10 hits
def parse_blast(filename):
    #list to hold top 10 hits
    x = []
    #open csv file and split lines into list rows
    blast_results = open(filename, 'r')
    rows = blast_results.read().splitlines()
    #for the first 10 rows
    for i in range(10):
        #replace commas w/a tab delimiter
        row = rows[i].replace(',', '\t')
        #add to list x
        x.append(row)
    blast_results.close()
    return x

#Step 8
def blast(contigs):
    #get longest contig from list
    max = 0
    for contig in contigs:
        if len(contig.seq) > max:
            longest = contig
    #write longest contig to fasta file
    SeqIO.write(longest, 'max_contig.fasta', 'fasta')

    #run a blastn on longest contig
    cmd = 'blastn -query max_contig.fasta -db hcmv_db -out hcmv_blastn_results.csv -outfmt "10 sacc pident length qstart qend sstart send bitscore evalue stitle" '
    os.system(cmd)

    #get top 10 hits froom results into list 'rows'
    rows = parse_blast('hcmv_blastn_results.csv')
    headers = 'sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle'
    #list for output
    output = []
    output.append(headers)
    output.extend(rows)

    #write to log file 
    outfile = open('miniProject.log', 'a')
    outfile.write('\n'.join(output))
    outfile.close()


'''
Run Functions
'''
#Set up argparser
parser = argparse.ArgumentParser(description = 'Run python wrapper.')
#Flag to specify what dataset to use (response either 'full' or 'test')
parser.add_argument('dataset', type = str, help = 'The test set you want to run: full or test')
args = parser.parse_args()
print(args.dataset)

#if testing full dataset
if args.dataset == 'full':
    #get current working directory
    orig_cwd = os.getcwd()
    dir = 'miniProject_Annie_Novak'
    #change working directory into mini proj folder
    cwd = orig_cwd + '/' + dir
    os.chdir(cwd)
    #create output log file
    outfile = open('miniProject.log', 'w')
    outfile.close()

    #run functions
    kmers = '55,77,99,127'
    get_fastq(orig_cwd)
    get_cds()
    run_kall_sleuth(cwd)
    run_bowtie()
    run_spades(kmers)
    contigs = get_contigs()
    blast(contigs)

#else run test dataset
elif args.dataset == 'test':
    #get current working directory
    cwd = os.getcwd()
    dir = 'test_outputs'
    #change working directory into test folder
    cwd = cwd + '/' + dir
    os.chdir(cwd)
    #create output log file
    outfile = open('miniProject.log', 'w')
    outfile.close()

    #run functions
    kmers = '127'
    get_cds()
    run_kall_sleuth(cwd)
    run_bowtie()
    run_spades(kmers)
    contigs = get_contigs()
    blast(contigs)
'''
End
'''
