import gzip, sys, getopt, shlex, subprocess, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.patches as mpatches
import os
import pandas as pd
import seaborn as sns



#Download required softwares and files########################################
subprocess.call(['mkdir' ,'annotation'])
subprocess.call(['mkdir' ,'alignment'])
subprocess.call(['mkdir' ,'counts'])
subprocess.call(['mkdir' ,'report'])
subprocess.call(['mkdir' ,'splice'])
subprocess.call(['wget' ,'https://github.com/SiddharthTomar/ICARuS2/raw/master/softwares/BBMap.tar.gz'])
stream = os.popen("tar xvzf BBMap.tar.gz")
subprocess.call(['wget' ,'https://github.com/SiddharthTomar/ICARuS2/raw/master/softwares/subread.tar.gz'])
stream = os.popen("tar xvzf subread.tar.gz")
subprocess.call(['wget' ,'https://github.com/SiddharthTomar/ICARuS2/raw/master/softwares/hisat.zip'])
subprocess.call(['unzip', 'hisat.zip'])
stream = os.popen("mv ./hisat2-2.1.0 ./hisat")
stream = os.popen("mv ./subread-1.5.2-Linux-x86_64 ./subread")
# #stream = os.popen("mkdir reads")
# #subprocess.call(['curl' ,'ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grcm38_snp_tran.tar.gz', '-o', 'index.tar.gz'])
# #stream = os.popen("mv ./grcm38_snp_tran.tar.gz ./reference/")
# #stream = os.popen("tar xvzf ./refrence/grcm38_snp_tran.tar.gz")
subprocess.call(['wget' , 'https://github.com/SiddharthTomar/ICARuS2/raw/master/softwares/samtools-1.5.tar.bz2'])
subprocess.call(['wget' , 'https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip'])
subprocess.call(['unzip', 'fastqc_v0.11.5.zip'])
subprocess.call(['wget' , 'https://raw.githubusercontent.com/SiddharthTomar/ICARuS2/master/Rscripts/postqc.r'])
subprocess.call(['wget' , 'https://raw.githubusercontent.com/SiddharthTomar/ICARuS2/master/Rscripts/preqc.r'])
subprocess.call(['wget' , 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M14/gencode.vM14.annotation.gtf.gz' ])
stream = os.popen("mv ./gencode.vM14.annotation.gtf.gz ./annotation/")
subprocess.call(['gzip' , '-d', './annotation/gencode.vM14.annotation.gtf.gz' ])

#----------------------------------------------------------------------------#

#Read configure file and store the parameters#################################
configename = "config.txt"
with open(configename) as f:
	content = f.readlines()
content = [x.strip() for x in content]
threads = content[0]
readtype = content [1]
#direction = content [2]
direction = "RF"
#clippingF = content [3]
clippingF = "13"
clippingR = content [4]
reads = content [5:]
strand = 2
#----------------------------------------------------------------------------#

#Read download file and get the reads#########################################
#configename = "url.txt"
#with open(configename) as f:
#	content = f.readlines()
#content = [x.strip() for x in content]
#path = "/reads/"
#for i in content:
#	name = i.split('/')[-1]
#	path = "./reads/"
#	name =  path + name
#	subprocess.call(['curl', i, '--output' , name])
#----------------------------------------------------------------------------#

#Generate data for plots######################################################																								
subprocess.call(['./FastQC/fastqc', '--threads', '7', './reads/*.fastq.gz'])#						
subprocess.call(['Rscript','preqc.r'])#										
#----------------------------------------------------------------------------#


#BBDUK IT!!!##################################################################
#Values defined for forward and reserve trim will be used#####################
for i in reads:
	path1 = "in1=" + i + "_1.fastq.gz"
	path2 = "in2=" + i + "_2.fastq.gz"
	out1 = "out1=./readProcessed/" + i.split('/')[-1] + "_1.fastq.gz"
	out2 = "out2=./readProcessed/" + i.split('/')[-1] + "_2.fastq.gz"
	ftl = "ftl=" + clippingF
	subprocess.call(['./bbmap/bbduk.sh', ftl, path1, path2, out1, out2])
	delete1 = "rm ./reads/" + i.split('/')[-1] + "_1.fastq.gz"
	delete2 = "rm ./reads/" + i.split('/')[-1] + "_2.fastq.gz"
	stream = os.popen(delete1)
	stream = os.popen(delete2)
	path1 = "in1=./readProcessed/" + i.split('/')[-1] + "_1.fastq.gz"
	path2 = "in2=./readProcessed/" + i.split('/')[-1] + "_2.fastq.gz"
	out1 = "out1=./readProcessed/" + i.split('/')[-1] + "_cleaned_1.fastq.gz"
	out2 = "out2=./readProcessed/" + i.split('/')[-1] + "_cleaned_2.fastq.gz"
	report = i + ".txt"
	subprocess.call(['./bbmap/bbduk.sh', path1, path2, out1, out2, "k=31", "ref=./reference/mrrna.fasta", ">>", report])
	delete1 = "rm ./readProcessed/" + i.split('/')[-1] + "_1.fastq.gz"
	delete2 = "rm ./readProcessed/" + i.split('/')[-1] + "_2.fastq.gz"
	stream = os.popen(delete1)
	stream = os.popen(delete2)
	out1 = "./readProcessed/" + i.split('/')[-1] + "_cleaned_1.fastq.gz"
	out2 = "./readProcessed/" + i.split('/')[-1] + "_cleaned_2.fastq.gz"
	subprocess.call(['./FastQC/fastqc', '--threads', '7', out1, out2])
#----------------------------------------------------------------------------#

#Generate data for plots######################################################
subprocess.call(['Rscript','postqc.r'])												
df = pd.read_csv('postq.csv')
#----------------------------------------------------------------------------#

#Align it using HISAT#########################################################
for i in reads:
	in1 = "./readProcessed/" + i.split('/')[-1] + "_cleaned_1.fastq.gz"
	in2 = "./readProcessed/" + i.split('/')[-1] + "_cleaned_2.fastq.gz"
	alignment = "./alignment/" + i.split('/')[-1]+ ".sam"
	alignmentbam = "./alignment/" + i.split('/')[-1]+ ".bam"
	sortedBam = "./alignment/" + i.split('/')[-1]+ "_sorted.bam"
	splicesite = "./splice/" + i.split('/')[-1]+ ".txt"
	filename = "./splice/" + i.split('/')[-1]+ ".log"
	file = open(filename, 'a')
	subprocess.call(['./hisat/hisat2','-f','-x','./reference/genome_snp_tran', '-q', '-1', in1, '-2', in2, '-S', alignment, '--no-discordant', '--rna-strandness', direction, '--dta-cufflinks', '--novel-splicesite-outfile', splicesite, '-p', '7'], stdout = file)
	file.close()
	subprocess.call(['samtools', 'view', '-b', '-S', alignment ,'-o', alignmentbam])
	remove = "rm " + "./alignment/" + i.split('/')[-1]+ ".sam"
	stream = os.popen(remove)
	subprocess.call(['samtools', 'sort', alignmentbam, '-o', sortedBam])
	remove = "rm " + "./alignment/" + i.split('/')[-1]+ ".bam"
	stream = os.popen(remove)
	delete1 = "rm ./readProcessed/" + i.split('/')[-1] + "_cleaned_1.fastq.gz"
	delete2 = "rm ./readProcessed/" + i.split('/')[-1] + "_cleaned_2.fastq.gz"
	stream = os.popen(delete1)
	stream = os.popen(delete2)
#----------------------------------------------------------------------------#

#Feature count################################################################
for i in reads:
	alignmentbam = "./alignment/" + i.split('/')[-1] + ".bam"
	countfile = './counts/' + i.split('/')[-1] + '_count.txt'
	subprocess.call(['./subread/bin/featureCounts', '-T', threads, '-s', strand,'-p', '-t', 'exon', '-g', 'gene_id', '-a', './annotation/gencode.vM14.annotation.gtf', '-R', '-o', countfile, alignmentbam])
#----------------------------------------------------------------------------#
