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
direction = content [2]
#This are non-functional
clippingF = content [3]
clippingR = content [4]
reads = content [5:]
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

#Pre cleanup QC plot############################################################
for i in reads:																	#
	path1 = i + "_1.fastq.gz"													#
	path2 = i + "_2.fastq.gz"													#
	subprocess.call(['./FastQC/fastqc', path1, path2])							#
subprocess.call(['Rscript','preqc.r'])											#
df = pd.read_csv('preq.csv')													#
df = df.groupby(['sample', 'tot.seq', 'module'])['status'].apply(', '.join).unstack().reset_index().rename_axis(None, axis=1)
mapping = {'PASS': 1, 'WARN': 2, 'FAIL': 3}										#
df = df.replace({'Basic Statistics': mapping, 'Per base sequence quality': mapping, 'Per tile sequence quality': mapping, 'Per sequence quality scores': mapping, 'Per base sequence content': mapping, 'Per sequence GC content': mapping, 'Per base N content': mapping, 'Sequence Length Distribution': mapping, 'Sequence Duplication Levels': mapping, 'Overrepresented sequences': mapping, 'Adapter Content': mapping,'Kmer Content': mapping})
preclean = df[['tot.seq','sample']]
df = df.drop('tot.seq', 1)														#
df = df.set_index('sample')														#
f, ax = plt.subplots()															#
hm = sns.heatmap(data = df, cmap="Pastel2", ax=ax, cbar=False)					#
box = ax.get_position()															#	
ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])					#	
legend_ax = f.add_axes([.7, .5, 1, .1])											#
legend_ax.axis('off')															#
colors = plt.cm.Pastel2(np.linspace(0, 1, len(mapping)))						#
patches = [mpatches.Patch(facecolor=c, edgecolor=c) for c in colors]			#
legend = legend_ax.legend(patches,												#
	{"PASS","WARN","FAIL"},														#
	handlelength=0.8, loc='lower left')											#
for t in legend.get_texts():													#	
	t.set_ha("left")															#
f.savefig('preqc.pdf', bbox_inches='tight')										#	
preclean = preclean.set_index('sample')
#------------------------------------------------------------------------------#


#BBDUK IT!!!##################################################################
#Values defined for forward and reserve trim will be used#####################
for i in reads:
	path1 = "in1=" + i + "_1.fastq.gz"
	path2 = "in2=" + i + "_2.fastq.gz"
	out1 = "out1=./readProcessed/" + i.split('/')[-1] + "_1.fastq.gz"
	out2 = "out2=./readProcessed/" + i.split('/')[-1] + "_2.fastq.gz"
	ftl = "ftl=" + clippingF
	subprocess.call(['./bbmap/bbduk.sh', ftl, path1, path2, out1, out2])
	path1 = "in1=./readProcessed/" + i.split('/')[-1] + "_1.fastq.gz"
	path2 = "in2=./readProcessed/" + i.split('/')[-1] + "_2.fastq.gz"
	out1 = "out1=./readProcessed/" + i.split('/')[-1] + "_cleaned_1.fastq.gz"
	out2 = "out2=./readProcessed/" + i.split('/')[-1] + "_cleaned_2.fastq.gz"
	report = i + ".txt"
	subprocess.call(['./bbmap/bbduk.sh', path1, path2, out1, out2, "k=31", "ref=./reference/mrrna.fasta", ">>", report])
	out1 = "./readProcessed/" + i.split('/')[-1] + "_cleaned_1.fastq.gz"
	out2 = "./readProcessed/" + i.split('/')[-1] + "_cleaned_2.fastq.gz"
	subprocess.call(['./FastQC/fastqc', out1, out2])
for i in reads:
 	out1 = "rm ./readProcessed/" + i.split('/')[-1] + "_1.fastq.gz"
	out2 = "rm ./readProcessed/" + i.split('/')[-1] + "_2.fastq.gz"
	stream = os.popen(out1)
	stream = os.popen(out1)
#----------------------------------------------------------------------------#


#Plot the heatmap from analysis###############################################
subprocess.call(['Rscript','postqc.r'])												
df = pd.read_csv('postq.csv')
df = df.groupby(['sample', 'tot.seq', 'module'])['status'].apply(', '.join).unstack().reset_index().rename_axis(None, axis=1)
mapping = {'PASS': 1, 'WARN': 2, 'FAIL': 3}
df = df.replace({'Basic Statistics': mapping, 'Per base sequence quality': mapping, 'Per tile sequence quality': mapping, 'Per sequence quality scores': mapping, 'Per base sequence content': mapping, 'Per sequence GC content': mapping, 'Per base N content': mapping, 'Sequence Length Distribution': mapping, 'Sequence Duplication Levels': mapping, 'Overrepresented sequences': mapping, 'Adapter Content': mapping,'Kmer Content': mapping})
postqclean = df[['tot.seq','sample']]
df = df.drop('tot.seq', 1)
df = df.set_index('sample')
f, ax = plt.subplots()
hm = sns.heatmap(data = df, cmap="Pastel2", ax=ax, cbar=False)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
legend_ax = f.add_axes([.7, .5, 1, .1])
legend_ax.axis('off')
colors = plt.cm.Pastel2(np.linspace(0, 1, len(mapping)))
patches = [mpatches.Patch(facecolor=c, edgecolor=c) for c in colors]
legend = legend_ax.legend(patches,
	{"PASS","WARN","FAIL"},
	handlelength=0.8, loc='lower left')
for t in legend.get_texts():
	t.set_ha("left")
f.savefig('postqc.pdf', bbox_inches='tight')
postqclean = postqclean.set_index('sample')
#----------------------------------------------------------------------------#

#Plot read charts`#############################################################
#Included in graphic core
#----------------------------------------------------------------------------#


#Align it using HISAT#########################################################
for i in reads:
	in1 = "./readProcessed/" + i.split('/')[-1] + "_cleaned_1.fastq.gz"
	in2 = "./readProcessed/" + i.split('/')[-1] + "_cleaned_2.fastq.gz"
	alignment = "./alignment/" + i.split('/')[-1]+ ".sam"
	alignmentbam = "./alignment/" + i.split('/')[-1]+ ".bam"
	splicesite = "./splice/" + i.split('/')[-1]+ ".txt"
	filename = "./splice/" + i.split('/')[-1]+ ".log"
	file = open(filename, 'a')
	subprocess.call(['./hisat/hisat2','-f','-x','./reference/genome_snp_tran', '-q', '-1', in1, '-2', in2, '-S', alignment, '--no-discordant', '--rna-strandness', direction, '--dta-cufflinks', '--novel-splicesite-outfile', splicesite, '-p', threads], stdout = file)
	file.close()
	subprocess.call(['samtools', 'view', '-b', '-S', alignment ,'-o', alignmentbam])
	remove = "rm " + "./alignment/" + i.split('/')[-1]+ ".sam"
	stream = os.popen(remove)
#----------------------------------------------------------------------------#

#Feature count################################################################
for i in reads:
	strand = 2
	alignmentbam = "./alignment/" + i.split('/')[-1] + ".bam"
	countfile = './counts/' + i.split('/')[-1] + '_count.txt'
	subprocess.call(['./subread/bin/featureCounts', '-T', threads, 's', strand,'-p', '-t', 'exon', '-g', 'gene_id', '-a', './annotation/gencode.vM14.annotation.gtf', '-R', '-o', countfile, alignmentbam])
#----------------------------------------------------------------------------#
