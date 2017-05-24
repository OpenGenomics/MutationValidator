#Take in a MAF File, and then the following Bams with flags

import sys
import os
import pdb

os.system(". /broad/tools/scripts/useuse; use .zlib-1.2.6")

import csv
import pysam

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

parser = ArgumentParser(description='Process input bam files')
parser.add_argument('--mafsnp', type=str, help='Maf file (snp)')
parser.add_argument('--mafindel', type=str, help='Maf file (indel)')
parser.add_argument('--wextumor',type=str, help='WEX tumor bam file path', default=None)
parser.add_argument('--wexnormal',type=str, help='WEX normal bam file path', default=None)
parser.add_argument('--wgstumor',type=str, help='WGS tumor bam file path', default=None)
parser.add_argument('--wgsnormal',type=str, help='WGS normal bam file path', default=None)
parser.add_argument('--rnatumor',type=str, help='RNA tumor bam file path', default=None)
parser.add_argument('--targetedtumor',type=str, help='targeted tumor bam file path', default=None)
parser.add_argument('--targetednormal',type=str, help='targeted normal bam file path', default=None)
parser.add_argument('--lowpasstumor',type=str, help='low pass tumor bam file path', default=None)
parser.add_argument('--lowpassnormal',type=str, help='low pass normal bam file path', default=None)
parser.add_argument('--othertumor',type=str, help="other validation tumor bam file path", default=None)
parser.add_argument('--othernormal',type=str, help="other validation normal bam file path", default=None)
parser.add_argument('--out', type=str, help='Output file', default='output_pileup_file' )

args = parser.parse_args()

##########################################
#DEFINE FUNCTIONS:
##########################################

#Convert maf file to a list of positions:
def maf_to_position(filename): #convert maf to a position file with Chromosome, Start_position, End_position.
 	columns_care = ['Chromosome','Start_position','End_position']
 	count=0
 	pos = []
 	with open(filename, 'rb') as csvfile:
 		fullfile = csv.reader(csvfile, delimiter='\t')
 		for line in fullfile:
 			if line[0][0]!='#':
 				if count==0:
 					header=line
 					ind = [int(i) for i, x in enumerate(header) if x in columns_care]
 					if len(ind)!=3:
 						print('Maf file must contain the following headers: Chromosome, Start_position, End_position') #FIX THIS TO THROW ERROR!
 						exit()
 					count = 1
 				else:
 					if line[ind[0]]=='M':
						print('Removing all mitochondria sites')
					else:
						pos.append([line[i] for i in ind])
  	return pos

#Gernate bam soft links from an input bam file, tumor normal status, and the type of sequencing method: This is used becase pysam requires that the index format is .bam.bai whereas picard by default uses .bai
def create_soft_links(bamfile, sequencing, type_bam):
	os.system("mkdir -p ./softlinked/")
	if not os.path.exists(bamfile):
		print("Bam file does not exist! Please check your input files.") #TURN THIS INTO AN ERROR MESSAGE!!!
		exit()
	if os.path.exists(bamfile + ".bai"):
		print("Bam file already indexed appropriately for pysam")
		os.system("ln -s " + bamfile + " ./softlinked/" + sequencing + "." + type_bam + ".bam")
		os.system("ln -s " + bamfile + ".bai" + " ./softlinked/" + sequencing +"." + type_bam + ".bam.bai")
	elif os.path.exists(bamfile):
		os.system("ln -s " + bamfile + " ./softlinked/" + sequencing +  "." + type_bam + ".bam")
		os.system("ln -s " + bamfile[:-4] + ".bai ./softlinked/" + sequencing + "." + type_bam + ".bam.bai")

#Generate list of bam files from those that had symbolic link:
def dir_to_bam_list(dir_list):
    files = os.listdir(dir_list)
    fls = []
    for f in files:
        if f.endswith(".bam"):
            fls.append(os.path.join(dir_list,f))
    return fls
 

#######################################
#Functions to pull down the mutations:
def position_to_ins_del(read_start, position, cigar):
    #''' [(0, 87), (2, 1), (0, 62), (4, 1)] '''
    result = 0
    cigar_final = ""
    for tpl in cigar:
        cigar_final += str(tpl[0])*tpl[1]
    return cigar_final[(position-read_start)-1]

#Actual pileup:
def validate_mutation(bam_file,chrom,position):
    for pileupcolumn in bam_file.pileup(chrom, position-1, position):
        if pileupcolumn.pos == position-1:
            bases = []
            indel_counts = 0
            ins_counts = 0
            del_counts = 0
            for pileupread in pileupcolumn.pileups:
                if pileupread.alignment.is_duplicate or pileupread.alignment.mapq <=5 or pileupread.alignment.is_qcfail or pileupread.alignment.is_secondary or pileupread.alignment.is_unmapped:
                    #print "read failed quality metrics and therefore not counting"
                    continue
                bases.append(pileupread.alignment.seq[pileupread.query_position])
                if pileupread.indel:
                    indel_counts += 1
                cnt = position_to_ins_del(pileupread.alignment.pos,position,pileupread.alignment.cigar)
                if str(cnt) == "1":
                    ins_counts += 1
                if str(cnt) == "2":
                    del_counts += 1
            return [pileupcolumn.n,bases.count("A"),bases.count("T"),bases.count("C"),bases.count("G"), ins_counts,del_counts]

#writes out the line in the preprocessing file based on the pileup for each position.
def validate_indel_mut_bam_position(sample,type_bam,bam,position_list,output_file):
    # pileup for indels and point mutations:
    for pos in position_list:
        samfile = pysam.Samfile(bam,"rb")
        if sample=='rna':
            chrom = 'chr' + pos[0]
            pileup = validate_mutation(samfile, chrom, int(pos[2]))
        else:
            pileup = validate_mutation(samfile, pos[0], int(pos[2]))
        if pileup==None:
            pileup = ['0']*7
        samfile.close()
        output_file.write(sample + "\t" + type_bam + "\t" + str(pos[0]) + "\t" + str(pos[1]) + "\t" + str(pos[2]) + "\t" + "\t".join([str(x) for x in pileup]) + "\n")

##########################################################################
###################################
#SCRIPT to EXECUTE on INPUT FILES:
####################################

#Convert input maf files to position lists:
if args.mafsnp != None and args.mafindel == None:
	positions = maf_to_position(args.mafsnp)
if args.mafindel == None and args.mafindel != None:
	positions = maf_to_position(arg.mafindel)
if args.mafsnp != None and args.mafindel !=None:
	posmaf =  maf_to_position(args.mafsnp)
	posindel =  maf_to_position(args.mafindel)
	positions = posmaf + posindel

#Generate symbolic links to the bam files (required due to pysam .bam.bai formatting requirement)
if args.wextumor != None:
	create_soft_links(args.wextumor,'wex','tumor')
if args.wexnormal != None:
	create_soft_links(args.wexnormal,'wex','normal')
if args.wgstumor != None:
	create_soft_links(args.wgstumor,'wgs','tumor')
if args.wgsnormal != None:
	create_soft_links(args.wgsnormal,'wgs','normal')
if args.rnatumor != None:
	create_soft_links(args.rnatumor,'rna','tumor')
if args.targetedtumor != None:
	create_soft_links(args.targetedtumor,'targeted','tumor')
if args.targetednormal != None:
	create_soft_links(args.targetednormal,'targeted','normal')
if args.lowpasstumor != None:
	create_soft_links(args.lowpasstumor, 'lowpass','tumor')
if args.lowpassnormal != None:
	create_soft_links(args.lowpassnormal,'lowpass','normal')
if args.othertumor !=None:
	create_soft_links(args.othertumor,'other','tumor')
if args.othernormal !=None:
	create_soft_links(args.othernormal,'other','normal')

#list to pull down info on the softlinked bam files (preprocessed earlier)
files = dir_to_bam_list("./softlinked/")
#if len(files) =<2:
#	print('There are no validation bams in this set. Therefore running this script seems useless')
#	exit()

#Create a preprocessing file to write out pileups, then determine pileup info for each bam
output_file = args.out + ".pileup_preprocessing.txt"
results = open(output_file,"w")
header = ["Sample","Type","Chromosome", "Start_position", "End_position", "ref", "A","T","C","G","ins","del"]
results.write("\t".join([str(x) for x in header]) + "\n")

if len(positions)==0:
    print('Warning: You are passing in an empty maf. Only a header will be present in the pileup. You might want to check your input files.')
else:
    for fl in files:
        validate_indel_mut_bam_position(os.path.basename(fl).strip().split(".")[0],os.path.basename(fl).strip().split(".")[1],fl,positions,results)

results.close()

os.system("rm -rf ./softlinked/")

#########################################
#EXAMPLE TO IMPLEMENT: THCA

#python validation_wrapper.py --mafsnp /xchip/cga/gdac-prod/cga/jobResults/filter_tsv_by_field/THCA-TCGA-BJ-A191-TP-NB-SM-1WK7Q-SM-1WK7S/4771644/THCA-TCGA-BJ-A191-TP-NB-SM-1WK7Q-SM-1WK7S.cross_center.SNP.maf --mafindel /xchip/cga/gdac-prod/genepattern/jobResults/6652540/THCA-TCGA-BJ-A191-Tumor-SM-1WK7Q.maf.annotated	--wextumor /seq/picard_aggregation/C499/TCGA-BJ-A191-01A-11D-A13W-08/v3/TCGA-BJ-A191-01A-11D-A13W-08.bam --wexnormal /seq/picard_aggregation/C499/TCGA-BJ-A191-10A-01D-A13W-08/v3/TCGA-BJ-A191-10A-01D-A13W-08.bam --wgstumor /seq/picard_aggregation/G32522/TCGA-BJ-A191-01A-11D-A13W-08/v1/TCGA-BJ-A191-01A-11D-A13W-08.bam --wgsnormal /seq/picard_aggregation/G32522/TCGA-BJ-A191-10A-01D-A13W-08/v1/TCGA-BJ-A191-10A-01D-A13W-08.bam --rnatumor /xchip/cga/gdac-prod/genepattern/jobResults/7105183/THCA-TCGA-BJ-A191-Tumor-SM-1WK7Q.recalibrated.bam


#python validation_wrapper.py -mafsnp /xchip/cga/gdac-prod/cga/jobResults/Oncotator_v1/THCA-TCGA-BJ-A0YZ-TP-NB-SM-1RAEJ-SM-1RAEM/8961784/THCA-TCGA-BJ-A0YZ-TP-NB-SM-1RAEJ-SM-1RAEM.snp.capture.maf.annotated --mafindel /xchip/cga/gdac-prod/cga/jobResults/Oncotator_v1/THCA-TCGA-BJ-A0YZ-TP-NB-SM-1RAEJ-SM-1RAEM/8961954/THCA-TCGA-BJ-A0YZ-TP-NB-SM-1RAEJ-SM-1RAEM.indel.capture.maf.annotated --wextumor /seq/picard_aggregation/C499/TCGA-BJ-A0YZ-01A-11D-A10S-08/v5/TCGA-BJ-A0YZ-01A-11D-A10S-08.bam  --wexnormal /seq/picard_aggregation/C499/TCGA-BJ-A0YZ-10A-01D-A10S-08/v5/TCGA-BJ-A0YZ-10A-01D-A10S-08.bam --rnatumor /cgaext/tcga/links/RNA-Seq/THCA/TCGA-BJ-A0YZ-01A-11R-A10U-07/UNC-LCCC_ILLUMINA_HG19/2012-06-09T13:00:24Z/5b6dd14c-538d-496e-b1b0-3c548a2f5eff/UNCID_1241222.1ba7877e-fd4e-452d-8eed-20fc9fb39698.sorted_genome_alignments.bam

##################################################
