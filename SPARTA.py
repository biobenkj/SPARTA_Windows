__author__ = 'benkjohnson'

import check_dependencies_linux
import optparse
import sys


cd = check_dependencies_linux.CheckDependencies()


optParser = optparse.OptionParser(
  usage="python SPARTA.py [options]",

  description="Simple Program for Automated reference-based bacterial RNA-seq Transcriptome Analysis (SPARTA)",

  epilog="Written by Benjamin K. Johnson (john3434@msu.edu), Michigan State University Department of " +
  "Microbiology and Molecular Genetics. (c) 2015")

optParser.add_option("--SE", help="Single-end read input. Default input choice is single-end if nothing is specified",
                    action="store_true", default="True", dest="seqtype")
optParser.add_option("--PE", help="Paired-end read input. Must have the exact same file name and end with _F for the forward read and _R for the reverse read",
                    action="store_false", default="False", dest="seqtype")
optParser.add_option("--cleanup", help="Clean up the intermediate files to save space. Default action is to retain the intermediate files. Usage: --cleanup=True",
                    action="store", default="False", dest="cleanup")


trim = optparse.OptionGroup(optParser, 'Trimmomatic options', "The order the options will be run are: ILLUMINACLIP, LEADING, TRAILING, SLIDINGWINDOW, MINLEN")
trim.add_option("--clip", help="ILLUMINACLIP options. MiSeq & HiSeq usually TruSeq3.fa; GAII usually TruSeq2.fa. Default is ILLUMINACLIP:TruSeq3-SE.fa:2:30:10. Usage: --clip=<adapterseqs>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>",
                  action="store", default="TruSeq3-SE.fa:2:30:10", dest="illuminaclip")
trim.add_option("--lead", help="Set the minimun quality required to keep a base. Default is LEADING=3. Usage: --lead=<quality>",
                  action="store", default=3, dest="leading")
trim.add_option("--trail", help="Set the minimum quality required to keep a base. Default is TRAILING=3. Usage: --trail=<quality>",
                  action="store", type="int", default=3, dest="trailing")
trim.add_option("--slidewin", help="SLIDINGWINDOW options. Default is SLIDINGWINDOW:4:15. Usage: --slidewin=<window_size>:<required_quality>",
                  action="store", default="4:15", dest="slidingwindow")

#bowtie = optparse.OptionGroup(optParser, 'Bowtie options for single-end reads')
#bowtie.add_option("")

htseq = optparse.OptionGroup(optParser, 'HTSeq options')
htseq.add_option("--stranded", help="Stranded options: yes, no, reverse. Default is --stranded=reverse. Usage: --stranded=yes/no/reverse",
                   action="store", default="reverse", dest="stranded")
htseq.add_option("--order", help="Order options: name, pos. Usage: --order=name/pos.", action="store", dest="order")
htseq.add_option("--minqual", help="Skip all reads with quality lower than the given value. Default is --minqual=10. Usage: --minqual=<value>",
                   action="store", type="int", default=10, dest="minqual")
htseq.add_option("--idattr", help="Feature ID from the GTF file to identify counts in the output table Default is --idattr=gene_id. Usage: --idattr=<id attribute>",
                   action="store", default="gene_id", dest="idattr")
htseq.add_option("--mode", help="Mode to handle reads overlapping more than one feature. Default is --mode=union. Usage: --mode=union/intersection-strict/intersection-nonempty",
                   action="store", default="union", dest="mode")

optParser.add_option_group(trim)
# optParser.add_option_group(bowtie)
optParser.add_option_group(htseq)
(options, args) = optParser.parse_args()

# args = parser.parse_args()
#
# if args.clip:
#     print args.clip

#Welcome the user to the software and check dependencies

print "Welcome to SPARTA!"
print "Let's make sure we have everything we need to get started..."
print "Now checking dependencies..."

#Check for Java, R, and NumPy

javacheck = cd.checkjava()
Rcheck = cd.checkR()
numpycheck = cd.checknumpy()
# matplotlibcheck = cd.checkmatplotlib()

#If NumPy can't be found, SPARTA will attempt to download and install it

if numpycheck == False:
    cd.installdependencies()
    answer = cd.getanswerstate()
    if answer:
        if numpycheck == False:
            cd.getNumPy()
            cd.installNumPy()
        # elif numpycheck == False and matplotlibcheck == True:
        #     cd.getNumPy()
        #     cd.installNumPy()
        # elif numpycheck == True and matplotlibcheck == False:
        #     cd.getmatplotlib()
        #     cd.installmatplotlib()

print "Everything appears to be fine. Moving on.\n"

import qc_analysis
import mapping_and_counting
import differential_expression

qc = qc_analysis.QC_analysis()
mac = mapping_and_counting.Mapping_and_Counting()
de = differential_expression.DifferentialExpression()

#Create a folder called RNAseq_Data with which to store all of the data analysis

subfolderpath = qc.create_folder()

#Find the RNAseq data folder location

rawdatapath = qc.finddata()

#Read in the data and sort out the genome feature file and genome sequence files

qc.findreferencefiles(rawdatapath)

#Trimmomatic

qc.trimmomatic(rawdatapath, subfolderpath)

#FastQC

qc.fastqc(rawdatapath, subfolderpath)

#Bowtie

mac.bowtie(rawdatapath, subfolderpath)

#HTSeq

mac.htseq(subfolderpath)

#edgeR

de.de_analysis(subfolderpath)











