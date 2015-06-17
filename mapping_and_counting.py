__author__ = 'benkjohnson'


import os
import subprocess
import glob
from shutil import copy
import check_dependencies_windows
import qc_analysis

class Mapping_and_Counting(object):
    def __init__(self):

        return

    def bowtie(self, datalocation, analysislocation):
        """Run Bowtie for SE reads less than 50 bp in length.
        Will add the ability to run Bowtie2 for PE and SE with
        reads greater than 50 bp."""

        cd = check_dependencies_windows.CheckDependencies()
        qc = qc_analysis.QC_analysis()
        gff, genref = qc.findreferencefiles(datalocation)
        copy(genref, os.path.join(analysislocation, 'Bowtie'))
        copy(gff, os.path.join(analysislocation, 'HTSeq'))
        # subprocess.Popen("cp " + genref + " " + analysislocation + "/Bowtie", shell=True).wait()
        # subprocess.Popen("cp " + gff + " " + analysislocation + "/HTSeq", shell=True).wait()
        os.chdir(os.path.join(cd.getSPARTAdir(), "Mapping_and_counting"))
        # os.chdir(cd.getSPARTAdir() + "/Mapping_and_counting")
        # if not os.path.lexists(os.path.join(cd.getSPARTAdir(), "Mapping_and_counting", "bowtie-1.1.1")):
        #     #This will be a problem for Windows users. Distribute with unzipped binaries?
        #     subprocess.call(["unzip", "bowtie-1.1.1-linux-x86_64.zip"], stdout=open(os.devnull, 'wb'))
        os.chdir(os.path.join(cd.getpwd(), "bowtie-1.1.1"))
        for file in os.listdir(os.path.join(analysislocation, "QC")):
            extension = file.split(".")[-1]
            if extension == "gz":
                subprocess.Popen("gzip -dc " + os.path.join(analysislocation, "QC", file) + " > " + os.path.join(analysislocation, "Bowtie", os.path.splitext(file)[0]), shell=True).wait()
        print "Building the Bowtie index from the reference genome"
        #This could be a problem... Have the user specify if the index is large or small
        subprocess.Popen(".\\bowtie-build-s.exe -q -f " + glob.glob(os.path.join(analysislocation, "Bowtie") + "/*.fa*")[0] + " " + glob.glob(os.path.join(analysislocation, "Bowtie") + "/*.fa*")[0].split(".")[0], shell=True).wait()
        allebwtfiles = glob.glob("*.ebwt")[:]
        for ebwtfile in allebwtfiles:
            copy(ebwtfile, os.path.join(analysislocation, "Bowtie"))
            # subprocess.Popen("cp " + ebwtfile + " " + analysislocation + "/Bowtie/", shell=True).wait()
        print "Mapping reads to the reference genome with Bowtie"
        for file in os.listdir(os.path.join(analysislocation, "Bowtie")):
            extension = os.path.splitext(file)[1]
            if extension == ".fq" or extension == ".fastq":
                fname = os.path.splitext(file)[0]
                strippedfile = fname[len('trimmed'):]
                #Again this could be a problem, have the user define whether the index is large or small
                subprocess.Popen(".\\bowtie-align-s.exe -S " + glob.glob(os.path.join(analysislocation, "Bowtie") + "/*.fa*")[0].split(".")[0] + " " + os.path.join(analysislocation, "Bowtie", file) + " > " + os.path.join(analysislocation, "Bowtie", "align" + strippedfile + ".sam"), shell=True).wait()
            # elif extension == ".fastq":
            #     fname = os.path.splitext(file)[0]
            #     strippedfile = fname[len('trimmed'):]
            #     subprocess.Popen("./bowtie -S " + analysislocation + "/Bowtie/trimmedMtbCDC1551" + " " + analysislocation + "/Bowtie/" + file + " > " + analysislocation + "/Bowtie/align" + strippedfile + ".sam", shell=True).wait()


    def htseq(self, analysislocation):
        """Run htseq-count to count gene features post-Bowtie mapping"""

        cd = check_dependencies_windows.CheckDependencies()
        os.chdir(os.path.join(cd.getSPARTAdir(), "Mapping_and_counting"))
        # if not os.path.lexists(os.path.join(cd.getSPARTAdir(), "Mapping_and_counting", "HTSeq-0.6.1")):
        #     subprocess.Popen("tar -zxf HTSeq-0.6.1.tar.gz", stdout=open(os.devnull, 'w'), shell=True).wait()
        # htseqcheck = cd.checkhtseq()
        # if htseqcheck == False:
        # os.chdir(os.path.join(cd.getpwd(), "HTSeq-0.6.1"))
        # subprocess.Popen("python setup.py build install --user", shell=True).wait()
        # os.chdir(os.path.join(cd.getpwd(), "build", "scripts-2.7"))
        gff = glob.glob(os.path.join(analysislocation, "HTSeq") + "/*.g*")[0]
        print "Counting gene features with HTSeq"
        for mapfile in os.listdir(os.path.join(analysislocation, "Bowtie")):
            extension = os.path.splitext(mapfile)[1]
            if extension == ".sam":
                fname = os.path.splitext(mapfile)[0]
                strippedmapfile = fname[len('align'):]
                subprocess.Popen("python -m HTSeq.scripts.count -m intersection-nonempty --stranded=yes " + os.path.join(analysislocation, "Bowtie", mapfile) + " " + gff + " > " + os.path.join(analysislocation, "HTSeq", "map" + strippedmapfile + ".sam"), shell=True).wait()

                # if htseqcheck == True:
                #     try:
                #         #Need to add the ability to change options/flags
                #         #Also, there is a bug here... Need a way to catch the output and check if /bin/sh htseq-count doesn't exist
                #         subprocess.Popen("htseq-count -m intersection-nonempty --stranded=yes " + os.path.join(analysislocation, "Bowtie", mapfile) + " " + gff + " > " + os.path.join(analysislocation, "HTSeq", "map" + strippedmapfile + ".sam"), shell=True).wait()
                #     except Exception:
                #         print "HTSeq appears to be installed, but htseq-count is not executable."
                #         print "Attempting to build HTSeq from source and execute htseq-count from there."
                #         os.chdir(os.path.join(cd.getpwd(), "HTSeq-0.6.1"))
                #         subprocess.Popen("python setup.py build", shell=True).wait()
                #         os.chdir(os.path.join(cd.getpwd(), "build", "scripts-2.7"))
                #         subprocess.Popen("./htseq-count -m intersection-nonempty --stranded=yes " + os.path.join(analysislocation, "Bowtie", mapfile) + " " + gff + " > " + os.path.join(analysislocation, "HTSeq", "map" + strippedmapfile + ".sam"), shell=True).wait()
                # else:
                #     subprocess.Popen("./htseq-count -m intersection-nonempty --stranded=yes " + os.path.join(analysislocation, "Bowtie", mapfile) + " " + gff + " > " + os.path.join(analysislocation, "HTSeq", "map" + strippedmapfile + ".sam"), shell=True).wait()

        return