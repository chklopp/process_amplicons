#!env python
# -*- coding: utf-8 -*-
#
# process_amplicons.py
# Copyright (C) 2021 INRA
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

__author__ = 'Christophe Klopp'
__copyright__ = 'Copyright (C) 2021 INRA'
__license__ = 'GNU General Public License'
__version__ = '0.1'
__email__ = 'christophe.klopp@inrae.fr'
__status__ = 'beta'

# from ConfigParser import ConfigParser
from optparse import *
import os, sys, gzip
#import subprocess

def version_string ():
    """
    Return process_amplicons.py version
    """
    return "process_amplicons.py " + __version__

def check_output_directory(od) :
    if not os.path.isdir(od) :
        os.mkdir(od)

def load_tags(tagfile) :
    fw = []
    rv = []
    tags = {}
    if os.path.exists(tagfile):
        with open(tagfile, 'r') as fi:
            try:
                for l in fi :
                    #print(l)
                    b = l[:-1].split("\t")
                    if len(b[1]) > 1 and len(b[2]) > 1 and len(b[0]) > 1 :
                        fw.append(b[1])
                        rv.append(b[2])
                        tags[b[1]+"_"+b[2]] = b[0]
            except :
                sys.stderr.write("#### can not open "+tagfile+" file\n")
    # print(fw); print(rv); print(tags)
    return(fw,rv,tags)

def rev_comp(s) :
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return("".join(complement.get(base, base) for base in reversed(s)))
    
def find_primers(s, p):
    # print("find primer :",s)
    fleft = s.find(p[0])
    fright = s.find(p[1])
    rleft = rev_comp(s).find(p[0])
    rright = rev_comp(s).find(p[1])
    #print(fleft+len(p[0])+1, fright, rleft+len(p[0])+1, rright)
    if fleft > 0 and fright > 0 :
            # extract 16 region 
            sr = s[fleft+len(p[0])+1:fright]
            # extract both flancking region to search tags
            lsr = s[0:fleft]
            rsr = s[fright:len(s)] # rsr = s[fright+len(p[1]):len(s)]
            return(("+", fleft+len(p[0])+1, fright, sr, lsr, rsr))
    elif rleft > 0 and rright > 0 :
            # extract 16 region
            lsr = rev_comp(s)[0:rleft]
            # extract both flancking region to search tags
            rsr = rev_comp(s)[rright:len(s)] # rsr = rev_comp(s)[rright+len(p[1]):len(s)]
            sr = rev_comp(s)[rleft+len(p[0])+1:rright]
            return(("-", rleft+len(p[0])+1, rright, sr, lsr, rsr))
    else :
        return("=",0,0,s,"","")


def find_tags(l, strand, lseq, rseq, fw, rv, tags):
    # print(l, strand, lseq, rseq)
    lfound = ""
    lcount = 0
    rfound = ""
    rcount = 0
    for i in fw :
        if lseq.find(i) > 0 :
            lcount += 1
            lfound = i
    for j in rv :
        if rseq.find(j) > 0 :
            rcount += 1
            rfound = j
    # print("find tags ",lfound,rfound,lcount,rcount)
    if lcount >= 1 and rcount >= 1 :
        if lfound+"_"+rfound in tags.keys() :
            return(tags[lfound+"_"+rfound])
        else :
            return("None")
    else :
        return("None")

    
def print_seq_to_file(outdir, sname, strand, left, right, seq, sample, fhandler):
    #out = outdir+"/"+sample+".fasta"
    #f = open(out, 'a')
    fhandler[sample].write(">"+sname+" "+strand+" "+str(left)+" "+str(right)+"\n")
    fhandler[sample].write(seq+"\n")
    #f.close()

def  open_output_files(outdir, tags) :
    # create a dictionnary of file handlers
    ftable = {}
    for k,v in tags.items() :
        ftable[v] = open(outdir+"/"+v+".fasta","w")    
    return(ftable)
        
def  close_output_files(ftable) :
    for k,v in ftable.items() :
        v.close()
    
def process_file(outdir, f, fw, rv, tags, fhandler):
    extension = os.path.splitext(f)[1]
    if extension != ".gz" and extension != ".fastq" and extension != ".fq" :
        sys.stderr.write("#### file extension not recognized "+extension+" \n")
        exit()
    else :
        line_count=0
        seq_name = "" 
        if os.path.exists(f):
            if extension == ".gz" :
            # print("file : "+f)
                fi = gzip.open(f, 'rt')
            else :
                fi = open(f, 'rt')
                
            for l in fi :
                # print(l)
                line_count += 1
                # get seqname
                if line_count%4 == 1:
                    seq_name = l[1:-1].split(" ")[0]
                    # print("seq name : "+seq_name)
                if line_count%4 == 2:
                    #print(l)
                    (strand, left, right, seq, lseq, rseq) = find_primers(l, primers)
                    # print(strand, left, right, seq)
                    # the sequence was found then it should be attributed to the correct sample 
                    if strand != "=" :
                        sample = "None"
                        sample = find_tags(seq, strand, lseq, rseq, fw, rv, tags)
                        # print("sample : "+sample )
                        # the sample found the seqeunces should be printed 
                        if sample != "None" :
                            print_seq_to_file(outdir, seq_name, strand, left, right, seq, sample, fhandler)
                    else :
                        continue
             
    sys.stderr.write("#### "+str(int(line_count/4))+" sequences processed\n")
             
    return None

def generate_stats_and_best_seq_files(outd):
    # read all fasta files 
    files = [f for f in os.listdir(outd) if f.endswith('.fasta')]
    fhso = open("general_statistics.txt","w")
    samples = {} # data structure to validate samples (sample_name + sample version in one char)
    for f in files :
        sample_name = os.path.splitext(f)[0]
        species_name = sample_name[0:-1]
        # print("sample name ",sample_name)
        seq = {}
        fh = open(outd+"/"+f,"r")
        # fill dictonnary of sequences 
        if os.stat(outd+"/"+f).st_size != 0 :
            fho = open(outd+"/"+f+".best_seq","w")
            nbseq = 0
            for l in fh :
                #print(l)general_statistics.txt
                if l[0:1] != ">" :
                    nbseq += 1
                    nuc = l[0:-1]
                    if nuc in seq.keys() :
                        seq[nuc] += 1
                    else :
                        seq[nuc] = 1
            # find element most represented in dictionnary
            #print(seq)
            sorted_seq = sorted(seq.items(), key=lambda kv: kv[1])
            # print(f, sorted_seq)
            fho.write(">ref_seq_"+f+"\n")
            fho.write(str(sorted_seq[-1][0])+"\n")
            fho.close()
            
            stat = ""
            # create stat line with all the counts for the different sequences found in the sample
            for k,v in sorted_seq :
                stat = str(v)+" "+stat
                
            # writing the stat line 
            if len(sorted_seq) > 1 :
                # print(sorted_seq, len(sorted_seq), sorted_seq[len(sorted_seq)][1])
                stat = f+"\t"+str(nbseq)+"\t"+str(sorted_seq[-1][1])+"\t"+str(sorted_seq[-2][1])+"\t"+stat+"\n"
            else :
                stat = f+"\t"+str(nbseq)+"\t"+str(sorted_seq[-1][1])+"\tO\t"+stat+"\n"
            fhso.write(stat)
            
            # writing the sample information
            if species_name in samples :
                #print("append ",species_name, sample_name, samples[species_name])
                samples[species_name].append([sample_name, str(sorted_seq[-1][0]),sorted_seq[-1][1]])
            else :
                #print("create ", species_name, sample_name)
                samples[species_name] = [[sample_name, str(sorted_seq[-1][0]),sorted_seq[-1][1]]]
                #print("after create ",sample_name, samples[species_name])
                
    fhso.close()
    return(samples)
   
def check_samples(samples):
    for k,v in samples.items() :
        #print(k,v)
        # Check sequence correspondence 
        seq = []
        for i,s in enumerate(v) :
            # print(i, s)
            seq.append(s[1])
            
        sseq = set(seq)
        if len(sseq) == 1 :
            print("OK", k, str(i+1), ",".join(sseq), sep="\t")
        else :
            for i,s in enumerate(v) :
                print("KO", k, str(i+1), s[0], s[1],  str(s[2]), sep="\t")
            
    

if __name__ == "__main__":

    parser = OptionParser(usage="Usage: process_amplicons.py -r READS -t TAGS -o OUTDIR")

    usage = "usage: %prog -r reads -t tags -o outdir"
    desc = "produces fasta files for the reads corresponding to each samples presented in the tag file, plus a file for each most represented sequence in each sample and a statitics file\n"\
           "ex : process_amplicons.py -r reads -t tags"
    parser = OptionParser(usage = usage, version = version_string(), description = desc)
    
    ogroup = OptionGroup(parser, "Options","")
    ogroup.add_option("-f", "--read-file", dest="readFile",
                      help="Indicate input read file in fastq (fq, fastq) gzipped or not",  type="string")
    ogroup.add_option("-t", "--tags", dest="tagFile",
                      help="Indicate input tag text file with sample name, left tag and reversed right tag separate by tabs",  type="string")
    ogroup.add_option("-o", "--output-dir", dest="outputDirectory",
                      help="Indicate the name of the output directory in which all the samples files will be written",  type="string", default="Samples")
    ogroup.add_option("-l", "--left-primer", dest="leftPrimer",
                      help="Indicate the sequence of the left primer which will be searched (exact) default : GACGAGAAGACCCTATA",  type="string", default="GACGAGAAGACCCTATA")
    ogroup.add_option("-r", "--right-primer", dest="rightPrimer",
                      help="Indicate the sequence of the reversed right primer which will be searched (exact) default : GACCTCGATGTTGGATTAAGA",  type="string", default="GACCTCGATGTTGGATTAAGA")
    parser.add_option_group(ogroup)
    
    (options, args) = parser.parse_args()
    if options.readFile == None and  options.tagFile == None :
        parser.print_help()
        sys.exit(1)
    else:
        # checks if output dir exists and creates it if not 
        primers = (options.leftPrimer,options.rightPrimer)
        check_output_directory(options.outputDirectory)
        # loading tags in a dictionnary
        (foward, reverse, tags) = load_tags(options.tagFile)
        fh = open_output_files(options.outputDirectory, tags)
        #print("fg1 :", fg1)
        # checking second genome file extension if exists
        process_file(options.outputDirectory, options.readFile, foward, reverse, tags, fh)
        close_output_files(fh)
        samples = generate_stats_and_best_seq_files(options.outputDirectory)
        check_samples(samples)
