# process_amplicons

# Aim

This script checks flash overlapped amplicon reads (fastq format) the presence of both end primers and separate sequences using two tags (forward and reverse) then splits the sequences in Fasta files using the sample names found in the tag file. Then it extracts the most represented sequence and produces a statitics file.

# Input file formats :

fq, fastq, fq.gz or fastq.gz for the read file.

Text tab separated tag file with 3 columns : 
- sample name
- left tag (forward)
- right tag (reverse)

example :

<pre>
Lamal1	ACACACAC	GTGTGTGT
Hasi1	ACAGCACA	GTGTGTGT
Lamar1	GTGTACAT	GTGTGTGT
Euni1	TATGTCAG	GTGTGTGT
Lamo1	TAGTCGCA	GTGTGTGT
Lapu1	TACTATAC	GTGTGTGT
Anni1	ACTAGATC	GTGTGTGT
Anli1	GATCGCGA	GTGTGTGT
Lamal2	ACACACAC	TGTGCTGT
Hasi2	ACAGCACA	TGTGCTGT
</pre>

# running the script

<pre>
process_amplicons.py                                                                             
Usage: process_amplicons.py -r reads -t tags -o outdir

produces fasta files for the reads corresponding to each samples presented in
the tag file, plus a file for each most represented sequence in each sample
and a statitics file ex : process_amplicons.py -r reads -t tags

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit

  Options:
    -f READFILE, --read-file=READFILE
                        Indicate input read file in fastq (fq, fastq) gzipped
                        or not
    -t TAGFILE, --tags=TAGFILE
                        Indicate input tag text file with sample name, left
                        tag and reversed right tag separate by tabs
    -o OUTPUTDIRECTORY, --output-dir=OUTPUTDIRECTORY
                        Indicate the name of the output directory in which all
                        the samples files will be written
    -l LEFTPRIMER, --left-primer=LEFTPRIMER
                        Indicate the sequence of the left primer which will be
                        searched (exact) default : GACGAGAAGACCCTATA
    -r RIGHTPRIMER, --right-primer=RIGHTPRIMER
                        Indicate the sequence of the reversed right primer
                        which will be searched (exact) default :
                        GACCTCGATGTTGGATTAAGA

</pre>

# outputs 

The read fasta files are created in the output directory.
The best_seq read files correspond to the most represented sequence in the samples and are also created in the output directory
The "general_statistics.txt" file is created in the currrent directory and contains information of the reads in each sample

Text tab separated tag file with 3 columns : 
- sample name
- sample read count
- most represented sequence read count
- second most represented sequence read count
- read count histogram 

<pre>
Xnit1.fasta	2	1	1	1 1 
Xnit2.fasta	76	62	2	62 2 1 1 1 1 1 1 1 1 1 1 1 1 
Xnit3.fasta	195	147	2	147 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
Xsy1.fasta	186	128	2	128 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
Xsy2.fasta	201	149	2	149 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
Xsy3.fasta	175	124	2	124 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
Xsyl.fasta	89	58	8	58 8 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
Xte1.fasta	228	177	2	177 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
Xte2.fasta	100	69	2	69 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 

</pre>

# known issues 

If you have more than you ulimit number of samples you will have to change your ulimit (bash) using for example 
<pre>
ulimit -n 2048
</pre>
