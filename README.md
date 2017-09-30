# biopython_orf_find
Python script utilised for identifying open reading frames in a study of Calliactis polypus regeneration (Transcriptomic investigation of wound healing and regeneration in the cnidarian Calliactis polypus, Scientific Reports, doi: 10.1038/srep41458) (http://www.nature.com/articles/srep41458)

# Dependencies
This script was designed to work with Python 3, and utilises the ‘Biopython’ package (https://github.com/biopython/biopython.github.io/). This script has been tested on Windows and Linux/SUSE, but it should work anywhere that Python does.

# Description of script logic
This script was designed to be usable by those unfamiliar with command-line operations as well as by experienced users familiar with Python operations. Thus, the script can accept arguments on the command-line, or it can (on a Windows environment) be double-clicked to launch an interactive console window with text prompts which specify to the user what commands are required at each point, with checks in place to ensure the user inputs the correct values. The order of this is to...

1. Specify the name of the fasta file which contains the nucleotide sequences from which ORFs will be extracted.
2. Specify the output file name containing the extracted ORFs.
3. Specify the minimum ORF length you wish to consider. 
4. Specify the maximum ORF length you wish to consider (can be unlimited). 
5. Specify the number of ORFs you wish to obtain from each nucleotide sequence which meet this length requirement.
6. Specify two stringency values which will determine the weighting with which we will consider ORFs with alternative (i.e., TTG, GTG, CTG) or no-codon (i.e., fragmented sequence) starts as opposed to traditional.
7. Specify what format ORFs should be presented as (i.e., protein translated, nucleotide CDS, or both).
8. Optional ability to replace alternative start codons with methionine ('M', only relevant if protein translated ORFs are being obtained).

Before delving into the specifics of how the stringencies (step 6) work, it should first be mentioned that this script works on the basis of identifying regions in-between stop codons. Thus, to this script, an ORF is any region uninterrupted by stop codons. Returning to the stringency values, these values have defaults which I recommend the script runs with, but if shorter peptides (such as those of 10-50AA length) which may commonly have alternative start codons are sought, then changing the stringency of these default values manually is a valid option.

###Alternative codon###
The default alternative codon parameter is 49. This means that the script will consider an ORF that starts with an alternative codon as "better" than one that starts with a methionine **only** if it is greater than 49 AA longer. 

###No-codon###
The default no-codon parameter is 99. This means that the script will consider an ORF that does not start with a codon as "better" than one that starts with **any** codon (traditional or alternative) **only** if it is greater than 99 AA longer.

These two stringencies affect the internal sorting process of the script, and determine what order ORFs are presented in the output files. Internally, when looking at any individual ORF, it will decide if a traditional start codon, an alternative start codon, or no codon best fits the ORF. The implicit assumption is that a no-codon start is a fragmentary ORF, and this is why it should be weighted against most heavily. In most scenarios, a methionine codon will be present in an ORF, so the two stringency values help to decide whether the ORF should start at the first methionine, or if it should start earlier. With regards to the output, the script will rank all ORFs obtained from a sequence using the two stringencies. Thus, the first ORF for each sequence is considered most likely to be the "best." Each subsequent ORF will be a bit shorter, or it may have an alternative or no-codon start which is weighted against.

# File in- and output
This script will read in fasta-formatted files containing nucleotide sequences. The output will be fasta-formatted file(s) containing protein translated ORFs, nucleotide CDS sequences, or both forms of output can be generated. The original sequence identifiers will be modified in this output to contain the ORF number as determined from this script. For example, if an original nucleotide sequence is titled ‘>contig1’, depending on the number of ORFs identified in this sequence, the output file will have entries titled ‘>contig1_ORF1’ and ‘>contig1_ORF2’, etc.

# Additional notes
This script does not require much RAM, and thus should be suitable for use on all types of computers. Unless your computer's processor is very weak, this script should be capable of processing files with hundreds of thousands of sequences in time spans of less than 10 minutes (approximately), though depending on certain parameter configurations this time can vary to some degree. As this script provides a progress bar, it can be roughly gauged how long the script should take to complete.

More complex ORF finders may often consider things such as GC content and the presence of Kozak consensus sequences among other features. While this script does not offer this, operating solely on the basis of ORF length, through personal testing I have found it to provide results which are more reliable than NCBI’s ORF Finder. Due to the ability to determine the strictness with which we consider alternative starts, the script is designed to be suitable for finding novel ORFs wherein assumptions of GC content and other sequence features may not hold. Additionally, as this script is capable of pulling many ORFs out of a sequence, it is also intended for performing analyses such as the one in this study, wherein multiple transcriptomes had potential ORFs extracted and compared via BLAST to identify conserved regions. Subsequently, as mentioned, this script is designed primarily with novel ORF identification in mind. If you intend to use this for yourself, you may want to consider what your goals are, as this script is not necessarily designed to find the most biologically valid start codon of conserved genes which typically demonstrate certain sequence features (although it tries). If you do find this script useful in any studies you perform, I’d appreciate if you cite the publication this script is associated with, and feel free to contact me if you have any questions.

Finally, the script has detailed usage details when called on the command-line. This is presented below. 

```
usage: biopython_orf_find.py [-h] [-i FILENAME] [-o OUTPUTFILENAME]
                             [-min MINPROLEN] [-max MAXPROLEN]
                             [-num HITSTOPULL] [-alt ALTCODONSTRINGENCY]
                             [-no NOCODONSTRINGENCY]
                             [-st {prot,nucl,both,PROT,NUCL,BOTH}]
                             [-r {y,n,Y,N}] [-f {y,n,Y,N}]

biopython_orf_find.py reads in a fasta formatted file containing nucleotide
sequences and, following user-specified parameters, produces an output fasta
file containing potential open reading frames (ORFs) as nucleotides/protein
translations/both.

optional arguments:
  -h, --help            show this help message and exit
  -i FILENAME, -input FILENAME
                        Input fasta file name
  -o OUTPUTFILENAME, -output OUTPUTFILENAME
                        Output fasta file name
  -min MINPROLEN, -minimum MINPROLEN
                        Minimum ORF amino acid length. Default == 30.
  -max MAXPROLEN, -maximum MAXPROLEN
                        Optional specification for maximum ORF amino acid
                        length. Default == 0, which means there is no upper
                        limit.
  -num HITSTOPULL, -numhits HITSTOPULL
                        Specify the number of ORFs you wish to extract from
                        each sequence. Default == 3.
  -alt ALTCODONSTRINGENCY, -altcodon ALTCODONSTRINGENCY
                        Control the stringency with which alternative start
                        codon ORFs are accepted. Recommended not to change
                        unless you understand the influence this has. Default
                        == 49.
  -no NOCODONSTRINGENCY, -nocodon NOCODONSTRINGENCY
                        Control the stringency with which fragmentary ORFs are
                        accepted (fragmentary means there is no traditional or
                        common alternative start codon in the sequence).
                        Recommended not to change unless you understand the
                        influence this has. Default == 99.
  -st {prot,nucl,both,PROT,NUCL,BOTH}, -seqtype {prot,nucl,both,PROT,NUCL,BOTH}
                        Specify the type of output you want to generate (i.e.,
                        protein translated ORF, nucleotide CDS, or both). If
                        you specify 'both', two outputs with '_prot' and
                        '_nucl' suffix will be generated. Default == 'prot'.
  -r {y,n,Y,N}, -replace {y,n,Y,N}
                        Optional ability to replace alternative starting
                        position with a methionine (M) [only relevant if
                        obtaining proteins]. Default == 'n'.
  -f {y,n,Y,N}, -force {y,n,Y,N}
                        Default == 'n', which means the program will not
                        overwrite existing files. Specify 'y' to allow this
                        behaviour at your own risk.
```
