#! python3
# Biopython based ORF Finder - SignalP version
# This script will obtain open reading frames from a fasta-formatted file containing nucleotide transcripts.
# Constraints can be altered to vary the strictness with which we accept or reject alternative start codons. 
# The output is a fasta file containing any number of ORFs that the user has specified based upon a minimum and maximum length also specified.
# Addition: This script will seek out ORFs with signal peptide start sites

# Load packages
import re, os, argparse, platform, subprocess, time, hashlib
from Bio import SeqIO
from pathlib import Path

# Define functions for later use [Rewriting this script would make it so much prettier/easy to manage...]
def validate_args(args):
        # Validate input file location
        if not os.path.isfile(args.fileName):
                print('The specified input file does not exist (i.e., "' + args.fileName + ')".')
                print('Make sure you entered this correctly and try again.')
                quit()
        # Validate parameters
        if args.minProLen < 1:
                print("minProLen must be >= 1. Fix your input and try again.")
                quit()
        if args.maxProLen < 0:
                print("maxProLen must be a positive number. Fix your input and try again.")
                quit()
        if args.maxProLen != 0 and args.minProLen > args.maxProLen:
                print("minProLen must be <= maxProLen (if maxProLen is not 0). Fix your input and try again.")
                quit()
        if args.unresolvedCodon < 0:
                print("unresolvedCodon must be >= 1. Fix your input and try again.")
                quit()
        if args.translationTable < 1 or args.translationTable > 31:
                print('-t translationTable value ranging from 1 to 31 inclusive. Fix this and try again.')
                quit()
        # Validate accessory program locations
        if args.cygwindir == None:
                args.cygwindir = ''
        if args.signalpdir == None:
                args.signalpdir = ''
        if args.signalpdir == None:
                print('signalpdir argument was not specified when -signalp was provided; fix your input and try again.')
                quit()
        if not os.path.isfile(os.path.join(args.signalpdir, 'signalp')):
                print('I am unable to locate the signalp execution file "signalp" within specified directory (' + args.signalpdir + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if platform.system() == 'Windows':
                program_execution_check(os.path.join(args.cygwindir, 'bash.exe --version'))
                cygwin_program_execution_check(os.path.dirname(args.outputFileName), args.cygwindir, args.signalpdir, 'signalp -h')
        else:
                program_execution_check(os.path.join(args.signalpdir, 'signalp -h'))
        # Validate output file location
        args.outputFileName = os.path.abspath(args.outputFileName)
        if not os.path.isdir(os.path.dirname(args.outputFileName)):
                print('The specified output directory does not exist (i.e., "' + os.path.dirname(args.outputFileName) + ')".')
                print('Create this directory first and then try again.')
                quit()
        else:
                # Check for file overwrites
                if args.sequenceType.lower() != 'both':
                        if os.path.isfile(args.outputFileName) and args.force == False:
                                print('There is already a file at "' + args.outputFileName + '".')
                                print('Either specify a new file name, delete this older file, or provide the -force argument')
                                quit()
                        elif os.path.isfile(args.outputFileName) and args.force == True:
                                os.remove(args.outputFileName)
                        args.protOutName = None # For consistency, since the below else statement generates these values
                        args.nuclOutName = None
                else:
                        # Derive output file names
                        outPrefix = args.outputFileName.rsplit('.', maxsplit=1)
                        if len(outPrefix) == 1:
                                outPrefix.append('.fasta')
                        args.protOutName = outPrefix[0] + '_prot.' + outPrefix[1]
                        args.nuclOutName = outPrefix[0] + '_nucl.' + outPrefix[1]

                        if os.path.isfile(args.protOutName) and args.force == False:
                                print('There is already a file at "' + args.protOutName + '". Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
                                quit()
                        elif os.path.isfile(args.protOutName) and args.force == True:
                                os.remove(args.protOutName)

                        if os.path.isfile(args.nuclOutName) and args.force == False:
                                print('There is already a file at "' + args.nuclOutName + '". Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
                                quit()
                        elif os.path.isfile(args.nuclOutName) and args.force == True:
                                os.remove(args.nuclOutName)
        return args

def program_execution_check(cmd):
        run_cmd = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
        cmdout, cmderr = run_cmd.communicate()
        if cmderr.decode("utf-8") != '' and not cmderr.decode("utf-8").startswith('Usage'):
                print('Failed to execute program "' + cmd + '". Is this executable in the location specified/discoverable in your PATH, or does the executable even exist? I won\'t be able to run properly if I can\'t execute this program.')
                print('---')
                print('stderr is below for debugging purposes.')
                print(cmderr.decode("utf-8"))
                print('Program closing now.')
                quit()

def cygwin_program_execution_check(outDir, cygwinDir, exeDir, exeFile):
        # Format script for cygwin execution
        scriptText = Path(exeDir, exeFile).as_posix()
        scriptFile = tmp_file_name_gen('tmpscript', '.sh', scriptText)
        with open(os.path.join(outDir, scriptFile), 'w') as fileOut:
                fileOut.write(scriptText)
        # Format cmd for execution
        cmd = os.path.join(cygwinDir, 'bash') + ' -l -c ' + os.path.join(outDir, scriptFile).replace('\\', '/')
        run_cmd = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
        cmdout, cmderr = run_cmd.communicate()
        os.remove(os.path.join(outDir, scriptFile))   # Clean up temporary file
        if cmderr.decode("utf-8") != '' and not 'perl: warning: falling back to the standard locale' in cmderr.decode("utf-8").lower():
                '''Need the above extra check for signalP since, on Windows at least, you can receive perl warnings which don't impact
                program operations. I think if that 'falling back' line is in stderr, nothing more serious will be present in stderr -
                this isn't completely tested, however.'''
                print('Failed to execute ' + exeFile + ' program via Cygwin using "' + cmd + '". Is this executable in the location specified/discoverable in your PATH, or does the executable even exist? I won\'t be able to run properly if I can\'t execute this program.')
                print('---')
                print('stderr is below for debugging purposes.')
                print(cmderr.decode("utf-8"))
                print('Program closing now.')
                quit()

def run_signalp_sequence(signalpdir, cygwindir, organism, tmpDir, seqID, protString):
        # Determine whether seqId and protString values are the proper type
        if not type(seqID) == str and not type(protString) == str:
                if not type(seqID) == list and not type(protString) == list:
                        print('run_signalp_sequence: seqID and protString inputs should both be str or both be list; this isn\'t true here, so I cannot procede.')
                        print('Fix the code leading up to this function call.')
                        quit()
                # If they are lists, ensure they have the same length
                else:
                        if len(seqID) != len(protString):
                                print('run_signalp_sequence: seqID and protString inputs are lists of nonequivalent length; I cannot procede unless this is true.')
                                print('Fix the code leading up to this function call.')
                                quit()
        # Generate temporary file for sequence
        if type(seqID) == list:
                tmpFileName = tmp_file_name_gen(os.path.join(tmpDir, 'tmp_sigpInput_' + ''.join([sid[0:5] for sid in seqID])[0:25] + '_'), '.fasta', ''.join([prot[0:10] for prot in protString]))
        else:
                tmpFileName = tmp_file_name_gen(os.path.join(tmpDir, 'tmp_sigpInput_' + seqID + '_'), '.fasta', protString)
        with open(tmpFileName, 'w') as fileOut:
                if type(seqID) == list:
                        for i in range(len(seqID)):
                                fileOut.write('>' + seqID[i].lstrip('>') + '\n' + protString[i] + '\n')      # lstrip any > characters just in case they're already present
                else:
                        fileOut.write('>' + seqID.lstrip('>') + '\n' + protString + '\n')
        # Run signalP
        if type(seqID) == list:
                sigpResultFile = tmp_file_name_gen(os.path.join(tmpDir, 'tmp_sigpResults_' + ''.join([sid[0:5] for sid in seqID])[0:25] + '_'), '.txt', ''.join([prot[0:10] for prot in protString]))
        else:
                sigpResultFile = tmp_file_name_gen(os.path.join(tmpDir, 'tmp_sigpResults_' + seqID + '_'), '.txt', protString)
        signalp_unthreaded(signalpdir, cygwindir, organism, tmpDir, tmpFileName, sigpResultFile)
        # Join and parse signalP results files
        sigPredictions = {}
        with open(sigpResultFile, 'r') as fileIn:
                for line in fileIn:
                        if line.startswith('#'):
                                continue
                        sl = line.split('\t')
                        sigPredictions[sl[0]] = [int(sl[3]), int(sl[4]), float(sl[5])] # [start, stop, score]
        # Clean up temporary 
        os.remove(tmpFileName)
        os.remove(sigpResultFile)
        # Return signalP prediction dictionary
        return sigPredictions

def signalp_unthreaded(signalpdir, cygwindir, organism, tmpDir, fastaFile, sigpResultFile):
        # Get the full fasta file location
        fastaFile = os.path.abspath(fastaFile)
        # Format signalP script text
        scriptText = '"' + os.path.join(signalpdir, 'signalp') + '" -t ' + organism + ' -f short -n "' + sigpResultFile + '" "' + fastaFile + '"'
        # Generate a script for use with cygwin (if on Windows)
        if platform.system() == 'Windows':
                sigpScriptFile = os.path.join(tmpDir, tmp_file_name_gen('tmp_sigpScript_' + os.path.basename(fastaFile), '.sh', scriptText))
                with open(sigpScriptFile, 'w') as fileOut:
                        fileOut.write(scriptText.replace('\\', '/'))
        # Run signalP depending on operating system
        if platform.system() == 'Windows':
                cmd = os.path.join(cygwindir, 'bash') + ' -l -c "' + sigpScriptFile.replace('\\', '/') + '"'
                runsigP = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
                sigpout, sigperr = runsigP.communicate()
                os.remove(sigpScriptFile)       # Clean up temporary file
        else:
                os.putenv("PYTHONPATH",os.pathsep.join([os.getenv("PYTHONPATH",""),signalpdir]))
                runsigP = subprocess.Popen(scriptText, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
                sigpout, sigperr = runsigP.communicate()
        # Process output
        okayLines = ['is an unknown amino amino acid', 'perl: warning:', 'LC_ALL =', 'LANG =', 'are supported and installed on your system']
        for line in sigperr.decode("utf-8").split('\n'):
                # If sigperr indicates null result, create an output file we can skip later
                if line.rstrip('\n') == '# No sequences predicted with a signal peptide':
                        with open(sigpResultFile, 'w') as fileOut:
                                fileOut.write(line)
                        break
                # Check if this line has something present within okayLines
                okay = 'n'
                for entry in okayLines:
                        if entry in line or line == '':
                                okay = 'y'
                                break
                if okay == 'y':
                        continue
                # If nothing matches the okayLines list, we have a potentially true error
                else:
                        raise Exception('SignalP error occurred when processing file name ' + fastaFile + '. Error text below\n' + sigperr.decode("utf-8"))

def tmp_file_name_gen(prefix, suffix, hashString):
        # Main function
        tmpHash = hashlib.md5(bytes(str(hashString) + str(time.time()), 'utf-8') ).hexdigest()       # This should always give us something unique even if the string for hashString is the same across different runs
        while True:
                if os.path.isfile(prefix + tmpHash + suffix):
                        tmpHash += 'X'
                else:
                        return prefix + tmpHash + suffix

def output_func(outputProt, outputNucl, ongoingCount, outputFileName, sequenceType, status, protOutName=None, nuclOutName=None):
        if sequenceType.lower() == 'prot':
                if (ongoingCount%10000 == 0 or status == 'final') and os.path.isfile(os.path.join(os.getcwd(), outputFileName)) == False:
                        with open(outputFileName, 'w') as output:
                                output.write('\n'.join(outputProt))
                        outputProt = []
                elif (ongoingCount%10000 == 0 or status == 'final') and os.path.isfile(os.path.join(os.getcwd(), outputFileName)) == True:
                        with open(outputFileName, 'a') as output:
                                output.write('\n')
                                output.write('\n'.join(outputProt))
                        outputProt = []
        elif sequenceType.lower() == 'nucl':
                if (ongoingCount%10000 == 0 or status == 'final') and os.path.isfile(os.path.join(os.getcwd(), outputFileName)) == False:
                        with open(outputFileName, 'w') as output:
                                output.write('\n'.join(outputNucl))
                        outputNucl = []
                elif (ongoingCount%10000 == 0 or status == 'final') and os.path.isfile(os.path.join(os.getcwd(), outputFileName)) == True:
                        with open(outputFileName, 'a') as output:
                                output.write('\n')
                                output.write('\n'.join(outputNucl))
                        outputNucl = []
        else:
                if (ongoingCount%10000 == 0 or status == 'final') and os.path.isfile(os.path.join(os.getcwd(), protOutName)) == False:       # Doesn't matter if we check for protOutName or nuclOutName. Theoretically, if one of the files is deleted while the program is running it could be problematic, but I mean, what can I really do about that without child-proofing the script excessively?
                        with open(protOutName, 'w') as protFile, open(nuclOutName, 'w') as nuclFile:
                                protFile.write('\n'.join(outputProt))
                                nuclFile.write('\n'.join(outputNucl))
                        outputProt = []
                        outputNucl = []
                elif (ongoingCount%10000 == 0 or status == 'final') and os.path.isfile(os.path.join(os.getcwd(), protOutName)) == True:
                        with open(protOutName, 'a') as protFile, open(nuclOutName, 'a') as nuclFile:
                                protFile.write('\n')
                                nuclFile.write('\n')
                                protFile.write('\n'.join(outputProt))
                                nuclFile.write('\n'.join(outputNucl))
                        outputProt = []
                        outputNucl = []
        return outputProt, outputNucl

### USER INPUT
usage = """%(prog)s reads in a fasta formatted file containing nucleotide sequences and, following user-specified parameters,
produces an output fasta file containing potential open reading frames (ORFs) as nucleotides/protein translations/both.
"""
# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-input", dest="fileName",
                   help="Input fasta file name")
p.add_argument("-o", "-output", dest="outputFileName",
                   help="Output fasta file name")
# Opts
p.add_argument("-min", "--minimum", type=int, dest="minProLen",
                   help="Minimum ORF amino acid length. Default == 30.", default=30)
p.add_argument("-max", "--maximum", type=int, dest="maxProLen",
                   help="Optional specification for maximum ORF amino acid length. Default == 0, which means there is no upper limit.", default=0)
p.add_argument("-st", "--seqtype", dest="sequenceType", choices = ['prot', 'nucl', 'both', 'PROT', 'NUCL', 'BOTH'],
                   help="Specify the type of output you want to generate (i.e., protein translated ORF, nucleotide CDS, or both). If you specify 'both', two outputs with '_prot' and '_nucl' suffix will be generated. Default == 'prot'.", default="prot")
p.add_argument("-f", "--force", dest="force", action="store_true", default=False,
                   help="Allow files to be overwritten at your own risk.")
p.add_argument("-u", "--unresolved", dest="unresolvedCodon", type=int,
                   help="Default == 0, which means the program will not discover ORFs with unresolved codons. If you want to risk chimeric ORF formation, you can change this value. You MUST validate any ORFs with unresolved portions. Recommended for this value to be less than 5.", default=0)
p.add_argument("-n", "--no_orf_num", dest="noOrfNum", action='store_true',
                   help="Provide this argument to prevent ORF nums being appended to sequence IDs. This can be useful when obtaining 1 ORF per transcript.", default=False)
p.add_argument("-t", "--translation", dest="translationTable", type=int, default=1,
                   help="Optionally specify the NCBI numeric genetic code to utilise for CDS translation (if relevant); this should be an integer from 1 to 31 (default == 1 i.e., Standard Code)")
# SignalP opts
p.add_argument("-sigpdir", "-signalpdir", dest="signalpdir", type=str,
               help="""Specify the directory where signalp executables are located.""")
p.add_argument("-sigporg", dest="signalporg", type = str, choices = ['euk', 'gram-', 'gram+'], default='euk',
               help="""Specify the type of organism for SignalP from the available
               options. Refer to the SignalP manual if unsure what these mean (default == 'euk').""")
p.add_argument("-c", "-cygwindir", dest="cygwindir", type=str, default="",
               help="""Cygwin is required since you are running this program on a Windows computer.
               Specify the location of the bin directory here or, if this is already in your PATH, you can leave this blank."""
               if platform.system() == 'Windows' else argparse.SUPPRESS)

## HARD-CODED TESTING
#args = p.parse_args()
#args = validate_args(args)
fileName = r"F:\toxins_annot\analysis\sigp_orf_testing\Telmatactis.fa"
translationTable = 1
unresolvedCodon = 0
minProLen = 30
maxProLen = 0
seqType = "both"
force = False
noORfNum = False
signalpdir = r"D:\Bioinformatics\Protein_analysis\signalp-4.1f.CYGWIN\signalp-4.1"
signalporg = "euk"
cygwindir = r""
## END

xRegex = re.compile(r'X+')                                              # Regex used to find start and stop positions of unresolved regions that are shorter than the cut-off
                
### RATIONALE FOR UNRESOLVED REGIONS ###
# I had to decide how to handle unresolved regions. I believe there are two valid approaches to this. The first is to replace any unresolved regions with stop codons and let the rest of the script process it like normal.
# This appears to be what NCBI does for their ORF Finder. The benefits of this approach is that you can't accidentally form chimeras. However, I think there is a second approach that has merit. If you are working with a genome that has
# short and rare occurrences of unresolved regions, you might not want to split up large ORFs on the basis of them having a very short stretch of unresolved codons. In this case, we can choose to not replace unresolved regions with stop codons
# if the region is shorter than an arbitrary and small limit. This isn't exactly perfect since even very short unresolved regions might hide stop codons. Thus, this option should be OFF by default, but we can provide the users an
# on/off switch so they can decide whether they want to risk discovering chimeric ORFs - providing a text prompt if the option is switched ON to verify any ORFs with unresolved regions would wash my hands of any mistake.

#if args.unresolvedCodon != 0:
if unresolvedCodon != 0:
        print('Program has noted that you are allowing the discovery of ORFs with unresolved codon regions. This is risky behaviour, since this program cannot guarantee that an unresolved region does not contain a stop codon. Subsequently, you can have chimeras form from two separate ORFs. YOU MUST VERIFY ANY ORFS WITH UNRESOLVED REGIONS! The best way to do this is with BLAST against homologous proteins. You have been warned.')

# Load the fasta file as a generator object, get the total number of sequences in the file, then re-load it for the upcoming loop
#records = SeqIO.parse(open(args.fileName, 'r'), 'fasta')
#records = SeqIO.parse(open(fileName, 'r'), 'fasta')
#totalCount = 0
#for record in records:
#        totalCount += 1
#records = SeqIO.parse(open(args.fileName, 'r'), 'fasta')
records = SeqIO.parse(open(fileName, 'r'), 'fasta')

### CORE PROCESSING LOOP
print('Starting the core processing of this script now. Progress bar is displayed below. Every 10,000 sequences, current progress will be saved to the output file(s) to reduce memory footprint.')

# Declare overall values needed before loop start
ongoingCount = 0
outputProt = []                                 # These get reset whenever we output to file
outputNucl = []

# Get the nucleotide (record) out of our generator (records) and grab them ORFs!
for record in records:
        tempOverallProt = []
        tempOverallNucl = []
        # Parental loop
        for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
                for frame in range(3):
                        nuc = nuc.upper() # Just in case?
                        length = 3 * ((len(record)-frame) // 3)
                        frameNuc = str(nuc[frame:frame+length])
                        #frameProt = str(nuc[frame:frame+length].translate(table=args.translationTable))
                        frameProt = str(nuc[frame:frame+length].translate(table=translationTable))
                        # Split protein/nucleotide into corresponding ORFs
                        ongoingLength = 0                                       # The ongoingLength will track where we are along the unresolvedProt sequence for getting the nucleotide sequence
                        splitNucleotide = []
                        splitProtein = []
                        frameProt = frameProt.split('*')
                        for i in range(len(frameProt)):
                                if len(frameProt) == 1 or i + 1 == len(frameProt):    # This means the splitProtein has no stop codons or we're looking at the last ORF which also has no stop codon
                                        splitProtein.append(frameProt[i])
                                        splitNucleotide.append(frameNuc[ongoingLength:ongoingLength+len(frameProt[i])*3])      # This will grab the corresponding nucleotide region
                                        ongoingLength += len(frameProt[i])*3
                                else:
                                        splitProtein.append(frameProt[i] + '*')
                                        splitNucleotide.append(frameNuc[ongoingLength:ongoingLength+len(frameProt[i] + '*')*3])       
                                        ongoingLength += (len(frameProt[i]) + 1)*3
                        
                        # Fix unresolved regions
                        resolvedProt = []
                        resolvedNuc = []
                        indicesForDel = []                                                              # We'll hold onto the indices of any splitProtein components that have X's in them. We'll loop through this later to delete them from splitProtein/Nucleotide, then we'll add the resolved segments to splitProtein/Nucleotide
                        for i in range(len(splitProtein)):
                                if 'X' in splitProtein[i]:
                                        posProt = []
                                        for x in re.finditer(xRegex, splitProtein[i]):
                                                #if x.end() - x.start() > args.unresolvedCodon:
                                                if x.end() - x.start() > unresolvedCodon:
                                                        posProt += [x.start(), x.end()]
                                        if posProt == []:                                               # If posProt still == [], that means we didn't find any unresolved regions that exceed our cut-off
                                                continue
                                        indicesForDel.insert(0, i)                                      # Insert it at 0 so we don't need to sort it at the end [we need to loop through a reversed list so we can directly delete the indices without messing up the order of splitProtein/Nucleotide]
                                        # Pull out resolved regions
                                        resolvedProt.append(splitProtein[i][:posProt[0]])               # We loop through our posProt by first grabbing everything before our first unresolved region
                                        resolvedNuc.append(splitNucleotide[i][:posProt[0]*3])
                                        for x in range(1, len(posProt)-1, 2):                           # We now, by skipping the first and last entry in posProt, can compare every coordinate pair that corresponds to a resolved region
                                                start = posProt[x]
                                                end = posProt[x+1]
                                                resolvedProt.append(splitProtein[i][start:end])
                                                resolvedNuc.append(splitNucleotide[i][start*3:end*3])
                                        resolvedProt.append(splitProtein[i][posProt[-1]:])              # We can now grab everything after our last unresolved region. If there was only one unresolved region, we never enter the above loop and just use the coordinate pair to get our start and end sequences
                                        resolvedNuc.append(splitNucleotide[i][posProt[-1]*3:])
                        # Delete old entries and add resolved entries
                        for index in indicesForDel:
                                del splitProtein[index]
                                del splitNucleotide[index]
                        splitProtein += resolvedProt                                                    # If we don't find any unresolved regions we wanted to delete, resolvedProt will be empty so nothing happens
                        splitNucleotide += resolvedNuc

                        # Enter the main processing loop with our resolved regions
                        prots = []
                        nucs = []
                        for i in range(len(splitProtein)):                                                      # Note that I have done a 'for i in range...' loop rather than a 'for value in splitProtein' loop which would have been simpler for a reason explained below on the 'elif i + 1 ==' line
                                # Declare blank values needed for each potential ORF region so we can tell which things were 'found'
                                noneCodonContingency = None
                                # Process sequences to determine whether we're ignoring this, or adding an asterisk for length counts
                                #if len(splitProtein[i]) < args.minProLen:                    # Disregard sequences that won't meet the size requirement without further processing
                                if len(splitProtein[i]) < minProLen:                    # Disregard sequences that won't meet the size requirement without further processing
                                        continue
                                #elif args.maxProLen != 0 and len(splitProtein[i]) > args.maxProLen:
                                elif maxProLen != 0 and len(splitProtein[i]) > maxProLen:
                                        continue
                                acceptedPro = str(splitProtein[i])
                                # Obtain possible start sites
                                ## Canonical start coding
                                ALLOWED_SHORT_RATIO = 0.50 # Arbitrary; we want to prevent a sequence from being shortened excessively
                                ABSOLUTE_SHORTEN_LEN = 50 # Arbitrary; we want to prevent a sequence being shortened more than is "realistic" in search of a start codon
                                allowedShortenLen = int(round(len(acceptedPro)*ALLOWED_SHORT_RATIO, 0))
                                
                                startIndices = []
                                for x in range(0, len(acceptedPro)):
                                        if acceptedPro[x].lower() == "m" and x <= allowedShortenLen:
                                                startIndices.append([x, 0]) # 0 is the "best" and represents a canonical start
                                ## Alternative start coding
                                nucSeqOfProt = splitNucleotide[i]                       # Don't need to do it, but old version of script extensively uses this value and cbf changing it
                                codons = re.findall('..?.?', nucSeqOfProt)              # Pulls out a list of codons from the nucleotide
                                for codon in codons:                                    # Cycle through this list of codons to find the first alternative start of the normal class (GTG and TTG) and the rare class (CTG)
                                        if codon == 'GTG' or codon == 'TTG':
                                                index = codons.index(codon)
                                                if index <= allowedShortenLen:
                                                        startIndices.append([index, 1]) # This will save the position of the first GTG or TTG encountered. Note that by breaking after this,  we stop looking for CTG as it is irrelevant after this
                                                break
                                        elif codon == 'CTG':
                                                if noneCodonContingency == None:        # noneCodonContingency is set to None at the end of each loop. Thus, this line of code will 'capture' the position of the first CTG in a sequence if a GTG or TTG was not encountered first
                                                        index = codons.index(codon)
                                                        if index <= allowedShortenLen:
                                                                startIndices.append([index, 2])
                                startIndices.sort(key = lambda x: x[0])
                                # Skip if no hits are obtained
                                if startIndices == []:
                                        continue
                                # Skip the hit if it doesn't meet our minimum length requirement anymore
                                #if len(acceptedPro[startIndices[0][0]:]) < args.minProLen: # Culling is necessary since we will have shortened the sequence somewhat unless we accepted the topHit as being a no codon start ORF. Note that we will here consider a stop codon in the length of the protein, such that a protein with 99AAs and a stop will pass a minimum 100AA length test. I think this is fair since not all regions here have a stop codon which allows weight to be added to these cases, especially since a stop codon is still conserved as part of an ORF.
                                if len(acceptedPro[startIndices[0][0]:]) < minProLen: # Culling is necessary since we will have shortened the sequence somewhat unless we accepted the topHit as being a no codon start ORF. Note that we will here consider a stop codon in the length of the protein, such that a protein with 99AAs and a stop will pass a minimum 100AA length test. I think this is fair since not all regions here have a stop codon which allows weight to be added to these cases, especially since a stop codon is still conserved as part of an ORF.
                                        continue
                                # Cull individual start sites that don't meet our minimum length requirement
                                exitLoop = False
                                while True:
                                        if exitLoop == True:
                                                break
                                        exitLoop = True
                                        for x in range(len(startIndices)):
                                                if startIndices[x][0] > allowedShortenLen:
                                                        del startIndices[x]
                                                        exitLoop = False
                                                        break
                                # Skip if we've culled all our hits
                                if startIndices == []:
                                        continue
                                # Assess remaining start sites for signal peptide presence
                                seqIDs = []
                                indexProts = []
                                for x in range(len(startIndices)):
                                        seqIDs.append(str(x))
                                        indexProts.append(acceptedPro[startIndices[x][0]:])
                                # Run signalP prediction and associate relevant results
                                #sigpPredictions = run_signalp_sequence(str(args.signalpdir), args.cygwindir, args.signalporg, str(args.signalpdir), seqIDs, indexProts)
                                sigpPredictions = run_signalp_sequence(str(signalpdir), cygwindir, signalporg, str(signalpdir), seqIDs, indexProts)
                                if sigpPredictions == {}:
                                        continue
                                sigpIndices = []
                                for x in range(len(startIndices)):
                                        if str(x) in sigpPredictions:
                                                startIndices[x].append(sigpPredictions[str(x)][2]) # SignalP score; closest to 1 is best
                                                sigpIndices.append(startIndices[x])
                                # If no signal peptides were found, skip sequence
                                if sigpIndices == []:
                                        continue
                                # Sort candidates on the basis of evidence signalP > canonical > startSite
                                sigpIndices.sort(key = lambda x: (x[1], x[0])) # i.e., [canonical(0 is best), startSite(smaller is better)]
                                # Identify the best candidate with a signalP prediction
                                bestCandidate = sigpIndices[0]
                                EXTENSION_FOR_BETTER_SIGP_SCORE = 7 # Arbitrary; this should help to address situations where two M's are close together, and one has a better signalP prediction score
                                STEPPING_STONE_EXTENSION_LIMITER = 3 # Arbitrary; after extending for a better sigp score once, we enforce a stricter mechanism for extensions afterwards
                                for sigpIndex in sigpIndices: ## TBD: Include extension
                                        # If new candidate has a better start site than bestCandidate...
                                        if sigpIndex[1] < bestCandidate[1]: # this measures canonical (0 is best)
                                                bestCandidate = sigpIndex
                                        # If new candidate has an equivalent start site...
                                        elif sigpIndex[1] == bestCandidate[1]:
                                                # ...and it is much longer than bestCandidate...
                                                if sigpIndex[0] + EXTENSION_FOR_BETTER_SIGP_SCORE < bestCandidate[0]:
                                                        bestCandidate = sigpIndex
                                                # ...and it isn't much shorter than bestCandidate...
                                                elif sigpIndex[0] - EXTENSION_FOR_BETTER_SIGP_SCORE <= bestCandidate[0]:
                                                        # ...and its signalP score is better than bestCandidate...
                                                        if sigpIndex[2] > bestCandidate[2]: # this measures signalpep (1 is best)
                                                                bestCandidate = sigpIndex
                                                                EXTENSION_FOR_BETTER_SIGP_SCORE = STEPPING_STONE_EXTENSION_LIMITER # Enforce stricter limit now that we've extended once
                                # Extract the relevant protein and nucleotide to add to our output list
                                protein = acceptedPro[bestCandidate[0]:]
                                nucleotide = nucSeqOfProt[bestCandidate[0]*3:]
                                prots.append(protein)
                                nucs.append(nucleotide)


        # Format and produce the output of this script
        tempOutputProt = []
        tempOutputNucl = []
        if args.sequenceType.lower() == 'prot' or args.sequenceType.lower() == 'both':
                for i in range(0, args.hitsToPull):                          
                        if tempOverallProt[i] == '-':                   # Because we made sure all the tempM/Alt/NoneLists had '-' added to pad out the list to have a length equal to the value of hitsToPull, when we cycle through our tempOverallList, we will often encounter '-' characters which signify the end of relevant ORFs identified in this sequence  
                                break                                   # Break out of this loop once we've fasta formatted all relevant ORF hits
                        if args.noOrfNum == False:
                                tempOutputProt.append('>' + record.id + '_ORF' + str(i+1) + '\n' + tempOverallProt[i])
                        else:
                                tempOutputProt.append('>' + record.id + '\n' + tempOverallProt[i])
                if len(tempOutputProt) > 0:    # Need this check for sequences that don't get any hits
                        outputProt.append('\n'.join(tempOutputProt))

        if args.sequenceType.lower() == 'nucl' or args.sequenceType.lower() == 'both':
                for i in range(0, args.hitsToPull):                          
                        if tempOverallNucl[i] == '-':                   # Because we made sure all the tempM/Alt/NoneLists had '-' added to pad out the list to have a length equal to the value of hitsToPull, when we cycle through our tempOverallList, we will often encounter '-' characters which signify the end of relevant ORFs identified in this sequence  
                                break                                   # Break out of this loop once we've fasta formatted all relevant ORF hits
                        if args.noOrfNum == False:
                                tempOutputNucl.append('>' + record.id + '_ORF' + str(i+1) + '\n' + tempOverallNucl[i])
                        else:
                                tempOutputNucl.append('>' + record.id + '\n' + tempOverallNucl[i])
                if len(tempOutputNucl) > 0:
                        outputNucl.append('\n'.join(tempOutputNucl))

        ongoingCount += 1
        
        # Save backup if ongoingCount == 10,000. There are two sections here, the first will create the file on the first loop, the second will add to the file on subsequent loops
        outputProt, outputNucl = output_func(outputProt, outputNucl, ongoingCount, args.outputFileName, args.sequenceType, 'processing', args.protOutName, args.nuclOutName)

# Dump the last few results after the script has finished, or create the output if there were less than 10,000 sequences
outputProt, outputNucl = output_func(outputProt, outputNucl, ongoingCount, args.outputFileName, args.sequenceType, 'final', args.protOutName, args.nuclOutName)
records.close()

#### SCRIPT ALL DONE
print('Program completed successfully!')
