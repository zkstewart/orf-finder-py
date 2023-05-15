#! python3
# Biopython based ORF Finder - SignalP version
# This script will obtain open reading frames from a fasta-formatted file containing nucleotide transcripts.
# Constraints can be altered to vary the strictness with which we accept or reject alternative start codons. 
# The output is a fasta file containing any number of ORFs that the user has specified based upon a minimum and maximum length also specified.
# Addition: This script will seek out ORFs with signal peptide start sites

# Load packages
import re, os, argparse, platform, subprocess, time, hashlib, queue, random, shutil
from Bio import SeqIO
from pathlib import Path
from threading import Thread

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
    if args.threads < 0:
        print("threads must be >= 1. Fix your input and try again.")
        quit()
    if args.translationTable < 1 or args.translationTable > 31:
        print('-t translationTable value ranging from 1 to 31 inclusive. Fix this and try again.')
        quit()
    # Validate accessory program locations
    if args.cygwindir == None:
        args.cygwindir = ''
    if args.signalpdir == None:
        print('signalpdir argument was not specified; fix your input and try again.')
        quit()
    if not os.path.isfile(os.path.join(args.signalpdir, 'signalp')):
        print('I am unable to locate the signalp execution file "signalp" within specified directory (' + args.signalpdir + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if platform.system() == 'Windows':
        program_execution_check(os.path.join(args.cygwindir, 'bash.exe --version'))
        cygwin_program_execution_check(os.path.dirname(os.path.abspath(args.outputFileName)), args.cygwindir, args.signalpdir, 'signalp -h')
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
        while True:
            tmpFileName = tmp_file_name_gen(os.path.join(tmpDir, 'tmp_sigpInput_' + ''.join([sid[0:5] for sid in seqID])[0:25] + '_'), '.fasta', ''.join([prot[0:10] for prot in protString]) + str(time.time()))
            if os.path.isfile(tmpFileName):
                continue
            else:
                break

    else:
        while True:
            tmpFileName = tmp_file_name_gen(os.path.join(tmpDir, 'tmp_sigpInput_' + seqID + '_'), '.fasta', protString + str(time.time()))
            if os.path.isfile(tmpFileName):
                continue
            else:
                break
    with open(tmpFileName, 'w') as fileOut:
        if type(seqID) == list:
            for i in range(len(seqID)):
                fileOut.write('>' + seqID[i].lstrip('>') + '\n' + protString[i] + '\n')      # lstrip any > characters just in case they're already present
        else:
            fileOut.write('>' + seqID.lstrip('>') + '\n' + protString + '\n')
    # Run signalP
    if type(seqID) == list:
        sigpResultFile = tmp_file_name_gen(os.path.join(tmpDir, 'tmp_sigpResults_' + ''.join([sid[0:5] for sid in seqID])[0:25] + '_'), '.txt', ''.join([prot[0:10] for prot in protString]) + str(time.time()))
    else:
        sigpResultFile = tmp_file_name_gen(os.path.join(tmpDir, 'tmp_sigpResults_' + seqID + '_'), '.txt', protString + str(time.time()))
    signalp_unthreaded(signalpdir, cygwindir, organism, tmpDir, tmpFileName, sigpResultFile)
    # Join and parse signalP results files
    sigPredictions = {}
    with open(sigpResultFile, 'r') as fileIn:
        for line in fileIn:
            if line.startswith('#'):
                continue
            sl = line.split('\t')
            sigPredictions[sl[0]] = [int(sl[3]), int(sl[4]), float(sl[5])] # [start, stop, score]
    # Clean up temporary files
    os.remove(tmpFileName)
    os.remove(sigpResultFile)
    # Return signalP prediction dictionary
    return sigPredictions

def signalp_unthreaded(signalpdir, cygwindir, organism, tmpDir, fastaFile, sigpResultFile):
    # Get the full fasta file location
    fastaFile = os.path.abspath(fastaFile)
    # Format signalP script text
    sigpTmpDir = tmp_file_name_gen(os.path.join(signalpdir, 'tmp_sigp_run_'), '', str(time.time()) + sigpResultFile)
    scriptText = '"{0}" -t {1} -f short -n "{2}" -T "{3}" "{4}"'.format(os.path.join(signalpdir, 'signalp'), organism, sigpResultFile, sigpTmpDir, fastaFile)
    # Generate a script for use with cygwin (if on Windows)
    if platform.system() == 'Windows':
        sigpScriptFile = os.path.join(tmpDir, tmp_file_name_gen('tmp_sigpScript_', '.sh', scriptText + str(time.time())))
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
    okayLines = ['is an unknown amino amino acid', 'perl: warning:', 'LC_ALL =', 'LANG =', 'are supported and installed on your system', '# temporary directory will not be removed']
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
    # Clean up tmp dir
    if os.path.isdir(sigpTmpDir):
        shutil.rmtree(sigpTmpDir)

def tmp_file_name_gen(prefix, suffix, hashString):
    # Main function
    tmpHash = hashlib.sha256(bytes(str(hashString) + str(time.time()) + str(random.randint(0, 100000)), 'utf-8') ).hexdigest()       # This should always give us something unique even if the string for hashString is the same across different runs
    while True:
        if os.path.isfile(prefix + tmpHash + suffix):
            tmpHash += 'X'
        else:
            return prefix + tmpHash + suffix

def orf_find_from_record(record, translationTable, unresolvedCodonLen, minProLen, maxProLen, signalpdir, signalporg, cygwindir):
    xRegex = re.compile(r'X+')                          # Regex used to find start and stop positions of unresolved regions that are shorter than the cut-off
    prots = []
    nucs = []
    for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
        for frame in range(3):
            nuc = nuc.upper() # Just in case?
            length = 3 * ((len(record)-frame) // 3)
            frameNuc = str(nuc[frame:frame+length])
            frameProt = str(nuc[frame:frame+length].translate(table=translationTable))
            # Split protein/nucleotide into corresponding ORFs
            ongoingLength = 0                       # The ongoingLength will track where we are along the unresolvedProt sequence for getting the nucleotide sequence
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
            indicesForDel = []                                  # We'll hold onto the indices of any splitProtein components that have X's in them. We'll loop through this later to delete them from splitProtein/Nucleotide, then we'll add the resolved segments to splitProtein/Nucleotide
            for i in range(len(splitProtein)):
                if 'X' in splitProtein[i]:
                    posProt = []
                    for x in re.finditer(xRegex, splitProtein[i]):
                        if x.end() - x.start() > unresolvedCodonLen:
                            posProt += [x.start(), x.end()]
                    if posProt == []:                           # If posProt still == [], that means we didn't find any unresolved regions that exceed our cut-off
                        continue
                    indicesForDel.insert(0, i)                      # Insert it at 0 so we don't need to sort it at the end [we need to loop through a reversed list so we can directly delete the indices without messing up the order of splitProtein/Nucleotide]
                    # Pull out resolved regions
                    resolvedProt.append(splitProtein[i][:posProt[0]])           # We loop through our posProt by first grabbing everything before our first unresolved region
                    resolvedNuc.append(splitNucleotide[i][:posProt[0]*3])
                    for x in range(1, len(posProt)-1, 2):               # We now, by skipping the first and last entry in posProt, can compare every coordinate pair that corresponds to a resolved region
                        start = posProt[x]
                        end = posProt[x+1]
                        resolvedProt.append(splitProtein[i][start:end])
                        resolvedNuc.append(splitNucleotide[i][start*3:end*3])
                    resolvedProt.append(splitProtein[i][posProt[-1]:])          # We can now grab everything after our last unresolved region. If there was only one unresolved region, we never enter the above loop and just use the coordinate pair to get our start and end sequences
                    resolvedNuc.append(splitNucleotide[i][posProt[-1]*3:])
            # Delete old entries and add resolved entries
            for index in indicesForDel:
                del splitProtein[index]
                del splitNucleotide[index]
            splitProtein += resolvedProt                            # If we don't find any unresolved regions we wanted to delete, resolvedProt will be empty so nothing happens
            splitNucleotide += resolvedNuc

            # Enter the main processing loop with our resolved regions
            for i in range(len(splitProtein)):                              # Note that I have done a 'for i in range...' loop rather than a 'for value in splitProtein' loop which would have been simpler for a reason explained below on the 'elif i + 1 ==' line
                # Declare blank values needed for each potential ORF region so we can tell which things were 'found'
                noneCodonContingency = None
                # Process sequences to determine whether we're ignoring this, or adding an asterisk for length counts
                if len(splitProtein[i]) < minProLen:            # Disregard sequences that won't meet the size requirement without further processing
                    continue
                elif maxProLen != 0 and len(splitProtein[i]) > maxProLen:
                    continue
                acceptedPro = str(splitProtein[i])
                # Obtain possible start sites
                ## Canonical start coding
                ALLOWED_SHORT_RATIO = 0.50 # Arbitrary; we want to prevent a sequence from being shortened excessively
                ABSOLUTE_SHORTEN_LEN = 50 # Arbitrary; we want to prevent a sequence being shortened more than is "realistic" in search of a start codon
                allowedShortenLen = min([int(round(len(acceptedPro)*ALLOWED_SHORT_RATIO, 0)), ABSOLUTE_SHORTEN_LEN])
                
                startIndices = []
                for x in range(0, len(acceptedPro)):
                    if acceptedPro[x].lower() == "m" and x <= allowedShortenLen:
                        startIndices.append([x, 0]) # 0 is the "best" and represents a canonical start
                ## Alternative start coding
                nucSeqOfProt = splitNucleotide[i]               # Don't need to do it, but old version of script extensively uses this value and cbf changing it
                codons = re.findall('..?.?', nucSeqOfProt)          # Pulls out a list of codons from the nucleotide
                for x in range(len(codons)):                # Cycle through this list of codons to find the first alternative start of the normal class (GTG and TTG) and the rare class (CTG)
                    codon = codons[x]
                    if codon == 'GTG' or codon == 'TTG':
                        if x <= allowedShortenLen:
                            startIndices.append([x, 1])     # This will save the position of the first GTG or TTG encountered.
                    elif codon == 'CTG':
                        if noneCodonContingency == None:    # noneCodonContingency is set to None at the end of each loop. Thus, this line of code will 'capture' the position of the first CTG in a sequence if a GTG or TTG was not encountered first
                            if x <= allowedShortenLen:
                                startIndices.append([x, 2])
                                noneCodonContingency = True # Stop accepting CTGs after the first to prevent "contamination"
                startIndices.sort(key = lambda x: x[0])
                # Skip if no hits are obtained
                if startIndices == []:
                    continue
                # Skip the hit if it doesn't meet our minimum length requirement anymore
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
                tmpDir = os.path.join(signalpdir, tmp_file_name_gen("tmp_orf_find_", "", str(time.time())))
                os.mkdir(tmpDir)
                sigpPredictions = run_signalp_sequence(str(signalpdir), cygwindir, signalporg, tmpDir, seqIDs, indexProts)
                os.rmdir(tmpDir)
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
    return prots, nucs

def record_worker(recordQ, outputQ, translationTable, unresolvedCodon, minProLen, maxProLen, signalpdir, signalporg, cygwindir):
    NoneType = type(None)
    while True:
        # Continue condition
        if recordQ.empty():
            time.sleep(0.5)
            continue
        # Grabbing condition
        record = recordQ.get()
        # Exit condition
        if type(record) == NoneType: # SeqIO doesn't allow comparison, need to do it another way
            recordQ.task_done()
            break
        # Perform work
        prots, nucs = orf_find_from_record(record, translationTable, unresolvedCodon, minProLen, maxProLen, signalpdir, signalporg, cygwindir)
        outputQ.put([record.description, prots, nucs])
        # Mark work completion
        recordQ.task_done()

def output_worker(outputQ, totalCount, outputFileName, sequenceType, protOutName, nuclOutName):
    ongoingCount = 0
    rememberPrint = -1
    while True:
        # Continue condition
        if outputQ.empty():
            time.sleep(0.5)
            continue
        # Grabbing condition
        output = outputQ.get()
        # Exit condition
        if output == None:
            outputQ.task_done()
            break
        # Perform work
        output_func(output[0], output[1], output[2], outputFileName, sequenceType, protOutName, nuclOutName)
        # Update progress bar
        ongoingCount += 1
        progress = ((ongoingCount)/totalCount)*100
        if int(progress)%1==0 and rememberPrint != int(progress):
            print('|' + str(int(progress)) + '% progress|' + str(ongoingCount+1) + ' sequences scanned for ORFs', end = '\r')       # Need to +1 to ongoingCount to counteract 0-index
            rememberPrint = int(progress)
        # Mark work completion
        outputQ.task_done()

def output_func(seqID, outputProts, outputNucs, outputFileName, sequenceType, protOutName=None, nuclOutName=None):
    # Create file(s) if necessary
    if (sequenceType.lower() == "prot" or sequenceType.lower() == "nucl") and os.path.isfile(outputFileName) == False:
        with open(outputFileName, 'w') as output:
            pass
    elif os.path.isfile(protOutName) == False:
        with open(protOutName, 'w') as output:
            pass
    elif os.path.isfile(nuclOutName) == False:
        with open(nuclOutName, 'w') as output:
            pass
    # Derive sequence IDs
    seqIDs = []
    for i in range(len(outputProts)):
        seqIDs.append(">{0}_ORF{1}".format(seqID, str(i+1)))
    # Write to relevant file(s)
    if sequenceType.lower() == "prot":
        with open(outputFileName, "a") as output:
            for i in range(len(seqIDs)):
                output.write("{0}\n{1}\n".format(seqIDs[i], outputProts[i]))
    elif sequenceType.lower() == 'nucl':
        with open(outputFileName, "a") as output:
            for i in range(len(seqIDs)):
                output.write("{0}\n{1}\n".format(seqIDs[i], outputNucs[i]))
    else:
        with open(protOutName, "a") as protFile, open(nuclOutName, "a") as nuclFile:
            for i in range(len(seqIDs)):
                protFile.write("{0}\n{1}\n".format(seqIDs[i], outputProts[i]))
                nuclFile.write("{0}\n{1}\n".format(seqIDs[i], outputNucs[i]))

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
p.add_argument("-tr", "--translation", dest="translationTable", type=int, default=1,
           help="Optionally specify the NCBI numeric genetic code to utilise for CDS translation (if relevant); this should be an integer from 1 to 31 (default == 1 i.e., Standard Code)")
p.add_argument("-th", "--threads", type=int, dest="threads",
           help="Number of threads to run. Default == 1, recommended == n-1.", default=1)
# SignalP opts
p.add_argument("-sdir", "--signalpdir", dest="signalpdir", type=str,
           help="""Specify the directory where signalp executables are located.""")
p.add_argument("-sorg", "--signalporg", dest="signalporg", type = str, choices = ['euk', 'gram-', 'gram+'], default='euk',
           help="""Specify the type of organism for SignalP from the available
           options. Refer to the SignalP manual if unsure what these mean (default == 'euk').""")
p.add_argument("-c", "--cygwindir", dest="cygwindir", type=str, default="",
           help="""Cygwin is required since you are running this program on a Windows computer.
           Specify the location of the bin directory here or, if this is already in your PATH, you can leave this blank."""
           if platform.system() == 'Windows' else argparse.SUPPRESS)

args = p.parse_args()
args = validate_args(args)

### RATIONALE FOR UNRESOLVED REGIONS ###
# I had to decide how to handle unresolved regions. I believe there are two valid approaches to this. The first is to replace any unresolved regions with stop codons and let the rest of the script process it like normal.
# This appears to be what NCBI does for their ORF Finder. The benefits of this approach is that you can't accidentally form chimeras. However, I think there is a second approach that has merit. If you are working with a genome that has
# short and rare occurrences of unresolved regions, you might not want to split up large ORFs on the basis of them having a very short stretch of unresolved codons. In this case, we can choose to not replace unresolved regions with stop codons
# if the region is shorter than an arbitrary and small limit. This isn't exactly perfect since even very short unresolved regions might hide stop codons. Thus, this option should be OFF by default, but we can provide the users an
# on/off switch so they can decide whether they want to risk discovering chimeric ORFs - providing a text prompt if the option is switched ON to verify any ORFs with unresolved regions would wash my hands of any mistake.

if args.unresolvedCodon != 0:
    print('Program has noted that you are allowing the discovery of ORFs with unresolved codon regions. This is risky behaviour, since this program cannot guarantee that an unresolved region does not contain a stop codon. Subsequently, you can have chimeras form from two separate ORFs. YOU MUST VERIFY ANY ORFS WITH UNRESOLVED REGIONS! The best way to do this is with BLAST against homologous proteins. You have been warned.')

# Load the fasta file as a generator object, get the total number of sequences in the file, then re-load it for the upcoming loop
records = SeqIO.parse(open(args.fileName, 'r'), 'fasta')
totalCount = 0
for record in records:
    totalCount += 1
records = SeqIO.parse(open(args.fileName, 'r'), 'fasta')

### CORE PROCESSING LOOP

# Set up queues
recordQ = queue.Queue(maxsize=50) # Arbitrary size; attempt to limit memory usage
outputQ = queue.Queue(maxsize=10000) # Arbitary size; attempt to limit memory usage

# Start threads
for i in range(args.threads):
    worker = Thread(target= record_worker, args=(recordQ, outputQ, args.translationTable, args.unresolvedCodon, args.minProLen, args.maxProLen, args.signalpdir, args.signalporg, args.cygwindir))
    worker.setDaemon(True)
    worker.start()
outputWorker = Thread(target = output_worker, args=(outputQ, totalCount, args.outputFileName, args.sequenceType, args.protOutName, args.nuclOutName))
outputWorker.setDaemon(True)
outputWorker.start()

# Put record in queue for worker threads
print('Starting the core processing of this script now. Progress bar will display below in a moment.')
for record in records:
    recordQ.put(record)
records.close()

# Close up shop on the threading structures
for i in range(args.threads):
    recordQ.put(None) # Add marker for record_workers to end
recordQ.join()
outputQ.put(None) # Add marker for output_worker to end
outputQ.join()

#### SCRIPT ALL DONE
print('Program completed successfully!')
