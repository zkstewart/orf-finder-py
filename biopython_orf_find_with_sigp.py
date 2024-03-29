#! python3
# Biopython based ORF Finder - SignalP version
# This script will obtain open reading frames from a fasta-formatted file containing nucleotide transcripts.
# Constraints can be altered to vary the strictness with which we accept or reject alternative start codons. 
# The output is a fasta file containing any number of ORFs that the user has specified based upon a minimum and maximum length also specified.
# Addition: This script will seek out ORFs with signal peptide start sites

# Load packages
import re, os, argparse, time, hashlib, queue, random
from Bio import SeqIO
from threading import Thread

# Load package from submodule
from Various_scripts import ZS_SignalPIO, ZS_SeqIO

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
    if args.signalpExe == None:
        print('signalpExe argument was not specified; fix your input and try again.')
        quit()
    if not os.path.isfile(os.path.join(args.signalpExe)):
        print('I am unable to locate the signalp execution file at the specified location (' + args.signalpExe + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
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

def run_signalp_sequence(signalpExe, organism, seqID, protString):
    # Coerce inputs into list
    if isinstance(seqID, str):
        seqID = [seqID]
    if isinstance(protString, str):
        protString = [protString]
    assert len(seqID) == len(protString), \
        "IDs and sequences lists don't match!"
    
    # Generate temporary file for sequence
    while True:
        tmpHash = ZS_SeqIO.Conversion.get_hash_for_input_sequences("".join([pStr[0:50] for pStr in protString]))
        tmpFileName = tmp_file_name_gen('tmp_sigpInput_', '.fasta', tmpHash, fullHashing=False)
        if os.path.isfile(tmpFileName):
            continue
        else:
            break
    
    with open(tmpFileName, 'w') as fileOut:
        for i in range(len(seqID)):
            fileOut.write('>' + seqID[i].lstrip('>') + '\n' + protString[i] + '\n')      # lstrip any > characters just in case they're already present
    
    # Run signalP
    sigpRunner = ZS_SignalPIO.SignalP(tmpFileName, signalpExe)
    sigpRunner.organism = organism
    sigPredictions = sigpRunner.signalp(withScore=True)
    
    # Clean up temporary files
    os.remove(tmpFileName)
    
    # Return signalP prediction dictionary
    return sigPredictions

def tmp_file_name_gen(prefix, suffix, hashString, fullHashing=True):
    # Main function
    if fullHashing:
        tmpHash = hashlib.sha256(bytes(str(hashString) + str(time.time()) + str(random.randint(0, 100000)), 'utf-8') ).hexdigest()       # This should always give us something unique even if the string for hashString is the same across different runs
    else:
        tmpHash = hashString
    while True:
        if os.path.isfile(prefix + tmpHash + suffix):
            tmpHash += 'X'
        else:
            return prefix + tmpHash + suffix

def orf_find_from_record(record, translationTable, unresolvedCodonLen, minProLen, maxProLen, signalpExe, signalporg):
    xRegex = re.compile(r'X+')                          # Regex used to find start and stop positions of unresolved regions that are shorter than the cut-off
    
    # Grab all "good" start sites for checking with signalP
    indices = []
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
            for i in range(len(splitProtein)):
                
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
                
                # Store start sites for batched signal peptide presence prediction
                for x in range(len(startIndices)):
                    indices.append(startIndices[x])
                    prots.append(acceptedPro[startIndices[x][0]:])
                    nucs.append(nucSeqOfProt[startIndices[x][0]*3:])
    
    # Run signalP prediction
    ids = [ str(i) for i in range(len(prots)) ]
    sigpPredictions = run_signalp_sequence(signalpExe, signalporg, ids, prots)
    if sigpPredictions == {}:
        return None, None
    
    # Reduce our indices, prots, and nucs lists to just those with hits
    indices = [ indices[int(id)] for id in sigpPredictions ]
    prots = [ prots[int(id)] for id in sigpPredictions ]
    nucs = [ nucs[int(id)] for id in sigpPredictions ]
    scores = [ sigpPredictions[id][2] for id in sigpPredictions  ] # also get the sigp prediction score
    
    # Immediately return value if there's only one hit
    if len(prots) == 1:
        return prots, nucs
    
    # Combine our index, prot, nuc, and score values
    indexProtNucs = list(zip(indices, prots, nucs, scores))
    
    # Sort candidates on the basis of evidence canonical > startSite > sigp score
    indexProtNucs.sort(key = lambda x: (x[0][1], -len(x[1]), -x[3])) # i.e., [canonical(0 is best), sequence length(longer is better), score (higher is better)]
    
    # Identify the best candidate with a signalP prediction
    bestCandidate = indexProtNucs[0]
    EXTENSION_FOR_BETTER_SIGP_SCORE = 7 # Arbitrary; this should help to address situations where two M's are close together, and one has a better signalP prediction score
    STEPPING_STONE_EXTENSION_LIMITER = 3 # Arbitrary; after extending for a better sigp score once, we enforce a stricter mechanism for extensions afterwards
    for indexPN in indexProtNucs: ## TBD: Include extension
        # If new candidate has a better start site than bestCandidate...
        if indexPN[0][1] < bestCandidate[0][1]: # this measures canonical (0 is best)
            bestCandidate = indexPN
        # If new candidate has an equivalent start site...
        elif indexPN[0][1] == bestCandidate[0][1]:
            # ...and it is much longer than bestCandidate...
            if (len(bestCandidate[1]) + EXTENSION_FOR_BETTER_SIGP_SCORE) < len(indexPN[1]):
                bestCandidate = indexPN
            # ...or it isn't much shorter than bestCandidate...
            elif (len(bestCandidate[1]) - EXTENSION_FOR_BETTER_SIGP_SCORE) <= len(indexPN[1]):
                # ...and its signalP score is better than bestCandidate...
                if indexPN[3] > bestCandidate[3]: # this measures signalpep (1 is best)
                    bestCandidate = indexPN
                    EXTENSION_FOR_BETTER_SIGP_SCORE = STEPPING_STONE_EXTENSION_LIMITER # Enforce stricter limit now that we've extended once
    
    return [bestCandidate[1]], [bestCandidate[2]] # function is expected to return a list

def record_worker(recordQ, outputQ, translationTable, unresolvedCodon, minProLen, maxProLen, signalpExe, signalporg):
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
        prots, nucs = orf_find_from_record(record, translationTable, unresolvedCodon, minProLen, maxProLen, signalpExe, signalporg)
        if prots != None:
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

def main():

    ### USER INPUT
    usage = """%(prog)s reads in a fasta formatted file containing nucleotide sequences and, following user-specified parameters,
    produces an output fasta file containing potential open reading frames (ORFs) as nucleotides/protein translations/both.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="fileName",
                required=True,
                help="Input fasta file name")
    p.add_argument("-o", dest="outputFileName",
                required=True,
                help="Output fasta file name")
    p.add_argument("-s", dest="signalpExe",
                required=True,
                type=str,
                help="Specify the location of the signalp executable.")
    # Opts
    p.add_argument("-min", "--minimum", dest="minProLen",
                required=False,
                type=int,
                help="Minimum ORF amino acid length. Default == 30.",
                default=30)
    p.add_argument("-max", "--maximum", dest="maxProLen",
                required=False,
                type=int,
                help="""Optional specification for maximum ORF amino acid length. Default == 0,
                which means there is no upper limit.""",
                default=0)
    p.add_argument("-st", "--seqtype", dest="sequenceType",
                required=False,
                choices = ['prot', 'nucl', 'both', 'PROT', 'NUCL', 'BOTH'],
                help="""Specify the type of output you want to generate (i.e., protein translated ORF,
                nucleotide CDS, or both). If you specify 'both', two outputs with '_prot' and '_nucl'
                suffix will be generated. Default == 'prot'.""",
                default="prot")
    p.add_argument("-f", "--force", dest="force",
                required=False,
                action="store_true",
                help="Allow files to be overwritten at your own risk.",
                default=False)
    p.add_argument("-u", "--unresolved", dest="unresolvedCodon",
                required=False,
                type=int,
                help="""Default == 0, which means the program will not discover ORFs with unresolved codons.
                If you want to risk chimeric ORF formation, you can change this value. You MUST validate any
                ORFs with unresolved portions. Recommended for this value to be less than 5.""",
                default=0)
    p.add_argument("-tr", "--translation", dest="translationTable",
                required=False,
                type=int,
                help="""Optionally specify the NCBI numeric genetic code to utilise for CDS translation (if relevant);
                this should be an integer from 1 to 31 (default == 1 i.e., Standard Code)""",
                default=1)
    p.add_argument("-th", "--threads", dest="threads",
                required=False,
                type=int, 
                help="Number of threads to run. Default == 1, recommended == n-1.",
                default=1)
    # SignalP opts
    p.add_argument("-sorg", "--signalporg", dest="signalporg",
                required=False,
                type=str,
                choices=['euk', 'gram-', 'gram+'],
                help="""Specify the type of organism for SignalP from the available
                options. Refer to the SignalP manual if unsure what these mean (default == 'euk').""",
                default='euk')
    
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
    
    # Load the fasta file as a generator object, get the total number of sequences in the file, and validate that they're nucleotides
    with open(args.fileName, 'r') as fileIn:
        records = SeqIO.parse(fileIn, 'fasta')
        totalCount = 0
        for record in records:
            totalCount += 1
    
    ### CORE PROCESSING LOOP
    # Set up queues
    recordQ = queue.Queue(maxsize=50) # Arbitrary size; attempt to limit memory usage
    outputQ = queue.Queue(maxsize=10000) # Arbitary size; attempt to limit memory usage
    
    # Start threads
    for i in range(args.threads):
        worker = Thread(target=record_worker,
                        args=(recordQ, outputQ, args.translationTable,
                            args.unresolvedCodon, args.minProLen, args.maxProLen,
                            args.signalpExe, args.signalporg))
        worker.setDaemon(True)
        worker.start()
    
    outputWorker = Thread(target=output_worker,
                        args=(outputQ, totalCount, args.outputFileName,
                                args.sequenceType, args.protOutName, args.nuclOutName))
    outputWorker.setDaemon(True)
    outputWorker.start()
    
    # Put records in queue for worker threads
    with open(args.fileName, 'r') as fileIn:
        records = SeqIO.parse(fileIn, 'fasta')
        print('Starting the core processing of this script now. Progress bar will display below in a moment.')
        for record in records:
            recordQ.put(record)
    
    # Close up shop on the threading structures
    for i in range(args.threads):
        recordQ.put(None) # Add marker for record_workers to end
    recordQ.join()
    
    outputQ.put(None) # Add marker for output_worker to end
    outputQ.join()
    
    #### SCRIPT ALL DONE
    print('Program completed successfully!')

if __name__ == "__main__":
    main()
