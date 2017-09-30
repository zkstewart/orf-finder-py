#! python3
# Biopython based ORF Finder
# This script will obtain open reading frames from a fasta-formatted file containing nucleotide transcripts.
# Constraints can be altered to vary the strictness with which we accept or reject alternative start codons. 
# The output is a .fasta file containing any number of ORFs that the user has specified based upon a minimum and maximum length also specified. 

# Load packages
import re, os, argparse
from Bio import SeqIO

### USER INPUT
usage = """%(prog)s reads in a fasta formatted file containing nucleotide sequences and, following user-specified parameters,
produces an output fasta file containing potential open reading frames (ORFs) as nucleotides/protein translations/both.
"""
# Reqs
p = argparse.ArgumentParser(description=usage)
#p.add_argument("input", type = str, help="Input fasta file name")
#p.add_argument("output", type = str, help="Output fasta file name")
p.add_argument("-i", "-input", dest="fileName",
                   help="Input fasta file name")
p.add_argument("-o", "-output", dest="outputFileName",
                   help="Output fasta file name")
# Opts
p.add_argument("-min", "-minimum", type=int, dest="minProLen",
                   help="Minimum ORF amino acid length. Default == 30.", default=30)
p.add_argument("-max", "-maximum", type=int, dest="maxProLen",
                   help="Optional specification for maximum ORF amino acid length. Default == 0, which means there is no upper limit.", default=0)
p.add_argument("-num", "-numhits", type=int, dest="hitsToPull",
                   help="Specify the number of ORFs you wish to extract from each sequence. Default == 3.", default=3)
p.add_argument("-alt", "-altcodon", type=int, dest="altCodonStringency",
                   help="Control the stringency with which alternative start codon ORFs are accepted. Recommended not to change unless you understand the influence this has. Default == 49.", default=49)
p.add_argument("-no", "-nocodon", type=int, dest="noCodonStringency",
                   help="Control the stringency with which fragmentary ORFs are accepted (fragmentary means there is no traditional or common alternative start codon in the sequence). Recommended not to change unless you understand the influence this has. Default == 99.", default=99)
p.add_argument("-st", "-seqtype", dest="sequenceType", choices = ['prot', 'nucl', 'both', 'PROT', 'NUCL', 'BOTH'],
                   help="Specify the type of output you want to generate (i.e., protein translated ORF, nucleotide CDS, or both). If you specify 'both', two outputs with '_prot' and '_nucl' suffix will be generated. Default == 'prot'.", default="prot")
p.add_argument("-r", "-replace", dest="replace", choices = ['y', 'n', 'Y', 'N'],
                   help="Optional ability to replace alternative starting position with a methionine (M) [only relevant if obtaining proteins]. Default == 'n'.", default='n')
p.add_argument("-f", "-force", dest="force", choices = ['y', 'n', 'Y', 'N'],
                   help="Default == 'n', which means the program will not overwrite existing files. Specify 'y' to allow this behaviour at your own risk.", default='n')

args = p.parse_args()

fileName = args.fileName
outputFileName = args.outputFileName
minProLen = args.minProLen
maxProLen = args.maxProLen
hitsToPull = args.hitsToPull
altCodonStringency = args.altCodonStringency
noCodonStringency = args.noCodonStringency
sequenceType = args.sequenceType
replace = args.replace
force = args.force

# Check if we should be overwriting files
if outputFileName != None:
        if os.path.isfile(outputFileName) and force.lower() != 'y':
                print('There is already a file named ' + outputFileName + '. Either specify a new file name, delete this older file, or provide the -force argument either "Y" or "y"')
                quit()
        elif os.path.isfile(outputFileName) and force.lower() == 'y':
                os.remove(outputFileName)

if fileName == None or outputFileName == None:
        # Locate our file of interest
        print('Enter the name of the fasta-formatted file you wish to extract ORFs from. This should be a nucleotide sequence file. Include file extension (e.g., ".fasta").')
        while True:
                try:
                        fileName = input()
                        if os.path.isfile(fileName) == False:
                                raise Exception
                        print('Fasta file located successfully')
                        break
                except KeyboardInterrupt:
                        quit()
                except:
                        print('Fasta file failed to load. If you misspelt the name, try again. If this script isn\'t in the same directory as the fasta file, move it there then try again.')
                        continue
        print('')

        # Allow user to determine output file name
        print('Enter the name which you want the output fasta file to be called. Include the file extension, and make sure not to use illegal characters (i.e. \\/:?"<>|).')
        while True:
                try:
                        illegalCharacters = '\/:?"<>|'
                        outputFileName = input()
                        if outputFileName == '':
                                print('You didn\'t name this file anything. You need to have at least one character in your output file name. Try again.')
                                continue
                        for character in illegalCharacters:
                             if character in outputFileName:
                                raise Exception
                        if os.path.isfile(outputFileName):
                                print('This is already a file. Try again.')
                                continue
                        break
                except KeyboardInterrupt:
                        quit()
                except:
                        print('You used an illegal character (i.e. \\/:?"<>|). Try to name your file without these characters again.')
                        continue
        print('')

        # Allow user to determine minimum protein length
        print('Enter the minimum amino acid length that you want to accept as a legitimate ORF. A setting of 100AA (equivalent to 300bp nucleotide) is recommended if searching for \'normal\' proteins, but you can reduce this number if you want to find smaller proteins.')
        while True:
                minProLen = input()
                try:
                        minProLen = int(minProLen)
                        break
                except KeyboardInterrupt:
                        quit()
                except:
                        print('You seem to have not typed a number here. Try again.')
        print('')

        # Allow user to determine maximum protein length
        print('Enter the maximum amino acid length that you want to accept as a legitimate ORF. A setting of 0AA will result in no upper limit.')
        while True:
                maxProLen = input()
                try:
                        maxProLen = int(maxProLen)
                        break
                except KeyboardInterrupt:
                        quit()
                except:
                        print('You seem to have not typed a number here. Try again.')
        print('')

        # Allow user to determine number of ORFs to pull
        print('Enter the number of ORFs you wish to extract from each sequence. Recommended number is from 1-3, but I won\'t force you to do that. If you want to pull all ORFs, just put in a large number (not too large though as it may reduce script speed).')
        while True:
                hitsToPull = input()
                try:
                        hitsToPull = int(hitsToPull)
                        break
                except KeyboardInterrupt:
                        quit()
                except:
                        print('You seem to have not typed a number here. Try again.')
        print('')

        # Allow user to customise script for different performance
        print('Do you want to run with the default settings of this program?')
        print('If you wish to run default (the recommended settings), type \'yes\', otherwise type \'no\' to customise these settings.')
        print('')
        print('Default means that we won\'t utilise a alternative start codon (i.e., TTG or GTG commonly, or a CTG rarely) unless it improves ORF length by 49 AA when compared to ORFs that start with a methionine (ATG), or an ORF which lacks a start codon entirely unless it improves ORF length by 99AA when compared to those that do have start codons.')
        print('')
        while True:
                yesOrNo = input()
                if yesOrNo.lower() == 'yes':
                        altCodonStringency = 49
                        noCodonStringency = 99
                        print('Good choice! Default is the best way to obtain results that are similar to or better than established ORF finders such as NCBI\'s ORF Finder.')
                        print('')
                        break
                if yesOrNo.lower() == 'no':
                        print('Have it your way! What do you want the length requirement to be before we accept a common alternative start codon (i.e. TTG or GTG).')
                        print('Note that, because of the way this script works, you should minus one from the number you want to use (i.e. 100 becomes 99).') # Numbers should be - 1 since, if we compare a protein to a no hit (aka '-' value) in the "Run a final size comparison to choose the best ORF" section below, using a minProLen of 100 instead of 99 will result in a protein with a length of 100A not being properly accepted in fringe cases where the length of the protein is exactly 100AA (our default minimum length). Reasoning: if 100 > len('-') + minProLen(99)... if minProLen were to be 100, then our protein of length 100 would not be > 101
                        while True:
                                altCodonStringency = input()
                                try:
                                        altCodonStringency = int(altCodonStringency)
                                        break
                                except KeyboardInterrupt:
                                        quit()
                                except:
                                        print('You seem to have not typed a number here. Try again.')
                        
                        print('What do you want the length requirement to be before we accept potentially incompletely sequenced ORFs with no start codon (i.e. at the 5\' and 3\' ends of the nucleotide sequence).')
                        while True:
                                noCodonStringency = input()
                                try:
                                        noCodonStringency = int(noCodonStringency)
                                        break
                                except KeyboardInterrupt:
                                        quit()
                                except:
                                        print('You seem to have not typed a number here. Try again.')
                        break
                else:
                        print('You didn\'t seem to type \'yes\' or \'no\', try again.')

        # Allow user to determine what kind of output they want to have
        print('Do you want to produce a translated protein ORF file, a nucleotide CDS file, or both? Enter \'prot\', \'nucl\', or \'both\'')
        seqtypeChoices = ['prot', 'nucl', 'both', 'PROT', 'NUCL', 'BOTH']
        while True:
                try:
                        sequenceType = input()
                        if sequenceType.lower() not in seqtypeChoices:
                                raise Exception
                        break
                except KeyboardInterrupt:
                        quit()
                except:
                        print('You didn\'t type \'prot\', \'nucl\', or \'both\'. Try again.')
        print('')

        # Allow user to determine whether alternative starts should be coded as M
        print('Do you want to replace alternative start positions with a methionine (M)? Enter \'y\' or \'n\'')
        metChoices = ['y', 'n', 'Y', 'N']
        while True:
                try:
                        replace = input()
                        if replace.lower() not in metChoices:
                                raise Exception
                        break
                except KeyboardInterrupt:
                        quit()
                except:
                        print('You didn\'t type \'y\' or \'n\'. Try again.')
        print('')

# Load the fasta file as a generator object
records = SeqIO.parse(open(fileName, 'rU'), 'fasta')

# Sort output file names if outputting both prot and nucl
if sequenceType.lower() == 'both':
        outPrefix = outputFileName.rsplit('.', maxsplit=1)
        protOutName = outPrefix[0] + '_prot.' + outPrefix[1]
        nuclOutName = outPrefix[0] + '_nucl.' + outPrefix[1]

### CORE PROCESSING LOOP
print('Starting the core processing of this script now. Notifications will be presented in this text area.')

# Declare overall values needed before loop start
startCodon = re.compile(r'^.*?(M.*)')           # Regex to pull out just the sequence starting with a Methionine (or ATG)
ongoingCount = 0
outputProt = []                                 # These get reset whenever we output to file
outputNucl = []

# Get the nucleotide (record) out of our list of nucleotides (records) and grab the ORFs!
for record in records:
        # Declare output holding values that should reset for each transcript/record
        tempOverallProt = []
        tempOverallNucl = []
        tempMProt = []
        tempMNucl = []
        tempAltProt = []
        tempAltNucl = []
        tempNoneProt = []
        tempNoneNucl = []
        for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
                for frame in range(3):
                        length = 3 * ((len(record)-frame) // 3)
                        ongoingLength = 0                                                                       # The ongoingLength will track where we are along the splitProtein sequence to determine the protPosition
                        splitProtein = nuc[frame:frame+length].translate(table=1).split("*")                    # An asterisk represents a stop codon, thus by splitting by asterisk we obtain a list of all the ORFs in this frame
                        for i in range(len(splitProtein)):                                                      # Note that I have done a 'for i in range...' loop rather than a 'for value in splitProtein' loop which would have been simpler for a reason explained below on the 'elif i + 1 ==' line
                                # Declare blank values needed for each potential ORF region so we can tell which things were 'found'
                                mPro = ''
                                altPro = ''
                                nonePro = ''
                                topHit = ''
                                codonIndex = None
                                noneCodonContingency = None
                                # Determine ongoingLength before we continue with the actual processing 
                                if len(splitProtein) == 1:                                      # This means the splitProtein has no stop codons
                                        ongoingLength += len(splitProtein[i])
                                elif i + 1 == len(splitProtein):                                # This means we're looking at the last ORF in the splitProtein, which means it will have no stop codon
                                        ongoingLength += len(splitProtein[i])
                                else:                                                           # Add +1 here since the protein does have a stop codon after it. I could ignore this step, but I think a stop codon should be considered in the protein length since it's a conserved feature. Mainly, this weights an ORF with a stop codon above one without assuming the amino-acid length is identical.
                                        ongoingLength += len(splitProtein[i]) + 1
                                # Process sequences to determine whether we're ignoring this, or adding an asterisk for length counts
                                if len(splitProtein[i]) < minProLen:                            # Disregard sequences that won't meet the size requirement without further processing
                                        continue
                                elif maxProLen != 0 and len(splitProtein[i]) > maxProLen:       # Disregard sequences that won't meet the size requirement without further processing
                                        continue
                                elif len(splitProtein) == 1:                            # If the length of splitProtein is 1, that means there was no stop codon in this frame
                                        acceptedPro = str(splitProtein[i])
                                elif i + 1 == len(splitProtein):                        # [i + 1 will equal the length of splitProtein on the final loop] Since we know there are at least 2 values in the splitProtein list, we know that the last entry will not have a stop codon. This line is the reason I used the 'for i in range' strategy rather than 'for value in splitProtein', since, depending on how short the minProLen has been specified as, it was possible doing a 'if value == splitProtein[-1]' would not accurately determine whether the ORF region we're looking at is actually the final region or is simply identical to the final region
                                        acceptedPro = str(splitProtein[i])
                                else:                                                   # Since we know there are at least 2 values in the splitProtein list, and we know that this is not the last entry, we know that this entry has a stop codon after it
                                        acceptedPro = str(splitProtein[i]) + '*'
                                # Alternative start coding      
                                preProtLength = (ongoingLength - len(acceptedPro))*3        # Get the length of the region leading up to the protein in nucleotides (hence times 3). We also minus the length of the current splitProtein since we have added this to the ongoingLength value already
                                protLenAsNucl = len(acceptedPro)*3
                                nucSeqOfProt = nuc[frame:frame+length][preProtLength:preProtLength+protLenAsNucl]  # Pull out the nucleotide sequence that corresponds to the ORF region
                                codons = re.findall('..?.?', str(nucSeqOfProt))         # Pulls out a list of codons from the nucleotide
                                for codon in codons:                                    # Cycle through this list of codons to find the first alternative start of the normal class (GTG and TTG) and the rare class (CTG)
                                        if codon == 'GTG' or codon == 'TTG':
                                                codonIndex = codons.index(codon)        # This will save the position of the first GTG or TTG encountered. Note that by breaking after this,  we stop looking for CTG as it is irrelevant after this
                                                break
                                        elif codon == 'CTG':
                                                if noneCodonContingency == None:        # noneCodonContingency is set to None at the end of each loop. Thus, this line of code will 'capture' the position of the first CTG in a sequence if a GTG or TTG was not encountered first
                                                        noneCodonContingency = codons.index(codon)

                                # Get the three ORF versions from each region inbetween stop codons
                                if 'M' in str(acceptedPro):                             # Obtains a traditional methionine initiated ORF starting from the first methionine if there is one in the sequence
                                        mPro = startCodon.search(str(acceptedPro)).groups()[0]  # Note that startCodon was declared at the start of this file         

                                if codonIndex != None:                                  # Gets the start position of the protein if we found a likely alternative start (aka a 'GTG' or 'TTG')
                                        altPro = acceptedPro[codonIndex:]
                                        if replace.lower() == 'y':
                                                altPro = 'M' + altPro[1:]               # This script assumes that the alternative start will be substituted with a methionine post-transcription
                                elif noneCodonContingency != None:                      # This will match an alternative start to 'CTG' only if 'TTG' or 'GTG' are not present
                                        altPro = acceptedPro[codonIndex:]
                                        if replace.lower() == 'y':
                                                altPro = 'M' + altPro[1:]

                                if splitProtein[i] == splitProtein[0]:                  # nonePro makes an assumption that the start of the ORF was not assembled properly resulting in the real start codon being cut off. Our stringency values will assess the likelihood of this hypothesis.
                                        nonePro = acceptedPro                           # Additionally, by only obtaining a 'nonePro' when it is in a protein fragment at the start of a frame (i.e., splitProtein[0]), we also make the (reasonable) assumption that any ORF inbetween two stop codons should itself have a start codon. This doesn't always hold true due to transcript assembly errors, but it must be assumed for the purpose of this script.                          

                                # Pull out the top hit from this protein fragment based upon how strict we want to be with accepting a traditional, alternative, or no codon start                                
                                if len(nonePro) > len(altPro) + noCodonStringency and len(nonePro) > len(mPro) + noCodonStringency:             # By adding on the stringency values declared earlier to the length of the protein, we can determine whether we want to consider an ORF without a start codon as legitimate
                                        topHit = nonePro
                                elif len(altPro) > len(mPro) + altCodonStringency:                                                              # Adding the stringency values here allows us to determine whether we can increase ORF length significantly by assuming an alternative start rather than a methionine start
                                        topHit = altPro
                                else:
                                        topHit = mPro                                                                                           # This is the default position unless either an alternative start or a no codon start outweighs the stringency values

                                # Cull the top hit if it doesn't meet our minimum length requirement anymore, or add it to the temporary list of ORF hits from this nucleotide sequence
                                if len(topHit) < minProLen:                                             # Culling is necessary since we will have shortened the sequence somewhat unless we accepted the topHit as being a no codon start ORF. Note that we will here consider a stop codon in the length of the protein, such that a protein with 99AAs and a stop will pass a minimum 100AA length test. I think this is fair since not all regions here have a stop codon which allows weight to be added to these cases, especially since a stop codon is still conserved as part of an ORF.
                                        doNothing = ''
                                elif maxProLen!= 0 and len(topHit) > maxProLen:                         # Culling should not be necessary here in almost all scenarios, but who knows?
                                        doNothing = ''
                                elif topHit == mPro:                                                    # These temp lists will be populated with potential ORFs from a single nucleotide sequence before being processed in the next major chunk of code starting with 'if len(tempMList + tempAltList + tempNoneList) >= 1:'
                                        if sequenceType.lower() == 'prot':
                                                tempMProt.append(topHit)
                                        elif sequenceType.lower() == 'nucl':
                                                newStartPosition = acceptedPro.find(mPro)
                                                tempMNucl.append(str(nucSeqOfProt[newStartPosition*3:]))
                                        else:
                                                tempMProt.append(topHit)
                                                newStartPosition = acceptedPro.find(mPro)
                                                tempMNucl.append(str(nucSeqOfProt[newStartPosition*3:]))
                                elif topHit == altPro:
                                        if sequenceType.lower() == 'prot':
                                                tempAltProt.append(topHit)
                                        elif sequenceType.lower() == 'nucl':
                                                newStartPosition = acceptedPro.find(altPro[1:]) - 1     # - 1 since we're looking at the second character in our altPro (just in case we're replacing with M)
                                                tempAltNucl.append(str(nucSeqOfProt[newStartPosition*3:]))
                                        else:
                                                tempAltProt.append(topHit)
                                                newStartPosition = acceptedPro.find(altPro[1:]) - 1
                                                tempAltNucl.append(str(nucSeqOfProt[newStartPosition*3:]))
                                elif topHit == nonePro:
                                        if sequenceType.lower() == 'prot':
                                                tempNoneProt.append(topHit)
                                        elif sequenceType.lower() == 'nucl':
                                                tempNoneNucl.append(str(nucSeqOfProt))
                                        else:
                                                tempNoneProt.append(topHit)
                                                tempNoneNucl.append(str(nucSeqOfProt))

        # Sort our top hits from each inter-stop codon fragment by size and category (i.e. mPro or altPro?) and select the top X hits
        if len(tempMProt + tempAltProt + tempNoneProt) >= 1 or len(tempMNucl + tempAltNucl + tempNoneNucl) >= 1:
                # Append '-' entries to lists which have less entries than we want to pull to allow the below 'for' loops to run without exceptions
                        # Prot list     [If we are only looking at nucleotides, then prot lists will be populated with hyphens which, realistically, won't impact memory consumption
                for i in range(0, hitsToPull-len(tempMProt)):
                        tempMProt.append('-')
                for i in range(0, hitsToPull-len(tempAltProt)):
                        tempAltProt.append('-')
                for i in range(0, hitsToPull-len(tempNoneProt)):
                        tempNoneProt.append('-')
                        # Nucl list
                for i in range(0, hitsToPull-len(tempMNucl)):
                        tempMNucl.append('-')
                for i in range(0, hitsToPull-len(tempAltNucl)):
                        tempAltNucl.append('-')
                for i in range(0, hitsToPull-len(tempNoneNucl)):
                        tempNoneNucl.append('-')
                # Sort the lists by size (largest on the bottom to allow the .pop() method to remove a hit when accepted) [as above, the prot or nucl variants might just be lists of hyphens. Running the sort twice shouldn't realistically impact time in that case.
                        # Prot
                tempSortedMProt = sorted(tempMProt, key=len)
                tempSortedAltProt = sorted(tempAltProt, key=len)
                tempSortedNoneProt = sorted(tempNoneProt, key=len)
                        # Nucl
                tempSortedMNucl = sorted(tempMNucl, key=len)
                tempSortedAltNucl = sorted(tempAltNucl, key=len)
                tempSortedNoneNucl = sorted(tempNoneNucl, key=len)
                # Run a final size comparison to choose the best ORF(s). We need to split this into two separate statements since the stringency values need to be *3 for nucls
                if sequenceType.lower() == 'prot' or sequenceType.lower() == 'both':
                        for i in range(0, hitsToPull):
                                if len(tempSortedNoneProt[-1]) > len(tempSortedAltProt[-1]) + noCodonStringency and len(tempSortedNoneProt[-1]) > len(tempSortedMProt[-1]) + noCodonStringency:         # Again, we add the stringency values to help with determining priority of ORF ordering. Since this script will often be returning either 1, 3, or 5 potential ORFs, it is important that we order these in the most logical way
                                        tempOverallProt.append(tempSortedNoneProt[-1])
                                        tempSortedNoneProt.pop()
                                elif len(tempSortedAltProt[-1]) > len(tempSortedMProt[-1]) + altCodonStringency:
                                        tempOverallProt.append(tempSortedAltProt[-1])
                                        tempSortedAltProt.pop()
                                else:
                                        tempOverallProt.append(tempSortedMProt[-1])                                                                                                                     # By using this as the 'else' position, sequences with methionine starts will be selected in the majority of situations as the default stringency settings ensure that it is rare an alternative start is used instead of a methionine start
                                        tempSortedMProt.pop()
                if sequenceType.lower() == 'nucl' or sequenceType.lower() == 'both':
                        for i in range(0, hitsToPull):
                                if len(tempSortedNoneNucl[-1]) > len(tempSortedAltNucl[-1]) + noCodonStringency*3 and len(tempSortedNoneNucl[-1]) > len(tempSortedMNucl[-1]) + noCodonStringency*3:         # Again, we add the stringency values to help with determining priority of ORF ordering. Since this script will often be returning either 1, 3, or 5 potential ORFs, it is important that we order these in the most logical way
                                        tempOverallNucl.append(tempSortedNoneNucl[-1])
                                        tempSortedNoneNucl.pop()
                                elif len(tempSortedAltNucl[-1]) > len(tempSortedMNucl[-1]) + altCodonStringency*3:
                                        tempOverallNucl.append(tempSortedAltNucl[-1])
                                        tempSortedAltNucl.pop()
                                else:
                                        tempOverallNucl.append(tempSortedMNucl[-1])                                                                                                                     # By using this as the 'else' position, sequences with methionine starts will be selected in the majority of situations as the default stringency settings ensure that it is rare an alternative start is used instead of a methionine start
                                        tempSortedMNucl.pop()
                # Format and produce the output of this script
                tempOutputProt = []
                tempOutputNucl = []
                if sequenceType.lower() == 'prot' or sequenceType.lower() == 'both':
                        for i in range(0, hitsToPull):                          
                                if tempOverallProt[i] == '-':                   # Because we made sure all the tempM/Alt/NoneLists had '-' added to pad out the list to have a length equal to the value of hitsToPull, when we cycle through our tempOverallList, we will often encounter '-' characters which signify the end of relevant ORFs identified in this sequence  
                                        break                                   # Break out of this loop once we've fasta formatted all relevant ORF hits
                                tempOutputProt.append('>' + record.id + '_ORF' + str(i+1) + '\n' + tempOverallProt[i])
                        if len(tempOutputProt) > 0:    # Need this check for sequences that don't get any hits
                                outputProt.append('\n'.join(tempOutputProt))

                if sequenceType.lower() == 'nucl' or sequenceType.lower() == 'both':
                        for i in range(0, hitsToPull):                          
                                if tempOverallNucl[i] == '-':                   # Because we made sure all the tempM/Alt/NoneLists had '-' added to pad out the list to have a length equal to the value of hitsToPull, when we cycle through our tempOverallList, we will often encounter '-' characters which signify the end of relevant ORFs identified in this sequence  
                                        break                                   # Break out of this loop once we've fasta formatted all relevant ORF hits
                                tempOutputNucl.append('>' + record.id + '_ORF' + str(i+1) + '\n' + tempOverallNucl[i])
                        if len(tempOutputNucl) > 0:
                                outputNucl.append('\n'.join(tempOutputNucl))
        else:
                doNothing = ''                # We don't need to do anything if no hits were found that pass the minimum length threshold

        ongoingCount += 1
        
        # Save backup if ongoingCount == 10,000. There are two sections here, the first will create the file on the first loop, the second will add to the file on subsequent loops
        if sequenceType.lower() == 'prot':
                if ongoingCount%10000 == 0 and os.path.isfile(os.getcwd() + '\\' + outputFileName) == False:
                        with open(outputFileName, 'w') as output:
                                output.write('\n'.join(outputProt))
                        print(str(ongoingCount) + ' sequences scanned for ORFs.')
                        outputProt = []
                elif ongoingCount%10000 == 0 and os.path.isfile(os.getcwd() + '\\' + outputFileName) == True:
                        with open(outputFileName, 'a') as output:
                                output.write('\n')
                                output.write('\n'.join(outputProt))
                        print(str(ongoingCount) + ' sequences scanned for ORFs.')
                        outputProt = []
        elif sequenceType.lower() == 'nucl':
                if ongoingCount%10000 == 0 and os.path.isfile(os.getcwd() + '\\' + outputFileName) == False:
                        with open(outputFileName, 'w') as output:
                                output.write('\n'.join(outputNucl))
                        print(str(ongoingCount) + ' sequences scanned for ORFs.')
                        outputNucl = []
                elif ongoingCount%10000 == 0 and os.path.isfile(os.getcwd() + '\\' + outputFileName) == True:
                        with open(outputFileName, 'a') as output:
                                output.write('\n')
                                output.write('\n'.join(outputNucl))
                        print(str(ongoingCount) + ' sequences scanned for ORFs.')
                        outputNucl = []
        else:
                if ongoingCount%10000 == 0 and os.path.isfile(os.getcwd() + '\\' + protOutName) == False:       # Doesn't matter if we check for protOutName or nuclOutName. Theoretically, if one of the files is deleted while the program is running it could be problematic, but I mean, what can I really do about that without child-proofing the script excessively?
                        with open(protOutName, 'w') as protFile, open(nuclOutName, 'w') as nuclFile:
                                protFile.write('\n'.join(outputProt))
                                nuclFile.write('\n'.join(outputNucl))
                        outputProt = []
                        outputNucl = []
                        print(str(ongoingCount) + ' sequences scanned for ORFs.')
                elif ongoingCount%10000 == 0 and os.path.isfile(os.getcwd() + '\\' + protOutName) == True:
                        with open(protOutName, 'a') as protFile, open(nuclOutName, 'a') as nuclFile:
                                protFile.write('\n')
                                nuclFile.write('\n')
                                protFile.write('\n'.join(outputProt))
                                nuclFile.write('\n'.join(outputNucl))
                        print(str(ongoingCount) + ' sequences scanned for ORFs.')
                        outputProt = []
                        outputNucl = []

# Dump the last few results after the script has finished, or create the output if there were less than 10,000 sequences
if sequenceType.lower() == 'prot':
        if os.path.isfile(os.getcwd() + '\\' + outputFileName) == False:
                with open(outputFileName, 'w') as output:
                        output.write('\n'.join(outputProt))
                print(str(ongoingCount) + ' sequences scanned for ORFs.')
        elif os.path.isfile(os.getcwd() + '\\' + outputFileName) == True:
                with open(outputFileName, 'a') as output:
                        output.write('\n')
                        output.write('\n'.join(outputProt))
                print(str(ongoingCount) + ' sequences scanned for ORFs.')
elif sequenceType.lower() == 'nucl':
        if os.path.isfile(os.getcwd() + '\\' + outputFileName) == False:
                with open(outputFileName, 'w') as output:
                        output.write('\n'.join(outputNucl))
                print(str(ongoingCount) + ' sequences scanned for ORFs.')
        elif os.path.isfile(os.getcwd() + '\\' + outputFileName) == True:
                with open(outputFileName, 'a') as output:
                        output.write('\n')
                        output.write('\n'.join(outputNucl))
                print(str(ongoingCount) + ' sequences scanned for ORFs.')
else:
        if os.path.isfile(os.getcwd() + '\\' + protOutName) == False:       # Doesn't matter if we check for protOutName or nuclOutName. Theoretically, if one of the files is deleted while the program is running it could be problematic, but I mean, what can I really do about that without child-proofing the script excessively?
                with open(protOutName, 'w') as protFile, open(nuclOutName, 'w') as nuclFile:
                        protFile.write('\n'.join(outputProt))
                        nuclFile.write('\n'.join(outputNucl))
                print(str(ongoingCount) + ' sequences scanned for ORFs.')
        elif os.path.isfile(os.getcwd() + '\\' + protOutName) == True:
                with open(protOutName, 'a') as protFile, open(nuclOutName, 'a') as nuclFile:
                        protFile.write('\n')
                        nuclFile.write('\n')
                        protFile.write('\n'.join(outputProt))
                        nuclFile.write('\n'.join(outputNucl))
                print(str(ongoingCount) + ' sequences scanned for ORFs.')

records.close()
#### SCRIPT ALL DONE
