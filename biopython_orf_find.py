#! python3
# Biopython based ORF Finder
# This script will obtain open reading frames from a fasta-formatted file containing nucleotide transcripts.
# Constraints can be altered to vary the strictness with which we accept or reject alternative start codons. 
# The output is a .fasta file containing any number of translated protein ORFs that the user has specified based upon a minimum length also specified. 
# Note that I am not formally trained in coding with any languages and am self-taught. You should not expect this code to be written in a 'professional' manner, though you should expect it to work.

# Load packages
import re, os
from Bio import SeqIO

# Declaration of values that are required before loop start. These are reset at the end of each loop as well. Note that I know writing this script from the ground up to rely on globals is 'bad practice', this script is simply intended to be a standalone file which will read a file in and return the wanted results.
mPro = ''
altPro = ''
nonePro = ''
topHit = ''

outputText = ''
tempOutputText = ''

tempOverallList = []
tempMList = []
tempAltList = []
tempNoneList = []
codonIndex = None
noneCodonContingency = None
startCodon = re.compile(r'^.*?(M.*)')           # Regex to pull out just the sequence starting with a Methionine (or ATG)

ongoingCount = 0

### PRE-CORE LOOP SET-UP

# Locate our file of interest
print('Enter the name of the fasta-formatted file you wish to extract ORFs from. This should be a nucleotide sequence file with .fa or .fasta extension. Do not write the file extension here.')
while True:
        try:
                fileName = input()
                if os.path.isfile(fileName + '.fasta') == False and os.path.isfile(fileName + '.fa') == False:
                        raise Exception
                print('Fasta file located successfully')
                break
        except:
                print('Fasta file failed to load. If you misspelt the name, try again. If this script isn\'t in the same directory as the fasta file, move it there then try again.')
                continue
print('')

# Allow user to determine output file name
print('Enter the name which you want the output .fasta file to be called. Do not write the file extension, and make sure not to use illegal characters (i.e. \\/:?"<>|).')
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
                if os.path.isfile(outputFileName + '.fasta'):
                        print('This is already a file. Try again.')
                        continue
                break
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
                        except:
                                print('You seem to have not typed a number here. Try again.')
                
                print('What do you want the length requirement to be before we accept potentially incompletely sequenced ORFs with no start codon (i.e. at the 5\' and 3\' ends of the nucleotide sequence).')
                while True:
                        noCodonStringency = input()
                        try:
                                noCodonStringency = int(noCodonStringency)
                                break
                        except:
                                print('You seem to have not typed a number here. Try again.')
                break
        else:
                print('You didn\'t seem to type \'yes\' or \'no\', try again.')

# Load the fasta file as a generator object
if os.path.isfile(fileName + '.fasta') == True:
        records = SeqIO.parse(open(fileName + '.fasta', 'rU'), 'fasta')
else:
        records = SeqIO.parse(open(fileName + '.fa', 'rU'), 'fasta')

### CORE PROCESSING LOOP
print('Starting the core processing of this script now. Regular backups will be saved and notifications will be presented in this text area.')

# Get the nucleotide (record) out of our list of nucleotides (records) and grab the ORFs!      # Note that I didn't write the code up to the len(splitProtein[i]) < minProLen line which means it hasn't been commented in any significant degree. This was taken from the biopython user guide.
for record in records:
        for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
                for frame in range(3):
                        length = 3 * ((len(record)-frame) // 3)                                                 # Multiple of three (Not my note, whoever made this section thought it necessary to write this)
                        splitProtein = nuc[frame:frame+length].translate(table=1).split("*")                    # An asterisk represents a stop codon, thus by splitting by asterisk we obtain a list of all the ORFs in this frame
                        for i in range(len(splitProtein)):                                                      # Note that I have done a 'for i in range...' loop rather than a 'for value in splitProtein' loop which would have been simpler for a reason explained below on the 'elif i + 1 ==' line
                                if len(splitProtein[i]) < minProLen:                    # Disregard sequences that won't meet the size requirement without further processing
                                        continue
                                elif len(splitProtein) == 1:                            # If the length of splitProtein is 1, that means there was no stop codon in this frame          to make it possible to identify the last entry in splitProtein when adding asterisks (stop codons) into sequences
                                        acceptedPro = str(splitProtein[i])
                                elif i + 1 == len(splitProtein):                        # [i + 1 will equal the length of splitProtein on the final loop] Since we know there are at least 2 values in the splitProtein list, we know that the last entry will not have a stop codon. This line is the reason I used the 'for i in range' strategy rather than 'for value in splitProtein', since, depending on how short the minProLen has been specified as, it was possible doing a 'if value == splitProtein[-1]' would not accurately determine whether the ORF region we're looking at is actually the final region or is simply identical to the final region
                                        acceptedPro = str(splitProtein[i])
                                else:                                                   # Since we know there are at least 2 values in the splitProtein list, and we know that this is not the last entry, we know that this entry has a stop codon after it
                                        acceptedPro = str(splitProtein[i]) + '*'
                                # Alternative start coding      
                                protPosition = re.findall(r'(.*)' + str(splitProtein[i]) + r'(.*)', str(nuc[frame:frame+length].translate(table=1)))       # Find where the ORF region is positioned in relation to the whole nucleotide sequence
                                preProtLength = len(protPosition[0][0])*3               # Get the length of the region leading up to the protein in nucleotides (hence times 3)
                                nucSeqOfProt = nuc[frame:frame+length][preProtLength:]  # Pull out the nucleotide sequence that corresponds to the ORF region
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
                                        altPro = 'M' + altPro[1:]                       # This script assumes that the alternative start will be substituted with a methionine post-transcription
                                elif noneCodonContingency != None:                      # This will match an alternative start to 'CTG' only if 'TTG' or 'GTG' are not present
                                        altPro = acceptedPro[codonIndex:]
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
                                elif topHit == mPro:
                                        tempMList.append(topHit)                                        # These temp lists will be populated with potential ORFs from a single nucleotide sequence before being processed in the next major chunk of code starting with 'if len(tempMList + tempAltList + tempNoneList) >= 1:'
                                elif topHit == altPro:
                                        tempAltList.append(topHit)
                                elif topHit == nonePro:
                                        tempNoneList.append(topHit)

                                # Reset values required for assessing the next protein fragment
                                topHit = ''
                                mPro = ''
                                altPro = ''
                                nonePro = ''
                                codonIndex = None
                                noneCodonContingency = None

        # Sort our top hits from each inter-stop codon fragment by size and category (i.e. mPro or altPro?) and select the top X hits
        if len(tempMList + tempAltList + tempNoneList) >= 1:
                # Append '-' entries to lists which have less entries than we want to pull to allow the below 'for' loops to run without exceptions
                for i in range(0, hitsToPull-len(tempMList)):                           # hitsToPull was declared earlier in the user entry part of this script
                        tempMList.append('-')
                for i in range(0, hitsToPull-len(tempAltList)):
                        tempAltList.append('-')
                for i in range(0, hitsToPull-len(tempNoneList)):
                        tempNoneList.append('-')
                # Sort the lists by size (largest on the bottom to allow the .pop() method to remove a hit when accepted)
                tempSortedMList = sorted(tempMList, key=len)
                tempSortedAltList = sorted(tempAltList, key=len)
                tempSortedNoneList = sorted(tempNoneList, key=len)
                # Run a final size comparison to choose the best ORF(s)
                for i in range(0, hitsToPull):
                        if len(tempSortedNoneList[-1]) > len(tempSortedAltList[-1]) + noCodonStringency and len(tempSortedNoneList[-1]) > len(tempSortedMList[-1]) + noCodonStringency:         # Again, we add the stringency values to help with determining priority of ORF ordering. Since this script will often be returning either 1, 3, or 5 potential ORFs, it is important that we order these in the most logical way
                                tempOverallList.append(tempSortedNoneList[-1])
                                tempSortedNoneList.pop()
                        elif len(tempSortedAltList[-1]) > len(tempSortedMList[-1]) + altCodonStringency:
                                tempOverallList.append(tempSortedAltList[-1])
                                tempSortedAltList.pop()
                        else:
                                tempOverallList.append(tempSortedMList[-1])                                                                                                                     # By using this as the 'else' position, sequences with methionine starts will be selected in the majority of situations as the default stringency settings ensure that it is rare an alternative start is used instead of a methionine start
                                tempSortedMList.pop()
                # Format and produce the output of this script
                for i in range(0, hitsToPull):                          
                        if tempOverallList[i] == '-':                   # Because we made sure all the tempM/Alt/NoneLists had '-' added to pad out the list to have a length equal to the value of hitsToPull, when we cycle through our tempOverallList, we will often encounter '-' characters which signify the end of relevant ORFs identified in this sequence  
                                break                                   # Break out of this loop once we've fasta formatted all relevant ORF hits
                        tempOutputText += '>' + record.id + '_ORF' + str(i+1) + '\n' + tempOverallList[i] + '\n'
                outputText += tempOutputText
                
        else:
                doNothing = ''                # We don't need to do anything if no hits were found that pass the minimum length threshold

        

        # Reset temporary values required for assessing nucleotides as a whole before we run through this all again
        tempOutputText = ''
        tempOverallList = []
        tempMList = []
        tempAltList = []
        tempNoneList = []

        ongoingCount += 1       # Keep a count so we can save backups

        # Save backup if ongoingCount == 10,000. There are two sections here, the first will create the file on the first loop, the second will add to the file on subsequent loops
        if ongoingCount%10000 == 0 and os.path.isfile(os.getcwd() + '\\' + outputFileName + '.fasta') == False:
            output = open(outputFileName + '.fasta', 'w')
            output.write(outputText)
            output.close()
            print('Backup made after ' + str(ongoingCount) + ' sequences scanned for ORFs.')
            outputText = ''
        elif ongoingCount%10000 == 0 and os.path.isfile(os.getcwd() + '\\' + outputFileName + '.fasta') == True:
            output = open(outputFileName + '.fasta', 'a')
            output.write(outputText)
            output.close()
            print('Backup made after ' + str(ongoingCount) + ' sequences scanned for ORFs.')
            outputText = ''

# Dump the last few results after the script has finished, or create the output if there were less than 10,000 sequences
if os.path.isfile(os.getcwd() + '\\' + outputFileName + '.fasta') == False:
        output = open(outputFileName + '.fasta', 'w')
        output.write(outputText)
        output.close()
        print('Final save made after ' + str(ongoingCount) + ' sequences scanned for ORFs.')
elif os.path.isfile(os.getcwd() + '\\' + outputFileName + '.fasta') == True:
        output = open(outputFileName + '.fasta', 'a')
        output.write(outputText)
        output.close()
        print('Final save made after ' + str(ongoingCount) + ' sequences scanned for ORFs!')

records.close()
#### SCRIPT ALL DONE
