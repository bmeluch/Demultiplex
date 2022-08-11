#!/usr/bin/env python
# Author: Bea Meluch

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
Updated 02 Aug 2022 for Demux assignment.'''

__version__ = "0.6"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = set('ATGCNatgcn')
RNA_bases = set('AUGCNaugcn')
DNA_pairing: dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
RNA_pairing: dict = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

def convert_phred(letter: str) -> int:
    """Converts a single character into a Phred+33 score"""
    # make sure input is a single character phred score
    if len(letter) != 1:
        raise ValueError("Input must be a single character.")
    return(ord(letter)-33)

def qual_score(phred_score: str) -> float:
    """Takes a Phred score string and returns the average quality score of the string as a float."""
    # initialize running total
    total_score: int = 0
    # loop through each character of string
    for ch in range(len(phred_score)):
        # add integer score output of convert_phred on each character to the running total
        total_score += convert_phred(phred_score[ch])
    # calculate the average
    return total_score/len(phred_score)

def qual_pass(seq: str, cutoff: int) -> bool:
    """Takes a Phred score string and returns True if all bases pass a quality threshold, False if any single base does not pass"""
    for base in seq:
        if convert_phred(base) < cutoff:
            return False
    return True

def validate_base_seq(seq: str,RNAflag=False) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNA_bases if RNAflag else DNA_bases)

def gc_content(seq: str) -> float:
    '''Given a string, returns GC content as a float'''
    seq = seq.upper()
    if not validate_base_seq(seq):
        raise ValueError("Not a valid nucleic acid sequence.")
    gc: int = 0
    for ch in seq:
        if ch == "G" or ch == "C":
            gc += 1
    return gc/len(seq.strip())

def oneline_fasta(input_file: str, output_file: str):
    '''Given a FASTA file, creates a new FASTA file with newlines
    removed from sequences'''

    lcount: int = 1
    header: str = ""
    seq: str = ""

    # open input and output files
    # iterate over input file lines
    with open(input_file, 'r') as in_fa, open(output_file, 'w') as out_fa:
        for line in in_fa:
            # header line test
            if line[0]==">":
                # SAVE EVERYTHING FROM LAST TIME
                # clean up
                seq = seq.replace("\n", "")

                # write to file IF lines are not blank (to avoid writing blank lines at the beginning)
                if header and seq:
                    out_fa.write(header+'\n')
                    out_fa.write(seq+'\n')
                
                # collect next record
                # get current header
                header = line.strip()
                # blank out sequence again
                seq = ""
            # not a header line: add sequence to sequence string
            else:
                seq += line
        # finish writing last record
        seq = seq.replace("\n", "")
        out_fa.write(header+'\n')
        out_fa.write(seq+'\n')

    return(output_file+" completed")

def reverse_complement(seq: str, RNAflag=False) -> str:
    '''Takes a string, if the string is a valid base sequence, returns the reverse complement. DNA or RNA, case insensitive.'''
    seq = seq.upper()
    rev_comp: str = ""
    if not validate_base_seq(seq):
        raise ValueError("Input is not a valid base sequence.")
    for base in seq:
        rev_comp += DNA_pairing[base]
    rev_comp = rev_comp[::-1]
    return rev_comp

if __name__ == "__main__":
    # write tests for functions above

    # autograder tests, do not change
    # Check that convert_phred returns the correct value for several different inputs"""
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")

    #autograder tests, do not change
    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    print("Passed DNA and RNA tests")

    #autograder tests, do not change
    #initializing variables, be sure to run this cell and do not change
    phred_score: str = "FFHHHHHJJJJIJIJJJIJJJJJJIIIJJJEHJJJJJJJIJIDGEHIJJFIGGGHFGHGFFF@EEDE@C??DDDDDDD@CDDDDBBDDDBDBDD@"
    assert qual_score(phred_score) == 37.62105263157895, "wrong average phred score"
    print("You calcluated the correct average phred score")

    #test gc_content
    assert gc_content("ATGATCACGACTGCTACGACTACGACTACG") == 0.5, "wrong GC content"
    print("Correct GC content")

    # test reverse complement
    assert reverse_complement("GATTACA") == "TGTAATC", "wrong rev comp"
    assert reverse_complement("augc") == "GCAU", "wrong rev comp"
    print("Correct reverse complementing")