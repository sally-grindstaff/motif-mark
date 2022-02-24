#!/usr/bin/env python

import argparse
import cairo
import re

def get_args():
    parser = argparse.ArgumentParser(description='Motif mark: Takes sequence and motif information and returns a png image that displays motif binding sites, to scale, in each exon and its flanking introns. Can handle multiple sequences and multiple motifs.')
    parser.add_argument('-f', '--fasta', help='Path to fasta file containing DNA sequence information. One header per sequence. Introns should be depicted in lower case nucleotide letters and exons should be depicted in upper case nucleotide letters. Note that the output file will have the same prefix as the input fasta file.', required=True, type=str)
    parser.add_argument('-m', '--motif', help='Path to txt file containing one motif per line. Motifs may contain any IUPAC degenerate base symbol.')
    return parser.parse_args()

args = get_args()

symbol_dict = {
    'a' : 'a',
    'c': 'c',
    'g' : 'g',
    't' : 't',
    'w' : ['a','t'],
    's' : ['c','g'],
    'm' : ['a','c'],
    'k' : ['g','t'],
    'r' : ['a','g'],
    'y' : ['c','t'],
    'b' : ['c','g','t'],
    'd' : ['a','g','t'],
    'h' : ['a','c','t'],
    'v' : ['a','c','g'],
    'n' : ['a','c','g','t']
}
# keys are IUPAC degenerate base symbols, values are represented bases

class Sequence:
    def __init__(self,seq,header):
        self.seq = seq
        self.header = header
        #self.exon_coords = get_exon_coords()
    #methods
    def get_exon_coords(self):
        pass #reference self.seq and use upper/lower notation to return tuple of start and end coordinates

class Motif:
    def __init__(self,seq,color):
        #data
        self.seq = seq.lower()
        self.color = color
    #methods
    def find_motifs(self,Sequence,Context):
        '''Interacts with Sequence object and Pycairo Context object; finds instances of motif in sequence and draws them'''
        # Compile regex search pattern from motif
        pattern = ''
        for i in self.seq:
            pattern = pattern + '('
            count = 0
            for j in symbol_dict[i]:
                if count > 0:
                    pattern = pattern + '|' # if base symbol is degenerate, use 'or' symbol between the corresponding bases
                pattern = pattern + j
                count += 1
            pattern = pattern + ')'
        p = re.compile(pattern)
        # search for all instances of motif pattern in Sequence and draw each one on Pycairo context object
            # note, using re.search in a loop rather than re.findall in order to allow for overlapping instances of motif
        while True:
            match = p.search(Sequence.seq)
            if match == None:
                break
            return match.group(0)

test = Motif('yaaan','red')

testseq = Sequence('caaaa','seq1')

#create sequence objects and something to keep track of them by parsing fasta file
def fasta_parser(fasta: str):
    '''Reads fasta file, creates object for each sequence'''
    seq = ''
    seq_obs = {}
    with open(fasta, 'r') as fh:
        for line in fh:
            line = line.strip()
            if line[0] == '>' and seq != '':
                #create sequence object from header and seq
                seq_obs[header] = Sequence(seq,header) 
                seq = ''
            if line[0] == '>':
                header = line[1:]
            else:
                seq += line
        #create final sequence object from header and seq
        seq_obs[header] = Sequence(seq,header)
    return seq_obs

seq_obs = fasta_parser(args.fasta)

#define cairo surface and context
width = 500 
height = 100 * len(keys(seq_obs))
surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
context = cairo.Context(surface)
context.set_line_width(1)
