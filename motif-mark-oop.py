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
    'u' : 't',
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
    #methods
    def get_exon_coords(self):
        m = re.search(r'[ATCG]+',self.seq)
        return [m.start(),m.end()]

class Motif:
    def __init__(self,seq,color):
        #data
        self.seq = seq.lower()
        self.color = color
    #methods
    def change_color(self,Context):
        Context.set_source_rgba(self.color[0],self.color[1],self.color[2],0.5)
    def find_motifs(self,Sequence,Context,base_coords,sequence_count):
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
        match = 'not_matched_yet'
        self.change_color(Context)
        while True:
            if match == 'not_matched_yet':
                match = p.search(Sequence.seq.lower())
            else:
                match = p.search(Sequence.seq.lower(), pos=match.start() + 1)
            if match == None:
                break
            Context.move_to(base_coords[0] + match.start(), base_coords[1] + (100 * sequence_count))
            Context.line_to(base_coords[0] + match.start() + len(self.seq), base_coords[1] + (100* sequence_count))
            Context.stroke()

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
width = 1100 
height = (100 * len(seq_obs.keys())) + 300
surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
context = cairo.Context(surface)

#draw sequence introns and exons
count = 0
base_coords = [50,100]
context.set_font_size(15)
for seq_ob in seq_obs:
    exon_coords = seq_obs[seq_ob].get_exon_coords()
    context.move_to(base_coords[0], base_coords[1]+(100*count)-40)
    context.show_text(seq_obs[seq_ob].header)
    context.move_to(base_coords[0], base_coords[1]+(100*count))
    context.set_line_width(10)
    context.line_to(base_coords[0]+len(seq_obs[seq_ob].seq), base_coords[1]+(100*count))
    context.stroke()
    context.set_line_width(40)
    context.move_to(base_coords[0]+exon_coords[0], base_coords[1]+(100*count))
    context.line_to(base_coords[0]+exon_coords[1], base_coords[1]+(100*count))
    context.stroke()
    count += 1

# main loop: for each Sequence object in seq_obs:
# - label sequence and draw introns and exons with pycairo
# - for each sequence, find and draw all occurrences of each motif on sequence

#dictionary containing rgb values for up to 5 different motifs. 
#keys are counts (corresponding to count below while looping through motifs file)
#values are lists of rgb values
#this dictionary is used to assign a color to each motif object.
color_dict = {
    1 : [0.2, 0.23, 0.9],
    2 : [0.9, 0.1, 0.1],
    3 : [0.4, 0.9, 0.4],
    4 : [0.9, 0.5, 0.1],
    5 : [0.5, 0.1, 0.5]
}

# test drawing lines
#with open(args.motif, 'r') as fh:
    #count = 0
    #coords = [50,50]
    #for line in fh:
        #count += 1
        #line = line.strip()
        #motif_ob = Motif(line,color_dict[count])
        #motif_ob.change_color(context)
        #context.move_to(coords[0],coords[1])
        #context.line_to(coords[0]+50, coords[1])
        #context.stroke()
        #coords[1] = coords[1] + 50
    #surface.write_to_png('test.png')
        

fh = open(args.motif, 'r')
color_count = 0
for line in fh:
    color_count += 1
    line = line.strip()
    motif_ob = Motif(line,color_dict[color_count])
    count = 0
    context.set_line_width(50)
    for seq_ob in seq_obs:
        motif_ob.find_motifs(seq_obs[seq_ob],context,base_coords,count)
        count += 1
    #add motif to legend
    context.set_line_width(30)
    context.move_to(base_coords[0], base_coords[1] + 100*len(seq_obs) + 40*color_count + 10)
    context.line_to(base_coords[0]+30, base_coords[1] + 100*len(seq_obs) + 40*color_count + 10)
    context.stroke()
    context.set_font_size(15)
    context.set_source_rgb(0,0,0)
    context.move_to(base_coords[0]+50, base_coords[1] + 100*len(seq_obs) + 40*color_count + 15)
    context.show_text(motif_ob.seq)
    # draw introns and exons
    # parse motif file
        # for each motif, draw all instances of motif
fh.close()

#add intron to legend
context.set_source_rgb(0,0,0)
context.set_line_width(10)
context.move_to(base_coords[0], base_coords[1] + 100*len(seq_obs) - 40)
context.line_to(base_coords[0]+30, base_coords[1] + 100*len(seq_obs) - 40)
context.stroke()
context.move_to(base_coords[0]+50, base_coords[1] + 100*len(seq_obs) - 35)
context.show_text('Intron')
#add exon to legend
context.set_line_width(50)
context.move_to(base_coords[0], base_coords[1] + 100*len(seq_obs))
context.line_to(base_coords[0]+30, base_coords[1] + 100*len(seq_obs))
context.stroke()
context.move_to(base_coords[0]+50, base_coords[1] + 100*len(seq_obs) + 5)
context.show_text('Exon')

out_name = re.sub(r'(.+)\.fa(sta)*',r'\1',args.fasta)
surface.write_to_png(out_name+'.png')

# make legend
# deal with reverse complements?
# why is the motif thing longer than the sequence