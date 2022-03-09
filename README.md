# motif-mark  
  
This is an object-oriented Python script that creates a .png image denoting binding sites for up to 5 motifs (of less than or equal to 10 bases each) on up to 10 sequences (of less than or equal to 1000 bases each). The image distinguishes between introns and exons and shows motif, sequence, intron, and exon sizes to scale. The prefix of the resulting .png file will be the same prefix used for the FASTA file.  
  
Input requirements:  
 - A FASTA file containing less than or equal to 10 sequences of less than or equal to 1000 bases each. Introns should be denoted by lower case nucleotides and exons should be denoted by upper case nucleotides.  
 - A txt file containing one motif per line. Motifs may be denoted in upper or lower case, and may contain any IUPAC degenerate base symbol.  

