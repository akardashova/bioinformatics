# bioinformatics
Several projects and works from the bioinformatics course   
* global_aligners: <br>Programs calculate the optimal global alignment of two sequences using the Needleman-Wunsch algorithm, using a DNAfull penalty matrix and an affine or linear gap penalty.<br>
* mapper: <br>The mapper based on the construction of a hash table for the reference genome. The program accepts a genome in FASTA format, the number k, creates a hash table and saves it to a file. After that, it reads this file and reads in FASTQ format and performs mapping of the reads (without gaps).<br>
*(a file with a hash table weighs a lot, so a part of it is attached as an example)*<br>
* human_dna_project: <br>The program determines at what minimum number of k there is a nucleotide sequence of length k that does not occur in the human genome and counts the number of such sequences.
