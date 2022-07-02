The mapper based on the construction of a hash table for the reference genome. 
The program accepts a genome in FASTA format, the number k, creates a hash table and saves it to a file. 
After that, it reads this file and reads in FASTQ format and performs mapping of the reads (without gaps).
(a file with a hash table weighs a lot, so a part of it is attached as an example)