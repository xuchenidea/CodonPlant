README
One program works for codon optimazation. This program is writen by the standard c language(c99). When you give a codon sequence, it can make an optimazation codon sequence. It use the Genetic Algorithm to get this goal.

How to use
codon [seqence file] [codon file] [codon optimazation file] [codon special] [population number] [cut factor length] [iteration num] [mution ration]
There are eight parameters.
[seqence file] 
the codon seqence file which needs to be optimazation
[codon file] 
Synonymous codon file. For example,
	{ "TTT", "TTC" },
	{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" },
	{ "ATT", "ATC", "ATA" },
	{ "ATG" },
	{ "GTT", "GTC", "GTA", "GTG" },
	{ "AGT", "AGC", "TCT", "TCC", "TCA", "TCG" },
	{ "CCT", "CCC", "CCA", "CCG" },
	{ "ACT", "ACC", "ACA", "ACG" },
	{ "GCT", "GCC", "GCA", "GCG" },
......
As you see, the same Synonymous codons should be on one line.
[codon optimazation file]
 the corresponding codon numbers as the order of [codon file], for example:
	{ { 4844 ,1.04 },{ 4508 ,0.96 } },  //{ "TTT", "TTC" },
	//{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" },
	{ { 1532, 0.47 },{ 4816, 1.48 },{ 5533, 1.70 },{ 3422, 1.05 },{ 1551, 0.48 },{ 2646 ,0.81 } },
The first is the codon numbers of the whole codon database. The second is the RSCU of this codon.
[codon special] 
These codons shoulb not be visibale in the final sequence. For example:
"AATAAA", "ATAAT", "AATTAA", "AACCAA"
Our goal is to vanish them.

the following 4 paramaters are our Genetic Algorithm needed.
[population number] 
The population number
[cut factor length]
We don't want the final sequence has the duplication string. So we fomulate the cut factor length. It should be 7-12.
[iteration number]
the iteration number of our algorithom. 
[mution ration]
In the each step, the mution ration of the whole sequence. It should be a percentage number.

Make a test
make sure these files in the same directory:
codon seq.txt codon.txt opt_codon.txt  special_codon.txt
then type:
codon seq.txt codon.txt opt_codon.txt  special_codon.txt 20 10 100 5
