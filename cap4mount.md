C H A P T E R <u><mark>4</mark></u> 

# **Multiple Sequence Alignment** 

##### **INTRODUCTION, 140** 

**Genome sequencing, 142 Uses of multiple sequence alignments, 142 Relationship of multiple sequence alignment to phylogenetic analysis, 143 METHODS, 144** 

**Multiple sequence alignment as an extension of sequence pair alignment by dynamic programming, 145 Scoring multiple sequence alignments, 151 Progressive methods of multiple sequence alignment, 152** CLUSTALW, 153 PILEUP, 155 Problems with progressive alignment, 155 **Iterative methods of multiple sequence alignment, 157** Genetic algorithm, 157 Hidden Markov models of multiple sequence alignment, 160 **Other programs and methods for multiple sequence alignment, 160 Localized alignments in sequences, 161** Profile analysis, 161 Block analysis, 165 Extraction of blocks from a global or local multiple sequence alignment, 165 Pattern searching, 170 Blocks produced by the BLOCKS server from unaligned sequences, 171 The eMOTIF method of motif analysis, 173 

**Statistical methods for aiding alignment, 173** Expectation maximization algorithm, 173 Multiple EM for motif elicitation (MEME), 177 The Gibbs sampler, 177 Hidden Markov models, 185 Motif-based hidden Markov models, 190 **Position-specific scoring matrices, 192** Information content of the PSSM, 195 Sequence logos, 196 **Multiple sequence alignment editors and formatters, 198** Sequence editors, 199 Sequence formatters, 200 **REFERENCES, 200** 

**139** 

**140** I CHAPTER 4 

## **INTRODUCTION** 

O NE OF THE MOST IMPORTANT CONTRIBUTIONS of molecular biology to evolutionary analysis is the discovery that the DNA sequences of different organisms are often related. Similar genes are conserved across widely divergent species, often performing a similar or even identical function, and at other times, mutating or rearranging to perform an altered function through the forces of natural selection. Thus, many genes are represented in highly conserved forms in organisms. Through simultaneous alignment of the sequences of these genes, sequence patterns that have been subject to alteration may be analyzed. 

Because the potential for learning about the structure and function of molecules by multiple sequence alignment (msa) is so great, computational methods have received a great deal of attention. In msa, sequences are aligned optimally by bringing the greatest number of similar characters into register in the same column of the alignment, just as described in Chapter 3 for the alignment of two sequences. Computationally, msa presents several difficult challenges. First, finding an optimal alignment of more than two sequences that includes matches, mismatches, and gaps, and that takes into account the degree of variation in all of the sequences at the same time poses a very difficult challenge. The dynamic programming algorithm used for optimal alignment of pairs of sequences can be extended to three sequences, but for more than three sequences, only a small number of relatively short sequences may be analyzed. Thus, approximate methods are used, including (1) a progressive global alignment of the sequences starting with an alignment of the most alike sequences and then building an alignment by adding more sequences, (2) iterative methods that make an initial alignment of groups of sequences and then revise the alignment to achieve a more reasonable result, (3) alignments based on locally conserved patterns found in the same order in the sequences, and (4) use of statistical methods and probabilistic models of the sequences. A second computational challenge is identifying a reasonable method of obtaining a cumulative score for the substitutions in the column of an msa. Finally, the placement and scoring of gaps in the various sequences of an msa presents an additional challenge. 

The msa of a set of sequences may also be viewed as an evolutionary history of the sequences. If the sequences in the msa align very well, they are likely to be recently derived from a common ancestor sequence. Conversely, a group of poorly aligned sequences share a more complex and distant evolutionary relationship. The task of aligning a set of sequences, some more closely and others less closely related, is identical to that of discovering the evolutionary relationships among the sequences. 

As with aligning a pair of sequences, the difficulty in aligning a group of sequences varies considerably with sequence similarity. On the one hand, if the amount of sequence variation is minimal, it is quite straightforward to align the sequences, even without the assistance of a computer program. On the other hand, if the amount of sequence variation is great, it may be very difficult to find an optimal alignment of the sequences because so many combinations of substitutions, insertions, and deletions, each predicting a different alignment, are possible. 

The availability of a subset of the many multiple sequence alignment programs is shown in Table 4.1. A flowchart illustrating the considerations to be made in choosing an alignment method is shown on page 144. 

When dealing with a sequence of unknown function, the presence of similar domains in several similar sequences implies a similar biochemical function or structural fold that may become the basis of further experimental investigation. A group of similar sequences may 

MULTIPLE SEQUENCE ALIGNMENT I **141** 

**Table 4.1.** _Web sites and program sources for multiple sequence alignment_ 

|**Name**|**Source**|**Reference**|
|---|---|---|
|**_Global alignments including progressive_**|<br>||
|CLUSTALW or CLUSTALX (latter has<br>graphical interface)|FTP to ftp.ebi.ac.uk/pub/software<sup>a,d</sup>|Thompson et al. (1994a, 1997); Higgins<br>et al. (1996)|
|MSA|http://www.psc.edu/<sup>b</sup><br>http://www.ibc.wustl.edu/ibc/msa.html<sup>c</sup><br>FTP to fastlink.nih.gov/pub/msa|Lipman et al. (1989);<br>Gupta et al. (1995)|
|PRALINE|http://mathbio.nimr.mrc.ac.uk/~jhering/<br>praline|Heringa (1999)|
|**_Iterative and other methods_**|||
|DIALIGN segment alignment|http://www.gsf.de/biodv/dialign.html|Morgenstern et al. (1996)|
|MultAlin|http://protein.toulouse.inra.fr/multalin.<br>html|Corpet (1988)|
|PRRP progressive global alignment<br>(randomly or doubly nested)|ftp.genome.ad.jp/pub/genome/saitama-<br>cc|Gotoh (1996)|
|SAGA genetic algorithm|http://igs-server.cnrs-mrs.fr/�<br>cnotred/Projects_home_page/saga_<br>home_page.html|Notredame and Higgins (1996)|
|**_Local alignments of proteins_**|||
|<br>Aligned Segment Statistical Evaluation<br>Tool (Asset)|FTP to ncbi.nlm.nih.gov/pub/neuwald/<br>asset|Neuwald and Green (1994)|
|BLOCKS Web site|http://blocks.fhcrc.org/blocks/|Henikoff and Henikoff (1991, 1992)|
|eMOTIF Web server|http://dna.Stanford.EDU/emotif/|Nevill-Manning et al. (1998)|
|GIBBS, the Gibbs sampler statistical<br>method|FTP to ncbi.nlm.nih.gov/pub/neuwald/<br>gibbs9_95/|Lawrence et al. (1993); Liu et al. (1995);<br>Neuwald et al. (1995)|
|HMMER hidden Markov model software|http://hmmer.wustl.edu/|Eddy (1998)|
|MACAW, a workbench for multiple<br>alignment construction and analysis|FTP to ncbi.nlm.nih.gov/pub/macaw|Schuler et al. (1991)|
|<br>MEME Web site, expectation<br>maximization method|http://meme.sdsc.edu/meme/website/|Bailey and Elkan (1995);<br>Grundy et al. (1996, 1997); Bailey<br>and Gribskov (1998)|
|Profile analysis at UCSD<sup>a,e</sup>|http://www.sdsc.edu/projects/profile/|<br>Gribskov and Veretnik (1996)|
|SAM hidden Markov model Web site|http://www.cse.ucsc.edu/research/comp<br>bio/sam.html|Krogh et al. (1994); Hughey and Krogh<br>(1996)|



> a Lists of additional Web sites for msa are maintained at: http://www.ebi.ac.uk/biocat/, http://www.hgmp.mrc.ac.uk/Registered/Menu/prot-mult.html, http://www.hum-molgen.de/BioLinks/Biocomp.html, http://biocenter.helsinki.fi/bi/rnd/biocomp/. Reviews on the performance of msa software are given in McClure et al. (1994; progressive alignment methods), Gotoh (1996) and Thompson et al. (1999), a review of Web sites is given in Briffeuil et al. (1998) and a review on iterative algorithms is given in Hirosawa et al. (1995) and Gotoh (1999). The performance of msa programs is commonly assessed by comparing the computed msa with a structural alignment of the proteins and by other objective methods (Notredame et al. 1998). Many of these programs are computationally complex and must be set up on a local site. 

b The Biomedical Supercomputing facility at the University of Pittsburgh Supercomputing Facility provides accounts (see http://www.psc.edu/biomed/seqanal/grants.html) that provide access to several different versions of MSA and profile analysis. MSA 50 150 will align no more than 50 sequences each less than 150 residues long, MSA 25 500 will align no more than 25 sequences each less than 200 residues long, and MSA10 1000 will align no more than 10 sequences each less than 1000 long. 

c The MSA server at the University of Washington will take up to 8 sequences, each less than 500 long. d CLUSTALW is also available as freeware that runs on PCs and Macintosh computers from the same FTP site. e Profile generating programs are available by FTP from ftp.sdsc.edu/pub/sdsc/biology and are included in the Genetics Computer Group suite of programs (http://www.gcg.com/), although the most recent features of Gribskov and Veretnik (1996) are not included. 

define a protein family that may share a common biochemical function or evolutionary origin. Similar proteins have been organized into databases of protein families that are described in Chapter 9. 

**142** I CHAPTER 4 

### **<u>GENOME SEQUENCING</u>** 

One application of multiple sequence alignment algorithms is in genome sequencing projects discussed in Chapter 2. Instead of cloning and arranging a very large number of fragments of a large DNA molecule, and then moving along the molecule and sequencing the fragments in order, random fragments of the large molecule are sequenced, and those that overlap are found by a msa program. This approach enables automated assembly of large sequences. Bacterial genomes have been quite readily sequenced by this method, and it has also been used to assemble portions of the _Drosophila_ and human genomes at Celera Genomics (Weber and Myers 1997 and see Chapter 10). 

The requirements for a msa program for genome projects differ in several respects from those for general sequence analysis. First, the sequences are fragments of the same large sequence molecule, and the sequences of overlapping fragments should be the same except for sequence copying and reading errors, which may introduce the equivalent of substitutions and insertions/deletions between the compared fragments. Thus, there should be one correct alignment that corresponds to that of the genome sequence instead of a range of possibilities. Second, the sequences may be from one DNA strand or the other and hence the complements of each sequence must also be compared. Third, sequence fragments will usually overlap, but by an unknown amount, and, in some cases, one sequence may be included within another. Finally, all of the overlapping pairs of sequence fragments must be assembled into a large, composite genome sequence, taking into account any redundant or inconsistent information. Interested readers may wish to consult a description of the type of methodology (Myers 1995 and see Chapter 10) and a comparison of the methods, including several commercial packages that are useful for managing the sequence data from laboratory sequencing projects (Miller and Powell 1994). The Institutue of Genome Research (http://www.tigr.org/) has also developed and made available software and methods for genome assembly and analysis. 

### **<u>USES OF MULTIPLE SEQUENCE ALIGNMENTS</u>** 

Just as the alignment of a pair of nucleic acid or protein sequences can reveal whether or not there is an evolutionary relationship between the sequences, so can the alignment of three or more sequences reveal relationships among multiple sequences. Multiple sequence alignment of a set of sequences can provide information as to the most alike regions in the set. In proteins, such regions may represent conserved functional or structural domains. 

If the structure of one or more members of the alignment is known, it may be possible to predict which amino acids occupy the same spatial relationship in other proteins in the alignment. In nucleic acids, such alignments also reveal structural and functional relationships. For example, aligned promoters of a set of similarly regulated genes may reveal consensus binding sites for regulatory proteins. Methods for finding such sites in nucleic acid sequences are discussed in Chapter 8. 

Another use for consensus information retrieved from a multiple sequence alignment is for the prediction of specific probes for other members of the same group or family of similar sequences in the same or other organisms. There are both computer and molecular biology applications. Once a consensus pattern has been found, database searching pro- 

MULTIPLE SEQUENCE ALIGNMENT I **143** 

grams (Chapter 7) may be used to find other sequences with a similar pattern. In the laboratory, a reasonable consensus of such patterns may be used to design polymerase chain reaction (PCR) primers for amplification of related sequences. 

### **~~RELATIONSHIP OF MULTIPLE SEQUENCE ALIGNMENT TO PHYLOGENETIC ANALYSIS~~** 

Once the msa has been found, the number or types of changes in the aligned sequence residues may be used for a phylogenetic analysis. The alignment provides a prediction as to which sequence characters correspond. Each column in the alignment predicts the mutations that occurred at one site during the evolution of the sequence family, as illustrated in Figure 4.1. Within the column are original characters that were present early, as well as other derived characters that appeared later in evolutionary time. In some cases, the position is so important for function that mutational changes are not observed. It is these conserved positions that are useful for producing an alignment. In other cases, the position is less important, and substitutions are observed. Deletions and insertions may also be present in some regions of the alignment. Thus, starting with the alignment, one can hope to dissect the order of appearance of the sequences during evolution. 



<!-- Start of picture text -->
seqA N • F L S<br>seqB N • F – S<br>seqC N K Y L S<br>seqD N • Y L S<br>N Y L S N K Y L S N F S N F L S<br>+K –L<br>Y to F<br><!-- End of picture text -->

**Figure 4.1.** The close relationship between msa and evolutionary tree construction. Shown is a short section of one msa of four protein sequences including conserved and substituted positions, an insertion (of K) and a deletion (of L). Below is a hypothetical evolutionary tree that could have generated these sequence changes. Each outer “branch” in the tree represents one of the sequences. The outer branches are also referred to as “leaves.” The deepest, oldest branch is that of sequence D, followed by A, then by B and C. The optimal alignment of several sequences can thereby be thought of as minimizing the number of mutational steps in an evolutionary tree for which the sequences are the outer branches or leaves. The mathematical solution to this problem was first outlined by Sankoff (1975). Fast multiple sequence alignment programs that are tree-based have since been developed (Ravi and Kececioglu 1998). However, such an approach depends on knowing the evolutionary tree to perform an alignment, and often this is not the case. Usually, pair-wise alignments are generated first and then used to predict the tree. In this example, the alignment could be explained by several different trees, including the one shown, following one of several types of analyses described in Chapter 6. The sequences then become the outer leaves of the tree, and the inner branches are constructed by this analysis. 

**144** I CHAPTER 4 

## **METHODS** 



<!-- Start of picture text -->
Choose Are the Yes Perform Is a<br>three or sequences global convincing<br>more protein alignment of alignment<br>sequences. 1 sequences? sequences. 2 produced? 4 Yes<br>No<br>Are the<br>Yes Translate Are there a<br>sequences into protein large<br>cDNA<br>sequences. number of<br>sequences? 3<br>sequences? 5<br>No<br>No<br>No<br>Are the sequences Yes Predict<br>genomic sequences Make a<br>that encode related gene profile<br>structure.<br>proteins? or PSSM<br>representation<br>of the<br>No<br>alignment.<br>Analyze for No Do the Analyze promoter<br>patterns, sequences regions, intron-exon Produce a<br>repeats, etc., encode RNA boundaries, etc., as hidden<br>as described molecules? described in Chapter 8. Markov<br>Yes<br>in Chapters model.<br>2 and 10.<br>Yes<br>Analyze for secondary Search for<br>structure as described blocks. 6<br>in Chapter 5.<br><!-- End of picture text -->

1. The sequence chosen for analysis may already be known to be similar on the basis of pair-wise alignments (Chapter 2), but sequences related by other criteria may also be used. Complex features of the sequences, including repeated or low-complexity regions that interfere with alignments, can be analyzed as described in Chapters 2 and 7. The flowchart describes the production of four classes of multiple sequence alignment. 

   - a. A global alignment includes the entire range of each sequence in the alignment, and is usually produced by extensions to the dynamic programming global alignment algorithm that is used for aligning pairs of sequences, but other methods are also used. 

   - b. A sequence block is an alignment of common patterns in protein sequences that includes matches and mismatches in each column found by using pattern-finding algorithms, but no gaps (insertions and deletions) are present. 

   - c. An alignment of common patterns in protein sequences that includes matches, mismatches, insertions, and deletions may be used to make a type of scoring matrix called a profile. 

   - d. A hidden Markov model is a probabilistic model of a global alignment of protein sequences or of a conserved local region (similar to a sequence profile) in those sequences that includes matches, mismatches, insertions, and deletions. The model is “trained” to represent the set of sequences. 

MULTIPLE SEQUENCE ALIGNMENT I **145** 

Methods for finding common patterns in DNA sequences are discussed in Chapter 8. 

2. Examples of global alignment, as well as other programs from which to choose, are given in the global alignments and iterative and other methods sections of Table 4.1. 

3. cDNA sequences of the same gene from a group of organisms may be multiply aligned by a global method so that synonymous (i.e., change the amino acid) and nonsynonymous (i.e., do not change the amino acid) sequences may be analyzed, as described in Chapter 6 (see also note 2). 

4. A convincing alignment should include a series of columns in which a majority of the sequences have the same amino acid or an amino acid that is a conservative substitution for that amino acid, with relatively few examples of other substitutions or gaps in these columns. These columns of alike amino acids should be found throughout the alignment, often clustered into domains. There may also be variable regions in the alignment that represent sequences that diverged more during the evolution of the protein family. 

5. This decision rests on whether or not there are enough sequences on which to build a hidden Markov model of the entire alignment or of a well-defined region in the alignment (a profile hidden Markov model). For sequences that are related but show considerable variations in many columns, as many as 100 sequences may be needed to produce a hidden Markov model of the alignment. This number is reduced to approximately 25–50 if there is less variation among the sequences. A scoring matrix representing the sequence variation found in each column of the alignment may also be made. These matrices may accommodate gaps in the alignment (a profile or HMM profile) or may not include gaps (position-specific scoring matrix). 

6. For finding patterns common to the sequences, pattern-searching algorithms and statistical methods are used. The former search for a set of matched sequence characters that are present in the sequences. The latter perform an exhaustive analysis of sequence “windows” in the sequences to find the most alike amino acid patterns by the expectation maximization (EM) or Gibbs sampling algorithms. These methods are described in the text. 

### **MULTIPLE SEQUENCE ALIGNMENT AS AN EXTENSION OF SEQUENCE PAIR** **~~ALIGNMENT BY DYNAMIC PROGRAMMING~~** 

The dynamic programming algorithm described in Chapter 2 provides an optimal alignment of two sequences. In the program MSA (Lipman et al. 1989), application of the global alignment algorithm has been extended to provide an optimal alignment of a small number of sequences greater than two. Gupta et al. (1995) have shown, however, that MSA rarely produces a provable optimal alignment. The number of sequences that can be aligned is limited because the number of computational steps and the amount of memory required grow exponentially with the number of sequences to be analyzed. This limitation means that the program has somewhat limited application to a small number of sequences. 

Recall that the dynamic programming method of sequence alignment between two sequences builds a scoring matrix where each position provides the best alignment up to that point in the sequence comparison. The number of comparisons that must be made to fill this matrix without using any short cuts and excluding gaps is the product of the length of the two sequences. Imagine extending this analysis to three or more sequences. For three sequences, instead of the two-dimensional matrix for two sequences, think of the lattice of a cube that is to be filled with calculated dynamic programming scores. Scoring positions on three surfaces of the cube will represent the alignment values between a pair of the sequences, ignoring the third sequence, as illustrated in Figure 4.2. In MSA, positions inside the lattice of the cube are given values based on the sum of the initial scores of the three pairs of sequences. 

**146** I CHAPTER 4 

For three protein sequences each 300 amino acids in length and excluding gaps, the number of comparisons to be made by dynamic programming is equal to 300<sup>3</sup> � 2.7 � 10<sup>7</sup> , whereas only 300<sup>2</sup> � 9 � 10<sup>4</sup> is required for two sequences of this length. This number is sufficiently small that alignment of three sequences by this method is practical. For alignment of more than three sequences, one has to imagine filling an _N_ -dimensional space or hypercube. The number of steps and memory required for a 300-amino-acid sequence (300<sup>_N_</sup> , where _N_ is the number of sequences) then becomes too large for most practical purposes, and it is necessary to find a way to reduce the number of comparisons that must be made without compromising the attempt to find an optimal alignment. Fortunately, Carrillo and Lipman (1988) found such a method, called the sum of pairs, or SP method. Since the publication of the MSA program, Gupta et al. (1995) have substantially reduced the memory requirements and number of steps required. The enhanced version of MSA is available by anonymous FTP from fastlink.nih.gov/pub/msa. 

The basic idea is that a multiple sequence alignment imposes an alignment on each of the pairs of sequences. The heavy arrow in Figure 4.2 represents the path followed in the cube to find a msa for three sequences, but the msa can be projected on to the sides of the cube, thus defining an alignment for each pair of sequences. The alignments found for each pair of sequences likewise impose bounds on the location of the msa within the cube, and thus defines the number of positions within the cube that have to be evaluated. Pair-wise alignments are first computed between each pair of sequences. Next, a trial msa is produced by first predicting a phylogenetic tree for the sequences (Saitou and Nei 1987; see Chapter 6 for the neighbor-joining method of tree construction), and the sequences are 



<!-- Start of picture text -->
sequence A<br>A-C<br>A-B<br>sequence C<br>B-C<br>sequence B<br><!-- End of picture text -->

**Figure 4.2.** Alignment of three sequences by dynamic programming. Arrows on the surfaces of the cube indicate the direction for filling in the scoring matrix for pairs of sequences, A with B, etc., performed as previously described. The alignment of all three sequences requires filling in the lattice of the cube space with optimal alignment scores following the same algorithm. The best score at each interior position requires a consideration of all possible moves within the cube up to that point in the alignment. The trace-back matrix will align positions in all three sequences including gaps. 

MULTIPLE SEQUENCE ALIGNMENT I **147** 

then multiply aligned in the order of their relationship on the tree. This method is used by other programs described below (e.g., PILEUP, CLUSTALW) and provides a heuristic alignment that is not guaranteed to be optimal. However, the alignment serves to provide a limit to the space within the cube within which optimal alignments are likely to be found. In Figure 4.3, the green area on the left surface of the cube is bounded by the optimal alignment of sequences B and C and a projection of the heuristic alignment for all three sequences. The orange and blue areas are similarly defined for other sequence pairs. The dark gray volume within the cube is bounded by projections from each of the three surface areas. For more sequences, a similar type of analysis of bounds may be performed in the corresponding higher-order space. 

In practice, MSA calculates the multiple alignment score within the cube lattice by adding the scores of the corresponding pair-wise alignments in the msa. This measure is known as the SP measure (for sum of pairs), and the optimal alignment is based on obtaining the best SP score. These scores may or may not be weighted so as to reduce the influence of more closely related sequences in the msa. The Dayhoff PAM250 matrix and an associated gap penalty are used by MSA for aligning protein sequences. MSA uses a constant penalty for any size of gap and scores gaps according to the scheme illustrated in Figure 4.4 (Altschul 1989; Lipman et al. 1989). MSA calculates a value � for each pair of sequences that provides an idea of how much of a role the alignment of those two sequences plays in the msa. � for a given sequence pair is the difference between the score of the alignment of that pair in the msa and the score of the optimal pair-wise alignment. The bigger the value of �, the more divergent the msa from the pair-wise alignment and the smaller the contribution of that alignment to the msa. For example, if an extra copy of one 



**Figure 4.3.** Bounds within which an optimal alignment will be found by MSA for three sequences. For MSA to find an optimal alignment among three sequences by the DP algorithm, it is only neccessary to calculate optimal alignment scores within the gray volume. This volume is bounded on the one side by the optimal alignments found for each pair of sequences, and on the other by a heuristic multiple alignment of the sequences. The colored areas on each cube surface are twodimensional projections of the gray volume. 

**148** I CHAPTER 4 

||||||Natural|gap cost|Quasi-natural gap cost|
|---|---|---|---|---|---|---|---|
|sequence 1|x|–|–|–|x|||
|sequence 2|x|x|–|x|x|3|4|
|sequence 3|x|x|x|x|x|||



**Figure 4.4.** Method of scoring gap penalties by the msa program MSA. _x_ indicates aligned residues, which may be a match or a mismatch, and – indicates a gap. In this example, each gap cost is 1, regardless of length. The “natural” gap cost is the sum of the number of gaps in all pair-wise combinations (sequences 1 and 2, 1 and 3, and 2 and 3). Note that the alignment of a gap of three in sequence 1 with a gap of length one in sequence 2 scores as gap of 1 because the gap in sequence 1 is longer. The quasi-natural gap cost is the natural cost for the gap plus an additional value for any gap that begins and ends within another. In this example, there is an additional penalty score for the presence of a single gap in sequence 2 that falls within a larger gap in sequence 1. The inclusion of this extra cost for a gap has little effect on the alignments produced but provides an enormous reduction in the amount of information that must be maintained in the DP scoring matrix (Altschul 1989), thus making possible the simultaneous alignment of more sequences by MSA. 

of the sequences is added to the alignment project, then � for sequence pairs that do not include that sequence will increase, indicating a lesser role because the contributions of that pair have been out-voted by the alike sequences (Altschul et al. 1989). Weighting the sequence pairs is designed to get around the common difficulty that some pairs in most sets of sequences are similar. Another score � is the sum of the �s and gives an indication of the degree of divergence among the sequences—closely related sequences will have low �s and �s and distantly related sequences will have high �s and �s. 

The MSA program avoids the bias in an alignment due to alike sequences by weighting the pair-wise scores before they are added to give the SP score. These weights are determined by using the predicted tree of the sequences discussed above. The pair-wise scores between all sequence pairs are adjusted to reduce the influence of the more unlike sequence pairs that occupy more distant “leaves” on the evolutionary tree (i.e., by sequences that are joined by more branches) based on the argument that these sequence pairs provide less useful information for computing the msa. This scheme is different from that used by other msa programs (see below), which generally increase the weight of scores from more distant sequences because these sequences represent greater divergence in the evolutionary tree (see Vingron and Sibbald 1993). 

In using MSA, several additional practical considerations should be considered (described on MSA Web sites given in Table 4.1): (1) MSA is a heavy user of machine resources and is limited to a small number of sequences of relatively short lengths. (2) In the UNIX command line mode of the program, there are options that allow users to specify gap costs, force the alignment of certain residues, specify maximum values for �, and tune the program in other ways. (3) When the output shows that some � are greater than the respective maximum �, a better alignment usually can be found by increasing the maximum � in question. However, increasing � also increases the computational time. (4) If the program bogs down, try dividing the problem into several smaller ones. 

Below is an example from http://www.psc.edu of using MSA to align a group of phospholipase a2 proteins. Note that the program uses the FASTA sequence format. The following steps are used: 

1. Calculate all pair-wise alignment scores (alignment costs). 

2. Use the scores (costs) to predict a tree. 

3. Calculate pair weights based on the tree. 

MULTIPLE SEQUENCE ALIGNMENT I **149** 

4. Produce a heuristic msa based on the tree. 

5. Calculate the maximum � € for each sequence pair. 

6. Determine the spatial positions that must be calculated to obtain the optimal alignment. 

7. Perform the optimal alignment. 

8. Report the � € found compared to the maximum � € . 

##### **Example of MSA** 



<!-- Start of picture text -->
MSA release 2.1 (PSC revision b) started on Thu Jun 19 14:55:31 1997<br>Sequence file format is Fasta.<br>Calculating pairwise alignments.<br>RE<br>Calculating weights.<br>~--------------- Tree given from ancestor ----------------<br>On the left: Internal Node Distance to parent = 278.83<br>On the left: Internal Node Distance to parent = 23.63<br>On the left: Internal Node Distance to parent = 118.62<br>On the left: SEQ#O1 Distance to parent = 230.50<br>On the right: SEQ#04 Distance to parent = 205.50<br>On the right: SEQ#05 Distance to parent = 238.37<br>On the right: SEQ#02 Distance to parent = 256.17<br>On the right: SEQ#03 Distance to parent = 0.00<br>Calculating epsilons.<br>Sequence ID Description<br>1 SEQ#O1 P1;1POA Phospholipase a2 (EC 3.1.1.4) - Chinese cobra<br>2 SEQ#O2 P1;1POD Phospholipase a2 (EC 3.1.1.4) - human<br>3 SEQ#O3 P1;1PPA Phospholipase a2 (EC 3.1.1.4) lys 49 variant<br>4 SEQ#O4 P1;1BPQ phospholipase A2 (EC 3.1.1.4) mutant (K56M) -<br>5 SEQ#O5 P1;1PP2R phospholipase A2 (EC 3.1.1.4) (calcium-free)<br>*** Heuristic Multiple Alignment ***<br>wee RQI5AT JOHIOU IOS UUOO ISO UOOU rE EEHAH35274 #435214 nee<br>NLYQFKNMIQCTVPSR-NLVNFHRMIK -LTTGKEAALSYGFYGCHCGVGGRGSPKDATDRCCVTHDCCYKRLEK  SWWDFADYGCYCGRGGSGTPVDDLDRCCQVHDNCYNEAEKISGC-RGC-- - -~- WPYGTK F KTLS Y SYKF<br>SVLELGKMIL-QETGKNAITSYGSYGCNCGWGHRGQPKDATDRCCFVHKCCYKKLT-~-DC--~--NHKTDRYSY<br>ALWQFNGMIKCKIPSSEPLLDFNNYGCYCGLGGSGTPVDDLDRCCQTHDNCYKQAMKLDSCKVLVDNPYTNNYSY<br>SLVQFETLIM-KIAGRSGLLWY SAYGCYCGWGGHGLPQDATDRCCFVHDCCYGKAT ~~ -DC- ----NPKTVSYTY<br><!-- End of picture text -->

**150** I CHAPTER 4 



<!-- Start of picture text -->
HERRERHARISI1A FERRERRRR ERE ERT AIIS<br>ECSQGTLTCKGGNNACAAAVCDCDRLAA<br>ICFAG-~ -APYNDNDYNINLKARC--~-~-~<br>SNSGSRITC-AKQDSCRSQLCECDKAAATCFARNKTTYNKKYQYYS-NKHCRGSTPRC<br>SWKNKAIIC-EEKNPCLKEMCECDKAVAICLRENLDTYNKKYKAYF<br>-KLKCKKPDT-C<br>SCSNNEITCSSENNACEAF ICNCDRNAAICFSK--VPYNKEHKNLD-KKNC----~-~<br>SEENGEIIC-GGDDPCGTQICECDKAAA<br>ICFRDNI PSYDNKYWLF P-PKDCREEPEPC<br>Calculating pairwise projection costs.<br>reerereery<br>Calculating multiple alignment.<br>coccihsoccBococlooscllococBo00cBo00sToc00bes0Mooca®)<br>JOSS IIIS S SII SIIS I IUISIS GIGI JUDO UIC TOIT III<br>*** Optimal Multiple Alignment ***<br>NLYQFKNMIQCTVPSR- SWWDFADYGCYCGRGGSGTPVDDLDRCCQVHDNCYNEAEKISGC<br>~~ ---WPYFKTYSY<br>NLVNFHRMIK-LTTGKEAALSYGFYGCHCGVGGRGSPKDATDRCCVTHDCCYKRLEK-RGC-~---GTKFLSYKF<br>SVLELGKMIL-QETGKNATTSYGSYGCNCGWGHRGQPKDATDRCCFVHKCCYKKL-~-TDC--~--NHKTDRYSY<br>ALWQFNGMIKCK1PSSEPLLDFNNYGCYCGLGGSGTPVDDLDRCCOTHDNCYKQAMKLDSCKVLVDNPYTNNYSY<br>SLVQFETLIM-KIAGRSGLLWYSAYGCYCGWGGHGLPQDATDRCCF<br>VHDCCYGKA- - - TDC---~-NPKTVSYTY<br>ECSQGTLTCKGGNNACAAAVCDCDRLAAICFAG--APYNDNDYNINLKARC---~----<br>SNSGSRITC-AKQDSCRSQLCECDKAAATCFARNKTTYNKKYQYYS-NKHCRGSTPRC<br>SWKNKAIIC-EEKNPCLKEMCECDKAVAICLRENLDTYNKKYKAYF<br>-KLKCK-KPDTC<br>SCSNNEITCSSENNACEAF<br>ICNCDRNAAICFSK - -VPYNKEHKNLD-KKNC-~~----~<br>SEENGEIIC-GGDDPCGTQICECDKAAAICFRDNI<br>PSYDNKYWLF P-PKDCREEPEPC<br>End gaps not penalized.<br>Costfile: pam250<br>Alignment cost: 35132 Lower bound: 34945<br>Delta: 187 Max. Delta: 285<br>Sequences Proj. Cost Pair. Cost Epsilon Max. Epsi. Weight Weight*Cost<br>1 2 1864 1825 39 39 1 1864<br>1 3 1891 1843 48 57 1 1891<br>1 4 1654 1653 1 5 4 6616<br>1 5 1814 1787 27 28 2 3628<br>2 3 1735 1733 2 8 4 6940<br>2 4 1876 1866 10 10 1 1876<br>2 5 1713 1712 1 8 2 3426<br>3.04 1901 1889 12 21 1 1901<br>3.5 1648 1648 0 11 2 3296<br>4 5 1847 1842 5 6 2 3694<br>Elapsed time = 0.895<br>Tree<br>A tree is given for the heuristic alianment (not shown).<br><!-- End of picture text -->

MULTIPLE SEQUENCE ALIGNMENT I **151** 

### **~~SCORING MULTIPLE SEQUENCE ALIGNMENTS~~** 

As discussed above, the SP method provides a way to score the msa by summing the scores of all possible combinations of amino acid pairs in a column of a msa. The method assumes a model for evolutionary change in which any of the sequences could be the ancestor of the others, as illustrated in Figure 4.5. This figure also illustrates a difficulty with the SP method when a substitution table of log odds scores such as BLOSUM62 is used for protein sequences (see Durbin et al. 1998, pp. 139–140). Shown is the effect of adding a small number of amino acid subsitutions to a column that initially has all matching amino acids. Scores in the msa column decrease rapidly as the number of mismatched residue pairs increases. For a larger number of sequences than five with all N, or with one or two C substitutions, these decreases should be greater because there will be more N-N matched pairs relative to mismatched N-C pairs. However, the reverse is true with the SP method of scoring. For _n_ sequences, the number of combinations of pairs in a column is 



<!-- Start of picture text -->
Sequence Column A Column B Column C<br>1 ......N................N...............N<br>2 ......N................N...............N<br>3 ......N................N...............N<br>4 ......N................N...............C<br>5 ......N................C...............C<br>N N N<br>N N N N N C<br>N N N C N C<br>Column A Column B Column C<br>No. of N–N matched pairs (each scores 6):<br>10 6 4<br>No. of N–C matched pairs (each scores –3):<br>0 4 6<br>BLOSUM62 score :<br>60 24 6<br><!-- End of picture text -->

**Figure 4.5.** The SP model for scoring a msa. This model represents one method for optimizing the msa by maximizing the number of matched pairs (or minimizing the cost or number of mismatched pairs) summed over all columns in the msa. Shown first are three columns of a five-sequence msa with all matched (A), four matched and one mismatched (B), or three matched and two mismatched (C) sequence characters. The SP method of calculating the cumulative scores for columns of a msa is then illustrated by a graph with the five sequences as vertices and representing the ten possible sequence pair-wise sequence comparisons. Solid lines represent a matched pair and dotted lines a mismatched pair. Shown are the BLOSUM62 scores for each column calculated by the SP method. (Adapted from Altschul 1989.) 

**152** I CHAPTER 4 



<!-- Start of picture text -->
A. B.<br>N N<br>N C<br>N C<br>N C N<br>N<br>N C N C<br><!-- End of picture text -->

**Figure 4.6.** Alternative methods for scoring a column in the msa (Altschul 1989b). The variations in column C of Fig. 4.5 are shown modeled by a phylogenetic tree ( _A_ ) and a simplified phylogenetic tree called a star phylogeny ( _B_ ) where one of the sequences is treated as the ancestor of all the others (instead of treating them as all equally possible ancestors as in the original sum of pairs scoring method). 

_n_ ( _n_ � 1)/2. If all are amino acid N, as in column A, then the BLOSUM62 score for the column is 6 � _n_ ( _n_ � 1)/2. If there is one C in the column, as in column B, then _n_ � 1 matched N-N pairs will be replaced by _n_ �1 mismatched N-C pairs, giving a score of 9( _n_ � 1) less. The score for one C in the column divided by that for zero Cs is 9( _n_ � 1)/[6 _n_ ( _n_ � 1)/2] � 3/ _n_ . For three sequences, the relative difference is 1, whereas for six sequences, the relative difference is 2. As more sequences are present in the column, the relative difference increases, not in agreement with expectation. Hence, the SP method is not providing a reasonable result when this type of scoring matrix is used. Two other methods for scoring a msa (Altschul 1989) have been described and are illustrated in Figure 4.6. The first is a treebased method. Because a phylogenetic tree describing the relationships among the sequences is found by the MSA program, the sum of the lengths of the tree branches can be calculated using the substitutions in the column of the msa. Alternatively, a simplified tree with one of the sequences as the ancestor of all of the others (a star phylogeny) can also be used (see Chapter 6). msa programs using these methods have not been implemented. Other scoring methods include information content (see p. 195) and a graph-based method called the trace method (Kececioglu 1993). A novel branch-and-cut algorithm for msa has been developed based on the trace method (Kececioglu et al. 2000). Other methods of scoring and producing an alignment guided by a tree are described below. 

### **~~PROGRESSIVE METHODS OF MULTIPLE SEQUENCE ALIGNMENT~~** 

The MSA program described above for obtaining an optimal alignment of multiple sequences is limited to three sequences or to a small number (six to eight) of relatively short sequences. Progressive alignment methods use the dynamic programming method to build a msa starting with the most related sequences and then progressively adding lessrelated sequences or groups of sequences to the initial alignment (Waterman and Perlwitz 1984; Feng and Doolittle 1987, 1996; Thompson et al. 1994a; Higgins et al. 1996). Relationships among the sequences are modeled by an evolutionary tree in which the outer branches or leaves are the sequences (Fig. 4.7). The tree is based on pair-wise comparisons of the sequences using one of the phylogenetic methods described in Chapter 6. Progenitor sequences represented by the inner branches of the tree are derived by alignment of the outermost sequences. These inner branches will have uncertainties where positions in the 

MULTIPLE SEQUENCE ALIGNMENT I **153** 



<!-- Start of picture text -->
N Y L S N K Y L S N F S N F L S<br>N  K / – Y  L  S N F L / – S<br>N  K / – Y/F  L / – S<br><!-- End of picture text -->

**Figure 4.7.** Progressive sequence alignment. Sequences are represented as the outermost branches (leaves) on an evolutionary tree. The most closely related sequences are first aligned by dynamic programming, providing a representation of ancestor sequences in deeper branches with uncertainties where amino acids have been substituted or positioned opposite a gap. These sequences are the same as those shown in EVMSA. The challenge to the msa method is to utilize an appropriate combination of sequence weighting, scoring matrix, and gap penalties so that the correct series of evolutionary changes may be found. 

outermost sequences are dissimilar, as illustrated in Figure 4.7. Two examples of programs that use progressive methods are CLUSTALW and the Genetics Computer Group program PILEUP. 

### **CLUSTALW** 

CLUSTAL has been around for more than 10 years, and the authors have done much to support and improve the program (Higgins and Sharp 1988; Thompson et al. 1994a; Higgins et al. 1996). CLUSTALW is a more recent version of CLUSTAL with the W standing for “weighting” to represent the ability of the program to provide weights to the sequence and program parameters, and CLUSTALX provides a graphic interface (see Table 4.1). These changes provide more realistic alignments that should reflect the evolutionary changes in the aligned sequences and the more appropriate distribution of gaps between conserved domains. 

CLUSTAL performs a global-multiple sequence alignment by a different method than MSA, although the initial heuristic alignment obtained by MSA is calculated the same way. The steps include: (1) Perform pair-wise alignments of all of the sequences; (2) use the alignment scores to produce a phylogenetic tree (for an explanation of the neighbor-joining method that is used, see Chapter 6); and (3) align the sequences sequentially, guided by the phylogenetic relationships indicated by the tree. Thus, the most closely related sequences are aligned first, and then additional sequences and groups of sequences are added, guided by the initial alignments to produce a msa showing in each column the sequence variations among the sequences. The initial alignments used to produce the guide tree may be obtained by a fast _k_ -tuple or pattern-finding approach similar to FASTA that is useful for many sequences, or a slower, full dynamic programming method may be used. An enhanced dynamic programming alignment algorithm (Myers and Miller 1988; see book Web site) is used to obtain optimal alignment scores. For producing a phylogenetic tree, genetic distances between the sequences are required. The genetic distance is the number of mismatched positions in an alignment divided by the total number of matched positions (positions opposite a gap are not scored). 

As with MSA, sequence contributions to the msa are weighted according to their relationships on the predicted evolutionary tree. A rooted tree with known branch lengths of which the sequences are outer branches (leaves) is examined (see Chapter 6). Weights are 

**154** I CHAPTER 4 

based on the distance of each sequence from the root, as illustrated in Figure 4.8. The alignment scores between two positions in the msa are then calculated using the resulting weights as multiplication factors. 

The scoring of gaps in a msa has to be performed in a different manner from scoring gaps in a pair-wise alignment. As more sequences are added to a profile of an existing msa, gaps accumulate and influence the alignment of further sequences (Thompson et al. 1994b; Taylor 1996). CLUSTALW calculates gaps in a novel way designed to place them between conserved domains. When Pascarella and Argos (1992; see book Web site **)** aligned sequences of structurally related proteins, the gaps were preferentially found between secondary structural elements. These authors also prepared a table of the observed frequency 



<!-- Start of picture text -->
A. Calculation of sequence weights<br>Weighting factor<br>0.2<br>A   0.2 + 0.3/2 = 0.35<br>0.3<br>0.1<br>B   0.1 + 0.3/2 = 0.25<br>0.5<br>C   0.5<br>B. Use of sequence weights<br>Column in alignment 1<br>Sequence A (weight a) ………K………<br>Sequence B (weight b) ………I………<br>Column in alignment 2<br>Sequence C (weight c) ………L………<br>Sequence D (weight d) ………V………<br>Score for matching these two column in an msa =<br>[ a x c x score (K,L) +<br>  a x d x score (K,V) +<br>  b x c x score (I,L) +<br>  b x d x score (I,V) ] / 4<br><!-- End of picture text -->

**Figure 4.8.** Weighting scheme used by CLUSTALW (Higgins et al. 1996). ( _A_ ) Sequences that arise from a unique branch deep in the tree receive a weighting factor equal to the distance from the root. Other sequences that arise from branches shared with other sequences receive a weighting factor that is less than the sum of the branch lengths from the root. For example, the length of a branch common to two sequences will only contribute one-half of that length to each sequence. Once the specific weighting factors for each sequence have been calculated, they are normalized so that the largest weight is 1. As CLUSTALW aligns sequences or groups of sequences, these fractional weights are used as multiplication factors in the calculation of alignment scores. ( _B_ ) Illustration of using sequence weights for aligning two columns in two separate alignments. Note that this sequence weighting scheme is the opposite to that used by MSA, because the more distant a sequence from the others, the higher the weight given. For a comparison of additional weighting schemes, see Vingron and Sibbald (1993). 

MULTIPLE SEQUENCE ALIGNMENT I **155** 

of gaps next to each amino acid in these regions. CLUSTALW uses the information in this table and also attempts to locate what may be the corresponding domains by appropriate gap placement in the msa. Like other alignment programs, CLUSTAL uses a penalty for opening a gap in a sequence alignment and an additional penalty for extending the gap by one residue. These penalties are user-defined (defaults are available). Gaps found in the initial alignments remain fixed. New gaps introduced as more sequences are added also receive this same gap penalty, even when they occur within an existing gap, but the gap penalties for an alignment are then modified according to the average match value in the substitution matrix, the percent identity between the sequences, and the sequence lengths (Higgins et al. 1996). These changes are attempts to compensate for the scoring matrix, expected number of gaps (alignment with more identities should have fewer gaps), and differences in sequence length (should limit placement of gaps if one sequence shorter). Tables of gaps are then calculated for each group of sequences to be aligned to confine them to less conserved regions in the alignment. Gap penalties are decreased where gaps already occur (another method for achieving this same result is to enhance the scores of more closely matching regions on the alignment as described in Taylor 1996), increased in regions adjacent to already gapped regions, decreased within stretches of hydrophilic regions (amino acids DEGKNQPRS), and increased or decreased according to the table in Pascarella and Argos (1992). These rules are most useful when a correct alignment of some of the sequences is already known. The CLUSTALW algorithm and the results of using the above sequence weighting gap adjustment method are illustrated in Figure 4.9. 

CLUSTALW also has options for adding one or more additional sequences with weights or an alignment to a existing alignment (Higgins et al. 1996). Once an alignment has been made, a phylogenetic tree may be made by the neighbor-joining method, with corrections for possible multiple changes at each counted position in the alignment (see Chapter 6). The predicted trees may also be displayed by various programs described in Chapter 6. 

### **PILEUP** 

PILEUP is the msa program that is a part of the Genetics Computer Group package of sequence analysis programs, owned since 1997 by Oxford Communications, and is widely used due to the popularity and availability of this package. PILEUP uses a method for msa that is very similar to CLUSTALW. The sequences are aligned pair-wise using the Needleman-Wunsch dynamic programming algorithm, and the scores are used to produce a tree by the unweighted pair-group method using arithmetic averages (UPGMA; Sneath and Sokal 1973 and see Chapter 6). The resulting tree is then used to guide the alignment of the most closely related sequences and groups of sequences. The resulting alignment is a global alignment produced by the Needleman-Wunsch algorithm. Standard scoring matrices and gap opening/extension penalties are used. Unfortunately, there have not been any recent enhancements of this program such as gap modifications or sequence weighting comparable to those introduced for CLUSTALW. As with other progressive alignment msa programs, PILEUP does not guarantee an optimal alignment. 

### **Problems with Progressive Alignment** 

The major problem with progressive alignment programs such as CLUSTAL and PILEUP is the dependence of the ultimate msa on the initial pair-wise sequence alignments. The very first sequences to be aligned are the most closely related on the sequence tree. If these sequences align very well, there will be few errors in the initial alignments. However, the more distantly related these sequences, the more errors will be made, and these errors will 

**156** I CHAPTER 4 



<!-- Start of picture text -->
Hbb_Human 1<br>Hbb_Horse 2 7<br>HbaHuman 3 2 Pairwise alignment:<br>HbaHorse 4[ 59 59 13 - Calculate dictance matrix<br>MygPhyca 5 > UB culate distance ma<br>Glb5_Petma 6 81 82 73 74 80 -<br>Lgb2_Luplu 7 87-86 868893<br>1 2 3 4 5 6<br>Hbb_Human<br>Hbb_Horse<br>Hba_Human<br>Hba_Horse Rooted Neighbortree (guide4, tree)  Joining<br>Myg_Phyca<br>Glb5_Petma<br>Lgb2_Luplu |<br>o-7o077-VHLAPEEKSAVIRLWGRIN- -{(DEVOGERLGRLEW MTOR FESFGDLST — [-y -<br>~_-----voLgSEEKAAVLALWDKWN--EEEVGGEALGALLVWyHwToRPFDSFGpLsn | Progressive<br>~==------VLGPADKTNVKAAWGKYGAHAGEYGAEALERMFLSPHITK VFPuEDLS-- | alignment: —<br>=== ==----VLGAADKTNVEAANS! AGEYGAEALERMFLGPHT<qyrPHrDs-- | Align following<br>---------VLGEGEWOLVLHVWAKVEATVAGHGQDILIRLFKSHFETLEKFDRFKELKT | the guide tree<br>PIVDTGSVAPLGAAEKTKIRSAWAPW'S1WETSGVDILVKFFTSTHAAQHEF PKFKGLTT<br>~~-~----GALTESQAALVKSSWEEENANPKETER TRARAROLFSFLKGTSE<br>PDAVMGNERVRARGRRVLGAF SDGIBHLD-<br>- ---NLRGTFA ACDKLEVOPENFRE<br>PGAVNGNPVKAHGKKVLHSFGEGVHHLD-<br>~~ --NLKGTFAALSELHCOKLHVEPENERI,<br>~~ --HGSROVKGHGKKVADALTNAVA )- - -- - DMPNALSALSDLHAHKLR' VNEKL<br>~~~ -HGSROVKAHGKEVGDALTLAVGHLD-----DLPGALSNLSDLHAHKLRVGPVNEKE<br>EAEMKASEDLKKHGVTVLTALGAILKKKG-~<br>~~ -HHEAELKPLAQSHATKHKT AIKYLEF<br>VP--QNNPELOAHAGKVFKLVYEAAMTADQLKKSPDVRWHAERT INAVNDAVASMDDT- OLOVTGVVVTDATLEN- EKMEMK L RDLSGKHAMSQVEPQYFKVGSVHVSKG-Vi Vy<br>TGNVLVCVLAH GKEF TPPVOARYORUURGVANALAHKYH-~--~-<br>LGNVLVVVLARHFGKDFT PELQASYQKVVAGVANALAHKYH- ----~<br>LSHCLLVTLAAHLPAEFTPAVHABLDKFLASVSTVLTSKYR-----~<br>LSHCLLSTLAVHLPNDFTPAVHASLDKFLSSVSTVLIBRYR-~-~--<br>ISEATTHVLSHHPGDFGADAQGAMNKALELFRKDIAAKYKELGYOG<br>LAAVIADTVAAG. ~~ --~-~-DAGFEKLMSMICILLABAY-----~-<br>VKEAILRTIKEWWGAKWSEELNSAWTIAYDELAIVIKKEMNDAA-—~<br><!-- End of picture text -->

**Figure 4.9.** A msa of seven globins by CLUSTALW. The protein identifiers are from the SwissProt database. The amino acid subsitution matrix was the Dayhoff PAM250 matrix, and gap penalties were varied to emphasize conserved ungapped regions. The approximate and known locations of seven � @ -helices in the structure of this group are shown in boxes. (Reprinted, with permission, from Higgins et al. 1996 [copyright Academic Press].) 

be propagated to the msa. There is no simple way to circumvent this problem. A second problem with the progressive alignment method is the choice of suitable scoring matrices and gap penalties that apply to the set of sequences (Higgins et al. 1996). 

For the difficult task of aligning more distantly related sequences, using Bayesian methods such as hidden Markov models (HMMs) may be useful. For more closely related 

MULTIPLE SEQUENCE ALIGNMENT I **157** 

sequences, CLUSTALW is designed to provide an adequate alignment of a large number of sequences and provide a very good indication of the domain structure of those sequences. 

### **~~ITERATIVE METHODS OF MULTIPLE SEQUENCE ALIGNMENT~~** 

The major problem with the progressive alignment method described above is that errors in the initial alignments of the most closely related sequences are propagated to the msa. This problem is more acute when the starting alignments are between more distantly related sequences. Iterative methods attempt to correct for this problem by repeatedly realigning subgroups of the sequences and then by aligning these subgroups into a global alignment of all of the sequences. The objective is to improve the overall alignment score, such as a sum of pairs score. Selection of these groups may be based on the ordering of the sequences on a phylogenetic tree predicted in a manner similar to that of progressive alignment, separation of one or two of the sequences from the rest, or a random selection of the groups. These methods are compared in Hirosawa et al. (1995). 

MultAlin (Corpet 1988) recalculates pair-wise scores during the production of a progressive alignment and uses these scores to recalculate the tree, which is then used to refine the alignment in an effort to improve the score. The program PRRP (Table 4.1) uses iterative methods to produce an alignment. An initial pair-wise alignment is made to predict a tree, the tree is used to produce weights for making alignments in the same manner as MSA except that the sequences are analyzed for the presence of aligned regions that include gaps rather than being globally aligned, and these regions are iteratively recalculated to improve the alignment score. The best scoring alignment is then used in a new cycle of calculations to predict a new tree, new weights, and new alignments, as illustrated in Figure 4.10. The process is repeated until there is no further increase in the alignment score (Gotoh 1994, 1995, 1996). 

The program DIALIGN (see Table 4.1) finds an alignment by a different iterative method. Pairs of sequences are aligned to locate aligned regions that do not include gaps, much like continuous diagonals in a dot matrix plot. Diagonals of various lengths are identified. A consistent collection of weighted diagonals that provides an alignment which is a maximum sum of weights is then found. 

Additional methods that use iterative procedures are described below. 

### **Genetic Algorithm** 

The genetic algorithm is a general type of machine-learning algorithm that has no direct relationship to biology and that was invented by computer scientists. The method has been recently adapted for msa by Notredame and Higgins (1996) in a computer program package called SAGA (Sequence Alignment by Genetic Algorithm; see Table 4.1). Zhang and Wong (1997) have developed a similar program. The method is of considerable interest because the algorithm can find high-scoring alignments as good as those found by other methods. Similar genetic algorithms have been used for RNA sequence alignment (Notredame et al. 1997) and for prediction of RNA secondary structure (Shapiro and Navetta 1994). Although the method is relatively new and not used extensively, it likely represents the first of a series of sequence analysis programs that produce alignments by attempted simulation of the evolutionary changes in sequences. 

The basic idea behind this method is to try to generate many different msas by rearrangements that simulate gap insertion and recombination events during replication in 

**158** I CHAPTER 4 



<!-- Start of picture text -->
Initial alignment<br>\ 1 ot<br>Phylogenetic<br>2 tree<br>©<br>Weights to.<br>wij Sequence Pairs<br>1sicgjsN<br>| wsP= >ea wij Siz<br>84Ba DGMASGLDKDYLKPDDSR-"DONAEOY SOaR-VTA  { tterativeie  Refinementea  of<br>44 ---------DSIPSGVDASKI- Muttiple Alignment<br>No<br>Converged?<br>Yes<br><!-- End of picture text -->

**Figure 4.10.** The iterative procedures used by PRRP to compute a multiple sequence alignment. (Reprinted, with permission, from Gotoh 1996 [copyright Academic Press].) 

order to generate a higher and higher score for the msa. The alignments are not guaranteed to be optimal or to be the highest scoring that is achievable (optimal alignment). Although SAGA can generate alignments for many sequences, the program is slow for more than about 20 sequences. 

A similar approach for obtaining a higher-scoring msa by rearranging an existing alignment uses a probability approach called simulated annealing (Kim et al. 1994). The program MSASA (Multiple Sequence Alignment by Simulated Annealing) starts with a heuristic msa and then changes the alignment by following an algorithm designed to identify changes that increase the alignment score. 

The success of the genetic algorithm may be attributed to the steps used to rearrange sequences, many of which might be expected to have occurred during the evolution of the protein family. The steps in the algorithm are as follows: 

1. The sequences to be aligned (up to �20 in number) are written in rows, as on a page, except that they are made to overlap by a random amount of sequence, up to 50 residues long for sequences about 200 in length. The ends are then padded with gaps. A typical population of 100 of these msas is made, although other numbers may be set. 

```
xxxxxxxxxx----
---xxxxxxxxxxx
-xxxxxxxxxx---
```

MULTIPLE SEQUENCE ALIGNMENT I **159** 

Shown is an initial msa for the genetic algorithm (1 of �100 in number). 

2. The 100 initial msas are scored by the sum of pairs method, except that both natural and quasi-natural gap-scoring schemes (Fig. 4.4) are used. Recall that the best SSP score for a msa is the minimum one and the one that is closest to the sum of the pair-wise sequence alignment. Standard amino acid scoring matrices and gap opening and extension penalties are used. 

3. These initial msas are now replicated to give another generation of msas. The half with the lowest SSP scores are sent to the next generation unchanged. The remaining half for the next generation are selectively chosen by lot, like picking marbles from a bag, except that the chance for a particular choice is inversely proportional to the msa score (the lower the score, the better the msa, therefore gives that one a greater chance of replicating). These latter one-half of the choices for the next generation are now subject to mutation, as described in step 4 below, to produce the children of the next generation. All members of the next-generation msas undergo recombination to make new child msas derived from the two parents, as described in step 5 below. The relative probabilities of these separate events are governed by program parameters. These parameters are also adjusted dynamically as the program is running to favor those processes that have been most useful for improving msa scores. 

4. In the mutation process, the sequence is not changed (else it would no longer be an alignment), but gaps are inserted and rearranged in an attempt to create a better-scoring msa. In the gap insertion process, the sequences in a given msa are divided into two groups based on an estimated phylogenetic tree, and gaps of random length are inserted into random positions in the alignment. Alternatively, in a “hill-climbing” version of the procedure, the position is so chosen as to provide the best possible score following the change. 

```
xxxxxxxxxx                   xxx--xxxxxxx
xxxxxxxxxx                   xxx--xxxxxxx
xxxxxxxxxx        ---->      xxxxxxxxx--x
xxxxxxxxxx                   xxxxxxxxx--x
xxxxxxxxxx                   xxxxxxxxx--x
```

Shown above are random gap insertions into phylogenetically related sequences. The first two and last three sequences comprise the two related groups in this example. x indicates any sequence character. 

Another mutational process is to move common blocks of sequence (overlapping ungapped regions) delineated by a gap, or blocks of gaps (overlapping gaps). Some of the possible moves are illustrated below. These moves may also be tailored to improve the alignment score. 

|`xxx--xxxxx`|`xx--xxxxxx`|`xxx--xxxxx`|`xxxxx--xxx`|
|---|---|---|---|
|`xxxxxxxxxx`|`xxxxxxxxxx`|`xxxxxxxxxx`|`xxxxxxxxxx`|
|`xx--xxxxxx`|`x--xxxxxxx`|`xxx--xxxxx`|`xx-xx-xxxx`|
|`xxxxxxxxxx`|`xxxxxxxxxx`|`xxxxxxxxxx`|`xxxxxxxxxx`|
|Starting block|Whole block<br>move|Split block<br>horizontally<br>(guided by<br>phylogenetic gro|Split block<br>vertically<br>uping)|



**160** I CHAPTER 4 

5. Recombination among next-generation parent msas is accomplished by one of two mechanisms. The first is not homology-driven. One msa is cut vertically through, and the other msa is cut in a staggered manner that does not lose any sequence after the fragments are spliced. The higher scoring of the two reciprocal recombinants is kept. The second, illustrated below, is recombination between msas driven by conserved sequence positions. It is driven by homology expressed as a vertical column of the same residue and is very like standard homologous recombination. 

```
xxGxxxxDxxxxGxx-xDxxxxGxx-xDxx
xxGx-xxDxxxxGxxxxDxxxxGxxxxDxx
xxGxx-xDxxxxGxxxxDxxxxGxxxxDxx
xxGxxxxDxxxxGx-xxDxxxxGx-xxDxx
```

**Parent A Parent B Child** alignment alignment alignment 

6. The next generation, an overlapping one of the previous one-half of the best-scoring parental msas and the mutated children, is now evaluated as in step 2, and the cycle of steps 2–5 is typically repeated as much as 100 times, although as many as 1000 generations can be run. The best-scoring msa is then obtained. 

7. The entire process of producing a set of msas for replication and mutation is repeated several times to obtain several possible msas, and the best scoring is chosen. 

### **Hidden Markov Models of Multiple Sequence Alignment** 

The HMM is a statistical model that considers all possible combinations of matches, mismatches, and gaps to generate an alignment of a set of sequences. A localized region of similarity, including insertions and deletions, may also be modeled by an HMM. Analysis of sequences by an HMM is discussed on page 185 along with other statistical methods. 

### **<u>OTHER PROGRAMS AND METHODS FOR MULTIPLE SEQUENCE ALIGNMENT</u>** 

The msa method often used, especially for 10 or more sequences, is to first determine sequence similarity between all pairs of sequences in the set. On the basis of these similarities, various methods are used to cluster the sequences into the most related groups or into a phylogenetic tree. 

In the group approach, a consensus is produced for each group and then used to make further alignments between groups. Two examples of programs using the group approach are the program PIMA (Smith and Smith 1992), which uses several novel alignment techniques, and the program MULTAL described by Taylor (1990, 1996; see Table 4.1). 

The tree method uses the distance method of phylogenetic analysis to arrange the sequences. The two closest sequences are then aligned, and the resulting consensus alignment is aligned with the next best sequence or cluster of sequences, and so on, until an alignment is obtained that includes all of the sequences. The programs PILEUP and CLUSTALW discussed above are examples. The ALIGN set of programs (Feng and Doolittle 1996) and the MS-DOS program by Corpet (1988) use this method. Additional programs for msa are also described in Barton (1994), Kim et al. (1994), and Morgenstern et al. (1996). 

Another program (Vingron and Argos 1991) aligns all possible pairs of sequences to cre- 

MULTIPLE SEQUENCE ALIGNMENT I **161** 

ate a set of dot matrices, and the matrices are then filtered sequentially to find motifs that provide a starting point for sequence alignment. A set of programs for interactive msa by dot matrix analysis and other alignment techniques has also been developed (Boguski et al. 1992). 

The program TREEALIGN takes the approach that multiple sequence alignments should be done in a fashion that simultaneously minimizes the number of changes needed during evolution to generate the observed sequence variation (Hein 1990). TREEALIGN (also named ALIGN in the program versions) has a method for performing the alignment and the most parsimonious tree construction at the same time. The initial steps are similar to other multiple sequence alignment methods, except for the use of a distance scale: i.e., the sequences are aligned pair-wise and the resulting distance scores are used sequentially to produce a tree, which is rearranged as more sequences are added. The sequences are then realigned so that the same tree can be produced by maximum parsimony. Finally, the tree is rearranged to maximize parsimony. The advantage to this method is the increased use of phylogenetic analysis to improve the multiple sequence alignment. 

### **~~LOCALIZED ALIGNMENTS IN SEQUENCES~~** 

Multiple sequence alignment programs based on the methods discussed above report a global alignment of the sequences, including all parts of all sequences. A portion of the alignment that is highly conserved may then be identified and a type of scoring matrix called a profile may be produced. A profile includes scores for amino acid substitutions and gaps in each column of the conserved region so that an alignment of the region to a new sequence can be determined. Alternatively, the alignment may be scanned for regions that include only substituted regions without gaps, called blocks, and these blocks may then be used in sequence alignments. 

There is also a third method for finding a localized region of sequence similarity in a set of sequences without first having to produce an alignment. In this method, the sequences are analyzed by pattern-searching or statistical methods. All of these methods for finding localized sequence similarity are discussed below. 

### **Profile Analysis** 

Profiles are found by performing the global msa of a group of sequences and then removing the more highly conserved regions in the alignment into a smaller msa. A scoring matrix for the msa, called a profile, is then made. The profile is composed of columns much like a mini-msa and may include matches, mismatches, insertions, and deletions. A tutorial on preparing profiles by the first method, prepared by M. Gribskov, is at Web address http://www.sdsc.edu/projects/profile/profile_tutorial.html, and the Web site at http://www.sdsc.edu/projects/profile/ will perform a motif analysis on the University of California at San Diego Supercomputer Center. The program Profilemake can be used to produce a profile from a msa (Gribskov et al. 1987, 1990; Gribskov and Veretnik 1996). A version of the Profilesearch program, which performs a database search for matches to a profile, is available at the University of Pittsburgh Supercomputer Center (http://www.psc.edu/general/software/packages/profiless/profiless.html). A special grant application may be needed to use this facility. Profile-generating programs are available by FTP from ftp.sdsc.edu/pub/sdsc/biology and are included in the Genetics Computer Group suite of programs (http://www.gcg.com/), although the more recent features (Gribskov and Veretnik 1996) are not included in GCG, v. 9.1. 

**162** I CHAPTER 4 

Once produced, the profile is used to search a target sequence for possible matches to the profile using the scores in the table to evaluate the likelihood at each position. For example, the table value for a profile that is 25 amino acids long will have 25 rows of 20 scores, each score in a row for matching one of the amino acids at the corresponding position in the profile. If a sequence 100 amino acids in length is to be searched, each 25amino-acid-long stretch of sequence will be examined, 1–25, 2–26, . . . . 76–100. The first 25-amino-acid-long stretch will be evaluated using the profile scores for the amino acids in that sequence, then the next 25-long stretch, and so on. The highest-scoring sections will be the most similar to the profile. 

The disadvantage of this method of profile extraction from an msa is that the profile produced is only as representative of the variation in the family of sequences as the msa itself. If several sequences in the msa are similar, the msa and the derived profile will be biased in favor of those sequences. Methods have been devised for partially circumventing this problem with the profile (Gribskov and Veretnik 1996), but the difficulty with the msa itself is not easily reconciled, as discussed at the beginning of this section. Sequence weighting is based on the production of a simple phylogenetic tree by distance methods; more closely related sequences then receive a reduced weight in the profile. Another problem is that some amino acids may not be represented in a particular column because not enough sequences have been included. Athough absence of an amino acid may mean that the amino acid may not occur at that position in the protein family, adding counts to such positions generally increases the usefulness of the profile. This feature is built into the profile method discussed below. 

An example of the generation of a profile and the matrix representation of this profile for a set of heat shock proteins is illustrated in Figure 4.11. The profile is similar to the log odds form of the amino acid substitution table, such as the PAM250 and BLOSUM62 



<!-- Start of picture text -->
cons A B C D E FOG HI K L MN P Q RS T V WwW YX Z Gap Len<br>I 8 3-2 5 4 5 5 -4 24 O 1513 1 1 1 -7 2 22 21-18 -6 4 100 100<br>T 13:19 -5 24 18 -18 19 71 7 -7-41411 10 -1 9 29 3-28 -14 15 100 100<br>L 5 5 -5 3 4 13 4 2 8 -4 1412 8-5 0-10 0 1010 -1 5 2 22 22<br>s 17:14:17 13 10-12 29 -5 -5 6 -14 -9 12 10 0 -2 34 19 1 -8 -15 4 100 100<br>T 15 3 22 0 -1 -5 12 -2 7 -3 -8 -6 5 7 -8 -7 16 29 9 =-22 6 -4 100 100<br>T 8-112 -2 0 5 6 -419 -4 8 5-1 2 -8 -8 7 22 19 -15 4 -3 100 100<br>c 17 024 -1 -3 11 8 -1 7 -10 1-2 1-3 -8 -14 8 5 9-5 14 -7 100 100<br>v 11 018 -1 -2 214 -10 26 -4 9 7-3 7 -7 -7 21 #10 31-19 -5 -5 100 100<br>c 10-8 15 -11 -11 6 8 -7 11 -10 4 3-7 0 -11 -4 11 5 15 -22 14 -11 100 100<br>v 7 7-3 8 8 -3 11 120 -1 1410 4 2 8-5 0 5 26 -24 -6 8 100 100<br><!-- End of picture text -->

**Figure 4.11.** Pattern identification by the profile method. A set of heat shock 70 (hsp70) proteins from a diverse group of organisms were aligned by the Genetics Computer Group msa program PILEUP. A profile was then made from one region in the alignment with the Genetics Computer Group program Profilemake. The profile represents the specific motif pattern found for the chosen location shown for this set of hsp70 proteins. The first column gives the consensus amino acid at each position in the profile. Thus, the consensus pattern is ITLSTTCVCV. This profile is used to search a target sequence for matches to the profile. The table values are a log odds score of giving the probability of finding the amino acid in the target sequence at that position in the profile divided by the probability of aligning the two amino acids by random chance. If a gap must be placed in the target sequence to align the sequence with the profile, then the penalties for opening a gap and extending the gap, respectively, are subtracted. The profile itself may include gaps, in which case the penalty is reduced, as seen for example in the row 3 of the profile table. The method of producing the substitution scores shown in the table is described in the text. 

MULTIPLE SEQUENCE ALIGNMENT I **163** 

matrices used for sequence alignments. The matrix is 23 columns wide, one column for each of the 20 amino acids, plus one column for an unknown amino acid **_z_** and two columns for a gap opening and extension penalty. There is one row for each column in the msa. The consensus sequence, derived from the most common amino acid in each column of the msa, is listed down the left-hand column. The scores on each row reflect the number of occurrences of each amino acid in the aligned sequences. For example, in the first row, I, T, and V were found, with I being the majority amino acid. The highest positive score on each row (underlined) is in the column corresponding to the consensus amino acid, the most negative score for an amino acid not expected at that position. These values are derived from the log odds amino acid substitution matrix that was used to produce the alignment, such as the log odds form of the Dayhoff PAM250 matrix. Two methods are used to produce profile tables, the average method and the evolutionary method. The evolutionary method seems somewhat better for finding family members. 

In the average method, the profile matrix values are weighted by the proportion of each amino acid in each column of the msa. For example, if column 1 in the msa has 5 Ile (I), 3 Thr (T), and 2 Val (V), then the frequency of each amino acid in this column is 0.5 I, 0.3 T, and 0.2 V. These amino acids are considered to have arisen with equal probability from any of the 20 amino acids as ancestors. In the example in Figure 4.11, the I, T, and V in column 1 could have arisen from any of the 20 amino acids by mutation. Suppose that they arose from an Ile (I). The profile values in the Ile (I) column of the corresponding row in the profile matrix would then use the amino acid scoring matrix values for I-I, I-T, and I-V, which are log odds scores of 5, 0, and 4 in the Dayhoff PAM250 matrix. Then the profile value for the I column is the frequency-weighted value, or 0.5 � 5 � 0.3 � 0 � 0.2 � 4 � 3.3. 

The profile table also includes penalties for matching a gap in the target sequence, shown in the two right columns. All of these table values are multiplied by a constant for convenience so that only the value of a score with one sequence relative to the score with another sequence matters. Once a profile table has been obtained, the table may be used in database searches for additional sequences with the same pattern (program Profilesearch) or as a scoring matrix for aligning sequences (program Profilegap). If several profiles characteristic of a protein family can be identified, the chance of a positive identification of additional family members is greatly increased (Bailey and Gribskov 1998; also see http://www.sdsc.edu/MEME). 

The evolutionary method for producing a profile table is based on the Dayhoff model of protein evolution (Chapter 2) (Gribskov and Veretnik 1996). The amino acids in each column of the msa are assumed to be evolving at a different rate, as reflected in the amount of amino acid variation that is observed. As with the average model, the object is to consider each of the 20 amino acids as a possible ancestor of the pattern of each column. In the evolutionary model, the evolutionary distance in PAM units that would be required to give the observed amino acid distribution in each column is determined. Recall that each PAM unit represents an overall probability of 1% change in a sequence position. For example, in the original Dayhoff PAM1 matrix for an evolutionary distance of 1 PAM unit (very roughly 10 my), the probability of an I not changing is 0.9872, and the probabilities for changing to a T or a V are 0.0011 and 0.0057, respectively. All of the probabilities of changing I to any other amino acid add up to 1.0000, for a combined probability of change of 1% for I. For an evolutionary distance of _n_ PAM, the PAM1 matrix is multiplied by itself _n_ times to give the expected changes at that distance. At a distance of 250 PAMs, the above three probabilities of an I not changing or of changing to a T or V are 0.10, 0.06, and 0.15, respectively, representing a much greater degree of change than for a shorter time, as might be expected (Dayhoff 1978). 

**164** I CHAPTER 4 

_Do not confuse these probabilities of one amino acid changing to another in the original Dayoff PAM250 matrix with scores from the log odds form of the PAM250 matrix, which have been used up to now. The log odds scores are derived from the original Dayhoff matrix by dividing each probability of change with the probability of a chance matching of the amino acids in a sequence alignment; i.e., that the one amino acid is not an ancestor of the other. These ratios are then converted to logarithms._ 

Thus, for the example of the msa column 1 with 5 Ile (I), 3 Thr (T), and 2 Val (V), the object is to find what amount of PAM distance from each of the 20 amino acids as possible ancestors will generate this much diversity. This amount can be found by a formula giving the amount of information (entropy) of the observed column variation given the expected variation in the evolutionary model, 



where _fa_ is the observed proportion of each amino acid _a_ in the msa column and _pa_ is the expected frequency of the amino acid when derived from a given ancestor amino acid. For a given column in the msa, _H_ is calculated for each 20 ancestor amino acids and for a large number of evolutionary distances (PAM1, PAM2, PAM4, . . . . ). The distance that gives the minimum value for _H_ for each column-possible ancestor combination is the best estimate of the distance that generates the column diversity from that ancestor. This analysis provides 20 possible models ( _Ma_ for _a_ � 1,2,3, . . . 20) as to how the amino acid frequencies in a column ( _F_ ) may have originated. The next step in the evolutionary profile construction determines the extent to which each _Ma_ predicts _F_ by the now-familiar Bayes conditional probability analysis. 



where the prior distribution _P_ ( _Ma_ ) is the given by the background amino acid frequencies and 



i.e., the product of the expected amino acid frequencies in _Ma_ raised to the power of the fraction observed for each amino acid in the msa column, as defined above. From _P_ ( _Ma_ � _F_ ), the weights for each of the 20 possible distributions that give rise to the msa column diversity are calculated as follows: 



where _Wa_ is the weight given to _Ma_ and _P_ ( _M_ random � _F_ ) is calculated as above using the background amino acid distribution. 

The log odds scores for the profile (Profile _ij_ ) are given by: 



where _Wai_ is the weight of an ancestral amino acid _a_ at row _i_ in the profile, _paij_ is the frequency of amino acid _j_ in the PAM amino acid distribution that best matches at row _i_ , and 

MULTIPLE SEQUENCE ALIGNMENT I **165** 

_p_ random _j_ is the background frequency of amino acid _j._ An example of a profile matrix for the ATP-dependent RNA helicase (“DEAD” box family) from the M. Gribskov laboratory is given in Figure 4.12. 

The usefulness of the evolutionary profile is demonstrated by the following: A profile for the 4Fe-4S ferredoxin family was prepared from six sequences. This profile was then used to search the SwissProt database for family members. Success was measured by the socalled receiver operating characteristic test (ROC) plot. The fraction of scores equal to or greater than a certain value is plotted for the true positive matches (a correct family member identified) on the _y_ axis and for the true negatives (unrelated sequences) on the _x_ axis. The area under the curve and the _x_ axis gives the probability of correct identification. The ROC50 is the area under the curve when it is truncated to the first 50 incorrect sequences, and can be used as a standard for success in a database search (Gribskov and Veretnik 1996). For the ferredoxin family search, the ROC50, 95.6 � 0.6% of the known family members, was identified in a search of SwissProt by an evolutionary profile, whereas 93.0 � 2.0% was identified by the average profile method (Gribskov and Veretnik 1996). The success rate was increased 0.4–0.6% by using 12 training sequences and 2–3% by using 134 training sequences. 

### **Block Analysis** 

Like profiles, blocks represent a conserved region in the msa. Blocks differ from profiles in lacking insert and delete positions in the sequences. Instead, every column includes only matches and mismatches. Like profiles, blocks may be made by searching for a section of an msa alignment that is highly conserved. However, aligned regions may also be found by searching each sequence in turn for similar patterns of the same length. These patterns may include a region with one or a few matching characters followed by a short spacer region of unmatched characters and then by another set of a few matching characters, and so on, until the sequences start to be different. These patterns are all of the same length, and when they are aligned, the matching sequence characters will appear in columns. The first alignments of this type were performed by computer programs that searched for patterns in sequences (Henikoff and Henikoff 1991; Neuwald and Green 1994). Several blocks located in different regions in a set of sequences may be used to produce a msa (Zhang et al. 1994), and blocks may be constructed from a set of aligned sequence pairs (Miller et al. 1994). Statistical and Bayesian statistical methods are also used to locate the most alike regions of sequences (Lawrence et al. 1993; Lawrence and Reilly 1990). Web sites that perform some of these types of analyses are discussed below and also given in Table 4.1. Finally, the information content of these tables can be displayed by a sequence logo (see p. 195). Note that few of these types of analyses presently provide a method for phylogenetic estimates of the sequence relationships so that sequence weighting can be used to make the changes more reflective of the phylogenetic histories among the sequences. Additionally, except where noted, these methods do not use substitution matrices such as the PAM and BLOSUM matrices to score matches. Rather, they are based on finding exact matches that have the same spacing in at least some of the input sequences, and that may be repeated in a given sequence. 

### **Extraction of Blocks from a Global or Local Multiple Sequence Alignment** 

A global msa of related protein sequences usually includes regions that have been aligned without gaps in any of the sequences. These ungapped patterns may be extracted from these aligned regions and used to produce blocks. Blocks found in this manner are 

**166** I CHAPTER 4 



<!-- Start of picture text -->
A. The multiple sequence alignment.<br>ghle_ecoli GVDVLVATPG RLLDLEWQWA ....VKLOQV EILVLDUADR MLDMGPINDI<br>dbp2_schpo GVEILCIATPG NLLOMLDSNE ....THLERVY TYLVLDEADR MLDMGFEPQI<br>dbp2_yeast GSEIVIATIO MLIDMLEIGR ....THLKRV TYLVLORADR MLDMGFEPOL<br>dbpa_ecoli APHIIVATPO ALLDALQKGT ....VSLDAL WILVMDEADR MLOMGPSDAL<br>rm62_dreme GCEIVIATPG RLIDFLSAGS ....THLKAC TYLVLOEADR MLDMGFEPOT<br>pob_human GVEICIATPG MLIDFLECGE ....THLERT TYLVLDEADR MLDMGFEPQL<br>Fhlb_ecoli GVDILIGITS MLIDYAKQNH ....INLGAT QVVVLDEADR MYDLGFIKDT<br>yr2l_cage] RPHIIVATSS ALVDHLENTE ...GPNLKAL KPLINDEADR ILMMDFEVEL<br>yhmS_yeast KPHIIIATPO RLMDALENTE ...GPSLRKL KFLVMDEADR LLOMEFGPVL<br>me31_crome KVQLIIATPG RILDLMDERV ....ADMSHC RILVLDEADK LLSLDFOGHL<br>drei_yeast RPDIVIATPG RFIDHIRNSA ...SFNVDSV EILVMDEADR MLEEGFQDEL<br>if‘a_rabit APHTIVGIPG AVFOMLNRRY ....LGPRYI EMFVLOEADE MLSRGFKDQI<br>if41_human APHTIVGTPG AVFDMLNRRY ....LSPKYI KMFVLDEADE MLSRGFKDOI<br>vasa_drome GCHVVIATPO RLLDFVDRIF ....ITFEDT REVVLDEADR MLOMGFSEDM<br>srmb_ecoli WQODIVVATTS ALLOYIKEEN ....FDCRAV ETLILDRADR MLOMGFAQDT<br>DEAD _ecoli GPQIVVGTPS RLLDHLERGT ..+-LDLSEL SGLVLORADE MLRNGFIEDY<br>if@a_oryaa GVIVVVGTPG AVFDMLRRQS ...-LRPDYI KMFVLDEADE MLSRGFKDQI<br>DEAD kiepn GPQIVVGTPG HLLDELKRGT ....LDLSKL SGLVLDEADE MLEMGFIEDV<br>pll0_mouse GCHLLVATPO HLVDOMMERGK ....IGLOFC KYLVLDEADR MLDMGFEPOI<br>pod_human TVHVVIATPO RILDLIKKGY ....ARVDHV QMIVLDEADK LLSQDFVQIM<br>if4a_drome GCHVVVGTPG RVYDMIHREL .....RTQ¥L KLFVLOEADE MLSRGFRKOQT<br><!-- End of picture text -->



<!-- Start of picture text -->
dedl_yeast GCDLLVATPG ALNDLLERGK ....ISLANV KYLVLOEADR MLDMGFEPQT<br>msl6_yeast RPHIVIATPG BLIDVLEEYS ...NKFFRFY DYEVLOEADR LLEIGFRODL<br>pr2a_yeast GCDILVATPG BLIDSLENHL ....LVMEQV ETLVLDEADK MYDLGFEDQVY<br>if4n_human GQHVVAGTPS RVPDMIARRS ...-LETRAT EMLVLOEADE MLNKGPREQI<br>ani_xenla GCHLLVATPG HLVOMMERGE ....IGLDFC KYLVLDEADR MUDMGFEPOQI<br>dbpl_yeast GCDLLVATPG RLNDLLERGE ...,VSLANI KYLVLDEADR MLDMGFEPOI<br>if4a_yeast DAQIVVGTPG EVFDNIQRRR ....FRTDKI FMFILDEADE MLSSGFREQI<br>Spb4_yeast RPQILIGTPS SVLDPLOMPA ....VKTSAC SMVVMDEADR LLDMSPIKDT<br>iff~a_caeel GIHVVVGIPG NVGDMINKNA ....LDTSRI KMFVLDEADE MUSRGFEDQI<br>Pros yeast GTEIVVATPC RFIDILTLND .GKLLSTKRI TFVVMDEADR LEDLGFEPQI<br>if42_mouse APHIVVGTPG EVFDMLNERY ....LSPKWI KMFVLDEADE MLSRGFRDQI<br>dhhl_yeast TVATLVGTPO RVLOLASRRV ....ADLSDC SLFIMDEADK MLSRDFRTIT<br>d@b73_drome KADIVVTTPG RLVDHLHATK ...GFCLESL KFLVIDEADR IMDAVPQNWL<br>yk04_yeast GCNFIIGIPG EVLDHLQNTR VIKEQLSQSL RYIVLDEGDR LMELGFDETI<br>ybz2_yeast SGQTVIATPG RFLELLEKDN .TLIKRPSEV WILILDEADR LLAQOGHFDEF<br>yhw9_yeast KPHFILATPO RLAHHIMSSG DDTVGGLMRA KYLVLDEADI LUTSTPADEL<br>gihl_cacel GATIIVGTVG RIKMFCEEGT ....IKLDEC RPFVLOEADR MIDAMGFGTD<br><!-- End of picture text -->

**Figure 4.12.** msa and the derived evolutionary profile. 



<!-- Start of picture text -->
B. The evolutionary profile. Note the location of red conserved regions in the alignment in the corresponding profileof these<br>sequences.<br>CcmsA EB C OB FE F G@ BB Tt KF FT © Ww P @ RB 8 T VW W F £ Gap ton<br>a 17) 18 «60 628 ta 8-1 3B 11 100 100<br>er ee se a i a a a eS a a eT ee ee ECL)<br># 5 24 -12 29 25 -20 @ 32 -9 & -10 -9 2 7 30 19 0 4& -8 -20 -7 27 100 100<br>Ioowl -120 6 <1 911 93 -12 12 63 -11 40 29 15 -3 -14 -15 -6 7 50 -17 B -11 100 100<br><!-- End of picture text -->



<!-- Start of picture text -->
vod -l 1 =I -9 (22 -3 -11 46 -9 37 3D -190 -0) -8 -130 -§ 066 650-19 2-H oo<br>v 8 =9 9 = 28 19 22 213 87 29 38 26 =19° <2 s1L =12 =4 0 69 58 92900 0) 29 100. 100<br>A 54 15 12 20 17 -2@ 44 -6 -4 -1 -1]1 -5 12 19 # -13 21 19 «9 -39 -20 10 100 10D<br>T 40 20 » 26 2 6-30 40-16 zo 25 «6-10 6 2 mM -l6 -16 ao 6186 20-60 «=30 io 100 «100<br>Pu & #F & & wil 2s 12 -8 € -16 -11 8 iT 17 2 8-5 48) 12 100 100<br>a aa oo 20 TH 50 =<60 380 =-20 -30 =-10 -50 <30 ag Ba mam =< 60 40 20-100 =70 »* 106 «Loo<br>14a<br>BR =30 19 =-30 o o -50 -30 so -20 eo -40 20 19 Mm a 159 19 -10 = 140 -650 20 100 #190<br>Lo o-2 -17 -18 -18 -12 22 -13 -9 38 -12 49 39 -15 <9 -9 -15 -11 6 38 & 42 -10 190 100<br>L OQ -12 -15 -14 -5 38 -12 -7 320 +7 8410 38 -12 +5 +6 He +8 o 2 6 & <7 100 100<br>D 15 58 -27 7H S54 -52 35 27 -12 46 -26 -21 38 6 41 3 9 10 -12 -57 -25 50 100 100<br>bo o-5 -<5 =? -8 =-¢ 24 -12 43 13 =-6 #25 17 -2 -7 ® -2 -8 -3 10 21 47 =-2 100 1090<br>bo o3 -19 -19 -13 -@ 31 -Il -8@ 34 -8 41 36 -12 -7 -§ -13 -§@ 2 2 -l @ -6 LOO 100<br>& 6 19 =15 23 27 -21 9 48 =8 48 <8 -1 16 6 22 12 6 5S -=6 -15 -16 25 190 100<br>K 3 Wd -12 2 12 -l6 2 40 -5 23 -7 # 15 6 18 2 8 F -5 7 -15 14 100 100<br>6 ll 27 @ 16 44 -16 29 5 =-6§ AE -11 -5 16 9 & & 4 25 -1 -12 -14 11 100 100<br>7 12 a -1 7 T <a a 2 4 a o a 10 5 4 3 5 42 7 <a <8 S 100 100<br><!-- End of picture text -->

**Figure 4.12.** _Continued._ 



<!-- Start of picture text -->
1a.<br>D t 2 ®© &@ ¥F = 4£ @ 12 ooo » 1 0 1 o oO ae ee oe | 1 92 22<br>T 2 © © 3 F -2 2 B9 F B09 GO O09 2@d4a =a er ee ee oe a ae ae|)<br>Ls a 1 <3 o 1 a o a 1 4 a 3 1 o 1 1 0 i 1 0 2 1 22 22<br>G a 3 o 4 4 -1 & -1 a o 1 1 3 1 1 2 4 a 5 -6 -3 2 a2 22<br>L Ss 6 -4 7 -4 16 2 4 21 «4 623 a7 = -4 -4 8 =-2 4 a9 o € -4 220 22<br>o 5 46 - GU -lb Ww 6 -2 WwW -8 -l i 4 * 1 7 -@ -3 -il 1g 109 «6100<br>L Lo =i} -12 -14 0-9 27 = =T Oe BG O20 57 HD = 7 2 6 9 --f 100 190<br>o 7 3 6-7 a 17-22 7 tT - WU -ll <3 4 a ub 4 7 6 <5 <-5 <168 16 100 100<br>KF ll WwW -3 1 9 -12 5 9 -4 #16 -6 o 19 6 ul 2 iW 4 -4 a -@ «610 200 6100<br><!-- End of picture text -->



<!-- Start of picture text -->
v 7 =L0 11 <11 -16 44 4 8 41-11 «619 16-10 Oo =-10 ~12 z @ 4 -22 9 -10 100 100<br>1 31<br>x 4 16-11 a6 17 -21 6 9 7 a7 -12 [) 16 7 16 20 iz a -7 o -18 17 100 «100<br>F 3 +10 a -14 -10 2o +10 =3 ao -9 26 Pt} <7 =-12 -14 -10 <6 4s 14 10 20 4-11 100 100<br>L -S5 -18 -19 -Z1 -14 az -16 -8 34-12 51 40 -15 -12 -10 -15 -It -3 aa 16 17-11 «100 «100<br>¥ 410 -1l1 tl -Li ~-31 15 8-17 66 -ll 46 a4 -17 4-12-17 -6 a 83-43 -4 -11 too too<br>BG =f -23 -37 <24 =14 53 =23 =11 9-1) 679 GT =1lS =14 “4 -l4 -19 (-4 aa iF Ik -# 100 100<br>oD 36 110 =80 186 106 -100 7000 «46-960 30 Bd 74 6 @& 20 -<25-116 -56 8 100 106<br>Se 36 WW -60 100 150 -70 SO 40 -20 Mf -30 -20 50 10 70 0 3 20 «20-110 -50 110 100 100<br>A 121 #17 #@ 28 @5 -41 Se -8 o o -3 6 WW 4b 16 -4 32 33 18-86-25) 16 100 too<br>Dp 30 110 -50 160 100-100 70 40 -20 30 -50 -—40 7 10) 70 0 2 20 -20-110 -50 #0 100 100<br>R =7 9 =16 7 10-2] =T 17) =f) 66 dd T a 10 18 49 5 0 -10 a5 25. 14 100 «100<br><!-- End of picture text -->

**Figure 4.12.** _Continued_ . 



<!-- Start of picture text -->
b 42<br>Woo-t -16 -28 -29 -19 49 -16 -22 31 1 60 62 -18 -10 -2 -1 -14 -1 31 -2 2 -6 100 100<br>L -6 -25 -35 -26 -16 61 -26 -10 41 -15 74 64 -20 -17 -7 -20 -20 -5 33 24 16 -12 i604 166<br>Bb 12 32 3 41 a2 -30 23 1 7 120-17) -12 a3 7 279 3 ig 1a =-7 =-27 -19 27) «100 )«(100<br>4 5 “2 = -15 -3 L a -3 -3 13 6 25 30 =i 0 5 6 -1 2 13 -1i “6 3 100 100<br>a 26 26 4 a1 23. 25. 54 =5 ~8 “1 -16 -3 17 41 10-10 23 17 6 -41 -27 15 100 «(100<br>F =«31 -44 =-6 <«63 -44 $6 -36 4 44-45 7 SL o=31 =$45 =50 =32 =-19 =19 13 Bz 90 -44 100 100<br>E 7 4a -7 ko 13 -7 5 4 7 10 a 5 6 2 10 a 6 5 5 => -7 12 100 «100<br>Bb lz 24-12 31 ay -27 as lz 6 lz -12 -8 17 1a 21 4 & a =-3 =30 -15 25 #100 «#106<br>DB a 18 -14 2 22 -18 il 13 o b -3 -1 13 5 24 5 3 a O -24 -12 74 «loo 100<br>Ir o -8 5 -8 -7 27 > «-10 => ai =7 36 29 -10 8 -? -li -6 6 35 —# 5 -6 #100 #100<br>1 $1<br>EB & 15 =5 16 21-11 7 it =2 a =-7 =-3 12 2 12 6 6 6 -3 8 -7 16 100 100<br># 3 47-4 mo al -19) 5 ay -7 18 -9 -3 15 5 2 2s 4 @& =7 =10 -12 21 100 100<br>I 2 -9 9 -19 -9 27 -9 -13 60 -9 34 @ -13 -€ -13 -13 -2 42 49 -23 4 -9 100 100<br><!-- End of picture text -->



<!-- Start of picture text -->
L 1 -10 -5 -11 -8 26 -9 -8 33 -7 34 28 -10 -7 -8 -9 5 2 28 -2 8 -7 100 100<br>K 8 9 -4 9 9 -13 11 1 0 16 -4 4 8 7 8 il 13 12 3 “2 -15 8 100 100<br>L 3 4 -9 3 6 3 -2 8 9 7 10 10 5 0 8 3 0 5 7 -2 0 7 100 100<br>L 1-13 -13 -13 -9 32-11 -7 32 -9 42 36-12 -7 -6 -13 -9 3 33 2 8 -7 100 100<br>* 99 tt) 25 208 120 94 137 44 181 105 256 94 41 62 64 144 59 99 162 3 35 CY)<br><!-- End of picture text -->

|**Figure 4.12.**_Continued._|
|---|



**170** I CHAPTER 4 

only as good as the msa from which they are derived. Using the BLOCKS (http://www.blocks.fhcrc.org/blocks/process_blocks.html), blocks of width 10–55 are extracted from a protein msa of up to 400 sequences (Henikoff and Henikoff 1991, 1992). The program accepts FASTA, CLUSTAL, or MSF formats, or manually reformatted msas. Several types of analyses may be performed with such extracted blocks. The BLOCKS server primarily generates blocks from unaligned sequences. The eMOTIFs server at http://dna.stanford.edu/emotif/ (Nevill-Manning et al. 1998) similarly extracts motifs from msas in several msa formats and provides a formatter for additional msa formats. These types of analyses are discussed below in greater detail. 

### **Pattern Searching** 

This type of analysis was performed on groups of related proteins, and the amino acid patterns that were located may be found in the Prosite catalog (Bairoch 1991). This catalog groups proteins that have similar biochemical functions on the basis of amino acid patterns such as those in the active site. Subsequently, these families were searched for amino acid patterns by the MOTIF program (Smith et al. 1990), which finds patterns of the type aa1 d1 aa2 d2 aa3, where aa1 and aa2 are conserved amino acids and d1 and d2 are stretches of intervening sequence up to 24 amino acids long. These initial patterns are then organized into blocks between 3 and 60 amino acids long by the Henikoff PROTOMAT program (Henikoff and Henikoff 1991, 1992). The BLOCKS database can be accessed at http://www.blocks.fhcrc.org/, and the server may also be used to produce new blocks by the original pattern-finding method or other methods described below. 

Although used successfully for making the BLOCKS database, the MOTIF program is limited in the pattern sizes that can be found. The MOTIF program distinguishes true motifs from random background patterns by requiring that motifs occur in a number of the input sequences and tend not to be internally repeated in any one sequence. As the length of the motif increases, there are many possible combinations of patterns of a given length where only a few characters match, e.g., �10<sup>9</sup> possible patterns for a 15-amino-acidlong pattern with only five matches. The MOTIF program always provides a motif, even for random sequences, thus making it difficult to decide how significant the found motif really is. This problem has been circumvented by combining the analysis performed by MOTIF with that of the Gibbs sampler (discussed on p. 177), which is based on sound statistical principles. A rigorous searching algorithm called Aligned Segment Statistical Evaluation Tool (ASSET) has been devised (Neuwald and Green 1994) that can find patterns in sequence up to 50 amino acids long, group them, and provide a measure of the statistical significance of the patterns. These patterns may also include certain pairs, the 26 positive scoring pairs in the BLOSUM62 scoring matrix. Consideration of all BLOSUM pairs is not possible because this would greatly increase the complexity of the analysis. 

The efficiency of ASSET is achieved by a combination of an efficient pattern search strategy called the depth-first method, which assures searching for the same patterns only once, and the use of formulas for efficiently organizing the patterns. Low-complexity regions with high proportions of the same residue and use of sequences, some of which are more similar than the others, can interfere with the ability of the method to find a range of patterns. ASSET removes low-complexity regions and redundant sequences from consideration. The program was easily able to find subtle motifs in the DNA methylase, reverse transcriptase, and tRNA ligase families, and previously identified by the MOTIF program. In addition, however, ASSET gave these motifs an expect score, the probability that these are random matches of unrelated sequences, of �0.001. The program also found motifs in 

MULTIPLE SEQUENCE ALIGNMENT I **171** 

families with only a fraction of the sequences sharing a motif (the acyltransferase family) and in a set of distantly related sequences sharing the helix-turn-helix motif. Finally, the program found several repeat sequences in a prenyltransferase and ankyrin-like repeats in an _E. coli_ protein. This source code of the program is available by anonymous FTP from ncbi.nlm.nih.gov/pub/neuwald/asset. The European Bioinformatics Institute has a Web page for another complex pattern-finding program (PRATT) at http://www2.ebi.ac. uk/pratt/ (Jonassen et al. 1995). 

### **Blocks Produced by the BLOCKS Server from Unaligned Sequences** 

As described above, the BLOCKS server can extract a conserved, ungapped region from a msa to produce a sequence block. This same server can also find blocks in a set of unaligned, input sequences and maintains a large database of blocks based on an analysis of proteins in the Prosite catalog. Blocks are found by the Protomat program (Henikoff and Henikoff 1991). Blocks are found in two steps: First, the program MOTIF (Smith et al. 1990) described on the previous page is used to locate spaced patterns. The second step takes the best and most consistent patterns found in step 1 and uses the program MOTOMAT to merge overlapping triplets and extend them, orders the resulting blocks, and chooses those that are in the largest subset of sequences. Since 1993, the Gibbs sampler (see below) has been used as an additional tool for finding the initial set of short patterns also by specifying that the sampler search for short motifs. This program is based on a statistical analysis of the sequences and can identify the most significant common patterns in a set of sequences. 

An example of BlockMaker output using an example from Lawrence et al. (1993) is shown below. The program first searches for blocks using either the MOTIFS or Gibbs sampler program to identify patterns, then the Protomat program to consolidate the patterns into meaningful blocks. The results of both types of analyses are reported. 



<!-- Start of picture text -->
A. Motif analysis<br>LipocalA, width = 15 LipocalB, width = 11<br>BBP_PIEBR 16 NFDWSNYHGKWWEVA ( 70) 101 VLSTDNKNYII<br>ICYA_MANSE 17 DFDLSAFAGAWHEIA ( 73) 105 VLATDYKNYAI<br>LACB_BOVIN 25 GLDIQKVAGTWYSLA ( 70) 110 VLDTDYKKYLL<br>MUP2_MOUSE 27 NFNVEKINGEWHTII ( 101) 143 DLSSDIKERFA<br>RETB_BOVIN 14 NFDKARFAGTWYAMA ( 77) 106 IIDTDYETFAV<br>B. Gibbs sampler analysis<br>LipocalA, width = 15 LipocalB, width = 11<br>BBP_PIEBR 16 NFDWSNYHGKWWEVA ( 70) 101 VLSTDNKNYIT<br>ICYA_MANSE 17 DFDLSAFAGAWHEIA ( 73) 105 VLATDYKNYAI<br>LACB_BOVIN 25 GLDIQKVAGTWYSLA ( 70) 110 VLDTDYKKYLL<br>MUP2_MOUSE 27 NFNVEKINGEWHTII ( 68) 110 IPKTDYDNFLM<br>RETB_BOVIN 14 NFDKARFAGTWYAMA ( 77) 106 IIDTDYETFAV<br><!-- End of picture text -->

**172** I CHAPTER 4 

In the above example, two blocks identified as Lipocal A and B are reported using both the MOTIF and Gibbs sampler programs for step 1, the initial pattern-finding step. The MOTIF program is based on a heuristic method that will always find motifs, even in random sequences, whereas the Gibbs sampler discriminates found motifs based on sound statistical methods. These blocks are identical to those determined from analysis of threedimensional structures. Note that MOTIF aligned MUP2_MOUSE incorrectly in the B 



<!-- Start of picture text -->
@ MFRRKAFLHWYTGEGMDEMEF TEAESNMNDPVAEYQQY<br>MFKRKAFLHWYTGEGMDEMEF TEAESNMNDLVSEYQQY<br>MF KRKAF LHWY TGEGMDEMEF TEVRANMNDLVAEYOQQY<br>MFKRKAF LHWY TSEGMDELEF SEAESNMNDLVSEYQQY<br>MFKRKGFMFRRKAF LHWYTGEGMEPVEFSEAQSDLEDL LHWF TGEGMDEMEF TEAESNMNDLVSEYQQY | LEYQQY<br>MFRRKAF LHWY TGEGMDEMEF SEAEGNTNDLVSEYQQY<br>MFRRKAF LHWYTGEGMDEMEF TEAESNMNDLMSEYQQY<br>MFRRKAF LHWY TGEGMDEMEF TEAESNMNDL VAEYQQY<br>MFRRKAFLHWYTGEGMDEMEF TEAESNMNDLVHEYQQY<br>MFRRKAFLHWYTGEGMDEMEF TEAESNMNDLVSEYQQY<br>MFRRKAFLHWY<br>TGEGMDEMEF TEAESNMNELVSEYQQY<br>MFRRKAFLHWYTLEGMEELEF TEAESNMNDLVYEYQQY<br>MFRRKAFLHWY<br>TNEGMD I TEFAEAESNMNDLVSEYQQY<br>MFRRKAFLHWYTSEGMDEMEF TEAESNMNDLVSEYQQY<br>MFRRKRFLHWY<br>MFRRNAF LHWYTGEGMDEMEFTGEGMDEMEF TEAESNMNDLVSEYQQYTEAESNMNDLVSEYQQY<br>MFRRQAFLHWYTSEGMDEMEF TEAESNMNDLVSEYQQY<br>MFSRKAFLHWYTGEGMEEGDFAEADNNVSDLLSEYQQY<br>MFGKRAFVHHYVGEGMEENEF TDARQDLYELEVDYANL.<br>MFKKRAFMFVKRAF VHWYVGEGMEEGEFVHWYVGEGME EGEFAT EAR ENDDLLALEKDYESVI AVLERDFEEV<br>MYAKRAF VHWYVGEGMEEGEFAEAREDLAALEKDYEEV<br>MYAKRAFMYAKRAF VHWYVGEGMEEGEFS VHWYVGEGMEEGEFS EAR EDEDLAALEKDFEEV I AALEKDYEEV<br>MYAKRAF VHWYVGEGMEEGEFSEAREDLAALEKDYEEV<br>MYAKRAF VHWYVGEGMEEE<br>MYAKRAF VHWYVGEGMEGE F SEAREDMAALEKDY E EV<br>GE F SEVREDL.AAL EKDYEEV<br>MYAKRAFVHWYVGEGMEE<br>MYAKRAF VHWYVGEGMEGE F TEAREDLAALEKDY E EV<br>MYAKRAF VHWYVGEGMEEVEFSEAREDGE F TEARED LAAL ERDYEKDYEEV| EV<br>MYAKRAF VHWYVSEGMEEGEFAEAREDLAALEKDYDEV<br>MYSKRAF VHWYVGEGMEEGEFSEAREDLAALEKDYEEV<br>MYSKRAFVHWYVGEGMEEGEFSEAREDLAALERDYEEV<br>b MF.K. .FVH.F..EGMQ. .QFPQ...Q......QF.<br>YR LY N NAN N NY<br>Ww ow —E EGE E EW<br>D OSD D D<br>T<br>C MF.KR.FLHWFT.EGMQ. .QFPE...Q. .DLI.DYQQY<br>R Y N NA N L<br>Ww EEG E M<br>DoS dD Vv<br>T<br><!-- End of picture text -->

**Figure 4.13.** Aligned block of 34 tubulin proteins. ( _a_ ) The sequences are divided into two groups based on the occurrence of R or L in the fourth position and Y in the last position. ( _b_ ) Specific substitution groups found in the columns of the block. If a group cannot be found, then the position is ambiguous and a dot is printed at the position. ( _c_ ) If only the first group of sequences is used, a more specific motif may be found because sequences in this group are more closely related to each other. (Reprinted, with permission, from Nevill-Manning et al. 1998 [copyright National Academy of Sciences].) 

MULTIPLE SEQUENCE ALIGNMENT I **173** 

block. The Gibbs sampler results may differ when the same sequences are submitted repeatedly with a different initial alignment (see below). 

### **The eMOTIF Method of Motif Analysis** 

Another somewhat different but extemely useful method of identifying motifs in protein sequences has been described (Nevill-Manning et al. 1998). From the BLOCKS database (derived from msa of proteins in the Prosite catalog) and the HSSP database (derived from msa of proteins based on predicted structural similarities), a set of amino acid substitution groups characteristic of each column in all of the alignments was found. These patterns reflect the higher log odds scores in the amino acid substitution matrices. A statistical analysis was performed to identify amino acids that are found together in the same msa column as opposed to amino acids that are found in different columns at the 0.01 level of significance. Thirty and 51 substitution groups that met this criterion were found in the BLOCKS and HSSP msas, respectively. For example, the chemically aromatic group of amino acids F, W, and Y were found to define a group often located in the same column of the msa. 

From the msa for a particular group of proteins, each column is examined to see whether these groups are represented in the column, as illustrated in Figure 4.13. In column 1, M is always present, and because M is one group, M is used in column 1 of the motif, as shown in part _b_ . Similarly for column 2, Y and F, which are members of the group FYW, are found, and hence this group is used as column 2 in the motif. The final motif shown in _b_ describes the variation in all the sequences. Instead, a motif may be made for only the first group of 19 sequences, and is shown in _c_ . This second motif ( _c_ ) has less variability and greater specificity for the first 19 sequences and thus would be more likely to find those sequences in a database search (i.e., it is a more sensitive motif for those sequences) than motif _b_ . 

The probability of each motif is estimated from the frequencies of the individual amino acids in the SwissProt database. The probability of the motif _b_ above is given by the product of the probability sums in each column, or _p_ (Motif) � _p_ ( _M_ ) � 1 � [ _p_ ( _F_ )� _p_ ( _W_ )� _p_ ( _y_ )] � [ _p_ ( _Y_ )� _p_ ( _R_ )] x . . . This value has been found to provide a good estimate of false positives, or of the selectivity of the motif, in a database search. Both the sensitivity and selectivity of a given motif must be taken into account in using the motif for a database search. Ideally, a motif can find all of the sequences used to generate the motif but none other. In practice, eMOTIF produces a large set of motifs, some more and some less sensitive for the set of aligned sequences. The more sensitive ones, which are also the most selective based on the value of _p_ (Motif), are then chosen. Some are useful for specifying subfamilies of a protein superfamily. A database of such motifs called Identify is a useful resource for discovering the function of a gene (Nevill-Manning et al. 1998; http://dna.stanford.edu/emotif/). 

### **~~STATISTICAL METHODS FOR AIDING ALIGNMENT~~** 

### **Expectation Maximization Algorithm** 

This algorithm has been used to identify both conserved domains in unaligned proteins and protein-binding sites in unaligned DNA sequences (Lawrence and Reilly 1990), including sites that may include gaps (Cardon and Stormo 1992). Given are a set of sequences that are expected to have a common sequence pattern and may not be easily recognizable by eye. An initial guess is made as to the location and size of the site of interest 

**174** I CHAPTER 4 

in each of the sequences, and these parts of the sequence are aligned. The alignment provides an estimate of the base or amino acid composition of each column in the site. The EM algorithm then consists of two steps, which are repeated consecutively. In step 1, the expectation step, the column-by-column composition of the site already available is used to estimate the probability of finding the site at any position in each of the sequences. These probabilities are used in turn to provide new information as to the expected base or amino acid distribution for each column in the site. In step 2, the maximization step, the new counts of bases or amino acids for each position in the site found in step 1 are substituted for the previous set. Step 1 is then repeated using these new counts. The cycle is repeated until the algorithm converges on a solution and does not change with further cycles. At that time, the best location of the site in each sequence and the best estimate of the residue composition of each column in the site will be available. 

As an example, suppose that there are 10 DNA sequences having very little similarity with each other, each about 100 nucleotides long and thought to contain a binding site near the middle 20 residues, based on biochemical and genetic evidence. As we will later see when examining the EM program MEME, the size and number of binding sites, the location in each sequence, and whether or not the site is present in each sequence do not necessarily have to be known. For the present example, the following steps would be used by the EM algorithm to find the most probable location of the binding sites in each of the 10 sequences. 

##### **The Initial Setup of the Algorithm** 

The 20-residue-long binding motif patterns in each sequence are aligned as an initial guess of the motif. The base composition of each column in the aligned patterns is then determined. The composition of the flanking sequence on each side of the site provides the surrounding base or amino acid composition for comparison, as illustrated below. For illustration purposes, each sequence is assumed to be the same length and to be aligned by the ends, and each character in the alignment represents five sequence positions (o, not in motif; x, in motif). 

ooooooooxxxxoooooooo ooooooooxxxxoooooooo ooooooooxxxxoooooooo ooooooooxxxxoooooooo ooooooooxxxxoooooooo ooooooooxxxxoooooooo ooooooooxxxxoooooooo ooooooooxxxxoooooooo ooooooooxxxxoooooooo ooooooooxxxxoooooooo 

**Columns defined by a preliminary alignment of the sequences provide initial estimates of frequencies of amino acids in each motif column** 

**Columns not in motif provide background frequencies** 

MULTIPLE SEQUENCE ALIGNMENT I **175** 

The number of each base in each column is determined and then converted to fractions. Suppose, for example, that there are four Gs in the first column of the 10 sequences, then the frequency of G in the first column of the site, _fs_ G � 4/10 � 0.4. This procedure is repeated for each base and each column. For the rest of the sequences not included in the sites, the background frequency of each base is calculated. For example, let one of these four values for the background frequency, the frequency of G, be _fb_ G � 224/800 � 0.28. These values are now placed in a 5 � 20 matrix of values, the first column for the background frequencies, and the next 20 columns for the base frequencies in each successive column in the sites. Thus, the counts in the first three columns of the matrix may appear as shown in Table 4.2. 

The following calculations are performed in the expectation step of the EM algorithm: 

1. The above estimates provide an initial estimate of the composition of the site and the location in each sequence. The object of this step is to improve this estimate by discriminating to the greatest possible extent between sequence within and sequence not within the site. Using the above estimates of base frequencies for (1) background sequences that are not within the site and (2) each column within the site, each sequence is scanned for all possible locations for the site to find the most probable location of the site. For the 10-residue DNA sequence example, there are 100 �20 � 1 possible starting sites for a 20-residue-long site, the first one being at position 1 in the sequence ending at 20 and the last beginnning at position 81 and ending at 100 (there is not enough sequence for a 20-residue-long site beyond position 81). 

## Sequence 1 xxxxoooooooooooooooo xxxx 



<!-- Start of picture text -->
A<br>oxxxxooooooooooooooo<br>xxxx<br>B<br>ooxxxxoooooooooooooo<br>xxxx<br>C<br>Use previous<br>...background<br>estimates of amino<br>frequencies in the<br>acid frequencies for<br>remaining positions.<br>each column in the<br>motif to calculate<br>probability of motif in<br>this position, and<br>multiply by...<br><!-- End of picture text -->

The resulting score gives the likelihood that the motif matches positions (a) 1-20, (b) 6-25, or (c) 11-30 in sequence 1. Repeat for all other positions and find most likely locator. Then repeat for the remaining sequences. 

**176** I CHAPTER 4 

**Table 4.2.** _Column frequencies of each base in the example given_ 

||**Background**|**Site column 1**|**Site column 2**|**. . .**|
|---|---|---|---|---|
|G|0.27|0.4|0.1|. . .|
|C|0.25|0.4|0.1|. . .|
|A|0.25|0.2|0.1|. . .|
|T|0.23|0.2|0.7|. . .|
||1.00|1.0|1.0||



The first column gives the background frequencies in the flanking sequence. Subsequent columns give base frequencies within the site given in the above example. 

For each possible site location, the probability that the site starts is just the product of the probabilities given by Table 4.2. For example, suppose that the site starts in column 1 and that the first two positions in sequence 1 are A and T, respectively. The site will then end at position 20 and the first two nonsites, flanking background sequence positions, are 21 and 22. Suppose that these positions have an A and a T, respectively. Then the probability of this location of the site in sequence 1 is given by _P_ site1,sequence1 � 0.2 (for A in position 1) � 0.7 (for T in position 2) � _P_ s for next 18 positions in site � 0.25 (for A in first flanking position) � 0.23 (for T in second flanking position) � _P_ s for next 78 flanking positions. Similar probabilities for _P_ site2, sequence1 to _P_ site78, sequence1 are then calculated, thus providing a comparative set of probabilities for the site location. The probability of this best location in sequence 1, say at site _k_ , is the ratio of the site probability at _k_ divided by the sum of all the other site probabilities _P_ (site _k_ in sequence 1) � _P_ site _k_ , sequence 1 / ( _P_ site 1, sequence 1 � _P_ site 2, sequence 1 � . . . . . � _P_ site 78, sequence 1). The probability of the site location in each sequence is then calculated in this manner. 

2. The above site probabilities for each sequence are then used to provide a new table of expected values for base counts for each of the site positions using the site probabilities as weights. For example, suppose that _P_ (site 1 in sequence 1) � 0.01 and that _P_ (site 2 in sequence 1) � 0.02. In the above example, the first base in site 1 is an A and the first base for site 2 is a T. Then 0.01 As and 0.02 Ts are added to the accumulated list of bases at site column 1. This procedure is repeated for every other 76 possible first columns in sequence 1. Similarly, site column 2 in the new table of expected values is augmented by counts from the 78 possible column 2 positions in sequence 1, the first, for example, being 0.01 Ts. The weighted sequence data from the remaining sequences are also added to the new table, resulting finally in a new estimate of the expected number of each base at each site position and providing a new version of Table 4.2. 

In this maximization step, the base frequencies found in the expectation step are used as an updated estimate of the site residue composition. In this case, the data are more complete than the initial estimate because all possible sites in each of the sequences have been evaluated. The expectation and maximization steps are repeated until the estimates of the base frequencies do not change. 

##### **An Alternative Method of Calculating Site Probabilities by the EM Algorithm** 

The example shown above uses the frequencies of each base in the trial alignment and background base frequencies to calculate the probabilities of each possible location in each sequence. An alternative method is to produce an odds scoring matrix calculated 

MULTIPLE SEQUENCE ALIGNMENT I **177** 

by dividing each base frequency by the background frequency of that base. The probability of each location is then found by multiplying the odds scores from each column. An even simpler method is to use log odds scores in the matrix. The column scores are then simply added. In this case, the log odds scores must be converted to odds scores before position probabilities are calculated. 

### **Multiple EM for Motif Elicitation (MEME)** 

A Web resource for performing local msas by the above expectation maximization method is the program Multiple EM for Motif Elicitation (MEME) developed at the University of California at San Diego Supercomputing Center. The Web page for two versions of MEME, ParaMEME, a Web program that searches for blocks by an EM algorithm (described below), and a similar program MetaMEME (which searches for profiles using HMMs, described below) is found at http://www.sdsc.edu/MEME/meme/website/ meme.html. The Motif Alignment and Search Tool (MAST) for searching through databases for matches to motifs may also be found at http://www.sdsc.edu/MEME/ meme/website/mast.html. 

MEME will locate one or more ungapped patterns in a single DNA or protein sequence or in a series of DNA or protein sequences. A search is conducted for a range of possible motif widths, and the most likely width for each profile is chosen on the basis of the loglikelihood score after one iteration of the EM algorithm. The EM algorithm then iterates to find the best EM estimate for that width. Three types of possible motif models may be chosen. The OOPS model is for one expected occurrence of a motif per sequence, the ZOOPS model is for zero or one occurrence per sequence, and the TCM model is for a motif to appear any number of times in a sequence. These models are reflected in the choices on the Web page (Fig. 4.14). The current version of MEME can use prior knowledge about a motif being present in all or only some of the sequences, the length of the motif and whether it is a palindrome (DNA sequences), and the expected patterns in individual motif positions (Dirichlet mixtures, see section on HMMs, p. 189) that provide information as to which amino acids are likely to be interchangeable in a motif (Bailey and Elkan 1995). Once a motif has been found, the motif and its position are effectively erased to prevent finding the same one twice. An example of the output from a ParaMEME analysis is given in Figure 4.15. 

### **The Gibbs Sampler** 

Another statistical method for finding motifs in sequences is the Gibbs sampler. The method is similar in principle to the EM method described above, but the algorithm is different. Like the EM method, given a set of sequences, the Gibbs sampler searches for the statistically most probable motifs and can find the optimal width and number of these motifs in each sequence (Lawrence et al. 1993; Liu et al. 1995; Neuwald et al. 1995). The source code of the program code is available by anonymous FTP from ncbi.nlm.nih.gov/pub/neuwald/gibbs9-95. A combinatorial approach of the Gibbs sampler and MOTIF may be used to make blocks at the BLOCKS Web site (http:// www.blocks.fhcrc.org/). The expected number of blocks in the search is one block for approximately each 40 residues of sequence. The Gibbs sampler is also an option of the msa block-alignment and editing program MACAW (Schuler et al. 1991), which runs on MS-DOS, Macintosh, and other computer platforms and is available by anonymous FTP from ncbi.nlm.nih.gov/pub/schuler/macaw. 



<!-- Start of picture text -->
Wednesday, June 28, 000 MEME - Submission Form (Aovanced Version) Pages<br>MEME -- Multiple EM for Motif Elicitation: Version 2.2<br>Motif discovery tool<br>Data Submission Form - Advanced Version es |<br>Your data will be processed on the Cray-T3E supercomputer at the SupercomputerSan Diego Center and the results will be sent to you by e-mail.<br>Please enter the e-mail address where you would like your results sent:<br>[Optional] Please enter a brief description of your sequences.<br>numberPlease enterof formats.the sequencesPlease whichenter either: you believe share one or more motifs. The sequences may contain no more than 100,000 characters total in any ofa large<br>1, the name of a file containing the| sequences here:<br>How do you think theistril occurrences of 7 MEMEmotifs.any motif Wide can choose is recommended the width ofif there  each are motif fewer favoring than  short 10 occurrences  or wide of<br>sea‘Oneee? per owsequenceMistributed among the Howike tomanylook differentfor motifs would you jeportedbutton below, toin have your or that enter sequences. width. SelectChoosing the width a number you want will with cause the all select motifs<br>Zenoor ont nor sequence "short", "wide" or any numbera width in from the text 2 to window. 300. (If you Legal enter choices something are  in<br>Any numberof repetitions the text window, it will override what is shown on the select button.)<br>Brief output format:<br>ADVANCED OPTIONS<br>‘Shuffle letters in input sequences:<br>DNA-ONLY OPTIONS<br>complementary strand, 5' to 3° (inverse complement)<br>Additional strands/directionsto search:<br>complementary strand, 3' to 5"<br>ClickReturn hereto MEME for informationmoreSYSTEM introduction.on MEME.<br>‘You might be interested in trying other motif-making programs such as BLOCK MAKER at the Fred Hutchinson Cancer Research Center .<br>Please send comments and questions to: thailey@sdsc.edu .<br><!-- End of picture text -->

**Figure 4.14.** The MEME Web page. The MEME program finds ungapped motifs (blocks) in unaligned protein or DNA sequences. As indicated, the program can be directed to search for the size and expected number of motifs or can predict motifs based on a statistical analysis based on the EM algorithm described in the text. 

**A. Summary line** MOTIF  1 width = 9 sites = 29.5 **B. Letter-probability matrix** Simplified A **: :   : : : :   :** 1 8 motif letterC **: : : : : : : : :** probability D **:   : : : : : : :** 8 matrix E **: : : : : : : : :** F **: : : : : : : : :** G **: :   : : : : :** 1 9 H **: : : : : : : : :** I 2 **:       : : : :** 212 K **: : : : : : : : :** L 3 **:     : : : : :** 18 M **: : : : : : : : :** N **: : : : :     : :** 89 P **: : : : : : : : :** Q **: : : : : : : : :** R **: : : : : : : : :** S **: : : : : : : : :** T **: : : : : : : : :** V 3 **:   :   : : : :** 3 7 W **: : : : : : : : :** Y **: : : : : : : : :** 

**C. Information content of the profile** Information bits 6.2 content 5.6 ( 22.0 bits ) 5.0 4.4 3.7 3.1 *** * *** 2.5 *** * * *** 1.9 *** * * * * * *** 1.2 *** * * * * * * *** 0.6 *** * * * * * * * *** 0.0 **---------** 

###### **D. The multilevel consensus sequence** 

Multilevel VDVLVNNAG consensus L sequence 

**Figure 4.15.** Results produced by a MEME analysis of sequences for motifs. The output diagrams are discussed in the text. ( _A_ ) Summary line giving the number of the next motif found in order of statistical significance, width, and expected number of occurrences in the given sequences. ( _B_ ) Simplified motif letter-probability matrix showing the frequency of each amino acid in each column of the matrix. The columns are the columns of the motif. For easier reading, the numbers shown are frequencies rounded to the nearest one-tenth and multiplied by 10, and zeros are shown as colons. ( _C_ ) The information content of the profile is given in a diagram. Basically, the diagram shows the degree of amino acid variation in each column of the profile: the lower the value, the greater the variation. The scale is logarithmic to the base 2 (bits). The total of all columns is also shown. The subject of information content is discussed in greater detail below under position-specific scoring matrices. ( _D_ ) The multilevel consensus sequence shows all letters in each column of the motif that occur with a frequency of �0.2. _Continued_ . 

**180** I CHAPTER 4 

###### **E. The next motif** 

|Motif 1 in<br>|BLOC<br>|KS forma<br>|t||
|---|---|---|---|---|
|BL MOTIF 1;  width =|9;   se|qs = 33|||
|2BHD_STREX|(|81)|VDGLVNNAG|1|
|3BHD_COMTE|(|81)|LNVLVNNAG|1|
|ADH_DROME|(|86)|VDVLINGAG|1|
|AP27_MOUSE|(|77)|VDLLVNNAA|1|
|BA72_EUBSP|(|86)|LDVMINNAG|1|
|BDH_HUMAN|(|138)|MWGLVNNAG|1|
|BPHB_PSEPS|(|79)|IDTLIPNAG|1|
|BUCD_KLETE|(|80)|FNVIVNNAG|1|
|DHES_HUMAN|(|84)|VDVLVCNAG|1|
|DHGB_BACME|(|87)|LDVMINNAG|1|
|DHMA_FLAS1|(|198)|VDVTGNNTG|1|
|ENTA_ECOLI|(|73)|LDALVNAAG|1|
|FIXR_BRAJA|(|112)|LHALVNNAG|1|
|GUTD_ECOLI|(|82)|VDLLVYSAG|1|
|HDE_CANTR|(|396)|IDILVNNAG|1|
|HDHA_ECOLI|(|89)|VDILVNNAG|1|
|NODG_RHIME|(|81)|VDILVNNAG|1|
|RIDH_KLEAE|(|89)|LDIFHANAG|1|
|YINL_LISMO|(|83)|VDAIFLNAG|1|
|YRTP_BACSU|(|84)|IDILINNAG|1|
|CSGA_MTXXA|(|13)|VDVLINNAG|1|
|DHB2_HUMAN|(|161)|LWAVINNAG|1|
|DHB3_HUMAN|(|125)|IGILVNNVG|1|
|DHCA_HUMAN|(|83)|LDVLVNNAG|1|
|FVT1_HUMAN|(|115)|VDMLVNCAG|1|
|HMTR_LEIMA|(|103)|CDVLVNNAS|1|
|MAS1_AGRRA|(|320)|IDGLVNNAG|1|
|PCR_PEA|(|165)|LDVLINNAA|1|
|YURA_MYXXA<br>//|(|90)|LDLVVANAG|1|



**Figure 4.15.** _Continued_ . ( _E_ ) Possible examples of the motif in the training set are shown. This list is based on using a position-dependent scoring matrix (log-odds matrix) to search each sequence. The threshold score for displaying a site is chosen such that the expected number of incorrect assignments will equal the expected number of missed but correct assignments. Positions before and after the motif are also shown. C _ontinued_ . 

MULTIPLE SEQUENCE ALIGNMENT I **181** 

###### **F. Possible examples of motif 1 in the training set** 

|Sequence name|Start|Score||Site||
|---|---|---|---|---|---|
|-------------------|------|------||------||
|2BHD_STREX|81|28.80|VAYAREEFGS|VDGLVNNAG|ISTGMFLETE|
|3BHD_COMTE|81|25.99|MAAVQRRLGT|LNVLVNNAG|ILLPGDMETG|
|ADH_DROME|86|22.33|LKTIFAQLKT|VDVLINGAG|ILDDHQIERT|
|AP27_MOUSE|77|24.36|TEKALGGIGP|VDLLVNNAA|LVIMQPFLEV|
|BA72_EUBSP|86|26.39|VGQVAQKYGR|LDVMINNAG|ITSNNVFSRV|
|BDH_HUMAN|138|23.46|PFEPEGPEKG|MWGLVNNAG|ISTFGEVEFT|
|BPHB_PSEPS|79|18.60|ASRCVARFGK|IDTLIPNAG|IWDYSTALVD|
|BUDC_KLETE|80|20.97|VEQARKALGG|FNVIVNNAG|IAPSTPIESI|
|DHES_HUMAN|84|25.67|AARERVTEGR|VDVLVCNAG|LGLLGPLEAL|
|DHGB_BACME|87|26.39|VQSAIKEFGK|LDVMINNAG|MENPVSSHEM|
|DHMA_FLAS1|198|16.36|ILVNMIAPGP|VDVTGNNTG|YSEPRLAEQV|
|ENTA_ECOLI|73|21.90|CQRLLAETER|LDALVNAAG|ILRMGATDQL|
|FIXR_BRAJA|112|23.67|EVKKRLAGAP|LHALVNNAG|VSPKTPTGDR|
|GUTD_ECOLI|82|17.17|SRGVDEIFGR|VDLLVYSAG|IAKAAFISDF|
|HDE_CANTR|92|20.90|VETAVKNFGT|VHVIINNAG|ILRDASMKKM|
|HDE_CANTR|396|29.32|IKNVIDKYGT|IDILVNNAG|ILRDRSFAKN|
|HDHA_ECOLI|89|30.18|ADFAISKLGK|VDILVNNAG|GGGPKPFDMP|
|NODG_RHIME|81|30.18|GQRAEADLEG|VDILVNNAG|ITKDGLFLHM|
|RIDH_KLEAE|89|16.02|LQGILQLTGR|LDIFHANAG|AYIGGPVAEG|
|YINL_LISMO|83|14.65|VELAIERYGK|VDAIFLNAG|IMPNSPLSAL|
|YRTP_BACSU|84|27.41|VAQVKEQLGD|IDILINNAG|ISKFGGFLDL|
|CSGA_MYXXA|13|28.94|AFATNVCTGP|VDVLINNAG|VSGLWCALGD|
|DHB2_HUMAN|161|19.62|KVAAMLQDRG|LWAVINNAG|VLGFPTDGEL|
|DHB3_HUMAN|125|18.63|HIKEKLAGLE|IGILVNNVG|MLPNLLPSHF|
|DHCA_HUMAN|83|30.23|RDFLRKEYGG|LDVLVNNAG|IAFKVADPTP|
|FVT1_HUMAN|115|24.21|IKQAQEKLGP|VDMLVNCAG|MAVSGKFEDL|
|HMTR_LEIMA|103|24.02|VAACYTHWGR|CDVLVNNAS|SFYPTPLLRN|
|MAS1_AGRRA|320|27.93|VTAAVEKFGR|IDGLVNNAG|YGEPVNLDKH|
|PCR_PEA|165|23.97|VDNFRRSEMP|LDVLINNAA|VYFPTAKEPS|
|YURA_MYXXA|90|18.59|IRALDAEAGG|LDLVVANAG|VGGTTNAKRL|



**Figure 4.15.** _Continued_ . ( _F_ ) The next motif is given in the format used for the BLOCKS database (http://www. blocks.fhcrc.org/blocks). The predicted locations of this motif in each sequence and the probability that the motif starts at that location are shown. The sites reported depend on the motif search model used: (1) OOPS, the most probable location in each sequence is given; (2) ZOOPS, the most probable location in each sequence is reported but only probabilities greater than 0.5 (a significant level for Bayesian statistics); TCM, all positions in each sequence with probabilities � 0.5 are shown. _Continued_ . 

**182** I CHAPTER 4 

###### **G. Position-specific scoring matrix** 

Log-odds matrix: alength = 20 w = 9 n = 9732 bayes = 8.36118 

|-2.725|0.818|-5.204|-4.539|-0.082|-4.432|-3.515|1.560|-4.218|1.814|0.701|-4.126|-3.146|-3.848 .|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|-3.441|-3.841|-4.023|-1.204|-4.313|-2.395|-0.889|-4.226|-4.009|-4.571|-3.882|-0.220|-4.682|.<br>-3.547|
|-0.768|-2.342|-4.756<br>|-4.189|-2.319|0.376|-3.154|1.757|-3.870|0.288|0.918|-3.149|-4.229|.<br>-3.492|
|-3.379|-2.600|-5.066|-4.331|-0.586|-5.089|-3.668|-0.081|-4.098|3.045|1.107|-4.393|-4.287|.<br>-3.383|
|-1.373|-1.895|-3.823|-3.574|-1.086|-1.952|-0.466|1.480|-3.565|-2.234|-1.834|-3.701|-3.612|.<br>-3.536|
|-1.879|-0.980|-2.231|-4.187|-3.807|-3.562|-0.892|-3.306|-3.238|-2.753|-3.337|4.193|-2.276|.<br>-2.750|
|-2.460|-0.912|-2.252<br>-|4.176|-3.833|-2.391|-0.968|-3.339|-3.262|-4.256|-3.364|4.217|-4.026|.<br>-2.768|
|-3.475|-1.137|-3.874|-3.535|-3.304|-2.080|-2.080|-2.826|-3.544|-3.127|-2.263|-3.592|-4.599|.<br>-3.533|
|-0.693|-3.833|-3.137|-3.879|-4.963|3.663|-3.647|-3.364|-3.716|-5.287|-4.212|-2.849|-4.518|.<br>-4.155|



###### **H. Motif letter-frequency matrix** 

|Letter-p|robabilit|y matrix|: alength|= 20 w|= 9 n = 9|732||||||
|---|---|---|---|---|---|---|---|---|---|---|---|
|0.011063|0.032022|0.001403|0.002682|0.038055|0.003212|0.001962|0.165990|0.003143|0.322510|0.037503|0.011063|
|0.006738|0.001268|0.841023|0.027061|0.002026|0.013178|0.012108|0.003008|0.003632|0.003860|0.001564|0.011063|
|0.124630|0.003583|0.001915|0.003418|0.008070|0.089951|0.002520|0.190255|0.004000|0.112000|0.043590|0.011063|
|0.007032|0.002996|0.001544|0.003098|0.026845|0.002037|0.001765|0.053213|0.003415|0.756853|0.049683|0.011063|
|0.028238|0.004883|0.003655|0.005236|0.018977|0.017917|0.016240|0.156947|0.004942|0.019499|0.006470|0.011063|
|0.019895|0.009211|0.011023|0.003422|0.002878|0.005871|0.012089|0.005691|0.006199|0.013606|0.002282|0.011063|
|0.013301|0.009656|0.010865|0.003449|0.002827|0.013217|0.011467|0.005564|0.006098|0.004800|0.002240|0.011063|
|0.813801|0.008259|0.003529|0.005378|0.004079|0.016396|0.005304|0.007937|0.005014|0.010499|0.004806|0.011063|
|0.045249|0.001275|0.005879|0.004237|0.001291|0.878064|0.001790|0.005467|0.004450|0.002354|0.001244|0.011063|



**Figure 4.15.** _Continued_ . ( _G_ ) Position-specific scoring matrix. This matrix is a log-odds matrix calculated by taking the log (base 2) of the ratio of the observed to expected counts for each amino acid in each column of the profile. Columns and rows in the matrix correspond to the amino acids in each column and positions of the motif, respectively. The counts for each column may have additional pseudocounts added to compensate for zero occurrences of an amino acid in a column or for a small number of sequences, as discussed below for this type of matrix. ( _H_ ) Motif letter-frequency matrix is given, showing the frequency of amino acid found in each column of the profile. Columns and rows correspond to the amino acids in each column and rows to columns in the motif, respectively. Shown also are the numbers of types of residues, the width of the motif, and number of characters in the sequences. Only portions of the output are shown. 

To understand the algorithm, consider a simple example using the Gibbs sampler algorithm to locate a single 20-residue-long motif in 10 sequences, each 200 residues long, as was done above to illustrate the EM algorithm. The method iterates through two steps. In the first step, the predictive update step, a random start position for the motif is chosen for all sequences but for one that is chosen at random or in a specified order. So let us choose sequence 1 as the outlier and use the other 9 to find an initial guess of the motif. These other 9 sequences are aligned with random overlaps. The following figure illustrates how this initial motif is located (an x equals 20 sequence positions, M indicates the random location of the motif chosen for each sequence, and � the 20 initially aligned motif positions). 

The objective is to find the most probable pattern common to all of the sequences by sliding them back and forth until the ratio of the motif probability to the background probability is a maximum. This is accomplished by first using the initial alignment shown above to estimate the residue frequencies in each column of the motif, and the sequence residues 

MULTIPLE SEQUENCE ALIGNMENT I **183** 



<!-- Start of picture text -->
Steps of the Gibbs sampler algorithm.<br>A. Estimate the amino acid frequencies in the motif columns of all but<br>1 sequence. Also obtain background<br><!-- End of picture text -->

#### **Motif** 

xxxMxxxxxx xxxMxxxxxx xxxxxxMxxx xxxxxxMxxx xxxxxMxxxx xxxxxMxxxx xMxxxxxxxx xMxxxxxxxx xxxxxxxxxM xxxxxxxxxM Mxxxxxxxxx Mxxxxxxxxx xxxxMxxxxx xxxxMxxxxx xMxxxxxxxx xMxxxxxxxx xxxxxxxxMx xxxxxxxxMx 

Random start Location of motif in each positions chosen sequence provides first estimate of motif composition 

```
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
M ->M ->M ->M ->M ->
```

```
xxxxxxxMxx
```

that are not included in the motif to estimate the background residue frequencies. For example, if these sequences are DNA sequences and the first column of the estimated motif in the 10 sequences includes 3 Gs, then the value for _fg_ , column1 � 3/9 � 0.33. Similarly, let _ft_ , column2 � 1/9 � 0.11 for illustration. These frequencies are determined for each of the 20 columns in our example. Similarly, if there are 240 Gs among the 10 � x 80 � = 800 sequence positions not within the estimated motif, then _fg_ , background � 240/800 � 0.30. Also let _ft_ , background � 180/800 � 0.225. If the first two positions in sequence 1 are G and T in that order, then the probability of the motif starting at position 1, Q1, is calculated as 0.33 � x 0.11 � x . . . . . . _x f_ last base, column20. The background probability of this first possible motif, _P_ 1, is also calculated as 0.30 � x 0.225 � x . . . .. _x f_ last base, background. 

**184** I CHAPTER 4 

_Note the difference between the Gibbs sampler method and the EM method, which calculates the probability of the entire sequence using the motif column frequencies within the motif and the background frequencies elsewhere._ 

The ratio _Q_ 1/ _P_ 1 is designated as weight _A_ 1 for motif position 1 in sequence 1. A1s are then calculated for all other 100 � 20 � 1 � 81 possible locations of the 20-residue-long motif in sequence 1. These weights are then normalized by dividing each weight by their sum to give a probability for each motif position. From this probability distribution, a random start position is chosen for position 1. In so doing, the chance of choosing a particular position is proportional to the weight of that position so that a higher scoring position is more likely to be chosen. (You can think of a bag with 81 kinds of balls, with the number of each ball proportional to the weight or probability of that kind. Drawing a random ball will favor the more prevalent ones.) This position in the left-out sequence is then used as an estimate of the location for the motif in sequence 1. The procedure is then repeated. Select the next sequence to be scanned, align the motifs in the other 9 sequences with sequence 1 now using the estimated location found above, and so on. This process is repeated until the residue frequencies in each column of the motif do not change. For different starting alignments, the number of iterations needed may range from several hundred to several thousand. 

As the above cycles are repeated, the more accurate the initial estimate of the motif in the aligned sequences, the more accurate the pattern location in the outlier sequence. The second step in the algorithm tends to move the sequence alignments in a direction that favors a better score but also has a random element to search for other possible better locations. When correct start positions have been selected in several sequences by chance, the compositions of the motif columns begin to reflect a pattern that the algorithm can search for in the other sequences, and the method converges on the optimal motif and the probability distribution of the motif location in each sequence. 

Several additional procedures are used to improve the performance of the algorithm. 

1. For a correct Bayesian statistical analysis, the amino acid counts in the motif and the background in the outlier sequence are estimated and added to the counts in the remaining aligned sequences. This step is the equivalent of combining prior and updated information to improve the estimation of the motif. These counts may be estimated by Dirichlet mixtures (see discussion of HMMs, p. 189), which give frequencies expected based on prior information from amino acid distributions (Liu et al. 1995). The missing background counts for each residue _bi_ are estimated by the formula _bi_ � _fi x_ , _B_ where _B_ is chosen based on experience with the method as √ _N_ , the number of sequences in the motif, and _fi_ is the frequency of residue _i_ in the sequences (Lawrence et al. 1993). 

2. Another feature is a procedure to prevent the algorithm from getting locked in a suboptimal solution. In the HMM method (see below), noise is introduced for this purpose. In the Gibbs sampler, after a certain number of iterations, the current alignments are shifted a certain number of positions to the right and left, and the scores from these shifted positions are found. A probability distribution of these scores is then used as a basis for choosing a new random alignment. 

3. The results of a range of motif widths can be investigated. The major difficulty in exploring motif width is to arrive at a criterion for comparing the resulting scores. One suitable measure is to optimize the average information (see below) per free parameter in the motif, a value that can be calculated (Lawrence et al. 1993; Liu et al. 1995). The number of free parameters for proteins is 20 � 1 � 19, and for DNA, 4 � 1 � 3, times the model width. 

4. The method can be readily extended to search for multiple motifs in the same set of sequences. 

5. The method has been extended to seek a pattern in only a fraction of the input sequences. 

MULTIPLE SEQUENCE ALIGNMENT I **185** 

The Gibbs sampler was used to align 30 helix-turn-helix DNA-binding domains showing very little sequence similarity. The information per parameter criterion was used to find the best motif width. Multiple motifs were found in lipocalins, a family with quite dissimilar motif sequences separated by variable spacer regions, and also in protein isoprenyltransferase subunits, which have very large numbers of repeats of several kinds (Lawrence et al. 1993). Thus, the method is widely applicable for discovering complex and variable motifs in proteins. 

### **Hidden Markov Models** 

The HMM is a statistical model that considers all possible combinations of matches, mismatches, and gaps to generate an alignment of a set of sequences (Fig. 4.16). A model of a sequence family is first produced and initialized with prior information about the sequences. A set of 20–100 sequences or more is then used as data to train the model. The trained model may then be used to produce the most probable msa as posterior information. Alternatively, the model may be used to search sequence databases to identify additional members of a sequence family. A different HMM is produced for each set of sequences. HMMs have been previously used very successfully for speech recognition, and an excellent review of the methodology is available (Rabiner 1989). In addition to their use in producing multiple sequence alignments (Baldi et al. 1994; Krogh et al. 1994; Eddy 1995, 1996), HMMs have also been used in sequence analysis to produce an HMM that represents a sequence profile (a profile HMM), to analyze sequence composition and patterns (Churchill 1989), to locate genes by predicting open reading frames (Chapter 8), and to produce protein structure predictions (Chapter 9). Pfam, a database of profiles that represent protein families, is based on profile HMMs (Sonhammer et al. 1997). 

HMMs often provide a msa as good as, if not better than, other methods. The approach also has a number of other strong features: It is well grounded in probability theory, no sequence ordering is required, insertion/deletion penalties are not needed, and experimentally derived information can be used. Two disadvantages to using HMMs are that at least 20 sequences and sometimes many more are required to accommodate the evolutionary history (see Mitchison and Durbin 1995). The HMM can be used to improve an existing heuristic alignment. The two HMM programs in common use are Sequence Alignment and Modeling Software System, or SAM (Krogh et al. 1994; Hughey and Krogh 1996), and HMMER (see Eddy 1998). The software is available at http://www.cse.ucsd.edu/ research/compbio/sam.html and http://hmmer.wustl.edu/. The algorithms used for producing HMMs are extensively discussed in Durbin et al. (1998). A comparison of HMMs with other methods is given at the end of this section. 

The HMM representation of a section of multiple sequence alignment that includes deletions and insertions was devised by Krogh et al. (1994) and is shown in Figure 4.6. This HMM generates sequences with various combinations of matches, mismatches, insertions, and deletions, and gives these a probability, depending on the values of the various parameters in the model. The object is to adjust the parameters so that the model represents the observed variation in a group of related protein sequences. A model trained in this manner will provide a statistically probable msa of the sequences. 

As illustrated in Figure 4.6, the object is to calculate the best HMM for a group of sequences by optimizing the transition probabilities between states and the amino acid compositions of each match state in the model. The sequences do not have to be aligned to use the method. Once a reasonable model length reflecting the expected length of the sequence alignment is chosen, the model is adjusted incrementally to predict the sequences. Several methods for training the model in this fashion have been described (Baldi et al. 1994; Krogh et al. 1994; Eddy et al. 1995; Eddy 1996; Hughey and Krogh 1996; 

**186** I CHAPTER 4 

#### **A. Sequence alignment** 

N • F L S N • F L S N K Y L T Q • W – T 

RED POSITION REPRESENTS ALIGNMENT IN COLUMN GREEN POSITION REPRESENTS INSERT IN COLUMN PURPLE POSITION REPRESENTS DELETE IN COLUMN 

#### **B. Hidden Markov model for sequence alignment** 



<!-- Start of picture text -->
D3 D2 D3 D4<br>I0 I1 I2 I3 I4<br>BEG M1 M2 M3 M4 END<br>match state insert state delete state transition probability<br><!-- End of picture text -->

**Figure 4.16.** Relationship between the sequence alignment and the hidden Markov model of the alignment (Krogh et al. 1994). This particular form for the HMM was chosen to represent the sequence, structural, and functional variation expected in proteins. The model accommodates the identities, mismatches, insertions, and deletions expected in a group of related proteins. ( _A_ ) A section of a multiple sequences alignment. The illustration shows the columns generated in a multiple sequence alignment. Each column may include matches and mismatches ( _red_ positions), insertions ( _green_ positions), and deletions ( _purple_ position). ( _B_ ) The HMM. Each column in the model represents the possibility of a match, insert, or delete in each column of the alignment in _A_ . The HMM is a probabilistic representation of a section of a msa. Sequences can be generated from the HMM by starting at the beginning state labeled BEG and then by following any one of many pathways from one type of sequence variation to another (states) along the state transition arrows and terminating in the ending state labeled END. Any sequence can be generated by the model and each pathway has a probability associated with it. Each square match state stores an amino acid distribution such that the probability of finding an amino acid depends on the frequency of that amino acid within that match state. Each diamond-shaped insert state produces random amino acid letters for insertions between aligned columns and each circular delete state produces a deletion in the alignment with probability 1. For example, one of many ways of generating the sequence **N K Y L T** in the above profile is by the sequence BEG → M1 → I1 → M2 → M3 → M4 → END. Each transition has an associated probability, and the sum of the probabilities of transitions leaving each state is 1. The average value of a transition would thus be 0.33, since there are three transitions from most states (there are only two from M4 and D4, hence the average from them is 0.5). For example, if a match state contains a uniform distribution across the 20 amino acids, the probability of any amino acid is 0.05. Using these average values of 0.33 or 0.5 for the transition values and 0.05 for the probability of each amino acid in each state, the probability of the above sequence **N K Y L T** is the product of all of the transition probabilities in the path BEG → M1 → I1 → M2 → M3 → M4 → END, and the probability that each state will produce the corresponding amino acid in the sequences, or 0.33 � 0.05 � 0.33 � 0.05 � 0.33 � 0.05 � 0.33 � 0.05 � 0.33 � 0.05 � 0.5 � 6.1 � 10<sup>�10</sup> . Since these probabilities are very small numbers, amino acid distributions and transition probabilities are converted to log odds scores, as done in other statistical methods (see pp. 176–177), and the logarithms are added to give the overall probability score. The secret of the HMM is to adjust the transition values and the distributions in each state by training the model with the sequences. The training involves finding every possible pathway through the model that can produce the sequences, counting the number of times each transition is used 

_Continued._ 

MULTIPLE SEQUENCE ALIGNMENT I **187** 

Durbin et al. 1998). For example, an EM algorithm from speech recognition methods known as the Baum-Welch algorithm is used as follows: 

1. The model is initialized with estimates of transition probabilities and amino acid composition for each match and insert date. If an initial alignment of the sequences is known, or some other kinds of data suggest which sequence positions are the same, then these data may be used in the model. For other cases, the initial distribution of amino acids to be used in each state is described below. The initial transition probabilities generally favor transitions from one match state to the next rather than favoring insert and delete states, which build more uncertainty into a sequence motif. 

2. All possible paths through the model for generating each sequence in turn are examined. There are many possible such paths for each sequence. This procedure would normally require a huge amount of time computationally. Fortunately, an algorithm, the forward-backward algorithm, reduces the number of computations to the number of steps in the model times the total length of the training sequences. This calculation provides a probability of the sequence, given all possible paths through the model, and, from this value, the probability of any particular path may be found. Another algorithm, the Baum-Welch algorithm, then counts the number of times a particular stateto-state transition is used and a particular amino acid is required by a particular match state to generate the corresponding sequence position. 

3. A new version of the HMM is produced that uses the results found in step 2 to generate new transition probabilities and match-insert state compositions. 

4. Steps 3 and 4 are repeated up to 10 more times until the parameters do not change significantly. 

5. The trained model is used to provide the most likely path for each sequence, as described in Figure 4.16. The algorithm used for this purpose, the Viterbi algorithm, does not have to go through all of the possible alignments of a given sequence to the HMM to find the most probable alignment, but instead can find the alignment by a dynamic programming technique very much like that used for the alignment of two sequences, discussed in Chapter 3. The collection of paths for the sequences provides a msa of the sequences with the corresponding match, insert, and delete states for each sequence. The columns in the msa are defined by the match states in the HMM such that amino acids from a particular match state are placed in the same column. For columns that do not correspond to a match state, a gap is added. 

6. The HMM may be used to search a sequence database for additional sequences that share the same sequence variation. In this case, the sum of the probabilities of all possible sequence alignments to the model is obtained. This probability is calculated by the forward component of the forward-backward algorithm described above. This analysis 

and which amino acids were required by each match and insert state to produce the sequences. This training procedure leaves a memory of the sequences in the model. As a consequence, the model will be able to give a better prediction of the sequences. Once the model has been adequately trained, of all the possible paths through the model that can generate the sequence **N K Y L T** , the most probable should be the match-insert-3 match combination (as opposed to any other combination of matches, inserts, and deletions). Likewise, the other sequences in the alignment would also be predicted with highest probability as they appear in the alignment; i.e., the last sequence would be predicted with highest probability by the path matchmatch-delete-match. In this fashion, the trained HMM provides a multiple sequence alignment, such as shown in _A_ . For each sequence, the objective is to infer the sequence of states in the model that generate the sequences. The generated sequence is a Markov chain because the next state is dependent on the current one. Because the actual sequence information is hidden within the model, the model is described as a hidden Markov model. 

**188** I CHAPTER 4 

gives a type of distance score of the sequence from the model, thus providing an indication of how well a new sequence fits the model and whether the sequence may be related to the sequences used to train the model. In later derivations of HMMs, the score was divided by the length of the sequence because it was found to be lengthdependent. A _z_ score giving the number of standard deviations of the sequence lengthcorrected score from the mean length-corrected score is therefore used (Durbin et al. 1998). 

Recall that for the Bayes block aligner, the initial or prior conditions were amino acid substitution matrices, block numbers, and alignments of the sequences. The sequences were then used as new data to examine the model by producing scores for every possible combination of prior conditions. By using Bayes’ rule, these data provided posterior probability distributions for all combinations of prior information. Similarly, the prior conditions of the HMM are the initial values given to the transition values and amino acid compositions. The sequences then provide new data for improving the model. Finally, the model provides a posterior probability distribution for the sequences and the maximum posterior probability for each sequence represented by a particular path through the model. This path provides the alignment of the sequence in the msa; i.e., the sequence plus matches, inserts, and deletes, as described in Figure 4.16. 

The success of the HMM method depends on having appropriate initial or prior conditions, i.e., a good prior model for the sequences and a sufficient number of sequences to train the model. The prior model should attempt to capture, for example, the expected amino acid frequencies found in various types of structural and functional domains in proteins. As the distributions are modified by adding amino acid counts from the training sequences, new distributions should begin to reflect common patterns as one moves through the model and along the sequences. It is important that the model reflect not only the patterns in the training sequences, but also pattern variations that might be present in other members of the same protein family. Otherwise, the model will only recognize the training sequences but not other family members. Thus, some smoothing of the amino acid frequencies is desirable, but not to the extent of suppressing highly conserved pattern information from the training sequences. Such problems are avoided by using a method called regularization to avoid overfitting the data to the model. Basically, the method involves using a carefully designed amino acid distribution as the prior condition and then modifying this distribution in a manner that uses training information in a complementary manner. 

Rather than using simple amino acid composition as a prior condition for the match states in the HMM, amino acid patterns that capture some of the important features of protein structure and function have been used with considerable success (Sjölander et al. 1996). Other prior conditions include using Dayhoff PAM or BLOSUM amino acid substitution matrices modified by adding additional counts (pseudocounts) to smooth the distributions (Tatusov et al. 1994; Eddy 1996; Henikoff and Henikoff 1996; Sonnhammer et al. 1997; and see Chapter 2). Sjölander et al. (1996) have prepared particularly useful amino acid distributions called Dirichlet mixtures to use as prior information in the match states of the HMM. These mixtures provide amino acid compositions that have proven to be useful for the detection of weak but significant sequence similarity. As an example, the amino acid frequencies that are characteristic of a particular set of nine blocks in the BLOCKS database have been determined. These blocks represent amino acid frequencies that are favored in certain chemical environments such as aromatic, neutral, and polar residues and are useful for detecting such environments in test sequences. The nine-component system has been used successfully for producing an HMM for globin sequences (Hughey and Krogh 1996). To use these frequencies as prior information, they are treated 

MULTIPLE SEQUENCE ALIGNMENT I **189** 

as possible posterior distributions that could have generated the given amino acid frequencies as posterior probabilities. The probability of a particular amino acid distribution given a known frequency distribution, i.e., 100 A, 67 G, 5 C, etc., where _p_ A is the probability of A given by the frequency of A, _p_ G the probability of G, etc., and _n_ is the total number of amino acids given by the multinomial distribution 



The prior distribution for the multinomial distribution is the Dirichlet distribution (Carlin and Louis 1996), whose formulation is similar to that given in Equation 6 with a similar set of parameters but with factorial and powers reduced by 1. The idea behind using this particular distribution is that if additional sequence data with a related pattern are added, then by the Bayesian procedure of multiplying prior probabilities with the likelihood of the new data to obtain the posterior distribution, the probability of finding the correct frequency of amino acids is favored statistically. Because the amino acid frequencies in the test sequences could be any one of several alternatives, a prior distribution that reflects these several choices is necessary. There is a method for weighting the prior distributions expected for several different multinomial distributions into a combined frequency distribution, the Dirichlet mixture. Calculation of these mixtures is a complex mathematical procedure (Sjölander et al. 1996). Dirichlet mixtures recommended for use in aligning proteins by the HMM method have been described previously (Karplus 1995) and are available from http://www.cse.ucsc.edu/research/compbio/. After the prior amino acid frequencies are in place in the match states of the model, these are modified by training the HMM with the sequences, as described in steps 2 and 3 above. For each match state in the model, a new frequency for each amino acid is calculated by dividing the sum of all new and prior counts for that amino acid by the new total of all amino acids. In this fashion, the new HMM (step 4 above) reflects a combination of expected distributions averaged over patterns in the Dirichlet mixture and patterns exhibited in the training sequences. A similar method is used to refashion the transition probabilities in the HMM during training following manual insertion of initial values. 

Another consideration in using HMMs is the number of sequences. If a good prior model such as the above Dirichlet distribution is used, it should be possible to train the HMM with as few as 20 sequences (SAM manual; Eddy 1996; Hughey and Krogh 1996). In general, the smaller the sequence number, the more important the prior conditions. If the number of sequences is �50, the initial conditions play a lesser role because the training step is more effective. As with any msa method, the more sequence diversity, the more challenging the task of aligning sequences with HMMs. HMMs are also more effective if methods to inject statistical noise into the model are used during the training procedure. As the model is refashioned to fit the sequence data, it sometimes goes into a form that provides locally optimal instead of globally optimal alignments of the sequences. One of several noise injection methods (Baldi et al. 1994; Krogh et al. 1994; Eddy et al. 1995; Eddy 1996; Hughey and Krogh 1996) may be used in the training procedure. One method called simulated annealing is used by SAM (Hughey and Krogh 1996). A user-defined number of sequences are generated from the model at each cycle and the counts so generated are added to those from the training sequences. The noise generated in this way is reduced as the cycle number is increased. Finally, the HMM program SAM has a built-in feature of model surgery during training. If a match state is used by fewer than half of the sequences, it is deleted. These same sequences then have to use an insert state in the revised model. Similarly, if an insert state is used by more than half of the sequences, a number of addi- 

**190** I CHAPTER 4 

tional match states equal to the average number of insertions is added, and the model has to be revised accordingly. These fractions may be varied in SAM to test the effect on the type of HMM model produced (Hughey and Krogh 1996). 

In trying to produce an HMM for a set of related sequences, the recommended procedure is to produce several models by varying the prior conditions. Using regularization by adding prior Dirichlet mixtures to the match states produces models that are more representative of the protein family from which the training sequences are derived. Varying the noise and model surgery levels is another way to vary the training procedure and the HMM model. The best HMM model is the one that predicts a family of related sequences with the lowest and most narrow distribution of NLL scores. An example of a portion of an HMM trained on a set of globin sequences is shown in Figure 4.17. 

### **Motif-based Hidden Markov Models** 

The program Meta-MEME uses the HMM method to find motifs (conserved sequence domains) in a set of related protein sequences and the spacer regions between them (Grundy et al. 1997) and is built in part on the HMM program HMMER (Eddy et al. 1995). A similar method was originally used to analyze prokaryotic promoters with two conserved patterns separated by a variable spacer region (Cardon and Stormo 1992). A Meta-MEME analysis may be performed at http://www.sdsc.edu/MEME using the University of California at San Diego Supercomputing Center. The use of hidden Markov models for producing a global msa is described in the above section. A problem with HMMs is that the training set has to be quite large (50 or more sequences) to produce a useful model for the sequences. For a smaller number of sequences, it is possible to obtain a model if suitable prior data are used, and an amino acid frequency that is a mixture of frequencies characteristic of certain structural domains (the Dirichlet mixture) is used as prior information of the match states of the model. This mixture is a reasonable guess of combinations of amino acid patterns that are likely to be found. A difficulty in training the HMM residues is that many different parameters must be found (the amino acid distributions, the number and positions of insert and delete states, and the state transition frequencies add up to thousands of parameters) to obtain a suitable model, and the purpose of the prior and training data is to find a suitable estimate for all of these parameters. When trying to make an alignment of short sequence fragments to produce a profile HMM, this problem is worsened because the amount of data for training the model is even further reduced. 

Two methods are used by Meta-MEME to circumvent this problem. First, another pattern-finding algorithm, the EM algorithm (discussed on p. 173), is used to locate ungapped regions that match in the majority of the sequences. Second, a simplified HMM with a much reduced number of parameters is produced. The model includes a series of match states that model the patterns located by MEME with transition probabilities of 1 between them and a single insert state between each of these patterns, as illustrated in Figure 4.18. As a result, fewer parameters need to be used, mostly for the amino acid frequencies in the match states. 

The most probable order and spacing of the patterns is next found. Another program (Motif Alignment and Search Tool, or MAST; Bailey and Gribskov 1997) is used for this purpose. MAST searches a sequence database for the patterns and reports the database sequences that have the statistically most significant matches. The order and spacing of the patterns found in the highest-scoring database sequences are then used by Meta-MEME as a basis for designing the number of match and insert states and the transition probabilities for the insert states. The match states are filled with modified Dirichlet mixtures (Baylor 

MULTIPLE SEQUENCE ALIGNMENT I **191** 



<!-- Start of picture text -->
PRRO<br>GOGO SHS OOO<br>® © OOO OOO ©-O @ OO<br>OO LOR OOM OORE<br><!-- End of picture text -->

**Figure 4.17.** HMM trained for recognition of globin sequences. Circles in the top row are delete states that include the position in the alignment; the diamonds in the second row are insert states showing the average length of the insertion, and the rectangles in the bottom row show the amino acid distribution in the match states: V is common at match position 1, L at 2, and so on. The width of each transition line joining these various states indicates the extent of use of that path in the training procedure, and dotted lines indicate a rarely used path. The most used paths are between the match states, but about one-half of the sequences use the delete states at model positions 56–60. Thus, for most of the sequences, the msa or profile will show the first two columns aligned with a V followed by an L, but at 56–60, about one-half of the sequences will have a 5-aminoacid deletion. (Reprinted, with permission, from Krogh et al. 1994 [copyright Academic Press].) 

and Gribskov 1996), and the model is trained by the motif models found by MEME. For the 4Fe-4S ferredoxins, a measure of the success of the HMM for database search, the ROC50 score (see p. 165), was approximately 0.6–0.8 for 4 to 8 training sequences, compared to 0.95–0.96 using an evolutionary profile of 6 to 12 sequences. However, this family was one of the most difficult ones to model, and other families produced an ROC50 of 0.9 or better when trained by 20 or more sequences. 

**192** I CHAPTER 4 



<!-- Start of picture text -->
ee a<br>Motif 1<br>Motif2<br><!-- End of picture text -->

**Figure 4.18.** The HMM used by Meta-MEME to estimate motifs in sequences. (Reprinted, with permission of Oxford University Press, from Grundy et al. 1997.) 

### **<u>POSITION-SPECIFIC SCORING MATRICES</u>** 

Analysis of msas for conserved blocks of sequence leads to production of the position-specific scoring matrix, or PSSM. An example of a PSSM produced by the MEME Web site is shown in Figure 4.15G. The PSSM may be used to search a sequence to obtain the most probable location or locations of the motif represented by the PSSM. Alternatively, the PSSM may be used to search an entire database to identify additional sequences that also have the same motif. Consequently, it is important to make the PSSM as representative of the expected sites as possible. The quality and quantity of information provided by the PSSM also varies for each column in the motif, and this variation profoundly influences the matches found with sequences. This situation can be accurately described by information theory, and the results can be displayed by a colored graph called a sequence logo (see Fig. 4.19). 

The PSSM is constructed by a simple logarithmic transformation of a matrix giving the frequency of each amino acid in the motif. Two considerations arise in trying to tune the PSSM so that it adequately represents the training sequences. First, if the number of sequences with the found motif is large and reasonably diverse, the sequences represent a good statistical sampling of all sequences that are ever likely to be found with that same motif. If a given column in 20 sequences has only isoleucine, it is not very likely that a different amino acid will be found in other sequences with that motif because the residue is probably important for function. In contrast, another column in the motif from the 20 sequences may have several amino acids, and some amino acids may not be represented at all. Even more variation may be expected at that position in other sequences, although the more abundant amino acids already found in that column would probably be favored. Thus, if a good sampling of sequences is available, the number of sequences is sufficiently large, and the motif structure is not too complex, it should, in principle, be possible to obtain frequencies highly representative of the same motif in other sequences also (Henikoff and Henikoff 1996; Sjölander et al. 1996). 

However, the number of sequences for producing the motif may be small, highly diverse, or complex, giving rise to a second level of consideration. If the data set is small, then unless the motif has almost identical amino acids in each column, the column frequencies in the motif may not be highly representative of all other occurrences of the motif. In such cases, 

MULTIPLE SEQUENCE ALIGNMENT I **193** 

it is desirable to improve the estimates of the amino acid frequencies by adding extra amino acid counts, called pseudocounts, to obtain a more reasonable distribution of amino acid frequencies in the column. Knowing how many counts to add is a difficult but fortunately solvable problem. On the one hand, if too many pseudocounts are added in comparison to real sequence counts, the pseudocounts will become the dominant influence in the amino acid frequencies, and searches using the motif will not work. On the other hand, if there are relatively few real counts, many amino acid variations may not be present because of the small sample of sequences. The resulting matrix would then only be useful for finding the sequences used to produce the motif. In such a case, the pseudocounts will broaden the evolutionary reach of the profile to variations in other sequences. Even in this case, the pseudocounts should not drown out but serve to augment the influence of the real counts. In summary, relatively few pseudocounts should be added when there is a good sampling of sequences, and more should be added when the data are more sparse. 

The goal of adding pseudocounts is to obtain an improved estimate of the probability _pca_ that amino acid _a_ is in column _c_ in all occurrences of the blocks, and not just the ones in the present sample. The current estimate of _pca_ is _fca,_ the frequency of counts in the data. A simplified Bayesian prediction improves the estimate of _pca_ by adding prior information in the form of pseudocounts (Henikoff and Henikoff 1996): 



where _nca_ and _bca_ are the real counts and pseudocounts, respectively, of amino acid _a_ in column _c_ , _Nc_ and _Bc_ are the total number of real counts and pseudocounts, respectively, in the column, and _fca_ � _nca_ / _Nc_ . It is obvious that as _bca_ becomes larger, the pseudocounts will have a greater infuence on _pca_ . Furthermore, not only the types of pseudocounts but also the total number added to the column ( _Bc_ ) will influence _pca_ . Finally, fractions such as _pca_ are used to produce the log odds form of the motif matrix, the PSSM, which is the most suitable representation of the data for sequence comparisons. A count and probability of zero for an amino acid _a_ in a given column, which is quite common in blocks, may not be converted to logarithms. Addition of a small number of _bca_ will correct this problem without producing a major change in the PSSM values. An equation similar to Equation 7 is used in the Gibbs sampler (p. 177), except that the number of sequences is _N_ � 1. 

Pseudocounts are added based on simple formulas or on the previous variations seen in aligned sequences. The amino acid substitution matrices, including the Dayhoff PAM and BLOSUM matrices, provide one source of information on amino acid variation. Another source is the Dirichlet mixtures derived as a posterior probability distribution from the amino acid substitutions observed in the BLOCKS database (see HMMs; Sjölander et al. 1996). 

One simple formula that has worked well in some studies is to make _B_ in Equation 7 equal to √ _N_ , where _N_ is the number of sequences, and to allot these counts to the amino acids in proportion to their frequencies in the sequences (Lawrence et al. 1993; Tatusov et al. 1997). As _N_ increases, the influence of pseudocounts will decrease because √ _N_ will increase more slowly. The main difficulties with this method are that it does not take into account known substitutions of amino acids in alignments and the observed amino acid variations from one column in the motif to the next, and it does not add enough pseudocounts when the number of sequences is small. 

The information in scoring matrices may be used to produce an average sequence profile, as illustrated in Figure 4.12. Rather than count amino acids, the scoring table values are averaged between each possible 20 amino acids and those amino acids found in the col- 

**194** I CHAPTER 4 

umn of the scoring matrix. Zero counts in a column are not a problem because amino acids not present are not used in the calculations. Because these averaging methods do not take into account the number of sequences in the block, they do not have the desirable effect of a reduced influence when there is a large number of sequences. 

Another method of using the information from amino acid substitution matrices is to base pseudocounts on these matrices. Recall the log odds form of the matrices is derived by taking the logarithm of the frequency of substitution _qia_ of amino acid _i_ for amino acid _a_ divided by the frequency of occurrence of amino acid _a_ , _pa_ . Then, _bca_ may be estimated from the total number of pseudocounts in the column by (Henikoff and Henikoff 1996), 



_bca_ in column c can also be made to depend on the observed data in that column (Tatusov et al. 1997), which is given by multiplying _Bc_ by the following conditional probabilities. 



where _nci_ is the real count of amino acid _i_ in column _c_ . 

The total number of pseudocounts in each column needs also to be estimated. As described above, one estimate is to make _Bc_ for each column equal to √ _N_ , where _N_ is the number of sequences, but this method does not take into account the differences between columns and, for a small number of sequences, the total number of pseudocounts is not sufficient. Allowing _Bc_ to be a constant that can exceed _Nc_ overcomes this limitation but still does not take into account variations in amino acid frequencies between columns, such that a column with conserved amino acids should receive fewer pseudocounts. Using the number of different amino acids in column _c_ , _Rc_ , as an indicator, _Bc_ has been estimated by the formula (Henikoff and Henikoff 1996) 



where _m_ is a positive number derived from trial database searches and _m_ � _m_ � _Bc_ � min ( _m_ � _Nc_ , _m_ /20) (the latter term meaning the minimum of the two given values). By this formula and a given value of _m_ , when _Nc_ � _m_ � 20, the total number of pseudocounts _Bc_ is greater, and when _Nc_ � _m_ � 20 , _Bc_ is smaller than the total number of real counts, _Nc_ , regardless of the value of _Rc_ . The number of pseudocounts is also reduced when _Rc_ �1. In a test search of the SwissProt and Prosite catalogs with various values of _m_ , a value of 5–6 for _m_ produced the most efficient PSSMs for finding known family members. Of the several methods for making PSSMs discussed above, the one with pseudocounts derived by Equations 9 and 10 was most successful. This search was performed with PSSMs derived from blocks with amino acid counts also weighted to account for redundancy (Henikoff and Henikoff 1996). However, pseudocounts added from Dirichlet mixtures, which also 

MULTIPLE SEQUENCE ALIGNMENT I **195** 

vary in each column of the scoring matrix, are also very effective (Henikoff and Henikoff 1996; Tatusov et al. 1997). 

Once pseudocounts have been added to real counts of amino acids in each column of the motif, the PSSM may be calculated. The PSSM has one column (or row) for each position in the motif and one row (or column) for each amino acid, and the entries are log odds entries. Each entry is derived by taking the logarithm to the base 2 (bit units, but sometimes also natural logarithms in nat units are used) of the total of the real counts plus pseudocounts for each amino acid, divided by the probability of that amino acid ( _bca_ / _Nc_ ). An example of a PSSM produced by MEME is shown in Figure 4.15G. 

As a sequence is searched with the PSSM, the value of the first amino acid in the sequence is looked up in the first column of the PSSM, then the value of the second amino acid in the matrix, and so on, until the length scanned is the same as the motif width represented by the matrix. All the log odds scores are added to produce a summed score for start position 1 in the sequence. The process is repeated starting at the second position in the sequence, and so on, until there is not enough sequence left. The highest log odds scoring sequence positions have the closest match statistically to the PSSM. Adding logarithms in this manner is the equivalent of mutiplying the probabilities of the amino acids at each sequence position. To convert each summed log odds score (S) to a likelihood or odds score of the sequence matching the PSSM, use the formula odds score � 2<sup>S</sup> . These odds scores may be summed and each individual score divided by the sum to normalize them and to thereby produce a probability of the motif at each sequence location. 

The above description and example are of using a PSSM to define motifs in protein families. PSSM are also used to define DNA sequence patterns that define regulatory sites, such as promoters or exon–intron junctions in genomic sequences. These topics are discussed in Chapter 8. 

### **Information Content of the PSSM** 

The usefulness of a PSSM in distinguishing real sequence patterns from background may be measured. The unit of measure is the information content in bits. The PSSM described above gives the log odds score for finding a particular matching amino acid in a target sequence corresponding to each motif position. Variations in the scores found in each column of the table are an indication of the amino acid variation in the original training sequences that were used to produce the motif. In some columns, only one amino acid may have been present, whereas in others several may have been present. The columns with highly conserved positions have more information than do the variable columns and will be more definitive for locating matches in target sequences. There is a formal method known as information theory for describing the amount of information in each column that is useful for evaluating each PSSM. The information content of a given amino acid substitution matrix was previously introduced (p. 83) and is discussed in greater detail here. T. Schneider has prepared a Web site that gives excellent tutorials and a review on the topic of information theory, along with methods to produce sequence logos (Schneider and Stephens 1990) at http://www-lmmb.ncifcrf.gov/�toms/sequencelogo. html. 

To illustrate the concepts of information and uncertainty (see above Web site), consider 64 cups in a row with an object hidden under one of them. The goal is to find the object with as few questions as possible. The solution is quite simple. First, ask whether the object is hidden under the first or second half of the cups. If the answer is the first 32, then ask which half of that 32, the first 16 or the second 16, and so on. The sequential questions reduce the possibilities from 64-32-16-8-4-2-1, and six questions will therefore suffice to locate the object. This number is also a measure of the amount of uncertainty in the data 

**196** I CHAPTER 4 

because this number of questions must be asked to find the object. After the first question has been asked, uncertainty has been reduced by 1, so that only five questions then need to be asked to find the object. The uncertainty is zero when the object is found. 

A method to calculate uncertainty (the number of questions to be asked) may be derived from the probability of finding the object under a given cup [ _p_ (object) � 1/64]. Uncertainty is found by taking the negative logarithm to the base 2 of 1/64 [�log2(1/64) � 6 bits]. A situation similar to the hidden object example is found with amino acids in the columns of a PSSM. Here, the interest is to find which amino acid belongs at a particular column in the motif. When we have no information at all, since there are 20 possible choices in all, the amount of uncertainty is log220 � 4.32. 

The data from the PSSM provide information that reduces this uncertainty. If only one amino acid is observed in a column of the PSSM, the uncertainty is zero because there are no other possibilities. If two amino acids are observed with equal frequency, there is still uncertainty as to which one it is, and one question must be asked to find the answer, or uncertainty � 1. The formula for finding the uncertainty in this example is the sum of the fractional information provided by each amino acid, or � [0.5 � log20.5 � 0.5 � log20.5] � 1. In general, the average amount of uncertainty ( _Hc_ ) in bits per symbol for column c of the PSSM is given by 



where _pic_ is the frequency of amino acid _i_ in column _c_ and is estimated by the frequency of occurrence of each amino acid ( _bca_ / _Nc_ ) and log2( _pic_ ) is the log odds score for each amino acid in column _c_ . Uncertainty for the entire PSSM may then be calculated as 



_H_ is also known as the entropy of the PSSM position in information theory because the higher the value, the greater the uncertainty. The lower the value of the uncertainty _H_ for the PSSM, the greater the ability of the PSSM to distinguish real occurrences of the motif from random matches. Conversely, the higher the information content, calculated as shown below, the more useful the PSSM. 

### **Sequence Logos** 

Sequence logos are graphs that illustrate the amount of information in each column of a motif. The logo is derived from sequence information in the PSSM described above. Conserved patterns in both protein and DNA sequences can be represented by sequence logos. A program for producing logos, along with several examples, is available from http://wwwlmmb.ncifcrf.gov/�toms/sequencelogo.html. The Web site of S.E. Brenner at http://www.bio.cam.ac.uk/seqlogo/ will produce sequence logos from an input alignment using the Gibbs sampler method, and an implementation of an extension of the logo method for structural RNA alignment (Gorodkin et al. 1997) is at http://www.cbs.dtu.dk/ gorodkin/appl/plogo.html. A logo representation for the BLOCKS database has been implemented (Henikoff et al. 1995) and may be viewed when the information on a partic- 

MULTIPLE SEQUENCE ALIGNMENT I **197** 

ular block is retrieved from the BLOCKS Web server (http://www.blocks.fhcrc.org/). An example of a Block logo is shown in Figure 4.19. Another example of a simple graph of information content is given in Figure 4.15C. In this case, the information for the entire motif has been calculated by the MEME server by summing the values in each column to a total value of 22 bits. Although logos are primarily used with ungapped motifs and sequence patterns, logos of alignments that include gaps in some sequence positions may also be made. If such is the case, then the height of the column with gaps is reduced by the proportion of sequence positions that are not gaps. 



<!-- Start of picture text -->
4<br>a2<br>A)<br>KAMFNOKOMHOKNM<br>ee eeTHORee<br>Logo of Gibbs Block D (Tcl) 9 sequences<br><!-- End of picture text -->

**Figure 4.19.** A sequence logo. The logo represents the amount of information in each column of a motif corresponding to the values in PSSM of the motif discussed above. The horizontal scale represents sequential positions in the motif. The height of each column gives the decrease in uncertainty provided by the information in that column. The higher the column, the more useful that position for finding matches in sequences. In each column are shown symbols of the amino acids found at the corresponding position in the motif, with the height of the amino acid proportional to the frequency of that amino acid in the column, and the amino acids shown in decreasing order of abundance from the top of the column. From each logo, the following information may thus be found: The consensus may be read across the columns as the top amino acid in each column, the relative frequency of each amino acid in each column of the motif is given by the size of the letters in each column, and the total height of the column provides a measure of how useful that column is for reducing the level of uncertainty in a sequence matching experiment. Note that the highest values are for columns with less diversity. 

**198** I CHAPTER 4 

The height of each logo position is calculated as the amount by which uncertainty has been decreased by the available data; in this case, the amino acid frequencies in each column of the motif. The relative heights of each amino acid within each column are calculated by determining how much each amino acid has contributed to that decrease. The uncertainty at column _c_ is given by Equation 11. Because the maximum uncertainty at a position/column when no information is available is log220 � 4.32, as more information about the motif is obtained by new data, the decrease in uncertainty (or increase in the amount of information) _Rc_ is 



where _Hc_ is given by Equation 11 and � _n_ is a correction factor for a small sequence number _n_ . _Rc_ is used as the total height of the logo column. The height of amino acid _a_ at position _c_ in the motif logo is then given by _fac_ � _Rc_ . 

The above description applies to protein sequences. Sequence logos are also produced for DNA sequences. The methodology is very similar to the above except that there are only four possible choices for each logo location. Hence, the maximum amount of uncertainty is log24 � 2. The above method assumes that the sequence pattern is less random than the background or expected sequence variation, and this assumption limits the ability of the method to locate subtle patterns in sequences. 

An improved method for finding more subtle patterns in sequences is called the relative entropy method (Durbin et al. 1998). In this case, differences between the observed frequencies and background frequencies are used (Gorodkin et al. 1997), and the decrease in uncertainty from background to observed (or amount of information) in bits is given by 



where _bi_ is the background frequency of residue _i_ in the organism and the maximum uncertainty in column _c_ is given by ��all _i_ [ _pic_ log2(1/ _bi_ )]. When background frequencies are taken into account, and the column frequency is less than the background frequency, it is possible for the information given by a particular residue in a logo column to be negative. To accommodate this change, the corresponding sequence character is inverted in the logo to indicate a less than expected frequency. There are also two ways used to illustrate the contribution of each character through the height of the symbol. The first method is described above. The second method is to display symbol heights in proportion to the ratio of the observed to the expected frequency, i.e., by the fraction ( _pic_ / _bi_ ) / (�all _i pic_ / _bi_ ) for each symbol _i_ . Gaps are included in the analysis by using _p_ gap � 1 and, as a result, will always give a negative contribution to the information (Gorodkin et al. 1997). 

### **~~MULTIPLE SEQUENCE ALIGNMENT EDITORS AND FORMATTERS~~** 

Once a multiple sequence alignment has been obtained by the global msa program, it may be necessary to edit the sequence manually to obtain a more reasonable or expected alignment. Several considerations must be kept in mind when choosing a sequence editor, which should include as many of the following features as possible: (1) provision for displaying the sequence on a color monitor with residue colors to aid in a clear visual repre- 

MULTIPLE SEQUENCE ALIGNMENT I **199** 

sentation of the alignment, (2) recognition of the multiple sequence format that was output by the msa program and maintenance of the alignment in a suitable format when the editing is completed, (3) provision of a suitable windows interface, allowing use of the mouse to add, delete, or move sequence followed by an updated display of the alignment. In addition, there are other types of editing that are commonly performed on msas such as, for example, shading conserved residues in the alignment. 

The large number of multiple sequence alignment formats that are in use were discussed in Chapter 2. Two commonly encountered examples are the Genetics Computer Group’s MSF format and the CLUSTALW ALN format. Because these formats follow a precise outline, one may be readily converted to another by computer programs. READSEQ by D.G. Gilbert at Indiana University at Bloomington is one such program. This program will run on almost any computer platform and may be obtained by anonymous FTP from ftp.bio.indiana.edu/molbio/readseq. There is also a Web-based interface for READSEQ from Baylor College of Medicine at http://dot.imgen.bcm.tmc.edu:9331/seq-util/sequtil.html/. A software package SEQIO, which provides C program modules for conversion of sequence files from one format to another, is available by anonymous FTP from ftp.pasteur.fr/pub/GenSoft/unix/programming/seqio-1.2.tar.gz; documentation is available at http://bioweb.pasteur.fr/docs/doc-gensoft/seqio/. 

A short list of the many available programs that have or exceed the above-listed features is discussed below. For a more comprehensive list, visit the catalog of software page at Web address http://www.ebi.ac.uk/biocat/. 

### **Sequence Editors** 

1. CINEMA (Colour Interactive Editor for Multiple Alignments) at http://www.biochem. ucl.ac.uk/bsm/dbbrowser/CINEMA2.02/kit.html is a broadly functional program for sequence editing and analysis, including dot matrix analysis. It features drag-and-drop editing, sequence shifting to left or right, viewing of different parts of an alignment using the split-screen option, multiple motif selection and manipulation, and a number of added features such as viewing of protein structures. CINEMA was developed by A.W.R. Payne, D.J. Parry-Smith, A.D. Michie, and T.K. Attwood. CINEMA is an applet that runs under a Web browser and therefore will run on almost any computer platform. 

2. GDE (Genetic Data Environment) provides a general interface on UNIX machines for sequence analysis, sequence alignment editing, and display (Smith et al. 1994) and is available from several anonymous FTP sites including ftp.ebi.ac.uk/pub/software/unix. GDE is described at http://bimas.dcrt.nih.gov/gde_sw.html, and http://www.tigr.org/ �jeisen/GDE/GDE.html. GDE features are incorporated into the Seqlab interface for the GCG software, vers. 9. This interface requires communication with a host UNIX machine running the Genetics Computer Group software. Interface with MS-DOS or Macintosh is possible if the computer is equipped with the appropriate X-Windows client software. 

3. GeneDoc is an alignment editing and display editor by K. Nicholas and H. Nicholas of the Pittsburgh Supercomputing Center for MSF-formatted msas. It can also import files in other formats. GeneDoc can move residues by inserting or deleting gap, and features drag-and-drop editing. As the alignment is edited, a new alignment score is calculated by sum of pairs method or based on a phylogenetic tree. GeneDoc is available from http://www.psc.edu/biomed/genedoc/ and runs under MS Windows. 

4. MACAW is both a local multiple sequence alignment program and a sequence editing tool (Schuler et al. 1991). Given a set of sequences, the program finds ungapped blocks in the sequences and gives their statistical significance. Later versions of the program 

