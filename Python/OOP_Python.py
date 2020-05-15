## RENAME this file YourLastName_OOP_FinalProject_2020.py

##Assignment: Add to the constructor and methods of a parent class and child classes
##            which inherit the base class properties.

## Begin with the parent Seq class and the child DNA class we created in lecture below.
## 

### Seq Class
#
#  Constructor:
#  (1) Use the string functions upper and string to clean up self.sequence.
#  (2) Add a variable self.kmers to the constructor and make it equal to an empty list.

#  Methods:
#  (1) Add a method called make_kmers that makes kmers of a given length from self.sequence
#      appends these to self.kmers. Default kmer parameter=3.
#  (2) Add a method called fasta that returns a fasta formatted string like this:
#      >species gene
#      AGATTGATAGATAGATAT


### DNA Class: INHERITS Seq class
#   
#  Constructor:
#  Use re.sub to change any non nucleotide characters in self.sequence into an 'N'.
#      re.sub('[^ATGCU]','N',sequence) will change any character that is not a
#      capital A, T, G, C or U into an N. (Seq already uppercases and strips.)

#  Methods:
#  Add a method called print_info that is like print_record, but adds geneid and an
#      empty space to the beginning of the string.

### RNA Class:  INHERITS DNA class
#  
#  Construtor:
#  Use the super() function (see DNA Class example).
#  (1) Automatically change all Ts to Us in self.sequence. 
#  (2) Add self.codons equals to an empty list

#  Methods:
#  (1) Add make_codons which breaks the self.sequence into 3 letter codons
#      and appends these codons to self.codons unless they are less than 3 letters long.
#  (2) Add translate which uses the Global Variable standard_code below to
#      translate the codons in self.codons and returns a protein sequence.

### Protein Class: INHERITS Seq class
#
#  Construtor:
#  Use the super() function (see DNA Class example).
#  (1) Use re.sub to change any non LETTER characters in self.sequence into an 'X'.
#  (2) self.aa_counts has already been added - it is a dictionary with all 21
#      amino acids as keys and values set to 0

#  Methods:
#  (1) Add tabulate_amino_acids, which counts the amino acids in self.sequence
#      and puts these in the dictionary self.counts
#  The next 2 methods use the kyte_doolittle dictionary you will create below.
#  (2) Add total_hydro, which return the sum of the total hydrophobicity of a sequence
#  (3) Add hydro_scan, which returns list of hydrophobicity values for self.kmers
#      Each element of the list will be the average for the kmer (sum/total length of the kmer)
#      See https://kelleybioinfo.org/algorithms/default.php?o=2 for details on how to do this.
#      If self.kmers is empty [], use make_kmers(5) to make a kmer list.

import re

standard_code = {
     "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "UCU": "S",
     "UCC": "S", "UCA": "S", "UCG": "S", "UAU": "Y", "UAC": "Y",
     "UAA": "*", "UAG": "*", "UGA": "*", "UGU": "C", "UGC": "C",
     "UGG": "W", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
     "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAU": "H",
     "CAC": "H", "CAA": "Q", "CAG": "Q", "CGU": "R", "CGC": "R",
     "CGA": "R", "CGG": "R", "AUU": "I", "AUC": "I", "AUA": "I",
     "AUG": "M", "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
     "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGU": "S",
     "AGC": "S", "AGA": "R", "AGG": "R", "GUU": "V", "GUC": "V",
     "GUA": "V", "GUG": "V", "GCU": "A", "GCC": "A", "GCA": "A",
     "GCG": "A", "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
     "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

#Make a dictionary called kyte_doolittle which has the hydrophobicity values for
# each amino acid from the Kyte-Doolittle scale: https://kelleybioinfo.org/algorithms/default.php?o=2
# For unknown ('X') amino acids, the value should be zero.

kyte_doolittle = {
"A":1.8, "L":3.8,"R":-4.5,"K":-3.9,"N":-3.5,"M":1.9,
"D":-3.5,"F":2.8,"C":2.5,"P":-1.6,"Q":-3.5,"S":-0.8,
"E":-3.5,"T":-0.7,"G":-0.4,"W":-0.9,"H":-3.2,"Y":-1.3,
"I":4.5,"V":4.2}  

class Seq:

    def __init__(self,sequence,gene,species):
        self.sequence=str(sequence.upper().strip())
        self.gene=gene
        self.species=species
        self.kmers = []

    def __str__(self):
        return self.sequence

    def print_record(self):
        print(self.species + " " + self.gene + ": " + self.sequence)

    def make_kmers(self, k=3):
        num_kmers= len(self.sequence) - k +1
        for i in range(num_kmers):
            self.kmers.append(self.sequence[i:i+k])
        return self.kmers

    def fasta(self):
        return ">"+self.species+" "+self.gene+"\n"+self.sequence
    
class DNA(Seq):

    def __init__(self,sequence,gene,species,geneid,**kwargs):
        super().__init__(sequence,gene,species)
        self.sequence=re.sub('[^ATGCU]','N',self.sequence)
        self.geneid=geneid
 
    def analysis(self):
        gc=len(re.findall('G',self.sequence) + re.findall('C',self.sequence))
        return gc

    def print_info(self):
        print(self.geneid+" "+self.species + " " + self.gene + ": " + self.sequence)

class RNA(DNA):

    def __init__(self,sequence,gene,species,geneid):
        super().__init__(sequence,gene,species,geneid)
        self.sequence=re.sub('[T]','U',self.sequence)
        self.codons=[]
        
    def make_codons(self):
        protein_seq=[]
        for n in range(0, len(self.sequence), 3):
            codon_seq = self.sequence[n:n+3]
            if len(codon_seq) == 3 :
                self.codons.append(codon_seq)
        return self.codons
    
    def translate(self):
        protein_seq=""
        for x in self.codons:
            try:
                protein_seq+=standard_code[x]
            except:
                protein_seq+='X'

        return protein_seq


class Protein(Seq):

    def __init__(self,sequence,gene,species,geneid):
        super().__init__(sequence,gene,species)
        self.sequence=re.sub('[^a-zA-Z]','X',self.sequence)
        self.aa_counts={'A':0,'C':0,'D':0,'E':0,'F':0,'G':0,'H':0,'I':0,'K':0,'L':0,'M':0,'P':0,'Q':0,'R':0,'S':0,'T':0,'U':0,'V':0,'W':0,'X':0,'Y' :0}

    def tabulate_amino_acids(self):
        for i in self.sequence:
            try:
                self.aa_counts[i]+=1
            except:
                self.aa_counts[i]=1
        return self.aa_counts

    def total_hydro(self):
        total_hydro=0
        for i in self.sequence:
            try:
                h_score = kyte_doolittle[i]
            except:
                h_score=0
            total_hydro = total_hydro + h_score
        return total_hydro

    def hydro_scan(self):
        hydroscan=[]
       
        if not self.kmers:
            self.make_kmers(5)
        
        for kmer in self.kmers:
            total_hydro=0
            for i in kmer:
                try:
                    h_score = kyte_doolittle[i]
                except:
                    h_score=0
                total_hydro = total_hydro + h_score  
            hydroscan.append((total_hydro/len(kmer)))  
        return hydroscan



    







