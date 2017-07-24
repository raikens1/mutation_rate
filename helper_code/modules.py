import os
import string
import copy
import subprocess
import math
import random
import pdb
from function_wrapper import *

class transition_mat:
    #It saves the information regarding the transition matrix of the kmers
    #Also saves the information regarding the count of particular kmers 
    #Public variables a)kmer (kmer like 3 or 4 or 5) b) mat(count of the kmer transitions) c) kmerseq(the list of possible kmers) 
    #d) kmerseq_count (dictionary which saves the count of a particular kmer)
    #e) each kmer is (before,change,after) i.e. if before is 2 and after is 1 then we consider changes like XYMZ to XYmZ only	
    def __init__(self,before,after): #before and after have to be >=1
	self.kmer = before+after+1
	self.before = before
	self.after = after
	ntd = ["A","C","G","T"]
	self.ntd = ntd
	keys = ntd
	k = before+after+1
	while(k>1): #This generates the possible kmers 
	    temp_keys=[]
	    k = k-1
	    for i in range(len(keys)):
	        for j in range(len(ntd)):
		    key = keys[i]+ntd[j] #Concatenate the previous key to nucleotide
		    temp_keys.append(key)
	    keys=temp_keys
	self.kmerseq = keys
	kmerseq_count = {}
	for i in self.kmerseq: #The count of kmerseq has been intialized
	    kmerseq_count[i] = 0
	self.kmerseq_count = kmerseq_count
	mat = {}
	for i in self.kmerseq:
	    for j in ntd:
		alpha = i
		beta = alpha[:self.before] + j + alpha[self.before+1:]
		orig = alpha[self.before]
		if not(orig==j):     
		    mat[(alpha,beta)] = {'total':0}    #Each transition is a different 
	self.mat = mat

    def increment(self,key,annot): #Increments the mat for a particular key (Kmer going to another kmer for a particular annotation)
	if key in self.mat:
	    if annot in self.mat[key]:
	        self.mat[key][annot] += 1
	    else:
		self.mat[key][annot] = 1
	    self.mat[key]['total'] += 1
    def kmer_increment(self,seq):
	#Increment the all possible kmers in a string. 123, 234, 345, 456
	#Sequence is in the list format
	for i in range(len(seq)-self.kmer+1): 
	    kmer_str = seq[i: (i+self.kmer)]
	    kmer_str = ''.join(kmer_str)
	    if kmer_str in self.kmerseq:
	        self.kmerseq_count[kmer_str] +=1

class transition_mat_exons:
    #It saves the information regarding the transition matrix of the codons in the exons
    #Also saves the information regarding the count of codons 
    #This differs from the previous because here all the transitions of a codon are allowed. Not just the middle position changes.
    def __init__(self):
        kmer = 3
        before = 1
        after = 1
        ntd = ["A","C","G","T"]
        keys = ntd
        k = before+after+1
        while(k>1): #This generates the possible kmers 
            temp_keys=[]
            k = k-1
            for i in range(len(keys)):
                for j in range(len(ntd)):
                    key = keys[i]+ntd[j] #Concatenate the previous key to nucleotide
                    temp_keys.append(key)
            keys=temp_keys
        self.kmerseq = keys
        kmerseq_count = {}
        for i in self.kmerseq: #The count of kmerseq has been intialized
            kmerseq_count[i] = 0
        self.kmerseq_count = kmerseq_count
        mat = {}
        for alpha in self.kmerseq:
            for beta in self.kmerseq:
                if not(alpha==beta):
                    mat[(alpha,beta)] = {'total':0}    #Each transition is different 
        self.mat = mat
    def increment(self,key,annot): #Increments the mat for a particular key (Kmer going to another kmer for a particular annotation)
        if key in self.mat:
            if annot in self.mat[key]:
                self.mat[key][annot] += 1
            else:
                self.mat[key][annot] = 1
            self.mat[key]['total'] += 1
    def kmer_increment_individual(self,kmer_str,count):
	#Increments the count of kmer_str by count
	self.kmerseq_count[kmer_str] += count

class file_handles:
    #This class fetches a before+after  sequence around a variant
    #This class also fetches a fasta sequence between two positions on a chromosome
    #Public varianbles a)chr_name (The allowed chromosome names) b) handles (The dictionary containing the handles of all the fasta files
    def __init__(self,directory):
	#k is the kmer length which calculates the str_len
    	self.chr_name = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]
	handles = {}
	for i in self.chr_name:
	    name = "handle"+i
	    filename = directory + "/chr" + i + ".fa"
	    handles[i] = file(filename,"r")
	self.handles = handles
    def pos_fasta(self,pos,init_offset): #This is an internal function which gives the coordinate in the fasta file corresponding to the position of the variant. It also needs the length of the first line (header of the fasta file)
	#init_offset is the length of the first line in the fasta file
	val = init_offset +  51*(pos/50) + (pos%50) - 1
	if(pos%50 == 0): #End of line of fasta file character so this special -1
	    val = val - 1
	return(val)
    def fetch_str_var(self,chrom,pos,before,after):
	#This fetches the before+pos+after length sequence around a particular variant on a chromsome
	#Chrom should be a char type Pos should be a integer
	if chrom in self.chr_name:
	    cur_handle = self.handles[chrom]
  	    #size of the first line in the fasta file
	    cur_handle.seek(0) 
	    first_line_offset = len(cur_handle.readline())
            offsets = range(-1*before,after+1) #We need before+after+1 length around a variant and the range function needs to be adjusted accordingly
	    str_read = []
	    for i in offsets:
	        cur_handle.seek( self.pos_fasta(pos+i,first_line_offset) )
	        str_read.append(string.upper(cur_handle.read(1)))
	    return(str_read)
	else:
	    print "Invalid Chromosome \n"
    def get_string(self,chrom,start,end,before,after): # returns a string from chromosome between start and end position
	#Chrom should be a string
        try:
	 if chrom in self.chr_name:
	    cur_handle = self.handles[chrom]
	    cur_handle.seek(0)
	    first_line_offset = len(cur_handle.readline())
	    k = before+after+1
	    str_read = []
	    for i in range(start-before, end+after+1): 
		cur_handle.seek( self.pos_fasta(i,first_line_offset) )
		str_read.append(string.upper(cur_handle.read(1)))
	 return(str_read)
	except:
	 print "Bad_Here"
	 pdb.set_trace()

class exon_file_properties:
#This is a class with has the information about the codons.
#It gets all the info from the fasta file. Public variables are pos_codon which maps each position in the sequence to a codon_index and its base in it.
#dictionary index_codon maps each codon_index to a codon sequence. nuc_amino maps amino_sequence to name
#Nuc_amino maps the codon_nucleotide code to amino acid name.
    def __init__(self,filename_fasta,reverse_or_not=False):
        #Add now. If reverse_or_not is 'reverse' then we make the index from the other hand. We did the same in save_fasta_seq. The seq has been saved approrpriatelt. Onlt need to correct this.
	self.filename = filename_fasta
	self.pos_codon = {}   #Dictionary in which each genomic position is mapped to (codon_index,position in that codon) or (splice, splice)
	self.index_codon = {} #Dictionary in which each codon index is mapped to its 3letter codon word
	self.pos_sequence = {} #Each index is mapped to 7mer context around it
	self.nuc_amino = {}
	self.reverse_or_not = reverse_or_not
	handle = file("/project/voight_subrate/avarun/Research/mutation_rate/codon_aa_table","r")
        content = handle.readlines()
	for entry in content:
	    entry = entry.rstrip('\n').split('\t')
	    acids = entry[1].split(',')
	    for aa in acids:
		self.nuc_amino[aa] = entry[0]
	self.total_codon = 0
    def update_pos_index(self): #This function reads the fasta sequence and updates the data structures
	handle = file(self.filename,"r")
	content = handle.readlines()
	final_sequence = []
	codon_no = 1
	codon_pos_in = 1
	if self.reverse_or_not == False:
	    exon_num = 1
	    for entry in content:
	        entry = entry.rstrip('\n').split('\t')
		if exon_num == 1:
		  if len(content)>1:
  		    sequence = entry[3][3:-5]
		    start = int(entry[1])
		    end = int(entry[2])-2
		    self.pos_codon[int(entry[2])-1] = ('splice_donor','splice')
                    self.pos_codon[int(entry[2])] = ('splice_donor','splice')
		  else:
                    sequence = entry[3][3:-3]
                    start = int(entry[1])
                    end = int(entry[2])
		if exon_num >1 and exon_num == len(content):
		    sequence = entry[3][5:-3]
		    start = int(entry[1])+2
		    end = int(entry[2])
                    self.pos_codon[int(entry[1])+1] = ('splice_acceptor','splice')
                    self.pos_codon[int(entry[1])] = ('splice_acceptor','splice')
		if not exon_num == 1 and not exon_num == len(content):
		    sequence = entry[3][5:-5]
		    start = int(entry[1])+2
		    end = int(entry[2])-2
                    self.pos_codon[int(entry[1])+1] = ('splice_acceptor','splice')
                    self.pos_codon[int(entry[1])] = ('splice_acceptor','splice')
                    self.pos_codon[int(entry[2])-1] = ('splice_donor','splice')
                    self.pos_codon[int(entry[2])] = ('splice_donor','splice')
		exon_num += 1
		final_sequence += sequence
	        for index in range(start,end+1): #Because range goes till upper_limit-1
		    self.pos_codon[index] = (codon_no,codon_pos_in)
		    if (codon_pos_in ==3):
		        codon_pos_in = 1
		        codon_no += 1
		    else:
		        codon_pos_in += 1
		counter = 0
		for index in range(int(entry[1]),int(entry[2])+1):
		    self.pos_sequence[index] = entry[3][counter:counter+7]
		    counter += 1		    
        if self.reverse_or_not == 'reverse':
            exon_num = 1
            for entry in content:
                entry = entry.rstrip('\n').split('\t')
                if exon_num == 1:
		  if len(content)>1:
                    sequence = entry[3][3:-5]
                    start = int(entry[1])+2
                    end = int(entry[2])
                    self.pos_codon[int(entry[1])+1] = ('splice_donor','splice')
                    self.pos_codon[int(entry[1])] = ('splice_donor','splice')
		  else:
                    sequence = entry[3][3:-3]
                    start = int(entry[1])
                    end = int(entry[2])
                if exon_num >1 and exon_num == len(content):
                    sequence = entry[3][5:-3]
                    start = int(entry[1])
                    end = int(entry[2])-2
                    self.pos_codon[int(entry[2])-1] = ('splice_acceptor','splice')
                    self.pos_codon[int(entry[2])] = ('splice_acceptor','splice')
                if not exon_num == 1 and not exon_num == len(content):
                    sequence = entry[3][5:-5]
                    start = int(entry[1])+2
                    end = int(entry[2])-2
                    self.pos_codon[int(entry[1])+1] = ('splice_acceptor','splice')
                    self.pos_codon[int(entry[1])] = ('splice_acceptor','splice')
                    self.pos_codon[int(entry[2])-1] = ('splice_donor','splice')
                    self.pos_codon[int(entry[2])] = ('splice_donor','splice')
                exon_num += 1
                final_sequence += sequence
		index = end
		while(index>=start):
                    self.pos_codon[index] = (codon_no,codon_pos_in)
                    if (codon_pos_in ==3):
                        codon_pos_in = 1
                        codon_no += 1
                    else:
                        codon_pos_in += 1
		    index -= 1
                counter = 0
		index = int(entry[2])
		while(index>=int(entry[1])):
                    self.pos_sequence[index] = entry[3][counter:counter+7]
                    counter += 1
		    index -= 1
	total_codon = len(final_sequence)/3
	self.total_codon = total_codon
	for index in range(1,total_codon+1):
	     self.index_codon[index] = [final_sequence[index*3-3],final_sequence[index*3-2],final_sequence[index*3-1]]	

class total_length_per_chr:
    #USGAE: MAKE SURE THAT save_fasta_sequence is called before hand.
    #This files takes an interval file and calculates the length of each chromosome present in the file.
    #It also saves the fasta sequence for each interval in a file. And then it updates the kmer count for the objects
    #Public variables a)chr_len (Directory containing the length for each chromosome) b) filename(the address of the filename containing the intervals) c) chr_name (list of acceptable chromosomes)
    def __init__(self,filename):
	self.filename = filename
	self.filename_seq = filename+"_fasta_sequence"
        self.chr_name = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
	chr_len = {}
	for i in self.chr_name:
	    chr_len[i] = 0	
	self.chr_len = chr_len
	self.total_len = 0.0   #Total length across all chromosomes (including the X chromosome) 
    def find_len_per_chr(self): #Find the length of all the chromosomes 
	handle = file(self.filename,"r")
	content = handle.readlines()
	for entry in content:
	    entry = entry.rstrip('\n').split('\t')
	    if entry[0] in self.chr_name:
		self.chr_len[entry[0]] = self.chr_len[entry[0]] + int(entry[2]) - int(entry[1]) +1
		self.total_len = self.total_len+int(entry[2]) - int(entry[1])+1
	print "total_length across all chromosomes in KBs " + str(self.total_len/1000)
    def update_kmer_count(self,dir_object_trans,before,after,strand=False,offset=False): #object_trans is the the directory containing objects of the class trans_mat
        #Updates the directory of the objects of class trans_mat with all possible kmers
	#If the chromosome did not have the object for class transition_mat then creates one too
	#If offset = True (or a number), then we update the kmer count by starting from a specific position and then increasing the scale by 3
	val = subprocess.call(["ls",self.filename_seq],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	if not(val==0):
	    print "Something is wrong. You have not yet saved the fasta sequence anywhere"+"\n"
	else:
	    handle = file(self.filename_seq,"r")
	    content = handle.readline()   #While saving the sequence we made sure that we only save for the right chromosome. So no need to check its validity
	    if offset==False: #Normal stuff... However if offset = True like in coding regions where I want to save the kmer based on middle position in the sequence. There I do different stuff
	      while(content):
	    	content = content.rstrip('\n').split('\t')
		if (strand == False):
	            seq = list(content[3])   #Because entry 3 is the sequence. Else entry 4 is the sequence 
		else:
		    seq = list(content[4])	    
		#######################################################NEW STUFF#############################################################################
		#ASSUME SEQ HAS BEEN CALLED +3 on either sides 
		#SO if I find kmer counts for 3mer or 5mer or 1mer or null then I need to not consider the sequence ends
		if before <3:
		    seq = seq[3-before:-1*(3-before)]
		#############################################################################################################################################
		if not(content[0] in dir_object_trans):
		    dir_object_trans[content[0]] = transition_mat(before,after)
	        dir_object_trans[content[0]].kmer_increment(seq)
	        content = handle.readline()
	      handle.close()
	    else:
	      #We need to save the sequence and then 
	      offset = offset - 1 #offset should be substracted 
	      seq = ''
	      while(content):	   
		content = content.rstrip('\n').split('\t')
		chr_name = content[0]
                if (strand == False):
                    seq = seq+(content[3].rstrip('\n'))   #Because entry 3 is the sequence. Else entry 4 is the sequence 
                else:
                    seq = seq+(content[4].rstrip('\n'))
                if not (content[0] in dir_object_trans):
                    dir_object_trans[chr_name] = transition_mat(before,after)
                content = handle.readline()
              handle.close()
	      seq = list(seq)
	      kmer = before+after+1	
	      #Now increment the dir_object_trans[chr_name]							    
	      for i in range(len(seq)-kmer+1):
                  kmer_str = seq[i: (i+kmer)]
		  if (i+(kmer/2))%3 == offset:
                      kmer_str = ''.join(kmer_str)
                      if kmer_str in dir_object_trans[chr_name].kmerseq:
                	dir_object_trans[chr_name].kmerseq_count[kmer_str] +=1
		      else:
			print "soemething weird as kmer not in object"
	#		pdb.set_trace()

          
    def save_fasta_sequence(self,object_files,before,after,strand=False,do_reverse_or_not=False): #This will save the fasta sequence for the file in a new file called (filename+fasta). Each line will have the coordinates and fasta sequence for that interval
	#object_files is the object of the class files_handles
        #Do_reverse_or_not is only needed because I forgot to save the strand information in coding data. If do_reverse_or_not='reverse' then we save the reverse_complement of seq

	handle = file(self.filename,"r")
	val = subprocess.call(["ls",self.filename_seq],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	if(val == 0):
	    pass
	else:
	    handle_w = file(self.filename_seq,"w") 
	    content = handle.readlines()
	    for entry in content:
                entry = entry.rstrip('\n').split('\t')
                if entry[0] in self.chr_name:
                    seq = ''.join(object_files.get_string(entry[0], int(entry[1]), int(entry[2]),before,after))
		    if strand == False:
			if do_reverse_or_not == 'reverse':
		             handle_w.write(entry[0]+'\t'+entry[1]+'\t'+entry[2]+'\t'+reverse_complement(seq)+'\n')       
			if do_reverse_or_not == False:
                             handle_w.write(entry[0]+'\t'+entry[1]+'\t'+entry[2]+'\t'+seq+'\n')			
		    if strand == True:
			if entry[3] == '-':
			    handle_w.write(entry[0]+'\t'+entry[1]+'\t'+entry[2]+'\t'+entry[3]+'\t'+reverse_complement(seq)+'\n')		        		
			if entry[3] == '+':
                            handle_w.write(entry[0]+'\t'+entry[1]+'\t'+entry[2]+'\t'+entry[3]+'\t'+seq+'\n')                                			
   	    handle_w.close()
        handle.close()
	
def analyze_var_file(filename,object,object_fa,before,after,filter="FALSE",filter_freq=0.0,file_exons=False,reverse_or_not=False): 
#This function opens the annotated variant file and updates the transition_mat object. 
#Object is the directory for the objects of the class transition_mat
#Object_fa which is the object for the class file_handles
#file_exons is the object of file exon_file_properties. It is used to get the offset of each position. Aka if a position is at the 1st, 2nd or 3rd position of the codon
#If reverse_or_not='reverse' then we mutated stuff is the reverse complement of original stuff
    handle = file(filename,"r")
    content = handle.readlines()
    for entry in content:
       entry = entry.rstrip("\n").split("\t")
       freq = float(entry[5])
       flag_con = "TRUE"
       if(filter=="TRUE"):
	    if(freq>=filter_freq and freq<=(1-filter_freq)):
		flag_con = "TRUE"
	    else:
		flag_con = "FALSE"
       if(flag_con=="TRUE"):		
	orig = object_fa.fetch_str_var(entry[0],int(entry[1]),before,after) #Extract before+1+after sequence around the variant
	anc = entry[2]  
	der = entry[3]  #The mutated nucleotide
	if reverse_or_not == 'reverse':
	    anc = reverse_complement(anc)
	    der = reverse_complement(der)
	    orig = reverse_complement(orig)
	    orig = list(orig)
	if(anc == der):
	   print "something_wrong as ancestral is same as derived"
	   print entry
	mutated = copy.deepcopy(orig)
	mutated_pos = before #because in python numbering starts from 0
	orig[mutated_pos] = anc
	mutated[mutated_pos] = der
	key = (''.join(orig),''.join(mutated)) 
	annot = entry[4]

	#Check if we also want to increment based on the position in the exons
	if file_exons==False:
	    if not(entry[0] in object):
	        object[entry[0]] = transition_mat(before,after)
            object[entry[0]].increment(key,annot)
	else:
	    offset = str(file_exons.pos_codon[int(entry[1])][1])
            if not(entry[0] in object[offset]):
                object[offset][entry[0]] = transition_mat(before,after)
            object[offset][entry[0]].increment(key,annot)	
    handle.close()

def analyze_var_file_ancestral(filename,object,object_fa,before,after):
#This function opens the annotated variant file and updates the transition_mat object. 
#Object is the directory for the objects of the class transition_mat
#Object_fa which is the object for the class file_handles
    handle = file(filename,"r")
    content = handle.readlines()
    for entry in content:
        entry = entry.rstrip("\n").split("\t")
        orig = object_fa.fetch_str_var(entry[0],int(entry[1]),before,after) #Extract before+1+after sequence around the variant
        ref = entry[2]
        alt = entry[3]  #The mutated nucleotide
        anc = entry[6].upper()
        if anc == 'A' or anc == 'G' or anc == 'C' or anc == 'T':
	    if anc == alt: #Ancestral does not match reference (so change the ref and alt) 
	        temp = ref
		ref = alt
		alt = temp
        mutated = copy.deepcopy(orig)
        mutated_pos = before #because in python numbering starts from 0
        orig[mutated_pos] = ref
        mutated[mutated_pos] = alt
        key = (''.join(orig),''.join(mutated))
        annot = entry[4]
        #Check if we also want to increment based on the position in the exons
        if not(entry[0] in object):
            object[entry[0]] = transition_mat(before,after)
        object[entry[0]].increment(key,annot)
    handle.close()


def analyze_var_exons(filename,object_exon_file,object,freq_flag="FALSE",filename_freq="FALSE",filter="FALSE",filter_freq=0.0,reverse_or_not=False):
#This function opens the annotated variant file and updates the transition_mat_exons object and also the count of each codon that is seen
#Object is an object of the class transition_mat
#Object_exon_file is the object of the class exon_file_properties
#It also generates a new file by the name of filename_freq which contains the amino acid change and its frequency
#If reverse_or_not='reverse' then we mutated stuff is the reverse complement of original stuff
    handle = file(filename,"r")
    content = handle.readlines()
    flag_con = "TRUE"
    if(freq_flag=="TRUE"):
	handle_w = file(filename_freq,"w")
    for entry in content:
     entry = entry.rstrip('\n').split('\t')
     if(filter=="TRUE"):
	flag_con = "FALSE"
     if (float(entry[5])>=filter_freq and float(entry[5])<=(1.0-filter_freq)):
        flag_con = "TRUE"
     if(flag_con == "TRUE"):
	codon_index = object_exon_file.pos_codon[int(entry[1])][0]
	codon_inside = object_exon_file.pos_codon[int(entry[1])][1]
	if(codon_index<= object_exon_file.total_codon):
	 codon_seq = copy.deepcopy(object_exon_file.index_codon[codon_index])
	 anc = entry[2]
	 der = entry[3]
         if reverse_or_not == 'reverse':
            anc = reverse_complement(anc)
            der = reverse_complement(der)
	 codon_seq[codon_inside-1] = anc
	 mutated = copy.deepcopy(codon_seq)
	 mutated[codon_inside-1] = der
         key = (''.join(codon_seq),''.join(mutated))
         annot = entry[4]
         object.increment(key,annot)
	 if (freq_flag=='TRUE'):
	    str_write = ''.join(codon_seq) + '\t' + ''.join(mutated) + '\t' + entry[5] +'\n'
	    handle_w.write(str_write)
    if(freq_flag=='TRUE'):
	handle_w.close()
    all_codons = object_exon_file.index_codon.values()
    for key in all_codons:
	kmer_str = (''.join(key))
	object.kmer_increment_individual(kmer_str,1)

def reverse_complement(string):
    #This function returns the reverse complement of a string. For example CAT will return ATG
    dict_c = {'A':'T', 'C':'G', 'G':'C', 'T':'A','N':'N'}
    rev_c = []
    for char in string:
        rev_c.insert(0,dict_c[char])
    return(''.join(rev_c))

def sum_by_1(k):
    #This finds the sum of the numbers 1/i where i varies from 1 to k-1
    count = 0.0
    for i in range(1,k):
	count += (1.0/i)
    return(count)

def save_objects_transition_mat(objects,filename_save,object_fa,acceptable_chr = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]):
    #This function takes a dictioary of objects of class transition matrix and saves the mutation rate for each chromosome and combined across all autosomes to filename_save
    #object_fa is the object of class total_length_per_chr.It saves the length for each chromosome
    #We save non normalized mutation rates.
    #We also save the Segregating sites per chromosome and also the length of each chromosome
    #We ASSUME THAT WE HAVE SEEN ALL AUTOMOSOMES AND CHROMOSOME X in the objects DIRECTORY, as we only use this for the TRAINING stuff only
    handle = file(filename_save,"w")
#    acceptable_chr = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]
    handle.write("key_alpha"+'\t'+"key_beta"+'\t')
    for entry in acceptable_chr:
	handle.write(entry+'\t')
    handle.write("combined"+'\n')
    chr_seen = objects.keys()
    keys = objects[chr_seen[0]].mat.keys()
    for key in keys:
	count_key = 0.0
	count_kmerseq = 0.0
	handle.write(key[0]+'\t'+key[1]+'\t')
	for chrom in acceptable_chr:
	    if not(chrom == 'X'):
	        count_key += objects[chrom].mat[key]['total']
	        count_kmerseq += objects[chrom].kmerseq_count[key[0]]
	    rate = 0.0
	    if (objects[chrom].kmerseq_count[key[0]]>0):
		rate = objects[chrom].mat[key]['total']/(1.0*objects[chrom].kmerseq_count[key[0]])
	    handle.write(str(round(rate,5))+'\t')
	rate_a = 0.0
	if (count_kmerseq > 0):
	    rate_a = count_key/count_kmerseq
	handle.write(str(round(rate_a,5))+'\n')
    #Save the number of segregating sites per chromosome
    handle.write("seg_sites"+'\t'+'NA'+'\t')
    for chrom in acceptable_chr:
	count = 0 #Save the total number of segregating sites per chromosome here
	for key in keys:
	     count = count + objects[chrom].mat[key]['total']
	handle.write(str(count)+'\t')
    handle.write('NA'+'\n') # Combined column does not have entry
    handle.write("length_chr"+'\t'+'NA'+'\t')
    for chrom in acceptable_chr:
	handle.write(str(object_fa.chr_len[chrom])+'\t')
    handle.write('NA'+'\n')
    handle.close()

def save_objects_transition_mat_all(objects,filename_save,acceptable_chr = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]):
    #This function takes a dictioary of objects of class transition matrix and saves the kmer1-kmer2 count and kmer1 count for each chromosome
    handle = file(filename_save,"w")
#    acceptable_chr = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]
    handle.write("key_alpha"+'\t'+"key_beta"+'\t')
    for entry in acceptable_chr:
        handle.write(entry+'\t'+entry+'\t')
    handle.write("\n")
    chr_seen = objects.keys()
    keys = objects[chr_seen[0]].mat.keys()
    for key in keys:
        handle.write(key[0]+'\t'+key[1])
        for chrom in acceptable_chr:
	    handle.write('\t'+str(objects[chrom].mat[key]['total'])+'\t'+str(objects[chrom].kmerseq_count[key[0]]))
        handle.write('\n')

def find_entropy(learned_rate_file):
    #This functions finds the entropy of the learned rate file. It actually only works on the combined rate across the autosomes only.
    #TODO: We can actually find separate ENTROPY FOR THE SEX CHROMOSOME AND SEE IF IT IS BETTER TO HAVE BIGGER CONTEXTS FOR THAT
    handle = file(learned_rate_file,"r")
    dict_rate = {}
    content = handle.readlines()
    for entry in content[1:len(content)-2]:
        entry = entry.split('\t')
        dict_rate[(entry[0],entry[1])] = float(entry[25])
    rates = dict_rate.values()
    sum_r = sum(rates)
    entropy = 0.0
    for entry in rates:
	if(entry >0):
	    entropy = entropy + -1*(entry/sum_r)*math.log((entry/sum_r),2) 
    print "Entropy is " + str(round(entropy,4))	    
    num_keys = len(rates)
    worst_case_entropy = num_keys*( -1*(1.0/num_keys)*math.log((1.0/num_keys),2))
    print "No model entropy would be "+ str(round(worst_case_entropy,4))
    print "Entropy Reduction is "+str(round(worst_case_entropy-entropy,4))
class interval_file_analysis:
#This class will find all the variants present in each line of the interval, and save them to another file 
#This will also then save the fasta sequence.
#Then we will save another files (akin to train_non_coding) where by we will save all the kmer counts, and kmer transition_info.
    def __init__(self,filename):
        self.filename = filename #Save the filename on which the whole analysis has to be done
    def save_variants(self,pop):
        #This file will save variants which are delim separated for each line in the original file
        dummy_annotations_interval_saved_delimiter(self.filename,self.filename+"_"+pop+"_delim_annot","all_var_"+pop+"_chr_loc")
    def save_sequence(self):
        #This function saves the sequence corresponding to each line in the interval file
        files = file_handles("/project/voight_datasets/hg19/")
        chr_len_obj = total_length_per_chr(self.filename)
        chr_len_obj.save_fasta_sequence(files,0,0,True)
    def update_kmer_count(self,chr_seen,dir_obj_trans_mat1,dir_obj_trans_mat2):
        #This function updates the kmerseq count of trans_mat1 object for the chromosome chr_seen, with kmerseq count of dir_obj_trans_mat2 for chr_seen
        if not(chr_seen in dir_obj_trans_mat1):
            dir_obj_trans_mat1[chr_seen] = transition_mat(2,2)
        dir_obj_trans_mat1[chr_seen].kmerseq_count = dir_obj_trans_mat2[chr_seen].kmerseq_count
    def write_string(self,chr,start,end,kmers_n,keys_n,obj_trans_mat1):
        #This function creates a string that is supposed to be writen to the results file. We save the chr, start, end, all kmerseq counts and transition counts in the string, separated by comma
        string = chr+','+start+','+end+','
        for kmer in kmers_n:
            string += str(obj_trans_mat1.kmerseq_count[kmer])+','
        for key in keys_n:
            string += str(obj_trans_mat1.mat[key]['total'])+','
        string = string.rstrip(',') + '\n'
        return(string)

    def train_interval_file(self,pop):
        #This function will save the kmer count and transition info
	files = file_handles("project/voight_datasets/hg19/")
        handle_all = file(self.filename,"r")
        handle_fasta = file(self.filename+"_fasta_sequence","r")
        handle_pop_annot = file(self.filename+"_"+pop+"_delim_annot","r")

        handle_results_pop_7 = file(self.filename+"_before3_after3_results_"+pop,"w")
        handle_results_pop_5 = file(self.filename+"_before2_after2_results_"+pop,"w")
        handle_results_pop_3 = file(self.filename+"_before1_after1_results_"+pop,"w")
        handle_results_pop_1 = file(self.filename+"_before0_after0_results_"+pop,"w")

        handle_keys_7 = file(self.filename+"_before3_after3_keys","w")
        handle_keys_5 = file(self.filename+"_before2_after2_keys","w")
        handle_keys_3 = file(self.filename+"_before1_after1_keys","w")
        handle_keys_1 = file(self.filename+"_before0_after0_keys","w")
        content_all = handle_all.readlines()
        length_all = len(content_all)
        #Saving each fasta sequence and the variants in a separate file. Do not want the program to get confused to add this random number on to them.
        random_nos = ''.join(random.sample(string.lowercase,15)) #Now more secure. Earlier, we just used to save a number but now a random string of length 10
        kmers = None
        keys  = None
        for index in range(length_all):
            content_fasta = handle_fasta.readline()
            entry = content_all[index].rstrip('\n').split('\t')
            chr_seen = entry[0]
            start_p = entry[1]
            end_p = entry[2]
            #Save fasta sequence for that interval in a separate file
            handle = file("junk_delete.interval_"+str(random_nos)+"_fasta_sequence","w")
            handle.write(content_fasta)

            #Save african variants in a separate file
            handle = file('junk_delete_'+pop+'_annot_'+str(random_nos),'w')
            flag = "TRUE"
            while(flag=="TRUE"):
                entry = handle_pop_annot.readline()
                if not(entry[0] == '#'):
                    handle.write(entry)
                else:
                    handle.close()
                    flag = 'FALSE'
	    test_dir_7_pop = {}
            test_dir_5_pop = {}
            test_dir_3_pop = {}
            test_dir_1_pop = {}

            #Directory containing objects of class trans_mat being updated for 5 mer only. For 3 and 1 mer we will get from this.
            analyze_var_file("junk_delete_"+pop+"_annot_"+str(random_nos),test_dir_7_pop,files,3,3)

            #Updating the counts of 1,3 and 5 and 7 mer
            chr_len_obj = total_length_per_chr("junk_delete.interval_"+str(random_nos))
            chr_len_obj.update_kmer_count(test_dir_7_pop,3,3,strand=True)
            chr_len_obj = total_length_per_chr("junk_delete.interval_"+str(random_nos))
            chr_len_obj.update_kmer_count(test_dir_5_pop,2,2,strand=True)
            chr_len_obj = total_length_per_chr("junk_delete.interval_"+str(random_nos))
            chr_len_obj.update_kmer_count(test_dir_3_pop,1,1,strand=True)
            chr_len_obj = total_length_per_chr("junk_delete.interval_"+str(random_nos))
            chr_len_obj.update_kmer_count(test_dir_1_pop,0,0,strand=True)


            if(kmers==None):
		kmers_7 = test_dir_7_pop[chr_seen].kmerseq #The possible kmers for 7 mer. Only save once
                kmers_5 = test_dir_5_pop[chr_seen].kmerseq #The possible kmers for 5 mer. Only save once
                kmers_3 = test_dir_3_pop[chr_seen].kmerseq #The possible kmers for 3 mer. Only save once
                kmers_1 = test_dir_1_pop[chr_seen].kmerseq #The possible kmers for 1 mer. Only save once
            if(keys==None):
		keys_7 = test_dir_7_pop[chr_seen].mat.keys()
                keys_5 = test_dir_5_pop[chr_seen].mat.keys()
                keys_3 = test_dir_3_pop[chr_seen].mat.keys()
                keys_1 = test_dir_1_pop[chr_seen].mat.keys()

            test_dir_5_pop[chr_seen].mat = update_trans_count(test_dir_7_pop[chr_seen],2,2)
            test_dir_3_pop[chr_seen].mat = update_trans_count(test_dir_7_pop[chr_seen],1,1)
            test_dir_1_pop[chr_seen].mat = update_trans_count(test_dir_7_pop[chr_seen],0,0)

            str_pop_7 = self.write_string(chr_seen,start_p,end_p,kmers_7,keys_7,test_dir_7_pop[chr_seen])
            str_pop_5 = self.write_string(chr_seen,start_p,end_p,kmers_5,keys_5,test_dir_5_pop[chr_seen])
            str_pop_3 = self.write_string(chr_seen,start_p,end_p,kmers_3,keys_3,test_dir_3_pop[chr_seen])
            str_pop_1 = self.write_string(chr_seen,start_p,end_p,kmers_1,keys_1,test_dir_1_pop[chr_seen])

            handle_results_pop_7.write(str_pop_7)
            handle_results_pop_5.write(str_pop_5)
            handle_results_pop_3.write(str_pop_3)
            handle_results_pop_1.write(str_pop_1)

        str_7 = 'chr_name'+'\t'+'start_pos'+'\t'+'end_pos'+'\n'
        for kmer in kmers_7:
            str_7 += str(kmer)+'\t'
        str_7 = str_7.rstrip('\t')+'\n'
        for key in keys_7:
            str_7 += str(key)+'\t'
        str_7 = str_7.rstrip('\t') + '\n'
        handle_keys_7.write(str_7)

        str_5 = 'chr_name'+'\t'+'start_pos'+'\t'+'end_pos'+'\n'
        for kmer in kmers_5:
            str_5 += str(kmer)+'\t'
        str_5 = str_5.rstrip('\t')+'\n'
        for key in keys_5:
            str_5 += str(key)+'\t'
        str_5 = str_5.rstrip('\t') + '\n'
        handle_keys_5.write(str_5)

        str_3 = 'chr_name'+'\t'+'start_pos'+'\t'+'end_pos'+'\n'
        for kmer in kmers_3:
            str_3 += str(kmer)+'\t'
        str_3 = str_3.rstrip('\t') + '\n'
        for key in keys_3:
            str_3 += str(key)+'\t'
        str_3 = str_3.rstrip('\t') + '\n'
        handle_keys_3.write(str_3)

        str_1 = 'chr_name'+'\t'+'start_pos'+'\t'+'end_pos'+'\n'
        for kmer in kmers_1:
            str_1 += str(kmer)+'\t'
        str_1 = str_1.rstrip('\t') + '\n'
        for key in keys_1:
            str_1 += str(key)+'\t'
        str_1 = str_1.rstrip('\t') + '\n'
        handle_keys_1.write(str_1)
    def save_mutated_sequence(self,pop):
	#This function saves takes the fasta sequence file, and the delim annotated variant files, and then creates a file with extra column of the mutated sequence
        handle_fasta = file(self.filename+"_fasta_sequence","r")
        handle_pop_annot = file(self.filename+"_"+pop+"_delim_annot","r")
	content_fasta = handle_fasta.readline()
	handle_mutated_fasta_pop = file(self.filename+"_fasta_sequence_mutated_"+pop,"w")
	counter = 0
	while(content_fasta):

	    var_pop = {}
            flag = "TRUE"
            while(flag=="TRUE"):
                entry = handle_pop_annot.readline()
		entry = entry.rstrip('\n').split('\t')
                if not(entry[0] == '#'):
		    var_pop[(entry[0],entry[1])] = (entry[2],entry[3])
            	else:
                    flag = 'FALSE'

	    entry = content_fasta.rstrip('\n').split('\t')
	    if len(var_pop) == 0:
                mut_sequence = list(copy.deepcopy(entry[4]))
	        handle_mutated_fasta_pop.write(entry[0]+'\t'+entry[1]+'\t'+entry[2]+'\t'+entry[3]+'\t'+entry[4]+'\t'+"".join(mut_sequence)+'\n')
	    else:
		strand = entry[3]
	        if strand == '+':
		    mut_sequence = list(copy.deepcopy(entry[4]))
	            for key in var_pop:
		        position_mutate = int(key[1]) - int(entry[1])
		        if (mut_sequence[position_mutate] == var_pop[key][0]):
		            mut_sequence[position_mutate] = var_pop[key][1]
		            counter += 1
		        elif mut_sequence[position_mutate] == reverse_complement(var_pop[key][0]):  #Had variant corresponding to reverse 
                            mut_sequence[position_mutate] = reverse_complement(var_pop[key][1])
                            counter += 1
		        else:
			    print "Something intertesting going on. Maybe the sequence is wrong"
			    pdb.set_trace()
		if strand == '-': #We are working with reverse complement so everything is broken. Take reverse complement again and then save the reverse complement of the final mutated thing
		    mut_sequence = list(reverse_complement(copy.deepcopy(entry[4])))
		    for key in var_pop:
                        position_mutate = int(key[1]) - int(entry[1])
                        if (mut_sequence[position_mutate] == var_pop[key][0]):
                            mut_sequence[position_mutate] = var_pop[key][1]
                            counter += 1
                        elif mut_sequence[position_mutate] == reverse_complement(var_pop[key][0]):  #Had variant corresponding to reverse 
                            mut_sequence[position_mutate] = reverse_complement(var_pop[key][1])
                            counter += 1
                        else:
                            print "Something intertesting going on. Maybe the sequence is wrong"
                            pdb.set_trace()
                    mut_sequence = list(reverse_complement(mut_sequence))
	        handle_mutated_fasta_pop.write(entry[0]+'\t'+entry[1]+'\t'+entry[2]+'\t'+entry[3]+'\t'+entry[4]+'\t'+"".join(mut_sequence)+'\n')
            content_fasta = handle_fasta.readline()	        
        print counter

    def generate_mapping(self,filename_keys):
        #This function takes a filename which contains information on all the keys etc, and creates a mapping out of it.
        handle = file(filename_keys,"r")
        content_keys = handle.readlines()
        index = 3   #So that we know which kmer (or kmer transition) has which index. Also the first 3 are just chr name, start and end so they do not count.
        mapping = {} #Dictionary that contains mapping from a name to a number(here column number) in the results file
        kmers = []  #Basically all possible kmers
        entry = content_keys[1].rstrip('\n').split('\t')
        for elem in entry:
            kmers.append(elem)
            mapping[index] = elem
            index += 1
        keys = [] #Basically kmer1->kmer2
        entry = content_keys[2].rstrip('\n').split('\t')
        for elem in entry:
            keys.append(eval(elem))
            mapping[index] = eval(elem)
            index += 1
        return([mapping,kmers,keys])

    def get_combined_results(self,handle_file,mapping):
        #This function takes the handle for the results file and creates a dictionary which contains the combined result of the entire results file
        dict_res = {}
        start_from = 3 #As the first 3 positions in the results file are just chromosome name, start and end
        content_results = handle_file.readline()
        while(content_results):
            entry = content_results.rstrip('\n').split(',')
            for index in range(len(entry)):
                if(index>=start_from):
                    if not(mapping[index] in dict_res):
                        dict_res[mapping[index]] = 0
                    dict_res[mapping[index]] += int(entry[index])
            content_results = handle_file.readline()
        return(dict_res)

    def write_combined_results(self,dict_res,filename):
        #This function takes the dict_res dictionary file and then saves it results in the file named filename
        handle = file(filename,"w")
        for entry in dict_res:
            if len(entry)==2:
                handle.write(entry[0]+'\t'+entry[1]+'\t'+str(dict_res[entry])+'\n')
            else:
                handle.write(entry+'\t'+str(dict_res[entry])+'\n')
        handle.close()

    def analyze_results(self,pop):
        #This function will save the rate files and the count,transition files. We will be doing this together, and not for each chromosome separately

        [mapping,kmers,keys] = self.generate_mapping(self.filename+"_before0_after0_keys")
        handle_results_pop = file(self.filename+"_before0_after0_results_"+pop,"r")
        dict_pop = self.get_combined_results(handle_results_pop,mapping)
        self.write_combined_results(dict_pop,self.filename+"_"+pop+"_before0_after0")

	#3mer
        [mapping,kmers,keys] = self.generate_mapping(self.filename+"_before1_after1_keys")
        handle_results_pop = file(self.filename+"_before1_after1_results_"+pop,"r")
        dict_pop = self.get_combined_results(handle_results_pop,mapping)
        self.write_combined_results(dict_pop,self.filename+"_"+pop+"_before1_after1")
	#5mer
        [mapping,kmers,keys] = self.generate_mapping(self.filename+"_before2_after2_keys")
        handle_results_pop = file(self.filename+"_before2_after2_results_"+pop,"r")
        dict_pop = self.get_combined_results(handle_results_pop,mapping)
        self.write_combined_results(dict_pop,self.filename+"_"+pop+"_before2_after2")
        #7mer
        [mapping,kmers,keys] = self.generate_mapping(self.filename+"_before3_after3_keys")
        handle_results_pop = file(self.filename+"_before3_after3_results_"+pop,"r")
        dict_pop = self.get_combined_results(handle_results_pop,mapping)
        self.write_combined_results(dict_pop,self.filename+"_"+pop+"_before3_after3")
	
    def find_waterson_theta_tf(self,handle,norm):
    #This function finds the waterson theta for each transciption factor and total mutable sites and total length
        content = handle.readline()
        total_len = 0
        total_mut = 0
        while(content):
            entry = content.rstrip('\n').split('\t')
            if len(entry)==3:
                total_mut += int(entry[2])
            else:
                total_len += int(entry[1])
            content = handle.readline()
	return total_mut*1.0/(norm*total_len),total_mut,total_len

    def find_waterson_theta_tf_pos(self,filename,norm):
        #This function finds the waterson theta for each position in the transcription factor
	handle = file(filename,"r")
        content = handle.readline()
	occur = 1   #total binding sites
        while(content):
	    content = content.rstrip('\n').split('\t')
            orig = content[4]
	    mut  = content[5]
	    length_site = len(orig)
	    if not('dict_mucount' in locals()):	
    	        dict_mucount = {}
	        for index in range(length_site):
		    dict_mucount[index] = 0
	    change = [i for i,(a1,a2)  in enumerate(zip(orig,mut)) if a1!=a2]
	    for index in change:
	        dict_mucount[index] += 1
            content = handle.readline()
	    occur +=  1
        #each position occurs = total binding sites time
        waterson_pos = [0]*length_site
	mut_pos = [0]*length_site
	len_pos = [0]*length_site
        for key in dict_mucount:
	    waterson_pos[key] = str(round(dict_mucount[key]*1.0/(norm*occur),5))
	    mut_pos[key]  = str(dict_mucount[key]*1.0)
	    len_pos[key] = str(occur)
        return ','.join(waterson_pos),','.join(mut_pos),','.join(len_pos)
    

	
