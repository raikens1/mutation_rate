import os
import string
import random
import numpy
import scipy
from scipy.stats import *
import pdb
import copy
from modules import *
from numpy import array, log, exp, mean
from scipy.special import gammaln
import operator


def dummy_annotations_interval_saved_delimiter(interval_file,annot_file,all_var_file,saved_index="F"):
    #This function finds the variants in the regions from the interval_file and saves there dummy annotations in the variant file. 
    #After finding the variants and their annotation for a particular line in the interval file, it also adds # at the next line.
    #The all_var_file is the mapping of all the variants for a particular population.
    #Format of all_var file is Chr_name Pos Ref Alt Anc DbSNP_id Freq
    #Format of annot file will be Chr_name Pos Ref Alt Dummy Freq Anc
    #This function would be useful if we a big file of intervals and instead of saving all variants for it in one file, we save variants for each interval in one file but separated by a delimiter
    if(saved_index=="F"):
        handle = file(all_var_file,"r")
        dict_all_var = {}
        content = handle.readline()
        while(content):
            content = content.rstrip('\n').split('\t')
            if content[0] in dict_all_var:
                dict_all_var[content[0]][int(content[1])] = [content[2],content[3],'Dummy',content[6],content[4]]
            else:
                dict_all_var[content[0]] = {}
                dict_all_var[content[0]][int(content[1])] = [content[2],content[3],'Dummy',content[6],content[4]]
            content = handle.readline()
	handle.close()
    else:
	dict_all_var = saved_index
    handle = file(interval_file,"r")
    handle_w = file(annot_file,"w")
    content = handle.readlines()
    for entry in content:
        entry = entry.rstrip('\t').split('\t')
        chr_name = entry[0]
        start = int(entry[1])
        end = int(entry[2])
        if chr_name in dict_all_var:
            for pos in range(start,(end+1)):
                if(pos in dict_all_var[chr_name]):
                    handle_w.write(chr_name+'\t'+str(pos)+'\t'+'\t'.join(dict_all_var[chr_name][pos])+'\n')
	handle_w.write('#'+'\n')
    handle.close()
    handle_w.close()


def update_trans_count(object,before,after):
    #This function updates the transition matrix count from the object counts. Object is the object of the class transition_mat. before and after are sizes for new reduced stuff
    dict_return = {}
    for key in object.mat.keys():
	kmer_size = len(key[0])		
	alpha = key[0][(kmer_size/2)-before:(kmer_size/2)+after+1]
	beta = key[1][(kmer_size/2)-before:(kmer_size/2)+after+1]
	new_key = (alpha,beta)
	if not(new_key in dict_return):
	    dict_return[new_key] = {}
	    dict_return[new_key]['total'] = 0 
	dict_return[new_key]['total'] += object.mat[key]['total']
    return(dict_return)

def get_dict_info_file(filename,acceptable_chr):
#This function returns the dictionary containing the count of transitions and occurances from the file _all
    handle = file(filename,"r")
    dict_rate = {} #Contains the rate for all pairs for all chromosomes
    for chrom in acceptable_chr:
        dict_rate[chrom] = {}
    content = handle.readlines()
    handle.close()
    for entry in content[1:len(content)]:
        entry = entry.split('\t')
        for index in range(2,len(acceptable_chr)*2+2): #From 2 till 24 we have rates for chromosomes 1 till X
            if(index%2 ==0):
                dict_rate[acceptable_chr[(index/2)-1]][(entry[0],entry[1])] = float(entry[index])
            else:
                dict_rate[acceptable_chr[(index/2)-1]][entry[0]] = float(entry[index])
        count_trans = 0
        count_key = 0
        for autosome in ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"]:
            count_trans += dict_rate[autosome][(entry[0],entry[1])]
            count_key += dict_rate[autosome][entry[0]]
        if not 'combined' in dict_rate:
            dict_rate['combined'] = {}
        dict_rate['combined'][(entry[0],entry[1])] = count_trans
        dict_rate['combined'][entry[0]] = count_key
    return(dict_rate)


def get_dict_rate_file(filename,bottom=True):
#This function returns the dictionary containing the transition rates from the filename
    acceptable_chr = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","combined"]
    handle = file(filename,"r")
    dict_rate = {} #Contains the rate for all pairs for all chromosomes
    for chrom in acceptable_chr:
        dict_rate[chrom] = {}
    content = handle.readlines()
    if bottom==True:
	lines_read = len(content)-2
    else:
	lines_read = len(content)
    for entry in content[1:lines_read]:
        entry = entry.split('\t')
	total_res_file = len(entry)
        for index in range(2,total_res_file): #From 2 till 25 we have rates for chromosomes 1 till combined
            dict_rate[acceptable_chr[index-2]][(entry[0],entry[1])] = float(entry[index])
    handle.close()
    return(dict_rate)

def get_bayesian_dict_rate_file(filename):
    #This will return the IGR rate with uniform prior
    acceptable_chr = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]
    dict_count = get_dict_info_file(filename,acceptable_chr)
    acceptable_chr.append("combined")
    dict_rate = {}
    for chrom in acceptable_chr:
        if not chrom in dict_rate:
	    dict_rate[chrom] = {}
        for entry in dict_count['combined']:
	    if len(entry)==2:
		dict_rate[chrom][entry] = (dict_count[chrom][entry]+1)/(dict_count[chrom][entry[0]]+2) 
    return dict_rate

 
def get_null_rate_file(filename):
#This function returns the rate under the Null Poisson model
    acceptable_chr = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]
    handle = file(filename,"r")
    dict_rate = {} #Contains the rate for all pairs for all chromosomes
    for chrom in acceptable_chr:
        dict_rate[chrom] = {}
    content = handle.readlines()
    sites = content[len(content)-2].rstrip('\n').split('\t')
    len_c = content[len(content)-1].rstrip('\n').split('\t')
    for entry in content[1:len(content)-2]:
        entry = entry.split('\t')
        for index in range(2,25): #From 2 till 24 we have rates for chromosomes 1 till X
            dict_rate[acceptable_chr[index-2]][(entry[0],entry[1])] = round(float(sites[index])/(float(len_c[index])*3),5)
    handle.close()
    return(dict_rate)


def get_counts_kmer_transitions(data_file,keys_file,folded=False):
#Keys_file is the file which contains information on which column name stands for what
#if folded true then reverse compliment results are returned
    handle_keys = file(keys_file,"r")
    content_keys = handle_keys.readlines()
    index = 0
    mapping = {}
    entry = content_keys[0].rstrip('\n').split('\t')
    for elem in entry:
        mapping[index] = elem
        index += 1
    entry = content_keys[1].rstrip('\n').split('\t')
    for elem in entry:
        mapping[index] = elem
        index += 1
    entry = content_keys[2].rstrip('\n').split('\t')
    for elem in entry:
        mapping[index] = eval(elem)
        index += 1
    handle_keys.close()
    handle_results = file(data_file,"r")
    content_results = handle_results.readlines()
    dict_data = {}
    for entry in content_results:
	entry = entry.rstrip('\n').split(',')
	for index in range(len(entry)):
          try:
            if(index>=3):	
                if not(mapping[index] in dict_data):
                    dict_data[mapping[index]] = 0
                dict_data[mapping[index]] += int(float(entry[index]))
	  except:
	    pdb.set_trace()
    handle_results.close()
    if(folded=='T'):
	dict_data_folded = {}
        folded_trans = set()
	for key in dict_data:
	    if(len(key)==2):
		(alpha,beta) = key
		before_after = len(alpha)/2
	        if (   alpha[before_after]=='A' or alpha[before_after]=='C'   ):
            	    folded_trans.add( (alpha,beta) )
	for (alpha,beta) in folded_trans:
            dict_data_folded[alpha] = dict_data[alpha]+dict_data[reverse_complement(alpha)]
            dict_data_folded[(alpha,beta)] = dict_data[(alpha,beta)]+dict_data[(reverse_complement(alpha),reverse_complement(beta))]
	return dict_data_folded
    else:
	return(dict_data)

def get_codon_info(data_file,mapping):
#This function returns the dictionary containing information about the count of all codons and codon to codon transitions that are present in the file
#mapping if the mapping of the entries in the file to the keys
    handle_results = file(data_file,"r")
    content_results = handle_results.readlines()
    dict_data = {}
    for entry in content_results:
        entry = entry.rstrip('\n').split(',')
        for index in range(len(entry)):
	    chrom = entry[0]
            if(index>=3):
		if(len(mapping[index]) ==2):
                    if not((mapping[index][0],mapping[index][1]) in dict_data):
                       dict_data[(mapping[index][0],mapping[index][1])] = 0
		    dict_data[(mapping[index][0],mapping[index][1])] += int(entry[index])
	        else:
		    if not(mapping[index] in dict_data):
			dict_data[mapping[index]] = 0
		    dict_data[mapping[index]] += int(entry[index])
    handle_results.close()
    return dict_data

def update_frequency_table(filename,dict_freq,amino_mapping='FALSE'):     
#This function takes the frequency table file and updates the dict_freq data structure. Amino mapping is the mapping from nucleotide to Amino acid name
    handle = file(filename,"r")
    content = handle.readlines()
    for entry in content:
	entry = entry.rstrip('\n').split('\t')
	if(amino_mapping=='FALSE'):
	    key = (entry[0],entry[1])
	else:
	    key = (amino_mapping[entry[0]],amino_mapping[entry[1]])
	if not(key in dict_freq):
	    dict_freq[key] = []
	dict_freq[key].append(float(entry[2]))
    handle.close()	

def keys_dict(transitions):
    #This function returns a dictionary where for each kmer it returns all the 4 possibilities.
    #Transitions are basically the set of all (alpha,beta) transitions
    keys_dict_return = {}
    nucleotide = ['A','C','G','T']
    for entry in transitions:
        if not(entry[0] in keys_dict_return):
            keys_dict_return[entry[0]] = []
    for key in keys_dict_return.keys():
        pos = len(key)/2
        for nuc in nucleotide:
            key_new = list(key)
            key_new[pos] = nuc
            keys_dict_return[key].append(''.join(key_new))
    return(keys_dict_return)

def rate_new_model(rate_dict_big,rate_dict_small,keys_more_detail):
    #This function makes a new rate model by only considering details of the transitions that fall in keys_more_detail. Otherwise it substitutes them by keeping transitions from a smaller kmer rate model
    #rate_dict_big is the rate for the bigger kmer, rate_dict_small is the rate of smaller kmers, keys_more_detail are the transitions for which we have to consider the larger context rates
    dict_return = {}
    chr_seen = rate_dict_big.keys()[0]
    size_small = len(rate_dict_small[chr_seen].keys()[0][0])
    size_large = len(rate_dict_big[chr_seen].keys()[0][0])
    start = (size_large/2)-(size_small/2)
    end = (size_large/2)+(size_small/2)+1
    chr_seen = rate_dict_small.keys()
    for chrom in chr_seen:
        if not(chrom in dict_return):
            dict_return[chrom] = {}
        for key in rate_dict_big[chrom].keys():
            if not( (key[0][start:end],key[1][start:end]) in keys_more_detail):
                dict_return[chrom][key] = rate_dict_small[chrom][(key[0][start:end],key[1][start:end])]
            else:
                dict_return[chrom][key] = rate_dict_big[chrom][key]
    return(dict_return)



def get_total_seg(keys,result):
#keys are the possible transitions. Result includes the count of the key in keys    
    count = 0
    for key in keys:
        count += result[key]
    return(count)

def get_dict_rate_from_region(chrom,res_1,res_3,res_5,keys_1,keys_3,keys_5):
    #This functions returns the MLE of all the rate matrices calculated from the best estimate from a region of the genome
    rate_0 = {}; rate_1 = {}; rate_3 = {}; rate_5 = {}
    rate_0[chrom] = {}; rate_1[chrom] = {}; rate_3[chrom] = {}; rate_5[chrom] = {};
    total_seg = 0
    total_len = res_1['A'] + res_1['C'] + res_1['G'] + res_1['T']
    for key in keys_1:
	total_seg += res_1[key]
    for key in keys_1:
	rate_0[chrom][key] = round(total_seg*1.0/(3*total_len),5)
	if(res_1[key[0]] > 0):
	    rate_1[chrom][key] = round(res_1[key]*1.0/res_1[key[0]],5)
	else:
	    rate_1[chrom][key] = 0.0
    for key in keys_3:
        if(res_3[key[0]] > 0):
            rate_3[chrom][key] = round(res_3[key]*1.0/res_3[key[0]],5)
        else:
            rate_3[chrom][key] = 0.0
    for key in keys_5:
        if(res_5[key[0]] > 0):
            rate_5[chrom][key] = round(res_5[key]*1.0/res_5[key[0]],5)
        else:
            rate_5[chrom][key] = 0.0
    return([rate_0,rate_1,rate_3,rate_5])

class test_model:
    #This class does a lot of thing. 
    # 1)It calculates the log likelihood of the observed data
    # 2)Samples from the data and returns the median value of each mutation type
    # 3)Samples and also returns the 95% confidence interval of each mutation type
    # 4)Also returns the categories for which observed falls into the simulated values
    # 5)Also returns the bias^2 + variance accuracy
    def __init__(self,chrom,dict_rate,results,keys_kmer):
        #chrom is the current chromosome, dict_rate is the dictionary indexed by each chromosome and has the rate for each transition. 
        #results has the counts for each kmer of particular type and also transition. Keys_kmer basically has for each kmer, all the possible kmers it may go to
        self.chrom = chrom
        self.dict_rate = dict_rate
        self.results = results
        self.keys_kmer = keys_kmer
        self.nsamp = 10000
    def get_prob(self,alpha):
        #This function returns the probability of transition of a kmer to all possible kmers. (input X, output [a,b,c,d] where a is prob of X mutating to A, b is prob of X mutating to C, respectively.
        probs = []                   #Finding probabilities for all 4 possible trasitions
        for beta in self.keys_kmer[alpha]:
            if not(alpha==beta):
                if(self.dict_rate[self.chrom][(alpha,beta)] == 0):
                    probs.append(self.dict_rate[self.chrom][(alpha,beta)]+.0000001)
                else:
                    probs.append(self.dict_rate[self.chrom][(alpha,beta)])
            else:
                prob = 0.0
                for gamma in self.keys_kmer[alpha]:
                    if not(alpha==gamma):
                        prob += self.dict_rate[self.chrom][(alpha,gamma)]
                if(prob >=1):
                    prob = 0.9999999
                probs.append(1-prob)
        prob_new = [x/sum(probs) for x in probs]
        return prob_new
    def log_factorial(self,x):
        """Returns the logarithm of x!
        Also accepts lists and NumPy arrays in place of x."""
        return gammaln(array(x)+1)

    def multinomial(self,n, xs, ps):
        #Returns log of pmf of mutlinomial distribution 
        xs, ps = array(xs), array(ps)
        result = self.log_factorial(n) - sum(self.log_factorial(xs)) + sum(xs * log(ps))
        return result

    def get_segsites(self,aa_cat=False,offset=False):
	#This function returns the total mutations present
	#If aa_cat is not False then it finds the counts in each individual category using offset
	if aa_cat == False:
	    total_seg = 0
            for alpha in self.keys_kmer.keys():
                for beta in self.keys_kmer[alpha]:
	    	    if not(alpha==beta):
	                total_seg += self.results[(alpha,beta)]
	    return total_seg
	else:
            amino_mapping  = returns_amino_mapping()
	    total_seg = {}
   	    for code in [1,2,3,4]:
	        total_seg[code] = 0
            for alpha in self.keys_kmer.keys():
                for beta in self.keys_kmer[alpha]:
                     if not(alpha==beta):
                         length = len(alpha)
                         extract = [(length/2)+1-offset, (length/2)+2-offset, (length/2)+3-offset]
                         aa1 = amino_mapping[''.join([ alpha[extract[0]], alpha[extract[1]], alpha[extract[2]]])]
                         aa2 = amino_mapping[''.join([ beta[extract[0]], beta[extract[1]], beta[extract[2]]])]
                         if(aa1==aa2):
                             code = 1
                         elif(aa2 == "Stop"):
                             code = 3
                         elif(aa1 == "Stop"):
                             code = 4
                         elif not(aa1==aa2):
                             code = 2
			 total_seg[code] += self.results[(alpha,beta)]
	    return total_seg
    def get_aa_segsites(self,aminoacid,offset):
        #This function returns the total mutations present for a AA
        #It finds the counts in each individual category using offset
        amino_mapping  = returns_amino_mapping()
        total_seg = 0
        for alpha in self.keys_kmer.keys():
            for beta in self.keys_kmer[alpha]:
                 if not(alpha==beta):
                     length = len(alpha)
                     extract = [(length/2)+1-offset, (length/2)+2-offset, (length/2)+3-offset]
                     aa1 = amino_mapping[''.join([ alpha[extract[0]], alpha[extract[1]], alpha[extract[2]]])]
                     aa2 = amino_mapping[''.join([ beta[extract[0]], beta[extract[1]], beta[extract[2]]])]
		     if len(aminoacid)==1:
                         if not(aa1==aa2) and (aa1==aminoacid):
                             total_seg += self.results[(alpha,beta)]
		     else:
                         if not(aa1==aa2) and (aa1 in aminoacid) and not(aa2 in aminoacid):
                             total_seg += self.results[(alpha,beta)]		
        return total_seg

    def get_ll(self):
        #This function returns the log likelihood of the observed data based on the model
        ll = 0.0
        for alpha in self.keys_kmer.keys():
            prob_new = self.get_prob(alpha)
            counts_seen = []

            for beta in self.keys_kmer[alpha]:
		if not(alpha==beta):
                    counts_seen.append(self.results[(alpha,beta)])
		else:
		    total_trans = 0
		    for gamma in self.keys_kmer[alpha]:
			if not(alpha==gamma):
			    total_trans += self.results[(alpha,gamma)]
		    counts_seen.append(self.results[alpha]-total_trans)
            if(self.results[alpha]>0):
                ll += self.multinomial(self.results[alpha],counts_seen,prob_new)
        return(ll)
    def multisamp(self,alpha):
        #This function takes the kmer and then returns the sampled stuff. If the original count is zero then it returns a zero sample stuff
        prob_new = self.get_prob(alpha)
        if(self.results[alpha]>0):
            return(numpy.random.multinomial(self.results[alpha],prob_new,size=self.nsamp))
        else:
            return(numpy.zeros((self.nsamp,len(self.keys_kmer[alpha]) ) ) )
    def get_accuracy(self,return_stuff=False):
        #This function calculated the bias^2 + variance and reports the result
        accuracy  = 0.0
	bias_sq = 0.0
	variance = 0.0
	counter = 0
        for alpha in self.keys_kmer.keys():
            possible_mut = self.multisamp(alpha)
            index = 0
            for beta in self.keys_kmer[alpha]:
                if not(alpha==beta):
		    counter += 1
                    exp = self.results[(alpha,beta)]		
                    obs = numpy.mean(possible_mut[:,index])
                    bias_sq  += math.pow(exp-obs,2)
                    variance += numpy.var(possible_mut[:,index])
                index += 1 # So that we can access the next sampled stuff
	accuracy = bias_sq + variance 
	accuracy = math.sqrt(accuracy)
	if return_stuff == False:
	    return(accuracy)
	else:
	    return[bias_sq,variance]
    def return_simulated_variants(self):
	#This function returns a list of variants that are simulated by the model. The list has repeated entries denoting the fact that some kmers have more than one variant present in them
	#Used in the simulation 
	list_variants = []
	random_sim = random.randrange(self.nsamp) #which simulation to return
        for alpha in self.keys_kmer.keys():
            possible_mut = self.multisamp(alpha)
            index = 0
            for beta in self.keys_kmer[alpha]:
                if not(alpha==beta):
		    for i in range(possible_mut[random_sim,index]>0):
			list_variants.append((alpha,beta))
                index += 1
        return(list_variants)
    def get_total_variant(self):
        #This function returns the total variants as predicted by model
        total_mut = numpy.zeros(self.nsamp) #We save the total variants here
        for alpha in self.keys_kmer.keys():
            possible_mut = self.multisamp(alpha)
            index = 0
            for beta in self.keys_kmer[alpha]:
                if not(alpha==beta):
                    total_mut += possible_mut[:,index]
		index += 1
        return(total_mut[random.randrange(self.nsamp)])
    def get_95conf(self):
        #This function returns the 95% of the total variants 
        total_mut = numpy.zeros(self.nsamp) #We save the total variants here
        for alpha in self.keys_kmer.keys():
            possible_mut = self.multisamp(alpha)
            index = 0
            for beta in self.keys_kmer[alpha]:
                if not(alpha==beta):
                     total_mut += possible_mut[:,index]
		index += 1
        total_mut.sort()
        return( ( int(total_mut[int(self.nsamp*0.025)]),int(total_mut[int(self.nsamp*0.975)]) ) )
    def get_zscore_total_mut(self):
       #This function will find zscore for total_mutations and see how it differs from the expected
       #Zscore less than 0 means that we find less than expected.
        total_seg = self.get_segsites()
        total_mut = numpy.zeros(self.nsamp) #We save the total variants here
        for alpha in self.keys_kmer.keys():
            possible_mut = self.multisamp(alpha)
            index = 0
            for beta in self.keys_kmer[alpha]:
                if not(alpha==beta):
                     total_mut += possible_mut[:,index]
                index += 1
        zscore_total_mut = (-numpy.mean(total_mut) + total_seg)/numpy.std(total_mut)
        return zscore_total_mut
    def get_zscore_tstv(self):
    #This function will find zscore for ts/tv ratio and see how it differs from expected
        ts_obs = 0.0
        tv_obs = 0.0
        for alpha in self.keys_kmer.keys():
	    kmer = len(alpha)
            for beta in self.keys_kmer[alpha]:
                if not(alpha==beta):
		    if alpha[kmer/2] == 'C' and beta[kmer/2] == 'T':
			ts_obs += self.results[(alpha,beta)]
		    elif alpha[kmer/2] == 'T' and beta[kmer/2] == 'C':
                        ts_obs += self.results[(alpha,beta)]
                    elif alpha[kmer/2] == 'A' and beta[kmer/2] == 'G':
                        ts_obs += self.results[(alpha,beta)]
                    elif alpha[kmer/2] == 'G' and beta[kmer/2] == 'A':
                        ts_obs += self.results[(alpha,beta)]
		    else:
			tv_obs += self.results[(alpha,beta)]
        ts_exp = numpy.zeros(self.nsamp)
        tv_exp = numpy.zeros(self.nsamp)
        for alpha in self.keys_kmer.keys():
            kmer = len(alpha)
            possible_mut = self.multisamp(alpha)
            index = 0
            for beta in self.keys_kmer[alpha]:
                if not(alpha==beta):
                    if alpha[kmer/2] == 'C' and beta[kmer/2] == 'T':
		        ts_exp += possible_mut[:,index]
		    elif alpha[kmer/2] == 'T' and beta[kmer/2] == 'C':
                        ts_exp += possible_mut[:,index]
                    elif alpha[kmer/2] == 'A' and beta[kmer/2] == 'G':
                        ts_exp += possible_mut[:,index]
                    elif alpha[kmer/2] == 'G' and beta[kmer/2] == 'A':
                        ts_exp += possible_mut[:,index]
		    else:
			tv_exp += possible_mut[:,index]
		index += 1
	if tv_obs == 0:
	    return 'NA'
	else:
	    zscore_tstv = ( (ts_obs/tv_obs) - numpy.mean(ts_exp/tv_exp))/numpy.std(ts_exp/tv_exp)
	    return zscore_tstv

    def get_conf_score(self,offset,sample,total_seg,return_stuff=False):
        #This function is used specifically in testing of coding transcripts. We find using the overall number of variants, the rank of observed variants in the expected distribution and also the 95%CI.
	#Since we use this for offset1,2,3 so we need to make sure that we sum over the previous samplings (aka for all offsets) 
	#sample is added to the total_mut and is returned, if return_stuff=True then we return the CI and rank
	#we also do individually for each offset, as now we can get results for each indidividual aa too
	amino_mapping  = returns_amino_mapping()
	if return_stuff == False:
	    total_mut = {}
	    for code in [1,2,3,4]:
	        total_mut[code] = numpy.zeros(self.nsamp) #We save the total variants here
	    for alpha in self.keys_kmer.keys():
                possible_mut = self.multisamp(alpha)
	        index = 0
        	for beta in self.keys_kmer[alpha]:
                    if not(alpha==beta):
	                length = len(alpha)
                	extract = [(length/2)+1-offset, (length/2)+2-offset, (length/2)+3-offset]
                	aa1 = amino_mapping[''.join([ alpha[extract[0]], alpha[extract[1]], alpha[extract[2]]])]
                	aa2 = amino_mapping[''.join([ beta[extract[0]], beta[extract[1]], beta[extract[2]]])]
			if(aa1==aa2):
		            code = 1
		        elif(aa2 == "Stop"):
            		    code = 3
		        elif(aa1 == "Stop"):
		            code = 4
		        elif not(aa1==aa2):
		            code = 2
                        total_mut[code] += possible_mut[:,index]
	            index += 1
	    if not sample == None:
		for code in [1,2,3,4]:
	            total_mut[code] += sample[code]
	    return total_mut
	else:
            #Find overall rank and individual category ranks				
	    #overall
	    overall = numpy.zeros(self.nsamp)
	    overall_seg = 0
            for code in [1,2,3,4]:
	        overall += sample[code]
		overall_seg += total_seg[code]
	    overall.sort()
	    CI0 = int(overall[int(self.nsamp*0.025)])
	    CI1 = int(overall[int(self.nsamp*0.975)]) 

            #Find rank all code [2,3,4] combined
            total_nonsyn = total_seg[2]+total_seg[3]+total_seg[4]
            total_sample = sample[2]+sample[3]+sample[4]
            z_nonsyn = (numpy.mean(total_sample) - total_nonsyn)/numpy.std(total_sample)

	    z_overall = (numpy.mean(overall) - overall_seg)/numpy.std(overall)
	    z_code = {}
	    for code in [1,2,3,4]:
		sample[code].sort()
		z_code[code] = (numpy.mean(sample[code]) - total_seg[code])/numpy.std(sample[code])
	    return [CI0,CI1,z_overall,z_code[1],z_code[2],z_code[3],z_code[4],z_nonsyn]

    def get_conf_zscore_denovo(self,total_seg,file_exons):
	#We use this to find the sample from distribution and find zscores for different categories of change.
        #We use the same to find overall confidence interval too.
	#Total_seg contains the total variants seen till now. #File_exons is the object of class exon_file_properties
        amino_mapping  = returns_amino_mapping()
	total_mut = {}
	for code in ['overall','syn','mis','nonsens','splice']:
	    total_mut[code] = numpy.zeros(self.nsamp) #We save simulated variants here
        for pos in file_exons.pos_codon:
            (codon_num,codon_pos) = file_exons.pos_codon[pos]
            origkmer = copy.deepcopy(file_exons.pos_sequence[pos])
            if origkmer in self.keys_kmer:
	        prob_new = self.get_prob(origkmer)
	        possible_mut = numpy.random.multinomial(1,prob_new,size=self.nsamp)
		index = 0
                for change in self.keys_kmer[origkmer]:
                    if not change ==  origkmer:
			total_mut['overall'] += possible_mut[:,index]
                        if codon_num == 'splice_donor' or codon_num == 'splice_acceptor':
                            total_mut['splice'] += possible_mut[:,index]
                        elif codon_num in file_exons.index_codon:
                            orig_codon = file_exons.index_codon[codon_num]
                            alpha = amino_mapping[''.join(orig_codon)]
                            mut_codon = copy.deepcopy(orig_codon)
                            mut_codon[codon_pos-1] = change[len(origkmer)/2]
                            beta = amino_mapping[''.join(mut_codon)]
                            if(alpha==beta):
                                total_mut['syn'] += possible_mut[:,index]
                            elif (beta == "Stop"):
	                        total_mut['nonsens'] += possible_mut[:,index]
                            elif not(alpha==beta):
                                total_mut['mis'] += possible_mut[:,index]
		    index += 1

	total_mut['overall'].sort()
        CI0 = total_mut['overall'][int(self.nsamp*0.025)]
        CI1 = total_mut['overall'][int(self.nsamp*0.975)]
	z_overall = (numpy.mean(total_mut['overall']) - total_seg['overall'])/numpy.std(total_mut['overall'])

        total_nonsyn = total_seg['mis']+total_seg['nonsens']+total_seg['splice']
        total_nonsynsample = total_mut['mis']+total_mut['nonsens']+total_mut['splice']
        z_nonsyn = (numpy.mean(total_nonsynsample) - total_nonsyn)/numpy.std(total_nonsynsample)

        total_LOF = total_seg['nonsens']+total_seg['splice']
        total_LOFsample = total_mut['nonsens']+total_mut['splice']
        z_LOF = (numpy.mean(total_LOFsample) - total_LOF)/numpy.std(total_LOFsample)
        return [CI0,CI1,z_overall,z_nonsyn,z_LOF]


    def get_aa_conf_score(self,aminoacid,offset,aa_sample,aa_total_seg,return_stuff=False):
        #This function is used specifically in testing of coding transcripts. We find using the overall number of variants for an AA, the rank of observed variants in the expected distribution
        #Since we use this for offset1,2,3 so we need to make sure that we sum over the previous samplings (aka for all offsets) 
        #sample is added to the total_mut and is returned, if return_stuff=True then we return the rank
        amino_mapping  = returns_amino_mapping()
        if return_stuff == False:
            total_mut = numpy.zeros(self.nsamp) #We save the total simulated variants here
            for alpha in self.keys_kmer.keys():
                possible_mut = self.multisamp(alpha)
                index = 0
                for beta in self.keys_kmer[alpha]:
                    if not(alpha==beta):
                        length = len(alpha)
                        extract = [(length/2)+1-offset, (length/2)+2-offset, (length/2)+3-offset]
                        aa1 = amino_mapping[''.join([ alpha[extract[0]], alpha[extract[1]], alpha[extract[2]]])]
                        aa2 = amino_mapping[''.join([ beta[extract[0]], beta[extract[1]], beta[extract[2]]])]
			if len(aminoacid)==1:
                            if not(aa1==aa2) and (aa1 == aminoacid):
                                total_mut += possible_mut[:,index]
			else:
                            if not(aa1==aa2) and (aa1 in aminoacid) and not(aa2 in aminoacid):
                                total_mut += possible_mut[:,index]
                    index += 1
            if not aa_sample == None:
                total_mut += aa_sample
            return total_mut
        else:
	    aa_sample.sort()
	    CI0 = int(aa_sample[int(self.nsamp*0.025)])
            CI1 = int(aa_sample[int(self.nsamp*0.975)])
            z_aa = (numpy.mean(aa_sample) - aa_total_seg)/numpy.std(aa_sample)
            aa_median = numpy.median(aa_sample)
	    return [CI0,CI1,aa_median,z_aa]
		
def returns_amino_mapping():
    #This function returns the mapping of codon to amino acids
    handle = file("/project/voight_subrate/avarun/Research/mutation_rate/codon_aa_table","r")
    content = handle.readlines()
    amino_mapping = {}
    for entry in content:
        entry = entry.rstrip('\n').split('\t')
        nucs = entry[1].split(',')
        for nuc in nucs:
            amino_mapping[nuc] = entry[0]
    return amino_mapping
def returns_codon_mapping():
    #This function returns the mapping of amino acid to codons
    handle = file("/project/voight_subrate/avarun/Research/mutation_rate/codon_aa_table","r")
    content = handle.readlines()
    amino_mapping = {}
    for entry in content:
        entry = entry.rstrip('\n').split('\t')
	codons = entry[1].split(',')
	if not entry[0] in amino_mapping:
	    amino_mapping[entry[0]] = []
	for codon in codons:
            amino_mapping[entry[0]].append(codon)
    return amino_mapping


def return_distributed_variants_igr(total_var,count,dict_rate):
    #This function will distribute total_var using the dict_rate and count of each kmer in IGR
    prob = {}
    for key in dict_rate['combined']:
	prob[key] = dict_rate['combined'][key]*count[key[0]]*1.0
    sum_prob = sum(prob.values())
    probability_mut = []
    keys_list = []
    for key in prob:
	prob[key] = prob[key]/sum_prob
	probability_mut.append(prob[key])
	keys_list.append(key)
    flag_correct =  False
    while flag_correct == False:
	flag_correct = True
        distribution = numpy.random.multinomial(total_var,probability_mut,1)[0]
        distribution_variants = []
        for index in range(len(distribution)):
	    variant = keys_list[index]
	    if distribution[index]>count[variant[0]]: #More variants distributed then the occurence itself.
		print "We are distributing variants once again. The simulated region is too small. Try giving a bigger sequence"
		flag_correct = False
	    else:
	        for i in range(distribution[index]):	
	            distribution_variants.append(variant)
    return distribution_variants

def get_denovo_prob(model,folded_st,only_folded="F"):
    #This function gets the denovo probability as a dictionary
    #folded_st then reversecomplement get the same rate as the folded one
    #onlyfolded returns only the folded rate (not the reverse complement)
    dict_rate = {}
    dict_rate['AFR'] = {}
    dict_rate['EUR'] = {}
    if model=="Constant":
	filename = "denovo_sub_prob_1mer_folded"
	handle = file(filename,"r")
	content = handle.readlines()
	for entry in content:
	    entry = entry.split('\t')
	    dict_rate['AFR'][(entry[0],entry[1])] = .4/100000000 #Because equal for the 3 changes
            dict_rate['EUR'][(entry[0],entry[1])] = .4/100000000
	    if only_folded == "F":
                dict_rate['AFR'][(reverse_complement(entry[0]),reverse_complement(entry[1]))] = .4/100000000 #Because equal for the 3 changes
                dict_rate['EUR'][(reverse_complement(entry[0]),reverse_complement(entry[1]))] = .4/100000000
    else:
	filename = "denovo_sub_prob_"+model+folded_st
	handle = file(filename,"r")
	content = handle.readlines()
        for entry in content:
            entry = entry.split('\t')
	    alpha = entry[0]
	    beta = entry[1]
            dict_rate['AFR'][(alpha,beta)] = float(entry[3])
            dict_rate['EUR'][(alpha,beta)] = float(entry[5])
	    if folded_st == "_folded":
		alpha_r = reverse_complement(alpha)
		beta_r = reverse_complement(beta)
		if only_folded == "F":
                    dict_rate['AFR'][(alpha_r,beta_r)] = dict_rate['AFR'][(alpha,beta)] 	
                    dict_rate['EUR'][(alpha_r,beta_r)] = dict_rate['EUR'][(alpha,beta)]
    return dict_rate

