import sys
sys.path.append("/project/voight_subrate/avarun/Research/mutation_rate")
from modules import *
from function_wrapper import *
import random
import pdb

def get_seq_context_variant(chr_name,position,before,after):
    #Chr_name: chromosome name 
    #position: position of variant
    #before after: nucleotides before or after the variant
    files = file_handles("/project/voight_datasets/hg19")
    chr_name = str(chr_name)
    position = int(position)
    before = int(before)
    after = int(after)
    str_around = files.fetch_str_var(chr_name,position,before,after) #string around the variant including it
    #print "sequence context "+str(before)+" basepairs before and "+str(after)+" after the variant including the middle position"+''.join(str_around)
    #return ''.join(str_around)
    return ''.join(str_around)

def get_seq_context_interval(chr_name,start_p,end_p,before_after):
    #Chr_name: chromosome name 
    #start_p: start position of interval
    #end_p: end postion of interval
    #before_after: nucleotides before or after the start/end (needed for padding)
    files = file_handles("/project/voight_datasets/hg19")
    chr_name = str(chr_name)
    start_p = int(start_p)
    end_p = int(end_p)
    before_after = int(before_after)
    acceptable_chr_name = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]
    if chr_name in acceptable_chr_name:
        seq = ''.join(files.get_string(chr_name, start_p, end_p, before_after, before_after))
        return seq
    else:
	print "Chromosome entered is not valid"

