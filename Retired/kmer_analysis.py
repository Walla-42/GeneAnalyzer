import numpy as np
import matplotlib.pyplot as plt

sequence = "atagagaggagagaggaatatatatattaagagagaggatgattgatgatgatagcgcgcgcgcgcgcgcgcgcgc"

def kmer_analyzer(sequence,k):
    kmer_dic = {}
    for i in range(0,len(sequence)):
        kmer = sequence[i:i+k]
        if kmer in kmer_dic:
            kmer_dic[kmer] += 1
        else:
            kmer_dic[kmer] = 1
    for kmer in kmer_dic:
        kmer_dic[kmer] = kmer_dic[kmer]/len(kmer_dic)
    return kmer_dic

def kmer_graph(kmer_dictionary):
    kmer_list = []
    kmer_frequency = []
    for kmer in kmer_dictionary:
        kmer_list.append(kmer)
        kmer_frequency.append(kmer_dictionary[kmer])
    plt.bar(kmer_list, kmer_frequency)
    plt.show()


kmer_graph(kmer_analyzer(sequence,3))



# Things to do:
# - Figure out how to make this useful - In genomics Kmers analysis is useful for sequencing data
# - Allow the user to analyze sequence data by kmer analysis
# - Sequence frequency vs number of times it appears in the sequenced data. 
# - In Sequence assembly, Kmers are the specific contigs 