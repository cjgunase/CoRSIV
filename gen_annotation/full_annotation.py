"""
Annotates a genome chromosome by chromosome
"""


import numpy as np
import sys
import math
import pandas as pd
from collections import defaultdict
import multiprocessing as mp



# Uses a binary search to find the exon start or end that has the lowest absolute distance to the bin_location
# This assumes that exon starts/ends are unique, ie we can use the index of that position to find the rest
# of the data in the corresponding row
def get_closest(bin_location):
   
        # Gene starts and ends, strand independent
        
        # Copying the list is an expensive operation in the algorithm
        # starts_copy = starts#list(starts)
        # ends_copy = ends#list(ends)
        starts_copy = list(starts)
        ends_copy = list(ends)
        forward_genes_ends = []
        forward_genes_starts = []
        reverse_genes_ends = []
        reverse_genes_starts = []

        true_tx_starts = []
        true_tx_ends = []

        #for i in range(10):
        while len(true_tx_starts) < 3 or len(true_tx_ends) < 3:
                # Relative start
                start = inexact_binary_search(starts_copy, bin_location)
                end = inexact_binary_search(ends_copy, bin_location)
                # Stop searching if we have gone too far
                if min(abs(start-bin_location),abs(end-bin_location)) > CUTOFF:
                        break
                # Find the closer of the nearest end locations of genes
                if abs(start-bin_location) < abs(end-bin_location):
                        starts_copy.remove(start)
                        # Determine if this is a true start or not
                        start_index = gene_data[gene_data["txStart"] == start].index[0]
                        gene = gene_data.iloc[start_index]
                        strand = gene["strand"]
                        gene_name = gene["name2"]
                        if strand == "+": # True start
                                if gene_name not in forward_genes_starts and len(true_tx_starts) < 3:
                                        forward_genes_starts.append(gene_name)
                                        true_tx_starts.append([start,strand, "starts", gene_name])
                        else: # Actually a gene end
                                if gene_name not in reverse_genes_ends and len(true_tx_ends) < 3:
                                        reverse_genes_ends.append(gene_name)
                                        true_tx_ends.append([start, strand, "starts", gene_name])
                else:
                        ends_copy.remove(end)
                        end_index = gene_data[gene_data["txEnd"] == end].index[0]
                        gene = gene_data.iloc[end_index]
                        strand = gene["strand"]
                        gene_name = gene["name2"]
                        if strand == "+": # True gene end
                                if gene_name not in forward_genes_ends and len(true_tx_ends) < 3:
                                        forward_genes_ends.append(gene_name)
                                        true_tx_ends.append([end, strand, "ends", gene_name])                   
                        else: # Actually a gene start
                                if gene_name not in reverse_genes_starts and len(true_tx_starts) < 3:
                                        reverse_genes_starts.append(gene_name)
                                        true_tx_starts.append([end,strand, "ends", gene_name])
        gs = []
        ge = []
        gs_names = []
        ge_names = []
        for gene in true_tx_starts:
                location = gene[0]
                strand = gene[1]
                name = gene[3]
                distance = float("inf")
                if strand == "+":
                        distance = bin_location - location
                else:
                        distance = location - bin_location
                gs_names.append(name)
                gs.append(distance)

        for gene in true_tx_ends:
                location = gene[0]
                strand = gene[1]
                name = gene[3]
                distance = float("inf")
                if strand == "+":
                        distance = bin_location - location
                else:
                        distance = location - bin_location
                ge_names.append(name)
                ge.append(distance)

   # Append None values to fill in the gaps, makes sure everything is the same length
        gs = gs + ([float("inf")]*(num_genes_to_find-len(gs)))
        ge = ge + ([float("inf")]*(num_genes_to_find-len(ge)))
        gs_names = gs_names + (["NA"]*(num_genes_to_find-len(gs_names)))
        ge_names = ge_names + (["NA"]*(num_genes_to_find-len(ge_names)))

        current_bin_data = {}
        current_bin_data["location"] = bin_location
        current_bin_data["GS1"] = gs[0]
        current_bin_data["GS2"] = gs[1]
        current_bin_data["GS3"] = gs[2]
        current_bin_data["GE1"] = ge[0]
        current_bin_data["GE2"] = ge[1]
        current_bin_data["GE3"] = ge[2]
        current_bin_data["GS1 Name"] = gs_names[0]
        current_bin_data["GS2 Name"] = gs_names[1]
        current_bin_data["GS3 Name"] = gs_names[2]
        current_bin_data["GE1 Name"] = ge_names[0]
        current_bin_data["GE2 Name"] = ge_names[1]
        current_bin_data["GE3 Name"] = ge_names[2]
        
        return current_bin_data

# Returns the element that is closest in abs value in alist (sorted ascending) to target
# TODO: modify this so if the element in the list is more than 10E6 Bp away we break and return None
def inexact_binary_search(alist, target):
    # alist - a sorted list, ascending order (left is low right is high)
    # target - the item we are looking for
    l = len(alist)
    # Base case
    if l == 1:
        return alist[0]
    # Recursive case
    else:
        # Split list in half
        L = alist[:l/2]
        R = alist[l/2:]
        if abs(L[-1] - target) < abs(R[0]-target):
            return inexact_binary_search(L,target)
        elif L[-1] == R[0]:
                if target <= L[-1]:
                        return inexact_binary_search(L,target)
                else:
                        return inexact_binary_search(R,target)
        else:
                        return inexact_binary_search(R,target)



# For parallelization
def get_overlapping_genes(bin_name):
    # Get bin coordinates
    right_edge = int(bin_name[bin_name.index("_")+1:])
    left_edge = right_edge-(bin_size-1)
    bin_overlapping_genes = set()
    
    # Look through each gene, this is slow, can we parallize it?
    for gene_index in range(len(gene_data)):
        # Get transcription coordinates
        tx_start = gene_data["txStart"].iloc[gene_index]
        tx_end = gene_data["txEnd"].iloc[gene_index]
        # Does the bin overlap with the gene?
        if overlaps(left_edge,right_edge,tx_start,tx_end):
            gene_name = gene_data["name2"].iloc[gene_index]
            bin_overlapping_genes.update({gene_name})
    return bin_overlapping_genes


# returns true if any part of list a overlaps with list b
# Input: leftmost element of a, rightmost element of a, same stuff for list b
def overlaps(a_L,a_R, b_L,b_R):
    # If the left index is in the list b or the right index is in the list b
    return ((a_L <= b_R) and (a_L >= b_L)) or ((a_R <= b_R) and (a_R >= b_L))


# # Computes any bins that should overlap a gene transcription, return as a list of bins
def compute_overlapping_bins(start, end, bin_size = 100):
    """
    start: transcription start
    end:   transcription end
    """
    # Note: a bin is defined by its rightmost coordinate
    first_bin =  start - (start%bin_size) + bin_size
    last_bin = end - ((end-1)%bin_size)
    return range(first_bin, last_bin + bin_size, bin_size)

def compute_overlap_data(gene_data):
    overlap_data = defaultdict(lambda:set()) # Mapping from bin name to set of overlapping genes
    for gene in list(gene_data.itertuples()):
        # Extract gene data
        gene_name = gene[13]
        tx_start = gene[5]
        tx_end = gene[6]

        # Compute any bins that should overlap this gene
        overlapping_bins = compute_overlapping_bins(tx_start, tx_end)

        # Update bin data
        for bin_ in overlapping_bins:
            overlap_data[bin_].update({gene_name})

    return overlap_data


# CONSTANTS
CUTOFF = 1000000 # Max distance to a gene to consider it in the primary analysis
num_genes_to_find=3 # Number of genes to find in primary analysis
CLOSE_CUTOFF = 2500 # Distance to a gene to consider it "close"
bin_size = 100 # Bin size used in binning analysis

# Parameters from command line
cpg_file = sys.argv[1] #e.g. CpGsInChrome1.csv
gene_file = sys.argv[2]#e.g. chrm1genes this data comes from UCSC
chrm_name = sys.argv[3] # e.g. chr1
chromosome_length = int(sys.argv[4])
#chromosome_length = 1234

# cpg_file = "all_bins_CpG_counts_multi_chr22.csv"
# gene_file = "chr22.gene"
# chrm_name = "chr22" chrY_ 159161600 58371970 chrY_58371800
# chrY_57216600
#      58371970
# Global
cpg_bins = pd.read_csv(cpg_file)# e.g. CpGsInChrome1.csv
gene_data = pd.read_table(gene_file) # e.g. chrm1genes


# Get transcription coordinates and sort
# Global
starts = list(gene_data["txStart"])
ends = list(gene_data["txEnd"])
starts.sort()
ends.sort()

# Create our new dataframe
annotated_bins = pd.DataFrame.from_dict({})
all_bin_names = []

# Add all chromosomal bins to the dataframe
for loc in range(bin_size,chromosome_length+bin_size, bin_size):
  bin_name = chrm_name + "_" + str(loc)
  all_bin_names.append(bin_name)

# Append CpG count data to the annotation
# First, create a mapping
bin_names = cpg_bins["Bin Name"]
CpG_counts = cpg_bins["CpG Count"]
bin_cpg__count_mapping = dict(zip(bin_names,CpG_counts))

all_bin_cpg_counts = []
for bin_ in all_bin_names:
    if bin_ in bin_cpg__count_mapping:
        all_bin_cpg_counts.append(bin_cpg__count_mapping[bin_]) # CpGs found
    else:
        all_bin_cpg_counts.append(0) # CpgS not found

print "Max is ",max(all_bin_cpg_counts)
annotated_bins["Bin Name"] = all_bin_names
annotated_bins["CpG Count"] = all_bin_cpg_counts

print annotated_bins["Bin Name"]
print "CpG Counts"
print annotated_bins["CpG Count"]



#Maps from a bin to another mapping that contains coordindates, gs, and ge data, global

location_data = []

for bin_ in list(annotated_bins["Bin Name"]):
    loc = int(bin_.split("_")[1])
    location_data.append(loc)

# Parallelization
pool = mp.Pool()
bin_data = pool.map(get_closest, location_data)
pool.close() 
pool.join() # close the pool and wait for the work to finish 


coordindates = []
g1_s = []
g1_s_names = []
g1_e = []
g1_e_names = []
g2_s = []
g2_s_names = []
g2_e = []
g2_e_names = []
g3_s = []
g3_s_names = []
g3_e = []
g3_e_names = []


for datem in bin_data:
        bin_name = datem["location"]
        bin_coordinates = chrm_name+":"+str(bin_name-(bin_size-1))+"-"+str(bin_name)
        
        #bin_size - 1 since geneticists don't like to index things starting at 0 :(
        coordindates.append(bin_coordinates)
        g1_s.append(datem["GS1"])
        g2_s.append(datem["GS2"])
        g3_s.append(datem["GS3"])

        g1_s_names.append(datem["GS1 Name"])
        g2_s_names.append(datem["GS2 Name"])
        g3_s_names.append(datem["GS3 Name"])    

        g1_e.append(datem["GE1"])
        g2_e.append(datem["GE2"])
        g3_e.append(datem["GE3"])

        g1_e_names.append(datem["GE1 Name"])
        g2_e_names.append(datem["GE2 Name"])
        g3_e_names.append(datem["GE3 Name"])


# Add data to the data frame
annotated_bins["UCSC Coordinates"] = coordindates + ["NA"]*(len(annotated_bins)-len(coordindates))

# # Starts
annotated_bins["GS1 Name"] = g1_s_names + ["NA"]*(len(annotated_bins)-len(g1_s_names))
annotated_bins["GS1 Distance"] = g1_s + ["NA"]*(len(annotated_bins)-len(g1_s))
annotated_bins["GS2 Name"] = g2_s_names + ["NA"]*(len(annotated_bins)-len(g2_s_names))
annotated_bins["GS2 Distance"] = g2_s + ["NA"]*(len(annotated_bins)-len(g2_s))
annotated_bins["GS3 Name"] = g3_s_names + ["NA"]*(len(annotated_bins)-len(g3_s_names))
annotated_bins["GS3 Distance"] = g3_s + ["NA"]*(len(annotated_bins)-len(g3_s))

# Ends
annotated_bins["GE1 Name"] = g1_e_names + ["NA"]*(len(annotated_bins)-len(g1_e_names))
annotated_bins["GE1 Distance"] = g1_e + ["NA"]*(len(annotated_bins)-len(g1_e))
annotated_bins["GE2 Name"] = g2_e_names + ["NA"]*(len(annotated_bins)-len(g2_e_names))
annotated_bins["GE2 Distance"] = g2_e + ["NA"]*(len(annotated_bins)-len(g2_e))
annotated_bins["GE3 Name"] = g3_e_names + ["NA"]*(len(annotated_bins)-len(g3_e_names))
annotated_bins["GE3 Distance"] = g3_e + ["NA"]*(len(annotated_bins)-len(g3_e))

#34063,34063,chr21_6812800,0,chr21_6812800,chr21:6812601-6812800,CH507-210P18.1,23209.0,CH507-210P18.4,-45738.0,MIR8069-1,-46370.0,CH507-210P18.1,503.0,MIR8069-1,-46456.0,CH507-210P18.4,-84463.0,NA,set([])
# Step 2
# Find "Close Genes"
close_genes = defaultdict(lambda:set())# mapping from bin index to list of close genes

# Find close genes
for idx in range(len(annotated_bins)):
    for i in range(1,num_genes_to_find+1):
        for x in ["GS","GE"]:
                distance = abs(float(annotated_bins[x+str(i)+ " Distance"].iloc[idx]))
                if distance <= CLOSE_CUTOFF:
                    gene_name = annotated_bins[x+str(i)+" Name"].iloc[idx]
                    close_genes[idx].update({gene_name})

# convert dictionary to a list and append it to the dataframe
sorted_keys = close_genes.keys()
sorted_keys.sort()
close_gene_list = []
for i in range(len(annotated_bins)):
    if len(close_genes[i]) == 0: # Empty set
        close_gene_list.append("NA")
    else:
        close_gene_list.append(",".join(list(close_genes[i])))

annotated_bins[" <= 2.5k bp Genes"] = close_gene_list + ["NA"]*(len(annotated_bins)-len(close_gene_list))

# Parallelization baby yeah

# pool = mp.Pool()
# overlap_data = pool.map(get_overlapping_genes, list(annotated_bins["bin_name"]))
# pool.close() 
# pool.join() 
# annotated_bins["Overlapping Genes"] =  overlap_data + ["NA"]*(len(annotated_bins)-len(overlap_data))


# Go through bins and created data for overlapping genes
overlapping_genes = []
overlap_data = compute_overlap_data(gene_data)
for bin_ in annotated_bins["Bin Name"]:
    bin_ = bin_.split("_")[1]
    overlapping_genes.append(",".join(list(overlap_data[int(bin_)])))

annotated_bins["Overlapping Genes"] = overlapping_genes



#Save the data
filename = chrm_name + "_annotated.csv"
annotated_bins.to_csv(filename)






