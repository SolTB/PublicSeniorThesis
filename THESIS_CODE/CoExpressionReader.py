'''
This program takes files of normalized counts from the DESeq run using gff files for pO157, pOSAKI, and the chromosome
Within each file there should be 7 rows: gene ID, Control1, Control2, Control3, CBD1, CBD2, and CBD3. (in the plasmid
txt files the three CBD files are first and the three control files are second). The program should read the three text files
make both a CBD and control gene->(expression1,expression2,expression3) dictionary. 


It will use these dictionaries to run a pearson correlation coefficient test on every gene x every gene. 
If the pearson correlation is above a certain value and has a significant (<0.01) p value, an edge will 
be created between those two genes. 

I'm going to try this two different ways: A) make a CBD network and a control network and compare them 
and B) somehow make 1 network that uses differential expression as the basis, but I haven't figured out how to make
this work with a pearson correlation coefficient. As a preliminary, I might just do like (ControlX-CBDX) for all three
but this seems statistically bad like there must be a better way to do this that somehow captures the difference between 
1 group vs another (Potentially ask Anna about this)

Class: THESIS
Author: Sol Taylor-Brill
Date of last edit: 12/11/19
'''
import math
import statistics
import numpy as np
import scipy
from scipy import stats
import time


def main():
	old_pO157_f = "pO157NormalizedCounts.txt" ## THIS FILE IS FORMATTED INCORRECTLY AND IS FIXED
	CBDDict_List = [] #list of all three TRIMMED CBD dictionaries (2 plasmids + 1 chromosome)
	ControlDict_List = [] #list of all three TRIMMED control dictionaries
	Master_Gene_Set = set() #all of the genes in the interactome

	## Normalized Counts Files ##
	chromosome_f = "ChromosomeNormalizedCounts.txt" ## CHROMOSOME FILE
	pOSAKI_f = "pOSAKI_Counts.txt" ## pOSAKI FILE
	pO157_f = "PO157_Counts.txt" ## FIXED pO157 FILE

	file_list = [chromosome_f, pOSAKI_f, pO157_f] #list of chromosome and plasmid files
	
	## Re-formatting and Reading Text Files ##
	#ReWrite(old_pO157_f, pO157_f) #reformatted plasmid file

	#For all three files: makes a gene->[expression1, expression2, expression3] dictionary for CBD and control
	for file in file_list:
		ControlDict, CBDDict, TrimmedControlDict, TrimmedCBDDict, gene_set = FileReader(file) #Trimmed files don't include genes will all the same values (ruins pearson)
		ControlDict_List.append(TrimmedControlDict) #adds TRIMMED dictionary to list of dictionaries
		CBDDict_List.append(TrimmedCBDDict) #adds TRIMMED dictionary to list of dictionaries

	MasterCBD_Dict = {**CBDDict_List[0], **CBDDict_List[1], **CBDDict_List[2]} #combines all 3 TRIMMED CBD Dicts
	MasterControl_Dict = {**ControlDict_List[0], **ControlDict_List[1], **ControlDict_List[2]} #combines all 3 TRIMMED Control Dicts

	## Correlation Analyses ##
	pear_cutoff = 0.3 #(0.1-0.3 = weak (exclude?), 0.3-0.5 = moderate, 0.5-1.0 = strong)
	p_val_cutoff = 0.01
	start = time.time()
	edge_set = Correlation_Testing(CBDDict_List[0], CBDDict_List[0], pear_cutoff, p_val_cutoff)
	end = time.time()
	change = end - start
	print("Pearson Correlation Time: ", change)

	file_Chrom2 = "ChromxChromCBD.txt.txt"
	Print_Edges(edge_set, file_Chrom2)

	return


## REFORMATTING AND READING TEXT FILES ##


#reformats the 2 plasmid normalized count files so they have the same format as the chromosome file
def ReWrite(file, outputf):
	Count = 0 #because header is formatted differently than the rest of the file
	Gene_List = []
	CBD1_List = []
	CBD2_List = []
	CBD3_List = []
	Control1_List = []
	Control2_List = []
	Control3_List = []

	f = open(file, 'r')
	print("opening ", file)

	for l in f:
		line = str(l)
		items = line.split()
		if  Count == 0:
			Gene_List.append("Gene Name")
			CBD1_List.append(items[0])
			CBD2_List.append(items[1])
			CBD3_List.append(items[2])
			Control1_List.append(items[3])
			Control2_List.append(items[4])
			Control3_List.append(items[5])
			Count = Count + 1 #add 1 to count so it runs through the rest of the file
		else:
			Gene_List.append(items[0])
			CBD1_List.append(items[1])
			CBD2_List.append(items[2])
			CBD3_List.append(items[3])
			Control1_List.append(items[4])
			Control2_List.append(items[5])
			Control3_List.append(items[6])

	f.close() #closes input file

	o_f = open(outputf, 'w') #opens output file
	print("opening", o_f)

	for i in range(len(Gene_List)): #writes info in the order that it is in the chromosome file
		output_line = str(Gene_List[i]) + '\t' + str(Control1_List[i]) + '\t' + str(Control2_List[i]) + '\t' + str(Control3_List[i]) + '\t' + str(CBD1_List[i]) + '\t' + str(CBD2_List[i]) + '\t' + str(CBD3_List[i]) + '\n'
		o_f.write(output_line) #writes line to file

	o_f.close() #output file closed
	print("Wrote to ", o_f) #function is finished printing

	return

#Input:Normalized counts text file (format = gene, control1, control2, control3, cbd1, cbd2, cbd3)
#Output: ControlDict (gene -> (control1, control2, control3)) dictionary & CBDDict (gene -> (cbd1, cbd2, cbd3)) dictionary
# & TrimmedCBDDict/TrimmedControlDict which don't include genes that have the same expression values (since it will return a pearson error)
def FileReader(file):
	count = 0 #this is just to make sure it doesn't add the header to the file
	CBDDict = {} # gene: [CBD1, CBD2, CBD3] dictionary
	TrimmedCBDDict = {} #this dictionary doesn't include genes that have the same expression for all three values to avoid pearson errors
	ControlDict = {} #gene: [Control1, Control2, Control3] dictionary
	TrimmedControlDict = {}
	gene_set = set() #set of genes in the file
	
	f = open(file, 'r')
	print("opening ", file)
	
	for l in f: #for every line in the file
		if count == 0:
			count = count + 1 #this line is just to get it to skip the header
		else: #so that the header isn't added
			line = str(l)
			items = line.split()
			gene = items[0] #the first item is gene
			control1 = items[1]
			control2 = items[2]
			control3 = items[3]
			ControlDict[gene] = (float(control1), float(control2), float(control3)) #adds tuple with all three control expression values to gene -> control dict
			if control1 != control2 or control2 != control3: #if at least two values are not the same
				TrimmedControlDict[gene] = (float(control1), float(control2), float(control3))
			
			cbd1 = items[4]
			cbd2 = items[5]
			cbd3 = items[6]
			CBDDict[gene] = (float(cbd1), float(cbd2), float(cbd3)) #adds tuple of all three cbd values to the gene -> expression dict
			if cbd1 != cbd2 or cbd2 != cbd3: #if at least two are different
				TrimmedCBDDict[gene] = (float(cbd1), float(cbd2), float(cbd3))
	return ControlDict, CBDDict, TrimmedControlDict, TrimmedCBDDict, gene_set



## TESTING COEXPRESSION BETWEEN GENES ##


#Input: EITHER MasterCBD_Dict or MasterControl_Dict ((gene->expressionlist) dictionary), p_cutoff (value below which pearson coefficient not considered important )
#Output: EITHER edges in CBD-only interactome or Control-only interactome
def Correlation_Testing(DictA, DictB, pear_cutoff, p_val_cutoff):
	edge_set = set() #set of all edge tuples (u,v) with a p-val<0.01
	seen_edge_set = set() #saves edges that we've already seen
	print("#Genes = ", len(DictA), len(DictB))

	for u in DictA: #all nodes in the dictionary
		a = DictA[u] #I'm defining a and u array here so that I only have to do it once for each node
		u_array = np.asarray(a, dtype=None, order=None)
		for v in DictB: #against all other nodes in the dictionary
			if u != v and (v,u) not in seen_edge_set: #don't test a gene by itself & don't go over a pair more than once
				b = DictB[v] #dictionary b -> expression list
				v_array = np.asarray(b, dtype=None, order=None) #turns list into an array
				
				pearson, p_val = scipy.stats.pearsonr(u_array, v_array) #compare expression list of one gene vs another

				seen_edge_set.add((v,u)) #adds inverse edge to seen edges so we don't go through it
				seen_edge_set.add((u,v)) #adds edge so we don't go through it 

				if p_val <= p_val_cutoff:
					edge_set.add((u,v, pearson)) #adds edge (node1,node2, pearson) as tuple to edgeset
	print("#Edges: ", len(edge_set))
	return edge_set


def Print_Edges(edge_set, file):
	f = open(file, 'w') #opens output file
	print("opening", file)

	for edge in edge_set:
		fileline = str(edge[0]) + '\t' + str(edge[1]) + '\t' + str(edge[2]) + '\n'
		f.write(fileline)

	print("Wrote to ", file)
	return

main()




