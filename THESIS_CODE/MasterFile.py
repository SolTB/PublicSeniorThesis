'''
SENIOR THESIS WORK
Author: Sol Taylor-Brill
Version: 10/13/19

This program combines several files: the original RNA-Seq data, a list of gene names in the order they appear in that file
and a name map that contains common names and gene locus IDs into one mega file to make the data easier to work with.
'''

def main():
	print("HELLO")
	RNA_Names = "RNA_NameList.txt"
	RNA_Excel = "VehicleCBDStd2.txt"
	Name_Map = "Combine_Name_Map.txt"
	OutputFile = ""

	RNA_list = ListMaker(RNA_File) #a list of genes in the same order as the excel file (to make formatting easier)
	IDlist, IDDict, NameDict = IDDict_maker(Name_Map) #makes name>ID/ID>name dictionaires
	return

#This function just goes through the RNA Name List file and makes a list of all the genes from the RNA Seq File
def ListMaker(file):
	f = open(file, 'r')
	print("opening ", file)
	gene_list = []
	for l in f:
		line = str(l)
		items = line.split()
		for gene in items:
			gene_list.append(gene)
	return gene_list

'''
This function goes through the name map files and makes ID>Name/Name>ID dictionaries and makes a list of gene
IDs for which the common name is known (ie gene appears in the original source file) 
'''
def IDDict_maker(MapFile):
	IDDict = {} #ID-> common name dictionary
	NameDict = {} #common name -> ID dictionary
	Known_IDset = set() #a list of all

	NameMap = open(MapFile,'r')
	print("opening ", MapFile)

	for l in NameMap:
		line = str(l) #line is equal to string
		itemsinline = line.split() #splits items up by tabs OR spaces and makes a list
		ID = itemsinline[0]
		name = itemsinline[1]
		if name != 'NA': #only includes genes for which the file has a known common name
			IDDict[ID] = name #adds to ID>name dictionary
			NameDict[name] = ID #adds to name>ID dictionary
			Known_IDset.add(ID) #adds gene to set of known genes
	Known_IDlist = list(Known_IDset) #makes the set a list
	return Known_IDlist, IDDict, NameDict

def excel_reader(RNA_Excel, RNA_list):
	basemean_Dict = {} #ID>basemean
	log2_Dict = {} #ID>log2
	lfcSE_Dict = {} #ID>lfcSE
	stat_Dict = {} #ID>stat
	pval_Dict = {} #ID>pval
	padj_Dict = {} #ID>padj

	f = open(RNA_Excel, 'r')
	print("opening ", RNA_Excel)
	listindex = 0 #just to keep track of what line we're to index RNA_name

	for l in f:
		line = str(l) #line is equal to string
		line_list = line.split() #splits items up by tabs OR spaces and makes a list
		GeneID = RNA_list[listindex]

		basemean = line_list[1]
		basemean_Dict[GeneID] = basemean

		log2 = line_list[2]
		log2_Dict[GeneID] = log2

		lfcSE = line_list[3]
		lfcSE_Dict[GeneID] = lfcSE

		stat = line_list[4]
		stat_Dict[GeneID] = stat

		pval = line_list[5]
		pval_Dict[GeneID]

		padj = line_list[6]
		padj_Dict[GeneID] = padj

		listindex = listindex + 1 #to track Gene ID for the excel file
	return basemean_Dict, log2_Dict, lfcSE_Dict, stat_Dict, pval_Dict, padj_Dict


main()