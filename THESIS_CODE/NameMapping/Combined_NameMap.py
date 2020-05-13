'''
SENIOR THESIS WORK
Author: Sol Taylor-Brill
Version: 10/12/19

This program takes input from ID-> common name info from BioCyc and STRING and combines them
'''

def main():
	BioCyc_Map = "BioCyc_NameMap.txt" #Name Map file made from BioCyc Data. Contains ID and common name
	STRING_Map = "STRING_NameMap.txt" #Name Map file made from STRING Data. Contains ID and common name

	BioCyc_Unknowns = "BioCyc_UnknownIDs.txt" #txt file containing all of the unknown gene IDs in the BioCyc file
	BC_UK_list = ListMaker(BioCyc_Unknowns)
	STRING_Unknowns = "STRINGUnknownIDs.txt" #txt file containing all the unknown gene IDs in the STRING file
	STRING_UK_list = ListMaker(STRING_Unknowns)

	RNA_File = "RNA_NameList.txt"
	RNA_list = ListMaker(RNA_File)

	BioCyc_IDlist, BioCyc_IDDict, BioCyc_NameDict = IDDict_maker(BioCyc_Map)
	STRING_IDlist, STRING_IDDict, STRING_NameDict = IDDict_maker(STRING_Map)

	OutputFile = "Combine_Name_Map.txt"
	Combined_FileMaker(RNA_list,BC_UK_list, STRING_UK_list, BioCyc_IDDict, STRING_IDDict, OutputFile)
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


#This function attempts to fill in the missing data from the BioCyc data with the STRING data
def Combined_FileMaker(RNA_list,BC_UK_list, STRING_UK_list, BioCyc_IDDict, STRING_IDDict, OutputFile):
	Combined_UK_set = set()
	numinSTRING = 0 #number of genes in STRING but not BioCyc
	f = open(OutputFile,'w')
	print("writing to ", OutputFile)

	for gene in RNA_list:
		if gene not in BC_UK_list: #if it's in the BioCyc data file
			name = BioCyc_IDDict[gene] 
			f.write(gene + '\t' + name + '\n') #write to output file with data from BC file
		elif gene in BC_UK_list and gene not in STRING_UK_list: #if its not in the BC data but it IS in the STRING data
			name = STRING_IDDict[gene] 
			f.write(gene + '\t' + name + '\n') #write to output file with data from STRING file
			numinSTRING = numinSTRING + 1 #counting number of genes in STRING but not BC
		else: #if not in either
			f.write(gene + '\t' + "UNKNOWN" + '\n') #say it's unknown
			Combined_UK_set.add(gene)

	Combined_UK_list = list(Combined_UK_set)
	UnknownNum = (len(Combined_UK_list)/len(RNA_list)) *100

	print("Number of genes in STRING but NOT BioCyc: ", str(numinSTRING))
	print("Number of genes in neither: " + str(len(Combined_UK_list)))
	print("Number of genes in RNA-Seq Output: " + str(len(RNA_list)))
	print("Percent still unknown: ", str(UnknownNum), "%")
	
	return



main()


