'''
SENIOR THESIS WORK
Author: Sol Taylor-Brill
Version: 02/03/20

This program takes an input file from NCBI containing a variety of info about proteins in E. coli O157 
parses through it and outputs a file containing just ID and common name for ease of use
'''

def main():
	file_name = "gene_result.txt"
	EXCEL_file = "DEDATA.txt"
	Output_file = "NameMapECS.txt"

	IDlist, IDlogDict, IDpvalDict = EXCEL_filereader(EXCEL_file)
	STRING_IDlist, IDDict, NameDict, NotesDict = STRING_filereader(file_name)
	final_filewriter(IDlist, STRING_IDlist, IDDict, NotesDict, Output_file)
	return


#Reads text file version of excel sheet and makes list of geneIDs (namelist) and dictionaries
# name-> logfold change (namelogDict) and name->p-value (namepvalDict)
#This function is from "CBDEdit.py"
#In this function name is not the common name!! (name = ID here)
def EXCEL_filereader(file):
	f=open(file,'r')
	print("Opening " + file)
	nameset = set()
	namelogDict = {}
	namepvalDict = {}
	for l in f:
		line = str(l)
		itemsinline = line.split()
		original_name = itemsinline[0]
		#cutname = original_name[5:]
		nameset.add(original_name)

		logval = itemsinline[2]
		if logval == "NA" :
			real_log = 0.0
		else:
			real_log = float(logval)
		namelogDict[original_name] = real_log
		
		pval = itemsinline[5]
		if pval == "NA":
			real_p = 1.0
		else:
			real_p = float(pval)
		namepvalDict[original_name] = real_p
	namelist = list(nameset)
	return namelist, namelogDict, namepvalDict

def STRING_filereader(file):
	IDset = set() #set containing all IDs
	IDDict = {} #ID->name Dictionary
	NameDict = {} #name -> ID Dictionary
	NotesDict = {} #ID-> notes (including protein name) dictionary
	f=open(file,'r')
	print("Opening " + file)
	nameset = set()
	for l in f:
		line = str(l)
		itemsinline = line.split()
		notes = str(itemsinline[11:-5])
		name = itemsinline[9] #common name of the gene
		ID = itemsinline[10] #ID of the gene
		IDset.add(ID) #adds ID to set of all IDs
		IDDict[ID] = name #adds ID->name
		NameDict[name] = ID #adds name->ID
		NotesDict[ID] = notes
	STRING_IDlist = list(IDset) #list of all IDs
	print(NotesDict)
	return STRING_IDlist, IDDict, NameDict, NotesDict


def final_filewriter(IDlist, STRING_IDlist, IDDict, NotesDict, Output_file):
	f=open(Output_file,'w')
	unknown_f = open("UnknownIDs.txt",'w')
	print("Writing to file " + Output_file)

	unknown_set = set()# all genes IDs from RNASeq File that do not have common names

	f.write("Gene ID" + '\t' + "Name" + '\t' + "Notes" + '\n')
	for gene in IDlist: #goes through all the genes from the RNASeq file by ID
		if gene in STRING_IDlist: #if the gene ID exist in the STRING file
			f.write(gene + '\t' + IDDict[gene] + '\t' + NotesDict[gene] + '\n') #write to file
		else:
			f.write(gene + '\t' + "NA" + '\n') #common name not known
			unknown_set.add(gene)
			unknown_f.write(gene + '\t')

	unknown_list = list(unknown_set)	
	unknown_num = (len(unknown_list)/len(IDlist)) * 100
	known_num = len(IDlist) - len(unknown_list)
	print("Number of unknowns: " + str(len(unknown_list)))
	print("Number of known genes: " + str(known_num))
	print("Total number of genes: " + str(len(IDlist)))
	print("Percent Unknown: " + str(unknown_num) + "%")
	return unknown_list


main()

