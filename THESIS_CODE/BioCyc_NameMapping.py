'''
SENIOR THESIS WORK
Author: Sol Taylor-Brill
Version: 10/12/19

This program takes an input file from BioCyc containing a variety of info about all of the genes in E. coli O157
(https://biocyc.org/group?id=biocyc14-41160-3779808068)
parses through it and outputs a file containing just ID and common name for ease of use
'''

def main():
	RNAOutput = "RNA_NameList.txt"
	BioCyc_file = "CommonNames.txt"
	EXCEL_file = "CBDCopy.txt"
	Output_file = "BioCyc_NameMap.txt"

	IDlist, IDlogDict, IDpvalDict = EXCEL_filereader(EXCEL_file, RNAOutput)
	BioCyc_IDlist, IDDict, NameDict = BioCyc_filereader(BioCyc_file)
	final_filewriter(IDlist, BioCyc_IDlist, IDDict, Output_file)
	return


#Reads text file version of excel sheet and makes list of geneIDs (namelist) and dictionaries
# name-> logfold change (namelogDict) and name->p-value (namepvalDict)
#This function is from "CBDEdit.py"
#In this function name is not the common name!! (name = ID here)
def EXCEL_filereader(file, RNAOutput):
	f=open(file,'r')
	print("Opening " + file)
	nameset = set()
	namelogDict = {}
	namepvalDict = {}

	RNA_f = open(RNAOutput, 'w')
	print("Writing to " + RNAOutput)
	for l in f:
		line = str(l)
		itemsinline = line.split()
		original_name = itemsinline[0]
		cutname = original_name[5:]
		nameset.add(cutname)
		RNA_f.write(cutname + '\t')

		logval = itemsinline[2]
		if logval == "NA":
			real_log = 0.0
		else:
			real_log = float(logval)
		namelogDict[cutname] = real_log
		
		pval = itemsinline[5]
		if pval == "NA":
			real_p = 1.0
		else:
			real_p = float(pval)
		namepvalDict[cutname] = real_p
	namelist = list(nameset)
	return namelist, namelogDict, namepvalDict


def BioCyc_filereader(file):
	IDset = set() #set containing all IDs
	IDDict = {} #ID->name Dictionary
	NameDict = {} #name -> ID Dictionary
	nameset = set()

	f=open(file,'r')
	print("Opening " + file)

	for l in f: #goes through every gene
		line = str(l) #line is equal to string
		itemsinline = line.split() #splits items up by tabs OR spaces and makes a list
		untrimmed_name = itemsinline[0] #common name of the gene
		name = untrimmed_name.replace('"','')
		untrimmed_ID = itemsinline[-1] #takes the last item in the list 
		ID = untrimmed_ID.replace('"','')
		IDset.add(ID) #adds ID to set of all IDs
		nameset.add(name)
		IDDict[ID] = name #adds ID->name
		NameDict[name] = ID #adds name->ID
	BioCyc_IDlist = list(IDset)
	return BioCyc_IDlist, IDDict, NameDict

def final_filewriter(IDlist, BioCyc_IDlist, IDDict, Output_file):
	f=open(Output_file,'w')
	unknown_f = open("BioCyc_UnknownIDs.txt",'w')
	print("Writing to file " + Output_file)

	unknown_set = set()# all genes IDs from RNASeq File that do not have common names

	for gene in IDlist: #goes through all the genes from the RNASeq file
		if gene in BioCyc_IDlist: #if the gene ID exist in the STRING file
			f.write(gene + '\t' + IDDict[gene] + '\n') #write to file
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