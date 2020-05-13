'''
SENIOR THESIS WORK
Author: Sol Taylor-Brill
Version: 10/12/19

This program takes an input file from STRING containing a variety of info about proteins in E. coli O157 
(https://string-db.org/cgi/download.pl?sessionId=Z6T87Go63l8T&species_text=Escherichia+coli+O157%3AH7+str.+EDL933)
parses through it and outputs a file containing just ID and common name for ease of use
'''

def main():
	file_name = "STRING_INFO.txt"
	EXCEL_file = "CBDCopy.txt"
	Output_file = "NameMap.txt"

	IDlist, IDlogDict, IDpvalDict = EXCEL_filereader(EXCEL_file)
	STRING_IDlist, IDDict, NameDict = STRING_filereader(file_name)
	final_filewriter(IDlist, STRING_IDlist, IDDict, Output_file)
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
		cutname = original_name[5:]
		nameset.add(cutname)

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

def STRING_filereader(file):
	IDset = set() #set containing all IDs
	IDDict = {} #ID->name Dictionary
	NameDict = {} #name -> ID Dictionary
	f=open(file,'r')
	print("Opening " + file)
	nameset = set()
	for l in f:
		line = str(l)
		itemsinline = line.split()
		description = itemsinline[3:] #long description of the gene including the locus ID
		name = itemsinline[1] #common name of the gene
		if 'locus_tag' in description: #if 'locus_tag' exists in description
			index = description.index('locus_tag')#index of 'locus_tag' which comes before the ID
			Z = index + 1 #index of ID
			ID = description[Z] #ID of the gene
			IDset.add(ID) #adds ID to set of all IDs
			IDDict[ID] = name #adds ID->name
			NameDict[name] = ID #adds name->ID
	STRING_IDlist = list(IDset)
	return STRING_IDlist, IDDict, NameDict


def final_filewriter(IDlist, STRING_IDlist, IDDict, Output_file):
	f=open(Output_file,'w')
	unknown_f = open("UnknownIDs.txt",'w')
	print("Writing to file " + Output_file)

	unknown_set = set()# all genes IDs from RNASeq File that do not have common names

	for gene in IDlist: #goes through all the genes from the RNASeq file
		if gene in STRING_IDlist: #if the gene ID exist in the STRING file
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

