'''
SENIOR THESIS WORK
Author: Sol Taylor-Brill
Version: 02/03/20

This program takes input from ID-> common name info from NCBI
'''

def main():
	UnmappedData = "DEDATA.txt" #txt file W/O common names
	NameMapFile = "NameMapECS.txt" #ID->common name/notes file
	OutputFile = "MappedDEDATA.txt" #DEDATA + common_name + protein notes

	Known_IDlist, IDDict, NameDict, NotesDict = IDDict_maker(NameMapFile)
	UnmappedIDLIST, IDStringDict = DEReader(UnmappedData)
	CombineOutput(UnmappedIDLIST, Known_IDlist, IDStringDict, IDDict, NotesDict, OutputFile)
	return


'''
This function goes through the name map files and makes ID>Name/Name>ID dictionaries and makes a list of gene
IDs for which the common name is known (ie gene appears in the original source file) 
'''
def IDDict_maker(MapFile):
	IDDict = {} #ID-> common name dictionary
	NameDict = {} #common name -> ID dictionary
	NotesDict = {} # ID-> protein notes dictionary
	Known_IDset = set() #a list of all

	NameMap = open(MapFile,'r')
	print("opening ", MapFile)

	for l in NameMap:
		line = str(l) #line is equal to string
		itemsinline = line.split() #splits items up by tabs OR spaces and makes a list
		ID = itemsinline[0] #ID
		name = itemsinline[1] #common name
		if name != 'NA': #only includes genes for which the file has a known common name
			noteslist = itemsinline[2:] #notes
			count = 0
			notes = []
			for item in noteslist:
				if count == 0:
					notes.append(item[2:-2])
					count = count + 1
				else:
					notes.append(item[1:-2])
			IDDict[ID] = name #adds to ID>name dictionary
			NameDict[name] = ID #adds to name>ID dictionary
			NotesDict[ID] = notes
			Known_IDset.add(ID) #adds gene to set of known genes
	Known_IDlist = list(Known_IDset) #makes the set a list
	return Known_IDlist, IDDict, NameDict, NotesDict


def DEReader(UnmappedData):
	IDStringDict = {} #Gene ID -> rest of the line
	UnmappedIDSet = set()

	Unmapped = open(UnmappedData,'r')
	print("opening ", Unmapped)

	for l in Unmapped:
		line = str(l) #line is equal to string
		itemsinline = line.split() #splits items up by tabs OR spaces and makes a list
		ID = itemsinline[0] #ID
		everything_else = itemsinline[1:] #everything but the ID
		IDStringDict[ID] = everything_else
		UnmappedIDSet.add(ID)

	UnmappedIDLIST = list(UnmappedIDSet)
	return UnmappedIDLIST, IDStringDict

def CombineOutput(UnmappedIDLIST, Known_IDlist, IDStringDict, IDDict, NotesDict, OutputFile):
	f = open(OutputFile, 'w')
	print("Opening", OutputFile)
	f.write("GeneID" + "\t" + "BaseMean" + "\t" + "Log2FC"+ "\t" + "StdErr" + "\t" + "Wald_Stats" + "\t" + "Pval" + "\t" + "P_adj"+ "\n")

	for gene in UnmappedIDLIST: #goes through all the genes in the original DEseq file
		original_items = IDStringDict[gene] #every item in a line
		item0 = original_items[0]
		item1 = original_items[1]
		item2 = original_items[2]
		item3 = original_items[3]
		item4 = original_items[4]
		item5 = original_items[5]
		if gene in Known_IDlist: #if the common name is in the file
			name = IDDict[gene]
			notelist = NotesDict[gene]
			notes = ' '.join(notelist)
			print(notes)
			f.write(gene + '\t' + item0 + '\t' + item1 + '\t' + item2 + '\t' + item3 + '\t' + item4 + '\t' + item5 + '\t' + name + '\t' + notes + "\n")
		else:
			f.write(gene + '\t' + item0 + '\t' + item1 + '\t' + item2 + '\t' + item3 + '\t' + item4 + '\t' + item5 + '\t' + "unknown" + '\t' + "unknown" + "\n")

	print("Finished writing")






main()