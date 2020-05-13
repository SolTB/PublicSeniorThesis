'''
This file is just to do some preliminary data analysis and clean up the gene IDs.
Ideally I'd like to figure out how to do these analyses in R instead.
'''

def main():
	file = "CBDCopy.txt"
	namelist, namelogDict, namepvalDict = filereader(file)
	basic_analysis(namelist, namelogDict, namepvalDict)
	return

#Reads text file version of excel sheet and makes list of geneIDs (namelist) and dictionaries
# name-> logfold change (namelogDict) and name->p-value (namepvalDict)
def filereader(file):
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

def basic_analysis(namelist, namelogDict, namepvalDict):
	UnderGeneList = []
	OverGeneList = []
	AllGeneList =[]
	for name in namelist:
		pval = namepvalDict[name]
		logval = namelogDict[name]
		if logval >5 and pval<=0.05:
			OverGeneList.append(name)
		elif logval<5 and pval<=0.05:
			UnderGeneList.append(name)
	AllGeneList = UnderGeneList + OverGeneList
	print("#Under Expressed Genes (p <=0.05): " + str(len(UnderGeneList)))
	print("#Over Expressed Genes (p <=0.05): " + str(len(OverGeneList)))
	print("#Differentially Expressed Genes (p <=0.05): " + str(len(AllGeneList)))
	return AllGeneList, UnderGeneList, OverGeneList



main()


