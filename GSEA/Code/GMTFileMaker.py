'''
This code reformats the GO terms file obtained from biocyc into a GMT format so that it can be used with fgsea
'''



def main():
	InputFile = "PathwayDataFIXED.txt"
	ClassFile = "Class_CommonNames.txt"
	OutputFile = "O157GOSuper.txt"
	PathwayDict = filereader(InputFile)
	ClassGeneDict = ClassReader(ClassFile, PathwayDict)
	OutputMaker(OutputFile, ClassGeneDict)
	return


#Input: Weirdly formatted BioCyc File
#Output: PathwayName -> pathway gene list
def filereader(file):
	PathwayDict = {}
	f=open(file,'r')
	print("Opening " + file)
	for l in f:
		genelist = []
		line = str(l)
		itemsinline = line.split('\t')
		pathway_name = itemsinline[0]
		count = 0
		for item in itemsinline:
			if count >0:
				if len(item) >0: #gets rid of empty spaces
					item = item.replace('"', '') #removes extra quotation makrs
					item = item.replace(" ", '')#removes spaces from genes (NOT PATHS)
					genelist.append(item)
			else:
				count = count + 1 #so that the pathway name is not included
		genelist.remove(genelist[-1]) #removes the "\n"
		PathwayDict[pathway_name] = genelist
	
	return PathwayDict

def ClassReader(ClassFile, PathwayDict):
	ClassPathDict = {} #superpathway -> pathway list
	PathClassDict = {} #path -> list of super pathways
	ClassGeneDict = {} #superpathway -> genes list
	ClassSet = set() #set of all classes
	f=open(ClassFile,'r')
	print("Opening " + ClassFile)
	for l in f:
		line = str(l)
		itemsinline = line.split('\t')
		path = itemsinline[0]
		path = path.replace('"', '') #removes extra quotation makrs
		classlist = itemsinline[1]
		classsplit = classlist.split('//') #split classes
		Edited_ClassList = []
		for c in classsplit:
			c = c.replace('"', '') #removes quotation marks
			Edited_ClassList.append(c)
			ClassSet.add(c)
		PathClassDict[path] = Edited_ClassList#path -> list of superpaths to which it belongs
	for c in ClassSet:
		ClassPathDict[c] = [] #initiates empty list for each class
	for p in PathClassDict: #goes through every path in dictionary
		for superpath in PathClassDict[p]: #goes through every super path/class that path belongs to
			splist = ClassPathDict[superpath] + [p] #adds path to that class's list of paths
			ClassPathDict[superpath] = splist #overwrites the dictionary entry for that class to include the path
	keyerrors = []
	for sp in ClassPathDict:#goes through all super paths
		spgenelist = [] #list of all genes in superpath
		for path in ClassPathDict[sp]: #for all the path paths in a super path
			if path in PathwayDict:
				genelist = PathwayDict[path] #all the genes of that path
			else:
				if path not in keyerrors:
					keyerrors.append(path)
			for gene in genelist:
				if gene not in spgenelist: #if it's not already in the list
					spgenelist.append(gene)
			ClassGeneDict[sp] = spgenelist
	print(len(keyerrors)) #15 paths result in key error		
	return ClassGeneDict


#writes output to simple txt file
def OutputMaker(Output_file, PathwayDict):
	f=open(Output_file,'w')
	print("Writing to file " + Output_file)

	for path in PathwayDict:
		substring = '' #makes a substring to put in output file with all gene names
		for gene in PathwayDict[path]:
			substring = substring + gene + '\t'
			inputstring = path + '\t' + "NA" + '\t' + substring + '\n'
			print(substring)
		f.write(path + '\t' + "NA" + '\t' + substring + '\n')


main()