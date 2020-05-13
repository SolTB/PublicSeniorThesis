
def main():
	Input = "EditedO157SuperInputText.txt"
	Output = "FINALGMTSP.txt"
	PathwayDict = filereader(Input)
	OutputMaker(Output,PathwayDict)
	return

def filereader(file):
	PathwayDict = {}
	f=open(file,'r')
	print("Opening " + file)
	for l in f:
		line = str(l)
		itemsinline = line.split('\t')
		pathway_name = itemsinline[0]
		genelist = itemsinline[1:]
		newgenelist = [] #list without empty spaces
		for item in genelist:
			if item != '' and item != '\n':
				newgenelist.append(item)
		PathwayDict[pathway_name] = newgenelist
	return PathwayDict


def OutputMaker(Output_file, PathwayDict):
	f=open(Output_file,'w')
	print("Writing to file " + Output_file)
	for path in PathwayDict:
		substring = '' #makes a substring to put in output file with all gene names
		for gene in PathwayDict[path]:
			substring = substring + gene + '\t'
		f.write(path + '\t' + substring+ '\n')

main()
