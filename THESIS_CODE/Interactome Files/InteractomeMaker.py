''''
SENIOR THESIS WORK
This program is used to make an E. coli K12 interactome from data from two studies.
Date: 04/05/20
Author: Sol Taylor-Brill
'''

from graphspace_python.api.client import GraphSpace
from graphspace_python.graphs.classes.gsgraph import GSGraph
graphspace = GraphSpace("soltb@reed.edu", "solTB") #Starting GraphSpace session

def main():
	NIHMf = "NIHM.txt"
	pbiof = "pbio.txt"
	RNAf = "RNASeq.txt"
	NIHM_nodeset, NIHM_edgeset, NIHM_JWDict = NIHM_reader(NIHMf) #reads NIHM file to create JW Dict, edgeset and nodeset
	pbio_nodeset, pbio_edgeset, pbio_JWDict = pbio_reader(pbiof) #reads pbio file to create JW Dict, edgeset, and nodeset
	NIHM_JWDict.update(pbio_JWDict) #combines two JW Dicts. NIHM Dict is now master dict

	NameExpressionDict = RNASeq_reader(RNAf) #reads gene expression file

	ALL_nodeset, ALL_edgeset = file_compare(NIHM_nodeset, NIHM_edgeset, pbio_nodeset, pbio_edgeset)

	candidates, neighbornodes, all_sub_nodes, subinteractome_edges = subinteractome(ALL_nodeset, ALL_edgeset)

	BigInteractomeNodes, BigInteractomeEdges = Big_SubInteractome(ALL_nodeset, ALL_edgeset, NameExpressionDict)
	#nodes_w_expression_list = subGraph_Maker(candidates, neighbornodes, all_sub_nodes, subinteractome_edges, NameExpressionDict, NIHM_JWDict)
	#fgsea_printer(nodes_w_expression_list)

	Graph_Maker(BigInteractomeNodes, BigInteractomeEdges, NameExpressionDict)
	return


def NIHM_reader(file):
	f=open(file,'r')
	print("Opening " + file)
	nodeset = set() #set of all the nodes in the NIHM file (Common names)
	edgeset = set() #set of all the edges in the NIHM file ((common_name1, common_name2))
	JWDict = {} #common name -> JWID dictionary
	
	for l in f:
		line = str(l)
		itemsinline = line.split()
		splitnames = itemsinline[0].split('__')
		commonA = splitnames[0] #common name for gene A
		nodeset.add(commonA)
		if len(splitnames) == 3:
			JWA = splitnames[2] #JW ID for gene A
			JWDict[commonA] = JWA
		splitnamesB = itemsinline[1].split('__') 
		commonB = splitnamesB[0] #common name for gene B
		nodeset.add(commonB)
		if len(splitnamesB) == 3:
			JWB = splitnamesB[2] #JW ID for gene B
			JWDict[commonB] = JWB
		edgeset.add((commonA,commonB))

	print("Number of Nodes in NIHM file:", len(nodeset))
	print("Number of Edges in NIHM file:", len(edgeset), '\n')

	return nodeset, edgeset, JWDict

def pbio_reader(file):
	f=open(file,'r')
	print("Opening " + file)
	nodeset = set()
	edgeset = set()
	JWDict = {} #common name -> JWID dictionary
	for l in f:
		line = str(l)
		itemsinline = line.split()
		if itemsinline[0] == 'Y':
			commonA = itemsinline[1]
			JWA = itemsinline[3]
			commonB = itemsinline[4]
			JWB = itemsinline[6]
		else:
			commonA = itemsinline[0]
			JWA = itemsinline[2]
			commonB = itemsinline[3]
			JWB = itemsinline[5]
		JWDict[commonA] = JWA
		JWDict[commonB] = JWB
		nodeset.add(commonA)
		nodeset.add(commonB)
		edgeset.add((commonA, commonB))

	print("Number of Nodes in pbio file", len(nodeset))
	print("Number of Edges in pbio file", len(edgeset), '\n')

	return nodeset, edgeset, JWDict

def file_compare(NIHM_nodeset, NIHM_edgeset, pbio_nodeset, pbio_edgeset):
	NIHM_nodes_only = 0
	NIHM_edges_only = 0
	nodes_inboth = 0
	edge_inboth = 0

	for node in NIHM_nodeset:
		if node in pbio_nodeset:
			nodes_inboth = nodes_inboth + 1 #adds 1 to the number of nodes that are in both files
		else:
			NIHM_nodes_only = NIHM_nodes_only + 1 #if it's not in the other file adds 1 to NIHM only sum

	for edge in NIHM_edgeset:
		if edge in pbio_edgeset:
			edge_inboth = edge_inboth + 1
		else:
			NIHM_edges_only = NIHM_edges_only + 1

	print("# nodes in both files:", nodes_inboth)
	print("# nodes in NIHM file only:", NIHM_nodes_only)
	pbio_nodes_only = len(pbio_nodeset) - nodes_inboth
	print("# nodes in pbio file only:", pbio_nodes_only, '\n')

	print("# edges in both files:", edge_inboth)
	print("# edges in NIHM file only:", NIHM_edges_only)
	pbio_edges_only = len(pbio_edgeset) - edge_inboth
	print("# edges in pbio file only:", pbio_edges_only, '\n')

	NIHM_nodeset.update(pbio_nodeset)
	ALL_nodeset = NIHM_nodeset
	print("Total # Nodes:", len(ALL_nodeset))

	NIHM_edgeset.update(pbio_edgeset)
	ALL_edgeset = NIHM_edgeset
	print("Total # Edges:", len(ALL_edgeset))

	return ALL_nodeset, ALL_edgeset

def RNASeq_reader(file):
	f=open(file,'r')
	print("Opening " + file)
	NameExpressionDict = {} #dictionary of gene expression for all genes with common names

	for l in f:
		line = str(l)
		itemsinline = line.split()
		if itemsinline[7] != itemsinline[0]: #basically if it has a common name
			NameExpressionDict[itemsinline[7]] = itemsinline[2]

	return NameExpressionDict

def subinteractome(ALL_nodeset, ALL_edgeset):
	candidates = ['sdhC', 'yecI', 'degP', 'fadL', 'fadR']
	subinteractome_edges = set()
	neighbornodes = set() #immediate neighbors
	newneighbornodes = set() #neighbors + neighbors of neighbors
	newsubinteractomeedges = set() #neighbors of neighbors
	for candidate in candidates:
		for edge in ALL_edgeset:
			if candidate in edge:
				subinteractome_edges.add(edge)
				neighbornodes.add(edge[0])
				neighbornodes.add(edge[1])

	for edge in ALL_edgeset:
		for node in neighbornodes:
			if node in edge:
				newsubinteractomeedges.add(edge)
				newneighbornodes.add(edge[0])
				newneighbornodes.add(edge[1])

	newneighbornodes.update(neighbornodes)
	return candidates, neighbornodes, newneighbornodes, newsubinteractomeedges

def subGraph_Maker(candidates, neighbornodes, all_sub_nodes, subinteractome_edges, NameExpressionDict, NIHM_JWDict):
	G = GSGraph()
	G.set_name('Ecoli_K12_SubInteractome')
	G.set_tags(['THESIS'])

	Nodes_w_expression = 0 #number of nodes in subgraph with expression values
	nodes_w_expression_list = []

	NameExpressionDict['yecI'] = '-5.380'

	for node in all_sub_nodes:
		if node in NameExpressionDict:
			link = 'https://www.uniprot.org/uniprot/?query=' + NIHM_JWDict[node]
			popupl = node + '\t' + NameExpressionDict[node] + '\n' + '<a href=' +'"' + link + '">' + "Uniprot" + "</a>"
			Nodes_w_expression = Nodes_w_expression + 1
			nodes_w_expression_list.append((node, float(NameExpressionDict[node])))
			if float(NameExpressionDict[node]) >= 0.5:
				ncolor = 'red'
			elif float(NameExpressionDict[node]) <= -0.5:
				ncolor = 'blue'
			else:
				ncolor = 'purple'
		else:
			ncolor = 'yellow'
			link = 'https://www.uniprot.org/uniprot/?query=' + NIHM_JWDict[node]
			popupl = node + '\n' + '<a href=' +'"' + link + '">' + "Uniprot" + "</a>"
		G.add_node(node, label=node, popup = popupl)
		
		if node in candidates:
			nshape = 'star'
			nsize = 90
		elif node in neighbornodes:
			nshape = 'rectangle'
			nsize = 50
		else:
			nshape = 'ellipse'
			nsize = 40
			
		G.add_node_style(node, color= ncolor, shape= nshape, height = nsize, width = nsize)

	for edge in subinteractome_edges:
		G.add_edge(edge[0], edge[1])

	print(nodes_w_expression_list, Nodes_w_expression)

	#graphspace.update_graph(G)
	print("Graph Updated")

	return nodes_w_expression_list

def fgsea_printer(nodelist):
	nodelist.sort(reverse = True, key = lambda x: x[1])

	f=open("Subinteractome_fgsea.txt",'w')
	print("Opening file")

	for node in nodelist:
		string = node[0] + '\t' + str(node[1]) + '\n'
		f.write(string)

	print("Finished Printing!")

	return

def Big_SubInteractome(ALL_nodeset, ALL_edgeset, NameExpressionDict):
	print("Starting to construct big sub interactome")
	BigInteractomeNodes = set()
	BigInteractomeEdges = set()
	for node in ALL_nodeset:
		if node in NameExpressionDict:
			BigInteractomeNodes.add(node)
	for edge in ALL_edgeset:
		if edge[0] in BigInteractomeNodes and edge[1] in BigInteractomeNodes:
			BigInteractomeEdges.add(edge)
	print("Number of Nodes:", len(BigInteractomeNodes)," Number of edges:", len(BigInteractomeEdges))
	return BigInteractomeNodes, BigInteractomeEdges

def Graph_Maker(BigInteractomeNodes,BigInteractomeEdges, NameExpressionDict):
	G = GSGraph()
	G.set_name('Ecoli_K12_w/_Expression')
	G.set_tags(['THESIS'])

	NameExpressionDict['yecI'] = '-5.380'

	for node in BigInteractomeNodes:
		popupl = node + '\t' + NameExpressionDict[node]
		if NameExpressionDict[node] != 'NA':
			if float(NameExpressionDict[node]) >= 0.5:
				ncolor = 'red'
			elif float(NameExpressionDict[node]) <= -0.5:
				ncolor = 'blue'
			else:
				ncolor = 'purple'
		elif NameExpressionDict[node] == 'NA':
			ncolor = 'yellow'
		G.add_node(node, label=node, popup = popupl)
		G.add_node_style(node, color= ncolor)

	for e in BigInteractomeEdges:
		G.add_edge(e[0],e[1])


	graphspace.post_graph(G)
	print("Graph Updated")
	return


main()