#....................................Takes input g1 and gm which are list of binaries ..................................................


def graycodes(g1,gm):   # return the graycode sequence starting at g1 and ending at gm 
	GCodes = [] 
	GCodes.append(g1)
	g=[]
	for i in range(len(g1)):
		g.append(g1[i])

	while(g!=gm):
		for i in range(len(g)):
			if(g[len(g)-i-1] != gm[len(g)-i-1]):
				g[len(g)-i-1] = gm[len(g)-i-1]
				break
		g_copy = [] 
		for i in range(len(g)):
			g_copy.append(g[i])
	
		#print(g)
		GCodes.append(g_copy) 
	return(GCodes)	


#.................................................DEBUGGED............................................................................
