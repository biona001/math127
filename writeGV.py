def writeGV(intermediateLengths, taxaOrdered):
	"""takes in final output of intermediateLengths and writes a .gv file
	   that can be turned into a visualization of the tree"""
	k = (len(intermediateLengths) - 1) / 2
	someList = [] # a list that looks like [0, 1, 2, ...] which counts the number of splits. aka they are intermediate nodes in a tree
	for i in range(k):
		someList.append(i + 1)

	with open('graph_test7.gv', 'a') as the_file:
		the_file.write('graph G {\n')
		#the_file.write('	rankdir=LR;\n')
		#the_file.write('	size="8,5"\n')
		#the_file.write('	node [shape = circle];\n')
		the_file.write('	' + str(taxaOrdered[0]) + ' -- split_' + str(someList[0]) + ' [ label = ' + str(intermediateLengths[0]) + '];\n')
		the_file.write('	' + str(taxaOrdered[1]) + ' -- split_' + str(someList[0]) + ' [ label = ' + str(intermediateLengths[1]) + '];\n')
		the_file.write('	split_' + str(someList[0]) + ' -- ' + 'split_' + str(someList[1]) + ' [ label = ' + str(intermediateLengths[2]) + '];\n')

		# auto assign every edge/node except the first 3 distance in given array
		numIteration = (len(intermediateLengths) - 3) / 2
		for i in range(0, numIteration): 
			if i == numIteration - 1:
				the_file.write('	split_' + str(someList[-1]) + ' -- ' + str(taxaOrdered[-2]) + ' [ label = ' + str(intermediateLengths[-2]) + ' ];\n')
				the_file.write('	split_' + str(someList[-1]) + ' -- ' + str(taxaOrdered[-1]) + ' [ label = ' + str(intermediateLengths[-1]) + ' ];\n')			
			else:
				the_file.write('	split_' + str(someList[i+1]) + ' -- ' + str(taxaOrdered[i+3]) + ' [ label = ' + str(intermediateLengths[i+3]) + ' ];\n')
				the_file.write('	split_' + str(someList[i+1]) + ' -- split_' + str(someList[i+2]) + ' [ label = ' + str(intermediateLengths[i+4]) + ' ];\n')
		
		the_file.write('}')