"""neighbor joining algorithm"""
import numpy

###### OVERALL LOGIC: #######
# Given a string of doctor's DNA sequence, and a lot of patient's DNA sequence,
# we can compute percentage of sites these sequences differ, denoted by p. Then we can 
# calculate estimated distance between doctor and patients based on d = 3/4 ln(1-4/3*p).
# Then the input matrix for neighbor joining is going to be these distances, where
# a = doctor, b = patient_1, .... 
###########################

###### HOW TO DO THIS?######
# This code basically follows the algorithm given by wiki here:
# 	https://en.wikipedia.org/wiki/Neighbor_joining#The_Q-matrix
# After running throught the alrithm, we output an array, called intermediateLengths, 
# which contains all the lengths we want. From this, we generate a .gv file 
# for graphviz, which outputs the tree.
##############################


# FOR TESTING PURPOSES:
# example matrix given in wiki:
#
#      a   b   c   d   e
# a [[ 0.  5.  9.  9.  8.]
# b  [ 5.  0.  10. 10. 9.]
# c  [ 9.  10. 0.  8.  7.]
# d  [ 9.  10. 8.  0.  3.]
# e  [ 8.  9.  7.  3.  0.]]
#
# From this, create a QMatrix, then join 2 vertices, which creates another matrix, 
# then create another QMatrix based on the new matrix... repeat until all 
# vertices are joined.

testMatrix = numpy.zeros((5, 5))

testMatrix[0, :] = numpy.array([0, 5, 9, 9, 8])
testMatrix[1, :] = numpy.array([5, 0, 10, 10, 9])
testMatrix[2, :] = numpy.array([9, 10, 0, 8, 7])
testMatrix[3, :] = numpy.array([9, 10, 8, 0, 3])
testMatrix[4, :] = numpy.array([8, 9, 7, 3, 0])


def calcRowSum(matrix, row):
	"""calculates the sum of a particular row, but works for columns too since these 
	matrices are symmetric."""

	sizeOfRow = matrix.shape[0]
	rowSum = 0
	colCounter = 0
	while colCounter < sizeOfRow:
		rowSum += matrix[row][colCounter]
		colCounter += 1
	return rowSum


def findMin(matrix):
	"""Returns the coordinate of the lowest entry in a matrix as a pair a number.
	Only searches the top triangular matrix."""
	currentMin = matrix[0][0]
	xCor, yCor = 0, 0
	for x in range(0, matrix.shape[0]):
		for y in range(x + 1, matrix.shape[1]):
			if matrix[x][y] < currentMin:
				currentMin = matrix[x][y]
				xCor, yCor = x, y
	return (xCor, yCor)	


# 
def QMatrix(distMatrix):
	"""Inputs a distance matrix d(i,j) and outputs a QMatrix Q(i,j).
	Use the following formula found on wiki:
	Q(i,j) = (n-2)*d(i,j) - (sum of row i) - (sum of column j)
	example: Q(a, b) = 3*5 - (5+9+9+8) - (5+10+10+9) = -50"""

	assert distMatrix.shape[0] == distMatrix.shape[1], 'matrix not square'
	assert distMatrix.shape[0] > 3, 'matrix size 3 or smaller'

	# first make a nxn zero matrix, and fill it with the correct numbers.
	sizeOfMatrix = distMatrix.shape[0]
	QMatrix = numpy.zeros((sizeOfMatrix, sizeOfMatrix))

	for x in range(0, sizeOfMatrix):
		rowSum = calcRowSum(distMatrix, x)
		for y in range(0, sizeOfMatrix):
			if x == y:
				continue # don't wanna calculate values along the diagonal
			colSum = calcRowSum(distMatrix, y)
			QMatrix[x][y] = (sizeOfMatrix - 2) * distMatrix[x][y] - colSum - rowSum

	return QMatrix


"""a list that stores the values of the 2 "new distance" created when makeNewMatrix
is ran. Need this to output a phylogenetic tree, but storing it this way might
not be the smartest thing to do."""
intermediateLengths = []


def makeNewMatrix(distMatrix, QMatrix):
	"""takes in the original distance matrix and a QMatrix, outputs a distance matrix whose 
	dimension is 1 smaller than the QMatrix. The new node will be on the first row and 
	column. """

	assert distMatrix.shape[0] == QMatrix.shape[0], 'Impossible: distance matrix and QMatrix have different size.'
	assert distMatrix.shape[1] == QMatrix.shape[1], 'Impossible: distance matrix and QMatrix have different size.'

	matrixSize = QMatrix.shape[0]
	minRowCoord, minColCoord = findMin(QMatrix)[0], findMin(QMatrix)[1]
	AtoBDist = distMatrix[minRowCoord][minColCoord] # distance between the closest 2 nodes
	randomConst = 1. / (2 * (matrixSize - 2)) # need . to turn 1 into a float. Fuck wasted 30 minutes here trying catch this bug. Fuck python 2.7

	#calculates the distance from Node A to the newly constructed Node.
	AtoNewNodeDist = 0.5 * AtoBDist + randomConst * (calcRowSum(distMatrix, minRowCoord) - calcRowSum(distMatrix, minColCoord)) 
	BtoNewNodeDist = AtoBDist - AtoNewNodeDist
	intermediateLengths.append(AtoNewNodeDist)
	intermediateLengths.append(BtoNewNodeDist)

	#makes a new matrix whose size is 1 smaller than the QMatrix
	newMatrix = numpy.zeros((matrixSize - 1, matrixSize - 1))

	# first search through the original matrix, and copies unaffected distances
	# to the lower right hand corner of the new matrix. The 2 integers keep tracks
	# of how many rows/columns have been "found" that needs to be copied, so the 
	# numbers can be transfered to the new matrix.
	rowCounter = 0
	for x in range(0, matrixSize):
		if (x == minRowCoord) or (x == minColCoord):
			continue # the new matrix should not include the original 2 nodes.
		rowCounter += 1
		colCounter = 0
		for y in range(0, matrixSize):
			if (y == minRowCoord) or (y == minColCoord):
				continue
			colCounter += 1
			newMatrix[rowCounter][colCounter] = distMatrix[x][y]

	# now search the original distance matrix, and locate the row of ANode (not BNode!)
	# and print that row into a list. From this list, subtract AtoNewNodeDistm from each 
	# entry that was not in ANode's row or BNode's row, so that distances from the new 
	# node and remaining nodes that wasn't changed can be updated. 
	newFirstRow = [row[minRowCoord] for row in distMatrix]
	newMatrixCounter = 1
	for i in range(0, len(newFirstRow)):
		if i == minRowCoord or i == minColCoord:
			continue
		newMatrix[0][newMatrixCounter] = newFirstRow[i] - AtoNewNodeDist
		newMatrix[newMatrixCounter][0] = newFirstRow[i] - AtoNewNodeDist
		newMatrixCounter += 1

	return newMatrix

def calc3x3Matrix(distMatrix):
	"""takes in the final 3x3 matrix, and outputs the final 3 distances.
	   here A is really node D, B is node E for visualization purposes, if you 
	   are looking at the wiki page. """
	AtoBDist = distMatrix[1][2]
	AtoNewNodeDist = 0.5 * AtoBDist + 0.5 * (calcRowSum(distMatrix, 1) - calcRowSum(distMatrix, 2)) 
	BtoNewNodeDist = AtoBDist - AtoNewNodeDist
	lastUnresolvedDist = distMatrix[0][1] - AtoNewNodeDist
	return (lastUnresolvedDist, AtoNewNodeDist, BtoNewNodeDist)

def writeGV(array):
	"""takes in final output of intermediateLengths and writes a .gv file
	   that can be turned into a visualization of the tree"""

	with open('graph.gv', 'a') as the_file:
		the_file.write('digraph finite_state_machine {\n')
		the_file.write('	rankdir=LR;\n')
		the_file.write('	size="8,5"\n')
		the_file.write('	node [shape = circle];\n')
		the_file.write('	Node_1 -> middle_1 [ label = ' + str(array[0]) + '];\n')
		the_file.write('	Node_2 -> middle_1 [ label = ' + str(array[1]) + '];\n')
		the_file.write('	middle_1 -> middle_2 [ label = ' + str(array[2]) + ' ];\n')

		# auto assign every edge/node except the first 3 distance in given array
		numIteration = (len(array) - 3) / 2
		for i in range(0, numIteration): 
			the_file.write('	middle_' + str(i+2) + ' -> Node_' + str(i+3) + ' [ label = ' + str(array[i+3]) + ' ];\n')
			if i == numIteration - 1:
				the_file.write('	middle_' + str(i+2) + ' -> Node_' + str(i+4) + ' [ label = ' + str(array[-1]) + ' ];\n')
			else:
				the_file.write('	middle_' + str(i+2) + ' -> middle_' + str(i+3) + ' [ label = ' + str(array[i+4]) + ' ];\n')
		
		the_file.write('}')


"""MAIN METHOD IS HERE: calling neightbor joining on testMatrix, returning an array
   that contains all distances. Hope this works."""
print testMatrix
while testMatrix.shape[0] > 3:
	matrix1 = QMatrix(testMatrix)
	matrix2 = makeNewMatrix(testMatrix, matrix1)
	print matrix2
	testMatrix = matrix2

last3Dist = calc3x3Matrix(testMatrix)

for i in range(3):
	intermediateLengths.append(last3Dist[i])

print intermediateLengths

writeGV(intermediateLengths)





		 



