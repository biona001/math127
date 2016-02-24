"""neighbor joining algorithm"""
import numpy

###### OVERALL LOGIC: #######
# Given a string of doctor's DNA sequence, and a lot of patient's DNA sequence,
# we can compute percentage of sites these sequences differ, denoted by p. Then we can 
# calculate estimated distance between doctor and patients based on d = 3/4 ln(1-4/3*p).
# Then the input matrix for neighbor joining is going to be these distances, where
# a = doctor, b = patient_1, .... 
###########################

###### HOW TO DO THIS?#######
# This code basically follows the algorithm given by wiki here:
# 	https://en.wikipedia.org/wiki/Neighbor_joining#The_Q-matrix
# After running throught the alrithm, we output an array, called intermediateLengths, 
# which contains all the lengths we want. From this, we generate a .gv file 
# for graphviz, which outputs the tree.
##############################

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
			#if matrix[x][y] - currentMin < 0.00000001 :
			if matrix[x][y] < currentMin:
				currentMin = matrix[x][y]
				xCor, yCor = x, y
	return (xCor, yCor)	

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


def saveMatrix(dist1, dist2, AtoBDist):
	"""makes a 3x3 distance matrix. The purpose of this matrix is to keep track of the distance 
	of the joined nodes. This mini distance matrix will be stored in matrixList, followed by a 
	single row containing the name of the rows."""
	mat = numpy.zeros((3, 3))
	mat[0][1], mat[1][0] = dist1, dist1
	mat[0][2], mat[2][0] = dist2, dist2
	mat[1][2], mat[2][1] = AtoBDist, AtoBDist
	return mat


"""a list of matrices that stores the values of the 2 "new distance" to the new node
created when makeNewMatrix is ran. """
matrixList = []
numSplits = 0


def makeNewMatrix(distMatrix, QMatrix, taxaList):
	"""takes in the original distance matrix, a QMatrix, and a list of taxas, 
	outputs a distance matrix whose dimension is 1 smaller than the QMatrix, and 
	a taxa list where the 2 joined nodes are removed. The newly created node will be on the 
	first row and column of the new matrix. """

	assert distMatrix.shape[0] == QMatrix.shape[0], 'Impossible: distance matrix and QMatrix have different size.'
	assert distMatrix.shape[1] == QMatrix.shape[1], 'Impossible: distance matrix and QMatrix have different size.'

	matrixSize = QMatrix.shape[0]
	minRowCoord, minColCoord = findMin(QMatrix)[0], findMin(QMatrix)[1]
	AtoBDist = distMatrix[minRowCoord][minColCoord] # distance between the closest 2 nodes
	randomConst = 1. / (2 * (matrixSize - 2)) 

	#calculates the distance from Node A to the newly constructed Node.
	AtoNewNodeDist = 0.5 * AtoBDist + randomConst * (calcRowSum(distMatrix, minRowCoord) - calcRowSum(distMatrix, minColCoord)) 
	BtoNewNodeDist = AtoBDist - AtoNewNodeDist

	# After connecting A and B to a new node split_0, put their distances into a 3x3 matrix. To 
	# keep track of their node names, put their names into an array. Then put [3x3matrix, node_name]
	# as a list into the matrixList. I know it's complicated but this is the best idea I have. 
	A3x3Matrix = saveMatrix(AtoNewNodeDist, BtoNewNodeDist, AtoBDist)
	global numSplits #need to declare global to modify a global variable
	newNode = "split_" + str(numSplits)
	numSplits += 1 
	smallMatrixNodeList = [newNode, taxaList[minRowCoord], taxaList[minColCoord]]
	matrixList.append([A3x3Matrix, smallMatrixNodeList])
	
	#deleting 2 taxa and adding the new taxa
	for index in sorted([minRowCoord, minColCoord], reverse=True): 
		del taxaList[index]
	taxaList.insert(0, newNode)

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

	# calculates the first row and colume of the new distance matrix. ie the distances
	# from the new node to all other unaffected node. 
	newMatrixCounter = 1
	for i in range(0, distMatrix.shape[0]):
		if i == minRowCoord or i == minColCoord:
			continue
		newMatrix[0][newMatrixCounter] = 0.5 * (distMatrix[minRowCoord][i] + distMatrix[minColCoord][i] - AtoBDist)
		newMatrix[newMatrixCounter][0] = 0.5 * (distMatrix[minRowCoord][i] + distMatrix[minColCoord][i] - AtoBDist)
		newMatrixCounter += 1

	return newMatrix

def calc3x3Matrix(distMatrix, taxaList):
	""" This is the 3-point-formula taught in lecture. Only exucuted when exactly 3 nodes
	are left to be joined. """
	AtoB = distMatrix[0][1]
	AtoC = distMatrix[0][2]
	BtoC = distMatrix[1][2]

	xDist = 0.5 * (AtoB + AtoC - BtoC)
	yDist = 0.5 * (AtoB + BtoC - AtoC)
	zDist = 0.5 * (AtoC + BtoC - AtoB)
	distList = [xDist, yDist, zDist]

	global numSplits 
	newNode = "split_" + str(numSplits)
	allNodes = [newNode, taxaList[0], taxaList[1], taxaList[2]]
	
	matrixList.append([distList, allNodes])
	return

def writeGV(matrixList):
	with open('final_test3.gv', 'a') as the_file:
		the_file.write('graph G {\n')

		for pair in matrixList: 
			if len(pair[1]) == 4: # if at the last pair
				the_file.write('	' + str(pair[1][0]) + ' -- ' + str(pair[1][1]) + '[ label = ' + str(pair[0][0]) + '];\n')
				the_file.write('	' + str(pair[1][0]) + ' -- ' + str(pair[1][2]) + '[ label = ' + str(pair[0][1]) + '];\n')
				the_file.write('	' + str(pair[1][0]) + ' -- ' + str(pair[1][3]) + '[ label = ' + str(pair[0][2]) + '];\n')

			else:
				the_file.write('	' + str(pair[1][0]) + ' -- ' + str(pair[1][1]) + '[ label = ' + str(pair[0][0][1]) + '];\n')
				the_file.write('	' + str(pair[1][0]) + ' -- ' + str(pair[1][2]) + '[ label = ' + str(pair[0][0][2]) + '];\n')

		the_file.write('}')


"""calling neighbor joining on testMatrix. Writes a .gv file."""
def main():
	#mat = numpy.zeros((5, 5))
	#mat[0, :] = numpy.array([0, 5, 9, 9, 8])
	#mat[1, :] = numpy.array([5, 0, 10, 10, 9])
	#mat[2, :] = numpy.array([9, 10, 0, 8, 7])
	#mat[3, :] = numpy.array([9, 10, 8, 0, 3])
	#mat[4, :] = numpy.array([8, 9, 7, 3, 0])
	#taxaList = ['a', 'b', 'c', 'd', 'e']

	#mat = numpy.zeros((5, 5))
	#mat[0, :] = numpy.array([0.0, 0.18, 0.11, 0.113, 0.215])
	#mat[1, :] = numpy.array([0.189, 0.0, 0.179, 0.192, 0.211])
	#mat[2, :] = numpy.array([0.11, 0.179, 0.0, 0.094, 0.205])
	#mat[3, :] = numpy.array([0.113, 0.192, 0.094, 0.0, 0.210])
	#mat[4, :] = numpy.array([0.215, 0.211, 0.205, 0.214, 0.0])
	#taxaList = ['gorrila', 'orangutan', 'human', 'chimp', 'gibbon']

	mat = numpy.zeros((6, 6))
	mat[0, :] = numpy.array([0, 5, 4, 7, 6, 8])
	mat[0, :] = numpy.array([5, 0, 7, 10, 9, 11])
	mat[0, :] = numpy.array([4, 7, 0, 7, 6, 8])
	mat[0, :] = numpy.array([7, 10, 7, 0, 5, 9])
	mat[0, :] = numpy.array([6, 9, 6, 5, 0, 8])
	mat[0, :] = numpy.array([8, 11, 8, 9, 8, 0])
	taxaList = ['a', 'b', 'c', 'd', 'e', 'f']

	while mat.shape[0] > 3:
		matrix1 = QMatrix(mat)
		matrix2 = makeNewMatrix(mat, matrix1, taxaList)
		mat = matrix2

	calc3x3Matrix(mat, taxaList)

	print matrixList

	writeGV(matrixList)
	print ".gv file written"


main()





		 



