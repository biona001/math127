"""neighbor joining algorithm"""
import numpy

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

testMatrix[0][1] = 5
testMatrix[0][2] = 9
testMatrix[0][3] = 9
testMatrix[0][4] = 8
testMatrix[1][0] = 5
testMatrix[1][2] = 10
testMatrix[1][3] = 10
testMatrix[1][4] = 9
testMatrix[2][0] = 9
testMatrix[2][1] = 10
testMatrix[2][3] = 8
testMatrix[2][4] = 7
testMatrix[3][0] = 9
testMatrix[3][1] = 10
testMatrix[3][2] = 8
testMatrix[3][4] = 3
testMatrix[4][0] = 8
testMatrix[4][1] = 9
testMatrix[4][2] = 7
testMatrix[4][3] = 3




# calculates the sum of a particular row, but works for columns too since these 
# matrix are symmetric
def calcRowSum(matrix, row):
	sizeOfRow = matrix.shape[0]
	rowSum = 0
	colCounter = 0
	while colCounter < sizeOfRow:
		rowSum += matrix[row][colCounter]
		colCounter += 1
	return rowSum


# Inputs a distance matrix d(i,j) and outputs a QMatrix Q(i,j).
# Use the following formula found on wiki:
# Q(i,j) = (n-2)*d(i,j) - (sum of row i) - (sum of column j)
# example: Q(a, b) = 3*5 - (5+9+9+8) - (5+10+10+9) = -50
def QMatrix(inputMatrix):
	assert inputMatrix.shape[0] == inputMatrix.shape[1], 'matrix not square'
	assert inputMatrix.shape[0] > 3, 'matrix size 3 or smaller'

	# first make a nxn zero matrix, and fill it with the correct numbers.
	sizeOfMatrix = inputMatrix.shape[0]
	QMatrix = numpy.zeros((sizeOfMatrix, sizeOfMatrix))

	for x in range(0, sizeOfMatrix):
		rowSum = calcRowSum(inputMatrix, x)
		for y in range(0, sizeOfMatrix):
			if x == y:
				continue # don't wanna calculate values along the diagonal
			colSum = calcRowSum(inputMatrix, y)
			QMatrix[x][y] = (sizeOfMatrix - 2) * inputMatrix[x][y] - colSum - rowSum

	return QMatrix


# Returns the coordinate of the lowest entry in a matrix as a pair a number.
# Only searches the top triangular matrix.
def findMin(matrix):
	currentMin = matrix[0][0]
	xCor, yCor = 0, 0
	for x in range(0, matrix.shape[0]):
		for y in range(x, matrix.shape[1]):
			if (x == y):
				continue

			if matrix[x][y] < currentMin:
				currentMin = matrix[x][y]
				xCor, yCor = x, y
	return (xCor, yCor)	


# takes in the original distance matrix and a QMatrix, outputs a distance matrix whose 
# dimension is 1 smaller than the QMatrix. The new node will be on the first row and 
# column. 
def makeNewMatrix(distMatrix, QMatrix):
	assert distMatrix.shape[0] == QMatrix.shape[0], 'Impossible: distance matrix and QMatrix have different size.'
	assert distMatrix.shape[1] == QMatrix.shape[1], 'Impossible: distance matrix and QMatrix have different size.'

	matrixSize = QMatrix.shape[0]
	minRowCoord, minColCoord = findMin(QMatrix)[0], findMin(QMatrix)[1]
	AtoBDist = distMatrix[minRowCoord][minColCoord] # distance between the closest 2 nodes
	randomConst = 1. / (2 * (matrixSize - 2)) # need . to turn 1 into a float. Fuck wasted 30 minutes here trying catch this bug. Fuck python 2.7

	#calculates the distance from Node A to the newly constructed Node.
	AtoNewNodeDist = 0.5 * AtoBDist + randomConst * (calcRowSum(distMatrix, minRowCoord) - calcRowSum(distMatrix, minColCoord)) 
	BtoNewNodeDist = AtoBDist - AtoNewNodeDist


	#makes a new matrix whose size is 1 smaller than the QMatrix
	newMatrix = numpy.zeros((matrixSize - 1, matrixSize - 1))

	#for x in range(0, matrixSize - 1):
		#newMatrix[0][x] = 


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

	return newMatrix



matrix1 = QMatrix(testMatrix)
print makeNewMatrix(testMatrix, matrix1)





		 



