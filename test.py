import numpy

daily_prices = numpy.array(
    [
        (4,3,3,1),
        (5,4,3,6),
        (6,3,2,7),
        (3,9,7,4),
        (8,4,6,3),
        (8,3,3,9)],
    dtype=[('MSFT','float'),('CSCO','float'),('GOOG','float'),('F','float') ]
    )

matrix_list = []
matrix_list.append(daily_prices)

print matrix_list[0][0]


def addHeader(matrix, taxaList):
	size = matrix.shape[0] + 1
	mat = numpy.zeros((size, size))

	# adds header based on taxaList
	for i in range(0, len(taxaList)):
		mat[0][i+1] = float(taxaList[i])
		mat[i+1][0] = float(taxaList[i])

	#copies matrix to the bottom right corner of larger matrix
	rowCounter = 0
	for x in range(0, size - 1):
		rowCounter += 1
		colCounter = 0
		for y in range(0, size - 1):
			colCounter += 1
			mat[rowCounter][colCounter] = matrix[x][y]
			print colCounter

	print mat
	return mat