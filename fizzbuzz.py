def fizzbuzz(int):
	for int in range(0, int + 1):
		if (int % 15 == 0):
			print "fizzbuzz"
		elif (int % 5 == 0):
			print "fizz"
		elif (int % 3 == 0):
			print "buzz"
		else:
			print int

fizzbuzz(100)