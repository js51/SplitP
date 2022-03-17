import itertools

class SquangleEvaluator:
	
	def __init__(self):
		import time
		self.time = time
		self.toPatternIndex = ['A', 'C', 'G', 'T']
		self.toIntIndex = {'A':0, 'C':1, 'G':2, 'T':3}
		self.toIntDict = {}
		self.toPatternDict = {}

	def transformed_prob_dist(self, prob_dist):
		state_space = ('A', 'C', 'G', 'T')
		banned = {('C','C'), ('G','G'), ('A','T')} | {(x, 'A') for x in state_space} | {('T', x) for x in state_space}
		q_dist = {}
		for pat in itertools.product(state_space, repeat=4): # quartets
			pattern = ''.join(pat)
			signed_sum = 0
			for table_pattern, value in prob_dist.items():
				product = 1
				for t in zip(pattern, table_pattern):
					if t not in banned:
						product *= -1
				signed_sum += product*value
			q_dist[pattern] = signed_sum
		return q_dist

	def get_polynomials(self):
		polynomials = []
		import os
		dir_path = os.path.dirname(os.path.realpath(__file__)) + '/polynomials/'
		for i in [1,2,3]:
			f = open(dir_path + 'sqg%s.txt' % str(i),'r')
			lines = f.readlines()
			f.close()
			polynomials.append([[int(number) for number in line.replace('\n','').split()] for line in lines])
		return polynomials[0], polynomials[1], polynomials[2]
	
	def read_polynomial(self, filename):
		f = open(filename,'r')
		lines = f.readlines()
		f.close()
		return [[int(number) for number in line.replace('\n','').split()] for line in lines]
	
	def to_int(self, pattern):
		if pattern in self.toIntDict:
			return self.toIntDict[pattern]
		else:
			index = self.toIntIndex
			integer = sum(4**(3-i) * index[pattern[i]] for i in range(4))
			self.toIntDict[pattern] = integer
			return integer
		
	def to_pattern(self, integer):
		try:
			return self.toPatternDict[integer]
		except KeyError:
			originalInteger = integer
			index = self.toPatternIndex
			pattern = ""
			for i in range(3,-1,-1):
				# print("Integer is", integer, "Dividing by", 4**i)
				temp = int(integer/4**(i))
				pattern += index[temp]
				integer -= temp * 4**(i)
			self.toPatternDict[originalInteger] = pattern
			return pattern
	
	def get_flattening(self, DTable):
		from tree_implementation import tree
		T = tree(-1, 4)
		return T.transformedFlattening(T.flattening("01|23", DTable))
	
	def evaluate_polynomial(self, poly, data):
		finalSum = 0
		for row in poly:
			product = row[0]
			for i in range(1,len(row)):
				pattern = self.to_pattern(int(row[i]) - 1)
				if pattern in data:
					product *= data[pattern]
				else: # Zero, not copied over earlier
					product = 0 # Don't bother multiplying by zero
					break # Don't bother continuing to multiply
			finalSum += product       
		return finalSum
	
	# Unit Tests
	def test_pattern_conversion(self):
		return all(self.toInt(self.to_pattern(i)) == i for i in range(256))