#%% Imports
import splitp
import matplotlib.pyplot as plt
import collections
import numpy as np
import pandas as pd

#%% Setup and settings
tree = splitp.NXTree("(A,B);", numStates=4) # Cherry
se = splitp.SquangleEvaluator()

file_location = "input/infile_mosquito.fa"
alignment = splitp.read_alignment_from_file(file_location, check = True)
taxa = [t for t in alignment]

windowSize = 10000
slideSize = 1000
unusedWindows = 0
sequenceLength = len(alignment[taxa[0]])
splits = ["01|23", "02|13", "03|12"]
colors = ['r', 'b', 'limegreen']

poly1, poly2, poly3 = se.get_polynomials()

#Fig1 = plt.figure("Flats")
#Fig2 = plt.figure("Subflats")
#Fig1 = plt.figure("Squangles")

numWindows = int((sequenceLength-windowSize)/slideSize)
linePos = 0

file = "output/sliding_window_mosquitos_test.out"

#%% Main loop
start = 0
try:
	with open(file, "r") as fle:
		for line in fle:
			try :
				start = int(line.split('#')[0].split('-')[1])
			except:
				pass
except FileNotFoundError:
	start = 0

f = open(file, "a")


for i in range(start, sequenceLength - 2, slideSize): # For each window

	print("Window {} of {}".format(int(i/slideSize), numWindows), end='')
	print("(" + str(round((int(i/slideSize)*100)/numWindows,1)) + "% Complete)")

	resultsFlat, resultsSubflat = [], []
	window = dict() # Keeps insertion order in python 3.5+

	for sequence in alignment:
		if (i + windowSize) < sequenceLength:
			window[sequence] = alignment[sequence][i : i + windowSize]
		else:
			window[sequence] = alignment[sequence][i : sequenceLength]
			i = sequenceLength

	counts, newSequenceLength = splitp.get_pattern_counts(window, True)
    print("New sequence length:", newSequenceLength)

	if newSequenceLength > 0.05*windowSize:

		# Pattern probs
		probs = pd.DataFrame(splitp.pattern_counts_to_probs(counts, newSequenceLength))
		probs = probs.set_index(probs[0]).astype({1: 'float64'})

		# Flats and subflats`
		for sp in splits:
			F = tree.sparse_flattening(sp, probs)
			SF = tree.signed_sum_subflattening(sp, probs)
			resultsFlat.append(tree.split_score(F))
			resultsSubflat.append(tree.split_score(SF))

		#plt.figure("Flats")
		flatBestFit = np.argmin(resultsFlat)
		#plt.axvline(linePos, color=colors[flatBestFit])
		#plt.figure("Subflats")
		subflatBestFit = np.argmin(resultsSubflat)
		#plt.axvline(linePos, color=colors[subflatBestFit])
		###

		# Squangles
		flat = tree.flattening("01|23", probs)
		flat = tree.transformed_flattening(flat)
		data = {}

		for r in flat:
			for c in flat:
				if flat.loc[r, c] != 0:
					data[r + c] = float(flat.loc[r, c])
		q1 = se.evaluate_polynomial(poly1, data)
		q2 = se.evaluate_polynomial(poly2, data)
		q3 = se.evaluate_polynomial(poly3, data)

		RSS = [0,0,0]
		for s in splits:
			if s == "01|23":
				if q2 <= q3:
					RSS[0] = (1/2)*(q1**2)
				else:
					RSS[0] = q2**2 + q3**2
			elif s == "02|13":
				if q3 <= q1:
					RSS[1] = (1/2)*(q2**2)
				else:
					RSS[1] =  q1**2 + q3**2
			elif s == "03|12":
				if q1 <= q2:
					RSS[2] = (1/2)*(q3**2)
				else:
					RSS[2] = q1**2 + q2**2

		#plt.figure("Squangles")
		squangleBestFit = np.argmin(RSS)
		#plt.axvline(linePos, color=colors[squangleBestFit])
		###

		linePos += 1

		# Write results to file
		f.write('\n' + str(i) + "-" + str(i+slideSize) + '#(' + str(splits[flatBestFit]) + ',' + str(splits[subflatBestFit]) + ',' + str(splits[squangleBestFit]) + ")#" + str(resultsFlat) + "#" + str(resultsSubflat) + "#" + str(RSS))

	else:
		print("Not using window of size ", newSequenceLength)
		unusedWindows += 1

print("Unused Windows: ", unusedWindows)