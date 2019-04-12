import sys
import csv
import scipy.stats as sp
import itertools as it
import random as ran

file_name = sys.argv[1] 
range_value1 = int(sys.argv[2])
sequenced_number = int(sys.argv[3])

csv.register_dialect('tab_delim', delimiter="\t", quoting=csv.QUOTE_NONE)

def read_lines(csv_reader, row_list):
	for row_number, row in enumerate(csv_reader):
		if row_number in row_list:
			yield row_number, row

with open(file_name, 'r') as File:
	reader = csv.reader(File, dialect='tab_delim')
	r = list(range(0, range_value1))

	# Generate all pairwise combinations of each window and group into tuples,
	# nested lists of individuals into large single list
	replicate_list = []
	for row_number, row in read_lines(reader, r):
		row_tuples = list(it.combinations(row, 2))
		replicate_list.append(row_tuples)

	# Groups window comparisons(tuples) into nested list, such that nested lists represent
	# each pairwise comparison made up of all individuals
	window_list = list(map(list, zip(*replicate_list)))
	
	#randomly selelcted pairwise windows for an X-A and A-A incompatibility
	window_inviableXA = window_list[25000]
	counted_window_inviable_XA = list(enumerate(window_inviableXA))
	window_inviableAA = window_list[79800]
	counted_window_inviable_AA = list(enumerate(window_inviableAA))
	
	inviable_list_XA = []
	inviable_list_AA = []

	for window in counted_window_inviable_XA:
		if window[1]==('0', '2'):
			prob = ran.randint(0, 3)
			if prob ==1:
				inviable_list_XA.append(int(window[0]))
			else:
				pass
		else:
			pass
	
	for window in counted_window_inviable_AA:
		if window[1]==('0','2'):
			prob = ran.randint(0, 3)
			if prob==1:
				inviable_list_AA.append(int(window[0]))
			else:
				pass
		else:
			pass

	combined_inviables = sorted(inviable_list_XA + list(set(inviable_list_AA)-set(inviable_list_XA)))
	
	#must iterate in reverse to avoid index errors	
	for value in reversed(combined_inviables):
		del replicate_list[value]

	replicate_list = replicate_list[0:sequenced_number]
	window_list = list(map(list, zip(*replicate_list)))

	fishergroup1 = []
	fishergroup2 = []
	fishergroup3 = []
	fishergroup4 = []

	for window in window_list[79800]:
		if window==('0', '2'):
			fishergroup1.append(1)
		elif window==('1', '2') or window==('2', '2'):
			fishergroup2.append(1)
		elif window==('0', '1') or window==('0', '0'):
			fishergroup3.append(1)
		elif window==('1', '1') or window==('1', '0'): 
			fishergroup4.append(1)
		elif window==('2', '0') or window==('2', '1'):
			fishergroup4.append(1)
		else:
			raise Exception('There is something wrong with your window values')

	Fisher1 = len(fishergroup1)
	Fisher2 = len(fishergroup2)
	Fisher3 = len(fishergroup3)
	Fisher4 = len(fishergroup4)

	odds, pvalue = sp.fisher_exact([[Fisher1, Fisher2], [Fisher3, Fisher4]])
	print(pvalue)

	File.close()


