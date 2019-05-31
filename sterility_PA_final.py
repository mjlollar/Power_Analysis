import sys
import csv
import scipy.stats as sp
import itertools as it
import random as ran


csv.register_dialect('tab_delim', delimiter="\t", quoting=csv.QUOTE_NONE)

file_name = sys.argv[1]
steriles_needed = int(sys.argv[2])
fertiles_needed = int(sys.argv[3])

# Function to enumerate replicates in ancestry file

def read_lines(csv_reader, row_list):
	for row_number, row in enumerate(csv_reader):
		if row_number in row_list:
			yield row_number, row

with open(file_name, 'r') as File:
	replicate_list = []
	reader = csv.reader(File, dialect='tab_delim')
	window1 = list(range(0, 2000))

	# Generate all pairwise combinations of window pairs

	for row_number, row in read_lines(reader, window1):
		row_tuples = list(it.combinations(row, 2))
		replicate_list.append(row_tuples)

	# Randomly picked X-A and A-A pair that are defined as our sterility loci
	
	window_list = list(map(list, zip(*replicate_list)))
	window_sterileXA = window_list[25000]
	counted_window_sterile_XA = list(enumerate(window_sterileXA))
	window_sterileAA = window_list[79800]
	counted_window_sterile_AA = list(enumerate(window_sterileAA))
	
	sterile_count_AA = []
	fertile_count_AA = []
	fertile_focal_count_AA = []

	# Count sterile A-A numbers in sequenced and fertile focal non-penatrent

	for window in counted_window_sterile_AA:
		if window[1] == ('0', '2'):
			prob = ran.randint(0, 3)
			if prob == 1:
				sterile_count_AA.append(window[0])
			else:
				fertile_focal_count_AA.append(window[0])
		else:
			fertile_count_AA.append(window[0])	
	
	sterile_count_XA = []
	fertile_count_XA = []

	# Count sterile X-A numbers in sequenced

	for window in counted_window_sterile_XA:
		if window[1] == ('0', '2'):
			prob = ran.randint(0, 3)
			if prob ==1:
				sterile_count_XA.append(window[0])
			else:
				fertile_count_XA.append(window[0])
		else:
			fertile_count_XA.append(window[0])

	
	
	int_sterile_XA_set = set(sterile_count_XA)
	int_sterile_AA_set = set(sterile_count_AA)
	final_sterile_XA_list = list(int_sterile_XA_set - int_sterile_AA_set)
	final_sterile_XA_set = set(final_sterile_XA_list)

	combined_sterile_sets = int_sterile_XA_set | int_sterile_AA_set
	ordered_combined_sterile_sets = sorted(combined_sterile_sets)
	sterile_calls = ordered_combined_sterile_sets[:steriles_needed]

	final_sterile_AA_count = []
	final_sterile_XA_count = []
	final_fertile_count = []
	final_fertile_focal_count = []

	for x in sterile_calls:
		if x in int_sterile_AA_set:
			final_sterile_AA_count.append(x)
		else:
			pass
	
	for x in sterile_calls:
		if x in final_sterile_XA_set:
			final_sterile_XA_count.append(x)
		else:
			pass

	int_fertile_focal_set = set(fertile_focal_count_AA)
	final_fertile_focal_set = int_fertile_focal_set - final_sterile_XA_set

	combined_fertile = set(fertile_count_XA + fertile_count_AA)
	int_fertile_list = combined_fertile - int_sterile_AA_set
	int_fertile_list2 = int_fertile_list - int_sterile_XA_set
	final_fertile_set = int_fertile_list2 - final_fertile_focal_set

	combined_fertile_sets = final_fertile_set | final_fertile_focal_set
	ordered_combined_fertile_sets = sorted(combined_fertile_sets)
	fertile_calls = ordered_combined_fertile_sets[:fertiles_needed]

	for x in fertile_calls:
		if x in final_fertile_focal_set:
			final_fertile_focal_count.append(x)
		else:
			pass
	
	for x in fertile_calls:
		if x in final_fertile_set:
			final_fertile_count.append(x)
		else:
			pass

	# Count groups for 2x2 fisher's exact test
	
	Fisher_SF = len(final_sterile_AA_count)
	Fisher_SNF = len(final_sterile_XA_count)
	Fisher_FF = len(final_fertile_focal_count)
	Fisher_FNF = len(final_fertile_count)

	odds, pvalue = sp.fisher_exact([[Fisher_SF, Fisher_SNF], [Fisher_FF, Fisher_FNF]])

	print(pvalue)

	File.close()
