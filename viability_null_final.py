# Python3.x compatibale 

import sys
import csv
import scipy.stats as sp
import itertools as it

csv.register_dialect('tab_delim', delimiter="\t", quoting=csv.QUOTE_NONE)

#Input arguments taken
file_name = sys.argv[1] 
number_sequenced = int(sys.argv[2])
number_sequenced2 = int(number_sequenced*2)

# function to enumerate rows
def read_lines(csv_reader, row_list):
	for row_number, row in enumerate(csv_reader):
		if row_number in row_list:
			yield row_number, row

#functions used in deletion of within-chromosome comparisons
def myrange1():
	for value in reversed(range(401,500)):
		yield value
def myrange2():
	for value in reversed(range(201,402)):
		yield value
def myrange3():
	for value in reversed(range(1,202)):
		yield value

# Use every other fly in text file, increases "randomness" of selection
# although I am not sure if this really matters. Does decrease null p-values in some cases.
def drop(mylist, n):
	del mylist[::n]


with open(file_name, 'r') as File:
	reader = csv.reader(File, dialect='tab_delim')
	r = list(range(0, number_sequenced2))
	
	# Generate tuples of all pairwise window combinations and add to master list "comparisons"
	# Nested list, nested lists are indidivduals 
	comparisons = []
	for row_number, row in read_lines(reader, r):
		row_tuples = list(it.combinations(row, 2))
		comparisons.append(row_tuples)

	drop(comparisons, 2)

	# Delete all within-chromosome comparisons
	# 500c2 yields 124750 pairwise comparisons. Tuple generation works
	# by comparing window1-window2[0], window1-window3[1], .... window1-window500[499], window2-window3[500], window2-window4[501], etc.
	# Determined the index ranges that corresponds to comparisons within each chromosome and put them in list
	# X-chromosome comparisons
	first_value = 0
	second_value = 98
	set1 = []
	for value in myrange1():
		removed_comps = list(range(first_value,second_value + 1))
		set1.append(removed_comps)
		first_value = first_value + value
		second_value = second_value + value - 1

	#2nd chromosome comparisons
	third_value = 44149
	fourth_value = 44349
	set2 = []
	for value in myrange2():
		removed_comps = list(range(third_value,fourth_value + 1))
		set2.append(removed_comps)
		third_value = third_value + value
		fourth_value = fourth_value + value - 1

	#3rd chromosome comparisons
	fifth_value = 104449
	sixth_value = 104649
	set3 = []
	for value in myrange3():
		removed_comps = list(range(fifth_value, sixth_value + 1))
		set3.append(removed_comps)
		fifth_value = fifth_value + value
		sixth_value = sixth_value + value - 1	

	#Combine list of all index values corresponding to within chromosome comparisons and reverse order
	#Reversal necessary so that iteration doesn't disrupt tuple index
	combined_sets = set1 + set2 + set3
	combined_sets = list((list(it.chain(*combined_sets))))
	combined_sets = list(set(combined_sets))
	combined_sets = sorted(combined_sets, reverse=True)

	#Iterate through nested list to remove all within chromosome comparisons
	for i in range(0,number_sequenced):
		for value in combined_sets:
			del comparisons[i][value]

	#Rearange nested list so that each new nested list contains all sampled individuals for a given window 
	window_list_comparisons = list(map(list, zip(*comparisons)))

	#calculate p-values for each window
	p_values = []
	for window in window_list_comparisons:
		fishergroup1 = []
		fishergroup2 = []
		fishergroup3 = []
		fishergroup4 = []
		for individual in window:
			if individual==('0', '2'):
				fishergroup1.append(1)
			elif individual==('1', '2'):
				fishergroup2.append(1) 
			elif individual==('2', '2'):
				fishergroup2.append(1)
			elif individual==('0', '0'):
				fishergroup3.append(1)
			elif individual==('0', '1'):
				fishergroup3.append(1)
			elif individual==('1', '0'):
				fishergroup4.append(1)
			elif individual==('1', '1'): 
				fishergroup4.append(1)
			elif individual==('2', '0'):
				fishergroup4.append(1)
			elif individual==('2', '1'):
				fishergroup4.append(1)
			else:
				raise Exception("There is something wrong with your window values")
		oddsratio, pvalue = sp.fisher_exact([[len(fishergroup1), len(fishergroup2)], [len(fishergroup3), len(fishergroup4)]])
		p_values.append(pvalue)

	#Take minimum pvalue from all windows and print it
	lowest_pvalue = min(p_values)
	print(lowest_pvalue)

	File.close()
