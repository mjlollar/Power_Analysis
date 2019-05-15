# Python3.x compatibale 

import sys
import csv
import collections
import scipy.stats as sp
import itertools as it

csv.register_dialect('tab_delim', delimiter="\t", quoting=csv.QUOTE_NONE)

file_name = sys.argv[1] 
number_sequenced = int(sys.argv[2])

def read_lines(csv_reader, row_list):
	for row_number, row in enumerate(csv_reader):
		if row_number in row_list:
			yield row_number, row

def flatten(my_list):
    for value in my_list:
        if isinstance(value, collections.Iterable) and not isinstance(value, (str, bytes)):
            yield from flatten(value)
        else:
            yield value

def myrange1():
	for value in range(0,100):
		yield value
def myrange2():
	for value in range(0,199):
		yield value
def myrange3():
	for value in range(0,200):
		yield value

with open(file_name, 'r') as File:
	reader = csv.reader(File, dialect='tab_delim')
	r = list(range(0, number_sequenced))
	
	comparisons = []
	for row_number, row in read_lines(reader,r):
		row_tuples = list(it.permutations(row, 2))
		comparisons.append(row_tuples)

	first_value = 99
	second_value = 498
	set1 = []
	for value in myrange1():
		perms = list(range(first_value,second_value + 1))
		set1.append(perms)
		first_value = first_value + 499
		second_value = second_value + 499	
	adds1 = range(49900, 50000)
	set1.append(adds1)
	set1 = flatten(set1)
	set1 = list(set1)

	third_value = 50199
	fourth_value = 50498
	set2 = []
	for value in myrange2():
		perms = list(range(third_value,fourth_value + 1))
		set2.append(perms)
		third_value = third_value + 499
		fourth_value = fourth_value + 499
	adds2 = range(149500, 149700)
	set2.append(adds2)
	set2 = flatten(set2)
	set2 = list(set2)
	set2.remove(149300)

	fifth_value = 149700
	sixth_value = 149999
	set3 = []
	for value in myrange3():
		perms = list(range(fifth_value, sixth_value + 1))
		set3.append(perms)
		fifth_value = fifth_value + 499
		sixth_value = sixth_value +	499
	set3 = flatten(set3)
	set3 = list(set3)

	combined_sets = set1 + set2 + set3

	int_comparisons = []
	for value in combined_sets:
		for row in comparisons:
			int_comparisons.append(row[value])
	
	final_comparisons = [int_comparisons[i:i+number_sequenced] for i in range(0, len(int_comparisons), number_sequenced)]
	
	#calculate p-values for each window
	p_values = []
	for window in final_comparisons:
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
