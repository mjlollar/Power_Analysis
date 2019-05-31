# Null power analysis for the detection of Sterility QTLs
# Model assumes two unique two-locus recessive sterility loci. One non-focal X-Autosome locus and one focal Autosome-Autosome locus
# Model also assumes 25% penetrance of sterile loci.
# Allows you to input the number of flies you wish to sequence, plus those used to scan, adding in additional steriles for sequencing
# Inputs include generated samples from SIBSAM (Pool 2015) and number of flies to be seuqenced/phenotypically scanned
# Input should be in the format: python3 25pen_PA.py "Ancestry file.txt" # of sequenced +1" "# of screened + # sequenced + 1"
# Example input for 200 sequenced flies and 1200 scanned flies: python3 25pen_Null_PA.py F2MaleAnc2000_1.txt 201 1401 

import sys
import csv
import scipy.stats as sp
import itertools as it

csv.register_dialect('tab_delim', delimiter="\t", quoting=csv.QUOTE_NONE)

file_name = sys.argv[1] 
range_value1 = int(sys.argv[2])
range_value2 = int(sys.argv[3])

# function to enumerate rows
def read_lines(csv_reader, row_list):
	for row_number, row in enumerate(csv_reader):
		if row_number in row_list:
			yield row_number, row


with open(file_name, 'r') as File:
	# generate various lists
	sterile_focal_counts = []
	sterile_nonfocal_counts = []
	fertile_focal_counts = []
	fertile_nonfocal_counts = []
	combined_list_sterile = []
	combined_list_fertile =[]
	
	# read in file
	reader = csv.reader(File, dialect='tab_delim')
	r = list(range(0, range_value1))
	r2 = list(range(range_value1, range_value2))
	
	# Generate tuples of all pairwise window combinations and add to master list
	for row_number, row in read_lines(reader, r):
		row_tuples = list(it.combinations(row, 2))
		combined_list_sterile.append(row_tuples)
	
	# Rearrange list of tuples for each row such that windows are grouped by position across rows
	window_list_sterile = map(list, zip(*combined_list_sterile))
	
	# calculate proportion of sterile focal and non-focal for each window (where sterile replicate is defined randomly as first range)
	# output of this loop produces two new lists of focal/non-focal counts by window
	for window in window_list_sterile:
		sterile_focal = []
		sterile_nonfocal =[]
		for pair in window:
			if pair==('0', '2'):
				sterile_focal.append(1)
			else:
				sterile_nonfocal.append(1)
	# Sum focal/nonfocal values per row in list for single value per window
		sumsf = sum(sterile_focal)
		sumnsf = sum(sterile_nonfocal)
		sterile_focal_counts.append(sumsf)
		sterile_nonfocal_counts.append(sumnsf)

	
	# Reset reader iterator
	File.seek(0)
	# repeat similar iteration for fertile simulated replicates	
	for row_number, row in read_lines(reader, r2):
		row_tuples2 = list(it.combinations(row, 2))
		combined_list_fertile.append(row_tuples2)


	window_list_fertile = map(list, zip(*combined_list_fertile))

	for window in window_list_fertile:
		fertile_focal =[]
		fertile_nonfocal = []
		for pair in window:
			if pair==('0', '2'):
				fertile_focal.append(1)
			else:
				fertile_nonfocal.append(1)

		sumff = sum(fertile_focal)
		sumfnf = sum(fertile_nonfocal)
		fertile_focal_counts.append(sumff)
		fertile_nonfocal_counts.append(sumfnf)

# combine focal and non-focal sums to two lists of two values ([#focal, #nonfocal] for each of sterile and fertile)
sterile_tuples = map(list, zip(sterile_focal_counts, sterile_nonfocal_counts))
fertile_tuples = map(list, zip(fertile_focal_counts, fertile_nonfocal_counts))


# make 2x2 matrices (contingency) for each window
# Contingency tables for each window take on format of [[sterile-focal, sterile-nonfocal], [fertile-focal, fertile-nonfocal]]
fisher_groups = map(list, zip(sterile_tuples, fertile_tuples))

p_values = []

# calculate Fisher's Exact for each window and obtain p-values
for groups in fisher_groups:
	oddsratio, pvalue = sp.fisher_exact(groups)
	p_values.append(pvalue)

# Determine the lowest window p-value and print it
lowest_pvalue = min(p_values)
print(lowest_pvalue)

File.close()
