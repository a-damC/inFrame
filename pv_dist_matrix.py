#Script that makes distance matrix from a csv file of values - outputs a nexus file to be used for SplitsTree


import csv

############################################################################################################################################################################################################################
#Opens csv file
with open("practice.csv") as csv_file:
	reader = csv.reader(csv_file)
	#print(reader)

#These 2 list_tuples make up the x and y of the distance matrix	
	list_tuples1 = []
	list_tuples2 = []
	for row in reader:
		list_tuples1.append(tuple(row))

############################################################################################################################################################################################################################		
#Create nexus file
nexus_file = open("distance_matrix.nex", "a")

nexus_file.write("#NEXUS \n[Distance matrix calculated by Adam's Basic PV Comparator] \n[Based on the BIGSdb Genome Comparator by Jolley & Maiden 2010 BMC Bioinformatics 11:595]\n\n\n")

nexus_file.write("BEGIN taxa;\n \tDIMENSIONS ntax = " + str(len(list_tuples1)) + ";")	

nexus_file.write("\n\nEND;")

nexus_file.write("\n\nBEGIN distances;\n \tDIMENSIONS ntax = " + str(len(list_tuples1)) + ";\n \tFORMAT\n\t\ttriangle=LOWER\n\t\tdiagonal\n\t\tlabels\n\t\tmissing=?\n\t\t;\nMATRIX")

############################################################################################################################################################################################################################
#makes second list as it goes to get the pyramid shape - always end with distance to self i.e. 0
#matrix_list only used to print what is being currently compared
for row1 in list_tuples1:
	list_tuples2.append(tuple(row1))
	nexus_file.write("\n" + str(row1[0]) + "\t")
	for row2 in list_tuples2:
		matrix_list = []
		matrix_list.append(row1)
		matrix_list.append(row2)
		#print(matrix_list)



#Only performs distance calculation if the index numbers are the same - which is the gene column - therefore only compares the same genes between isolates.
#T_distance keeps running total of distance. i & j count the index positions, with j being reset to 0 after each iteration.
#na_count keeps track of n/a for weighted calculation at the end
		T_distance = 0
		na_count = 0
		i = 0
		for state1 in row1[1:]:
			j = 0
			for state2 in row2[1:]:
				if i == j:
					i += 1
					j = len(row1)*2
					if state1 != 'n/a' and state2 != 'n/a':
						dist = int(state1) - int(state2)
						if dist<0:
							distance = dist*(-1)
						else:
							distance = dist
					else:
						na_count += 1
						distance = 0
						
					T_distance = T_distance + distance
				else:
					j += 1
				
		nexus_file.write(str(T_distance) + "\t")
		

nexus_file.write("\n\t ; \nEND;")

############################################################################################################################################################################################################################
