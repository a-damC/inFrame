#Code used to extract coordinates of gene from sequence file using BLAST search

import os
import csv
from Bio.Seq import Seq
from Bio import SeqIO
########################################################################################################################################################################################################################################

#collection of gene sequences from seq data - to be manipulated in jalview
NEIS1518_1_output = open('gene_output_NEIS1518_1.txt','a')
NEIS0774_1_output = open('gene_output_NEIS0774_1.txt','a')
NEIS1310_59_output = open('gene_output_NEIS1310_59.txt','a')
NEIS1194_24_output = open('gene_output_NEIS1194_24.txt','a')
NEIS1297_1_output = open('gene_output_NEIS1297_1.txt','a')
NEIS0354_1_output = open('gene_output_NEIS0354_1.txt','a')
NEIS2011_1_output = open('gene_output_NEIS2011_1.txt','a')
NEIS1778_1_output = open('gene_output_NEIS1778_1.txt','a')
NEIS0213_1_output = open('gene_output_NEIS0213_1.txt','a')
NEIS1587_1_output = open('gene_output_NEIS1587_1.txt','a')
error_gene_output = open('error_gene_output.txt','a')

count = 0
directory = r'SeqFilesTested'
files = os.listdir(directory)
for f in files:
	count += 1
	print(f + '  ..... file number:' + str(count))
	#output file
	gene_coord_file = open(f+"_gene_coord.txt","a")
	#blastn command to create csv file with useful info
	blastn_cmd = "blastn -query " + 'SeqFiles/'+f + " -gapextend 2 -db /home/ac814/Documents/IRP_2020/BIGSdb/PVgenes_ntd/combined_PV_genes -outfmt '10 qseqid length qstart qend sseqid slen sstart send' -out " +f+"_hits_blastn.txt"

	os.system(blastn_cmd)
	
	blastn_check_cmd = "blastn -query " + 'SeqFiles/'+f + " -gapextend 2 -db /home/ac814/Documents/IRP_2020/BIGSdb/PVgenes_ntd/combined_PV_genes -out " +f+"_checkHit_blastn.txt"
	os.system(blastn_check_cmd)
	
	aa_output = open(f[0:5]+'_aa_output', 'a')

########################################################################################################################################################################################################################################

	#opens file created by blastn command
	with open(str(f)+"_hits_blastn.txt") as hit_file:
		contents = csv.reader(hit_file)	

		#used to get list of gene names for later use in next block
		gene_name_list = []
		genelistCount = 0
		for i in contents:
			gene_name_list.append(i[4])
			genelistCount += 1
			if genelistCount == len(open(f+'_hits_blastn.txt').readlines()):
				break


	#re-open this file for next part of code
	with open(str(f)+"_hits_blastn.txt") as hit_file:
		contents = csv.reader(hit_file)
	
	
#gene_coord_list stores gene coordinates in the isolate sequence
#side list is used to compare genes with 2 hits
#count is used to initialise the loop and prevent out of range error on "if(gene_id not in gene_cord_list[len(gene_cord_list)-1][4]):"
#delete list for genes that have more than 2 hits - delete the rest of the hits, most likely to be short unuseful seqs

		gene_cord_list = [['initialiser','initialiser','initialiser','initialiser','initialiser']]
		side_list = []
		delete_list = [['initialiser']]
		count = 0
	
		for line in contents:
			isolate_id = line[0]
			alignment_length = line[1]
			iso_start = line[2]
			iso_end = line[3]
			gene_id = line[4]
			gene_length = line[5]
			gene_start = line[6]
			gene_end = line[7]
			#print(gene_id)
			#print(gene_name_list.count(gene_id))

			#Goes through here if gene hasn't already been added and not in the delete list
			if(gene_id not in gene_cord_list[len(gene_cord_list)-1][4] and gene_id not in delete_list[len(delete_list)-1][0]):
				#this is for if the alignment is greater than or equal to the gene length
				if(int(alignment_length) >= int(gene_length)):
					gene_cord_list.append(line)
					gene_coord_file.write(str(line))
					
				
				#is the alignment smaller than gene length and is this the only hit of this gene id (checks to gene list made from earlier)
				if(int(alignment_length) < int(gene_length) and gene_name_list.count(gene_id) == 1):
					if(int(gene_start) == 1 or int(gene_end) == 1 and int(gene_start) == int(gene_length) or int(gene_end) == int(gene_length)):
						gene_cord_list.append(line)
					else:
						positionCount = 0
						#does the sequence contain the first letters of the hit
						if(int(gene_start) == 1 or int(gene_end) == 1):
							positionCount += 1
							if(int(gene_start) == 1):
								isoStart = iso_start
							if(int(gene_end) == 1):
								isoStart = iso_end
								positionCount += 2
							
						#does the sequence contain the last letters of the hit
						if(int(gene_start) == int(gene_length) or int(gene_end) == int(gene_length)):
							positionCount += 2
							if(int(gene_start) == int(gene_length)):
								isoEnd = iso_start
								positionCount += 2
							if(int(gene_end) == int(gene_length)):
								isoEnd = iso_end
												
						#gene_end value == to gene_length, find value for isoStart which should be smallest isolate value minus smallest gene value
						if(positionCount == 2):
							smallestIsoValue = min(int(iso_start),int(iso_end))
							smallestGeneValue = min(int(gene_start), int(gene_end))
							safeSide = smallestGeneValue*2
							isoStart = smallestIsoValue - safeSide
						
						#gene_start value == to gene_length i.e. gene flipped relative to one in database need to find value for isoStart which should be largest isolate value plus smallest gene value
						if(positionCount == 4):
							highestIsoValue = max(int(iso_start),int(iso_end))
							smallestGeneValue = min(int(gene_start), int(gene_end))
							safe_Side = smallestGeneValue*2
							isoStart = highestIsoValue + safe_Side	
								
						#isoStart == iso_start, need find where end is therefore lookfor largest iso value and add the difference of (gene length - highest gene value)
						if(positionCount == 1):
							highestIsoValue = max(int(iso_start),int(iso_end))
							highestGeneValue = max(int(gene_start), int(gene_end))
							dist_to_geneLength = int(gene_length) - highestGeneValue
							safe_Side = dist_to_geneLength*2
							isoEnd = highestIsoValue + safe_Side
						#isoStart == iso_end, query seq is reversed, as such find end value by finding smallest isovalue and minusing (gene Length - highest gene value) 	
						if(positionCount == 3):
							smallestIsoValue = min(int(iso_start),int(iso_end))
							highestGeneValue = max(int(gene_start), int(gene_end))
							dist_to_geneLength = int(gene_length) - highestGeneValue
							safe_Side = dist_to_geneLength*2
							isoStart = smallestIsoValue - safeSide
						
						#adds adjusted hit to the final output list gene_cord_list
						adjusted_Hit = [isolate_id,alignment_length, str(isoStart), str(isoEnd),gene_id, gene_length, gene_start, gene_end]
						gene_cord_list.append(adjusted_Hit)
						gene_coord_file.write(str(adjusted_Hit))
			
			#does gene have 3 or more blastn hits? will have been added to delete list. Most likely hit to be something short and not useful
			if(gene_id in delete_list[len(delete_list)-1][0]):
				continue

			#goes through here if query hit is smaller than gene length and its not alread been added and there are more than one hit for this gene
			if(gene_id not in gene_cord_list[len(gene_cord_list)-1][4] and int(alignment_length)<int(gene_length) and gene_name_list.count(gene_id) > 1):
				side_list.append(line)
				#comes through here if there are 2 lines in the side_list 
				if(len(side_list) == 2):
					positionCount = 0
					#Is gene's position 1 included in one of the hits?
					if(int(side_list[0][6]) == 1 or int(side_list[0][7]) == 1 or int(side_list[1][6]) == 1 or int(side_list[1][7] == 1)):
						positionCount += 1
						if(int(side_list[0][6]) == 1):
							isoStart = side_list[0][2]
						elif(int(side_list[0][7]) == 1):
							isoStart = side_list[0][3]
							positionCount += 2
						elif(int(side_list[1][6]) == 1):
							isoStart = side_list[1][2]
						else:
							isoStart = side_list[1][3]
							positionCount += 2

		        	#Is gene's last position included in any of the hits?
					if(int(side_list[0][6]) == int(gene_length) or int(side_list[0][7]) == int(gene_length) or int(side_list[1][6]) == int(gene_length) or int(side_list[1][7] == int(gene_length))):
						positionCount += 2
						if(int(side_list[0][6]) == int(gene_length)):
							isoEnd = side_list[0][2]
							positionCount += 2
						elif(int(side_list[0][7]) ==int(gene_length)):
							isoEnd = side_list[0][3]
						elif(int(side_list[1][6]) == int(gene_length)):
							isoEnd = side_list[1][2]
							positionCount += 2
						else:
							isoEnd = side_list[1][3]
					
					#positionCount if even means we have isoEnd, if odd we have isoStart, if 0 then need both isoStart and isoEnd

					#have isoEnd, find isoStart by finding smallest iso value and minus smallest gene value
					if(positionCount == 2):
						smallestIsoValue = min(int(side_list[0][2]),int(side_list[0][3]),int(side_list[1][2]),int(side_list[1][3]))
						smallestGeneValue = min(int(side_list[0][6]),int(side_list[0][7]),int(side_list[1][6]),int(side_list[1][7]))
						safeSide = smallestGeneValue*2
						isoStart = str(smallestIsoValue - safeSide)
					#have isoEnd but gene is reversed, find isoStart by finding largest iso value and adding smallest gene value
					if(positionCount == 4):
						highestIsoValue = max(int(side_list[0][2]),int(side_list[0][3]),int(side_list[1][2]),int(side_list[1][3]))
						smallestGeneValue = min(int(side_list[0][6]),int(side_list[0][7]),int(side_list[1][6]),int(side_list[1][7]))
						safe_Side = smallestGeneValue*2
						isoStart = highestIsoValue + safe_Side
					
					#have isoStart, find isoEnd by looking for largest iso value and adding difference of (gene length - highest gene value)
					if(positionCount == 1):
						highestIsoValue = max(int(side_list[0][2]),int(side_list[0][3]),int(side_list[1][2]),int(side_list[1][3]))
						highestGeneValue = max(int(side_list[0][6]),int(side_list[0][7]),int(side_list[1][6]),int(side_list[1][7]))
						dist_to_geneLength = int(gene_length) - highestGeneValue
						safe_Side = dist_to_geneLength*2
						isoEnd = str(highestIsoValue + safe_Side)
					#have isoStart but gene is reversed, find isoEnd by looking for smallest iso value and minusing the difference of (gene length - highest gene value)
					if(positionCount == 3):
						smallestIsoValue = min(int(side_list[0][2]),int(side_list[0][3]),int(side_list[1][2]),int(side_list[1][3]))
						highestGeneValue = max(int(side_list[0][6]),int(side_list[0][7]),int(side_list[1][6]),int(side_list[1][7]))
						dist_to_geneLength = int(gene_length) - highestGeneValue
						safe_Side = dist_to_geneLength*2
						isoStart = smallestIsoValue - safe_Side
				
					#need to find both isoStart and isoEnd
				#first see if gene is reversed or not (assumption here is that if one hit is reversed, both will be as the hit is usually split on the repeatative tract as dust is turned on which masks these low complex regions)	
					if(positionCount == 0):
						if(int(side_list[0][2]) < int(side_list[0][3]) and int(side_list[0][6]) < int(side_list[0][7])):
							smallestIsoValue = min(int(side_list[0][2]),int(side_list[0][3]),int(side_list[1][2]),int(side_list[1][3]))
							smallestGeneValue = min(int(side_list[0][6]),int(side_list[0][7]),int(side_list[1][6]),int(side_list[1][7]))
							safeSide = smallestGeneValue*2
							isoStart = str(smallestIsoValue - safeSide)

							highestIsoValue = max(int(side_list[0][2]),int(side_list[0][3]),int(side_list[1][2]),int(side_list[1][3]))
							highestGeneValue = max(int(side_list[0][6]),int(side_list[0][7]),int(side_list[1][6]),int(side_list[1][7]))
							dist_to_geneLength = int(gene_length) - highestGeneValue
							safe_Side = dist_to_geneLength*2
							isoEnd = str(highestIsoValue + safe_Side)
						
						
						else:
							highestIsoValue = max(int(side_list[0][2]),int(side_list[0][3]),int(side_list[1][2]),int(side_list[1][3]))
							smallestGeneValue = min(int(side_list[0][6]),int(side_list[0][7]),int(side_list[1][6]),int(side_list[1][7]))
							safe_Side = smallestGeneValue*2
							isoStart = highestIsoValue + safe_Side

							smallestIsoValue = min(int(side_list[0][2]),int(side_list[0][3]),int(side_list[1][2]),int(side_list[1][3]))
							highestGeneValue = max(int(side_list[0][6]),int(side_list[0][7]),int(side_list[1][6]),int(side_list[1][7]))
							dist_to_geneLength = int(gene_length) - highestGeneValue
							safe_Side = dist_to_geneLength*2
							isoEnd = str(smallestIsoValue - safe_Side)

					combinedHit = [side_list[0][0],side_list[0][1],str(isoStart),str(isoEnd),side_list[0][4],side_list[0][5],side_list[0][6],side_list[0][7]]
					gene_cord_list.append(combinedHit)
					gene_coord_file.write(str(combinedHit))
					
					delete_list.append(side_list[0][4])
					side_list.clear()
						
					
	gene_coord_file = open(str(isolate_id)+"_gene_coord.txt","a")
	gene_coord_file.write(str(gene_cord_list[1:]))
	gene_coord_file.close
	
	print(gene_cord_list[1][0])
	for i in gene_cord_list[1:]:
		print(str(i[4]) + ' ' + str(int(i[2]) - int(i[3])) + ' ' + str(i[5]) + '--- isolate start and end ' + str(i[2]) + ' ' + str(i[3]))


#########################################################################################################################################################################################################################################
#Below codes the part that will extract the genes from the sequence files 
	
	for row in gene_cord_list[1:]:
		#print(row)

		isolate_id = row[0]
		alignment_length = row[1]
		iso_start = min(row[2],row[3])
		iso_end = max(row[3],row[2])
		gene_id = row[4]
		gene_length = row[5]
		gene_start = row[6]
		gene_end = row[7]	

		print(str(row[4]) + ' -Gene length in iso: ' + str((int(row[3])+1) - int(row[2])) + ' ,Gene length in database: ' + str(row[5]) + ' -------- ' + str(row[6]) + ' ' + str(row[7]))


		for seq_record in SeqIO.parse('SeqFiles/'+f, 'fasta'):
			iso_id = seq_record.id
			sequence = seq_record.seq
			
			gene_coord_contents = sequence[(int(iso_start)-1):int(iso_end)]
			if(int(gene_start) > int(gene_end)):
				actual_gene = gene_coord_contents.reverse_complement()
				
			else:
				actual_gene = gene_coord_contents
			
			
			if(gene_id == "NEIS1518_1"):
				NEIS1518_1_output.write(">"+isolate_id+"\n"+str(actual_gene)+"\n")

			elif(gene_id == "NEIS0774_1"):
				NEIS0774_1_output.write(">"+isolate_id+"\n"+str(actual_gene)+"\n")

			elif(gene_id == "NEIS1310_59"):
				NEIS1310_59_output.write(">"+isolate_id+"\n"+str(actual_gene)+"\n")

			elif(gene_id == "NEIS1194_24"):
				NEIS1194_24_output.write(">"+isolate_id+"\n"+str(actual_gene)+"\n")

			elif(gene_id == "NEIS1297_1"):
				NEIS1297_1_output.write(">"+isolate_id+"\n"+str(actual_gene)+"\n")

			elif(gene_id == "NEIS0354_1"):
				NEIS0354_1_output.write(">"+isolate_id+"\n"+str(actual_gene)+"\n")

			elif(gene_id == "NEIS2011_1"):
				NEIS2011_1_output.write(">"+isolate_id+"\n"+str(actual_gene)+"\n")
		
			elif(gene_id == "NEIS1778_1"):
				NEIS1778_1_output.write(">"+isolate_id+"\n"+str(actual_gene)+"\n")

			elif(gene_id == "NEIS0213_1"):
				NEIS0213_1_output.write(">"+isolate_id+"\n"+str(actual_gene)+"\n")

			elif(gene_id == "NEIS1587_1"):
				NEIS1587_1_output.write(">"+isolate_id+"\n"+str(actual_gene)+"\n")
	
			else:
				error_gene_output.write(">"+isolate_id+"\n"+str(actual_gene)+"\n")


#######################################################################################################################################################################################################################################
#Below codes part that take ntd seq and finds the amino acid seq
			
			table = 11
			min_prot_len = float(gene_length)*0.18
			
			for strand, nuc in [(+1, actual_gene)]:
				#the three reading frames for each direction
				for frame in range(3):
					#need a sequence that is divisible by 3. // is a divide that returns an integer and not a floating point
					length = 3 * ((len(actual_gene)-frame) // 3)
					#for loop that runs through strand in 3 reading frames and splits everytime it comes across a stop codon
					for prot in nuc[frame:frame+length].translate(table).split("*"):
						if len(prot) >= min_prot_len:
							print("%s...%s - length %i, strand %i, frame %i" %(prot[:30], prot[-3:], len(prot), strand, frame))
							aa_output.write(">"+ str(gene_id) + "\n" + str(prot) + "\n")

#######################################################################################################################################################################################################################################

	
