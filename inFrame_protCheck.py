#Code checks whether amino acid sequences are in frame or not

import os
import csv
from Bio.Seq import Seq
from Bio import SeqIO

############################################################################################################################################################################################################################

states_output = open("PVstates_output.csv","a")
states_output.write("isolate,NEIS1518,NEIS0774,NEIS1310,NEIS1194,NEIS1297,NEIS0354,NEIS2011,NEIS1778,NEIS0213,NEIS1587 \n")

directory = r'Amino_acid_iso_Output'
files = os.listdir(directory)
#query against protein database
for f in files:
	blastp_cmd = "blastp -query " + "Amino_acid_iso_Output/"+f + " -evalue 0.5 -use_sw_tback -seg no -db /home/ac814/Documents/IRP_2020/BIGSdb/PVprot_aa/combined_PV_prot -outfmt '10 qseqid length qstart qend pident nident sseqid slen sstart send' -out " +f[:5]+"_blastp.txt"
	
	os.system(blastp_cmd)

	blastp_check_cmd = "blastp -query " + "Amino_acid_iso_Output/"+f + " -evalue 0.5 -use_sw_tback -seg no -db /home/ac814/Documents/IRP_2020/BIGSdb/PVprot_aa/combined_PV_prot -out " + f[:5]+"_check_blastp.txt"

	os.system(blastp_check_cmd)
	
#State is 0 unless the query has a hit in the database
	NEIS1518_state = 0
	NEIS0774_state = 0
	NEIS1310_state = 0
	NEIS1194_state = 0
	NEIS1297_state = 0
	NEIS0354_state = 0
	NEIS2011_state = 0
	NEIS1778_state = 0
	NEIS0213_state = 0
	NEIS1587_state = 0


	with open(f[:5]+"_blastp.txt") as prot_hit_file:
		contents = csv.reader(prot_hit_file)
		for i in contents:
						
			geneID = i[0]
			queryLength = i[1]
			hitStart = i[2]
			hitEnd = i[3]
			pIdent = i[4]
			nIdent = i[5]
			dbProtID = i[6]
			protLength = i[7]
			protStart = i[8]
			protEnd = i[9]

			if(int(protStart) == 1 and int(protEnd) == int(protLength)):
				print(str(geneID) + " is ON!")
				if(str(geneID) == 'NEIS1518_1'):
					NEIS1518_state = 1
				if(str(geneID) == 'NEIS0774_1'):
					NEIS0774_state = 1
				if(str(geneID) == 'NEIS1310_59'):
					NEIS1310_state = 1
				if(str(geneID) == 'NEIS1194_24'):
					NEIS1194_state = 1
				if(str(geneID) == 'NEIS1297_1'):
					NEIS1297_state = 1
				if(str(geneID) == 'NEIS0354_1'):
					NEIS0354_state = 1
				if(str(geneID) == 'NEIS2011_1'):
					NEIS2011_state = 1 
				if(str(geneID) == 'NEIS1778_1'):
					NEIS1778_state = 1
				if(str(geneID) == 'NEIS0213_1'):
					NEIS0213_state = 1
				if(str(geneID) == 'NEIS1587_1'):
					NEIS1587_state = 1

#optional print statement below may be useful but is a copy of what is outputted
	#print(str(f[:5]) + ' ' + str(NEIS1518_state) + ' ' + str(NEIS0774_state) + ' ' + str(NEIS1310_state) + ' ' + str(NEIS1194_state) + ' ' + str(NEIS1297_state) + ' ' + str(NEIS0354_state) + ' ' + str(NEIS2011_state) + ' ' + str(NEIS1778_state) + ' ' + str(NEIS0213_state) + ' ' + str(NEIS1587_state))

	states_output.write(str(f[:5]) + ',' + str(NEIS1518_state) + ',' + str(NEIS0774_state) + ',' + str(NEIS1310_state) + ',' + str(NEIS1194_state) + ',' + str(NEIS1297_state) + ',' + str(NEIS0354_state) + ',' + str(NEIS2011_state) + ',' + str(NEIS1778_state) + ',' + str(NEIS0213_state) + ',' + str(NEIS1587_state) + '\n')


states_output.close()			







