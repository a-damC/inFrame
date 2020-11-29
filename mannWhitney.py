#Code that performs Mann-Whitney U test on the PV data

import pandas as pd

#######################################################################################################
#the 2 datasets I want to use
wholeGenome_df = pd.read_csv('wholeGenome_statsFormat.csv')
pvData_df = pd.read_csv('PVgenesData_statsFormat.csv')
#print(wholeGenome_df)

#######################################################################################################

#list of the isolates to iterate over
isolates = [54441, 54440, 54439, 54422, 54421, 54412, 54405, 54403, 51380, 51266, 51258, 51250,	51190, 51189, 51188, 51187, 51186, 51185, 51184, 51183, 51182, 51181, 51180, 51179, 51178, 51177, 51176, 51174, 51173, 51171, 51170, 51169, 51168, 51166, 51164, 51163, 51162, 51161, 51160, 51159, 51158, 51157, 51156, 51154, 51153, 51152, 51151, 51150, 51148, 51147, 51146, 51145]

########################################################################################################
#list which isolates are in which group

#GROUP_A
groupA_wg_inhabitants = [51145, 51147 , 51159 , 51160 , 51161 , 51162, 51163, 51164, 51169, 51180, 51187, 54403, 54405, 54439, 54440]
groupA_pv_inhabitants = ['51145_p', '51147_p' , '51159_p' , '51160_p' , '51161_p' , '51162_p' , '51163_p' , '51164_p' , '51169_p' , '51180_p' , '51187_p' , '54403_p' , '54405_p' , '54439_p' , '54440_p']

#GROUP_B
groupB_wg_inhabitants = [51168, 51171, 51182, 51189, 51258]
groupB_pv_inhabitants = ['51168_p', '51171_p', '51182_p', '51189_p', '51258_p']

#GROUP_C
groupC_wg_inhabitants = [51148, 51151, 51154, 51166, 51170, 51188]
groupC_pv_inhabitants = ['51148_p', '51151_p', '51154_p', '51166_p', '51170_p', '51188_p']

#GROUP_D
groupD_wg_inhabitants = [51146, 51156, 51157, 51179, 51183, 51184, 51250, 54412, 54422]
groupD_pv_inhabitants = ['51146_p', '51156_p', '51157_p', '51179_p', '51183_p', '51184_p', '51250_p', '54412_p','54422_p']

#GROUP_E
groupE_wg_inhabitants = [51158, 51181, 51185, 51186, 51266]
groupE_pv_inhabitants = ['51158_p', '51181_p', '51185_p', '51186_p', '51266_p']

#GROUP_F
groupF_wg_inhabitants = [51150, 51152, 51153, 51173, 51176, 51177, 51190, 51380]
groupF_pv_inhabitants = ['51150_p', '51152_p', '51153_p', '51173_p', '51176_p', '51177_p', '51190_p', '51380_p']

#GROUP_G
groupG_wg_inhabitants = [51174, 51178, 54421, 54441]
groupG_pv_inhabitants = ['51174_p', '51178_p', '54421_p', '54441_p']

########################################################################################################

#Output files

A_output = open('A_MannWhitney.csv','a')
B_output = open('B_MannWhitney.csv','a')
C_output = open('C_MannWhitney.csv','a')
D_output = open('D_MannWhitney.csv','a')
E_output = open('E_MannWhitney.csv','a')
F_output = open('F_MannWhitney.csv','a')
G_output = open('G_MannWhitney.csv','a')


########################################################################################################

#the for loop - to get U stat values for all 52 isolates
for i in isolates:
	
	print("----- " + str(i) + " -----")
	if(i in groupA_wg_inhabitants):
		wgGroup = groupA_wg_inhabitants
		pvGroup = groupA_pv_inhabitants
		output = A_output
		print(" group A used ")

	if(i in groupB_wg_inhabitants):
		wgGroup = groupB_wg_inhabitants
		pvGroup = groupB_pv_inhabitants
		output = B_output
		print(" group B used ")

	if(i in groupC_wg_inhabitants):
		wgGroup = groupC_wg_inhabitants
		pvGroup = groupC_pv_inhabitants
		output = C_output
		print(" group C used ")

	if(i in groupD_wg_inhabitants):
		wgGroup = groupD_wg_inhabitants
		pvGroup = groupD_pv_inhabitants
		output = D_output
		print(" group D used ")

	if(i in groupE_wg_inhabitants):
		wgGroup = groupE_wg_inhabitants
		pvGroup = groupE_pv_inhabitants
		output = E_output
		print(" group E used ")

	if(i in groupF_wg_inhabitants):
		wgGroup = groupF_wg_inhabitants
		pvGroup = groupF_pv_inhabitants
		output = F_output
		print(" group F used ")

	if(i in groupG_wg_inhabitants):
		wgGroup = groupG_wg_inhabitants
		pvGroup = groupG_pv_inhabitants
		output = G_output
		print(" group G used ")

	ID = str(i)

	info_output = open(ID +'_MWinfo.txt','a')
	
	#remove self from equation
	wgGroup.remove(i)
	pvGroup.remove(ID + '_p')
	info_output.write(ID + '\n' + str(wgGroup) + '\n' + str(pvGroup) + '\n')
########################################################################################################

	#Get the initial distance ranks for wg data  
	columnWG = wholeGenome_df.loc[0:,['isolate_ID',ID]]
	columnWG['distRank'] = columnWG[ID].rank(ascending = 1)
	#just get ranks for the specific group
	groupA_WGranks = columnWG.loc[columnWG['isolate_ID'].isin(wgGroup)]

	info_output.write('WG data distance rankings: ' + '\n' + str(groupA_WGranks) + '\n')


#########################################################################################################

	#Get the initial distance ranks for the pv data
	columnPV = pvData_df.loc[0:,['isolate_ID',ID]]
	columnPV['distRank'] = columnPV[ID].rank(ascending = 1)
	#just get ranks for the specific group
	groupA_PVranks = columnPV.loc[columnPV['isolate_ID'].isin(pvGroup)]

	info_output.write('PV data distance rankings: ' + '\n' + str(groupA_PVranks) + '\n')

########################################################################################################

	#creates table of PV and WG data, ready for the overall ranks
	concat_column = pd.concat([groupA_PVranks,groupA_WGranks])
	#ranks the ranks in new column
	concat_column['RankOfRank'] = concat_column['distRank'].rank(ascending = 1)

	info_output.write('The Rank of distance Ranks:  \n' + str(concat_column) + '\n')

########################################################################################################

	#Rank sum equation: n(n+1)/2 value
	n_plus1 = len(wgGroup) + 1
	n_times_above = len(wgGroup)*n_plus1
	equation_value = n_times_above/2
	print(equation_value)

	info_output.write('n(n+1)/2 = ' + str(equation_value) + '\n')

#########################################################################################################

	#Calcs U stat for wg data
	wgRanks = concat_column.loc[concat_column['isolate_ID'].isin(wgGroup)]
	wgRankSum = wgRanks['RankOfRank'].sum()
	print(wgRankSum)

	#(n(n+1))/2 (part of the mann whitney test)
	U_wg = wgRankSum - equation_value
	print('Uwg value = ' + str(U_wg))

########################################################################################################

	#Calcs U stat for pv data
	pvRanks = concat_column.loc[concat_column['isolate_ID'].isin(pvGroup)]
	pvRankSum = pvRanks['RankOfRank'].sum()
	print(pvRankSum)
	
	#(n(n+1))/2 (part of the mann whitney test)
	U_pv = pvRankSum - equation_value
	print('Upv value = ' + str(U_pv))

	info_output.write('Uwg value = ' + str(U_wg) + '\n' + 'Upv value = ' + str(U_pv) + '\n')

########################################################################################################
	
	#Ustat = smallest number of either U_wg and U_pv
	if(U_wg < U_pv):
		output.write(ID + ',' + str(U_wg) + '\n')
	else:
		output.write(ID + ',' + str(U_pv) + '\n')

########################################################################################################
	#re-add the self groups for the other iterations
	wgGroup.append(i)
	pvGroup.append(ID + '_p')

