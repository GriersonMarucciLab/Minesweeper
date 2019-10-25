#!/usr/bin/env python

# Author: Joshua Rees, joshua.rees@bristol.ac.uk
# Affiliation: Life Sciences, University of Bristol 
# Last Updated: 2019-06-04

"""
Step 4 to X: The largest deletion segment determines the remaining low / no essentiality genes that have not been deleted. 
According to ranked results from INPUT_script_4X/conquerko_red_0_endtimes.txt, repeat for yellow and blue.
These remaining genes are divided into eight groups, a powerset generated for these eight groups, and each combination from the powerset  
individually appended to the current largest deletion combination and simulated. 
If none of these simulations produces a dividing cell, the remaining genes are appended as single knockouts to the current largest deletion 
combination and simulated. The individual remaining genes that don’t produce a dividing cell are temporarily excluded and a reduced remaining gene list produced.
Expects: Place N conquerko_COLOUR_N_endtimes.txt in INPUT_script_4X folder.
Output: either bash scripts and exp / ko lists in OUTPUT_script_4X folder OR final results in OUTPUT_final 
"""

# Imports
from sys import exit
from itertools import chain, combinations
import re
from datetime import datetime
import math
import os
import glob

# Global variables
#JOB = "What is this experiment/run of simulations called? "
NAME = 'mine'
JOB = NAME + 'gapko{}'
SIMNAMESTART = len(NAME)
SIMNAMEEND = (len(JOB))-2
SIMNAME = JOB[SIMNAMESTART:SIMNAMEEND]

def splashscreen():
	""" Fancy splash screen for Lego scripts """
	pass

def alreadyRunCheck():
	''' Checks to see if input files are present and output files have been generated, checks with the user'''

	# What round are we on? Input from conquerko = 1, input from gapko = n + 1, should work as the number of files decreases (i.e. red, yellow finish, leaving blue)
	fileexists = "OUTPUT_script_4X/bashexpkofiles/{}gapko_blue_1.sh"
	filexistname = fileexists.format(NAME)
	if not os.path.isfile(filexistname):
		pass
	if os.path.isfile(filexistname):
		filen = glob.glob("INPUT_script_4X/*_endtimes.txt")
		filenumber = len(filen)
		if filenumber <= 3:
			roundnumber = 0
		elif filenumber > 3:
			files = glob.glob("INPUT_script_4X/*_endtimes.txt")
			# getctime (creation) vs getmtime (modification), oldest to newest top to bottom
			files.sort(key=os.path.getmtime)
			files.reverse()
			topfile = files[0]
			roundnumber = topfile[-14]
			roundnumber = int(roundnumber)

		filex = glob.glob("OUTPUT_script_4X/bashexpkofiles/*.sh")
		filex.sort(key=os.path.getmtime)
		filex.reverse()
		topfile = filex[0]
		roundnum = topfile[-4]
		roundnum = int(roundnum)

		roundnumberadded = roundnumber + 1

		if roundnumberadded == roundnum:
			print(f"The latest endtimes.txt is numbered {roundnumber}, and the latest bash file is numbered {roundnum}." )
			print("\nThis increment of 1 suggests you have already run this script on this input. If you run it again the order of the powersets will randomise." )
			print("\nIf you are mid simulation this make the interpretation of your results incorrect." )
			response = input("\nDo you want to run this script? yes / no\n> " )
			if response == 'no':
				exit(0)

def interpretResults():
	""" Match the simulation results of Stage 3 / Stage 4 to the specific powerset combination, create unsortedsimresults and matchedresults.txt. """

	# What round are we on? Input from stage 3 = 1, input from stage 4 = n + 1, should work as the number of files decreases (i.e. red, yellow finish, leaving blue)
	filen = glob.glob("INPUT_script_4X/*_endtimes.txt")
	filenumber = len(filen)
	if filenumber <= 3:
		previousroundnumber = '0'
		roundnumber = '1'
	elif filenumber > 3:
		files = glob.glob("INPUT_script_4X/*_endtimes.txt")
		# getctime (creation) vs getmtime (modification), oldest to newest top to bottom
		files.sort(key=os.path.getmtime)
		files.reverse()
		topfile = files[0]
		roundnumber = topfile[-14]
		roundnumber = int(roundnumber)
		roundnumber = roundnumber + 1
		previousroundnumber = roundnumber - 1
		previousroundnumber = str(previousroundnumber)
		roundnumber = str(roundnumber)

	# previous round number used throughout to process the results of the last round > producing dividing_COLOUR_N, matchedresults_COLOUR_N, unsorted.._COLOUR_N
	# roundnumber used in the creation of the next scripts and {}eightsandgenes.txt

	## have to iterate through names just with indexes not named variables (dividingcounter_list is similar example below)

	#search for round results, should work as number of files decreases
	filesearch = 'INPUT_script_4X/*' + previousroundnumber + '_endtimes.txt'
	filesearching = glob.glob(filesearch)
	filesearchingnumber =  len(filesearching)
	if filesearchingnumber == 3:
		names = []
		names.append('_red_' + previousroundnumber)
		names.append('_yellow_' + previousroundnumber)
		names.append('_blue_' + previousroundnumber)
	if filesearchingnumber < 3:
		names = []
		for file in filesearching:
			filecolourextract1 = file[-18:-15]
			filecolourextract2 = file[-19:-15]
			filecolourextract3 = file[-21:-15]
			if filecolourextract1 == 'red':
				names.append('_red_' + previousroundnumber)
			if filecolourextract2 == 'blue':
				names.append('_blue_' + previousroundnumber)
			if filecolourextract3 == 'yellow':
				names.append('_yellow_' + previousroundnumber)

	#input
	previousroundnumber = int(previousroundnumber)
	if previousroundnumber == 0:
		explist_path = "OUTPUT_script_3/{}conquer{}_exp.list"
		combolist = "OUTPUT_script_3/{}_combosandgenes.txt"
		endtimes = "INPUT_script_4X/conquerko{}_endtimes.txt"
	elif previousroundnumber > 0:
		explist_path = "OUTPUT_script_4X/bashexpkofiles/{}gapko{}_exp.list"
		combolist = "OUTPUT_script_4X/{}_eightsandgenes.txt" #transition from combosandgenes to eightsandgenes
		endtimes = "INPUT_script_4X/gapko{}_endtimes.txt"
	previousroundnumber = str(previousroundnumber)

	#output
	unsortedsimresults = "OUTPUT_script_4X/unsortedsimresults{}.txt"
	matchedresults = "OUTPUT_script_4X/matchedresults{}.txt"
	deletionlog = "OUTPUT_final/deletionlog.txt"

	for name in names:
		full_results_list = []
		nthendtimes = endtimes.format(name)

		# open endtimes results and solve results being on two lines
		firstinput = open(nthendtimes).read()
		firstinput = re.sub(r"(\d)\n^(\w)", r"\1\t\2", firstinput, flags=re.MULTILINE)
		
		nthunsortedsimresults = unsortedsimresults.format(name)
		firstoutput = open(nthunsortedsimresults,"w+", newline ='\n')
		firstoutput.truncate()
		firstoutput.write(firstinput)
		firstoutput.close()

		# open unsortedsimresults.txt, create list and a list of equal length, account for intended length, insert sim result in location, 
		# append lists together, save as file
		secondinput = open(nthunsortedsimresults, "r+")
		result_list = [line.rstrip() for line in secondinput.readlines()]

		# create an equal length list of default 'No_Result'
		result_list_copy = []
		for line in result_list:
			result_list_copy.append('No_Result')

		# add additional lines to account for crashed sims
		explist = explist_path.format(NAME, name)
		lenexplist = sum(1 for line in open(explist))
		lenendtimes = sum(1 for line in open(nthunsortedsimresults))
		difference = (lenexplist - lenendtimes)
		for iteration in range(0, difference):
			result_list_copy.append('No_Result')

		# remove the control sims for results (the last two sims in every *exp.list), and the blank combination result (duplicate wildtype) produced by powerset	
		duplicatewildtype = 1
		if lenexplist >= 250:
			wildtype = 3000
			mutant = 3005
		elif lenexplist < 250:
			wildtype = lenexplist-1
			mutant = lenexplist

		# copy actual results into prepared list, by deleting and inserting at sorted location (row = sim number (0 indexed))
		for line in result_list:
			results = line.split("\t")
			sim = results[0]
			sim = int(sim)
			if sim == duplicatewildtype:
				continue
			elif sim == mutant:
				continue
			elif sim == wildtype:
				continue
			else: 
				#account for 0 index
				sim = sim-1
				time = results[1]
				outcome = results[2]
				#print(result_list_copy[sim-1])
				#print(result_list_copy[sim])
				#print(result_list_copy[sim+1])
				del result_list_copy[sim]
				#print(result_list_copy[sim-1])
				#print(result_list_copy[sim])
				#print(result_list_copy[sim+1])
				new_result = time + '\t' + outcome
				result_list_copy.insert(sim, new_result)
				#print(result_list_copy[sim-1])
				#print(result_list_copy[sim])
				#print(result_list_copy[sim+1])
				#input("Press Enter to continue...")

		# remove the 'NoResult' for wildtype and mutant sims that were not copied
		result_list_copy = result_list_copy[:-2]

		# append the output of this loop / list to full results list
		for line in result_list_copy:
			full_results_list.append(line)
		
		# zip fullresults list and gene_list.txt
		# as combosandgenes and eightsandgenes are only files that use name as prefix, not insert / suffix, so drop the first _ from _COLOUR_N = COLOUR_N
		name_variant = name[1:]
		combolist_path = combolist.format(name_variant)
		thirdinput = open(combolist_path).read()
		thirdinput = re.sub(r"(^[^']*)[^\n]*", r"\1", thirdinput, flags=re.MULTILINE)
		combo_list = thirdinput.split('\n')
		# make sure combos on same line are joined by spaces, not by tabs - so aren't seperated by split
		for line in combo_list:
			line = ' '.join(line)
		#combo_list = [line.rstrip() for line in thirdinput]
		full_results_list = zip(combo_list, full_results_list)

		# save matchedresults in outputdir
		nthmatchedresults = matchedresults.format(name)
		secondoutput = open(nthmatchedresults,"w+", newline ='\n')
		secondoutput.truncate()
		# convert zip created list of tuples, line by line into string 
		for line in full_results_list:
		  secondoutput.write('\t'.join(str(s) for s in line) + '\n')
		secondoutput.close()
		
		# there's rogue double tabbing at the beginning of the matched results file - removing it
		nthmatchedresults = matchedresults.format(name)
		roguetabtext = open(nthmatchedresults).read()
		roguetabtext = re.sub(r"\t\t", r"\t", roguetabtext)
		secondoutput = open(nthmatchedresults,"w+", newline ='\n')
		secondoutput.truncate()
		secondoutput.write(roguetabtext)
		secondoutput.close()
	
	print(f'\n\nHave matched the simulation results for {names} and saved the results in OUTPUT_script_4X/matchedresults.txt.')
	log = open(deletionlog,"a+", newline ='\n')
	log.write(f"\n\nHave matched the simulation results for {names} and saved the results in OUTPUT_script_4X/matchedresults.txt.")
	log.close()

	return names, roundnumber, previousroundnumber

def createDividingTxt(names):
	'''Check powerset combinations for success (i.e division) using matchedresults.txt.
	    Used by endingDecision() to determine the next step, and remainingGenes() for locating remaining genes.
		Creates dividing_COLOUR_N.txt'''

	matchedresults = "OUTPUT_script_4X/matchedresults{}.txt"
	dividingresults = "OUTPUT_script_4X/dividing{}.txt"
	deletionlog = "OUTPUT_final/deletionlog.txt"
	dividingcounter_list = []

	for name in names:
		dividingcounter = 0
		nthmatchedresults = matchedresults.format(name)
		results = open(nthmatchedresults,"r+")
		results_list = [line.rstrip() for line in results.readlines()]
		dividing_results = []
		divided = 'Divided'

		for line in results_list:
			if line.endswith(divided): 
				dividing_results.append(line)
				dividingcounter = dividingcounter + 1

		# save non_essential_genes in output dir
		nthdividingresults = dividingresults.format(name)
		neresults = open(nthdividingresults,"w+", newline ='\n')
		neresults.truncate() 
		for line in dividing_results:
			line = line + '\n'
			neresults.write(line)
		neresults.close()

		print(f'\nHave filtered the dividing results and saved in OUTPUT_script_4X/{nthdividingresults}.')
		log = open(deletionlog,"a+", newline ='\n')
		log.write(f"\nHave filtered the dividing results and saved in OUTPUT_script_4X/{nthdividingresults}.")
		log.close()
		dividingcounter_list.append(dividingcounter)

	return dividingcounter_list

def remainingGenes(names, dividingcounter_list, previousroundnumber):
	'''Calculate remaining genes from dividing results. Three Options:
		No division :: Remaining genes = re-record prior round results
		Division + Single KO Append Flag :: prior round was an appended single KO round, Remaining genes = reduced set that divided (when singly appended)
		Division + No Flag :: Remaining genes = of those that divided, select smallest number remaining (code starts at line 487) 
		Passed to endingDecision() with dividingcounter_list to determine next step.
		'''

	
	dividingcombos = "OUTPUT_script_4X/dividing{}.txt"
	remaininggenestxt = "OUTPUT_script_4X/remaininggenes{}.txt"
	deletedgenestxt = "OUTPUT_script_4X/deletedgenes{}.txt"
	deletionlog = "OUTPUT_final/deletionlog.txt"
	negenes = "OUTPUT_script_2/nonessential.txt"
	#format with names (but just red, yellow, blue = name_variant2)
	singleko_trigger = 'OUTPUT_final/singleko{}.txt'

	previousroundnumber = int(previousroundnumber)
	if previousroundnumber == 0:
		matchingcomboandgenes = "OUTPUT_script_3/{}_combosandgenes.txt"
	elif previousroundnumber > 0:
		matchingcomboandgenes = "OUTPUT_script_4X/{}_eightsandgenes.txt" #transition from combosandgenes to eightsandgenes
	previousroundnumber = str(previousroundnumber)
	
	for name in names:
		x = names.index(name)
		counter = dividingcounter_list[x]

		genecodes = []
		combonames = []
		matchedlist = []
		bestresult = []
		negenes_list = []
		negenes_deleted = []
		
		if counter == 0:
			name_variant = name[1:]
			pastname_components = name_variant.split('_')
			pastname_name = pastname_components[0]
			pastname_number = pastname_components[1]
			pastname_number = int(pastname_number)
			pastname_number = pastname_number - 1
			pastname_number = str(pastname_number)
			pastname_full = '_' + pastname_name + '_' + pastname_number
			
			# As no new division, copy remaining / deleted genes from past file to current file
			nthremaininggenes = remaininggenestxt.format(pastname_full)
			remaininggenes_results = open(nthremaininggenes).read()
			copypastremaininggenes = remaininggenestxt.format(name)
			copypastremaininggenestxt = open(copypastremaininggenes,"w+", newline ='\n')
			copypastremaininggenestxt.write(remaininggenes_results)
			copypastremaininggenestxt.close()
			
			# deletedgenes_results = past results, copy past deleted genes = copying past results into current results
			nthdeletedgenes = deletedgenestxt.format(pastname_full)
			deletedgenes_results = open(nthdeletedgenes).read()
			copypastdeletedgenes = deletedgenestxt.format(name)
			copypastdeletedgenestxt = open(copypastdeletedgenes,"w+", newline ='\n')
			copypastdeletedgenestxt.write(deletedgenes_results)
			copypastdeletedgenestxt.close()

			log = open(deletionlog,"a+", newline ='\n')
			log.write(f"\n\nNo dividing in-silico cells were produced in this {name} round, re-recording {name} last round's results.\n")
			log.close()

			continue

		elif counter > 0:
			name_variant2 = name[1:-2]
			nthsingleko_trigger = singleko_trigger.format(name_variant2)
			name_variant = name[1:]
			
			if os.path.isfile(nthsingleko_trigger):
				# single gene knockouts appended to the largest deletion, have been run, and produced division
				# this happened because the prior eights simulation round produced no division
				# so create deleted genes from prior deleted genes (which were passed from eights to singleko to here) (i.e. what the single kos were appended to)
				pastname_components = name_variant.split('_')
				pastname_name = pastname_components[0]
				pastname_number = pastname_components[1]
				pastname_number = int(pastname_number)
				pastname_number = pastname_number - 1
				pastname_number = str(pastname_number)
				pastname_full = '_' + pastname_name + '_' + pastname_number
				# deletedgenes_results = past results, copy past deleted genes = copying past results into current results
				nthdeletedgenes = deletedgenestxt.format(pastname_full)
				deletedgenes_results = open(nthdeletedgenes).read()
				copypastdeletedgenes = deletedgenestxt.format(name)
				copypastdeletedgenestxt = open(copypastdeletedgenes,"w+", newline ='\n')
				copypastdeletedgenestxt.write(deletedgenes_results)
				copypastdeletedgenestxt.close()
				
				# and create a reduced set of remaining genes (i.e. combining the appended single kos) 
				# temporarily exclude genes that did not divide when appended singly to largest deletion
				# load last rounds (singlekos) eights and combos
				# split to genecodes and code names
				# match dividing sims to genecodes via code names

				matchingcombosandgenestxt = matchingcomboandgenes.format(name_variant) 
				matchingcombosandgenesallcolumns = open(matchingcombosandgenestxt, "r+")
				matchingcombosandgenesallcolumns_list = [line.rstrip() for line in matchingcombosandgenesallcolumns.readlines()]
				
				# open combosandgenes / eightsandgenes (get combination name and gene codes)
				for line in matchingcombosandgenesallcolumns_list:
					if not line.strip():
						continue
					else:
						combo, genecode = line.split("\t")
						combonames.append(combo)
						genecode = re.sub(r"\n", r"", genecode)
						genecodes.append(genecode)
				
				# open dividing_COLOUR_N (get combination names that divided)
				dividingtxt = dividingcombos.format(name)
				dividingallcolumns = open(dividingtxt, "r+")
				dividingallcolumns_list = [line.rstrip() for line in dividingallcolumns.readlines()]
				
				# match only combination names that divided, with gene codes, via lists produced from combosandgenes/eightsandgenes
				for line in dividingallcolumns_list:
					results = line.split("\t")
					combostring = results[0]
					genelocation = combonames.index(combostring)
					genestring = genecodes[genelocation]
					matchedresult = combostring + '\t' + genestring + '\n'
					matchedlist.append(matchedresult)
				
				# produce list of genes that did divide when appended to largest deletion
				reducedremaininggenes = []
				for line in matchedlist:
					line = line.split("\t")
					genestring = line[1]
					# re.sub the last gene code (i.e. the appended single gene ko)
					genestring = re.sub(r".*('\w\w\w\w\w\w',\n)", r"\1", genestring)
					generesult = re.sub(r"\n", "", genestring)
					reducedremaininggenes.append(generesult)

				# create a list of genes that did not divide when appended to largest deletion
				nondividingappended = [item for item in genecodes if item not in reducedremaininggenes]
				nondividingappendedn = len(nondividingappended)
				remaininggenesn = len(reducedremaininggenes)
				negenes_deletedn = sum(1 for line in open(nthdeletedgenes)) 

				log = open(deletionlog,"a+", newline ='\n')
				log.write(f"\n In {name_variant} there were {remaininggenesn} genes found that could be appended to largest deletion and still divide, and {nondividingappended} genes that couldn't.\n")
				log.write(f"The reduced remaining genes are: {remaininggenes}\n")
				log.write(f"The temporarily excluded genes (potential conditional essentials) are: {nondividingappended}\n")
				log.close()

				remaininggenes = reducedremaininggenes

				nthremaininggenes = remaininggenestxt.format(name)
				remaininggenes_results = open(nthremaininggenes,"w+", newline ='\n')
				remaininggenes_results.truncate() 
				for line in remaininggenes:
					line = line + ' '
					remaininggenes_results.write(line)
				remaininggenes_results.close()

				continue

			elif not os.path.isfile(nthsingleko_trigger):
				# best deletion is selected, written to new deleted gene list
				# remaining genes is generated by comparing non essential genes with the new deleted gene list
				matchingcombosandgenestxt = matchingcomboandgenes.format(name_variant) 
				matchingcombosandgenesallcolumns = open(matchingcombosandgenestxt, "r+")
				matchingcombosandgenesallcolumns_list = [line.rstrip() for line in matchingcombosandgenesallcolumns.readlines()]
				
				# open combosandgenes / eightsandgenes (get combination name and gene codes)
				for line in matchingcombosandgenesallcolumns_list:
					if not line.strip():
						continue
					else:
						combo, genecode = line.split("\t")
						combonames.append(combo)
						genecode = re.sub(r"\n", r"", genecode)
						genecodes.append(genecode)
				
				# open dividing_COLOUR_N (get combination names that divided)
				dividingtxt = dividingcombos.format(name)
				dividingallcolumns = open(dividingtxt, "r+")
				dividingallcolumns_list = [line.rstrip() for line in dividingallcolumns.readlines()]
				
				# match only combination names that divided, with gene codes, via lists produced from combosandgenes/eightsandgenes
				for line in dividingallcolumns_list:
					results = line.split("\t")
					combostring = results[0]
					genelocation = combonames.index(combostring)
					genestring = genecodes[genelocation]
					matchedresult = combostring + '\t' + genestring + '\n'
					matchedlist.append(matchedresult)
				
				# get best (i.e. most genes deleted) out of matched list
				longestlength = 0
				for line in matchedlist:
					line = line.split("\t")
					combostring = line[0]
					combostring = re.sub(r"\n", r"", combostring)
					genestring = line[1]
					genestring = re.sub(r"\n", r"", genestring)
					templist = genestring.split(' ')
					# covert list to set then list = list(set(t)) 
					# this removes order, but order should come from remaining genes anyway / negenes.txt
					templist = list(set(templist)) 
					lencount = len(templist)
					genestring = ' '.join(templist)
					if lencount > longestlength:
						bestresult = combostring + '\t' + genestring + '\n'
						longestlength = lencount
				
				# compare best result vs OUTPUT_script_1/gene_list.txt
				negenes_txt = open(negenes, "r+")
				full_negenes_list = [line.rstrip() for line in negenes_txt.readlines()]
				for line in full_negenes_list:
					line = line.split("\t")
					line = line[0]
					line = re.sub(r"\n", r"", line)
					negenes_list.append(line)
				bestgeneresult = bestresult.split("\t") 
				bestgeneresult = bestgeneresult[1]
				bestgeneresult = re.sub(r"\n", "", bestgeneresult)
				bestresult_list = bestgeneresult.split(' ')
				
				for gene in negenes_list:
					for unorderedgene in bestresult_list:
						if gene == unorderedgene:
							negenes_deleted.append(gene)
				
				remaininggenes = [item for item in negenes_list if item not in negenes_deleted]
				remaininggenesn = len(remaininggenes)
				negenes_deletedn = len(negenes_deleted)

				log = open(deletionlog,"a+", newline ='\n')
				log.write(f"\n\nA {name} combination deleted {negenes_deletedn} genes, leaving {remaininggenesn} remaining genes.\n")
				#log.write(f"The largest number of genes deleted was by combination: {combostring}\n")
				# cannot get to work > i.e. Blue_0 prints 12.5h 25b 12.5b (sim 26), but the genes are from 25b 12.5g 12.5f (sim 18)
				log.write(f"The genes deleted: {negenes_deleted}\n")
				log.write(f"The remaining genes are: {remaininggenes}\n")
				log.close()

				nthremaininggenes = remaininggenestxt.format(name)
				remaininggenes_results = open(nthremaininggenes,"w+", newline ='\n')
				remaininggenes_results.truncate() 
				for line in remaininggenes:
					line = line + ' '
					remaininggenes_results.write(line)
				remaininggenes_results.close()

				nthdeletedgenes = deletedgenestxt.format(name)
				deletedgenes_results = open(nthdeletedgenes,"w+", newline ='\n')
				deletedgenes_results.truncate() 
				for line in negenes_deleted:
					line = line + ' '
					deletedgenes_results.write(line)
				deletedgenes_results.close()
	
def powerset(iterable):
	''' Daughter function of outputToLists() 
		Generates a powerset, all possible unique combinations of a set, in this case our top 3 largest deletions (colour variants) and matching deletion segments
		e.g. [1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)
		From this StackOverflow thread: https://stackoverflow.com/questions/464864/how-to-get-all-possible-combinations-of-a-list-s-elements
		Which pulls from this python documentation: https://docs.python.org/3/library/itertools.html#itertools-recipes 
		'''

	s = list(set(iterable))  
	return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def outputToLists(eightsnames_list, deletedgenes, name_variant, roundnumbername):
	''' Daughter function of eightPanelGroupingsGeneration() 
		Uses powerset() to generate combinations
		Records the names of 8 groups in the combinations, and matches with the genes within group
		Outputs COLOUR_eightsandgenes.txt (complete powerset)
	'''

	# Have to code around naming change that occurs here (number change, and _ (underscore) to handle): 
	# eightssegments uses name_variant (== previousroundnumber)
	# eightsandgenes uses roundnumbername (== current round number)
	# combosandgenes / eightsandgenes / eightssegments are only files that use name as prefix
	# not insert / suffix, drop the first _ from _COLOUR_N = COLOUR_N using name_variant
	
	# produce power sets of eightsnames
	# use to generate deletion lines by opening eightssegment (containing genes) with associated eightsnames
	eightssegment = "OUTPUT_script_4X/eightssegments/{}remaining_{}.txt" #name_variant, eightsname
	eightsandgenes = "OUTPUT_script_4X/{}_eightsandgenes.txt" #roundnumbername
	deletedgenestxt = "OUTPUT_script_4X/deletedgenes{}.txt"
	deletionlog = "OUTPUT_final/deletionlog.txt"
	powersetcombos_key = []
	genesdeleted = []
	deletedgenes_list = []
	powersetcombos = []

	### see http://book.pythontips.com/en/latest/enumerate.html for use of enumerate
	### use it here to loop over the powerset of the eight groups, skipping the first entry in the powerset (always = 0 components combination)
	for i, combo in enumerate(powerset(eightsnames_list), 1):
		# for each combination of 8 groups, record the names of the groups, to give a conjoined name for the combination
		combo = ' '.join(combo)
		powersetcombos_key.append(combo)
		
		# then using the seperate names of the groups
		combosplit = combo.split(" ")
		genesdeleted = []
		for item in combosplit:

			# empty strings are considered FALSE normally. 
			# Flip It - Is the string empty ( item.strip() )? IF NOT makes Yes = TRUE.
			# == if the item is just an empty string, append a space, and go to next iteration of loop 
			if not item.strip():
				genesdeleted.append(" ")
				continue
			
			#look up genes from eightssegment, using eightsname
			else:
				eightssegment_name = eightssegment.format(name_variant, item)
				eightssegmenttxt = open(eightssegment_name).read()
				genesdeleted.append(eightssegmenttxt)
		
		genesdeleted_txt = ''.join(genesdeleted)
		powersetcombos.append(genesdeleted_txt)

	name = '_' + name_variant
	nthdeletedgenes = deletedgenestxt.format(name)
	with open(nthdeletedgenes, 'r') as deletedgenes_txt:
		deletedgenes_txtall = deletedgenes_txt.read().replace('\n', '')
	linecount = 1

	for line in powersetcombos:
		if linecount == 1:
			deletedgenes_list.append('')
			linecount = linecount + 1
		elif linecount > 1:
			for line in deletedgenes_txtall:
				deletedgenes_list.append(deletedgenes_txtall)
	
	# PRODUCE eightsandgenes using roundnumbername RATHER than name_variant
	# CHANGE to eightsandgenes
	eightsandgenes_list = zip(powersetcombos_key, deletedgenes_list, powersetcombos)
	eightsandgenes_name = eightsandgenes.format(roundnumbername)
	eightsandgenes_txt = open(eightsandgenes_name,"w+", newline ='\n')
	eightsandgenes_txt.truncate() 
	# convert zip created list of tuples, line by line into string 
	for line in eightsandgenes_list:
		eightsandgenes_txt.write('\t'.join(str(s) for s in line) + '\n')
	eightsandgenes_txt.close()
	
	tempscript = open(eightsandgenes_name).read()
	tempscript = re.sub(r",'", ", '", tempscript, flags=re.MULTILINE)
	tempscript = re.sub(r",[\s]*'", ", '", tempscript, flags=re.MULTILINE)
	eightsandgenes_txt = open(eightsandgenes_name, 'w+', newline ='\n')
	eightsandgenes_txt.truncate()
	eightsandgenes_txt.write(tempscript)
	
	print(f'\nCreated and saved {eightsandgenes_name}.')
	log = open(deletionlog,"a+", newline ='\n')
	log.write(f"\nCreated and saved {eightsandgenes_name}.")
	log.close()

def eightsegmentwrite(name_variant, eights_name, eightsegment_list):
	''' Daughter function of eightPanelGroupingsGeneration() 
		Outputs COLOUR_remaining_N.txt (each of the 8 groups individually)'''

	eightssegment = "OUTPUT_script_4X/eightssegments/{}remaining_{}.txt" #name_variant, eightsname
	eightssegmenttxt = eightssegment.format(name_variant, eights_name)
	output_eightsegment = open(eightssegmenttxt,"w+", newline ='\n')
	output_eightsegment.truncate() 
	for line in eightsegment_list:
		output_eightsegment.write(line)
	output_eightsegment.close()

def eightPanelGroupingsGeneration(remaininggenes, deletedgenes, name, roundnumbername): 
	''' Daughter function of endingDecision() 
		Depending on the number of remaining genes, different code is executed:
		15 genes or greater :: works out the number of genes in each of the 8 groups (minimum 2, apart from final/eighth group) and allocates them
		14 genes to 8 genes :: hardcoded number of genes selected, maintaining 8 groups by transitioning the minimum number of genes per group from 2 to 1
		7 genes to 1 genes :: number of groups equal to number of genes, with each group containing one hardcoded gene
		Calls eightsegmentwrite() to save each of the eight groups of genes individually
		Calls outputToLists() to save complete powerset.

		This is a long function, collapse the if / elseif statements if using something like Sublime to read the code.
	'''

	eightssegment = "OUTPUT_script_4X/eightssegments/{}remaining_{}.txt" #name_variant, eightsname
	deletionlog = "OUTPUT_final/deletionlog.txt"
	clonelist = remaininggenes
	name_variant = name[1:] 
	roundnumbername_variant = roundnumbername[1:]
	name_variant2 = name[1:-2]

	whole = len(remaininggenes) # math.trunc used to round down > prevent overlapping gene knockouts 

	### Decision based on number of remaining genes
	if whole >= 15:
		percent = 12.5
		hundred = 100
		zero = 0
		ngenes = (whole/hundred)*percent
		ngenes = math.trunc(ngenes)
		nofgroups = 8
		groupcounter = 1
		# number of groups and their names in a list for def(outputToLists)
		eightsname_list = ['1', '2', '3', '4', '5', '6', '7', '8']

		for x in range(1, nofgroups + 1):
			if groupcounter == nofgroups:
				eights_name = str(x)
				eights_start = ngenes*(x-1) + (x-1)  
				eights_end = whole
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				break
			elif x == 1:
				eights_name = str(x)
				eights_start = zero
				eights_end = ngenes + 1
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 2:
				eights_name = str(x)
				eights_start = ngenes + 1
				eights_end = ngenes*2 + 2
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 3:
				eights_name = str(x)
				eights_start = ngenes*2 + 2
				eights_end = ngenes*3 + 3 
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 4:
				eights_name = str(x)
				eights_start = ngenes*3 + 3 
				eights_end = ngenes*4 + 4 
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 5:
				eights_name = str(x)
				eights_start = ngenes*4 + 4 
				eights_end = ngenes*5 + 5
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 6:
				eights_name = str(x)
				eights_start = ngenes*5 + 5
				eights_end = ngenes*6 + 6
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 7:
				eights_name = str(x)
				eights_start = ngenes*6 + 6 
				eights_end = ngenes*7 + 7
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
	
	##############################################

	elif whole == 14:
		zero = 0
		nofgroups = 8
		groupcounter = 1
		# number of groups and their names in a list for def(outputToLists)
		eightsname_list = ['1', '2', '3', '4', '5', '6', '7', '8']

		for x in range(1, nofgroups + 1):
			if groupcounter == nofgroups:
				eights_name = str(x)
				eights_start = 13  
				eights_end = whole
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				break
			elif x == 1:
				eights_name = str(x)
				eights_start = zero
				eights_end = 2
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 2:
				eights_name = str(x)
				eights_start = 2
				eights_end = 4
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 3:
				eights_name = str(x)
				eights_start = 4
				eights_end = 6 
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 4:
				eights_name = str(x)
				eights_start = 6 
				eights_end = 8
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 5:
				eights_name = str(x)
				eights_start = 8
				eights_end = 10
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 6:
				eights_name = str(x)
				eights_start = 10
				eights_end = 12
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 7:
				eights_name = str(x)
				eights_start = 12
				eights_end = 13
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue

	elif whole == 13:
		zero = 0
		nofgroups = 8
		groupcounter = 1
		# number of groups and their names in a list for def(outputToLists)
		eightsname_list = ['1', '2', '3', '4', '5', '6', '7', '8']

		for x in range(1, nofgroups + 1):
			if groupcounter == nofgroups:
				eights_name = str(x)
				eights_start = 12  
				eights_end = whole
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				break
			elif x == 1:
				eights_name = str(x)
				eights_start = zero
				eights_end = 2
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 2:
				eights_name = str(x)
				eights_start = 2
				eights_end = 4
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 3:
				eights_name = str(x)
				eights_start = 4
				eights_end = 6 
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 4:
				eights_name = str(x)
				eights_start = 6 
				eights_end = 8
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 5:
				eights_name = str(x)
				eights_start = 8
				eights_end = 10
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 6:
				eights_name = str(x)
				eights_start = 10
				eights_end = 11
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 7:
				eights_name = str(x)
				eights_start = 11
				eights_end = 12
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue

	elif whole == 12:
		zero = 0
		nofgroups = 8
		groupcounter = 1
		# number of groups and their names in a list for def(outputToLists)
		eightsname_list = ['1', '2', '3', '4', '5', '6', '7', '8']

		for x in range(1, nofgroups + 1):
			if groupcounter == nofgroups:
				eights_name = str(x)
				eights_start = 11  
				eights_end = whole
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				break
			elif x == 1:
				eights_name = str(x)
				eights_start = zero
				eights_end = 2
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 2:
				eights_name = str(x)
				eights_start = 2
				eights_end = 4
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 3:
				eights_name = str(x)
				eights_start = 4
				eights_end = 6 
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 4:
				eights_name = str(x)
				eights_start = 6 
				eights_end = 8
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 5:
				eights_name = str(x)
				eights_start = 8
				eights_end = 9
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 6:
				eights_name = str(x)
				eights_start = 9
				eights_end = 10
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 7:
				eights_name = str(x)
				eights_start = 10
				eights_end = 11
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue

	elif whole == 11:
		zero = 0
		nofgroups = 8
		groupcounter = 1
		# number of groups and their names in a list for def(outputToLists)
		eightsname_list = ['1', '2', '3', '4', '5', '6', '7', '8']

		for x in range(1, nofgroups + 1):
			if groupcounter == nofgroups:
				eights_name = str(x)
				eights_start = 10  
				eights_end = whole
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				break
			elif x == 1:
				eights_name = str(x)
				eights_start = zero
				eights_end = 2
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 2:
				eights_name = str(x)
				eights_start = 2
				eights_end = 4
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 3:
				eights_name = str(x)
				eights_start = 4
				eights_end = 6 
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 4:
				eights_name = str(x)
				eights_start = 6 
				eights_end = 7
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 5:
				eights_name = str(x)
				eights_start = 7
				eights_end = 8
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 6:
				eights_name = str(x)
				eights_start = 8
				eights_end = 9
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 7:
				eights_name = str(x)
				eights_start = 9
				eights_end = 10
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue

	elif whole == 10:
		zero = 0
		nofgroups = 8
		groupcounter = 1
		# number of groups and their names in a list for def(outputToLists)
		eightsname_list = ['1', '2', '3', '4', '5', '6', '7', '8']

		for x in range(1, nofgroups + 1):
			if groupcounter == nofgroups:
				eights_name = str(x)
				eights_start = 9  
				eights_end = whole
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				break
			elif x == 1:
				eights_name = str(x)
				eights_start = zero
				eights_end = 2
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 2:
				eights_name = str(x)
				eights_start = 2
				eights_end = 4
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 3:
				eights_name = str(x)
				eights_start = 4
				eights_end = 5 
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 4:
				eights_name = str(x)
				eights_start = 5 
				eights_end = 6
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 5:
				eights_name = str(x)
				eights_start = 6
				eights_end = 7
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 6:
				eights_name = str(x)
				eights_start = 7
				eights_end = 8
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 7:
				eights_name = str(x)
				eights_start = 8
				eights_end = 9
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue

	elif whole == 9:
		zero = 0
		nofgroups = 8
		groupcounter = 1
		# number of groups and their names in a list for def(outputToLists)
		eightsname_list = ['1', '2', '3', '4', '5', '6', '7', '8']

		for x in range(1, nofgroups + 1):
			if groupcounter == nofgroups:
				eights_name = str(x)
				eights_start = 8  
				eights_end = whole
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				break
			elif x == 1:
				eights_name = str(x)
				eights_start = zero
				eights_end = 2
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 2:
				eights_name = str(x)
				eights_start = 2
				eights_end = 3
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 3:
				eights_name = str(x)
				eights_start = 3
				eights_end = 4 
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 4:
				eights_name = str(x)
				eights_start = 4 
				eights_end = 5
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 5:
				eights_name = str(x)
				eights_start = 5
				eights_end = 6
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 6:
				eights_name = str(x)
				eights_start = 6
				eights_end = 7
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 7:
				eights_name = str(x)
				eights_start = 7
				eights_end = 8
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue

	elif whole == 8:
		zero = 0
		nofgroups = 8
		groupcounter = 1
		# number of groups and their names in a list for def(outputToLists)
		eightsname_list = ['1', '2', '3', '4', '5', '6', '7', '8']

		for x in range(1, nofgroups + 1):
			if groupcounter == nofgroups:
				eights_name = str(x)
				eights_start = 7  
				eights_end = whole
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				break
			elif x == 1:
				eights_name = str(x)
				eights_start = zero
				eights_end = 1
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 2:
				eights_name = str(x)
				eights_start = 1
				eights_end = 2
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 3:
				eights_name = str(x)
				eights_start = 2
				eights_end = 3 
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 4:
				eights_name = str(x)
				eights_start = 3 
				eights_end = 4
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 5:
				eights_name = str(x)
				eights_start = 4
				eights_end = 5
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 6:
				eights_name = str(x)
				eights_start = 5
				eights_end = 6
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 7:
				eights_name = str(x)
				eights_start = 6
				eights_end = 7
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue

	##############################################

	elif whole == 7:
		zero = 0
		nofgroups = whole
		groupcounter = 1
		# number of groups and their names in a list for def(outputToLists)
		eightsname_list = ['1', '2', '3', '4', '5', '6', '7']

		for x in range(1, nofgroups + 1):
			if groupcounter == nofgroups:
				eights_name = str(x)
				eights_start = 6  
				eights_end = whole
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				break
			elif x == 1:
				eights_name = str(x)
				eights_start = zero
				eights_end = 1
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 2:
				eights_name = str(x)
				eights_start = 1
				eights_end = 2
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 3:
				eights_name = str(x)
				eights_start = 2
				eights_end = 3 
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 4:
				eights_name = str(x)
				eights_start = 3 
				eights_end = 4
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 5:
				eights_name = str(x)
				eights_start = 4
				eights_end = 5
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 6:
				eights_name = str(x)
				eights_start = 5
				eights_end = 6
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue

	elif whole == 6:
		zero = 0
		nofgroups = whole
		groupcounter = 1
		# number of groups and their names in a list for def(outputToLists)
		eightsname_list = ['1', '2', '3', '4', '5', '6']

		for x in range(1, nofgroups + 1):
			if groupcounter == nofgroups:
				eights_name = str(x)
				eights_start = 5  
				eights_end = whole
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				break
			elif x == 1:
				eights_name = str(x)
				eights_start = zero
				eights_end = 1
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 2:
				eights_name = str(x)
				eights_start = 1
				eights_end = 2
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 3:
				eights_name = str(x)
				eights_start = 2
				eights_end = 3 
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 4:
				eights_name = str(x)
				eights_start = 3 
				eights_end = 4
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 5:
				eights_name = str(x)
				eights_start = 4
				eights_end = 5
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue

	elif whole == 5:
		zero = 0
		nofgroups = whole
		groupcounter = 1
		# number of groups and their names in a list for def(outputToLists)
		eightsname_list = ['1', '2', '3', '4', '5']

		for x in range(1, nofgroups + 1):
			if groupcounter == nofgroups:
				eights_name = str(x)
				eights_start = 4  
				eights_end = whole
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				break
			elif x == 1:
				eights_name = str(x)
				eights_start = zero
				eights_end = 1
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 2:
				eights_name = str(x)
				eights_start = 1
				eights_end = 2
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 3:
				eights_name = str(x)
				eights_start = 2
				eights_end = 3 
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 4:
				eights_name = str(x)
				eights_start = 3 
				eights_end = 4
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue

	elif whole == 4:
		zero = 0
		nofgroups = whole
		groupcounter = 1
		# number of groups and their names in a list for def(outputToLists)
		eightsname_list = ['1', '2', '3', '4']

		for x in range(1, nofgroups + 1):
			if groupcounter == nofgroups:
				eights_name = str(x)
				eights_start = 3  
				eights_end = whole
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				break
			elif x == 1:
				eights_name = str(x)
				eights_start = zero
				eights_end = 1
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 2:
				eights_name = str(x)
				eights_start = 1
				eights_end = 2
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 3:
				eights_name = str(x)
				eights_start = 2
				eights_end = 3 
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue

	elif whole == 3:
		zero = 0
		nofgroups = whole
		groupcounter = 1
		# number of groups and their names in a list for def(outputToLists)
		eightsname_list = ['1', '2', '3']

		for x in range(1, nofgroups + 1):
			if groupcounter == nofgroups:
				eights_name = str(x)
				eights_start = 2  
				eights_end = whole
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				break
			elif x == 1:
				eights_name = str(x)
				eights_start = zero
				eights_end = 1
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue
			elif x == 2:
				eights_name = str(x)
				eights_start = 1
				eights_end = 2
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue

	elif whole == 2:
		zero = 0
		nofgroups = whole
		groupcounter = 1
		# number of groups and their names in a list for def(outputToLists)
		eightsname_list = ['1', '2']

		for x in range(1, nofgroups + 1):
			if groupcounter == nofgroups:
				eights_name = str(x)
				eights_start = 1  
				eights_end = whole
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				break
			elif x == 1:
				eights_name = str(x)
				eights_start = zero
				eights_end = 1
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				groupcounter = groupcounter + 1
				continue

	elif whole == 1:
		zero = 0
		nofgroups = whole
		groupcounter = 1
		# number of groups and their names in a list for def(outputToLists)
		eightsname_list = ['1']

		for x in range(1, nofgroups + 1):
			if groupcounter == nofgroups:
				eights_name = str(x)
				eights_start = 1  
				eights_end = whole
				eightsegment_list = clonelist[eights_start:eights_end]
				eightsegmentwrite(name_variant, eights_name, eightsegment_list)
				break

	##############################################
				
	lastcreated = eightssegment.format(name_variant, eights_name)
	print(f'\n\nLast created 1/8ths deletion segments of remaining genes = {lastcreated}.')
	log = open(deletionlog,"a+", newline ='\n')
	log.write(f"\n\nLast created 1/8ths deletion segments of remaining genes = {lastcreated}.")
	log.close()

	outputToLists(eightsname_list, deletedgenes, name_variant, roundnumbername_variant)

def runappendsingleKOS(remaininggenes_list, deletedgenes_list, name, roundnumbername):
	''' Daughter function of endingDecision() 
		single gene knockouts are appended to the largest deletion with a successful dividing cell, to be run
		this is due to the prior eights simulation round producing no division (see createDividingTxt() and remainingGenes())
		Ouputs a singlekos version of COLOUR_eightsandgenes.txt
	'''

	eightsandgenes = "OUTPUT_script_4X/{}_eightsandgenes.txt" #roundnumbername
	deletedgenestxt = "OUTPUT_script_4X/deletedgenes{}.txt"
	deletionlog = "OUTPUT_final/deletionlog.txt"
	deletedgenes = []
	remaininggenes = []
	codename = []

	# produce deleted genes as a text block per line (not single gene per line aka from the list) (with starting line set to blank to lineup with powerset format)
	nthdeletedgenes = deletedgenestxt.format(name)
	with open(nthdeletedgenes, 'r') as deletedgenes_txt:
		deletedgenes_txtall = deletedgenes_txt.read().replace('\n', '')
	
	remaininggenes_list_len = len(remaininggenes_list) + 1

	linecount = 1
	for line in range(0, remaininggenes_list_len):
		if linecount == 1:
			deletedgenes.append('')
			linecount = linecount + 1
		elif linecount > 1:
			for line in deletedgenes_txtall:
				deletedgenes.append(deletedgenes_txtall)

	# produce remaining genes list (with starting line set to blank to lineup with powerset format)
	linecount = 1
	linenumber = 0
	for line in range(0, remaininggenes_list_len):
		if linecount == 1:
			remaininggenes.append('')
			linecount = linecount + 1
		elif linecount > 1:
			gene = remaininggenes_list[linenumber]
			remaininggenes.append(gene)
			linenumber = linenumber + 1

	# produce names for each line
	linecount = 1
	linenumber = 0
	for line in remaininggenes:
		if linecount == 1:
			remaininggenes.append('')
			linecount = linecount + 1
		elif linecount > 1:
			gene = remaininggenes[linenumber]
			gene = re.sub(r"'", "", gene)
			gene = re.sub(r",", "", gene)
			codename.append(gene)
			linenumber = linenumber + 1

	# PRODUCE eightsandgenes using roundnumbername RATHER than name_variant
	# CHANGE to eightsandgenes
	roundnumbername_variant = roundnumbername[1:]
	eightsandgenes_list = zip(codename, deletedgenes, remaininggenes)
	eightsandgenes_name = eightsandgenes.format(roundnumbername_variant)
	eightsandgenes_txt = open(eightsandgenes_name,"w+", newline ='\n')
	eightsandgenes_txt.truncate() 
	# convert zip created list of tuples, line by line into string 
	for line in eightsandgenes_list:
		eightsandgenes_txt.write('\t'.join(str(s) for s in line) + '\n')
	eightsandgenes_txt.close()
	
	tempscript = open(eightsandgenes_name).read()
	tempscript = re.sub(r",'", ", '", tempscript, flags=re.MULTILINE)
	tempscript = re.sub(r",[\s]*'", ", '", tempscript, flags=re.MULTILINE)
	eightsandgenes_txt = open(eightsandgenes_name, 'w+', newline ='\n')
	eightsandgenes_txt.truncate()
	eightsandgenes_txt.write(tempscript)
	
	print(f'\nCreated and saved a (singlekos) version of an eightsandgenes: {eightsandgenes_name}.')
	log = open(deletionlog,"a+", newline ='\n')
	log.write(f"\nCreated and saved a (singlekos) version of an eightsandgenes: {eightsandgenes_name}.")
	log.close()

def createScripts(roundnumbername, job, simname):
	''' Daughter function of endingDecision() 
		Convert template script using user input, create exp and ko txt files for COLOUR_eightsandgenes.txt'''

	user_input_txt = "INPUT_script_1/user_input.txt"
	templatescript = "templatescript/TemplateScript.sh"
	eightsandgenes = "OUTPUT_script_4X/{}_eightsandgenes.txt"
	kolist = "OUTPUT_script_4X/bashexpkofiles/{}_ko.list"
	explist = "OUTPUT_script_4X/bashexpkofiles/{}_exp.list"
	experimentscript = "OUTPUT_script_4X/bashexpkofiles/{}.sh"
	deletionlog = "OUTPUT_final/deletionlog.txt"

	jobname = job.format(roundnumbername)
	roundnumbername_variant = roundnumbername[1:]

	#create knock out txt file 
	eightsandgenes_clone = eightsandgenes.format(roundnumbername_variant)
	kolist_clone = kolist.format(jobname) 
	ko = open(kolist_clone, 'w+', newline ='\n')
	kotempscript = open(eightsandgenes_clone).read()
	kotempscript = re.sub(r"^[^']*'", "'", kotempscript, flags=re.MULTILINE)
	# replace blank combo from beginning, just removed with regex
	kotempscript = " \n" + kotempscript
	ko.write(kotempscript)
	
	kocount = sum(1 for line in open(kolist_clone)) 
	#IN this script don't append control simulations to gene candidates if powerset of 8
	if kocount < 250:
		ko.write("%s\n" % "")
		ko.write("%s\n" % "'MG_006',")

	#create exp txt file
	explist_clone = explist.format(jobname)
	exp = open(explist_clone, 'w+', newline ='\n')
	stopper = kocount
	for gene in range(0,stopper):
		exp.write("%s\n" % simname) #gapko
	
	if kocount < 250:
		#append control simulations to gene candidates
		exp.write("%s\n" % "wildtype")
		exp.write("%s\n" % "mutant")

	#create simulation bash script
	user_input_values = open(user_input_txt, "r+")
	tempscript = open(templatescript).read()
	experimentscript_clone = experimentscript.format(jobname)
	expscript = open(experimentscript_clone, 'w+', newline ='\n')

	value_list = [line.rstrip() for line in user_input_values.readlines()]
	user_input_values.close()
	value_list = value_list[0:10]

	### List of templatescript KEYS that need to change
	### Imported from user_input.txt
	### NAME
	name_value = value_list[0]
	name_value = name_value.split("\t")
	name_value = name_value[1]
	tempscript = re.sub(r"NAME", name_value, tempscript, flags=re.MULTILINE)
	### EMAIL
	email_value = value_list[1]
	email_value = email_value.split("\t")
	email_value = email_value[1]
	tempscript = re.sub(r"EMAIL", email_value, tempscript, flags=re.MULTILINE)
	### DEPARTMENT
	department_value = value_list[2]
	department_value = department_value.split("\t")
	department_value = department_value[1]
	tempscript = re.sub(r"DEPARTMENT", department_value, tempscript, flags=re.MULTILINE)
	### UNI
	uni_value = value_list[3]
	uni_value = uni_value.split("\t")
	uni_value = uni_value[1]
	tempscript = re.sub(r"UNI", uni_value, tempscript, flags=re.MULTILINE)
	### SUPERPC
	superpc_value = value_list[4]
	superpc_value = superpc_value.split("\t")
	superpc_value = superpc_value[1]
	tempscript = re.sub(r"SUPERPC", superpc_value, tempscript, flags=re.MULTILINE)
	### SUPERQUEUE
	superqueue_value = value_list[5]
	superqueue_value = superqueue_value.split("\t")
	superqueue_value = superqueue_value[1]
	tempscript = re.sub(r"SUPERQUEUE", superqueue_value, tempscript, flags=re.MULTILINE)
	### GROUP
	GROUP_value = value_list[6]
	GROUP_value = GROUP_value.split("\t")
	GROUP_value = GROUP_value[1]
	tempscript = re.sub(r"GROUP", GROUP_value, tempscript, flags=re.MULTILINE)
	### USER
	USER_value = value_list[7]
	USER_value = USER_value.split("\t")
	USER_value = USER_value[1]
	tempscript = re.sub(r"USER", USER_value, tempscript, flags=re.MULTILINE)
	### PROJECTFOLDER
	projectfolder_value = value_list[8]
	projectfolder_value = projectfolder_value.split("\t")
	projectfolder_value = projectfolder_value[1]
	tempscript = re.sub(r"PROJECTFOLDER", projectfolder_value, tempscript, flags=re.MULTILINE)
	
	### Decided by script
	### JOB is determined globally, so that users can see folder structure
	job_clone = jobname
	tempscript = re.sub(r"JOBx", job_clone, tempscript, flags=re.MULTILINE)
	
	### ENDTIMENAME default is 'e', producing 'e' + 'ndtimes.txt'
	ENDTIMENAME = 'gapko' + roundnumbername + '_e'
	endtimename_clone = ENDTIMENAME
	tempscript = re.sub(r"'e'", endtimename_clone, tempscript, flags=re.MULTILINE)

	### Calculate: DATE, SIMNUMBER
	### DATE
	DATE = "{:%Y-%m-%d}".format(datetime.now())
	tempscript = re.sub(r"DATE", DATE, tempscript, flags=re.MULTILINE)
	
	### SIMNUMBER
	if kocount < 250:
		SIMNUMBER = kocount + 2
	elif kocount >= 250:
		SIMNUMBER = kocount
	SIMNUMBER = str(SIMNUMBER)
	tempscript = re.sub(r"SIMNUMBER", SIMNUMBER, tempscript, flags=re.MULTILINE)

	expscript.write(tempscript)

	ko.close()
	exp.close()
	expscript.close()
		
	print('\nCreated bash script (*.sh), *_exp.lists, and *_ko.lists in OUTPUT_script_4X/bashexpkofiles/\nYou only need /gapko, /pdfs, and /figs internal folders.\n')
	log = open(deletionlog,"a+", newline ='\n')
	log.write(f"\nCreated bash script (*.sh), *_exp.lists, and *_ko.lists in OUTPUT_script_4X/bashexpkofiles/\nYou only need /gapko, /pdfs, and /figs internal folders.\n")
	log.close()
	#print('\nOne of your expected folder structures is:')
	#print(f'- projects\n\t- {GROUP_value}\n\t\t- {USER_value}\n\t\t\t- output\n\t\t\t\t- {projectfolder_value}\n\t\t\t\t\t- {job}\n\t\t\t\t\t\t- {simname}\n\t\t\t\t\t\t- wildtype\n\t\t\t\t\t\t- mutant\n\t\t\t\t\t\t- pdfs\n\t\t\t\t\t\t- figs')
	
def endingDecision(names, dividingcounter_list, JOB, SIMNAME, previousroundnumber, roundnumber):
	''' Given the outcome of createDividingTxt() and remainingGenes() determine next step for remaining variants (e.g. Red, Yellow, Blue).
		Four Options per variant:
		If remaining genes =< 8 in last simulation round, triggering final round, record final results (if all three variants complete = end of Minesweeper)
		If conducted singlegene kos for remaining genes, individually appended to largest deletion, in last simulation round
			and none divided, record final results
		If conducted a normal eight group powerset deletion in last simulation round, and none divided, start an appended single kos round
		If conducted a normal eights group powerset deletion in last simulation round, and some divided, 
			continue to another round / final round of eights (depending on n of remaining genes)

		Call eightPanelGroupingsGeneration(), runandcombinesingleKOs(), createScripts() as needed to:
		Outputs results in OUTPUT_final and next simulation round files in /bashexpkofiles folder
		'''

	remaininggenestxt = "OUTPUT_script_4X/remaininggenes{}.txt"
	deletedgenestxt = "OUTPUT_script_4X/deletedgenes{}.txt"
	#format with names (but just red, yellow, blue = name_variant2)
	singleko_trigger = 'OUTPUT_final/singleko{}.txt'
	finalround_trigger = 'OUTPUT_final/finalround{}.txt'
	finaloutput = 'OUTPUT_final/{}_finalresult.txt'
	deletionlog = "OUTPUT_final/deletionlog.txt"

	# name == previousroundnumber_name
	roundnumber_names = []

	for name in names:
		colourextract1 = name[1:-2]
		colourextract2 = name[1:-2]
		colourextract3 = name[1:-2]
		if colourextract1 == 'red':
			roundnumber_names.append('_red_' + roundnumber)
		elif colourextract2 == 'yellow':
			roundnumber_names.append('_yellow_' + roundnumber)
		elif colourextract3 == 'blue':
			roundnumber_names.append('_blue_' + roundnumber)

	for name in names:
		x = names.index(name)
		counter = dividingcounter_list[x]

		name_variant2 = name[1:-2]
		### use name_variant2 > as colour only for triggers
		nthsingleko_trigger = singleko_trigger.format(name_variant2)
		nthfinalround_trigger = finalround_trigger.format(name_variant2)
		nthfinaloutput = finaloutput.format(name_variant2)
		
		#get final results from remaininggenes / deleted genes
		remaininggenes_list = []
		nthremaininggenes = remaininggenestxt.format(name)
		remaininggenes_txt = open(nthremaininggenes, "r+")
		for line in remaininggenes_txt:
			linecomponents = line.split(' ')
			for components in linecomponents:
				remaininggenes_list.append(components)
		remaininggenes_list = list(filter(None, remaininggenes_list))

		deletedgenes_list = []
		nthdeletedgenes = deletedgenestxt.format(name)
		deletedgenes_txt = open(nthdeletedgenes, "r+")
		for line in deletedgenes_txt:
			linecomponents = line.split(' ')
			for components in linecomponents:
				deletedgenes_list.append(components)
		deletedgenes_list = list(filter(None, deletedgenes_list))

		#remaininggenesn = len(remaininggenes_list)
		#print(f"Name: {name} Counter: {counter} remaininggenes: {remaininggenes_list} n of remaininggenes: {remaininggenesn}")

		## if reached remaining genes =< 8 in last simulation round, triggering final round, record final results
		if os.path.isfile(nthfinalround_trigger):
			finaloutputtxt = open(nthfinaloutput,"w+", newline ='\n')
			finaloutputtxt.write(f"{name_variant2} deleted {deletedgenes_list}\n")
			finaloutputtxt.write(f"{name_variant2} did not delete {remaininggenes_list}\n")
			finaloutputtxt.close()

			print(f'\n\n{name_variant2} has finished its final round. See OUTPUT_final/{name_variant2}_finalresult.txt for result.')
			log = open(deletionlog,"a+", newline ='\n')
			log.write(f"\n\n{name_variant2} has finished its final round. See OUTPUT_final/{name_variant2}_finalresult.txt for result.")
			log.close()

			continue

		## if conducted singlegene kos for remaining genes individually appended to largest deletion in last simulation round, and non divided, record final results
		elif os.path.isfile(nthsingleko_trigger) and counter == 0:
			finaloutputtxt = open(nthfinaloutput,"w+", newline ='\n')
			finaloutputtxt.write(f"{name_variant2} deleted {deletedgenes_list}\n")
			finaloutputtxt.write(f"{name_variant2} deleted {remaininggenes_list}\n")
			finaloutputtxt.close()

			print(f'\n\n{name_variant2} has completed appended single ko sims, and with no dividing, has finished its final round.\nSee OUTPUT_final/{name_variant2}_finalresult.txt for result.\n')
			log = open(deletionlog,"a+", newline ='\n')
			log.write(f"\n\n{name_variant2} has completed appended single ko sims, and with no dividing, has finished its final round.\nSee OUTPUT_final/{name_variant2}_finalresult.txt for result.\n")
			log.close()

			continue

		## if conducted a normal eight group powerset deletion in last simulation round, and none divided, start an appended single kos round
		## also write a single ko trigger to interpret the results differently in the following round
		elif not os.path.isfile(nthsingleko_trigger) and counter == 0:
			singleko_triggertxt = open(nthsingleko_trigger,"w+", newline ='\n')
			singleko_triggertxt.write("\n")
			singleko_triggertxt.close()

			### use roundnumbername > as future
			roundnumbername = roundnumber_names[x]
			#### use name > as past / current results
			runappendsingleKOS(remaininggenes_list, deletedgenes_list, name, roundnumbername)

			### use roundnumbername > as future
			#!!!! gets input from {}eightsandgenes.txt / {}singlekopowersets.txt OR list returned from prior function?
			createScripts(roundnumbername, JOB, SIMNAME)
			
			print(f'\n{name_variant2} produced 0 divisions and has more than 8 genes remaining, so progressing to appended single ko sims.\n')
			log = open(deletionlog,"a+", newline ='\n')
			log.write(f"\n{name_variant2} produced 0 divisions and has more than 8 genes remaining, so progressing to appended single ko sims.\n")
			log.close()

			continue

		# if conducted a normal eights deletion in last simulation round, and some divided, continue to another round / final round of eights (depending on n of remaining genes)
		elif counter > 0:
			remaininggenesn = len(remaininggenes_list)
			if remaininggenesn <= 8:
				#write finalround trigger into folder
				finalround_triggertxt = open(nthfinalround_trigger,"w+", newline ='\n')
				finalround_triggertxt.write("\n")
				finalround_triggertxt.close()
				
				#delete singleko_trigger if exists
				if os.path.isfile(nthsingleko_trigger):
					os.remove(nthsingleko_trigger)
				
				### use roundnumbername > as future
				roundnumbername = roundnumber_names[x]
				#### use name > as past / current results
				eightPanelGroupingsGeneration(remaininggenes_list, deletedgenes_list, name, roundnumbername)

				### use roundnumbername > as future
				#!!!! gets input from {}eightsandgenes.txt / {}singlekopowersets.txt OR list returned from prior function?
				createScripts(roundnumbername, JOB, SIMNAME)
				
				print(f'\n{name_variant2} has less than 8 genes remaining, so progressing to final round.')
				log = open(deletionlog,"a+", newline ='\n')
				log.write(f"\n{name_variant2} has less than 8 genes remaining, so progressing to final round.")
				log.close()

				continue	

			elif remaininggenesn > 8:
				#delete singleko_trigger if exists
				if os.path.isfile(nthsingleko_trigger):
					os.remove(nthsingleko_trigger)

				### use roundnumbername > as future
				roundnumbername = roundnumber_names[x]
				#### use name > as past / current results
				eightPanelGroupingsGeneration(remaininggenes_list, deletedgenes_list, name, roundnumbername)

				### use roundnumbername > as future
				#!!!! gets input from {}eightsandgenes.txt / {}singlekopowersets.txt OR list returned from prior function?
				createScripts(roundnumbername, JOB, SIMNAME)
				
				print(f'\n{name_variant2} has more than 8 genes remaining, so continuing onto next round.')
				log = open(deletionlog,"a+", newline ='\n')
				log.write(f"\n{name_variant2} has more than 8 genes remaining, so continuing onto next round.")
				log.close()

				continue

	# Completion message
	finalyellow = finaloutput.format('yellow')
	finalred = finaloutput.format('red')
	finalblue = finaloutput.format('blue')
	if os.path.isfile(finalyellow) and os.path.isfile(finalred) and os.path.isfile(finalblue):
		print(f'\n\nMinesweeper has finished :) see OUTPUT_final for results.')
		log = open(deletionlog,"a+", newline ='\n')
		log.write(f"\n\nMinesweeper has finished :) see OUTPUT_final for results.")
		log.close()

def main():
	splashscreen()
	# use previousroundnumber
	alreadyRunCheck()
	names, roundnumber, previousroundnumber = interpretResults()
	dividingcounter_list = createDividingTxt(names)
	remainingGenes(names, dividingcounter_list, previousroundnumber)
	endingDecision(names, dividingcounter_list, JOB, SIMNAME, previousroundnumber, roundnumber)
	#Called by endingDecision: eightPanelGroupingsGeneration()
	#Called by endingDecision: runandcombinesingleKOs()
	# use roundnumber > as producing new scripts
	#Called by endingDecision: createScripts()

main()