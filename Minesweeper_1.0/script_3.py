#!/usr/bin/env python

# Author: Joshua Rees, joshua.rees@bristol.ac.uk
# Affiliation: Life Sciences, University of Bristol 
# Last Updated: 2019-06-04

"""
Step 3: Results from Stage 2 are ranked, the largest deletion segment produced is matched with all other dividing, non-overlapping segments. 
This is controlled by logic in variantCombinations().
Powersets are used to generate combinations of the matched segements for simulation. 
Deletion segments that do not prevent division go to the Stage 4.
Expects: Place divideko_endtimes.txt in INPUT_script_3 folder.
Output: combination list, remaining genes list, multiple bash scripts and exp / ko lists, in OUTPUT_script_2 folder
"""

# Imports
from sys import exit
from itertools import chain, combinations
import re
from datetime import datetime
import os

# Global variables
#JOB = "What is this experiment/run of simulations called? "
name = 'mine'
JOB = name + 'conquer'
SIMNAMESTART = len(name)
SIMNAMEEND = (len(JOB))
SIMNAME = JOB[SIMNAMESTART:SIMNAMEEND]

def splashscreen():
	""" Fancy splash screen for Lego scripts """

	pass

def alreadyRunCheck():
	''' Checks to see if input files are present and output files have been generated, checks with the user'''

	# What round are we on? Input from conquerko = 1, input from gapko = n + 1, should work as the number of files decreases (i.e. red, yellow finish, leaving blue)
	inputfile = "INPUT_script_3/divideko_endtimes.txt"
	outputfile = "OUTPUT_script_3/{}conquer_blue_0.sh"
	outputfile = outputfile.format(name)

	if os.path.isfile(inputfile) and os.path.isfile(outputfile):
		print(f"You have expected files in your INPUT and OUTPUT folders for this script." )
		print("\nThis suggests you have already run this script on this input. If you run it again the order of the powersets will randomise." )
		print("\nIf you are mid simulation this make the interpretation of your results incorrect." )
		response = input("\nDo you want to run this script? yes / no\n> " )
		if response == 'no':
			exit(0)

def interpretResults():
	""" Match the simulation results of Stage 2 to the % board, create unsortedsimresults and matchedresults.txt. """

	explist = "OUTPUT_script_2/{}divide_exp.list"
	explist = explist.format(name)
	boardlist = "OUTPUT_script_2/alldivisionsegments_codes.txt"
	endtimes = "INPUT_script_3/divideko_endtimes.txt"
	unsortedsimresults = "OUTPUT_script_3/unsortedsimresults.txt"
	matchedresults = "OUTPUT_script_3/matchedresults.txt"
	full_results_list = []

	# open endtimes results and solve results being on two lines
	firstinput = open(endtimes).read()
	firstinput = re.sub(r"(\d)\n^(\w)", r"\1\t\2", firstinput, flags=re.MULTILINE)
		
	firstoutput = open(unsortedsimresults,"w+", newline ='\n')
	firstoutput.truncate()
	firstoutput.write(firstinput)
	firstoutput.close()

	# open unsortedsimresults.txt, create list and a list of equal length, account for intended length, insert sim result in location, append lists together, save as file
	secondinput = open(unsortedsimresults, "r+")
	result_list = [line.rstrip() for line in secondinput.readlines()]

	# create an equal length list of default 'No_Result'
	result_list_copy = []
	for line in result_list:
		result_list_copy.append('No_Result')

	# add additional lines to account for crashed sims
	lenexplist = sum(1 for line in open(explist))
	lenendtimes = sum(1 for line in open(unsortedsimresults))
	difference = (lenexplist - lenendtimes)
	for iteration in range(0, difference):
		result_list_copy.append('No_Result')
		
	# remove the control sims for results (the last two sims in every *exp.list)
	wildtype = lenexplist-1
	mutant = lenexplist

	# copy actual results into prepared list, by deleting and inserting at sorted location (row = sim number (0 indexed))
	for line in result_list:
		results = line.split("\t")
		sim = results[0]
		sim = int(sim)
		if sim == mutant:
			continue
		elif sim == wildtype:
			continue
		else: 
			#account for 0 index
			sim = sim-1
			time = results[1]
			outcome = results[2:]
			outcome = '\t'.join(outcome)
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
	thirdinput = open(boardlist, "r+")
	board_list = [line.rstrip() for line in thirdinput.readlines()]
	full_results_list = zip(board_list, full_results_list)

	# save matchedresults in outputdir
	secondoutput = open(matchedresults,"w+", newline ='\n')
	secondoutput.truncate()

	# convert zip created list of tuples, line by line into string 
	for line in full_results_list:
	  secondoutput.write('\t'.join(str(s) for s in line) + '\n')
	secondoutput.close()

	print('\nHave matched the simulation results with % boards and saved the results in OUTPUT_script_3/matchedresults.txt.')
	deletionlog = "OUTPUT_final/deletionlog.txt"
	log = open(deletionlog,"a+", newline ='\n')
	log.write("\nHave matched the simulation results with % boards and saved the results in OUTPUT_script_3/matchedresults.txt.\n")
	log.close()


def dividedAndProducedProteinRNA(line):
	linetocheck = line
	divided = 'Divided'
	outcome = 'UP'

	if linetocheck.endswith(divided) and linetocheck.count(outcome) == 2:
		return True
	else:
		return False

def successCheckAndVariants():
	'''Check boards for success (i.e division) using matchedresults.txt (ordered largest to smallest when generated) and save top 3 largest deletions.
		Top 3 are assigned to colours (red, yellow, blue) and used as three variants / avenues of deletion going forward.
		Creates successfulboards.txt + variants.txt.'''

	matchedresults = "OUTPUT_script_3/matchedresults.txt"
	variants = "OUTPUT_script_3/variants.txt"
	successfulboards = "OUTPUT_script_3/successfulboards.txt"
	divisionsegment_path = "OUTPUT_script_2/divisionsegment{}.txt"
	deletionlog = "OUTPUT_final/deletionlog.txt"
	
	# add additional names to names list if you want to increase number of variants
	# also add additional names in createScripts(), line ~367
	names = ['red_0', 'yellow_0', 'blue_0']

	boards = open(matchedresults, "r+")
	board_list = [line.rstrip() for line in boards.readlines()]
	variants_list = []
	success_list = []
	divided = 'Divided'
	counter = 0

	for line in board_list:
		if dividedAndProducedProteinRNA(line) is True:
			line = line.split("\t")
			line = line[0]
			success_list.append(line)
			if counter < 3: 
				variants_list.append(line)
				counter = counter + 1

	counter = 0
	for line in success_list:
		line = re.sub(r"(\d\d)\W(\d)", r"\1_\2", line)
		success_division = divisionsegment_path.format(line) 
		success_division_contents = open(success_division).read()
		success_division_contents = re.sub(r"\n", r"", success_division_contents)
		line = re.sub(r"(\d\d)_(\d)", r"\1.\2", line)
		del success_list[counter]
		new_line = line + '\t' + success_division_contents + '\n'
		success_list.insert(counter, new_line)
		counter = counter + 1

	# save succesful boards in output dir
	successfulboards_txt = open(successfulboards,"w+", newline ='\n')
	successfulboards_txt.truncate() 
	# convert zip created list of tuples, line by line into string 
	for line in success_list:
	  successfulboards_txt.write(line)
	successfulboards_txt.close()
	
	named_variants = zip(names, variants_list)
	# save variants in output dir
	variants_txt = open(variants,"w+", newline ='\n')
	log = open(deletionlog,"a+", newline ='\n')
	log.write('Largest initial deletion segments:\n')
	variants_txt.truncate() 
	# convert zip created list of tuples, line by line into string 
	for line in named_variants:
	  variants_txt.write('\t'.join(str(s) for s in line) + '\n')
	  log.write('\t'.join(str(s) for s in line) + '\n')
	variants_txt.close()
	log.close()

	successlist = success_list

	print('\nHave filtered the % boards results for successful division and largest 3 division segments, and saved in OUTPUT_script_3/successfulboards.txt + variants.txt.')
	deletionlog = "OUTPUT_final/deletionlog.txt"
	log = open(deletionlog,"a+", newline ='\n')
	log.write("\nHave filtered the % boards results for successful division and largest 3 division segments.\n")
	log.close()

	return successlist

def powerset(iterable):
	''' Daughter function of outputToLists() 
		Generates a powerset, all possible unique combinations of a set, in this case our top 3 largest deletions (colour variants) and matching deletion segments
		e.g. [1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)
		From this StackOverflow thread: https://stackoverflow.com/questions/464864/how-to-get-all-possible-combinations-of-a-list-s-elements
		Which pulls from this python documentation: https://docs.python.org/3/library/itertools.html#itertools-recipes 
		'''
	
	s = list(set(iterable))  # set option = no duplicate elements
	return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def outputToLists(components, boardnames, genecodes, red_yellow_blue):
	''' Daughter function of variantCombinations() 
		Uses powerset() to generate combinations
		Records the names of % boards in the combinations, and matches with the combined deleted genes
		Outputs COLOUR_combosandgenes.txt'''

	combos_genes = "OUTPUT_script_3/{}_combosandgenes.txt"
	successfulcomponents = []
	powersetcombos_key = []
	genesdeleted = []
	powersetcombos = []

	### This checks the list of segments that have been logic matched to the largest deletion, against a 
	###	a list of segments that produce a successfully dividing cell. Only those on both lists get kept
	for component in components:
		for board in boardnames:
			if component == board:
				successfulcomponents.append(component)
			else:
				continue

	### see http://book.pythontips.com/en/latest/enumerate.html for use of enumerate
	### use it here to loop over the powerset of the successful components, skipping the first entry in the powerset (always = 0 components combination)
	for i, combo in enumerate(powerset(successfulcomponents), 1):
		# for each combination of matched segments, record the names of the segments, to give a conjoined name for the combination
		combo = ' '.join(combo)
		powersetcombos_key.append(combo)
		
		# then using the seperate names of the segments 
		combosplit = combo.split(" ")
		genesdeleted = []
		for item in combosplit:
			
			# empty strings are considered FALSE normally. 
			# Flip It - Is the string empty ( item.strip() )? IF NOT makes Yes = TRUE.
			# == if the item is just an empty string, append a space, and go to next iteration of loop  
			if not item.strip():
				genesdeleted.append(" ")
				continue
			
			# use the name of the segment to find the index location (using the boardnames list)
			# use the index to recall the gene codes of that segment / board
			else:
				x = boardnames.index(item)
				genesdeleted.append(genecodes[x])
		
		# record the genes
		genesdeleted_txt = ''.join(genesdeleted)
		powersetcombos.append(genesdeleted_txt)

	# combine the recorded names and recorded genes
	combosandgenes = zip(powersetcombos_key, powersetcombos)
	combosandgenes_name = combos_genes.format(red_yellow_blue)
	combosandgenes_txt = open(combosandgenes_name,"w+", newline ='\n')
	combosandgenes_txt.truncate() 
	# convert zip created list of tuples, line by line into string 
	for line in combosandgenes:
	  combosandgenes_txt.write('\t'.join(str(s) for s in line) + '\n')
	combosandgenes_txt.close()

	print(f'\nCreated and saved {combosandgenes_name}.')
	deletionlog = "OUTPUT_final/deletionlog.txt"
	log = open(deletionlog,"a+", newline ='\n')
	log.write(f"\nCreated and saved {combosandgenes_name}.\n")
	log.close()

def variantCombinations(successlist):
	'''The Variants (three largest deletion segments) are matched with all other dividing, non-overlapping segments
		Using logic outlined below. Removing 100% of genes ends Minesweeper. Removing 90% of the genes ends this script, moving onto the next'''
	variants = "OUTPUT_script_3/variants.txt"
	top3 = open(variants, "r+")
	top3_list = [line.rstrip() for line in top3.readlines()]
	boardnames = []
	genecodes = []

	for line in successlist:
		boardname, genecode = line.split("\t")
		boardnames.append(boardname)
		genecode = re.sub(r"\n", r"", genecode)
		genecodes.append(genecode)

	for line in top3_list:
		line = line.split("\t")
		red_yellow_blue = line[0]
		boardcode = line[1]

		if boardcode == '100':
			print("Final results have been written in OUTPUT_final/finalresults.txt")
			finalresults = "OUTPUT_final/finalresults.txt"
			finalresults_txt = open(finalresults,"w+", newline ='\n')
			finalresults_txt.truncate() 
			for line in successlist:
			  finalresults_txt.write(line)
			  break
			finalresults_txt.close()

			deletionlog = "OUTPUT_final/deletionlog.txt"
			log = open(deletionlog,"a+", newline ='\n')
			log.write("\n\nMinesweeper is finished :) Final results have been written in OUTPUT_final/finalresults.txt. ")
			log.close()

			exit()
		
		elif boardcode == '90a':
			## If already removed 90% of the genes, going to create resulting files without running the code / simulations,
			## to spoof the next script into working.  
			explist_path = "OUTPUT_script_3/mineconquer_red_0_exp.list"
			combolist = "OUTPUT_script_3/red_0_combosandgenes.txt"
			endtimes = "INPUT_script_4X/conquerko_red_0_endtimes.txt"

			explist_path_txt = open(explist_path,"w+", newline ='\n')
			explist_path_txt.truncate() 
			explist_path_txt.write('conquer\n')
			explist_path_txt.write('wildtype\n')
			explist_path_txt.write('mutant\n')
			explist_path_txt.close()

			combolist_txt = open(combolist,"w+", newline ='\n')
			combolist_txt.truncate()
			combolist_txt.write('\n')
			# As it's the top result, we can just read the first line from successfulboards.txt 
			topline = open('OUTPUT_script_3/successfulboards.txt').readline()
			combolist_txt.write(topline)
			combolist_txt.close()

			endtimes_txt = open(endtimes,"w+", newline ='\n')
			endtimes_txt.truncate() 
			endtimes_txt.write('1\t11.22\nNon Essential Divided\n')
			endtimes_txt.write('2\t9.183\nNon Essential Divided\n')
			endtimes_txt.write('1\t13.89\nDNA/RNA/Protein/Metabolic NoDivision\n')
			endtimes_path_txt.close()

			print("This stage has finished. You should run the next stage.\nYou are carrying forward one variant (red), your results have been placed in INPUT_script_4X/comboresults.txt\n")
			deletionlog = "OUTPUT_final/deletionlog.txt"
			log = open(deletionlog,"a+", newline ='\n')
			log.write("\n\nThis stage has finished. You should run the next stage.\nYou are carrying forward one variant (red), your results have been placed in INPUT_script_4X/comboresults.txt\n")
			log.close()

			exit()
		
		elif boardcode == '90b':
			## If already removed 90% of the genes, going to create resulting files without running the code / simulations,
			## to spoof the next script into working.  
			explist_path = "OUTPUT_script_3/mineconquer_red_0_exp.list"
			combolist = "OUTPUT_script_3/red_0_combosandgenes.txt"
			endtimes = "INPUT_script_4X/conquerko_red_0_endtimes.txt"

			explist_path_txt = open(explist_path,"w+", newline ='\n')
			explist_path_txt.truncate() 
			explist_path_txt.write('conquer\n')
			explist_path_txt.write('wildtype\n')
			explist_path_txt.write('mutant\n')
			explist_path_txt.close()

			combolist_txt = open(combolist,"w+", newline ='\n')
			combolist_txt.truncate()
			combolist_txt.write('\n')
			# As it's the top result, we can just read the first line from successfulboards.txt 
			topline = open('OUTPUT_script_3/successfulboards.txt').readline()
			combolist_txt.write(topline)
			combolist_txt.close()

			endtimes_txt = open(endtimes,"w+", newline ='\n')
			endtimes_txt.truncate() 
			endtimes_txt.write('1\t11.22\nNon Essential Divided\n')
			endtimes_txt.write('2\t9.183\nNon Essential Divided\n')
			endtimes_txt.write('1\t13.89\nDNA/RNA/Protein/Metabolic NoDivision\n')
			endtimes_path_txt.close()

			print("This stage has finished. You should run the next stage.\nYou are carrying forward one variant (red), your results have been placed in INPUT_script_4X/comboresults.txt\n")
			deletionlog = "OUTPUT_final/deletionlog.txt"
			log = open(deletionlog,"a+", newline ='\n')
			log.write("\n\nThis stage has finished. You should run the next stage.\nYou are carrying forward one variant (red), your results have been placed in INPUT_script_4X/comboresults.txt\n")
			log.close()

			exit()
		
		### Logic of matching segments starts here
		elif boardcode == '80a':
			components = ['80a', '12.5g', '12.5h']
			outputToLists(components, boardnames, genecodes, red_yellow_blue)

		elif boardcode == '80b':
			components = ['80b', '12.5a', '12.5b']
			outputToLists(components, boardnames, genecodes, red_yellow_blue)

		elif boardcode == '70a':
			components = ['70a', '12.5f', '12.5g', '12.5h']
			outputToLists(components, boardnames, genecodes, red_yellow_blue)

		elif boardcode == '70b':
			components = ['70b', '12.5a', '12.5b', '12.5c']
			outputToLists(components, boardnames, genecodes, red_yellow_blue)

		elif boardcode == '60a':
			components = ['60a', '33c', '12.5e', '12.5f', '12.5g', '12.5h']
			outputToLists(components, boardnames, genecodes, red_yellow_blue)

		elif boardcode == '60b':
			components = ['60b', '33a', '12.5a', '12.5b', '12.5c', '12.5d']
			outputToLists(components, boardnames, genecodes, red_yellow_blue)

		elif boardcode == '50a':
			components = ['50a', '33c', '12.5e', '12.5f', '12.5g', '12.5h']
			outputToLists(components, boardnames, genecodes, red_yellow_blue)

		elif boardcode == '50b':
			components = ['50b', '33a', '12.5a', '12.5b', '12.5c', '12.5d']
			outputToLists(components, boardnames, genecodes, red_yellow_blue)

		elif boardcode == '33a':
			components = ['33a', '33b', '33c', '12.5d', '12.5e', '12.5f', '12.5g', '12.5h']
			outputToLists(components, boardnames, genecodes, red_yellow_blue)

		elif boardcode == '33b':
			components = ['33b', '33a', '33c', '12.5a', '12.5b', '12.5f', '12.5g', '12.5h']
			outputToLists(components, boardnames, genecodes, red_yellow_blue)

		elif boardcode == '33c':
			components = ['33c', '33a', '33b', '12.5a', '12.5b', '12.5c', '12.5d', '12.5e']
			outputToLists(components, boardnames, genecodes, red_yellow_blue)

		elif boardcode == '25a':
			components = ['25a', '12.5c', '12.5d', '12.5e', '12.5f', '12.5g', '12.5h']
			outputToLists(components, boardnames, genecodes, red_yellow_blue)

		elif boardcode == '25b':
			components = ['25b', '12.5a', '12.5b', '12.5e', '12.5f', '12.5g', '12.5h']
			outputToLists(components, boardnames, genecodes, red_yellow_blue)

		elif boardcode == '25c':
			components = ['25c', '12.5a', '12.5b', '12.5c', '12.5d', '12.5g', '12.5h']
			outputToLists(components, boardnames, genecodes, red_yellow_blue)

		elif boardcode == '25d':
			components = ['25d', '12.5a', '12.5b', '12.5c', '12.5d', '12.5e', '12.5f']
			outputToLists(components, boardnames, genecodes, red_yellow_blue)

		elif boardcode == '12.5a':
			components = ['12.5a', '12.5b', '12.5c', '12.5d', '12.5e', '12.5f', '12.5g', '12.5h']
			outputToLists(components, boardnames, genecodes, red_yellow_blue)

		elif boardcode == '12.5b':
			components = ['12.5a', '12.5b', '12.5c', '12.5d', '12.5e', '12.5f', '12.5g', '12.5h']
			outputToLists(components, boardnames, genecodes, red_yellow_blue)

		elif boardcode == '12.5c':
			components = ['12.5a', '12.5b', '12.5c', '12.5d', '12.5e', '12.5f', '12.5g', '12.5h']
			outputToLists(components, boardnames, genecodes, red_yellow_blue)

		elif boardcode == '12.5d':
			components = ['12.5a', '12.5b', '12.5c', '12.5d', '12.5e', '12.5f', '12.5g', '12.5h']
			outputToLists(components, boardnames, genecodes, red_yellow_blue)

		elif boardcode == '12.5e':
			components = ['12.5a', '12.5b', '12.5c', '12.5d', '12.5e', '12.5f', '12.5g', '12.5h']
			outputToLists(components, boardnames, genecodes, red_yellow_blue)

		elif boardcode == '12.5f':
			components = ['12.5a', '12.5b', '12.5c', '12.5d', '12.5e', '12.5f', '12.5g', '12.5h']
			outputToLists(components, boardnames, genecodes, red_yellow_blue)

		elif boardcode == '12.5g':
			components = ['12.5a', '12.5b', '12.5c', '12.5d', '12.5e', '12.5f', '12.5g', '12.5h']
			outputToLists(components, boardnames, genecodes, red_yellow_blue)

		elif boardcode == '12.5h':
			components = ['12.5a', '12.5b', '12.5c', '12.5d', '12.5e', '12.5f', '12.5g', '12.5h']
			outputToLists(components, boardnames, genecodes, red_yellow_blue)

def createScripts(job, simname):
	''' Convert template script using user input, create exp and ko txt files for Variants(.txt)'''

	names = ['red_0', 'yellow_0', 'blue_0']

	user_input_txt = "INPUT_script_1/user_input.txt"
	templatescript = "templatescript/TemplateScript.sh"
	variantcombos = "OUTPUT_script_3/{}_combosandgenes.txt"
	kolist = "OUTPUT_script_3/{}_ko.list"
	explist = "OUTPUT_script_3/{}_exp.list"
	experimentscript = "OUTPUT_script_3/{}.sh"

	for name in names:
		jobname = job + '_' + name
		
		#create knock out txt file 
		variantcombos_clone = variantcombos.format(name)
		kolist_clone = kolist.format(jobname) 
		ko = open(kolist_clone, 'w+', newline ='\n')
		kotempscript = open(variantcombos_clone).read()
		kotempscript = re.sub(r"^[^']*'", "'", kotempscript, flags=re.MULTILINE)
		# replace blank combo from beginning, just removed with regex
		kotempscript = " \n" + kotempscript
		ko.write(kotempscript)
		#append control simulations to gene candidates
		ko.write("%s\n" % "")
		ko.write("%s\n" % "'MG_006',")
		# use to create exp file of same length
		kocount = sum(1 for line in open(kolist_clone)) 

		#create exp txt file
		explist_clone = explist.format(jobname)
		with open(explist_clone, 'w+', newline ='\n') as exp:
			stopper = kocount
			for gene in range(0,stopper):
				exp.write("%s\n" % simname)
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
		ENDTIMENAME = 'conquerko' + '_' + name + '_e'
		endtimename_clone = ENDTIMENAME
		tempscript = re.sub(r"'e'", endtimename_clone, tempscript, flags=re.MULTILINE)

		### Calculate: DATE, SIMNUMBER
		### DATE
		DATE = "{:%Y-%m-%d}".format(datetime.now())
		tempscript = re.sub(r"DATE", DATE, tempscript, flags=re.MULTILINE)
		
		### SIMNUMBER
		SIMNUMBER = kocount + 2
		SIMNUMBER = str(SIMNUMBER)
		tempscript = re.sub(r"SIMNUMBER", SIMNUMBER, tempscript, flags=re.MULTILINE)

		expscript.write(tempscript)

		ko.close()
		exp.close()
		expscript.close()
		
	print('\nCreated bash scripts (*.sh), *_exp.lists, and *_ko.lists in OUTPUT_script_3.\n\n> See readme / lines 10-32 in bash script for what to do next.\n> See lines 128 - 154 in bash script for directories you need to create locally or on supercomputer.\n')
	print('\nOne of your expected folder structures is:\n')
	print(f'- projects\n\t- {GROUP_value}\n\t\t- {USER_value}\n\t\t\t- output\n\t\t\t\t- {projectfolder_value}\n\t\t\t\t\t- {job}\n\t\t\t\t\t\t- {simname}\n\t\t\t\t\t\t- wildtype\n\t\t\t\t\t\t- mutant\n\t\t\t\t\t\t- pdfs\n\t\t\t\t\t\t- figs\n')
	
	deletionlog = "OUTPUT_final/deletionlog.txt"
	log = open(deletionlog,"a+", newline ='\n')
	log.write("\nCreated bash scripts (*.sh), *_exp.lists, and *_ko.lists in OUTPUT_script_3.\n\n> See readme / lines 10-32 in bash script for what to do next.\n> See lines 128 - 154 in bash script for directories you need to create locally or on supercomputer.\n")
	log.write("\nOne of your expected folder structures is:\n")
	log.write(f'- projects\n\t- {GROUP_value}\n\t\t- {USER_value}\n\t\t\t- output\n\t\t\t\t- {projectfolder_value}\n\t\t\t\t\t- {job}\n\t\t\t\t\t\t- {simname}\n\t\t\t\t\t\t- wildtype\n\t\t\t\t\t\t- mutant\n\t\t\t\t\t\t- pdfs\n\t\t\t\t\t\t- figs\n\n')
	log.close()

def main():
	splashscreen()
	alreadyRunCheck()
	interpretResults()
	successlist = successCheckAndVariants()
	variantCombinations(successlist)
	createScripts(JOB, SIMNAME)

main()