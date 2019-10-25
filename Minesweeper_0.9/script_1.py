#!/usr/bin/env python

# Author: Joshua Rees, joshua.rees@bristol.ac.uk
# Affiliation: Life Sciences, University of Bristol 
# Last Updated: 2019-06-04

"""
Stage 1: single gene knockouts are conducted to identify low / no essentiality genes (whose knockout does not prevent cell division) for Stage 2. 
Creates single gene knock out simulation bash scripts / exp list / ko list from list of genes.
Expects: Geme file in INPUT_script_1 FOLDER in format MG_XXX, one per line
Output: bash scripts and / exp list / ko list and gene list text into OUTPUT_script_1
"""


# Imports
import re
from datetime import datetime
import math
import os

# Global variables
#JOB = "What is this experiment/run of simulations called? "
name = 'mine'
JOB = name + 'inputko{}'
SIMNAMESTART = len(name)
SIMNAMEEND = (len(JOB))-2
SIMNAME = JOB[SIMNAMESTART:SIMNAMEEND]

def splashscreen():
	""" Fancy splash screen for Minesweeper scripts """

	pass

def userInput(job, simname):
	""" Ask user for variables to modify template script """

	JOB = job
	SIMNAME = simname
	user_input_txt = "INPUT_script_1/user_input.txt"
	linecount = 0

	if os.path.isfile(user_input_txt):
		linecount = len(open(user_input_txt).readlines(  ))

	if linecount >= 10:
		print('\nYou have pre-existing values in INPUT_script_1/user_input.txt - you can update them there if needed.')
	else:
		input("Have you read the readme? " )
		input("Have you placed genes.txt (your list of candidate gene deletions) in INPUT_script_1 folder? " )
		input("The data you are about to enter is only stored locally. Press Enter to continue." )
		NAME = input("What is your name? ")
		EMAIL = input("What is your email? ")
		DEPARTMENT = input("What is your department? ")
		UNI = input("What is your uni / institute? ")
		SUPERPC = input("What computer/supercomputer will you be using? ")
		SUPERQUEUE = input("Will you be using a queuing system or running locally? ")
		print('\nIf you are not using the SLURM queuing system you will have to modify /templatescript/TemplateScript.sh.')
		print('\nBelow is the expected folder structure for the simulation output.')
		print('- projects\n\t- GROUP\n\t\t- USER\n\t\t\t- output\n\t\t\t\t- PROJECTFOLDER\n\t\t\t\t\t- JOB\n\t\t\t\t\t\t- inputko\n\t\t\t\t\t\t- wildtype\n\t\t\t\t\t\t- mutant\n\t\t\t\t\t\t- pdfs\n\t\t\t\t\t\t- figs')
		print('\nIf you are not using this folder structure you will have to modify /templatescript/TemplateScript.sh.')
		print('\nThe following questions ask you for the names that are CAPITALISED. Please give answers in lowercase :)\n')
		GROUP = input("What is your lab group's name / project folder on the supercomputer? ")
		USER = input("What is your USER? ")
		PROJECTFOLDER = input("What is this research project called? ")
		JOB = JOB.format("")
		print(f'\nThis experiment/run of simulations (referred to as JOB) is {JOB}.\n')
		print('\nYour expected folder structure is:')
		print(f'- projects\n\t- {GROUP}\n\t\t- {USER}\n\t\t\t- output\n\t\t\t\t- {PROJECTFOLDER}\n\t\t\t\t\t- {JOB}\n\t\t\t\t\t\t- {SIMNAME}\n\t\t\t\t\t\t- wildtype\n\t\t\t\t\t\t- mutant\n\t\t\t\t\t\t- pdfs\n\t\t\t\t\t\t- figs')
		
		target = open(user_input_txt, 'w', newline ='\n')
		target.truncate()

		line = "{}\t{}\n"                      
		NAME_line = line.format("NAME", NAME)
		target.write(NAME_line)

		EMAIL_line = line.format("EMAIL", EMAIL)
		target.write(EMAIL_line)

		DEPARTMENT_line = line.format("DEPARTMENT", DEPARTMENT)
		target.write(DEPARTMENT_line)

		UNI_line = line.format("UNI", UNI)
		target.write(UNI_line)

		SUPERPC_line = line.format("SUPERPC", SUPERPC)
		target.write(SUPERPC_line)

		SUPERQUEUE_line = line.format("SUPERQUEUE", SUPERQUEUE)
		target.write(SUPERQUEUE_line)	

		GROUP_line = line.format("GROUP", GROUP)
		target.write(GROUP_line)	

		USER_line = line.format("USER", USER)
		target.write(USER_line)	

		PROJECTFOLDER_line = line.format("PROJECTFOLDER", PROJECTFOLDER)
		target.write(PROJECTFOLDER_line)	

		JOB_line = line.format("JOB", JOB)
		target.write(JOB_line)

		target.write('\nYour expected folder structure is:')
		target.write(f'- projects\n\t- {GROUP}\n\t\t- {USER}\n\t\t\t- output\n\t\t\t\t- {PROJECTFOLDER}\n\t\t\t\t\t- {JOB}\n\t\t\t\t\t\t- inputko\n\t\t\t\t\t\t- wildtype\n\t\t\t\t\t\t- mutant\n\t\t\t\t\t\t- pdfs\n\t\t\t\t\t\t- figs')
		
		target.close()

		print('\nYour input values have been saved in INPUT_script_1/user_input.txt - you can update them there if needed.')
		
		deletionlog = "OUTPUT_final/deletionlog.txt"
		log = open(deletionlog,"a+", newline ='\n')
		log.write("\nYour input values have been saved in INPUT_script_1/user_input.txt - you can update them there if needed.\n")
		log.close()

def createGeneList():
	""" Convert given txt file into correct gene format, create a usable list
		Used by createScripts() """

	given_genes_txt = "INPUT_script_1/genes.txt"
	gene_list_txt = "OUTPUT_script_1/gene_list.txt"

	# convert to format 'MG_XXX' and save as text file in output folder
	# re.sub(pattern, repl, string, count=0, flags=0)
	# re.sub(r'^(\w\w\w\w\w\w)$', r''\1',', NAMEOFTEXTFILE, flags=re.MULTILINE)
	target = open(gene_list_txt,"w")
	target.truncate()
	input = open(given_genes_txt).read()
	target.write(re.sub(r"^(\w\w\w\w\w\w)$", r"'\1',", input, flags=re.MULTILINE))
	target.close()
	print('\nCreated gene_list.txt in OUTPUT_script_1.')

	aslist = open(gene_list_txt, "r+")
	gene_list = [line.rstrip() for line in aslist.readlines()]
	aslist.close()

	return gene_list
	# Access with LISTNAME[N] > N is 0 indexed (i.e. 0 = first entry in list)

def createScripts(gene_list,job, simname):
	""" Convert template script to bash script using user input, create exp and ko txt files """
	# Two pathways depending on number of genes: <200 - one script, >200 - multiple scripts

	JOB = job
	SIMNAME = simname

	user_input_txt = "INPUT_script_1/user_input.txt"
	templatescript = "templatescript/TemplateScript.sh"
	kolist = "OUTPUT_script_1/mineinputko{}_ko.list"
	explist = "OUTPUT_script_1/mineinputko{}_exp.list"
	experimentscript = "OUTPUT_script_1/mineinputko{}.sh"
	gene_list_retrieved = gene_list

	nofgenecandidates = len(gene_list_retrieved)

	if nofgenecandidates < 200:
		
		#create knock out txt file

		kolist = kolist.format("")
		with open(kolist, 'w+', newline ='\n') as ko:
			for gene in gene_list_retrieved:
				ko.write("%s\n" % gene)
			#append control simulations to gene candidates
			ko.write("%s\n" % "")
			ko.write("%s\n" % "'MG_006',")
		
		
		#create exp txt file

		explist = explist.format("")
		with open(explist, 'w+', newline ='\n') as exp:
			stopper = nofgenecandidates
			for gene in range(0,stopper):
				exp.write("%s\n" % SIMNAME)
			#append control simulations to gene candidates
			exp.write("%s\n" % "wildtype")
			exp.write("%s\n" % "mutant")	

		
		#create simulation bash script

		user_input_values = open(user_input_txt, "r+")
		tempscript = open(templatescript).read()
		experimentscript = experimentscript.format("")
		expscript = open(experimentscript, 'w+', newline ='\n')

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
		job_clone = JOB.format("")
		tempscript = re.sub(r"JOBx", job_clone, tempscript, flags=re.MULTILINE)
		
		### ENDTIMENAME default is 'e', producing 'e' + 'ndtimes.txt'
		ENDTIMENAME = 'inputko_e'
		endtimename_clone = ENDTIMENAME
		tempscript = re.sub(r"'e'", endtimename_clone, tempscript, flags=re.MULTILINE)

		### Calculate: DATE, SIMNUMBER
		### DATE
		DATE = "{:%Y-%m-%d}".format(datetime.now())
		tempscript = re.sub(r"DATE", DATE, tempscript, flags=re.MULTILINE)
		
		### SIMNUMBER
		SIMNUMBER = nofgenecandidates + 2
		SIMNUMBER = str(SIMNUMBER)
		tempscript = re.sub(r"SIMNUMBER", SIMNUMBER, tempscript, flags=re.MULTILINE)

		expscript.write(tempscript)

		ko.close()
		exp.close()
		expscript.close()
		
		print('\nCreated bash script (*.sh), *_exp.list, and *_ko.list in OUTPUT_script_1.\n\n> See readme / lines 10-32 in bash script for what to do next.\n> See lines 128 - 154 in bash script for directories you need to create locally or on supercomputer.')


	elif nofgenecandidates > 200:
		# division based on 200 simulations (+2 for controls), 
		# need X files which contain N simulations
		# divide into X new txt files > ko list, exp list in loop
		# modify X generic scripts 

		division = nofgenecandidates / 200
		divisions = math.ceil(division)
		penultimatedivision = divisions - 1

		divisioncounter = 1
		startinggene = 0
		scriptnameincrement = 1
		
		while divisioncounter < divisions:
			endinggene = startinggene + 200
			part_gene_list = gene_list_retrieved[startinggene:endinggene]
			scriptnameincrement_str = str(scriptnameincrement)
		
			#create knock out txt file
			tempko = kolist.format(scriptnameincrement_str)
			with open(tempko, 'w+', newline ='\n') as ko:
				for gene in part_gene_list:
					ko.write("%s\n" % gene)
				ko.write("%s\n" % "")
				ko.write("%s\n" % "'MG_006',")

			#create exp txt file
			tempexp = explist.format(scriptnameincrement_str)
			with open(tempexp, 'w+', newline ='\n') as exp:
				for gene in range(startinggene,endinggene):
					exp.write("%s\n" % SIMNAME)
				#append control simulations to gene candidates
				exp.write("%s\n" % "wildtype")
				exp.write("%s\n" % "mutant")

			#create simulation bash script
			user_input_values = open(user_input_txt, "r+")
			tempscript = open(templatescript).read()

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
			job_clone = JOB.format(scriptnameincrement_str)
			tempscript = re.sub(r"JOBx", job_clone, tempscript, flags=re.MULTILINE)

			### ENDTIMENAME default is 'e', producing 'e' + 'ndtimes.txt'
			ENDTIMENAME = 'inputko{}_e'
			endtimename_clone = ENDTIMENAME.format(scriptnameincrement_str)
			tempscript = re.sub(r"'e'", endtimename_clone, tempscript, flags=re.MULTILINE)

			### Calculate: DATE, SIMNUMBER
			### DATE
			DATE = "{:%Y-%m-%d}".format(datetime.now())
			tempscript = re.sub(r"DATE", DATE, tempscript, flags=re.MULTILINE)
			
			### SIMNUMBER
			SIMNUMBER = 202
			SIMNUMBER = str(SIMNUMBER)
			tempscript = re.sub(r"SIMNUMBER", SIMNUMBER, tempscript, flags=re.MULTILINE)

			tempexpscript = experimentscript.format(scriptnameincrement_str)
			expscript = open(tempexpscript, 'w+', newline ='\n')
			expscript.write(tempscript)

			startinggene = startinggene + 200
			scriptnameincrement = scriptnameincrement + 1
			divisioncounter = divisioncounter + 1
			ko.close()
			exp.close()
			expscript.close()

		while divisioncounter == divisions:
			endinggene = nofgenecandidates
			part_gene_list = gene_list_retrieved[startinggene:endinggene]
			scriptnameincrement_str = str(scriptnameincrement)

			#create knock out txt file
			tempko = kolist.format(scriptnameincrement_str)
			with open(tempko, 'w+', newline ='\n') as ko:
				for gene in part_gene_list:
					ko.write("%s\n" % gene)
				ko.write("%s\n" % "")
				ko.write("%s\n" % "'MG_006',")

			#create exp txt file
			tempexp = explist.format(scriptnameincrement_str)
			with open(tempexp, 'w+', newline ='\n') as exp:
				for gene in range(startinggene,endinggene):
					exp.write("%s\n" % SIMNAME)
				#append control simulations to gene candidates
				exp.write("%s\n" % "wildtype")
				exp.write("%s\n" % "mutant")

			#create simulation bash script
			user_input_values = open(user_input_txt, "r+")
			tempscript = open(templatescript).read()

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
			job_clone = JOB.format(scriptnameincrement_str)
			tempscript = re.sub(r"JOBx", job_clone, tempscript, flags=re.MULTILINE)

			### ENDTIMENAME default is 'e', producing 'e' + 'ndtimes.txt'
			ENDTIMENAME = 'inputko{}_e'
			endtimename_clone = ENDTIMENAME.format(scriptnameincrement_str)
			tempscript = re.sub(r"'e'", endtimename_clone, tempscript, flags=re.MULTILINE)

			### Calculate: DATE, SIMNUMBER
			### DATE
			DATE = "{:%Y-%m-%d}".format(datetime.now())
			tempscript = re.sub(r"DATE", DATE, tempscript, flags=re.MULTILINE)
			
			### SIMNUMBER
			SIMNUMBER = (nofgenecandidates - (penultimatedivision*200))+2
			SIMNUMBER = str(SIMNUMBER)
			tempscript = re.sub(r"SIMNUMBER", SIMNUMBER, tempscript, flags=re.MULTILINE)

			tempexpscript = experimentscript.format(scriptnameincrement_str)
			expscript = open(tempexpscript, 'w+', newline ='\n')
			expscript.write(tempscript)

			ko.close()
			exp.close()
			expscript.close()

			break

		print('\nCreated multiple bash scripts (*.sh), *_exp.lists, and *_ko.lists in OUTPUT_script_1.\n\n> See readme / lines 10-32 in bash script for what to do next.\n> See lines 128 - 154 in bash script for directories you need to create locally or on supercomputer.')

	deletionlog = "OUTPUT_final/deletionlog.txt"
	log = open(deletionlog,"a+", newline ='\n')
	log.write("\nCreated bash script/s (*.sh), *_exp.list/s, and *_ko.list/s in OUTPUT_script_1.\n\n> See readme / lines 10-32 in bash script for what to do next.\n> See lines 128 - 154 in bash script for directories you need to create locally or on supercomputer.\n")
	log.write("\nYour/one of your expected folder structures is:\n")
	log.write(f'- projects\n\t- {GROUP_value}\n\t\t- {USER_value}\n\t\t\t- output\n\t\t\t\t- {projectfolder_value}\n\t\t\t\t\t- {job}\n\t\t\t\t\t\t- {simname}\n\t\t\t\t\t\t- wildtype\n\t\t\t\t\t\t- mutant\n\t\t\t\t\t\t- pdfs\n\t\t\t\t\t\t- figs\n\n')
	log.close()

def main():
	splashscreen()
	userInput(JOB, SIMNAME)
	gene_list = createGeneList()
	createScripts(gene_list, JOB, SIMNAME)

main()