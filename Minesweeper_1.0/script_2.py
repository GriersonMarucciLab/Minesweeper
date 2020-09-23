#!/usr/bin/env python

# Author: Joshua Rees, joshua.rees@bristol.ac.uk
# Affiliation: Life Sciences, University of Bristol 
# Last Updated: 2019-06-04

"""
Stage 2: 26 deletion segments, ranging in size from 100% to 12.5% of the low / no essentiality genes, are generated. 
Deletion segments that do not prevent division go to Stage 3. 
Create deletion segments from inputko*_endtimes.txt, output segments bash scripts and text files.
Expects: Place N inputko*.txt files in INPUT_script_2 folder.
Output: deletion segments, deletion segments list, single bash script and exp / ko lists, in OUTPUT_script_2 folder
"""

# Imports
import os
import fnmatch
import re
import math
from datetime import datetime

# Global variables
#JOB = "What is this experiment/run of simulations called? "
name = 'mine'
JOB = name + 'divide'
SIMNAMESTART = len(name)
SIMNAMEEND = (len(JOB))
SIMNAME = JOB[SIMNAMESTART:SIMNAMEEND]

def splashscreen():
    """ Fancy splash screen for Lego scripts """

    pass

def interpretResults():
    """ Match the simulation results of Stage 1 to the knocked out gene, create unsortedsimresults and matchedgeneresults.txt. """

    explist_path = "OUTPUT_script_1/mineinputko{}_exp.list"
    genelist = "OUTPUT_script_1/gene_list.txt"
    endtimes = "INPUT_script_2/inputko{}_endtimes.txt"
    unsortedsimresults = "OUTPUT_script_2/unsortedsimresults{}.txt"
    matchedresults = "OUTPUT_script_2/matchedgeneresults.txt"
    full_results_list = []

    numberofendtimes = len(fnmatch.filter(os.listdir('INPUT_script_2/'), '*endtimes.txt'))
    for n in range(1, numberofendtimes + 1):
        iteration = n
        nthendtimes = endtimes.format(iteration)

        # open endtimes results and solve results being on two lines
        firstinput = open(nthendtimes).read()
        firstinput = re.sub(r"(\d)\n^(\w)", r"\1\t\2", firstinput, flags=re.MULTILINE)

        nthunsortedsimresults = unsortedsimresults.format(iteration)
        firstoutput = open(nthunsortedsimresults,"w+", newline ='\n')
        firstoutput.truncate()
        firstoutput.write(firstinput)
        firstoutput.close()

        # open unsortedsimresults.txt, create list and a list of equal length, account for intended length, insert sim result in location, append lists together, save as file
        secondinput = open(nthunsortedsimresults, "r+")
        result_list = [line.rstrip() for line in secondinput.readlines()]

        # create an equal length list of default 'No_Result'
        result_list_copy = []
        for line in result_list:
            result_list_copy.append('No_Result')

        # add additional lines to account for crashed sims
        explist = explist_path.format(iteration)
        lenexplist = sum(1 for line in open(explist))
        lenendtimes = sum(1 for line in open(nthunsortedsimresults))
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
    thirdinput = open(genelist, "r+")
    gene_list = [line.rstrip() for line in thirdinput.readlines()]
    full_results_list = zip(gene_list, full_results_list)

    # save matchedresults in outputdir
    secondoutput = open(matchedresults,"w+", newline ='\n')
    secondoutput.truncate()
    # convert zip created list of tuples, line by line into string
    for line in full_results_list:
      secondoutput.write('\t'.join(str(s) for s in line) + '\n')
    secondoutput.close()

    print('\nHave matched the simulation results with gene codes and saved the results in OUTPUT_script_2/matchedgeneresults.txt.\n')

    deletionlog = "OUTPUT_final/deletionlog.txt"
    log = open(deletionlog,"a+", newline ='\n')
    log.write("\nHave matched the simulation results with gene codes and saved the results in OUTPUT_script_2/matchedgeneresults.txt.\n")
    log.close()

def dividedAndProducedProteinRNA(line):
    linetocheck = line
    divided = 'Divided'
    outcome = 'UP'
    
    if linetocheck.endswith(divided) and linetocheck.count(outcome) == 2:
        return True
    else:
        return False
    

def createNEList():
    """ Filter the Stage 1 matched results, finding those that divided, while avoiding excluded genes.
        Creates nonessential.txt """

    matchedresults = "OUTPUT_script_2/matchedgeneresults.txt"
    nonessential = "OUTPUT_script_2/nonessential.txt"
    exclusionlist = "INPUT_script_2/exclusionlist.txt"

    generesults = open(matchedresults, "r+")
    generesults_list = [line.rstrip() for line in generesults.readlines()]
    non_essential_genes = []

    for line in generesults_list:
        if dividedAndProducedProteinRNA(line) is True:
            non_essential_genes.append(line)

    neresults_gene_list = []
    # the tab split and line[0] element selection produces a single line of output in divisionsegment{}.txt
    for line in non_essential_genes:
         line = line.split("\t")
         line = line[0]
         line = line + ' '
         neresults_gene_list.append(line)

    if os.path.isfile(exclusionlist):
        exclusiongenes = open(exclusionlist, "r+")
        exclusiongenes_list = [line.rstrip() for line in exclusiongenes.readlines()]
        for gene in exclusiongenes_list:
            genelocation = exclusiongenes_list.index(gene)
            del exclusiongenes_list[genelocation]
            gene = gene + ' '
            exclusiongenes_list.insert(genelocation, gene)

        for excludedgene in exclusiongenes_list:
            genelocation = neresults_gene_list.index(excludedgene)
            del neresults_gene_list[genelocation]

    # save non_essential_genes in output dir
    neresults = open(nonessential,"w+", newline ='\n')
    neresults.truncate()
    for line in neresults_gene_list:
        line = line + '\n'
        neresults.write(line)
    neresults.close()

    print('\nHave filtered the results (deletions that still produce dividing cells) and saved in OUTPUT_script_2/nonessential.txt.\n')
    deletionlog = "OUTPUT_final/deletionlog.txt"
    log = open(deletionlog,"a+", newline ='\n')
    log.write("\nHave filtered the results (deletions that still produce dividing cells) and saved in OUTPUT_script_2/nonessential.txt.\n")
    log.close()

def segmentGeneration(segment_n, segment_s, neresults_gene_list):
    """ Daughter function of createDivisionSegments()
        used to select which genes to include in the 26 deletion segments """

    divisionseg = "OUTPUT_script_2/divisionsegment{}.txt"

    # 26 set segments
    whole = len(neresults_gene_list)
    hundred = 100
    zero = 0
    # math.trunc = round down > prevent overlapping gene knockouts (! may result in genes being missed)

    clonelist = neresults_gene_list

    ### percentage calculation (from segment_name which is passed into the function)
    percent = segment_n

    ngenes = (whole/hundred)*percent
    ngenes = math.trunc(ngenes)
    inversengenes = whole - ngenes
    print(f"\n{percent}% gene deletion = {ngenes} genes deleted out of {whole} non-essential genes, leaving {inversengenes} not deleted.")

    ### gene selection depending on percentage size and starting gene list from left or right (/top or bottom)
    if segment_s == 'a' and percent > 49:
        ### a (left anchored) > arbitary in the case of 100
        leftstart = zero
        leftend = whole - inversengenes
        divisionsegment_list = clonelist[leftstart:leftend]
    elif segment_s == 'b' and percent > 49:
        ### b (right anchored)
        rightstart = zero + inversengenes
        rightend = whole
        divisionsegment_list = clonelist[rightstart:rightend]
    elif segment_s == 'a' and percent < 49 :
        start = zero
        end = ngenes + 1
        divisionsegment_list = clonelist[start:end]
    elif segment_s == 'b' and percent < 49 :
        start = ngenes + 1
        end = ngenes*2 + 2
        divisionsegment_list = clonelist[start:end]
    elif segment_s == 'c' and percent == 33:
        start = ngenes*2 + 2
        end = whole
        divisionsegment_list = clonelist[start:end]
    elif segment_s == 'c' and percent < 33:
        start = ngenes*2 + 2
        end = ngenes*3 + 3
        divisionsegment_list = clonelist[start:end]
    elif segment_s == 'd' and percent == 25:
        start = ngenes*3 + 3
        end = whole
        divisionsegment_list = clonelist[start:end]
    elif segment_s == 'd' and percent == 12.5:
        start = ngenes*3 + 3
        end = ngenes*4 + 4
        divisionsegment_list = clonelist[start:end]
    elif segment_s == 'e' and percent == 12.5:
        start = ngenes*4 + 4
        end = ngenes*5 + 5
        divisionsegment_list = clonelist[start:end]
    elif segment_s == 'f' and percent == 12.5:
        start = ngenes*5 + 5
        end = ngenes*6 + 6
        divisionsegment_list = clonelist[start:end]
    elif segment_s == 'g' and percent == 12.5:
        start = ngenes*6 + 6
        end = ngenes*7 + 7
        divisionsegment_list = clonelist[start:end]
    elif segment_s == 'h' and percent == 12.5:
        start = ngenes*7 + 7
        end = whole
        divisionsegment_list = clonelist[start:end]

    ### output to own text file
    if percent == 100:
        divisionsegment_name = str(segment_n)
        divisionsegment = divisionseg.format(divisionsegment_name)
    elif percent == 12.5:
        divisionsegment_name = '12_5' + segment_s
        divisionsegment = divisionseg.format(divisionsegment_name)
    else:
        divisionsegment_name = str(segment_n) + segment_s
        divisionsegment = divisionseg.format(divisionsegment_name)

    output_divisionsegment = open(divisionsegment,"w+", newline ='\n')
    output_divisionsegment.truncate()
    for line in divisionsegment_list:
        output_divisionsegment.write(line)
    output_divisionsegment.close()

    print(f'\nHave created the {percent}% deletion segment in OUTPUT_script_2/divisionsegment{divisionsegment_name}.txt')

    deletionlog = "OUTPUT_final/deletionlog.txt"
    log = open(deletionlog,"a+", newline ='\n')
    log.write(f"\nHave created the {percent}% deletion segment in OUTPUT_script_2/divisionsegment{divisionsegment_name}.txt\n")
    log.close()

    return divisionsegment_list

def outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, trigger):
    """ Daughter function of createDivisionSegments()
        used to save the 26 deletion segments as seperate files,
        and a combined file on the final segment / YES trigger """

    alldivisionsegments_txt = "OUTPUT_script_2/alldivisionsegments.txt"
    alldivisionsegments_codes_txt = "OUTPUT_script_2/alldivisionsegments_codes.txt"

    if segment_n == 100:
        ### output to alldivisionsegments_list and codes
        alldivisionsegments_codes.append(str(segment_n) + '\n')
        alldivisionsegments_list.append(str(segment_n) + '\n')
        for line in temp_list:
            alldivisionsegments_list.append(line)
        alldivisionsegments_list.append('\n')
    else:
        ### output to alldivisionsegments_list and codes
        alldivisionsegments_codes.append(str(segment_n) + segment_s + '\n')
        alldivisionsegments_list.append(str(segment_n) + segment_s + '\n')
        for line in temp_list:
            alldivisionsegments_list.append(line)
        alldivisionsegments_list.append('\n')

    if trigger == 'Yes':
        # OUTPUT alldivisionsegments_list to text
        output_alltxt = open(alldivisionsegments_txt,"w+", newline ='\n')
        output_alltxt.truncate()
        for line in alldivisionsegments_list:
            output_alltxt.write(line)
        output_alltxt.close()

        # OUTPUT alldivisionsegments_codes to text
        output_alltxt2 = open(alldivisionsegments_codes_txt,"w+", newline ='\n')
        output_alltxt2.truncate()
        for line in alldivisionsegments_codes:
            output_alltxt2.write(line)
        output_alltxt2.close()

        deletionlog = "OUTPUT_final/deletionlog.txt"
        log = open(deletionlog,"a+", newline ='\n')
        log.write("\nHave created the 26 deletion segments, see OUTPUT_script_2/alldivisionsegments.txt\n")
        log.close()

    else:
        pass

    return alldivisionsegments_codes, alldivisionsegments_list

def createDivisionSegments():
    """ Generates the 26 deletion segments """

    nonessential = "OUTPUT_script_2/nonessential.txt"
    alldivisionsegments_list = []
    alldivisionsegments_codes = []
    # divsionsegments_txt == minedivide_exp.list

    neresults = open(nonessential, "r+")
    neresults_list = [line.rstrip() for line in neresults.readlines()]
    neresults_gene_list = []
    for line in neresults_list:
         line = line + ' '
         neresults_gene_list.append(line)

    ### The 26 Segments

    ## 100
    segment_n = 100
    segment_s = 'a'
    temp_list = segmentGeneration(segment_n, segment_s, neresults_gene_list)
    alldivisionsegments_codes, alldivisionsegments_list = outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, 'No')

    # 90a
    segment_n = 90
    segment_s = 'a'
    temp_list = segmentGeneration(segment_n, segment_s, neresults_gene_list)
    alldivisionsegments_codes, alldivisionsegments_list = outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, 'No')
    # 90b
    segment_n = 90
    segment_s = 'b'
    temp_list = segmentGeneration(segment_n, segment_s, neresults_gene_list)
    alldivisionsegments_codes, alldivisionsegments_list = outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, 'No')

    # 80a
    segment_n = 80
    segment_s = 'a'
    temp_list = segmentGeneration(segment_n, segment_s, neresults_gene_list)
    alldivisionsegments_codes, alldivisionsegments_list = outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, 'No')
    # 80b
    segment_n = 80
    segment_s = 'b'
    temp_list = segmentGeneration(segment_n, segment_s, neresults_gene_list)
    alldivisionsegments_codes, alldivisionsegments_list = outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, 'No')

    # 70a
    segment_n = 70
    segment_s = 'a'
    temp_list = segmentGeneration(segment_n, segment_s, neresults_gene_list)
    alldivisionsegments_codes, alldivisionsegments_list = outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, 'No')
    # 70b
    segment_n = 70
    segment_s = 'b'
    temp_list = segmentGeneration(segment_n, segment_s, neresults_gene_list)
    alldivisionsegments_codes, alldivisionsegments_list = outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, 'No')

    # 60a
    segment_n = 60
    segment_s = 'a'
    temp_list = segmentGeneration(segment_n, segment_s, neresults_gene_list)
    alldivisionsegments_codes, alldivisionsegments_list = outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, 'No')
    # 60b
    segment_n = 60
    segment_s = 'b'
    temp_list = segmentGeneration(segment_n, segment_s, neresults_gene_list)
    alldivisionsegments_codes, alldivisionsegments_list = outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, 'No')

    # 50a
    segment_n = 50
    segment_s = 'a'
    temp_list = segmentGeneration(segment_n, segment_s, neresults_gene_list)
    alldivisionsegments_codes, alldivisionsegments_list = outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, 'No')
    # 50b
    segment_n = 50
    segment_s = 'b'
    temp_list = segmentGeneration(segment_n, segment_s, neresults_gene_list)
    alldivisionsegments_codes, alldivisionsegments_list = outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, 'No')

    # 33a
    segment_n = 33
    segment_s = 'a'
    temp_list = segmentGeneration(segment_n, segment_s, neresults_gene_list)
    alldivisionsegments_codes, alldivisionsegments_list = outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, 'No')
    # 33b
    segment_n = 33
    segment_s = 'b'
    temp_list = segmentGeneration(segment_n, segment_s, neresults_gene_list)
    alldivisionsegments_codes, alldivisionsegments_list = outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, 'No')
    # 33c
    segment_n = 33
    segment_s = 'c'
    temp_list = segmentGeneration(segment_n, segment_s, neresults_gene_list)
    alldivisionsegments_codes, alldivisionsegments_list = outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, 'No')

    # 25a
    segment_n = 25
    segment_s = 'a'
    temp_list = segmentGeneration(segment_n, segment_s, neresults_gene_list)
    alldivisionsegments_codes, alldivisionsegments_list = outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, 'No')
    # 25b
    segment_n = 25
    segment_s = 'b'
    temp_list = segmentGeneration(segment_n, segment_s, neresults_gene_list)
    alldivisionsegments_codes, alldivisionsegments_list = outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, 'No')
    # 25c
    segment_n = 25
    segment_s = 'c'
    temp_list = segmentGeneration(segment_n, segment_s, neresults_gene_list)
    alldivisionsegments_codes, alldivisionsegments_list = outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, 'No')
    # 25d
    segment_n = 25
    segment_s = 'd'
    temp_list = segmentGeneration(segment_n, segment_s, neresults_gene_list)
    alldivisionsegments_codes, alldivisionsegments_list = outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, 'No')

    # 12_5a
    segment_n = 12.5
    segment_s = 'a'
    temp_list = segmentGeneration(segment_n, segment_s, neresults_gene_list)
    alldivisionsegments_codes, alldivisionsegments_list = outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, 'No')
    # 12_5b
    segment_n = 12.5
    segment_s = 'b'
    temp_list = segmentGeneration(segment_n, segment_s, neresults_gene_list)
    alldivisionsegments_codes, alldivisionsegments_list = outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, 'No')
    # 12_5c
    segment_n = 12.5
    segment_s = 'c'
    temp_list = segmentGeneration(segment_n, segment_s, neresults_gene_list)
    alldivisionsegments_codes, alldivisionsegments_list = outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, 'No')
    # 12_5d
    segment_n = 12.5
    segment_s = 'd'
    temp_list = segmentGeneration(segment_n, segment_s, neresults_gene_list)
    alldivisionsegments_codes, alldivisionsegments_list = outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, 'No')
    # 12_5e
    segment_n = 12.5
    segment_s = 'e'
    temp_list = segmentGeneration(segment_n, segment_s, neresults_gene_list)
    alldivisionsegments_codes, alldivisionsegments_list = outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, 'No')
    # 12_5f
    segment_n = 12.5
    segment_s = 'f'
    temp_list = segmentGeneration(segment_n, segment_s, neresults_gene_list)
    alldivisionsegments_codes, alldivisionsegments_list = outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, 'No')
    # 12_5g
    segment_n = 12.5
    segment_s = 'g'
    temp_list = segmentGeneration(segment_n, segment_s, neresults_gene_list)
    alldivisionsegments_codes, alldivisionsegments_list = outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, 'No')
    # 12_5h
    segment_n = 12.5
    segment_s = 'h'
    temp_list = segmentGeneration(segment_n, segment_s, neresults_gene_list)
    alldivisionsegments_codes, alldivisionsegments_list = outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, 'Yes')

def createScripts(job, simname):
    """ Convert template script to bash script using user input, create exp and ko txt files """

    user_input_txt = "INPUT_script_1/user_input.txt"
    templatescript = "templatescript/TemplateScript.sh"
    divisions = "OUTPUT_script_2/alldivisionsegments.txt"
    kolist = "OUTPUT_script_2/{}_ko.list"
    explist = "OUTPUT_script_2/{}_exp.list"
    experimentscript = "OUTPUT_script_2/{}.sh"

    #create knock out txt file using the combined 26 segment file

    kolist_clone = kolist.format(job)
    ko = open(kolist_clone, 'w+', newline ='\n')
    kotempscript = open(divisions).read()
    kotempscript = re.sub(r"\d\d\w\n", "", kotempscript, flags=re.MULTILINE)
    kotempscript = re.sub(r"\d\d\.\d\w\n", "", kotempscript, flags=re.MULTILINE)
    kotempscript = re.sub(r"\d\d\_\d\n", "", kotempscript, flags=re.MULTILINE)
    kotempscript = re.sub(r"\d\d\_\d\d\n", "", kotempscript, flags=re.MULTILINE)
    kotempscript = re.sub(r"\d\d\_\d\d\d\n", "", kotempscript, flags=re.MULTILINE)
    ko.write(kotempscript)
    #append control simulations to gene candidates
    ko.write("%s\n" % "")
    ko.write("%s\n" % "'MG_006',")

    kocount = sum(1 for line in open(kolist_clone))

    #create exp txt file

    explist_clone = explist.format(job)
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
    experimentscript_clone = experimentscript.format(job)
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
    job_clone = JOB
    tempscript = re.sub(r"JOBx", job_clone, tempscript, flags=re.MULTILINE)

    ### ENDTIMENAME default is 'e', producing 'e' + 'ndtimes.txt'
    ENDTIMENAME = 'divideko_e'
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

    print('\nCreated bash script (*.sh), *_exp.list, and *_ko.list in OUTPUT_script_2.\n\n> See readme / lines 10-32 in bash script for what to do next.\n> See lines 128 - 154 in bash script for directories you need to create locally or on supercomputer.')
    print('\nYour expected folder structure is:')
    print(f'- projects\n\t- {GROUP_value}\n\t\t- {USER_value}\n\t\t\t- output\n\t\t\t\t- {projectfolder_value}\n\t\t\t\t\t- {job}\n\t\t\t\t\t\t- {simname}\n\t\t\t\t\t\t- wildtype\n\t\t\t\t\t\t- mutant\n\t\t\t\t\t\t- pdfs\n\t\t\t\t\t\t- figs')

    deletionlog = "OUTPUT_final/deletionlog.txt"
    log = open(deletionlog,"a+", newline ='\n')
    log.write("\nCreated bash script (*.sh), *_exp.list, and *_ko.list in OUTPUT_script_2.\n\n> See readme / lines 10-32 in bash script for what to do next.\n> See lines 128 - 154 in bash script for directories you need to create locally or on supercomputer.\n")
    log.write("\nYour expected folder structure is:\n")
    log.write(f'- projects\n\t- {GROUP_value}\n\t\t- {USER_value}\n\t\t\t- output\n\t\t\t\t- {projectfolder_value}\n\t\t\t\t\t- {job}\n\t\t\t\t\t\t- {simname}\n\t\t\t\t\t\t- wildtype\n\t\t\t\t\t\t- mutant\n\t\t\t\t\t\t- pdfs\n\t\t\t\t\t\t- figs\n\n')
    log.close()

def main():
    splashscreen()
    interpretResults()
    createNEList()
    createDivisionSegments()
    createScripts(JOB, SIMNAME)

main()