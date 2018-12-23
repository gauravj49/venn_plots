#!/usr/local/bin/python
"""
***********************************************
- PROGRAM: venn_plots.py
- CONTACT: Gaurav Jain(gauravj49@gmail.com)
***********************************************
"""
print (__doc__)

# Built in modules
import argparse
import os.path
import sys

# 3rd party modules
import textwrap
import re
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib as mp
#mp.use('Agg') # to use matplotlib without X11
import matplotlib.pyplot as plt
import subprocess
import binascii as bi
import scipy.stats as stats
from collections import *
#from scipy.stats.stats import nanmean
from numpy import nanmean

# for looping files in a dir
import glob

# user defined modules
from gjainLIB import *      # import all the functions from the Gaurav`s python library

### for color scale
from  matplotlib import colors
from itertools import cycle, islice # barplot colors

# Venn diagram for sets (2,3)
from itertools import combinations, chain
from matplotlib_venn import *

### for 5-set venn diagrams
import PIL
from PIL import ImageFont
from PIL import Image
from PIL import ImageDraw
from collections import defaultdict, OrderedDict

################ USER CONFIGURATION ###################
# 1)  Save matplotlib figure with text as vectorized text and not as a path
#     The matplotlib documentation in http://matplotlib.org/users/customizing.html states that the default value for the parameter svg.fonttype is 'path', which means that characters will be converted to paths when exporting to svg format.
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['pdf.fonttype'] = 42

# 2) Set print precision option
np.set_printoptions(precision=6)
#######################################################

def main():
    # Get input options
    args = check_options()

    # Store the variables
    input_file  = args.input_file
    output_file = args.output_file
    column_nums = args.column_nums
    output_dir  = args.output_dir
    caseIns     = args.caseIns

    if output_dir:
        output_dir = output_dir.rstrip('\/')

    # Get the column numbers to be plotted
    if not column_nums:
        with open(input_file, 'rU') as f: column_nums = ",".join([str(i) for i in xrange(int(len(f.readline().strip().split("\t"))))])
    column_numbers = map(int, column_nums.split(","))

    # Get the output file if not passed
    if not output_file:
        if output_dir:
            create_dir(output_dir)
        else:
            output_dir = get_file_info(input_file)[0]
        output_file = input_file
    else:
        output_dir = get_file_info(output_file)[0]

    # Rename output file to have column numbers added to the filename
    output_file = "{0}/{1}_columns_{2}.png".format(output_dir, get_file_info(output_file)[1],"_".join(map(str, column_numbers))) 

    # Get the data structure to pass it forward
    named_sets = defaultdict(set)
    group_list = list()
    named_sets, group_list = get_interactions_dict(input_file, column_numbers, caseIns)
    mtitle = "Number of unique and overlapping elements in :\n" + " vs. ".join(sorted(group_list))
    
    # Get the output peaks_log file
    profile_file = get_file_info(output_file)[0] + "/" + get_file_info(output_file)[1] + '_profile.txt'
    olf      = open(profile_file, 'w')
    summary_file = get_file_info(output_file)[0] + "/" + get_file_info(output_file)[1] + '_summary.txt'
    slf      = open(summary_file, 'w')

    # Create all possible combination of the groups and store the respective number peaks for each
    all_elements = defaultdict()
    peak_count = defaultdict()

    print "- Overlap:"
    for intersected, unioned, count, peaks_list in venn_count(named_sets):
        # Get the groups
        group = '|'.join(sorted(intersected))
        print "\t", group," = ", str(count)
        slf.write(group + " = " + str(count) + "\n")
        
        # Get the peaks for each groups
        plist = "\n".join([i for i in sorted(peaks_list)])
        #print str("> " + group + " = " + str(count) + "\n" + plist + "\n")

        # Write it to the text file
        olf.write(str("> " + group + " = " + str(count) +  "\n" + plist + "\n"))
        
        # Store it to draw the venn diagrams
        peak_count[group] = str(count)

        # Store the list to use it in the tabular format
        # all_elements['mir_191'] = 'CA1|Blood'
        for i in peaks_list:
            all_elements[i] = group
    olf.close()

    # draw the venn diagram and save the data in summary and profile(fasta) format
    if len(group_list) == 2:
        generate_venn2(group_list, peak_count, output_file, named_sets, 2, mtitle)
    elif len(group_list) == 3:
        generate_venn3(group_list, peak_count, output_file, named_sets, 3, mtitle)
    elif len(group_list) == 4:
        generate_venn4(group_list, peak_count, output_file, named_sets, mtitle)
    elif len(group_list) == 5:
        generate_venn5(group_list, peak_count, output_file, named_sets, mtitle)
        #Generate_venn5(group_list, peak_count, output_file, named_sets)

    # Save the data in the tabular format
    toutput_file = get_file_info(output_file)[0] + "/" + get_file_info(output_file)[1] + '_tabular.txt'
    save_tabular(all_elements, sorted(group_list), toutput_file)

    # Print output file locations
    print "- Your output profile file is: " + profile_file
    print "- Your output summary file is: " + summary_file


################ USER DEFINED FUNCTIONS ###################
def generate_venn2(group_list, peak_count, output_file, named_sets, ngroups, mtitle):
    # The functions venn2 and venn2_circles accept as their only required 
    # argument a 3-element list (A, B, AB) of subset sizes
    # e.g.: venn2(subsets = (3, 2, 1))
    # for colors: http://www.javascripter.net/faq/colornam.htm

    subset = ['A', 'B', 'A|B']
    subset_values, labels = get_venn_subset(group_list, peak_count, named_sets, ngroups, subset)
    venn2(subsets=(subset_values), set_labels=(labels), set_colors=('deepskyblue', 'orangered'), alpha=0.5, normalize_to=1)
    save_plot(plt, "{0}.png".format(get_file_info(output_file)[3]), mtitle)
    save_plot(plt, "{0}.pdf".format(get_file_info(output_file)[3]), mtitle)

def generate_venn3(group_list, peak_count, output_file, named_sets, ngroups, mtitle):
    # The functions venn3 and venn3_circles accept as their only required 
    # argument a 3-element list (A, B, AB, C, AC, BC, ABC) of subset sizes

    subset = ['A', 'B', 'A|B', 'C', 'A|C', 'B|C', 'A|B|C']
    subset_values, labels = get_venn_subset(group_list, peak_count, named_sets, ngroups, subset)

    venn3(subsets=(subset_values), set_labels = (labels),set_colors=('deepskyblue', 'orangered', 'g'), alpha=0.50, normalize_to=1)
    save_plot(plt, "{0}.png".format(get_file_info(output_file)[3]), mtitle)
    save_plot(plt, "{0}.pdf".format(get_file_info(output_file)[3]), mtitle)

def generate_venn4(group_list, peak_count, output_file, named_sets, mtitle, **kwds):
    ''' This function generates the venn diagram for 4 sets'''

    from matplotlib.patches import Circle, Ellipse
    alignment = {'horizontalalignment':'center', 'verticalalignment':'baseline'}

    # Define groups and their relative position
    groups = defaultdict()
    groups['A']       = (120, 200)
    groups['B']       = (280, 200)
    groups['C']       = (155, 250)
    groups['D']       = (245, 250)
    groups['A|B']     = (200, 115)
    groups['A|C']     = (140, 225)
    groups['A|D']     = (145, 155)
    groups['B|C']     = (255, 155)
    groups['B|D']     = (260, 225)
    groups['C|D']     = (200, 240)
    groups['B|C|D']   = (235, 205)
    groups['A|C|D']   = (165, 205)
    groups['A|B|D']   = (225, 135)
    groups['A|B|C']   = (175, 135)
    groups['A|B|C|D'] = (200, 175)

    # Set figure size
    if 'figsize' in kwds and len(kwds['figsize']) == 2:
        # if 'figsize' is in kwds, and it is a list or tuple with length of 2
        figsize = kwds['figsize']
    else: # default figure size
        figsize = (10, 11)

    # Set colors for different Circles or ellipses
    if 'colors' in kwds and isinstance(kwds['colors'], Iterable) and len(kwds['colors']) >= 4:
        colors = kwds['colors']
    else:
        colors = ['orange', 'r', 'b', 'c',]

    # Draw ellipse, the coordinates are hard coded in the rest of the function
    fig = plt.figure(figsize=figsize, frameon = False)   # set figure size
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    patches = []
    width, height = 170, 110  # width and height of the ellipses
    patches.append(Ellipse((170, 170), width, height, -45 , color=colors[0], alpha=0.4))
    patches.append(Ellipse((200, 200), width, height, -45 , color=colors[1], alpha=0.2))
    patches.append(Ellipse((200, 200), width, height, -135, color=colors[2], alpha=0.2))
    patches.append(Ellipse((230, 170), width, height, -135, color=colors[3], alpha=0.2))
    for e in patches:
        ax.add_patch(e)
    ax.set_xlim(80, 320); ax.set_ylim(80, 320)
    ax.set_xticks([]); ax.set_yticks([]);
    ax.set_aspect("equal")
    
    # Print the groups in the plot
    group_mapping = dict(zip(sorted(group_list),['A','B','C','D']))
    for g in get_groups_combination(sorted(group_list)):
        map_key = '|'.join([group_mapping[v] for v in g.split('|')])
        if map_key in groups:
            plt.text(groups[map_key][0], groups[map_key][1],peak_count[g], **alignment)

    # Annotate group names
    sorted_group_list = sorted(group_list)
    plt.text(110, 110,sorted_group_list[0] + "(" + str(len(named_sets[sorted_group_list[0]])) + ")", **alignment)
    plt.text(290, 110,sorted_group_list[1] + "(" + str(len(named_sets[sorted_group_list[1]])) + ")", **alignment)
    plt.text(130, 275,sorted_group_list[2] + "(" + str(len(named_sets[sorted_group_list[2]])) + ")", **alignment)
    plt.text(270, 275,sorted_group_list[3] + "(" + str(len(named_sets[sorted_group_list[3]])) + ")", **alignment)
    save_plot(plt, "{0}.png".format(get_file_info(output_file)[3]), mtitle)
    save_plot(plt, "{0}.pdf".format(get_file_info(output_file)[3]), mtitle)


def generate_venn5(group_list, peak_count, output_file, named_sets, mtitle, **kwds):
    ''' This will add the numbers to the 5-set venn diagram image '''
    from matplotlib.patches import Circle, Ellipse

    groups = defaultdict()
    groups['A']         = (200,270) 
    groups['B']         = (270,215) 
    groups['C']         = (215,115) 
    groups['D']         = (115,140) 
    groups['E']         = (100,225) 

    groups['A|B']       = (242,231) 
    groups['A|C']       = (186,130) 
    groups['A|D']       = (208,252) 
    groups['A|E']       = (166,253) 
    groups['B|C']       = (234,157) 
    groups['B|D']       = (126,172) 
    groups['B|E']       = (248,188) 
    groups['C|D']       = (156,146) 
    groups['C|E']       = (132,239) 
    groups['D|E']       = (117,195)

    groups['A|B|C']     = (204,152)  
    groups['A|B|D']     = (232,237)  
    groups['A|B|E']     = (234,196)  
    groups['A|C|D']     = (169,142)  
    groups['A|C|E']     = (163,245)  
    groups['A|D|E']     = (188,242)  
    groups['B|C|D']     = (145,166)  
    groups['B|C|E']     = (237,167)  
    groups['B|D|E']     = (124,186)  
    groups['C|D|E']     = (140,221)  

    groups['A|B|C|D']   = (170,156) 
    groups['A|B|C|E']   = (217,168) 
    groups['A|B|D|E']   = (215,225)
    groups['A|C|D|E']   = (163,237)
    groups['B|C|D|E']   = (133,200)
    groups['A|B|C|D|E'] = (180,195)

    # Set figure size
    if 'figsize' in kwds and len(kwds['figsize']) == 2:
        # if 'figsize' is in kwds, and it is a list or tuple with length of 2
        figsize = kwds['figsize']
    else: # default figure size
        figsize = (15, 16)

    # Set colors for different Circles or ellipses
    if 'colors' in kwds and isinstance(kwds['colors'], Iterable) and len(kwds['colors']) >= 4:
        colors = kwds['colors']
    else:
        colors = ['orange', 'r', 'b', 'c', 'g']

    # Draw ellipse, the coordinates are hard coded in the rest of the function
    fig = plt.figure(figsize=figsize, frameon = False)   # set figure size
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    patches = []
    width, height = 180, 90  # width and height of the ellipses
    patches.append(Ellipse((170, 190), width, height, 45  , color=colors[0], alpha=0.20)) # D
    patches.append(Ellipse((180, 210), width, height, 165 , color=colors[1], alpha=0.20)) # E
    patches.append(Ellipse((201, 213), width, height, 83  , color=colors[2], alpha=0.20)) # A
    patches.append(Ellipse((212, 195), width, height, 10  , color=colors[3], alpha=0.20)) # B
    patches.append(Ellipse((190, 175), width, height, 125 , color=colors[4], alpha=0.20)) # C
    for e in patches:
        ax.add_patch(e)
    ax.set_xlim(80, 320); ax.set_ylim(80, 320)
    ax.set_xticks([]); ax.set_yticks([]);
    ax.set_aspect("equal")

    # Print the groups in the plot
    alignment = {'horizontalalignment':'left', 'verticalalignment':'baseline'}
    group_mapping = dict(zip(sorted(group_list),['A','B','C','D','E']))
    for g in get_groups_combination(sorted(group_list)):
        map_key = '|'.join([group_mapping[v] for v in g.split('|')])
        if map_key in groups:
            plt.text(groups[map_key][0], groups[map_key][1],peak_count[g], fontsize=14, **alignment)

    # Annotate group names
    sorted_group_list = sorted(group_list)
    plt.text(100, 110, sorted_group_list[3] + "(" + str(len(named_sets[sorted_group_list[3]])) + ")", fontsize=14, **alignment) # D
    plt.text(100, 260, sorted_group_list[4] + "(" + str(len(named_sets[sorted_group_list[4]])) + ")", fontsize=14, **alignment) # E
    plt.text(220, 300, sorted_group_list[0] + "(" + str(len(named_sets[sorted_group_list[0]])) + ")", fontsize=14, **alignment) # A
    plt.text(255, 242, sorted_group_list[1] + "(" + str(len(named_sets[sorted_group_list[1]])) + ")", fontsize=14, **alignment) # B
    plt.text(255, 110, sorted_group_list[2] + "(" + str(len(named_sets[sorted_group_list[2]])) + ")", fontsize=14, **alignment) # C
    save_plot(plt, "{0}.png".format(get_file_info(output_file)[3]), mtitle)
    save_plot(plt, "{0}.pdf".format(get_file_info(output_file)[3]), mtitle)

def get_binary_string_for_combination(items):
    ''' 
    This gets an input a list of items. Then for each combination of items in the list, 
    it gets a binary number and returns the item->binary_number dictionary.
    
    For example: 
    items: [A,B,C]
    k=[list(compress(items,mask)) for mask in product(*[[0,1]]*len(items))]
    k
    [[], ['C'], ['B'], ['B', 'C'], ['A'], ['A', 'C'], ['A', 'B'], ['A', 'B', 'C']]

    ks = map(str,k)
    ks
    ['[]', "['C']", "['B']", "['B', 'C']", "['A']", "['A', 'C']", "['A', 'B']", "['A', 'B', 'C']"]

    v=[mask for mask in itertools.product(*[[0,1]]*len(items))]
    v=[(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1), (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1)]

    d=dict(itertools.izip(ks,v))
    d
    {"['A', 'B']": (1, 1, 0), "['B', 'C']": (0, 1, 1), "['B']": (0, 1, 0), '[]': (0, 0, 0), "['A', 'B', 'C']": (1, 1, 1), "['A']": (1, 0, 0), "['C']": (0, 0, 1), "['A', 'C']": (1, 0, 1)}

    '''
    from itertools import product, compress, izip
    # Get the keys in a list
    #k = [list(compress(items,mask)) for mask in product(*[[0,1]]*len(items))]
    k  = ['|'.join(list(compress(items,mask))) for mask in product(*[[0,1]]*len(items))]

    # Since you can not use a list as a dictionary key (why: please see the link below),
    # convert the list items of lists into strings
    # link: http://stackoverflow.com/questions/7257588/why-cant-i-use-a-list-as-a-dict-key-in-python
    #ks = map(str,k)
    
    # Get the binary values of the keys
    v=[mask for mask in itertools.product(*[[0,1]]*len(items))]
        
    # return the dictionary 
    #return dict(itertools.izip(ks,v))
    return dict(itertools.izip(k,v))

def get_venn_subset(group_list, peak_count, named_sets, ngroups, subset):
    ''' Get the correct subsets for venn function '''

    # annotate groups
    sorted_group_list = sorted(group_list) # ['ACC', 'CA1']

    # Generate lables with the counts of elements: ['ACC(78)', 'CA1(102)']
    labels = [g + "(" + str(len(named_sets[g])) + ")" for g in sorted_group_list]
    
    # Get all the combinations of the groups: ['A', 'B', 'A|B']
    groups = [g for g in get_groups_combination(sorted(['A','B','C','D','E'][0:int(ngroups)]))]

    # Map groups to the alphabets: {'ACC': 'A', 'CA1': 'B'}
    group_mapping = dict(zip(sorted(group_list),['A','B','C','D','E'][0:int(ngroups)]))

    # Get the pairing for venn functions
    venn_pairing = defaultdict()
    for g in get_groups_combination(sorted(group_list)):
        # For each g('ACC', 'CA1', 'ACC|CA1'), 
        # get the map_key({'ACC': 'A', 'CA1': 'B', 'ACC|CA1': 'A|B'})
        map_key = '|'.join([group_mapping[v] for v in g.split('|')])
        if map_key in groups:
            venn_pairing[map_key] = int(peak_count[g])
        
    pairing = [venn_pairing[v] for v in subset]
    return pairing, labels

def save_plot(plt, output_file, mtitle):
    ''' Save the plot with proper title and features '''

    # Add the title
    plt.title(mtitle)

    # print default png
    plt.savefig(output_file, dpi=300)
    
    # print pdf
    plt.savefig(get_file_info(output_file)[0] + "/" + get_file_info(output_file)[1] + ".pdf", dpi=300)
    #plt.show()
    
def venn_count(named_sets):
    ''' Finds the number of items that are found only in the intersection for each combination of sets '''

    names = set(named_sets)
    for i in range(1, len(named_sets)+1):
        # for each group of file/set 
        # ex. i = 1 is single (Caco2, Calu3, Capan1, etc.)
        # ex. i = 2 is double (Caco2/Calu3, Caco2/Capan1, etc.)
        for to_intersect in combinations(sorted(named_sets), i):
            #Loop through all the combinations of the files/set
            # ex. to_intersect = ('Caco2',) for i=1
            # ex. to_intersect = ('Caco2', 'Calu3') for i=2

            # Get the rest of the combination for this set
            # ex. others = set(['Capan1', 'HepG2', 'Calu3', 'GM12878']) for i=1
            # ex. others = set(['HepG2', 'GM12878', 'Capan1']) for i=2
            others = names.difference(to_intersect)

            # Get all the peaks that are the intersection of the sets in to_intersect
            # for example if to_intersect = ('Caco2'), then get all the peaks with Caco2
            # for example if to_intersect = ('Caco2', 'Calu3'), then get all the peaks that common between Caco2 and Calu3. They can also be there in other groups.
            intersected = set.intersection(*(named_sets[k] for k in to_intersect))
    
            # Get all the peaks that are there in 'others'
            unioned = set.union(*(named_sets[k] for k in others)) if others else set()

            # Return the peaks that are only there in the intersaction of the above ... means remove all the peaks from intersected that are there in unioned
            yield to_intersect, others, len(intersected - unioned), intersected-unioned

def get_groups_combination(symbols):
    ''' Generate all possible combinations of the groups from a sequence '''
    max_length = len(symbols)
    for length in xrange(1, max_length + 1):
        for group in map('|'.join, combinations(symbols, length)):
            yield group
            
def get_interactions_dict(input_file, column_nums, caseIns):
    ''' For the input file of the columns, get the data from the relevant colums  '''

    # create an empty set and dict
    group_list = list()
    interactions_set = defaultdict(set)
    with open(input_file, 'rU') as f:
        # ACC	CA1	CA3	DG	Blood							mmu-miR-423-3p		common signature
        # mmu-miR-434-5p	mmu-miR-6240	mmu-miR-378d	mmu-miR-10a-5p	mmu-miR-423-3p	Yes	mmu-miR-423-3p					mmu-miR-143-3p		
        # mmu-miR-143-3p	mmu-miR-1249-3p	mmu-miR-6240	mmu-miR-425-3p	mmu-miR-143-3p	YES	mmu-miR-143-3p					mmu-miR-27a-3p	
        
        # Get the headers
        headers = f.readline().strip().split("\t")
        for i in column_nums:
            group_list.append(headers[i])
    f.close()
    
    # Extract information of each column and then store it in the datastructre
    data = np.genfromtxt(input_file, comments='#', delimiter="\t", skip_header=1, usecols=column_nums, dtype=str)
    if caseIns:
        data = np.char.lower(data)
    
    # Remove empty strings from a list of strings
    # str_list = filter(None, str_list) # fastest
    
    # Get the entires for all the groups
    for i,c in enumerate(column_nums):
        interactions_set[headers[c]] = set(filter(None, data[:,i]))
    
    return interactions_set, group_list

def save_tabular(all_elements_info, sorted_group_list, output_file):
    ''' Saves the element distribution for the regions in a tabular format '''

    # Get the output file
    print "\n- Your output tabular file is: " + output_file
    of = open(output_file , 'w')
    
    # Write the header
    of.write("#element\t"+"\t".join(sorted_group_list) + "\n")

    # Get the categories for each featureid
    item_map = get_binary_string_for_combination(sorted_group_list)

    # Iterate all featureid`s and save it with its respective categories in the output file
    for key, val in all_elements_info.iteritems():
        line = key + "\t" + "\t".join(map(str,item_map[str(val)]))
        of.write(line)
        of.write("\n")

def print_help():
    ''' Print system help '''
    print >> sys.stderr, "\n ----------------- HELP ------------------\n", parser.print_help(), "\n"

def check_options():
    ''' Checks the options to the program '''

    # Create parser object
    global parser
    parser = argparse.ArgumentParser(add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=textwrap.dedent('''\
        ----------------- SAMPLE USAGE ------------------
        - python scripts/venn_plots.py -if=test/test_input.txt -of=test/01_two_sets_venn.txt   -cl=0,1
        - python scripts/venn_plots.py -if=test/test_input.txt -of=test/02_three_sets_venn.txt -cl=0,1,3
        - python scripts/venn_plots.py -if=test/test_input.txt -of=test/03_four_sets_venn.txt  -cl=0,1,4,6
        - python scripts/venn_plots.py -if=test/test_input.txt -of=test/04_five_sets_venn.txt  -cl=1,4,5,6,7
        - python scripts/venn_plots.py -if=test/test_input.txt -of=test/05_eight_sets_venn.txt
        -------------------------------------------------
        '''))

    # Add arguments 
    parser.add_argument("-if", metavar="--input_file"  , help="*Input file (tab separated) with headers", dest="input_file", type=str, required=True)
    parser.add_argument("-of", metavar="--output_file" , help=" Output file", dest="output_file", type=str)
    parser.add_argument("-od", metavar="--output_dir"  , help=" Output directory", dest="output_dir", type=str)
    parser.add_argument("-cl", metavar="--cols"        , help=" Comma separated list of columns you want to use to draw venn diagrams: \"col0,col1,col2,col5,col9\" \n for example:\n\t-cl=\"0,1,2,5,9\"\n Its ZERO based indexing i.e. Count the first column as zero ", dest="column_nums", type=str)
    parser.add_argument('-ci', "--caseInsensitive"     , help="if set, use case insensitive comparisions", action='store_true', default=False, dest="caseIns")

    # Print the help message only if no arguments are supplied
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    # Save the STDOUT output in a log file
    if parser.parse_args().output_file or parser.parse_args().output_dir:
        if parser.parse_args().output_file:
            logdir = get_file_info(parser.parse_args().output_file)[0]
            logfile = get_file_info(parser.parse_args().output_file)[3] + ".log"

        if parser.parse_args().output_dir:
            logdir = "{0}/logs".format(parser.parse_args().output_dir)
            logfile = "{0}/{1}.log".format(logdir, get_file_info(sys.argv[0])[1])
    elif parser.parse_args().input_file:
        logdir = get_file_info(parser.parse_args().input_file)[0]
        logfile = get_file_info(parser.parse_args().input_file)[3] + ".log"
    else:
        logdir  = "{0}/logs".format(os.getcwd())
        logfile = "{0}/{1}.log".format(logdir, get_file_info(sys.argv[0])[1])
    
    create_dir(logdir)
    logf = open(logfile, 'w')
    sys.stdout = Log(logf, sys.stdout)

    # Parse command line with parse_args and store it in an object
    args = parser.parse_args()
    print_initial_arguments(parser)
    return args

if __name__=="__main__":
      main()


