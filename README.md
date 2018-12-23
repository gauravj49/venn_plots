Venn diagram plotting script
====================================================

## Description: 
A script for plotting:
* Scaled two and three set venn diagrams
* Unscaled four and five set venn diagrams 
* Profile, summary and tabular output for all sets including (six sets or more)

## Dependencies:
* ``numpy``
* ``scipy``
* ``matplotlib``
* ``matplotlib-venn``

## Usage:
```
[gjain@gjvbx venn_plots ] $ python scripts/venn_plots.py 

***********************************************
- PROGRAM: venn_plots.py
- CONTACT: Gaurav Jain(gauravj49@gmail.com)
***********************************************

usage: venn_plots.py [-h] -if --input_file [-of --output_file]
                     [-od --output_dir] [-cl --cols] [-ci]

arguments:
-h, --help               show this help message and exit
-if --input_file        *Input file (tab separated) with headers
-of --output_file        Output file
-od --output_dir         Output directory
-cl --cols               Comma separated list of columns you want to use to draw venn diagrams: "col0,col1,col2 col5,col9" 
                         for example: -cl="0,1,2,5,9"
                         Its ZERO based indexing i.e. Count the first column as zero 
-ci, --caseInsensitive   if set, use case insensitive comparisions

----------------- SAMPLE USAGE ------------------
- python scripts/venn_plots.py -if=test/test_input.txt -of=test/01_two_sets_venn.txt   -cl=0,1
- python scripts/venn_plots.py -if=test/test_input.txt -of=test/02_three_sets_venn.txt -cl=0,1,3
- python scripts/venn_plots.py -if=test/test_input.txt -of=test/03_four_sets_venn.txt  -cl=0,1,4,6
- python scripts/venn_plots.py -if=test/test_input.txt -of=test/04_five_sets_venn.txt  -cl=1,4,5,6,7
- python scripts/venn_plots.py -if=test/test_input.txt -of=test/05_eight_sets_venn.txt
-------------------------------------------------
```
