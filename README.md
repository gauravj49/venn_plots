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

usage: venn_plots.py [-h] -if --input_file [-of --output_file] [-od --output_dir] [-cl --cols] [-ci]

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

### Input file (tab-delimited text file with headers):

| Set1 | Set2 | Set3 | Set4 | Set5 | Set6 | Set7 | Set8 |
| ---- |----- | ---- | ---- | ---- |----- | ---- | ---- |
|   34 |   50 |   49 |    4 |   35 |   11 |   44 |   20 |
|   45 |   14 |    5 |   43 |   40 |   23 |   15 |    0 |
|    7 |   27 |   27 |    1 |   25 |   16 |   16 |   48 |
|   28 |   31 |   42 |   11 |   17 |   40 |   35 |   37 |
|   26 |   13 |   44 |   27 |    9 |   36 |   23 |   26 |
|   15 |   28 |   34 |   22 |   29 |   12 |   38 |   18 |
|    8 |   39 |   25 |   13 |    4 |   13 |   17 |   22 |
|   23 |    5 |   15 |    6 |   38 |   46 |   24 |   12 |
|   18 |    9 |   16 |   50 |   48 |    3 |   39 |   24 |
|    1 |   42 |   14 |    3 |   32 |   14 |   46 |   25 |
|   49 |    3 |    1 |   34 |   33 |   48 |   31 |   49 |
|   12 |   22 |   37 |   35 |   50 |   50 |   22 |    1 |
|   38 |    4 |   11 |   20 |   27 |   42 |   33 |   19 |
|   37 |   15 |   31 |   44 |   12 |   34 |    4 |   21 |
|   42 |   33 |   50 |   17 |   19 |   47 |    5 |    2 |
|   33 |   26 |   20 |   25 |    0 |   19 |   50 |   42 |
|   40 |      |   19 |   12 |   36 |    7 |   14 |   23 |
|   22 |      |    0 |      |   11 |   22 |      |    8 |
|   20 |      |   41 |      |   22 |   24 |      |   31 |
|   21 |      |    4 |      |      |      |      |   28 |
|    4 |      |      |      |      |      |      |   10 |

### Output files:

1. Plots (png and pdf) for 2 to 5 sets  
 ![plots-01](https://user-images.githubusercontent.com/10153240/50388588-5f8b8000-071b-11e9-9a3f-fdcd9cafacef.jpg)

1. Summary file:
   ```
   - Set1 = 14
   - Set2 = 9
   - Set1|Set2 = 7
   ```

1. Profile file:
    ```
    > Set1 = 14
    1
    12
    18
    20
    21
    23
    34
    37
    38
    40
    45
    49
    7
    8
    > Set2 = 9
    13
    14
    27
    3
    31
    39
    5
    50
    9
    > Set1|Set2 = 7
    15
    22
    26
    28
    33
    4
    42
    ```

1. Tabular file (1=PRESENT, 0 = ABSENT):

    | #element | Set1 | Set2 |
    |----------|------|------|
    |       28 |    1 |    1 |
    |       26 |    1 |    1 |
    |       27 |    0 |    1 |
    |       20 |    1 |    0 |
    |       21 |    1 |    0 |
    |       22 |    1 |    1 |
    |       49 |    1 |    0 |
    |       45 |    1 |    0 |
    |       42 |    1 |    1 |
    |       40 |    1 |    0 |
    |        1 |    1 |    0 |
    |        3 |    0 |    1 |
    |        5 |    0 |    1 |
    |        4 |    1 |    1 |
    |        7 |    1 |    0 |
    |        9 |    0 |    1 |
    |        8 |    1 |    0 |
    |       13 |    0 |    1 |
    |       38 |    1 |    0 |
    |       39 |    0 |    1 |
    |       12 |    1 |    0 |
    |       15 |    1 |    1 |
    |       14 |    0 |    1 |
    |       33 |    1 |    1 |
    |       18 |    1 |    0 |
    |       31 |    0 |    1 |
    |       23 |    1 |    0 |
    |       37 |    1 |    0 |
    |       50 |    0 |    1 |
    |       34 |    1 |    0 |

Courtsey GitDiagram

![image](https://github.com/user-attachments/assets/b83dbcb5-85a4-4789-acd7-a26688ab0c68)
