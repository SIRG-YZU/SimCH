---
Package release: SimCH (version 0.1, April 2022)  
Authors: Gongming Wang, Lei Sun  
Email: sunlei(at)yzu.edu.cn  
Description: SimCH is a Python-based simulator for benchmarking and evaluating scRNA-seq data
             computational methods.   
License: SimCH is distributed under the GNU GENERAL PUBLIC (GPL) license (Version 3).  
Copyright (C) 2022 Yangzhou University
---

# Contents
 1. Introduction
 2. Package components
 3. Installation
 4. Input files
 5. Testing SimCH

# 1. Introduction
SimCH is a Python-based simulator for bechmarking and evaluating scRNA-seq data processing and analysis methods.
It depends on several state-of-art packages, e.g. Numpy, Scipy, pandas, scanpy, sympy, etc.

# 2. Package components
- README.md (instructions)  
- ChangeLog (records of version changes)  
- SimCH  
    - SimCH.py (the main program)  
- LICENSE  
    - GPL_v3.0.txt  (GPL license, version 3.0)

# 3. Installation
SimCH can be run via command line interface (CLI) of Windows/Linux/macOS as follows:  
(1) Install below dependancy packages manually.
  - Python 3.7 (https://www.python.org/downloads/)
  - Numpy (https://pypi.org/project/numpy/)  
    ```pip3 install numpy```
  - Scipy (https://scipy.org)  
    ```pip3 install scipy```
  - pandas (https://pandas.pydata.org/getting_started.html)  
    ```pip3 install pandas```
  - matplotlib (https://matplotlib.org/users/installing.html)  
    ```pip3 install matplotlib```
  - scanpy (https://scanpy.readthedocs.io/en/stable/index.html)  
    ```pip3 install scanpy```
  - sympy (https://www.sympy.org/en/index.html)  
    ```pip3 install sympy```

(2) SimCH installation  
  --unzip SimCH_v0.1.zip in an available directory.

****All python packages should be installed in Python3 directory and be run on Python3.***

# 4. Input files
- A count matrix file (genes by cells) of real scRNA-seq data in tab-seperated TXT file.  
The first row of the file records the cell IDs while the first column records the gene IDs. And each cell
records the numeric count of gene expression in each cell. The segments of each row are seperated by tabs.

- A TXT file annotating cell groups (optionally required for simulating on a tagged heterougenous dataset)  
The file contains only one row of cell group IDs that sequentially corresponds to the cell columns of the count matrix file.

# 5. Testing SimCH
## 5.1 Download test datasets
After SimCH installation, download test datasets from https://sourceforge.net/projects/crisprvi/files/data/dataset-1.zip .  
The unpacked folder contains 13 test datasets as follows.
- Camp_GSE75140 - Camp.txt
- Klein_GSE65525 - Klein.txt
- Zeisel_GSE60361 - Zeisel.txt, Zeisel_3.txt, cell_groups.txt
- Avraham_GSE65528 - Avraham.txt 
- Darmanis_GSE67835 - Darmanis.txt
- Engel_GSE74596 - Engel.txt 
- Tung_GSE77288 - Tung.txt 
- Li_GSE81861 - Li.txt, cell_groups.txt
- Ziegenhain_GSE75790 - Ziegenhain_SCRBseqB.txt
- Ziegenhain_GSE75790 - Ziegenhain_CELseq2B.txt
- Ziegenhain_GSE75790 - Ziegenhain_SmartSeqA.txt
- Ziegenhain_GSE75790 - Ziegenhain_SmartSeq2B.txt
- Haber_GSE92332 - GSE92332_FAE_UMIcounts.txt

## 5.2 Simulation of SimCH
Enter the folder 'SimCH' containing a Python script named 'SimCH.py' and run  

  ```$ python SimCH.py```  
  
And it will give you an interactive menu to perform simulation step by step, including 
selecting simulation mode, giving path to the count matrix of the real data, and parameter settings.  

In the top menu, there are 11 simulation sub-modes:  
```
1       simulation based on SimCH-flex-NB
2       simulation based on SimCH-flex-NBZI
3       simulation based on SimCH-fit-NB
4       simulation based on SimCH-fit-NBZI
5       simulation based on SimCH-copula-NB
6       simulation based on SimCH-copula-NBZI
7       simulation of independent multiple groups based on SimCH-copula-NB
8       simulation of independent multiple groups based on SimCH-fit-NB
9       simulation of independent multiple groups based on SimCH-copula-NBZI
10      simulation of independent multiple groups based on SimCH-fit-NBZI
11      simulation based on SimCH-copula-NB for evaluating imputation methods
        in terms of gene co-expression preservation
```
Then you can test above SimCH sub-modes on the test datasets as follows:  
- **Test 1** -- simulation on a homogeneous dataset, e.g. Ziegenhain_SmartSeq2B.txt  

  ```Step A: input '1' -->  give path to Ziegenhain_SmartSeq2.txt --> change seeds/parameters or not --> estimate parameters```  
  Waiting for parameters estimation ...  
  ```Step B: change simulation parameters or not --> give output directory --> run simulation```  

  After simulation, all output files will be stored in the given directory.  
  Likewise, you can choose other sub-modes,
  e.g. 2 to 6, to run the simulation on a homogeneous dataset or untagged heterogeneous one.
  

- **Test 2** -- simulation on a heterogeneous dataset, e.g. Li.txt  

  ```Step A: input '7' --> give path to Li.txt --> give the path to the cell groups file --> change seeds/parameters or not --> estimate parameters```  
  Waiting for parameters estimation ...  
  ```Step B: change simulation parameters or not --> give output directory --> run simulation```  

  After simulation, all output files will be stored in the given directory.  
  Likewise, you can choose other sub-modes,
  e.g. 8 to 10, to run the simulation on a tagged heterogeneous dataset.  
  
The resulting output files mainly include:
- Cell gene mean(bcv2).txt
- Cell gene mean(size factor).txt
- Counts.txt
- Gene dispersion.txt
- Gene mean.txt
- Library size.txt
- Simulation information.txt
- Size factor.txt

where 'Counts.txt' is the simulated data and 'Simulation information.txt' contains basic information of the simulation. Other files involved in the simulation would be useful for benchmarking and evaluating the scRNA-seq computational methods.

