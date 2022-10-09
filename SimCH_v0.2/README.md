---
Package release: SimCH (version 0.2, Oct 2022)  
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
- test
	- example.txt (the TXT file of an example count matrix)
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
  --unzip SimCH-main.zip in an available directory.

****All python packages should be installed in the Python3 directory and be run on Python3.***

# 4. Input files
- A count matrix file (genes by cells) of experimental scRNA-seq data in tab-seperated TXT file.  
The first row of the file records the cell IDs while the first column records the gene IDs. And each cell
records the numeric count of gene expression in each cell. The segments of each row are seperated by tabs.
Please see text/example.txt as reference.

- A TXT file annotating cell groups (optionally required for simulating on a tagged heterougenous dataset)  
The file contains only one row of cell group IDs that sequentially corresponds to the cell columns of the count matrix file.
Please see Li_GSE81861's data as reference.

# 5. Testing SimCH
## 5.1 Download test datasets
After SimCH installation, download test datasets from https://github.com/SIRG-YZU/SimCH/blob/main/datasets.zip .  
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
Enter the folder 'SimCH-main/SimCH_v0.1/SimCH' containing a Python script named 'SimCH.py' and run  

  ```$ python SimCH.py```  
  
And it will give you an interactive menu to perform simulation step by step, including 
selecting simulation mode, giving path to the count matrix of the experimental data, and parameter settings.  

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

****Please note that the options 1-6 can start parameter estimation for basic simulation (SimCH-flex-NB, SimCH-flex-NBZI, SimCH-fit-NB, SimCH-fit-NBZI, SimCH-copula-NB, or SimCH-copula-NBZI) as well as further extended simulation (SimCH-ext). The options 7-10 will start the independent muti-group simulation, which was used for the comparision between the simulators on tagged heterogeneous data and for the power evaluation as mentioned in the paper. And the option 11 will start simulation based on SimCH-copula-NB for evaluating imputation methods in terms of gene co-expression preservation.***

Then you can test above SimCH sub-modes on the test datasets as follows:  
- **Test 1** -- simulation on a homogeneous dataset, e.g. Ziegenhain_SmartSeq2B.txt  

  ```Step A: input '1' -->  give path to Ziegenhain_SmartSeq2.txt --> change seeds/parameters or not --> estimate parameters```  
  Waiting for parameters estimation ...  
  ```Step B: change simulation parameters or not --> run simulation --> give output directory```  

  After simulation, all output files will be stored in the given directory.  
  Likewise, you can choose other sub-modes,
  e.g. 2 to 6, to run the simulation on a homogeneous dataset or untagged heterogeneous one.
  

- **Test 2** -- simulation on a heterogeneous dataset with cell groups annotation, e.g. Li.txt  

  ```Step A: input '7' --> give path to Li.txt --> give the path to the cell groups file --> change seeds/parameters or not --> estimate parameters```  
  Waiting for parameters estimation ...  
  ```Step B: change simulation parameters or not --> run simulation (-s) --> give output directory```  

  After simulation, all output files will be stored in the given directory.  
  Likewise, you can choose other sub-modes,
  e.g. 8 to 10, to run the simulation on a tagged heterogeneous dataset.  
  
- **Test 3** -- estimate parameters from a homogeneous dataset (e.g. example.txt) and perform extended simulation based on SimCH-fit-NB

  ```Step A: input '3' -->  give path to example.txt --> change seeds/parameters or not --> estimate parameters```  
  Waiting for parameters estimation ...  After that, it will give following options:  

  ```
	-1      Group ratio (for extended simulation):[1]
	-2      Batch ratio (for extended simulation):[1]
	-3      Batch variation (for extended simulation):[0, 0.2]
	-4      Path (for extended simulation):No
	-5      DEgene ratio (for extended simulation):0.2
	-6      DEgene variation (for extended simulation):[0, 0.5]
	-7      Non-linear gene ratio (for extended simulation):0.01
	-8      Marker gene ratio (for extended simulation):0.01
	-9      Library magnification (for extended simulation):1
	-10     Cell number (for extended simulation):864
	
	If you want to modify the above parameters, please enter a number
	-v means view current parameters
	-s means start simulation
  ```
  You can set the parameters above to simulate multi-groups, batch effects and differentiation paths.  
- For simulate multi-groups, choose '-1' and set '[0.5,0.5]' (means two groups); choose '-5' and set '0.1' (DE genes ratio).  
- For simulate batch effects, choose '-2' and set '[0.5,0.5]' (means two batches); choose '-3' and set '[0,0.1]'.  
- For simulate differentiation path, choose '-1' and set '[0.3,[0.2,0.2],0.3]' (means three differentiation paths, and the second path has two sub-paths); choose '-4' and set 'Yes'.   

  ```Step B: change simulation parameters or not --> run simulation (-s) --> give output directory```

  The resulting output files mainly include:
- Cell gene mean(bcv2).txt
- Cell gene mean(size factor).txt
- Counts.txt
- Gene dispersion.txt
- Gene mean.txt
- Library size.txt
- Simulation information.txt
- Size factor.txt
- DEgene factor.txt (output of multi-group and differentiation path simulation)
- Group.txt (output of multi-group simulation)
- Marker gene.txt (output of multi-group and differentiation path simulation)
- Batch factor.txt (output of batch effects simulation)
- Batch.txt (output of batch effects simulation)
- Path.txt (output of differentiation path simulation)
- Nonlinear gene factor.txt (output of differentiation path simulation)

  where 'Counts.txt' is the simulated data and 'Simulation information.txt' contains basic parameter settings of the simulation. Other files involved in the simulation would be useful for benchmarking and evaluating the scRNA-seq computational methods.

