# ARIC

## Section 1: Introduction

**<font color=red>ARIC</font>** is a bioinfomatics software for bulk gene expression and DNA methylation data deconvolution. ARIC utilizes a novel two-step marker selection strategy, including **component-wise condition number-based feature collinearity elimination** and **adaptive outlier markers removal**. This strategy can systematically obtain effective markers that ensure a robust and precise **weighted υ-SVR-based** rare proportion prediction.

## Section 2: Installation Tutorial

### Section 2.1: System requirement
ARIC is implemented using python and can be install in windows, UNIX/LINUX and MAC OS. ARIC requires python version >= 3 and all the dependent packages will be installed using pip.

### Section 2.2: Installation
ARIC can be installed from pypi by the following command. The source code can also be downloaded from pypi.

```Shell
pip install ARIC
```

## Section 3: A Quick Tutorial for Demo data Deconvolution

In this section, we will demonstrate how to perform bulk data deconvolution using the demo data.

### Section 3.1: Quick Start

We provide a small demo data [here](https://github.com/XWangLabTHU/ARIC/tree/main/demo_data). 
There are two main files in csv format. One saves the mixture bulk data and another saves external reference data. Just put the file path to the function "**ARIC**", and the program will do every thing.

```python
from ARIC import *

ARIC(ref_path="ref.csv", mix_path="mix.csv")
```

### Section 3.2: Function Introduction

The main function in ARIC is **decipher**.

```Python
decipher(ref_path, mix_path, save_path='prop_predict.csv', 
         marker_path='', scale=0.1, delcol_factor=10, 
         iter_num=10, confidence=0.75, w_thresh=10, 
         unknown=False, is_markers=False, is_methylation=False)
```

+ **'ref_path'**: Path to reference csv file. The first row should be the name of cell types and the first column should be marker names. 
+ **'mix_path'**: Path to mixture csv file. The first row should be the samples id and the first column corresponds to the marker names.
+ **'save_path'**: Path to save the deconvolution results.
+ **'marker_path'**: Path to user provided markers. In csv file, one column. 
+ **'scale'**: Control the stop criterion for the SVR, smaller scale means earlier stop. We recommend scale=0.1 for RNA-Seq and scale=1 for methylation. default = 0.1.
+ **'delcol_factor'**: larger values means less loops in the eliminating collinearity part. It should be in the range from 2 to 20. default=10.
+ **'iter_num'**: The number of loops in the outlier detection process. default=10.
+ **'confidence'**: The ratio of markers treated as normal points in the first outlier detection loop. It should be in [0, 1]. default=0.75. 
+ **'w_thresh'**: The threshold used to cut weights which are too large. default=10.
+ **'unknown'**: If there is unknown component in the mixture data. default=False.
+ **'is_markers'**: False if users choose to select markers by ARIC. default=False
+ **'is_methylation'**: If the data is methylation data. default=False

