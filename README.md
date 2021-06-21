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

We provide a small demo data [here](https://github.com/XWangLabTHU/ARIC/tree/main/data/demo). 
There are two main files in csv format. One saves the mixture bulk data and another saves external reference data. Just put the file path to the function "**ARIC**", and the program will do every thing.

```python
from ARIC import *

ARIC(mix_path="mix.csv", ref_path="ref.csv")
```

### Section 3.2: Function Introduction

The main function in ARIC is **decipher**.

```Python
ARIC(mix_path, ref_path, save_path=None, marker_path=None, 
     selected_marker=False, scale=0.1, delcol_factor=10,
     iter_num=10, confidence=0.75, w_thresh=10, 
     unknown=False, is_methylation=False)
```

+ **'mix_path'**: Path to mixture data, must be an csv file with colnames and rownames.
+ **'ref_path'**: Path to reference data, must be an csv file with colnames and rownames.
+ **'save_path'**: Where to save the deconvolution results. Default: mix_path_prefix_prop.csv.
+ **'marker_path'**: Path to the user specificed markers. Must be an csv file.
+ **'selected_marker'**: Output selected marker for every sample. Marker files will be saved in a folder named "sample_marker.csv".
+ **'scale'**: Used for controlling the convergence of SVR. A smaller value makes the convergence much faster. Default: 0.1.
+ **'delcol_factor'**: Used for controlling the extent of removing collinearity. Default: 10.
+ **'iter_num'**: Iterative numbers of outliers detection. Default: 10.
+ **'confidence'**: Ratio of remained markers in each outlier detection loop. Default: 0.75.
+ **'w_thresh'**: Threshold to cut the weights designer. Default: 10.
+ **'unknown'**: Whether to estimate unknown content proportion.
+ **'is_methylation'**: Whether the data type belongs to methylation data. If true, preliminary marker selection will be performed.