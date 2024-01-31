# LineageFilter

LineageFilter is a Python library enabling the estimation of the taxonomic composition of complex samples using metaproteomics and machine learning.

## Install
After downloading the package file (.whl file), install it as follows:
```bash
#Unix/macOS
python3 -m pip install --force-reinstall .\lineagefilter-0.0.1-py3-none-any.whl
#Windows
py -m pip install --force-reinstall .\lineagefilter-0.0.1-py3-none-any.whl
```
## Data download
All the data necessary to test the package can be downloaded from the following link : https://drive.google.com/file/d/1trTl9YnHN0_WgnvESQ4UOEt7hWT8cbBh/view?usp=drive_link  
Moreover, mascot DAT files, which are necessary to obtain expect tables (see below) can be downloaded from 

## Step by step Demo
The following steps will demonstrate how to generate the LineageFilter's data for one sample : Candida_MiniMix100
Start by loading LineageFilter package, as well as others necessary python packages and data:
```python
import LineageFilter

#Import others necessary python packages
import pandas as pd
import os
import pickle

#Import necessary data
##Used to create Quantification & Expect tables
Lineage_tab = pd.read_csv("DATA\\Taxonomy2021_levelsRelevant.tsv", sep='\t')
Prot2Tax = pd.read_csv("DATA\\Prot2TaxID_NCBInrS2021.tsv", sep='\t')

##Used in LineageFilter process
Weight_tree = pd.read_csv("DATA\\Weight_Tree.tsv", sep='\t')
FULL_Lineage = pd.read_csv("DATA\\Taxonomy2021_FULL.tsv", sep='\t')
TaxID_ChangeLOG = pd.read_csv("DATA\\taxid-changelog__relevant.tsv", sep='\t')
```

The next step is to create the quantification's tables for this sample. This can be done using Unipept's Metaproteomics analysis output file as follows :
```python
all_tables_Quanti = LineageFilter.convertUNIPEPT_mpaOUT("DATA\\UnipeptMPA_results\\MPA_results__0.050\\Candida_MiniMix100_mpa.csv", sep='\t')
```

Then, expect tables can be obtained using the extractEXPECT function. To do so, the mascot DAT file, as well as a python dictionary (key=(spectrum's query number, number of the peptide matched), value=expect) are necessary :
```python
with open("DATA\\dict_queries\\Candida_MiniMix100.pkl", 'rb') as f:
    query_peprank2expect = pickle.load(f)
all_tables_Expect = LineageFilter.extractEXPECT(query_peprank2expect, path_to_DATfile, Lineage_tab, Prot2Tax)

```

## Simple Demo

```python
from github import Github

# Authentication is defined via github.Auth
from github import Auth

# using an access token
auth = Auth.Token("access_token")

# First create a Github instance:

# Public Web Github
g = Github(auth=auth)

# Github Enterprise with custom hostname
g = Github(base_url="https://{hostname}/api/v3", auth=auth)

# Then play with your Github objects:
for repo in g.get_user().get_repos():
    print(repo.name)

# To close connections after use
g.close()
```
