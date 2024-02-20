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
All the data necessary to test the package can be downloaded from the following link : https://drive.google.com/file/d/1uT9X9I-YLRJA79mpt0XKkRvRpXj_W-C7/view?usp=drive_link  
Moreover, mascot DAT files, which are necessary to obtain expect tables (see below) can be downloaded from the following link : https://drive.google.com/file/d/18CyChHcjhUHyuan3YNMDQDq8KKldZc-p/view?usp=drive_link  

The second dataset's datas (low signal) can be downloaded from the following link : https://drive.google.com/file/d/102mIL-4EiQGEPV4j5d7v3q75qsEQvS6j/view?usp=drive_link

## Step by step Demo
The following steps will demonstrate how to generate the LineageFilter's data for one sample : Candida_MiniMix100  
For a demo on all the samples, please see below (Simple Demo part).  
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
all_tables_Expect = LineageFilter.extractEXPECT(query_peprank2expect, "DATmascot\\Candida_MiniMix100.dat", Lineage_tab, Prot2Tax)
```

Finally, the LineageFilter's data for this sample can be obtained as follows :
```python
LF_dat = LineageFilter.getLFdat(all_tables_Quanti, all_tables_Expect, Weight_tree, FULL_Lineage, TaxID_ChangeLOG)
```

Moreover, the LineageFilter's data can be annotated (for the training/testing process, see below) as follows :
```python
#Get positive identifications for this sample
with open("DATA\\TaxID__TP\\Candida_MiniMix100.txt",'r') as f:
    list_TaxID_TP = [int(i.rstrip()) for i in f.readlines()]

#Annotate this sample's data as a test sample
Test_data = LineageFilter.Create_TRAIN_TEST_data(LF_dat, list_TaxID_TP, test=True)
#Annotate this sample's data as a train sample
Train_data = LineageFilter.Create_TRAIN_TEST_data(LF_dat, list_TaxID_TP, test=False)
```

## Getting Unipept's results with LF filtering
The following steps will demonstrate how to generate a table containing all the genera validated by the LF method (threshold = LF_threshold parameter of LineageFilter.get_FilteredMPA function, default 0.30), using data generated as previously described.  
This table contains the the taxID and name of each genus, and its lineage (superkingdom, phylum, class, order and family), as well as a quantification feature (quantification parameter of LineageFilter.get_FilteredMPA function, default : 'specific', corresponds to Unipept's default value (MPA table of the website, filtered by LF)).
```python
import LineageFilter

#Import necessary data
TaxID_ChangeLOG = pd.read_csv("DATA\\taxid-changelog__relevant.tsv", sep='\t')
path_to_MPA = "DATA\\UnipeptMPA_results\\MPA_results__0.050\\Candida_MiniMix100_mpa.csv"
path_saveRF = "DATA\TRAIN_TEST\\Unipept___MetaRF__genus_0.050.sav"

#Get predictions values
LF_dat = LineageFilter.predict_class(LF_dat, path_saveRF)

#Get Unipept's filtered table using LF method
MPA_save_table = LineageFilter.get_FilteredMPA(LF_dat, path_to_MPA, TaxID_ChangeLOG, LF_threshold=0.15, quantification='specific')
```

## Simple Demo
The following code demonstrate how this package can be used, by showing a full pipeline of its use : generate the train and test datas (at a given p-value), train the random forest model and test it.
```python
import LineageFilter

#Import others necessary python packages
import pandas as pd
import os
import pickle

#Import necessary data
Weight_tree = pd.read_csv("DATA\\Weight_Tree.tsv", sep='\t')
FULL_Lineage = pd.read_csv("DATA\\Taxonomy2021_FULL.tsv", sep='\t')
TaxID_ChangeLOG = pd.read_csv("DATA\\taxid-changelog__relevant.tsv", sep='\t')

#Set the names of the quantifications columns, the taxonomical levels considered, and the samples used to test the model
QnonSPE_colname = '# PEPS_shared'  ;  QSPE_colname = '# spePEPS'  ;  levels_LF = ['phylum', 'class', 'order', 'family', 'genus']
names_TEST_dat = ['Q28078_20210309_M21_D1.txt', 'Run1_C1_2000ng.txt', 'Run1_U1_2000ng.txt', 'Run2_P1_2000ng.txt', 'SIHUMI_4bandes.txt']

#Set the p-value considered, and the paths to the quantifications and expect tables
pval = '0.050'
path_QuantiTABs = "DATA\\Quanti_tables\\Unipept_"+pval+"\\"
path_ExpectTABs = "DATA\\Expect_tables\\Pvalue_"+pval+"\\"

List_DATAtrain = []  ;  List_DATAtest = []
for file in os.listdir("DATA\\TaxID__TP\\"):
    print('\n\t\t'+file.split('.txt')[0]+'\n')
    print('\tGetting Quantification tables')  ;  all_tables_Quanti = []
    for level in levels_LF:
        all_tables_Quanti.append(pd.read_csv(path_QuantiTABs+file.split('.txt')[0]+'__'+level+'.tsv', sep='\t'))
    
    print('\tGetting Expect tables')  ;  all_tables_Expect = []
    for level in levels_LF:
        all_tables_Expect.append(pd.read_csv(path_ExpectTABs+file.split('.txt')[0]+'__'+level+'.tsv', sep='\t'))
    
    print('Getting LF tab')
    LF_dat = LineageFilter.getLFdat(all_tables_Quanti, all_tables_Expect, Weight_tree, FULL_Lineage, TaxID_ChangeLOG,
                      QnonSPE_colname=QnonSPE_colname, QSPE_colname=QSPE_colname)
    
    with open("DATA\\TaxID__TP\\"+file,'r') as f:
        list_TaxID_TP = [int(i.rstrip()) for i in f.readlines()]
    if file in names_TEST_dat:
        tmp_test = LineageFilter.Create_TRAIN_TEST_data(LF_dat, list_TaxID_TP, test=True)  ;  tmp_test['sample_ID'] = file.split('.txt')[0]
        List_DATAtest.append(tmp_test)
    else:
        tmp_train = LineageFilter.Create_TRAIN_TEST_data(LF_dat, list_TaxID_TP)  ;  tmp_train['sample_ID'] = file.split('.txt')[0]
        List_DATAtrain.append(tmp_train)

TRAIN_data = pd.concat(List_DATAtrain)  ;  TEST_data = pd.concat(List_DATAtest)

#Set the path to save the random forest model
path_saveRF = "DATA\TRAIN_TEST\\Unipept___MetaRF__genus_0.050.sav"
#Create the model
LineageFilter.Create_RFmodel(TRAIN_data, path_saveRF)

#Test the model on the test data
TEST_data_pred = LineageFilter.predict_class(TEST_data, path_saveRF)
```
