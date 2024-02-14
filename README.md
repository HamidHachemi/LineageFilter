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

The second dataset's datas can be downloaded from the following links :  
LF's data : https://drive.google.com/file/d/1WzwSUzPO66xKu38piF1Dq-KgnQgmytUJ/view?usp=drive_link  
mascot DAT files : https://drive.google.com/file/d/1eUfAZH4uVD-ErYM6I-bsvXQ8QIJS2WFI/view?usp=drive_link

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
The following steps will demonstrate how to generate a table containing all the genera validated by the LF method (threshold = LF_threshold parameter of LineageFilter.get_FilteredMPA function, default 0.20), using data generated as previously described.  
This table contains the the taxID and name of each genus, and its lineage (superkingdom, phylum, class, order and family), as well as a quantification feature (quantification parameter of LineageFilter.get_FilteredMPA function, default : 'specific', corresponds to Unipept's default value (MPA table of the website, filtered by LF)).
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

## Evaluate the performances
Once the predictions have been done, the performances of the model can be assessed by plotting the ROC curve (TPR according to FPR) and the F1-score curve (F1-score according to FPR).  
This can be done, for the previously generated data, using the following code :
```python
#Import the necessary packages
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve
from matplotlib import rc
import numpy as np

#Assess Unipept performances
UNI_spepep = pd.DataFrame()
n_posTOT = sum(TEST_data_pred['Class']==True)  ;  n_negTOT = sum(TEST_data_pred['Class']==False)
for thr in list(range(int(max(TEST_data_pred[TEST_data_pred['Class']==1]['# spePEPS__genus']))))[::-1]:
    tmp = TEST_data_pred[TEST_data_pred['# spePEPS__genus']>=thr]  ;  TP = sum(tmp['Class']==True)  ;  FN = n_posTOT-TP
    FP = sum(tmp['Class']==False)  ;  TN = n_negTOT-FP
    if (TP!=0)|(FP!=0):
        UNI_spepep = pd.concat([UNI_spepep, pd.DataFrame({'# spePEPS':[thr], 'FPR':[FP/(FP+TN)], 'TPR':[TP/(TP+FN)],
                                                    'precision':[TP/(TP+FP)], 'F1':[2*TP/(2*TP+FP+FN)], 'thr':[thr]})])

#ROC curve
plt.figure(figsize=(20,20))  ;  linewidth = 6
##plot Unipept performances
plt.plot(UNI_spepep['FPR'], UNI_spepep['TPR'], label="Unipept", color='blue', linestyle='-', linewidth=linewidth)
##plot LineageFilter performances
true_y = TEST_data_pred['Class']  ;  y_prob = TEST_data_pred['pred']  ;  fpr, tpr, thresholds = roc_curve(true_y, y_prob)
thresholds = list(np.arange(0,1+1e-3,1e-3))  ;  thresholds = thresholds[::-1]  ;  FPRs = []  ;  TPRs = []
for i in range(len(thresholds)):
    TP = sum(true_y[y_prob>=thresholds[i]])  ;  FN = sum(true_y[y_prob<thresholds[i]])
    FP = sum([1 if j==0 else 0 for j in true_y[y_prob>=thresholds[i]]])  ;  TPRs.append(TP/(TP+FN))
    TN = sum([1 if i==0 else 0 for i in true_y[y_prob<thresholds[i]]])  ;  FPRs.append(FP/(FP+TN))

plt.plot(FPRs, TPRs, label="LineageFilter", color='green', linestyle='-', linewidth=linewidth)

plt.plot(np.arange(0,1,1e-5), np.arange(0,1,1e-5), color='black', linestyle='--', linewidth=linewidth/5)
plt.grid(visible=True)  ;  plt.legend(loc='lower right')
Xticks = np.array(list(np.arange(0,1e-4,1e-5))+list(np.arange(1e-4,1e-3,1e-4))+list(
    np.arange(1e-3,1e-2,1e-3))+list(np.arange(1e-2,1e-1,1e-2))+list(np.arange(1e-1,1,1e-1))+[1])
plt.gca().set_xscale('log')  ;  plt.xscale('log',base=10)  ;  plt.xticks(Xticks, minor=False)
plt.xlim(1e-5,1)  ;  plt.ylim(-1e-2,1)  ;  plt.yticks(np.arange(0,1.1,0.1), minor=False)

plt.xlabel('False Positive Rate')  ;  plt.ylabel('True Positive Rate')
rc('font', **{'family':'DejaVu Sans', 'weight':'bold', 'size':26})  ;  plt.rc('font', size=35)
plt.show()


#F1-score curve
plt.figure(figsize=(20,20))  ;  linewidth = 6
##plot Unipept performances
plt.plot(UNI_spepep['FPR'], UNI_spepep['F1'], label="Unipept", color='blue', linestyle='-', linewidth=linewidth)
##plot LineageFilter performances
true_y = TEST_data_pred['Class']  ;  y_prob = TEST_data_pred['pred']  ;  fpr, tpr, thresholds = roc_curve(true_y, y_prob)  ;  F1_scores = []
thresholds = list(np.arange(0,1+1e-3,1e-3))  ;  thresholds = thresholds[::-1]  ;  FPRs = []  ;  F1_scores = []
for i in range(len(thresholds)):
    TP = sum(true_y[y_prob>=thresholds[i]])  ;  FN = sum(true_y[y_prob<thresholds[i]])
    FP = sum([1 if j==0 else 0 for j in true_y[y_prob>=thresholds[i]]])  ;  F1_scores.append((2*TP)/(2*TP+FP+FN))
    TN = sum([1 if i==0 else 0 for i in true_y[y_prob<thresholds[i]]])  ;  FPRs.append(FP/(FP+TN))

plt.plot(FPRs, F1_scores, label="LineageFilter", color='green', linestyle='-', linewidth=linewidth)

plt.grid(visible=True)  ;  plt.legend(loc='lower right')
Xticks = np.array(list(np.arange(0,1e-4,1e-5))+list(np.arange(1e-4,1e-3,1e-4))+list(
    np.arange(1e-3,1e-2,1e-3))+list(np.arange(1e-2,1e-1,1e-2))+list(np.arange(1e-1,1,1e-1))+[1])
plt.gca().set_xscale('log')  ;  plt.xscale('log',base=10)  ;  plt.xticks(Xticks, minor=False)
plt.xlim(1e-5,1)  ;  plt.ylim(-1e-2,1)  ;  plt.yticks(np.arange(0,1.1,0.1), minor=False)

plt.xlabel('False Positive Rate')  ;  plt.ylabel('F1-score')
rc('font', **{'family':'DejaVu Sans', 'weight':'bold', 'size':26})  ;  plt.rc('font', size=35)
plt.show()
```
