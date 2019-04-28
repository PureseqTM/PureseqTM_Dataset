# Training and testing dataset
Here we provide two datasets for the prediction of transmembrane topology label:
```
   (i) 328 dataset for training and validation, and
   (ii) 39 dataset for testing.
```
Note that there are NO redundant sequences that share >25% identity between the two datasets.


# Data description
Below we describe the details in each dataset (use XX to indicate '328' and '39' in Train_dataset and Test_dataset, respectively)

```
1) data list:
   XX_list     -> the entry is in '1pdbA' format where the first 4-character indicate the PDB id and the 5th character indicates chain id

2) amino acid sequence as the original input:
   XX_fasta    -> the L*1 amino acide sequence

3) transmembrane topology label (9-state) from PDBTM annotation:
   XX_pdbtm    -> the L*1 ground-truth from PDBTM

4) transmembrane topology label (2-state) only considering alpha helix:
   XX_truth    -> the L*1 ground-truth

5) 2-state prediction results from 4 different methods:
   XX_phobius  -> results from Phobius
   XX_philius  -> results from Philius
   XX_topcons  -> results from TOPCONS2 (Web Server)
   XX_purestm  -> results from our method
```


Note that for Train_dataset, we additional provide the training list and validation list:
```
   164_train_list     -> the data for training our method
   164_validate_list  -> the data for validating our method during training
```

To transfer from 9-state PDBTM label to 2-state TM topology label, run the below command:
```
   ../util/pdbtm2binary.py 'input_pdbtm' > 'output_topology'
```

Type the below command to evaluate the performance of the 4 methods:
```
   ./evaluate_total
```


