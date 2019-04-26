# Discrinimation dataset
Here we provide two datasets for the purpose of discrmination of transmembrane (TM) and non-transmembrane (non-TM):
```
   (i)  440 dataset for TM , and
   (ii) 6418 dataset for non-TM.
```
Note that the TM dataset only contains alpha-helix TM proteins.


# Data description
Below we describe the details in each dataset (use XX to indicate '440_TM' and '6418_nonTM' in TM_dataset and nonTM_dataset, respectively)

```
1) data list:
   XX_list     -> the entry is in '1pdbA' format where the first 4-character indicate the PDB id and the 5th character indicates chain id

2) amino acid sequence as the original input:
   XX_fasta    -> the L*1 amino acide sequence

3) transmembrane topology label (2-state) only considering alpha helix:
   XX_truth    -> the L*1 ground-truth

4) 2-state prediction results from 4 different methods:
   XX_phobius  -> results from Phobius
   XX_philius  -> results from Philius
   XX_topcons  -> results from TOPCONS2 (Web Server)
   XX_purestm  -> results from our method
```


Note that for TM_dataset, we additionally provide the ground-truth label from PDBTM:
```
   440_TM_pdbtm       -> the 9-state TM topology label
```

To transfer 9-state from PDBTM to 2-state TM topology label, run the below command:
```
   ../util/pdbtm2binary.py 'input_pdbtm' > 'output_topology'
```

