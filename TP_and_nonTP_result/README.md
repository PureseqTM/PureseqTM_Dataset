# Discrinimation dataset
Here we provide two datasets for the purpose of discrmination of transmembrane proteins (TPs) and non-transmembrane proteins (non-TPs):
```
   (i)  440 dataset for TPs, and
   (ii) 6418 dataset for non-TPs.
```
Note that the TP dataset only contains alpha-helix TP proteins.


# Data description
Below we describe the details in each dataset (use XX to indicate '440_TP' and '6418_nonTP' in TP_dataset and nonTP_dataset, respectively)

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


Note that for TP_dataset, we additionally provide the ground-truth label from PDBTM:
```
   440_TP_pdbtm       -> the 9-state PDBTM topology label
```

To transfer from 9-state PDBTM label to 2-state TM topology label, run the below command:
```
   ../util/pdbtm2binary.py 'input_pdbtm' > 'output_topology'
```

Type the below command to evaluate the performance of the 4 methods:
```
   ./evaluate_total
```

