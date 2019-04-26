# Human proteome from UniProt
Here we provide the detailed step of extracting Reviewed transmembrane (TM) regions from Human proteome.

# Data generation
## Step 1: get Human proteome from UniProt
```
url="https://www.uniprot.org/uniprot/?query=organism:9606+AND+reviewed:yes&sort=score&format=tab&columns=id,length"
wget -q $url -O human_20416_len
tail -n+2 human_20416_len | awk '{print $1}' > human_20416_list
```

Note that we manually delete entry 'P0DPR3' as this protein only contains two residues.


## Step 2: download FASTA sequence files
```
mkdir -p human_20416_fasta
for i in `cat human_20416_list`;
do
    wget -q https://www.uniprot.org/uniprot/$i.fasta -O human_20416_fasta/$i.fasta;
done
```

## Step 3: download GFF annotation files
```
mkdir -p human_20416_gff
for i in `cat human_20416_list`;
do
    wget -q https://www.uniprot.org/uniprot/$i.gff -O human_20416_gff/$i.gff;
done
```

## Step 4: extract TM regions from GFF
```
mkdir -p human_20416_truth
for i in `cat human_20416_list`;
do
    ../util/GFF_to_TransMemb human_20416_fasta/$i.fasta human_20416_gff/$i.gff > human_20416_truth/$i.top;
done
```

Note that the 'truth' here comes from UniProt Reviewed annotations, instead from experimental verified results, such as PDBTM.


## Step 5: extract ALL reviewed Human transmembrne proteins (TPs)
```
url="https://www.uniprot.org/uniprot/?query=annotation:(type:transmem)&fil=reviewed:yes+AND+organism:9606&sort=score&format=tab&columns=id,length"
wget -q $url -O human_5238_TP_len
tail -n+2 human_5238_TP_len | awk '{print $1}' > human_5238_TP_list
```

# Prediction results
Below we provide the prediction results of 2-state TP from different methods, where XX indicates human_20416:

| Folder name   | Description   | 
| ------------- | ------------- |
|   XX_phobius  | results from Phobius   |
|   XX_philius  | results from Philius   |
|   XX_topcons  | results from TOPCONS2 (Web Server) |
|   XX_purestm  | results from PureseqTM |


# Result evaluation
To evaluate the prediction results, please type the following command:
```
xx=human_20416
list=human_5238_TP_list
for method in purestm phobius philius topcons
do
	rm -f 5238_TP_${method}_reso
	for i in `cat $list`
	do
		../util/TM2_to_Evaluation ${xx}_truth/$i.top ${xx}_${method}/$i.top 2 >> 5238_TP_${method}_reso
	done
done
```


