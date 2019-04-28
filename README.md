# Dataset for PureseqTM
The datasets for training and testing PureseqTM can be found here. In addition, the Human proteome data is provided. For each dataset, we also generated the prediction results from [Phobius](http://phobius.sbc.su.se/), [Philius](http://www.yeastrc.org/philius/pages/philius/runPhilius.jsp), [Topcons2](http://topcons.cbr.su.se/), and [PureseqTM](http://pureseqtm.predmp.com/).


| Folder name   | Description   |
| ------------- | ------------- |
| pdbtm_database        | ground-truth labeling from [PDBTM](http://pdbtm.enzim.hu/) database. |
| Train_and_Test_result | training and testing set, also include the prediction results from the four methods. |
| TP_and_nonTP_result   | dataset for discrinating Transmembrane Proteins (TPs) and non-TPs. |
| Human_proteome_result | human proteome dataset from [UniProt](https://www.uniprot.org/). |
| source_code,util      | source code for evaluation and label generation.  |


# UniProt Human membrane proteome
Given a UniProt ID (e.g., [Q9UMS5](https://www.uniprot.org/uniprot/Q9UMS5)), users may find the prediction result of PureseqTM by the following link:
http://pureseqtm.predmp.com/view.html?id=Q9UMS5_PureTM&name=Q9UMS5
![alt text](http://18.209.146.171:8081/site_media/data/Q9UMS5_PureTM/Q9UMS5.png "Q9UMS5 result")

