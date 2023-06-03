# Building a Profile Hidden Markov Model for the Kunitz-type protease inhibitor domain
This repository contains the files generated to complete the final project for module 2 of the Laboratory of Bioinformatics course of the Master's Degree in Bioinformatics of the Università di Bologna, course 2022-2023.

This project consists of building a Hidden Markov Model profile for the Kunitz-type protease inhibitor domain, starting from available structural information. Once the model is built, it is to be used to annotate Kunitz domain(s) in UniProtKB/Swiss-Prot proteins.

For details, see the written report (`project-lb1b-Laia_Torres_Masdeu.pdf`).

Written by Laia Torres Masdeu.

**Table of Contents**
-  [0. Requirements and databases](#0-requirements-and-databases)
- [1. Structure Selection](#1-structure-selection)
    - [1.1. PDB Search](#11-pdb-search)
    - [1.2. Protein Alignment](#12-protein-alignment)
- [2. HMM Generation](#2-hmm-generation)
    - [2.1. Train a HMM profile](#21-train-a-hmm-profile)
    - [2.2. Training set](#22-training-set)
- [3. Method Testing](#3-method-testing)
    - [3.1. Validation set](#31-validation-set)
    - [3.2. Cross-validation and Testing sets](#32-cross-validation-and-testing-sets)
    - [3.3. Model evaluation](#33-model-evaluation)

## 0. Requirements and used versions and releases
To be able to conduct this project, download the HMMER package from their [website](http://hmmer.org/). This project was conducted using HMMER v3.3.2 (Nov 2020) version. 

All UniProt searches were performed using Release 2023_02.
All PDB searches were conducted on 27/05/2023.
All PDBeFold submissions were issued under PDBe Fold v2.59.

## 1. Structure selection
The goal of this first section is to retrieve any available structures that have been annotated with the Kunitz domain from the [PDB database](https://www.rcsb.org/) and to generate a Mulitple Sequence Alignment (MSA) from these through [PDBeFold](https://www.ebi.ac.uk/msd-srv/ssm/) web-based tool.

### 1.1. PDB Search
First, access the [Advanced search](https://www.rcsb.org/search/advanced) section of the PDB website, and write the desired query. In my case, I opted for these parameters:
- Pfam identifier, SCOP2 lineage identifier, or CATH lineage identifier (`PF00014`, `4003337` or `4.10.410.10`): to select those structurally-resolved PDB structures that have been annotated as containing a Kunitz domain.
- Sequences with no mutations (`Polymer Entity Mutation Count=0`): wild-type versions of the protein, no mutants.
- Resolution (`<= 3`).
- Polymer Entity Sequence Length (`51 - 76` residues): size range of the Kunitz domains.

This query returned >100 structures. However, it is known that structures belonging to the same protein can be found in different PDB files. Therefore, group those structures that belong to the same protein, and get a representative for each group. To do this, select the `Return Polymer Entities` that are `grouped by Sequence Identity 95%` and `displaying as Representatives`. 

After retrieving these entities and before moving on to the protein alignment, generate a Custom Tabular Report. Select:
- `PDB ID`.
- `Auth Asym ID`: chain ID given by the author.
- `Accession code(s)`: UniProtIDs related to each PDB entry.

The first two will be necessary to conduct the protein alignment, while the last will be necessary to make the positive test set.

Once all of the desired data items are selected, download the tabular report in CSV format (`pdb_report.csv`). 

### 1.2. Protein Alignment
To be able to conduct the MSA, create a text file (`pdb_codes.txt`) that contains the PDB identifiers and desired chains in the PDBeFold format (3abc:A):
```
grep -v 'Identifier\|Entity\|,,' pdb_report.csv | cut -d ',' -f 2,3 | tr -d \" |tr "," ":" >pdb_codes.txt
```
- `grep -v`: select those lines in the input file that do <u>not</u> contain the subsequent terms. `-i` enables case insensitive matching.
- `cut`: cut out selected portions of each line of the input file, with the desired delimiter character `-d`, and output the selected fields `-f`
- `tr`: substitute (translate) first character by second character. `-d` option deletes that character from the input file
- `>`: redirect final output to the output file

Lastly, access the PDBeFold website:
1. Click on `Launch PDBeFold`
2. Select, inside `Multiple Submission Form`, `List of PDB codes` as the source
3. Upload the file (`pdb_codes.txt`) 
4. Once all entries are displayed, click on `Update List` (with Jmol as the viewer)
5. Submit the query

Once the request has been completed, download the MSA in FASTA format (`msa.fasta`).

&emsp;<u>Note</u>: To check that the MSA is correct, download the Multiple Alignment Results (`msa.dat`). All sequences should have an RMSD <2. If this is not the case, remove those entries (`grep -v`) from the PDB codes file (`pdb_codes_ref.txt`), re-run the MSA, and download the new file (`msa_ref.fasta`).

## 2. HMM Generation
The goal of this second section is to generate the HMM profile from the selected PDB entities and to verify its performance on a training set.

### 2.1. Train an HMM profile
Once all the MSA sequences have an RMSD <2, build the HMM profile. 

&emsp;<u>Optional</u>: If the MSA has poorly aligned or highly gapped regions that could influence the resulting HMM profile, trim these areas (being `pos1` and `pos2` the initial and final position of the residues to include, respectively):
```
grep . msa.fasta |awk '{if (substr($1,1,1)==">") {printf "\n%s ",$1} else {printf "%s",$0}}'|awk '{print $1; print toupper(substr($2,pos1,pos2))}'>kunitz_mod.fasta
```
- `grep .`: remove empty lines
- `awk`: used to clean the fasta file, so that it just contains the PDB ID and chain, followed by the sequence, in uppercases, from `pos1` to `pos2`.

&emsp;<u>Note</u>: After trimming, check that there are no duplicated sequences within the MSA. If there are, remove one of the sequences (duplicate sequences can introduce redundancy and potentially bias the model training process). To do so, upload the trimmed file (`kunitz_mod.fasta`) to the [Align](https://www.uniprot.org/aligng) tool at the UniProt website, and delete manually (using vim) the FASTA entry and sequence of one of the duplicated sequences.

Once the (trimmed) MSA is ready, build the HMM profile (where `kunitz_mod.hmm` will be the output HMM profile and `kunitz_mod.fasta` is the input MSA):
```
hmmbuild kunitz_mod.hmm kunitz_mod.fasta
```

### 2.2. Training set
To verify that the trained HMM is able to recognise the proteins in the dataset (consistency test), perform a `hmmsearch` with the sequences of the proteins used to generate the model as the sequence database (second argument). 

The FASTA sequences can be retrieved with the [Retrieve/ID Mapping](https://www.uniprot.org/id-mapping) tool at the UniProt website. The IDs to retrieve can be loaded from a text file. 

Therefore, create a text file (`training_set.idlist`) that contains the (unique) UniProt identifiers of the model proteins:
```
grep -v 'Identifier\|Entity\|,,' pdb_report.csv | cut -d ',' -f 4 | tr -d \" |sort |uniq >training_set.idlist
```
- `sort`: sort input file in ascending order
- `uniq`: filter out repeated lines in the input file (file must be sorted)

&emsp;<u>Note</u>: Remember to remove (if any) entries that scored RMSD <2 when performing the MSA, and that finally are not included in the model. 

Map the IDs in the UniProt website to retrieve their sequences, and download them, compressed, as FASTA (canonical) format (`training_set.fasta`). Finally, run `hmmsearch` (Usage: `hmmsearch [options] <hmmfile> <seqdb>`):
```
hmmsearch --max -o training_set.out kunitz_mod.hmm training_set.fasta
```
- `--max`: turn off all the heuristics for cutting off distantly related proteins
- `-o`: name of the output file

Check the output file on the training set, evaluating the E-value, alignment quality, consistency, and coverage. These should all be very high, as we are running the model over the same set of sequences used to generate it.

## 3. Method Testing
The goal of this third and last section is to retrieve a suitable validation set, consisting of known proteins that contain (positive set) or not (negative set) the Kunitz domain, search it against the trained model, and compute the scoring indexes to evaluate the HMM profile on cross-validation and testing sets.

### 3.1. Validation set
The FASTA sequences for the validation set can be retrieved with the Advanced search option at the [UniProt website](https://www.uniprot.org/), restricting the search to UniProtKB/Swiss-Prot entries (reviewed proteins). To retrieve those proteins that are known to contain a Kunitz domain, search for the Pfam identifier `PF00014` in the Advanced Search. This search will contain the proteins that will make up the **positive set** (`/UniProt_Swiss-Prot_data/uniprot_sp_kunitz.fasta`). All the remaining proteins in UniProtKB/Swiss-Prot will be used as the **negative set** (`/UniProt_Swiss-Prot_data/uniprot_sp_nonkunitz.fasta`).

Once both datasets are downloaded, create a file that contains the UniProtIDs for each set:
```
grep '^>' uniprot_sp_kunitz.fasta |cut -d '|' -f 2 >uniprot_sp_kunitz.idlist
grep '^>' uniprot_sp_nonkunitz.fasta |cut -d '|' -f 2 >uniprot_sp_nonkunitz.idlist
```

For a fair test on the model, the positive test set should exclude the training data. 
To remove those sequences that are part of the positive set and that have been used in the training test, create a file (`pos_set.idlist`) that doesn't include the UniProtIDs of the proteins used to build the model:
```
comm -13 <(sort training_set.idlist) <(sort ../UniProt_Swiss-Prot_data/uniprot_sp_kunitz.idlist) >pos_set.idlist
```
- `comm`: compare two files and show the lines that are unique or in common. The first column shows lines unique to the first file, the second column shows the lines unique to the second file and the third column shows lines common in both. `-13` returns only the unique lines found in the second file.

Once the UniprotIDs of the proteins that have <u>not</u> been used in the training set are obtained, go to the [Retrieve/ID Mapping](https://www.uniprot.org/id-mapping) tool at the UniProt website to retrieve and download the FASTA file (`pos_set.fasta`).

Next, create a FASTA file (`val_set.fasta`) that contains all of the sequences (in FASTA format) of the validation set (positive and negative set):
```
cat pos_set.fasta ../UniProt_Swiss-Prot_data/uniprot_sp_nonkunitz.fasta >val_set.fasta
```

Run `hmmsearch`, with the file that was just generated (`val_set.fasta`) as a sequence database:
```
hmmsearch --max --noali -o val_set.search kunitz_mod.hmm val_set.fasta 
```
- `--noali`: exclude the alignments from the output (to simplify and ease output visualisation)

Refine the output to just include the hits that contain the E-values, UniProtIDs, and other information (`final_pos` corresponds to the number of the last row containing this information):
```
head -n final_pos val_set.search |tail -n +18|grep -v inclusion >val_set.out
```
- `head -n [x]`: display the first x lines of the input file
- `tail -n +[x]`: display the lines starting from line x of the input file, till the end

Now that the hits are refined, generate a file for each sub-validating set () that will be used later for the testing, in the format `UniProtID E-value class`, where `E-value` is the E-value for the highest scoring hit (domain) and the `class` is the actual ("real") annotation, obtained from UniProt: 1 for proteins that contain a Kunitz domain and 0 for proteins that don't contain a Kunitz domain.

```
awk '{split($9,i,"|"); print i[2],$4}' val_set.out |grep -f pos_set.idlist |awk '{print $1, $2, 1}' >val_set.classp
awk '{split($9,i,"|"); print i[2],$4}' val_set.out |grep -v -f pos_set.idlist |awk '{print $1, $2, 0}' >val_set.classn
```

Finally, add all the proteins retrieved from UniProt as non-Kunitz to the class file. To do so, first, get the UniProtIDs that have been a hit against the model:
```
cut -d ' ' -f 1 val_set.classn >val_set.idn 
```

And lastly, retrieve those IDs that did not appear in the `hmmsearch` output and add them to the class file (with a high E-value, eg `100`, to evidence that they did not match the model):
```
comm -13 <(sort val_set.idn) <(sort ../UniProt_Swiss-Prot_data/uniprot_sp_nonkunitz.idlist) |awk '{print $0,100,0}' >> val_set.classn
```

### 3.2. Cross-validation and Testing sets
The testing procedure consists of the implementation of a 2-fold cross-validation test: the positive and negative sets are split into two subsets (`set1.txt` and `set2.txt`), optimising the classification (E-value) threshold on one subset (cross-validation set) and testing the performance on the other subset (testing set), and reversibly.

Before dividing the positive and training sets, shuffle the content of the files, and save it in a new file (`val_set.randomp` and `val_set.randomn`):

```
sort -R val_set.classp > val_set.randomp
sort -R val_set.classn > val_set.randomn
```
- `sort -R`: sort the file randomly

To divide the sets into two different files, count the number of entries/lines (`wc -l`) there are in each of the `val_set.randomp` and `val_set.randomn` files, add it to a variable (`n_p`, number of positives, `n_n`, number of negatives) and half these variables to obtain the number of entries in half of each set (`n_p_h`, half number of positives, `n_n_h`, half number of negatives):
```
n_p=$(wc -l val_set.randomp |awk '{print $1}')
n_p_h=$(($n_p/2))
n_n=$(wc -l val_set.randomn |awk '{print $1}')
n_n_h=$(($n_n/2))
```

Finally, select the first half of each of the `val_set.randomp` and `val_set.randomn` files, and add it to the first subset (`set1.txt`) and the second half of each file and add it to the second subset (`set2.txt`):
```
head -n $n_n_h val_set.randomn >set1.txt
head -n $n_p_h val_set.randomp >>set1.txt
tail -n +$(($n_n_h+1)) val_set.randomn >set2.txt
tail -n +$(($n_p_h+1)) val_set.randomp >>set2.txt
```

### 3.3. Model evaluation
Compute the scoring indexes for evaluating the profile HMM on the validation sets. This is done in 5 sub-steps:
1. Calculate the optimised E-value on the first dataset (`set1.txt` = cross-validation set)
2. Use this E-value on the second dataset and evaluate its performances (`set2.txt` = testing set)
3. Calculate the optimised E-value on the second dataset (`set2.txt` = cross-validation set)
4. Use this E-value on the first dataset and evaluate its performances (`set1.txt` = testing set)
5. Compute the average E-value and get the performances on the whole dataset (`val_set.txt` = testing set)

To find the optimised E-values and their performances - Matthews correlation coefficient (MCC), accuracy (ACC), and confusion matrix (CM) - run the python script `confusion_matrix.py`, which takes as first argument the file to be evaluated and as second argument the E-value threshold to use (where `eval2` and `eval1` are the optimised E-values obtained from the cross-validation sets 1 and 2, respectively). The output is saved in a file (`.cvs` stores the results of the cross-validation sets and `.tv` of the testing set):
```
for i in `seq 1 20`; do python ../confusion_matrix.py set1.txt 1e-$i; done >set1.cvs
for i in `seq 1 20`; do python ../confusion_matrix.py set2.txt 1e-$i; done >set2.cvs
python ../confusion_matrix.py set1.txt eval2 > set1.ts
python ../confusion_matrix.py set2.txt eval1 > set2.ts
```

Finally, create the text file that will contain the full validation set (`val_set.txt`), to compute the performance on the whole dataset: 
```
cat val_set.randomp val_set.randomn > val_set.txt
```

And run the Python script with the average E-value (`evalavg`), and save the output on a file:
```
python ../confusion_matrix.py val_set.txt evalavg >val_set.ts
```
