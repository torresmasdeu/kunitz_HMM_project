# 1. Structure selection
grep -v 'Identifier\|Entity\|,,' pdb_report.csv |cut -d ',' -f 2,3 |tr -d \" |tr "," ":" >pdb_codes.txt

# 2. HMM Generation
# 2.1. Train a HMM profile
grep . msa.fasta |awk '{if (substr($1,1,1)==">") {printf "\n%s ",$1} else {printf "%s",$0}}'|awk '{print $1; print toupper(substr($2,13,54))}'>kunitz_mod.fasta
vim kunitz_mod.fasta #to remove 1YC0 sequence
hmmbuild kunitz_mod.hmm kunitz_mod.fasta
# 2.2. Training set
grep -v -i 'Identifier\|Entity\|,,' pdb_report.csv |cut -d ',' -f 4 |tr -d \" |sort |uniq >training_set.idlist 
hmmsearch --max -o training_set.out kunitz_mod.hmm training_set.fasta

# 3. Method Testing
# 3.1. Validation dataset
comm -13 <(sort training_set.idlist) <(sort ../UniProt_Swiss-Prot_data/uniprot_sp_kunitz.idlist) >pos_set.idlist
cat pos_set.fasta ../UniProt_Swiss-Prot_data/uniprot_sp_nonkunitz.fasta >val_set.fasta
hmmsearch --max --noali -o val_set.search kunitz_mod.hmm val_set.fasta 
head -n 409 val_set.search |tail -n +18|grep -v inclusion >val_set.out
awk '{split($9,i,"|"); print i[2],$4}' val_set.out |grep -f pos_set.idlist |awk '{print $1, $2, 1}' >val_set.classp
awk '{split($9,i,"|"); print i[2],$4}' val_set.out |grep -v -f pos_set.idlist |awk '{print $1, $2, 0}' >val_set.classn
cut -d ' ' -f 1 val_set.classn >val_set.idn
comm -13 <(sort val_set.idn) <(sort ../UniProt_Swiss-Prot_data/uniprot_sp_nonkunitz.idlist) |awk '{print $0,100,0}' >> val_set.classn
# 3.2. Cross-validation and Testing sets
sort -R val_set.classp > val_set.randomp
sort -R val_set.classn > val_set.randomn
n_p=$(wc -l val_set.randomp |awk '{print $1}')
n_p_h=$(($n_p/2))
n_n=$(wc -l val_set.randomn |awk '{print $1}')
n_n_h=$(($n_n/2))
head -n $n_n_h val_set.randomn >set1.txt
head -n $n_p_h val_set.randomp >>set1.txt
tail -n +$(($n_n_h+1)) val_set.randomn >set2.txt
tail -n +$(($n_p_h+1)) val_set.randomp >>set2.txt
# 3.3. Model evaluation
for i in `seq 1 20`; do python ../confusion_matrix.py set1.txt 1e-$i; done >set1.cvs
for i in `seq 1 20`; do python ../confusion_matrix.py set2.txt 1e-$i; done >set2.cvs
python ../confusion_matrix.py set1.txt 1e-5 > set1.ts
python ../confusion_matrix.py set2.txt 1e-3 > set2.ts
cat val_set.randomp val_set.randomn > val_set.txt
python ../confusion_matrix.py val_set.txt 5.05e-4 >val_set.ts