# SEMQuant

This repository contains tools for SEMQuant developed by Biocomputing-Research-Group. These include Raxport, SiprosEnsemble and some other scripts.

[You can find the simple tutorial for download the yeast-ups1 2fmol on our wiki page](https://github.com/xyz1396/SiprosToolKits-Sipros4/wiki/13C-labeled-E.-coli-SIP-proteomic-search-tutorial)

### Install environment


This repository contains tools for SEMQuant developed by Biocomputing-Research-Group. These include Raxport, SiprosEnsemble, Sipros5 and some scripts. there are two pipeline for SEMQuant: SEMQUant-DDA and SEMQuant-DIA.

## SEMQuant-DDA
SEMQuant-DDA 

### Citation

under review

### Install environment

#### Raxport relies on .net. Some scripts rely on python2, python3, and R.

``` bash
conda create -n py2 scikit-learn python=2.7
conda create -n mono -c conda-forge mono
conda create -n py3 python=3.10
```
Raxport relies on .net. Some other scripts rely on python2, python3, and R.

### Make folder for the workflow

```bash
mkdir raw samples result

### Download raw file

``` bash
cd raw 
# Download raw file with 2fmol yeast-ups1
wget ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2020/01/PXD002099/*_2fmol*.raw
```
The raw folder should contain 3 samples.

### Convert raw to mzml

```bash
conda activate mono
# -j is the threads that you plan to use
# -f set the spectra output format to indexed mzML
# -L select MS level to be MS2
mono ThermoRawFileParser.exe -d=./raw -o=./raw -f=2 -L=2 
```
### Generate Reverse Sequences

``` bash
python sipros_prepare_protein_database.py -i work_dir/db.fasta -o db_rev.fasta -c work_dir/configs/SiprosConfig_db.cfg

#need to update the FASTA_Database path in .cfg file 
FASTA_Database = work_dir/db_rev.fasta
```
The step will generate a new database file with reverse sequences. Update the path of `FASTA_Database` in the configuration file.

### Running The Database-searching

There are two two to run Sipros_OpemMP for database searching: one for running on a single MS2 file via `-f` and another another processing multiple MS2 files via `-w workingdirectory`.


```bash
#!/bin/bash

# Single MS2 file
Sipros_OpemMP -o output_dir -f ms_data -c SiprosConfig.cfg

# Multiple MS2 files in a working directory
Sipros_OpemMP -o output_dir -w work_dir -c SiprosConfig.cfg

```

Results (`.Spe2Pep` files) will be saved on the output directory. if you have many configure files, specify `-g`, like `Sipros_OpemMP -o output_dir -w workingdirectory -g configurefiledirectory`. Use `./Sipros_OpemMP -h` for help information. 

```bash
conda activate py2

#prep sample for store results
for sample in result/E10/*Spe2Pep.txt; do path=${sample%.SE*}; name=${path#result/E10}; mkdir samples_E10/${name}; done
for sample in result/E10/*Spe2Pep.txt; do path=${sample%.SE*}; name=${path#result/E10/}; cp ${sample} samples_E5/${name}/${sample#result/E5/}; done

for folder in /home/UNT/fs0199/three-species/samples_E5/*/; do echo ${folder}; done

### filter PSMs
for folder in /home/UNT/fs0199/three-species/samples_E5/*/; do ./runSiprosFiltering.sh  -in ${folder} -c ../configs/SiprosConfig_yeast.cfg -o ${folder}; mv ${folder}/*.p*.txt ${folder}/peptide/;mv ${folder}/*.tab ${folder}/peptide/; done
for folder in /home/UNT/bz0053/Documents/Projects/MBR/dataset/yeast_ups/samples_2fmol/*/; do ./runSiprosFiltering.sh  -in ${folder} -c ../configs/SiprosConfig_yeast.cfg -o ${folder}; done

# filter for pep 1 %
mkdir PEP_1pct
cp *peptide.tsv PEP_1pct/
for folder in /home/UNT/fs0199/three-species/samples_E10/20*/; do cp ${folder}/*.SE.pep.txt ./PEP_1pct; done
ls -la PEP_1pct/

for folder in /home/UNT/fs0199/three-species/samples_E10/*_01/ ;do python /home/UNT/fs0199/sipros5-master/script3/sipros_ensemble_filtering.py -i ${folder} -c /home/UNT/fs0199/three-species/SiprosConfig.cfg -o ${folder}; done


for folder in /home/UNT/fs0199/three-species/samples_E10/*_01/; do python /home/UNT/fs0199/sipros5-master/script3/sipros_peptides_assembling.py  -c /home/UNT/fs0199/three-species/SiprosConfig.cfg -w ${folder}; donehttps://chatgpt.com/c/67d229a2-725c-8001-bd57-48f02bf9d8e6

```
### IonQuant for quantification

```bash
conda activate mono
java -Xmx21G -Dlibs.bruker.dir=/home/UNT/bz0053/Documents/fragpipe/tools/MSFragger-4.0/ext/bruker -Dlibs.thermo.dir=/home/UNT/bz0053/Documents/fragpipe/tools/MSFragger-4.0/ext/thermo -cp /home/UNT/bz0053/Documents/fragpipe/tools/jfreechart-1.5.3.jar:/home/UNT/bz0053/Documents/fragpipe/tools/batmass-io-1.30.0.jar:/home/UNT/bz0053/Documents/fragpipe/tools/IonQuant-1.10.12.jar ionquant.IonQuant --threads 23 --perform-ms1quant 1 --perform-isoquant 0 --isotol 20.0 --isolevel 2 --isotype tmt10 --ionmobility 0 --site-reports 1 --minexps 1 --mbr 1 --maxlfq 1 --requantify 1 --mztol 10 --imtol 0.05 --rttol 0.4 --mbrmincorr 0 --mbrrttol 1 --mbrimtol 0.05 --mbrtoprun 40 --ionfdr 0.01 --proteinfdr 0.01 --peptidefdr 0.01 --normalization 1 --minisotopes 2 --minscans 3 --writeindex 0 --tp 0 --minfreq 0 --minions 2 --locprob 0.75 --uniqueness 0 --multidir . --filelist /home/UNT/bz0053/Documents/Projects/MBR/dataset/Astral/fragpipe_MBR/E45/filelist_ionquant.txt --modlist /home/UNT/bz0053/Documents/Projects/MBR/dataset/Astral/fragpipe_MBR/E45/modmasses_ionquant.txt

```


### Citation

under review
