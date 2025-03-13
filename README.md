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
python sipros_prepare_protein_database.py -i /home/UNT/bz0053/Documents/Projects/MBR/dataset/yeast_ups/yeast_ups.fasta -o yeast_ups_rev.fasta -c ../configs/SiprosConfig_yeast.cfg

#need to update the FASTA_Database path in .cfg file 
FASTA_Database = /home/UNT/bz0053/Documents/Projects/MBR/dataset/mock_uneven/Mock_Comm_RefDB_V3_rev.fasta
```
The step will generate a new database file with reverse sequences. Update the path of `FASTA_Database` in the configuration file.

### Running The Database-searching

There are two two to run Sipros_OpemMP for database searching: one for running on a single MS2 file via `-f` and another another processing multiple MS2 files via `-w workingdirectory`.


### Citation

under review
```bash
#!/bin/bash

# Single MS2 file
Sipros_OpemMP -o output_dir -f ms_data -c SiprosConfig.cfg

# Multiple MS2 files in a working directory
Sipros_OpemMP -o output_dir -w workingdirectory -c SiprosConfig.cfg

```

Results (`.Spe2Pep` files) will be saved on the output directory. if you have many configure files, specify `-g`, like `Sipros_OpemMP -o output_dir -w workingdirectory -g configurefiledirectory`. Use `./Sipros_OpemMP -h` for help information. 

