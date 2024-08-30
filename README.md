### SEMQuant

This repository contains tools for SEMQuant developed by Biocomputing-Research-Group. These include Raxport, SiprosEnsemble and some other scripts.

[You can find the simple tutorial for download the yeast-ups1 2fmol on our wiki page](https://github.com/xyz1396/SiprosToolKits-Sipros4/wiki/13C-labeled-E.-coli-SIP-proteomic-search-tutorial)

### Install environment

```bash
conda create -n py2 scikit-learn python=2.7
conda create -n mono -c conda-forge mono
conda create -n py3 python=3.10
```
Raxport relies on .net. Some other scripts rely on python2, python3, and R.

### Make folder for the workflow

```bash
mkdir raw 
```

### Download raw file

```bash
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






### Citation

under review
