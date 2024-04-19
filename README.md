### SEMQuant

This repository contains tools for EMQuant developed by Biocomputing-Research-Group.  These include Raxport, SiprosEnsemble and some scripts.

[You can find the simple tutorial for download the yeast-ups1 2fmol on our wiki page](https://github.comxyz1396/SiprosToolKitswiki/13C-labeled-E.-coli-SIP-proteomic-search-tutorial)

### Citation

under review

### Install environment

# Raxport relies on .net. Some scripts rely on python2 and R.

```bash
conda create -n py2 scikit-learn python=2.7
conda create -n mono -c conda-forge mono
conda create -n r -c conda-forge -c bioconda r-base r-stringr r-tidyr bioconductor-biostrings
```

### Make folder for the workflow

mkdir raw 

### Download raw file

```bash
cd raw 
# Download raw file with 2fmol yeast-ups1
wget ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2020/01/PXD002099/*_2fmol*.raw
```

### Download Sipros-Ensamble 





