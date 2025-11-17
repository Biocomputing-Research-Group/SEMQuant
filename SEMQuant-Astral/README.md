# SEMQuant-Astral
SEMQuant-Astral

### Install environment

Please following the tutorial to install Sipors 5: https://github.com/xyz1396/sipros5
### Make folder for the workflow

```bash
#work_dir = YOUR_WORKING_DIRECTORY
export work_dir=/home/UNT/bz0053/Documents/Projects/test/
export tool_dir=/home/UNT/bz0053/Documents/Projects/SEMQuant/SEMQuant-Astral
# generate the following folder under your work_dir 
mkdir raw samples results 

### Download raw file

``` bash
cd raw 
# Download raw file with three-mixed species Astral data - E45 in your work_dir
wget ftp://www.ebi.ac.uk/pride/archive/projects/PXD046444/20230324_OLEP08_200ng_30min_E45H50Y5*.raw
```
The raw folder should contain 3 samples.

### Convert raw to FT2 use Raxport

```bash
conda activate mono
# -j is the threads that you plan to use
# -f set the spectra output format to indexed mzML
# -L select MS level to be MS2
mono ThermoRawFileParser.exe -d=./raw -o=./raw -f=2 -L=2 
```

### Running The Database-searching

There are two ways to run Sipros_OpemMP for database searching: one for running on a single MS2 file via `-f` and another another processing multiple MS2 files via `-w workingdirectory`.

```bash
#!/bin/bash

# Single MS2 file
Sipros -o output_dir -f ms_data -c SiprosConfig.cfg
$tool_dir/Sipros -o $work_dir/results -f $work_dir/raw/YOUR_MS2.FT2 -c $tool_dir/configs/SiprosConfig.cfg
# Multiple MS2 files in a working directory
Sipros -o output_dir -w work_dir -c SiprosConfig.cfg
$tool_dir/Sipros -o $work_dir/results -w $work_dir/raw -c $tool_dir/configs/SiprosConfig.cfg
```

Note that the source code for the Sipros-Astral mode is currently in-house and will be made available after publication via: https://github.com/xyz1396/sipros5.git. In the meantime, we provide a compiled binary named ```Sipros``` for immediate use.


### Use Rscript to generate PSM features for Percolator required input

### Percolator
```bash
Percolator [Work dir: /home/UNT/bz0053/Documents/Projects/MBR/dataset/Astral/fragpipe_MBR/E45/2]
$tool_dir/percolator --only-psms --no-terminate --post-processing-tdc --num-threads 23 --results-psms 20230324_OLEP08_200ng_30min_E45H50Y5_180K_2Th3p5ms_02_percolator_target_psms.tsv --decoy-results-psms 20230324_OLEP08_200ng_30min_E45H50Y5_180K_2Th3p5ms_02_percolator_decoy_psms.tsv --protein-decoy-pattern rev_ 20230324_OLEP08_200ng_30min_E45H50Y5_180K_2Th3p5ms_02.pin

```

### convert percolator output for protein inference input
```bash
python $tool_dir/Percolator2PeptideProphet.py ~/three-species/result_V2/SE_percolator_DIANN/E20/01/20230324_OLEP08_200ng_30min_E20H50Y30_180K_2Th3p5ms_01.target.Spe2Pep.txt ~/three-species/result_V2/SE_percolator_DIANN/E20/01/20230324_OLEP08_200ng_30min_E20H50Y30_180K_2Th3p5ms_01.pin ~/three-species/result_V2/SE_percolator_DIANN/E20/01/ 20230324_OLEP08_200ng_30min_E20H50Y30_180K_2Th3p5ms_01
```

### Protein inference
```bash
ProteinProphet [Work dir: /home/UNT/bz0053/Documents/Projects/MBR/dataset/Astral/fragpipe_MBR/E45]
$tool_dir/philosopher_v5.1.0_linux_amd64/philosopher proteinprophet --maxppmdiff 2000000 --output combined $work_dir/E45/filelist_proteinprophet.txt
```
### FDR control
```bash
#Process 'ProteinProphet' finished, exit code: 0
PhilosopherDbAnnotate [Work dir: F:\Astral\DIANN\E30]
$tool_dir\philosopher_v5.1.0_windows_amd64\philosopher.exe database --annotate F:\Astral\2024-09-22-decoys-mix_HYE.fasta.fas --prefix rev_
Process 'PhilosopherDbAnnotate' finished, exit code: 0

#PhilosopherFilter [Work dir: F:\Astral\DIANN\E30]
$tools_dir\philosopher_v5.1.0_windows_amd64\philosopher.exe filter --picked --prot 0.01 --minPepLen 8 --tag rev_ --pepxml F:\Astral\DIANN\E30 --protxml F:\Astral\DIANN\E30\combined.prot.xml --razor

```

### generate spec lib
```bash
SpecLibGen [Work dir: F:\Astral\DIANN\E30]
python -u $tool_dir\speclib\gen_con_spec_lib.py F:\Astral\2024-09-22-decoys-mix_HYE.fasta.fas F:\Astral\DIANN\E30 unused F:\Astral\DIANN\E30 True unused use_easypqp noiRT;noIM 16 "--unimod $tool_dir/unimod_old.xml --max_delta_unimod 0.02 --max_delta_ppm 15.0 --fragment_types [\'b\',\'y\',]" "--rt_lowess_fraction 0.0" delete_intermediate_files $work_dir\E45\filelist_speclibgen.txt

```

### DIANN for quantification
```bash
DIA-NN [Work dir: /home/UNT/bz0053/Documents/MBR/DIANN/E30]
$tool_dir/fragpipe/tools/diann/1.8.2_beta_8/linux/diann-1.8.1.8 --lib library.tsv --threads 15 --verbose 1 --out diann-output/report.tsv --qvalue 0.01 --matrix-qvalue 0.01 --matrices --no-prot-inf --smart-profiling --no-quant-files --peak-center --no-ifs-removal --report-lib-info --cfg $work_dir/E45/filelist_diann.txt--
DIA-NN 1.8.2 beta 8 (Data-Independent Acquisition by Neural Networks)
```
