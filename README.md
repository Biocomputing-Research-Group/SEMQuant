# SEMQuant

SEMQuant is a set of proteomics quantification pipelines developed by the
[Biocomputing Research Group](https://github.com/Biocomputing-Research-Group). The repository bundles
the tools used by the workflows — Raxport, Sipros Ensemble, Percolator, Philosopher, IonQuant and
DIA-NN — together with the helper scripts that glue them together.

Two pipelines are provided:

| Pipeline | Directory | Acquisition | Quantification |
| --- | --- | --- | --- |
| [SEMQuant-DDA](#semquant-dda) | `SEMQuant-DDA/` | Thermo DDA | IonQuant (MS1 LFQ + MBR) |
| [SEMQuant-Astral](#semquant-astral) | `SEMQuant-Astral/` | Orbitrap Astral | DIA-NN (spectral library) |

## Table of Contents

- [Repository layout](#repository-layout)
- [Requirements](#requirements)
- [Conventions](#conventions)
- [SEMQuant-DDA](#semquant-dda)
- [SEMQuant-Astral](#semquant-astral)
- [Citation](#citation)

## Repository layout

```
SEMQuant-DDA/
├── Raxport                     # Thermo .raw → FT2 converter (.NET / mono)
├── Sipros_OpenMP               # database search engine
├── SiprosEnsembleOMP           # Sipros Ensemble search engine
├── configs/                    # SiprosConfig*.cfg templates
├── Scripts/                    # Python/shell post-processing (filtering, assembly)
│   └── runSiprosFiltering.sh   # wrapper: tabulating → filtering → assembly → clustering
├── SE2IQ/                      # Sipros Ensemble → IonQuant/ProteinProphet converters
├── IonQuant-1.10.12.jar        # quantification (with jfreechart / batmass-io jars)
├── bruker/  thermo/            # vendor libraries required by IonQuant
├── yeast_ups.fasta             # example database (yeast + UPS1)
└── runSEMQuant.sh              # end-to-end example driver

SEMQuant-Astral/
├── percolator_3_6_4/           # PSM rescoring (linux/ and windows/ builds)
├── philosopher_v5.1.0_linux_amd64/   # protein inference + FDR control
├── linux_diann.zip             # DIA-NN 1.8.2 beta 8 (Git LFS; unzip before use)
├── script3/                    # Python 3 port of the Sipros post-processing scripts
└── SE_IQ_Astral/               # converters: Percolator/Sipros → ProteinProphet, IonQuant
```

## Requirements

Create the conda environments used throughout the workflows:

```bash
conda create -n py2  -c conda-forge scikit-learn python=2.7
conda create -n py3  -c conda-forge python=3.10
conda create -n mono -c conda-forge mono
```

- **Raxport** and **ThermoRawFileParser** require .NET / `mono`.
- **`SEMQuant-DDA/Scripts/`** targets Python 2.7 (`py2`); **`SEMQuant-Astral/script3/`** is the
  Python 3 port (`py3`).
- **IonQuant** requires a Java 11+ runtime.
- Some post-processing steps use **R**.
- `linux_diann.zip` is stored with **Git LFS** — run `git lfs install && git lfs pull` after cloning,
  then unzip it before running the DIA-NN step.

Tools **not** bundled in this repository and which must be obtained separately:
[ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser), and the FragPipe
`speclib/gen_con_spec_lib.py` spectral-library generator with its `unimod.xml`.

## Conventions

Every command below assumes these two variables are exported:

```bash
export work_dir=/PATH_TO_YOUR_WORK_DIR      # where raw data and results live
export tool_dir=/PATH_TO/SEMQuant           # this repository
```

Individual pipeline sections then use `$tool_dir/SEMQuant-DDA` or `$tool_dir/SEMQuant-Astral`.

---

## SEMQuant-DDA

DDA workflow: raw conversion → Sipros Ensemble database search → PSM filtering and protein
assembly → IonQuant label-free quantification with match-between-runs.

### 1. Set up the working directory

```bash
export work_dir=/PATH_TO_YOUR_WORK_DIR
export dda_dir=$tool_dir/SEMQuant-DDA

mkdir -p $work_dir/{raw,samples,results}
```

### 2. Download example raw files

```bash
cd $work_dir/raw
# yeast-UPS1 2 fmol replicates (PXD002099)
wget ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2020/01/PXD002099/*_2fmol*.raw
```

This yields three sample raw files.

### 3. Convert raw to FT2

```bash
conda activate mono
$dda_dir/Raxport -i $work_dir/raw -o $work_dir/raw
```

### 4. Generate the target–decoy database

```bash
conda activate py2
python $dda_dir/Scripts/sipros_prepare_protein_database.py \
    -i $dda_dir/yeast_ups.fasta \
    -o $work_dir/yeast_ups_rev.fasta \
    -c $dda_dir/configs/SiprosConfig_yeast.cfg
```

Replace the input/output FASTA with your own database, then point the configuration file at the
new decoy-appended database:

```ini
FASTA_Database = /PATH_TO_YOUR_WORK_DIR/yeast_ups_rev.fasta
```

### 5. Database search

`SiprosEnsembleOMP` searches either a single MS2 file (`-f`) or every MS2 file in a directory
(`-w`). Use `-g` instead of `-c` to sweep a directory of configuration files, and
`SiprosEnsembleOMP -h` for the full option list.

```bash
export OMP_NUM_THREADS=24

# single MS2 file
$dda_dir/SiprosEnsembleOMP -o $work_dir/results \
    -f $work_dir/raw/YOUR_MS2.FT2 \
    -c $dda_dir/configs/SiprosConfig_yeast.cfg

# all MS2 files in a directory
$dda_dir/SiprosEnsembleOMP -o $work_dir/results \
    -w $work_dir/raw \
    -c $dda_dir/configs/SiprosConfig_yeast.cfg
```

Search results are written to the output directory as `.Spe2Pep.txt` files.

### 6. Split results into per-sample folders

```bash
cd $work_dir
for sample in results/*Spe2Pep.txt; do
    name=$(basename "${sample%%.*}")
    mkdir -p samples/"$name"
    cp "$sample" samples/"$name"/
done
```

### 7. Filter PSMs and assemble proteins

`runSiprosFiltering.sh` chains PSM tabulating, ensemble filtering, protein assembly and SIP-mode
clustering:

```bash
conda activate py2
for folder in $work_dir/samples/*/; do
    $dda_dir/Scripts/runSiprosFiltering.sh \
        -in "$folder" \
        -c  $dda_dir/configs/SiprosConfig_yeast.cfg \
        -o  "$folder"
done
```

The individual steps can also be run directly:

```bash
for folder in $work_dir/samples/*/; do
    python $dda_dir/Scripts/sipros_ensemble_filtering.py \
        -i "$folder" -c $dda_dir/configs/SiprosConfig_yeast.cfg -o "$folder"
    python $dda_dir/Scripts/sipros_peptides_assembling.py \
        -w "$folder" -c $dda_dir/configs/SiprosConfig_yeast.cfg
done
```

### 8. Convert results for IonQuant

```bash
conda activate py3
for folder in $work_dir/samples/*/; do
    name=$(basename "$folder")
    python $dda_dir/SE2IQ/SE2prophet.py   "$folder/$name.psm.txt" "$folder" "$name"
    python $dda_dir/SE2IQ/SE2Ionquant.py  $work_dir/yeast_ups_rev.fasta "$folder"
done

# build the IonQuant file list, then copy the modification-mass list next to it
python $dda_dir/SE2IQ/generate_filelist.py
cp $dda_dir/modmasses_ionquant.txt $work_dir/samples/
```

`filelist_ionquant.txt` and `modmasses_ionquant.txt` in `SEMQuant-DDA/` are templates — edit
`--specdir` in the file list to point at your `raw/` directory.

### 9. Quantify with IonQuant

```bash
java -Xmx21G \
  -Dlibs.bruker.dir=$dda_dir/bruker \
  -Dlibs.thermo.dir=$dda_dir/thermo \
  -cp "$dda_dir/jfreechart-1.5.3.jar:$dda_dir/batmass-io-1.30.0.jar:$dda_dir/IonQuant-1.10.12.jar" \
  ionquant.IonQuant \
  --threads 23 --perform-ms1quant 1 --perform-isoquant 0 \
  --isotol 20.0 --isolevel 2 --isotype tmt10 --ionmobility 0 \
  --site-reports 1 --minexps 1 --mbr 1 --maxlfq 1 --requantify 1 \
  --mztol 10 --imtol 0.05 --rttol 0.4 \
  --mbrmincorr 0 --mbrrttol 1 --mbrimtol 0.05 --mbrtoprun 40 \
  --ionfdr 0.01 --proteinfdr 0.01 --peptidefdr 0.01 \
  --normalization 1 --minisotopes 2 --minscans 3 --minions 2 --minfreq 0 \
  --writeindex 0 --tp 0 --locprob 0.75 --uniqueness 0 --multidir . \
  --filelist $work_dir/samples/filelist_ionquant.txt \
  --modlist  $work_dir/samples/modmasses_ionquant.txt
```

`runSEMQuant.sh` in `SEMQuant-DDA/` shows the same sequence as a single driver script.

---

## SEMQuant-Astral

Astral workflow: raw conversion → Sipros database search → Percolator rescoring → ProteinProphet
inference → Philosopher FDR control → spectral library generation → DIA-NN quantification.

### 1. Set up the working directory

```bash
export work_dir=/PATH_TO_YOUR_WORK_DIR
export astral_dir=$tool_dir/SEMQuant-Astral

mkdir -p $work_dir/{raw,samples,results}
```

### 2. Download example raw files

```bash
cd $work_dir/raw
# three-species Astral benchmark, E45 condition (PXD046444)
wget ftp://ftp.ebi.ac.uk/pride/archive/projects/PXD046444/20230324_OLEP08_200ng_30min_E45H50Y5*.raw
```

This yields three sample raw files.

### 3. Convert raw to indexed mzML

```bash
conda activate mono
# -f=2 indexed mzML output, -L=2 keep MS2 only
mono ThermoRawFileParser.exe -d=$work_dir/raw -o=$work_dir/raw -f=2 -L=2
```

### 4. Database search

Uses the same search engine as the DDA pipeline (`$tool_dir/SEMQuant-DDA/Sipros_OpenMP`), on a
single MS2 file (`-f`) or a whole directory (`-w`):

```bash
export OMP_NUM_THREADS=24
export sipros=$tool_dir/SEMQuant-DDA/Sipros_OpenMP

# single MS2 file
$sipros -o $work_dir/results -f $work_dir/raw/YOUR_MS2.mzML \
        -c $tool_dir/SEMQuant-DDA/configs/SiprosConfig.cfg

# all MS2 files in a directory
$sipros -o $work_dir/results -w $work_dir/raw \
        -c $tool_dir/SEMQuant-DDA/configs/SiprosConfig.cfg
```

### 5. Generate PSM features for Percolator

An R script converts the Sipros PSM output into a Percolator `.pin` feature file. Run it once per
sample to produce `<sample>.pin`.

### 6. Rescore PSMs with Percolator

```bash
export sample=20230324_OLEP08_200ng_30min_E45H50Y5_180K_2Th3p5ms_02

$astral_dir/percolator_3_6_4/linux/percolator \
    --only-psms --no-terminate --post-processing-tdc --num-threads 23 \
    --results-psms       ${sample}_percolator_target_psms.tsv \
    --decoy-results-psms ${sample}_percolator_decoy_psms.tsv \
    --protein-decoy-pattern rev_ \
    ${sample}.pin
```

### 7. Convert Percolator output for protein inference

```bash
python $astral_dir/SE_IQ_Astral/Percolator2PeptideProphet.py \
    $work_dir/results/E45/01/${sample}.target.Spe2Pep.txt \
    $work_dir/results/E45/01/${sample}.pin \
    $work_dir/results/E45/01/ \
    ${sample}
```

This writes `interact-<sample>.pep.xml` files for ProteinProphet.

### 8. Protein inference

List the per-run `interact-*.pep.xml` files in a file list (see
`SEMQuant-Astral/filelist_proteinprophet.txt` for the expected format), then:

```bash
$astral_dir/philosopher_v5.1.0_linux_amd64/philosopher proteinprophet \
    --maxppmdiff 2000000 --output combined \
    $work_dir/E45/filelist_proteinprophet.txt
```

### 9. FDR control

```bash
export philosopher=$astral_dir/philosopher_v5.1.0_linux_amd64/philosopher
export db=$work_dir/2024-09-22-decoys-mix_HYE.fasta.fas

$philosopher database --annotate $db --prefix rev_

$philosopher filter --picked --razor --prot 0.01 --minPepLen 8 --tag rev_ \
    --pepxml  $work_dir/E45 \
    --protxml $work_dir/E45/combined.prot.xml
```

### 10. Generate the spectral library

Uses FragPipe's EasyPQP-based library builder (obtain `speclib/gen_con_spec_lib.py` and
`unimod.xml` from a FragPipe installation):

```bash
python -u $speclib_dir/gen_con_spec_lib.py \
    $db $work_dir/E45 unused $work_dir/E45 True unused use_easypqp noiRT;noIM 16 \
    "--unimod $speclib_dir/unimod.xml --max_delta_unimod 0.02 --max_delta_ppm 15.0 --fragment_types [\'b\',\'y\',]" \
    "--rt_lowess_fraction 0.0" \
    delete_intermediate_files \
    $work_dir/E45/filelist_speclibgen.txt
```

### 11. Quantify with DIA-NN

Unzip the bundled DIA-NN 1.8.2 beta 8 build first:

```bash
git lfs pull                                  # if not already fetched
unzip $astral_dir/linux_diann.zip -d $astral_dir/diann

$astral_dir/diann/diann-1.8.1.8 \
    --lib library.tsv --threads 15 --verbose 1 \
    --out diann-output/report.tsv \
    --qvalue 0.01 --matrix-qvalue 0.01 --matrices \
    --no-prot-inf --smart-profiling --no-quant-files \
    --peak-center --no-ifs-removal --report-lib-info \
    --cfg $work_dir/E45/filelist_diann.txt
```

---

## Citation

If you use SEMQuant in your research, please cite:

> Zhang B., Feng S., Xiong Y., Pan C., Guo X. SEMQuant: A computational platform for accurate and
> comprehensive quantitative metaproteomics analysis. *Journal of Computer Science and Technology*,
> 2026, 41(3): 1087–1100. https://doi.org/10.1007/s11390-026-5122-3

```bibtex
@article{zhang2026semquant,
  title   = {SEMQuant: A Computational Platform for Accurate and Comprehensive
             Quantitative Metaproteomics Analysis},
  author  = {Zhang, Bailu and Feng, Shichao and Xiong, Yi and Pan, Chongle and Guo, Xuan},
  journal = {Journal of Computer Science and Technology},
  volume  = {41},
  number  = {3},
  pages   = {1087--1100},
  year    = {2026},
  doi     = {10.1007/s11390-026-5122-3},
  url     = {https://doi.org/10.1007/s11390-026-5122-3}
}
```
