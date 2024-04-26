#!/bin/bash
eval "$(conda shell.bash hook)"

rawfile_folder=""
configfile=""
Working_Directory=""
Number_threads=""
FASTA_Database="" #yeast_ups1.fasta
Ionquant_dir=""

while [[ $# >0 ]]
do
	key="$1"
	case $key in
		-h|--help)
		echo -e "Usage:\n"
		echo -e "inference.sh [OPTION}...<PARAM>...\n\n"
		echo -e "	-r\t raw file folder for experimental mass spectrum\n"
		echo -e "	-w\t working directory for Sipros-Ensemble\n"
		echo -e "	-c\t Configuration file\n"
		echo -e "   	-t\t Number of threads\n"
		exit 1
		;;
		-r)
		rawfile_folder="$2"
		shift
		;;
		-c)
		config_file="$2"
		shift
		;;
		-w)
		Working_Directory="$2"
		shift
		;;
		-t)
		Number_threads="$2"
		shift
		;;
		-f)
		FASTA_Database="$2"
		shift
		;;
		*)
		echo "ERROR: Unidentified user variable $key"
		exit 1
		;;
	esac
	shift
done


conda activate mono
# generate FT1 FT2 files
mono Raxport.exe -i ${rawfile_folder} -o ${Working_Directory} -j ${Number_threads}
conda deactivate
# generate rev_fasta
conda activate py2
python ./Scripts/sipros_prepare_protein_database.py -i ${FASTA_Database} -o ${FASTA_Database%.*}_rev.fasta -c ${config_file}

#  change origional DB to rev_DB 
conda deactivate
conda activate py3
python SE2IQ/UpdateConfig.py ${config_file} ${FASTA_Database}
conda deactivate

export OMP_NUM_THREADS=${Number_threads}
mkdir ${Working_Directory}/result 
./SiprosEnsembleOMP -w ${Working_Directory} -c ${config_file} -o ${Working_Directory}/result

for sample in ${Working_Directory}/result/*Spe2Pep.txt; do path=${sample%.SE*}; name=${path#result/}; mkdir ${name}; done
for sample in ${Working_Directory}/result/*Spe2Pep.txt; do path=${sample%.SE*}; name=${path#result/}; mv ${sample} ${name}/${sample#${Working_Directory}/result/}; done

### start here
conda activate py2
for folder in ${Working_Directory}/result/*/; do ./Scripts/runSiprosFiltering.sh  -in ${folder} -c ${config_file} -o ${folder}; done
conda deactivate

# generate pep.xml file
for folder in ${Working_Directory}/result/*; do python SE2IQ/SE2prophet.py ${folder}/${folder##*result/}.SE.psm.txt ${folder} ${folder##*result/}; done   
# generate protein.tsv and psm.tsv
for folder in ${Working_Directory}/result/*; do python SE2IQ/SE2Ionquant.py /${FASTA_Database%.*}_rev.fasta ${folder}; done

##IQ
#create filelist use generate.py and copy modmasses file to sample file
conda activate py3
cd 
python SE2IQ/generate_filelist.py ${config_file} ${FASTA_Database}

# generate filelist_ionquant.txt
python ./SE2IQ/generate_filelist.py ${Working_Directory}
cp ./modmasses_ionquant.txt ${Working_Directory}/result
conda deactivate


# IonQuant
# Need to download IonQuant 1.10.12 from https://msfragger.arsci.com/ionquant/
conda activate mono


java -Xmx23G -Dlibs.bruker.dir=/home/UNT/bz0053/Documents/Projects/fragpipe/tools/MSFragger-4.0/ext/bruker \
     -Dlibs.thermo.dir=/home/UNT/bz0053/Documents/Projects/fragpipe/tools/MSFragger-4.0/ext/thermo \
     -cp /home/UNT/bz0053/Documents/Projects/fragpipe/tools/jfreechart-1.5.3.jar:/home/UNT/bz0053/Documents/Projects/fragpipe/tools/batmass-io-1.30.0.jar:/home/UNT/bz0053/Documents/Projects/fragpipe/tools/IonQuant-1.10.12.jar ionquant.IonQuant \
     --threads 23 --perform-ms1quant 1 --perform-isoquant 0 --isotol 20.0 --isolevel 2 --isotype tmt10 --ionmobility 1 --site-reports 1 --minexps 1 \
     --mbr 1 --maxlfq 1 --requantify 1 --mztol 10 --imtol 0.05 --rttol 0.4 --mbrmincorr 0 --mbrrttol 1 --mbrimtol 0.05 --mbrtoprun 3 --ionfdr 0.01 --proteinfdr 1 \
     --peptidefdr 1 --normalization 1 --minisotopes 2 --minscans 3 --writeindex 0 --tp 0 --minfreq 0 --minions 2 --locprob 0.75 --uniqueness 0 --multidir . \
     --filelist ${Working_Directory}/result/filelist_ionquant.txt \
     --modlist ${Working_Directory}/result/modmasses_ionquant.txt



#run in the SEMQuant folder
#wd=/home/UNT/bz0053/Documents/Projects/MBR/SEMQuant

rawfile_folder=/home/UNT/bz0053/Documents/Projects/test/raw
Working_Directory=/home/UNT/bz0053/Documents/Projects/test
Number_threads=23
FASTA_Database=/home/UNT/bz0053/Documents/Projects/test/yeast_ups.fasta
config_file=/home/UNT/bz0053/Documents/Projects/MBR/SEMQuant/SiprosConfigBenchmark.cfg



for folder in ${Working_Directory}/result/*; do echo ${folder}; done