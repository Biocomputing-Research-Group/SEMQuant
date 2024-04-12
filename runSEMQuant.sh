#!/bin/bash 

#-------------------------------------------------------------------#

# Get path of run script
if [ -L $0 ] ; then
    exePath=$(dirname $(readlink -f $0)) ;
else
    exePath=$(dirname $0) ;
fi;

SEMQuanFolder=""

SipFolder=""
ConfigureFile=""
OutputFolder=""
TabFile=""

# Read Parameter
while [[ $# > 0 ]]
do
key="$1"
case $key in
    -h|--help)              # help output
    echo -e "Usage:\n"
    echo -e "   runSiprosPostprocessing.sh [OPTION]...<PARAM>...\n\n"
    echo -e "<PARAMS>\n"
    echo -e "   -in\t SIP file directory (directory containing search results from Sipros Ensemble).\n"
    echo -e "   -c\t configuration file.\n"
    echo -e "   -o\t output folder (If not exist, will create one).\n"
    echo -e "<OPTIONS>\n"
    echo -e "   -h\t help.\n"
    exit 1
    ;;
    -o)			# output directory
    OutputFolder="$2"
    shift # past argument
    ;;
    -in)					# Forward paired end read file -- single file
    SipFolder="$2"
    shift # past argument
    ;;
    -t)
    TabFile="$2"
    shift
    ;;
    -c)			# Output file prefix
    ConfigureFile="$2"
    shift # past argument
    ;;
    *)
    echo "ERROR: Unidentified user variable $key"
    exit 1        				# unknown option
    ;;
esac
shift # past argument or value
done

if [ -z "$SipFolder" ] && [ -z "$ConfigureFile" ]  && [ -z "$OutputFolder" ] ; then
   echo "Input incomplete. Not all required parameters specified. Run with -h for help. Exiting..."
   exit 1
fi

if [ ! -d ${OutputFolder} ]; then
	mkdir -p ${OutputFolder}
fi


# Generate PSM table
if [ "$TabFile" == "" ]; then
TabFile=$(python2 ${exePath}/sipros_psm_tabulating.py -i ${SipFolder}/ -c ${ConfigureFile} -o ${OutputFolder}/)
fi

# PSM Filtering
# echo ${TabFile}
python2 ${exePath}/sipros_ensemble_filtering.py -i ${TabFile} -c ${ConfigureFile} -o ${OutputFolder}/

# Protein Assembly
python2 ${exePath}/sipros_peptides_assembling.py -w ${OutputFolder}/ -c ${ConfigureFile}

# Protein SIP mode clustering
python2 ${exePath}/ClusterSip.py -w ${OutputFolder}/ -c ${ConfigureFile}


#### need to organize
conda activate mono
mono ./ThermoRawFileParser.exe -d=/home/UNT/bz0053/Documents/Projects/MBR/dataset/yeast_ups/raw -L=2 -f=2
conda activate py2
export OMP_NUM_THREADS=24
echo ${OMP_NUM_THREADS}

./Sipros_OpenMP -o /home/UNT/bz0053/Documents/Projects/MBR/dataset/yeast_ups/ -f /home/UNT/bz0053/Documents/Projects/MBR/dataset/yeast_ups/raw/110714_yeast_ups1_2fmol_r1.mzML -c ../configs/SiprosConfig_yeast_Nterm.cfg 
 ./Sipros_OpenMP -o /home/UNT/bz0053/Documents/Projects/MBR/dataset/yeast_ups/result -w /home/UNT/bz0053/Documents/Projects/MBR/dataset/yeast_ups/raw -c ../configs/SiprosConfig_yeast.cfg 

python correctNterm.py
mkdir samples
for sample in result/*Spe2Pep.txt; do path=${sample%..*}; name=${path#result/}; mkdir samples/${name}; done
for sample in result/*Spe2Pep.txt; do path=${sample%..*}; name=${path#result/}; cp ${sample} samples/${name}/${sample#result/}; done
for sample in result_Nterm/correction/*Spe2Pep.txt; do path=${sample%..*}; name=${path#result_Nterm/correction/}; cp ${sample} samples/${name}/${name#result_Nterm/correction/}_Nterm..SE.Spe2Pep.txt; done
conda activate py2

for folder in /home/UNT/bz0053/Documents/Projects/MBR/dataset/yeast_ups/samples_2fmol/*/; do echo ${folder}; done
for folder in /home/UNT/bz0053/Documents/Projects/MBR/dataset/yeast_ups/samples_2fmol/*/; do ./runSiprosFiltering.sh  -in ${folder} -c ../configs/SiprosConfig_yeast.cfg -o ${folder}; mv ${folder}/*.p*.txt ${folder}/peptide/;mv ${folder}/*.tab ${folder}/peptide/; done
for folder in /home/UNT/bz0053/Documents/Projects/MBR/dataset/yeast_ups/samples_2fmol/*/; do ./runSiprosFiltering.sh  -in ${folder} -c ../configs/SiprosConfig_yeast.cfg -o ${folder}; done

# modify for 1% fdr
for folder in /home/UNT/bz0053/Documents/Projects/MBR/dataset/yeast_ups/samples_50fmol/*/; do ./runSiprosFiltering.sh  -in ${folder} -c ../configs/SiprosConfig_yeast.cfg -o ${folder}; done
#modify for single file
./runSiprosFiltering.sh  -in /home/UNT/bz0053/Documents/Projects/MBR/dataset/yeast_ups/samples_50fmol/110618_yeast_ups_50fmol_r1 -c ../configs/SiprosConfig_yeast.cfg -o /home/UNT/bz0053/Documents/Projects/MBR/dataset/yeast_ups/samples_50fmol/110618_yeast_ups_50fmol_r1

#SE2IQ
for folder in /home/UNT/bz0053/Documents/Projects/MBR/dataset/yeast_ups/samples_4fmol/110*; do python SE2prophet.py ${folder}/${folder##*samples_4fmol/}.psm.txt ${folder} ${folder##*samples_4fmol/}; done
for folder in /home/UNT/bz0053/Documents/Projects/MBR/dataset/yeast_ups/samples_4fmol/110**; do python SE2Ionquant.py /home/UNT/bz0053/Documents/Projects/MBR/dataset/yeast_ups/yeast_ups_rev.fasta ${folder}; done

##IQ
conda activate mono
#create filelist use generate.py and copy modmasses file to sample file

java -Xmx23G -Dlibs.bruker.dir=/home/UNT/bz0053/Documents/Projects/fragpipe/tools/MSFragger-4.0/ext/bruker -Dlibs.thermo.dir=/home/UNT/bz0053/Documents/Projects/fragpipe/tools/MSFragger-4.0/ext/thermo -cp /home/UNT/bz0053/Documents/Projects/fragpipe/tools/jfreechart-1.5.3.jar:/home/UNT/bz0053/Documents/Projects/fragpipe/tools/batmass-io-1.30.0.jar:/home/UNT/bz0053/Documents/Projects/fragpipe/tools/IonQuant-1.10.12.jar ionquant.IonQuant --threads 23 --perform-ms1quant 1 --perform-isoquant 0 --isotol 20.0 --isolevel 2 --isotype tmt10 --ionmobility 1 --site-reports 1 --minexps 1 --mbr 1 --maxlfq 1 --requantify 1 --mztol 10 --imtol 0.05 --rttol 0.4 --mbrmincorr 0 --mbrrttol 1 --mbrimtol 0.05 --mbrtoprun 3 --ionfdr 0.01 --proteinfdr 1 --peptidefdr 1 --normalization 1 --minisotopes 2 --minscans 3 --writeindex 0 --tp 0 --minfreq 0 --minions 2 --locprob 0.75 --uniqueness 0 --multidir . --filelist /home/UNT/bz0053/Documents/Projects/MBR/dataset/yeast_ups/samples_50fmol/filelist_ionquant.txt --modlist /home/UNT/bz0053/Documents/Projects/MBR/dataset/yeast_ups/samples_50fmol/modmasses_ionquant.txt
