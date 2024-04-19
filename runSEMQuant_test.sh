#!/bin/bash
eval "$(conda shell.bash hook)"
rawfile_folder=""
configfile=""
Working_Directory=""
Number_threads=""
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
		echo -e "       -t\t Number of threads\n"
		exit 1
		;;
		-r)
		rawfile_folder="$2"
		shift
		;;
		-c)
		configfile="$2"
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
		*)
		echo "ERROR: Unidentified user variable $key"
		exit 1
		;;
	esac
	shift
done


conda activate mono
mono Raxport.exe -i ${rawfile_folder} -o ${Working_Directory} -j ${Number_threads}
conda deactivate
export OMP_NUM_THREADS=${Number_threads}
./SiprosEnsembleOMP -w ${Working_Directory} -c ${configfile} -o ${Working_Directory}
