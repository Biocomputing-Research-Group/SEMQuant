# debugProxy.py
import os
import sys
import runpy

# 1. cd WORKDIR
os.chdir('test')

# 2A. python test.py 4 5
# args = 'python ../script3/sipros_psm_tabulating.py \
#         -i regular \
#         -c SiprosEnsembleConfig.cfg -o regular'

args = 'python ../script3/sipros_ensemble_filtering.py \
        -i regular \
        -c SiprosEnsembleConfig.cfg \
        -o regular'

# args = 'python ../script3/sipros_psm_tabulating.py \
#         -i sip \
#         -c SiprosEnsembleC13.cfg -o sip'

# args = 'python ../script3/sipros_prepare_protein_database.py \
#         -i fasta/EcoliWithCrapNodup.fasta \
#         -o fasta/DecoySIP.fasta \
#         -c test/SiprosEnsembleC13.cfg'

# args = 'python ../script3/sipros_ensemble_filtering.py \
#         -i sip \
#         -c SiprosEnsembleC13.cfg \
#         -o sip'

# args = 'python ../script3/sipros_peptides_assembling.py \
#         -c SiprosEnsembleC13.cfg \
#         -w sip'

# 2B. python -m mymodule.test 4 5
# args = 'python -m mymodule.test 4 5'

args = args.split()
if args[0] == 'python':
    """pop up the first in the args"""
    args.pop(0)
if args[0] == '-m':
    """pop up the first in the args"""
    args.pop(0)
    fun = runpy.run_module
else:
    fun = runpy.run_path
sys.argv.extend(args[1:])
fun(args[0], run_name='__main__')
