import glob
import sys

file_path=sys.argv[1]

files=glob.glob(file_path+"/20*")
f=open(file_path + '/filelist_ionquant.txt','w')
f.write('flag\tvalue\n')
for file in files:
	f.write('--psm\t'+file+'/psm.tsv\n')

#f.write('--specdir\t/home/UNT/bz0053/Documents/Projects/test/raw/\n')
f.write('--specdir\t'+sys.argv[2]+'\n')
f.close()


