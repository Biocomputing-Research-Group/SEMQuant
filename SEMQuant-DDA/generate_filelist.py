import glob
files=glob.glob('Run*')
f=open('filelist_ionquant.txt','w')
f.write('flag\tvalue\n')
for file in files:
	f.write('--psm\t'+file+'/psm.tsv\n')

f.write('--specdir\t/home/UNT/bz0053/Documents/Projects/MBR/dataset/mock_uneven/U1_2rep/raw/raw_file/\n')
#f.write('--specdir\t/home/UNT/bz0053/Documents/Projects/test/raw/\n')
f.close()
