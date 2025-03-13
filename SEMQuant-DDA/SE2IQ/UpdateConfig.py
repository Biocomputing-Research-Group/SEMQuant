import sys
filename=sys.argv[1]
fasta_path=sys.argv[2]
contents=[]
with open(filename) as f:
	for line_id, line in enumerate(f):
		if line_id==12:
			s=line.strip().split(' = ')
			s[1]=fasta_path.split('.')[0]+'_rev.fasta'
			news=' = '.join(s)
			contents.append(news+'\n')
		else:
			contents.append(line)
with open(filename,'w') as f:
	for line in contents:
		f.write(line)