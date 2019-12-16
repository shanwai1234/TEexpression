import sys

fh = open(sys.argv[1],'r')

for line in fh:
	if line.startswith('@'):
		print (line.strip())
	else:
		new = line.strip().split('\t')
		nlist = []
		for i,j in enumerate(new):
			if 'NH:i:' not in j:
				nlist.append(j)
			else:
				n = j.split(':')
				num = n[-1]
				x = 'RH:i:'+num
				nlist.append('NH:i:1')
				nlist.append(x)
		print ('\t'.join(nlist))

fh.close()
