import argparse
import subprocess
import os
import csv, sqlite3
import re

con = sqlite3.connect("./sequences/inject.db")
cur = con.cursor()

#type, pos, A, B, BOOL 

cur.execute("DROP TABLE IF EXISTS ref;") # use your column names here
cur.execute("DROP TABLE IF EXISTS aln;") # use your column names here

cur.execute("CREATE TABLE ref (snp INTEGER, seqname VARCHAR(40), pos INTEGER);") # use your column names here
cur.execute("CREATE TABLE aln (seqname VARCHAR(40), start INTEGER, PRIMARY KEY (seqname) );") # use your column names here

parser = argparse.ArgumentParser(description='make a mutation file.')

parser.add_argument('-s', metavar='--snps', type=str, default="polymap.txt",
                        help='polymap')
parser.add_argument('-a', metavar='--aln', type=str, default="mseqs.sort.bam",
                        help='align')
parser.add_argument('-r', metavar='--ref', type=str, default="merge.sort.bam",
                        help='reference')

#def get_pos(pos, cigar):
#	for num1, i_or_d, num2, m in re.findall('(\d+)([ID])(\d+)?([A-Za-z])?', a):	
#		print num1, i_or_d, num2, m
#	return pos+diff
#
#quit()

args = parser.parse_args()

os.system("echo \" PRAGMA synchronous=OFF;\" | sqlite3 inject.db")

proc = subprocess.Popen(['samtools', 'view', '-f', '0x40', '-F', '0x904', args.a], stdout=subprocess.PIPE)

inserts=[]
X=0

ref=""

for line in proc.stdout:
	line=line.split('\t')
	inserts.append([line[0],line[3] ])
	X+=1
	
	if ref=="":
		ref=line[2]

	if X==1000:
		cur.executemany("INSERT INTO aln VALUES (?, ?);", inserts )
		con.commit()
		inserts=[]
		X=0
quit()

snp_file=open(args.s)

inserts=[]
X=0

for line in snp_file:
	snp=line.split(' ')[1]
	X=int(snp)
	pos=ref+":"+str(X)+"-"+str(X+1)
	proc = subprocess.Popen(['samtools', 'view', '-f', '0x40', '-F', '0x904', args.r, pos], stdout=subprocess.PIPE)

	for line in proc.stdout:
		line=line.split('\t')
		inserts.append([str(snp), line[0], str(int(snp)-int(line[3])) ])
		X+=1

	if X>=1000:
		cur.executemany("INSERT INTO ref VALUES (?, ?, ?);", inserts )
		con.commit()
		inserts=[]
		X=0

snp_file.close()
		

if X!=0:
	cur.executemany("INSERT INTO ref VALUES (?, ?, ?);", inserts )
	con.commit()

cur.execute("SELECT ref.snp, aln.start+ref.pos, count(*) FROM ref INNER JOIN aln ON ref.seqname=aln.seqname GROUP BY ref.snp, aln.start+ref.pos")

#rows = cur.fetchall()

this_row=cur.fetchone()
print this_row 

snp_file=open(args.s)

for line in snp_file:
	line=line.strip('\n').split(' ')
	new_snp=[]


	if not (this_row is None):
		while int(line[1])>int(this_row[0]):
			this_row=cur.fetchone()
			if this_row is None:
				break

	if not (this_row is None):
		while int(line[1])==int(this_row[0]):
			new_snp.append(this_row[:])
			this_row=cur.fetchone()
			if this_row is None:
				break
		
	if len(new_snp)>0:
		new_snp.sort(key=lambda x: x[2], reverse=True)
		line[1]=str(new_snp[0][1])
		print ' '.join(line)
	else: 
		line[1]='NaN'
		print ' '.join(line)

con.close()
