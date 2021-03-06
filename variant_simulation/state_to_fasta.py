import argparse
import os
import csv, sqlite3


parser = argparse.ArgumentParser(description='make a mutation file.')
parser.add_argument('-N', metavar='--number', type=int, default=0,
                        help='column')
parser.add_argument('-s', metavar='--state', type=str, default="states.txt",
                        help='state file')
parser.add_argument('-m', metavar='--mutation', type=str, default="mutation.txt",
                        help='mutation file')
parser.add_argument('-d', metavar='--database', type=str, default="poly.db",
                        help='mutation file')

args = parser.parse_args()

state_file=open(args.s)
mutat_file=open(args.m)

N=args.N

states=[]

con = sqlite3.connect(args.d)
cur = con.cursor()

#type, pos, A, B, BOOL 
cur.execute("DROP TABLE IF EXISTS snps;")
cur.execute("CREATE TABLE snps (sample INTEGER, var VARCHAR(40) );") # use your column names here

os.system("echo \" PRAGMA synchronous=OFF;\" | sqlite3 poly.db")

for line in state_file:
	if line[0]!='@':
		states.append( line[:]  ) 
X=0
inserts=[]
for line in mutat_file:
	state=states.pop(0).strip('\n').split('\t')
	line=line.strip('\n')
	if X==1000:
		cur.executemany("INSERT INTO snps VALUES (?, ?);", inserts )
		con.commit()
		inserts=[]
		X=0
#	print state, len(state), N
	for n in range (0, N):
		if state[n][0]=='1':
			inserts.append([n*2+0, line])
		if state[n][1]=='1':
			inserts.append([n*2+1, line])
		
	
	X=X+1
cur.executemany("INSERT INTO snps VALUES (?, ?);", inserts )
con.commit()
cur.execute("CREATE INDEX all_snps ON snps (sample, var);") # use your column names here
con.close()
