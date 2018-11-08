import sys
import argparse
import os
import csv, sqlite3

con = sqlite3.connect(sys.argv[2])
cur = con.cursor()
cur.execute("CREATE TABLE mut (snpid INTEGER, type VARCHAR(10), locus INTEGER, arg1 INTEGER, arg2 INTEGER, arg3 VARCHAR(6), arg4 INTEGER, arg5 INTEGER, PRIMARY KEY (snpid) );") # use your column names here
os.system("echo \" PRAGMA synchronous=OFF;\" | sqlite3 "+sys.argv[2])
File=open(sys.argv[1])
X=0
SNPID=0
inserts=[]
for line in File:
	line=[str(SNPID)]+line.strip('\n').split(' ')
	SNPID+=1
	inserts.append(line)
	X+=1
	if X==1000:
		cur.executemany("INSERT INTO mut VALUES (?, ?, ?, ?, ?, ?, ?, ?);", inserts )
		con.commit()
		inserts=[]
		X=0
		print "committed ", SNPID


cur.executemany("INSERT INTO mut VALUES (?, ?, ?, ?, ?, ?, ?, ?);", inserts )
con.commit()

rows=cur.execute("SELECT * FROM mut WHERE locus IN (SELECT locus FROM mut GROUP BY locus HAVING COUNT (locus) > 1) ORDER BY locus;").fetchall()

group=""

print len(rows)

for row in rows:
	snpid=row[0]
	locus=int(row[2])
	if(locus!=group):
		group=locus
		L=group
		continue
	while (len(cur.execute("SELECT locus FROM mut WHERE locus = "+str(L) ).fetchall() ) != 0) :
		L+=1
	cur.execute("UPDATE mut SET locus = "+str(L)+" WHERE snpid = "+str(snpid) )
	print ("UPDATE mut SET locus = "+str(L)+" WHERE snpid = "+str(snpid) )
		
for row in cur.execute("SELECT * FROM mut ORDER BY locus;"):
	#line=row[0]
        print row[1], row[2], row[3], row[4], row[5], row[6], row[7]

con.close()

#File=open(sys.argv[1])
#out={}
#for line in File:
#	line=line.strip('\n').split(' ')
#	while(True):
#		try:
#			get=out[line[1]]
#			line[1]=str(int(line[1])+1)
#		except:
#			out[line[1]]=line
#			break
#keys=map(int, out.keys() )
#keys.sort()
#for k in keys:
#	print ' '.join(out[str(k)])

