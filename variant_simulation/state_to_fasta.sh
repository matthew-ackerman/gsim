ref=$1
states=$2
#SIZE="`head -1 $states | grep '	' -o | wc -l`"
SIZE=`head -2 ../sequences/states.txt | tail -n 1 | cut -d '	' -f 1 | cut -d ':' -f 2`
SIZE=$(($SIZE*2))
POLY=../sequences/polymorphisms.map
POLYDB=../sequences/poly.db
REFSIZE=$((`tail -n +2 $1 | wc -c`-`tail -n +2 $1 | wc -l`))

python mutation_simulation2.py -n $((`wc -l ../sequences/states.txt | cut -d ' ' -f 1` -3)) -l $REFSIZE -o True > ../sequences/temp
python mutation_sort2.py ../sequences/temp $POLYDB > $POLY

rm ../sequences/temp

python state_to_fasta.py -N $((SIZE/2)) -s $states -m $POLY -d $POLYDB

for N in `seq 0 2 $((SIZE-2))` 
do
	name=`printf %03d $((N/2))`
	echo $name

	echo "SELECT var FROM snps WHERE sample=$N;" | sqlite3 $POLYDB  | ./mutate -r $ref -n seq_$name  | gzip - > ../sequences/seq_$name.0.fa.gz
	echo "SELECT var FROM snps WHERE sample=$((N+1));" | sqlite3 $POLYDB  | ./mutate -r $ref -n seq_$name | gzip - > ../sequences/seq_$name.1.fa.gz
done

rm -f $POLYDB
