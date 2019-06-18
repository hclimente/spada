# Input variables:
#    - FASTA    Fasta file.
# Output file:
#    - *.tsv

\$(python <<EOF
with open('${FASTA}', 'r') as SEQS:
     for line in SEQS:
         if line[0] == '>':
             with open('{}.fa'.format(line[1:-1]), 'w') as OUT:
                 OUT.write(line)
EOF
)

for f in `ls *.fa`
do
    iupred2a.py \$f long | grep -v '#' | awk '\$3 > 0.5' | sed 's/\\t/ /g' >iupred.out
    tx=`head -n1 \$f | sed 's/>//'`

    touch \$f.tsv
    start=""
    end=""
    seq=""
    while read p; do
    pos=`echo \$p | cut -f1 -d' '`
    res=`echo \$p | cut -f2 -d' '`

    if [ "\$start" == "" ]; then
        start="\$pos"
    elif [ "\$pos" -ne "\$((end + 1))" ]; then
        if [ 5 -lt "\$((end - start))" ]; then
        echo -e "\$tx\tIDR\t\$seq\t\$start\t\$end" >>\$f.tsv
        fi
        seq="\$res"
        start="\$pos"
        end="\$pos"
    else
        end="\$pos"
        seq="\$seq\$res"
    fi
    done <iupred.out

    if [ 5 -lt "\$((end - start))" ]; then
    echo -e "\$tx\tIDR\t\$seq\t\$start\t\$end" >>\$f.tsv
    fi
done

cat *tsv >idr.tsv
