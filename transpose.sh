#!/bin/bash

while getopts "a:b:" option
do 
	case $option in
		a) f1=$OPTARG;;
		b) f2=$OPTARG;;
	esac
done

awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' $f1 > $f2
