#!/usr/bin/env bash
set -eux

if [ $# -ne 2 ]
then
	echo Usage $0 numssd gpuid && exit 1
fi

#when varying cacheline size, equally vary the num pages per cache (-p) to retain cache capacity.
P=4096
R=1 #1 random 0 sequential
G=2
B=64
D=1024
Q=135
PTR=0
G=$2

CTRLS=$1
for ((C=1; C<=$CTRLS; C++))
do
    echo "++++++++++++++++++ $C Controllers ++++++++++++++++++"
    for T in 1024
    do
        echo "------------------ $T Threads ------------------"
        ../../build/bin/nvm-cache-bench -k $C -p $((1024*1024*2)) -P $P -n 1 -t $T -b $B -d $D -q $Q -T $PTR -e 8589934592 --gpu=$G -r $R
    done
done
