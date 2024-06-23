#!/usr/bin/env bash
set -eux


P=512
R=1
B=1024
G=0
R=true
A=0
RT=50
for C in 1
do
    echo "++++++++++++++++++ $C Controllers ++++++++++++++++++"
    for T in 1024
    do
        echo "------------------ $T Threads ------------------"
        ../../build/bin/nvm-readwrite-bench --threads=$T --blk_size=$B --reqs=1 --pages=$T --queue_depth=1024 --num_queues=128 --page_size=$P --n_ctrls=$C --gpu=$G --access_type=$A --ratio=$RT --num_blks=$B --random=$R --input ./benchmark1.sh
    done
done