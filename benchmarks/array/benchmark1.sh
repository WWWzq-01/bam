#!/usr/bin/env bash
set -eux


P=4096
R=1
B=1024
G=0
for C in 1
do
    echo "++++++++++++++++++ $C Controllers ++++++++++++++++++"
    for T in 1024
    do
        echo "------------------ $T Threads ------------------"
        ../../build/bin/nvm-array-bench --threads=$T --blk_size=$B --reqs=$R --pages=$T --queue_depth=1024 --num_queues=128 --page_size=$P --n_ctrls=$C --gpu=$G | grep "IO"
    done

done