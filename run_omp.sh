#!/bin/sh

N_RUNS=5

args="2100 457 6300 0.4 1 452 20 2000 0.1 16 50 0.2 583 223 712"

res_ref=$(./wind_seq $args)
l_time=$(echo "$res_ref" | grep 'Time:')
l_expected=$(echo "$res_ref" | grep 'Result:')

echo "Reference result: ${l_expected#Result: }"
echo "Reference time:   ${l_time#Time: }"

for n_threads in 1 2 4 8 16 32; do
    echo -n "$n_threads"

    for j in $(seq $N_RUNS); do
        res=$(OMP_NUM_THREADS=$n_threads ./wind_omp $args)

        l_time=$(echo "$res" | grep 'Time:')
        l_result=$(echo "$res" | grep 'Result:')

        echo -n " ${l_time#Time: }"

        [ "$l_result" != "$l_expected" ] && echo "wrong result" && exit
    done
    echo
done
