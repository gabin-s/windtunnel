#!/bin/sh

N_RUNS=5
EXEC_PATH='./wind_omp'
n_threads=4



max_iter=1000
threshold=0 # threshold=0 to prevent early stopping
inlet_pos=0
inlet_size=100
rows=1000

export OMP_NUM_THREADS=$n_threads

args_t="$rows <cols> $max_iter $threshold $inlet_pos $inlet_size 0 0 0 0 0 0 123 456 789"  # no particles

echo "# running $EXEC_PATH, with $n_threads thread(s)"
echo "# args: $args_t"

for cols in $(seq 100 500 5100); do
    printf "%d" $cols

    args="$rows $cols $max_iter $threshold $inlet_pos $inlet_size \
        0 0 0 0 0 0 123 456 789"  # no particles

    for _ in $(seq $N_RUNS); do
        time_l=$( echo $args | xargs "$EXEC_PATH" | grep "Time: ")

        printf " ${time_l#Time: }"
    done
    echo
done
