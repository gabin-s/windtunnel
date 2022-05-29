#!/bin/sh

N_RUNS=5
EXEC_PATH='./wind_omp'

max_iter=1000
threshold=-1 # threshold=-1 to prevent early stopping
inlet_pos=0
inlet_size=100
rows=5000
cols=5000

args="$rows $cols $max_iter $threshold $inlet_pos $inlet_size 0 0 0 0 0 0 123 456 789"  # no particles

echo "# running $EXEC_PATH, with <nthreads> thread(s)"
echo "# args: $args"

for nthreads in $(echo 1 2 4 8 $(seq 16 16 256)); do
    printf "%d" $nthreads

    export OMP_NUM_THREADS=$nthreads

    for _ in $(seq $N_RUNS); do
        time_l=$( echo $args | xargs "$EXEC_PATH" | grep "Time: ")

        printf " ${time_l#Time: }"
    done
    echo
done
