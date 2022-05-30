#!/bin/sh

N_RUNS=1
EXEC_PATH='./wind_seq'
n_threads=8

max_iter=1000
threshold=0 # threshold=0 to prevent early stopping
inlet_pos=0
inlet_size=100
rows=1000
cols=1000

# moving random particles
partm_start=1              # first row (CANNOT START AT 0)
partm_size=$(($rows - $partm_start))  # number of rows 
partm_d=0.1                # density

# fixed fixed particles
partf_start=0            # first row (CANNOT START AT 0)
partf_size=0             # number of rows 
partf_d=0                # density

# compute the number of particles, and stores it in $n_particles
estimate_n_particles()
{
    # moving particles
    local d=$(printf "%.0f" ${partm_d}e3)
    local a=$(($partm_size * $cols * $d))
    n_particles=$(printf "%.0f" ${a}e-3)

    # fixed particles
    local d=$(printf "%.0f" ${partf_d}e3)
    local a=$(($partf_size * $cols * $d))
    n_particles=$(($n_particles + $(printf "%.0f" ${a}e-3)))
}

export OMP_NUM_THREADS=$n_threads
args="$rows $rows $max_iter $threshold $inlet_pos $inlet_size"
args="$args $partf_start $partf_size $partf_d $partm_start $partm_size <partm_d>" # particles
args="$args 123 456 789"

echo "# running $EXEC_PATH, with $n_threads thread(s)"
echo "# args: $args"

for partm_d in $(seq 0 0.1 1.0); do
    estimate_n_particles
    printf "%s %d" $partm_d $n_particles

    args="$rows $rows $max_iter $threshold $inlet_pos $inlet_size"
    args="$args $partf_start $partf_size $partf_d $partm_start $partm_size $partm_d" # particles
    args="$args 123 456 789"

    for j in $(seq $N_RUNS); do
        time_l=$( echo $args | xargs "$EXEC_PATH" | grep "Time: ")

        printf " ${time_l#Time: }"
    done
    echo
done
