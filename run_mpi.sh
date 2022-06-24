#!/bin/sh
#                                 thresh   inlet_size
#                                     | inlet_pos | fixed lines
#                       row col   T   |   |       | [fix] [mob] [short rnd]
mpiexec -n 4 ./wind_mpi 32 210000 118 0.1 0 2100000 0 0 0 0 0 0 673 3902 43