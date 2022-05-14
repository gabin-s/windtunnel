/*
 * Simplified simulation of air flow in a wind tunnel
 *
 * MPI version
 *
 * Computacion Paralela, Grado en Informatica (Universidad de Valladolid)
 * 2020/2021
 *
 * v1.4
 *
 * (c) 2021 Arturo Gonzalez Escribano
 */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<limits.h>
#include<sys/time.h>

#include <assert.h>

/* Headers for the MPI assignment versions */
#include<mpi.h>
#include<stddef.h>

#define	PRECISION	10000
#define	STEPS		8

#define TAG_PROPAGATE_RIGHT 123
#define TAG_PROPAGATE_LEFT 456
#define TAG_SYNC_FLOW 789

/*
 * Student: Comment these macro definitions lines to eliminate modules
 *	Module2: Activate effects of particles in the air pressure
 *	Module3: Activate moving particles
 */
#define MODULE2
#define MODULE3

/* Structure to represent a solid particle in the tunnel surface */
typedef struct {
	unsigned char extra;		// Extra field for student's usage
	int pos_row, pos_col;		// Position in the grid
	int mass;			// Particle mass
	int resistance;			// Resistance to air flow
	int speed_row, speed_col;	// Movement direction and speed
	int old_flow;			// To annotate the flow before applying effects
} Particle;


/* 
 * Function to get wall time
 */
double cp_Wtime(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + 1.0e-6 * tv.tv_usec;
}

/* 
 * Macro function to simplify accessing with two coordinates to a flattened array
 * 	This macro-function can be changed and/or optimized by the students
 *
 */
#define accessMat( arr, exp1, exp2 )	arr[ (int)(exp1) * columns + (int)(exp2) ]

/**
 * Macro to simplify access to local matrices to take into consideration the
 * ghost cells (j=0 and j=local_colums+1).
 * For example a column index of 0 will map to the actual column 1.
 */
#define localAccessMat( arr, exp1, exp2 )	arr[ (int)(exp1) * (local_columns + 2) + (int)(exp2) + 1]

/*
 * Function: Update flow in a matrix position
 * 	This function can be changed and/or optimized by the students
 */
int update_flow( int *flow, int *flow_copy, int *particle_locations, int row, int col, int local_columns, int skip_particles, int is_first, int is_last) {
	// Skip update in particle positions
	if ( skip_particles && localAccessMat( particle_locations, row, col ) != 0 ) return 0;

	// Update if border left
	if ( is_first && col == 0 ) {
		localAccessMat( flow, row, col ) = 
			( 
			localAccessMat( flow_copy, row, col ) + 
			localAccessMat( flow_copy, row-1, col ) * 2 + 
			localAccessMat( flow_copy, row-1, col+1 ) 
			) / 4;
	} 
	// Update if border right
	else if ( is_last && col == local_columns-1 ) {
		localAccessMat( flow, row, col ) = 
			( 
			localAccessMat( flow_copy, row, col ) + 
			localAccessMat( flow_copy, row-1, col ) * 2 + 
			localAccessMat( flow_copy, row-1, col-1 ) 
			) / 4;
	}
	// Update in central part
	else {
		localAccessMat( flow, row, col ) = 
			( 
			localAccessMat( flow_copy, row, col ) + 
			localAccessMat( flow_copy, row-1, col ) * 2 + 
			localAccessMat( flow_copy, row-1, col-1 ) +
			localAccessMat( flow_copy, row-1, col+1 ) 
			) / 5;
	}

	// Return flow variation at this position
	return abs( localAccessMat( flow_copy, row, col ) - localAccessMat( flow, row, col ) );
}

/*
 * Function: Move particle
 * 	This function can be changed and/or optimized by the students
 */
void move_particle( int *flow, Particle *particles, int particle, int rows, int local_columns, int columns, int j0 ) {
	// Compute movement for each step
	int step;
	for( step = 0; step < STEPS; step++ ) {
		// Highly simplified phisical model
		int row = particles[ particle ].pos_row / PRECISION;
		int col = particles[ particle ].pos_col / PRECISION;
		int local_col = col - j0;

		int pressure = localAccessMat( flow, row-1, local_col );
		int left, right;

		if ( col == 0 ) left = 0;
		else 			left = pressure - localAccessMat( flow, row-1, local_col-1 );
		if ( col == columns-1 ) right = 0;
		else 					right = pressure - localAccessMat( flow, row-1, local_col+1 );
				
		int flow_row = (int)( (float)pressure / particles[ particle ].mass * PRECISION );
		int flow_col = (int)( (float)(right - left) / particles[ particle ].mass * PRECISION );

		// Speed change	
		particles[ particle ].speed_row =
			( particles[ particle ].speed_row + flow_row ) / 2;
		particles[ particle ].speed_col =
			( particles[ particle ].speed_col + flow_col ) / 2;
	
		// Movement	
		particles[ particle ].pos_row =
			particles[ particle ].pos_row + particles[ particle ].speed_row / STEPS / 2;
		particles[ particle ].pos_col =
			particles[ particle ].pos_col + particles[ particle ].speed_col / STEPS / 2;

		// Control limits
		if ( particles[ particle ].pos_row >= PRECISION * rows )
			particles[ particle ].pos_row = PRECISION * rows - 1;
		if ( particles[ particle ].pos_col < 0 )
			particles[ particle ].pos_col = 0;
		if ( particles[ particle ].pos_col >= PRECISION * columns )
			particles[ particle ].pos_col = PRECISION * columns - 1;
	}
}

#ifdef DEBUG
/* 
 * Function: Print the current state of the simulation 
 */
void print_status_local( int iteration, int rows, int local_columns, int *flow, int num_particles, int *particle_locations, int max_var ) {
	/* 
	 * You don't need to optimize this function, it is only for pretty 
	 * printing and debugging purposes.
	 * It is not compiled in the production versions of the program.
	 * Thus, it is never used when measuring times in the leaderboard
	 */
	int i,j;
	printf("Iteration: %d, max_var: %f\n", 
		iteration, 
		(float)max_var / PRECISION
		);

	printf("  +");
	for( j=0; j<local_columns; j++ ) printf("---");
	printf("+\n");
	for( i=0; i<rows; i++ ) {
		if ( i % STEPS == iteration % STEPS )
			printf("->|");
		else
			printf("  |");

		for( j=0; j<local_columns; j++ ) {
			char symbol;
			if ( localAccessMat( flow, i, j )  >= 10 * PRECISION ) symbol = '*';
			else if ( localAccessMat( flow, i, j ) >= 1 * PRECISION ) symbol = '0' + localAccessMat( flow, i, j ) / PRECISION;
			else if ( localAccessMat( flow, i, j ) >= 0.5 * PRECISION ) symbol = '+';
			else if ( localAccessMat( flow, i, j ) > 0 ) symbol = '.';
			else symbol = ' ';

			if ( localAccessMat( particle_locations, i, j ) > 0 ) printf("[%c]", symbol );
			else printf(" %c ", symbol );
		}
		printf("|\n");
	}
	printf("  +");
	for( j=0; j<local_columns; j++ ) printf("---");
	printf("+\n\n");
}

/* 
 * Function: Print the current state of the simulation 
 */
void print_status( int iteration, int rows, int columns, int *flow, int num_particles, int *particle_locations, int max_var ) {
	/* 
	 * You don't need to optimize this function, it is only for pretty 
	 * printing and debugging purposes.
	 * It is not compiled in the production versions of the program.
	 * Thus, it is never used when measuring times in the leaderboard
	 */
	int i,j;
	printf("Iteration: %d, max_var: %f\n", 
		iteration, 
		(float)max_var / PRECISION
		);

	printf("  +");
	for( j=0; j<columns; j++ ) printf("---");
	printf("+\n");
	for( i=0; i<rows; i++ ) {
		if ( i % STEPS == iteration % STEPS )
			printf("->|");
		else
			printf("  |");

		for( j=0; j<columns; j++ ) {
			char symbol;
			if ( accessMat( flow, i, j )  >= 10 * PRECISION ) symbol = '*';
			else if ( accessMat( flow, i, j ) >= 1 * PRECISION ) symbol = '0' + accessMat( flow, i, j ) / PRECISION;
			else if ( accessMat( flow, i, j ) >= 0.5 * PRECISION ) symbol = '+';
			else if ( accessMat( flow, i, j ) > 0 ) symbol = '.';
			else symbol = ' ';

			if ( accessMat( particle_locations, i, j ) > 0 ) printf("[%c]", symbol );
			else printf(" %c ", symbol );
		}
		printf("|\n");
	}
	printf("  +");
	for( j=0; j<columns; j++ ) printf("---");
	printf("+\n\n");
}
#endif

/*
 * Function: Print usage line in stderr
 */
void show_usage( char *program_name ) {
	fprintf(stderr,"Usage: %s ", program_name );
	fprintf(stderr,"<rows> <columns> <maxIter> <threshold> <inlet_pos> <inlet_size> <fixed_particles_pos> <fixed_particles_size> <fixed_particles_density> <moving_particles_pos> <moving_particles_size> <moving_particles_density> <short_rnd1> <short_rnd2> <short_rnd3> [ <fixed_row> <fixed_col> <fixed_resistance> ... ]\n");
	fprintf(stderr,"\n");
}

/*
 * MAIN PROGRAM
 */
int main(int argc, char *argv[]) {
	int i,j;

	// Simulation data
	int max_iter;			// Maximum number of simulation steps
	int var_threshold;		// Threshold of variability to continue the simulation
	int rows, columns;		// Cultivation area sizes

	int *flow = NULL;		// Wind tunnel air-flow 
	int *flow_copy = NULL;		// Wind tunnel air-flow (ancillary copy)
	int *particle_locations = NULL;	// To quickly locate places with particles

	int inlet_pos;			// First position of the inlet
	int inlet_size;			// Inlet size
	int particles_f_band_pos;	// First position of the band where fixed particles start
	int particles_f_band_size;	// Size of the band where fixed particles start
	int particles_m_band_pos;	// First position of the band where moving particles start
	int particles_m_band_size;	// Size of the band where moving particles start
	float particles_f_density;	// Density of starting fixed particles
	float particles_m_density;	// Density of starting moving particles

	unsigned short random_seq[3];		// Status of the random sequence

	int		num_particles;		// Number of particles
	Particle	*particles;		// List to store cells information


	/* 0. Initialize MPI */
	MPI_Init( &argc, &argv );
	int rank;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );


	/* 1. Read simulation arguments */
	/* 1.1. Check minimum number of arguments */
	if (argc < 16) {
		fprintf(stderr, "-- Error: Not enough arguments when reading configuration from the command line\n\n");
		show_usage( argv[0] );
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}

	/* 1.2. Read simulation area sizes, maximum number of iterations and threshold */
	rows = atoi( argv[1] );
	columns = atoi( argv[2] );
	max_iter = atoi( argv[3] );
	var_threshold = (int)(atof( argv[4] ) * PRECISION);

	/* 1.3. Read inlet data and band of moving particles data */
	inlet_pos = atoi( argv[5] );
	inlet_size = atoi( argv[6] );
	particles_f_band_pos = atoi( argv[7] );
	particles_f_band_size = atoi( argv[8] );
	particles_f_density = atof( argv[9] );
	particles_m_band_pos = atoi( argv[10] );
	particles_m_band_size = atoi( argv[11] );
	particles_m_density = atof( argv[12] );

	/* 1.4. Read random sequences initializer */
	for( i=0; i<3; i++ ) {
		random_seq[i] = (unsigned short)atoi( argv[13+i] );
	}

	/* 1.5. Allocate particles */
	num_particles = 0;
	// Check correct number of parameters for fixed particles
	if (argc > 16 ) {
		if ( (argc - 16) % 3 != 0 ) {
			fprintf(stderr, "-- Error in number of fixed position particles\n\n");
			show_usage( argv[0] );
			MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
		}
		// Get number of fixed particles
		num_particles = (argc - 16) / 3;
	}
	// Add number of fixed and moving particles in the bands
	int num_particles_f_band = (int)( particles_f_band_size * columns * particles_f_density );
	int num_particles_m_band = (int)( particles_m_band_size * columns * particles_m_density );
	num_particles += num_particles_f_band;
	num_particles += num_particles_m_band;

	// Allocate space for particles
	if ( num_particles > 0 ) {
		particles = (Particle *)malloc( num_particles * sizeof( Particle ) );
		if ( particles == NULL ) {
			fprintf(stderr,"-- Error allocating particles structure for size: %d\n", num_particles );
			MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
		}
	}
	else particles = NULL;

	/* 1.6.1. Read fixed particles */
	int particle = 0;
	if (argc > 16 ) {
		int fixed_particles = (argc - 16) / 3;
		for (particle = 0; particle < fixed_particles; particle++) {
			particles[ particle ].pos_row = atoi( argv[ 16 + particle*3 ] ) * PRECISION;
			particles[ particle ].pos_col = atoi( argv[ 17 + particle*3 ] ) * PRECISION;
			particles[ particle ].mass = 0;
			particles[ particle ].resistance = (int)( atof( argv[ 18 + particle*3 ] ) * PRECISION);
			particles[ particle ].speed_row = 0;
			particles[ particle ].speed_col = 0;
		}
	}
	/* 1.6.2. Generate fixed particles in the band */
	for ( ; particle < num_particles-num_particles_m_band; particle++) {
		particles[ particle ].pos_row = (int)( PRECISION * ( particles_f_band_pos + particles_f_band_size * erand48( random_seq ) ) );
		particles[ particle ].pos_col = (int)( PRECISION * columns * erand48( random_seq ) );
		particles[ particle ].mass = 0;
		particles[ particle ].resistance = (int)( PRECISION * erand48( random_seq ) );
		particles[ particle ].speed_row = 0;
		particles[ particle ].speed_col = 0;
	}

	/* 1.7. Generate moving particles in the band */
	for ( ; particle < num_particles; particle++) {
		particles[ particle ].pos_row = (int)( PRECISION * ( particles_m_band_pos + particles_m_band_size * erand48( random_seq ) ) );
		particles[ particle ].pos_col = (int)( PRECISION * columns * erand48( random_seq ) );
		particles[ particle ].mass = (int)( PRECISION * ( 1 + 5 * erand48( random_seq ) ) );
		particles[ particle ].resistance = (int)( PRECISION * erand48( random_seq ) );
		particles[ particle ].speed_row = 0;
		particles[ particle ].speed_col = 0;
	}

#ifdef DEBUG
	// 1.8. Print arguments 
	if ( rank == 0 ) {
		printf("Arguments, Rows: %d, Columns: %d, max_iter: %d, threshold: %f\n", rows, columns, max_iter, (float)var_threshold/PRECISION );
		printf("Arguments, Inlet: %d, %d  Band of fixed particles: %d, %d, %f  Band of moving particles: %d, %d, %f\n", inlet_pos, inlet_size, particles_f_band_pos, particles_f_band_size, particles_f_density, particles_m_band_pos, particles_m_band_size, particles_m_density );
		printf("Arguments, Init Random Sequence: %hu,%hu,%hu\n", random_seq[0], random_seq[1], random_seq[2]);
		printf("Particles: %d\n", num_particles );
		for (int particle=0; particle<num_particles; particle++) {
			printf("Particle[%d] = { %d, %d, %d, %d, %d, %d }\n",
					particle,
					particles[particle].pos_row, 
					particles[particle].pos_col, 
					particles[particle].mass, 
					particles[particle].resistance, 
					particles[particle].speed_row, 
					particles[particle].speed_col
					);
		}
		printf("\n");
	}
#endif // DEBUG


	/* 2. Start global timer */
	MPI_Barrier( MPI_COMM_WORLD );
	double ttotal = cp_Wtime();


/*
 *
 * START HERE: DO NOT CHANGE THE CODE ABOVE THIS POINT
 *
 */

	// Declare variables used in later output
	int max_var = INT_MAX;
	int iter = 0;
	int resultsA[6];
	int resultsB[6];
	int resultsC[6];

	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int prev_rank = rank - 1;
	int next_rank = rank + 1;

	// indicate if it is the last rank (right) 
	int is_last  = (rank == size-1);
	int is_first = (rank == 0);

	// number of columns for each MPI rank
	int columns_per_task = columns / size;
	int local_columns = columns_per_task;
	if(is_last)
		local_columns += columns % size;	
	
	// create the column_t type
	MPI_Datatype column_t;
    MPI_Type_vector(rows, 1, local_columns + 2, MPI_INT, &column_t);
    MPI_Type_commit(&column_t);

	// absolute index of the first column
	int j0 = rank * columns_per_task;

	// array of which rank manage each particle:
	// "who_cares[ particle ]" is the rank which the particle is managed by.
	int* who_cares = (int *)malloc(num_particles * sizeof(int));

	for( particle = 0; particle < num_particles; particle++ ) {
		int r = (particles[ particle ].pos_col / PRECISION) / columns_per_task;
		if(r == size) r = size - 1;
		who_cares[ particle ] = r;

		assert(r >= 0 && r <= size - 1);
	}

	/* 3. Initialization */
	flow               = (int *)calloc(  (size_t)rows * (size_t)(local_columns + 2), sizeof(int) );
	flow_copy          = (int *)malloc(  (size_t)rows * (size_t)(local_columns + 2)* sizeof(int) );
	particle_locations = (int *)calloc(  (size_t)rows * (size_t)(local_columns + 2), sizeof(int) );

	if ( flow == NULL || flow_copy == NULL || particle_locations == NULL ) {
		fprintf(stderr,"-- Error allocating culture structures for size: %d x %d \n", rows, columns );
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}

	/* 4. Simulation */
	for( iter=1; iter<=max_iter && max_var > var_threshold; iter++) {
		// 4.1. Change inlet values each STEP iterations - MPI: IMPLEMENTED, TESTED
		if ( iter % STEPS == 1 ) {
			for ( j=inlet_pos; j<inlet_pos+inlet_size; j++ ) {
				// compute the noise *before* break/continue for predictable results
				double noise = 0.5 - erand48( random_seq );

				// only enter the loop if column j belongs to the current rank
				int j_local = j - j0;
				if(j < j0) continue;
				if(j > j0 + local_columns - 1) continue;

				// 4.1.1. Change the fans phase
				double phase = iter / STEPS * ( M_PI / 4 );
				double phase_step = M_PI / 2 / inlet_size;
				double pressure_level = 9 + 2 * sin( phase + (j-inlet_pos) * phase_step );
	
				// 4.1.2. Add some random noise
				// 4.1.3. Store level in the first row of the ancillary structure
				localAccessMat( flow, 0, j_local ) = (int)(PRECISION * (pressure_level + noise));
			}
		} // End inlet update

	    if(!is_last)
			// send right and receive right
			MPI_Sendrecv(flow + local_columns, 1, column_t, next_rank, TAG_PROPAGATE_RIGHT,
						 flow + local_columns + 1, 1, column_t, next_rank, TAG_PROPAGATE_LEFT,
						 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			
		if(!is_first)
			// send left and receive left
			MPI_Sendrecv(flow + 1, 1, column_t, prev_rank, TAG_PROPAGATE_LEFT,
						 flow, 1, column_t, prev_rank, TAG_PROPAGATE_RIGHT,
						 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

#ifdef MODULE2
#ifdef MODULE3
		// 4.2. Particles movement each STEPS iterations - MPI: IMPLEMENTED, NOT TESTED
		if ( iter % STEPS == 1 ) {
			// Clean particle positions
			memset(particle_locations, 0, sizeof(int) * rows * (local_columns + 2));

			int particle;
			for( particle = 0; particle < num_particles; particle++ ) {
				if(who_cares[ particle ] != rank) continue;

				int mass = particles[ particle ].mass;
				// Fixed particles
				if ( mass == 0 ) continue;
				// Movable particles
				move_particle( flow, particles, particle, rows, local_columns, columns, j0);
			} 

			// synchronize particle locations - IMPLEMENTED, NOT TESTED
			for(particle = 0; particle < num_particles; particle++) {
				Particle *p = &particles[ particle ];

				if(p->mass == 0) continue; // fixed particle

				// update pos_row, pos_col by broadcasting
				MPI_Bcast(&(p->pos_row), 2, MPI_INT, who_cares[ particle ], MPI_COMM_WORLD);

				// update speed_row, speed_col by broadcasting
				MPI_Bcast(&(p->speed_row), 2, MPI_INT, who_cares[ particle ], MPI_COMM_WORLD);

				// update who_cares
				int r = (p->pos_col / PRECISION) / columns_per_task;
				if(r == size) r = size - 1;
				who_cares[ particle ] = r;
				assert(r >= 0 && r < size);
			}

			// Annotate position
			for(particle = 0; particle < num_particles; particle++) {	
				if(who_cares[ particle ] != rank) continue;
				
				localAccessMat( particle_locations, 
					particles[ particle ].pos_row / PRECISION,
					particles[ particle ].pos_col / PRECISION - j0 ) += 1;
			}
		} // End particles movements
#endif // MODULE3

		// 4.3. Effects due to particles each STEPS iterations
		if ( iter % STEPS == 1 ) {
			int particle;
			for( particle = 0; particle < num_particles; particle++ ) {
				if(who_cares[ particle ] != rank) continue;

				int row = particles[ particle ].pos_row / PRECISION;
				int col = particles[ particle ].pos_col / PRECISION - j0;

				update_flow( flow, flow_copy, particle_locations, row, col, local_columns, 0, is_first, is_last);
				particles[ particle ].old_flow = localAccessMat( flow, row, col );
			}
			for( particle = 0; particle < num_particles; particle++ ) {
				if(who_cares[ particle ] != rank) continue;

				int row = particles[ particle ].pos_row / PRECISION;
				int col = particles[ particle ].pos_col / PRECISION - j0;
				int resistance = particles[ particle ].resistance;

				int back = (int)( (long)particles[ particle ].old_flow * resistance / PRECISION ) / localAccessMat( particle_locations, row, col );
				localAccessMat( flow, row, col ) -= back;
			
				localAccessMat( flow, row-1, col ) += back / 2;

				if(is_first && col <= 0)
					localAccessMat( flow, row-1, col ) += back / 4;
				else
					localAccessMat( flow, row-1, col-1 ) += back / 4;
				
				if(is_last && col >= local_columns - 1)
					localAccessMat( flow, row-1, col ) += back / 4;
				else
					localAccessMat( flow, row-1, col+1 ) += back / 4;					
			}

			// synchronize old_flow
			for(particle = 0; particle < num_particles; particle++) {
				MPI_Bcast(&(particles[ particle ].old_flow), 1, MPI_INT, who_cares[ particle ], MPI_COMM_WORLD);
			}
		} // End effects
#endif // MODULE2

		// 4.4. Copy data in the ancillary structure
		memcpy(flow_copy, flow, sizeof(int) * rows * (local_columns + 2));

		// 4.5. Propagation stage
		// 4.5.1. Initialize data to detect maximum variability
		if ( iter % STEPS == 1 ) max_var = 0;

		// 4.5.2. Execute propagation on the wave fronts
		int wave_front = iter % STEPS;
		if ( wave_front == 0 ) wave_front = STEPS;
		int wave;
		for( wave = wave_front; wave < rows; wave += STEPS) {
			if ( wave > iter ) continue;
			int col;
			for( col=0; col<local_columns; col++) {
				int var = update_flow( flow, flow_copy, particle_locations, wave, col, local_columns, 1, is_first, is_last );
				if ( var > max_var ) {
					max_var = var;
				}
			}
		} // End propagation

		// synchronize max_var among all processes
		// TODO: optimize: only stop ranks that converged
		MPI_Allreduce(&max_var, &max_var, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

#ifdef DEBUG
		// 4.7. DEBUG: Print the current state of the simulation at the end of each iteration 
		if (rank == 0)
			print_status_local( iter, rows, local_columns, flow, num_particles, particle_locations, max_var );
#endif

	} // End iterations

	// to store global results
	MPI_Datatype local_mat;
	
	if(rank != 0) {
		// defines a datatype to transfer (rows, local_columns + 2) into (rows, local_columns)
		MPI_Type_vector(rows, local_columns, local_columns + 2, MPI_INT, &local_mat);
		MPI_Type_commit(&local_mat);
		MPI_Send(flow + 1, 1, local_mat, 0, TAG_SYNC_FLOW, MPI_COMM_WORLD);
		MPI_Type_free(&local_mat);
	} 
	else {
		int* flow_global = (int*)malloc(sizeof(int) * rows * columns);

		// copy rank 0's flow to flow_global
		for(int i = 0; i < rows; i++) {
			int stride = (local_columns + 2);
			memcpy(flow_global + i * columns, flow + 1 + i * stride, sizeof(int) * (local_columns - 2));
		}

		// defines a datatype to transfer (rows, local_columns + 2) into (rows, local_columns)
		MPI_Type_vector(rows, columns_per_task, columns, MPI_INT, &local_mat);
		MPI_Type_commit(&local_mat);
		for(int r = 1; r < size - 1; r++) {
			MPI_Recv(flow_global + columns_per_task * r, 1, local_mat, r, 
					 TAG_SYNC_FLOW, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		MPI_Type_free(&local_mat);

		if(size != 1) {
			// receive from the last rank, which may have a different number of columns
			int n_col = columns - columns_per_task * (size - 1);
			int r = size - 1;
			MPI_Type_vector(rows, n_col, columns, MPI_INT, &local_mat);
			MPI_Type_commit(&local_mat);
			MPI_Recv(flow_global + columns_per_task * r, 1, local_mat, r, TAG_SYNC_FLOW, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Type_free(&local_mat);
		}

		// MPI: Fill result arrays used for later output
		int ind;
		for ( ind=0; ind<6; ind++ )
			resultsA[ ind ] = accessMat( flow_global, STEPS-1, ind * columns/6 );

		int res_row = ( iter-1 < rows-1 ) ? iter-1 : rows-1;
		for ( ind=0; ind<6; ind++ )
			resultsB[ ind ] = accessMat( flow_global, res_row/2, ind * columns/6 );

		for ( ind=0; ind<6; ind++ )
			resultsC[ ind ] = accessMat( flow_global, res_row, ind * columns/6 );

#ifdef DEBUG
		print_status( iter, rows, columns, flow_global, num_particles, particle_locations, max_var );
#endif

		free(flow_global);
	}

/*
 *
 * STOP HERE: DO NOT CHANGE THE CODE BELOW THIS POINT
 *
 */

	/* 5. Stop global timer */
	MPI_Barrier( MPI_COMM_WORLD );
	ttotal = cp_Wtime() - ttotal;

	/* 6. Output for leaderboard */
	if ( rank == 0 ) {
		printf("\n");
		/* 6.1. Total computation time */
		printf("Time: %lf\n", ttotal );

		/* 6.2. Results: Statistics */
		printf("Result: %d, %d", 
			iter-1, 
			max_var
			);
		int i;
		for ( i=0; i<6; i++ ) printf(", %d", resultsA[i] );
		for ( i=0; i<6; i++ ) printf(", %d", resultsB[i] );
		for ( i=0; i<6; i++ ) printf(", %d", resultsC[i] );
		printf("\n");
	}

	/* 7. Free resources */	
	free( flow );
	free( flow_copy );
	free( particle_locations );
	free( particles );

	/* 8. End */
	MPI_Finalize();
	return 0;
}
