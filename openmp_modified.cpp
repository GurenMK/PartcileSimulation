#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "omp.h"
#include <algorithm>
#include <vector>

//
//  benchmarking program
//
using namespace std; 
int main( int argc, char **argv )
{   
    int navg,nabsavg=0,numthreads; 
    double dmin, absmin=1.0,davg,absavg=0.0;
	
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" ); 
        printf( "-no turns off all correctness checks and particle output\n");   
        return 0;
    }

    int n = read_int( argc, argv, "-n", 1000 );
    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;      

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

    //plane size
    int sx = (int)ceil(sqrt((double)n));
    int sy = ((n+sx-1)/sx);
    //individual bin size
    double bin_x = (sx/sqrt(n))/11;
    double bin_y = (sy/sqrt(n))/11;
    //number of rows/columns in grid
    int s1 = int(sx/bin_x);
    int s2 = int(sy/bin_y);
    //in case plane is not a square
    int s = max(s1,s2);
    //padded vector row/column length
    //empty bins on all edges
    int size = s+2;
    
    //vector of bins
    vector<particle_t> *bins;
    //create bins
    bins = new vector<particle_t>[size*size];

    #pragma omp parallel private(dmin) shared(bins)
    {
    numthreads = omp_get_num_threads();
    for( int step = 0; step < 1000; step++ )
    {
        navg = 0;
        davg = 0.0;
	dmin = 1.0;   
        
        //assign particles to bins
        #pragma omp single
        for( int i = 0; i < n; i++ ) {
            //find i-th particle bin location
            int x = floor(particles[i].x/bin_x);
            int y = floor(particles[i].y/bin_y);
            //account for padding and add particle pointer to the bin
            bins[(y+1)*(s+2)+(x+1)].push_back(particles[i]);
        }
        //
        //  compute all forces
        //
        #pragma omp for reduction (+:navg) reduction(+:davg)
        for( int i = 0; i < n; i++ )
        {
            particles[i].ax = particles[i].ay = 0;
            //particle at i bin location
            int x = floor(particles[i].x/bin_x);
            int y = floor(particles[i].y/bin_y);
            //location of the i-th particle 
            int location_i = (y+1)*(s+2)+(x+1);
            //location of top leftmost neighboring bin of the i particle
            int first_neighbor = location_i-1-size;
            //check neighboring bins and particle binS
            for (int j = 0; j < 9; j++ ) {
                vector <particle_t> :: iterator it;
                //check surrounding bins and own bin
                for(it = bins[first_neighbor+((j/3)*size)+(j%3)].begin(); it != bins[first_neighbor+((j/3)*size)+(j%3)].end(); ++it) {
				    apply_force( particles[i], *it,&dmin,&davg,&navg);
                }
            }    
        }
		
        //
        //  move particles
        //
        #pragma omp for
        for( int i = 0; i < n; i++ ) 
            move( particles[i] );
  
        if( find_option( argc, argv, "-no" ) == -1 ) 
        {
          //
          //  compute statistical data
          //
          #pragma omp master
          if (navg) { 
            absavg += davg/navg;
            nabsavg++;
          }

          #pragma omp critical
	  if (dmin < absmin) absmin = dmin; 
		
          //
          //  save if necessary
          //
          #pragma omp master
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
        #pragma omp single
        for( int i = 0; i < size*size; i++ )
            bins[i].clear();
    }
}
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d,threads = %d, simulation time = %g seconds", n,numthreads, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    // 
    //  -the minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");
    
    //
    // Printing summary data
    //
    if( fsum)
        fprintf(fsum,"%d %d %g\n",n,numthreads,simulation_time);

    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );

    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
