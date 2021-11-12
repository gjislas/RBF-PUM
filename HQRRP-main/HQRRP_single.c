#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#define max( a, b ) ( (a) > (b) ? (a) : (b) )
#define min( a, b ) ( (a) < (b) ? (a) : (b) )

#define PRINT_DATA

// ============================================================================
// Declaration of local prototypes.

static void init_pvt( int n_p, int * buff_p );

static void set_pvt_to_zero( int n_p, int * buff_p );

// Read and write binary for floating point (single precision)
static void read_float_matrix(int m_A, int n_A, float * buff_M, char *filename);

static void write_float_matrix(int m_A, int n_A, float * buff_M, char *filename);

static void write_float_vector(int size, float * buff_V, char *filename);

// Write function for permutation vector
static void write_int_vector(int size, int * buff_V, char *filename);

// ============================================================================
int main( int argc, char *argv[] ) {
  int     nb_alg, pp, m_A, n_A, mn_A, ldim_A, ldim_Q, info, lwork;
  float  * buff_A, * buff_tau, * buff_Q, * buff_wk_qp4, * buff_wk_orgqr;
  int     * buff_p;

  // Time variable definitions
  struct timespec beginWall, beginCPU, endWall, endCPU; 
  long WallSeconds, WallNanoseconds;
  long CPUSeconds, CPUNanoseconds;
  double CPUTime, WallTime;

  if (argc < 4){
		printf("Please specify the rows, columns, input file, and output file.\n");
		exit(0); 
	}

	m_A  = atoi(argv[1]);
	n_A  = atoi(argv[2]);
	char fname[256], pname[256];	
  strcpy(fname, argv[3]);	
  strcpy(pname, argv[4]); // output file name for buff_p
  
  // Create matrix A, vector p, vector s, and matrix Q.
  mn_A     = min( m_A, n_A );
  buff_A   = ( float * ) malloc( m_A * n_A * sizeof( float ) );
  ldim_A   = max( 1, m_A );

  buff_p   = ( int * ) malloc( n_A * sizeof( int ) );

  buff_tau = ( float * ) malloc( n_A * sizeof( float ) );

  buff_Q   = ( float * ) malloc( m_A * mn_A * sizeof( float ) );
  ldim_Q   = max( 1, m_A );

  // reads a binary file of type double which contains a matrix of size m_A*n_A
  read_float_matrix(m_A, n_A, buff_A, fname);
  //write_float_matrix(m_A,n_A,buff_A,"ai.bin");
  

  // Initialize vector with pivots.
  set_pvt_to_zero( n_A, buff_p );
  // changing an index to 1 will fix the row
  // buff_p[ 0 ] = 1;
  // buff_p[ 1 ] = 1;
  // buff_p[ 2 ] = 1;
  // buff_p[ 3 ] = 1;
  // buff_p[ 4 ] = 1;

  // Create workspace.
  lwork       = max( 1, 128 * n_A );
  buff_wk_qp4 = ( float * ) malloc( lwork * sizeof( float ) );

  // Factorize matrix.
  printf( "%% Just before computing factorization.\n" );
  
  // Start Clock
  clock_gettime(CLOCK_REALTIME, &beginWall);
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &beginCPU);

  sgeqp4( & m_A, & n_A, buff_A, & ldim_A, buff_p, buff_tau, 
          buff_wk_qp4, & lwork, & info );
  
  // Stop Clock
  clock_gettime(CLOCK_REALTIME, &endWall);
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &endCPU);

  // Calculate CPU and Wall elapsed time
  WallSeconds = endWall.tv_sec - beginWall.tv_sec;
  WallNanoseconds = endWall.tv_nsec - beginWall.tv_nsec;

  CPUSeconds = endCPU.tv_sec - beginCPU.tv_sec;
  CPUNanoseconds = endCPU.tv_nsec - beginCPU.tv_nsec;

  printf( "%% Just after computing factorization.\n" );
  printf( "%% Wall time: %d.%09d seconds\n", WallSeconds, WallNanoseconds);
  printf( "%% CPU time: %d.%09d seconds\n", CPUSeconds, CPUNanoseconds);

  // Remove workspace.
  free( buff_wk_qp4 );

  // Build matrix Q.
  lwork     = max( 1, 128 * n_A );
  buff_wk_orgqr = ( float * ) malloc( lwork * sizeof( float ) );
  slacpy_( "All", & m_A, & mn_A, buff_A, & ldim_A, buff_Q, & ldim_Q );
  sorgqr_( & m_A, & mn_A, & mn_A, buff_Q, & ldim_Q, buff_tau,
           buff_wk_orgqr, & lwork, & info );
  if( info != 0 ) {
    fprintf( stderr, "Error in sorgqr: Info: %d\n", info );
  }
  free( buff_wk_orgqr );

  
  //........ Write results to individual files.............

  // write_float_matrix(m_A, n_A, buff_A, "af.bin");
  write_int_vector(n_A, buff_p, pname);
  // write_float_vector(n_A, buff_tau, "tauf.bin");
  // write_float_matrix(m_A, n_A, buff_Q, "qf.bin");

  // Free matrices and vectors.
  free( buff_A );
  free( buff_p );
  free( buff_tau );
  free( buff_Q );

  printf( "%% End of Program\n" );

  return 0;
}

// ============================================================================
static void read_float_matrix(int m_A, int n_A, float * buff_M, char *filename){
  FILE *rBinFile;
  rBinFile = fopen(filename,"r");

  fread(buff_M, sizeof (float), m_A*n_A, rBinFile);

  fclose(rBinFile);
}

// ============================================================================
static void write_float_matrix(int m_A, int n_A, float * buff_M, char *filename){
  FILE *wBinFile;
  wBinFile = fopen(filename,"w");

  fwrite(buff_M, sizeof (float), m_A*n_A, wBinFile);
  
  fclose(wBinFile);
}

// ============================================================================
static void write_float_vector(int size, float * buff_V, char *filename){
  FILE *wBinFile;
  wBinFile = fopen(filename,"w");

  fwrite(buff_V, sizeof (float), size, wBinFile);
  
  fclose(wBinFile);
}

// ============================================================================
static void write_int_vector(int size, int * buff_V, char *filename){
  FILE *wBinFile;
  wBinFile = fopen(filename,"w");

  fwrite(buff_V, sizeof (int), size, wBinFile);
  
  fclose(wBinFile);
}

// ============================================================================
static void init_pvt( int n_p, int * buff_p ) {
  int  i;

  for( i = 0; i < n_p; i++ ) {
    buff_p[ i ] = ( i + 1 );
  }
}

// ============================================================================
static void set_pvt_to_zero( int n_p, int * buff_p ) {
  int  i;

  for( i = 0; i < n_p; i++ ) {
    buff_p[ i ] = 0;
//  printf( "--%d/%d: %d\n",i,n_p,buff_p[ i ] );
  }
}

