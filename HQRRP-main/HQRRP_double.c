#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#define max( a, b ) ( (a) > (b) ? (a) : (b) )
#define min( a, b ) ( (a) < (b) ? (a) : (b) )

#define PRINT_DATA
// Turn on CSV_DEBUG to do...

// ============================================================================
// Declaration of local prototypes.

static void matrix_generate( int m_A, int n_A, double * buff_A, int ldim_A, char *filename );

static void print_double_matrix( char * name, int m_A, int n_A, 
                double * buff_A, int ldim_A );

static void print_double_vector( char * name, int n_v, double * buff_v );

static void print_int_vector( char * name, int n_v, int * buff_v );

static void init_pvt( int n_p, int * buff_p );

static void set_pvt_to_zero( int n_p, int * buff_p );

// Read and write binary for doubles
static void read_double_matrix(int m_A, int n_A, double * buff_M, char *filename);

static void write_double_matrix(int m_A, int n_A, double * buff_M, char *filename);

static void write_double_vector(int size, double * buff_V, char *filename);

// Write function for permutation vector
static void write_int_vector(int size, int * buff_V, char *filename);

// ============================================================================
int main( int argc, char *argv[] ) {
  int     nb_alg, pp, m_A, n_A, mn_A, ldim_A, ldim_Q, info, lwork;
  double  * buff_A, * buff_tau, * buff_Q, * buff_wk_qp4, * buff_wk_orgqr;
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
  buff_A   = ( double * ) malloc( m_A * n_A * sizeof( double ) );
  ldim_A   = max( 1, m_A );

  buff_p   = ( int * ) malloc( n_A * sizeof( int ) );

  buff_tau = ( double * ) malloc( n_A * sizeof( double ) );

  buff_Q   = ( double * ) malloc( m_A * mn_A * sizeof( double ) );
  ldim_Q   = max( 1, m_A );

  // reads a binary file of type double which contains a matrix of size m_A*n_A
  read_double_matrix(m_A, n_A, buff_A, fname);
  //write_double_matrix(m_A,n_A,buff_A,"ai.bin");

  #ifdef PRINT_DATA
    // print_double_matrix( "ai", m_A, n_A, buff_A, ldim_A );
    // print_double_vector( "taui", n_A, buff_tau );
  #endif
  
  printf("Post matrix gen\n");
  // Initialize vector with pivots.
  set_pvt_to_zero( n_A, buff_p );
  // changing an index to 1 will fix the row
  // buff_p[ 0 ] = 1;
  // buff_p[ 1 ] = 1;
  // buff_p[ 2 ] = 1;
  // buff_p[ 3 ] = 1;
  // buff_p[ 4 ] = 1;
  printf("Post pivot init\n");

  #ifdef PRINT_DATA
    // print_int_vector( "pi", n_A, buff_p );
  #endif

  // Create workspace.
  lwork       = max( 1, 128 * n_A );
  buff_wk_qp4 = ( double * ) malloc( lwork * sizeof( double ) );

  // Factorize matrix.
  printf( "%% Just before computing factorization.\n" );
  // New factorization.

  // Start Clock
  clock_gettime(CLOCK_REALTIME, &beginWall);
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &beginCPU);

  dgeqp4( & m_A, & n_A, buff_A, & ldim_A, buff_p, buff_tau, 
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
  buff_wk_orgqr = ( double * ) malloc( lwork * sizeof( double ) );
  dlacpy_( "All", & m_A, & mn_A, buff_A, & ldim_A, buff_Q, & ldim_Q );
  dorgqr_( & m_A, & mn_A, & mn_A, buff_Q, & ldim_Q, buff_tau,
           buff_wk_orgqr, & lwork, & info );
  if( info != 0 ) {
    fprintf( stderr, "Error in dorgqr: Info: %d\n", info );
  }
  free( buff_wk_orgqr );

  // // Print results.
  // #ifdef PRINT_DATA
  //   print_double_matrix( "af", m_A, n_A, buff_A, ldim_A );
  //   print_int_vector( "pf", n_A, buff_p );
  //   print_double_vector( "tauf", n_A, buff_tau );
  //   print_double_matrix( "qf", m_A, mn_A, buff_Q, ldim_Q );
  // #endif
  
  // Write results to individual files
  // write_double_matrix(m_A, n_A, buff_A, "af.bin");
  write_int_vector(n_A, buff_p, pname);
  // write_double_vector(n_A, buff_tau, "tauf.bin");
  // write_double_matrix(m_A, n_A, buff_Q, "qf.bin");

  // Free matrices and vectors.
  free( buff_A );
  free( buff_p );
  free( buff_tau );
  free( buff_Q );

  printf( "%% End of Program\n" );

  return 0;
}

static void read_double_matrix(int m_A, int n_A, double * buff_M, char *filename){
  FILE *rBinFile;
  rBinFile = fopen(filename,"r");

  fread(buff_M, sizeof (double), m_A*n_A, rBinFile);

  fclose(rBinFile);
}

static void write_double_matrix(int m_A, int n_A, double * buff_M, char *filename){
  FILE *wBinFile;
  wBinFile = fopen(filename,"w");

  fwrite(buff_M, sizeof (double), m_A*n_A, wBinFile);
  
  fclose(wBinFile);
}

static void write_double_vector(int size, double * buff_V, char *filename){
  FILE *wBinFile;
  wBinFile = fopen(filename,"w");

  fwrite(buff_V, sizeof (double), size, wBinFile);
  
  fclose(wBinFile);
}

static void write_int_vector(int size, int * buff_V, char *filename){
  FILE *wBinFile;
  wBinFile = fopen(filename,"w");

  fwrite(buff_V, sizeof (int), size, wBinFile);
  
  fclose(wBinFile);
}


// ============================================================================
static void print_double_matrix( char * name, int m_A, int n_A, 
                double * buff_A, int ldim_A ) {
  int  i, j;

  printf( "%s = [\n", name );
  for( i = 0; i < m_A; i++ ) {
    for( j = 0; j < n_A; j++ ) {
      printf( "%le ", buff_A[ i + j * ldim_A ] );
    }
    printf( "\n" );
  }
  printf( "];\n" );
}

// ============================================================================
static void print_double_vector( char * name, int n_v, double * buff_v ) {
  int  i;

  printf( "%s = [\n", name );
  for( i = 0; i < n_v; i++ ) {
    printf( "%le\n", buff_v[ i ] );
  }
  printf( "\n" );
  printf( "];\n" );
}

// ============================================================================
static void print_int_vector( char * name, int n_v, int * buff_v ) {
  int  i;

  printf( "%s = [\n", name );
  for( i = 0; i < n_v; i++ ) {
    printf( "%s -- %d/%d: %d\n",name,i,n_v-1,buff_v[i] );
  }
  printf( "];\n" );
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

