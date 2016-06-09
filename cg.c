//
// Simple conjugate gradient solver
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#include "mmio.h"

#if COO
#include "COO/common.h"
#elif CSR
#include "CSR/common.h"
#endif

struct
{
  int    num_blocks;
  int    max_itrs;       // max iterations to run
  double conv_threshold; // convergence threshold to stop CG
  const char *matrix_file;

  abft_mode mode;

  int    num_bit_flips;        // number of bits to flip in a matrix element
  int    bitflip_region_start; // start of region in matrix element to bit-flip
  int    bitflip_region_end;   // end of region in matrix element to bit-flip
} params;

double        get_timestamp();
void          parse_arguments(int argc, char *argv[]);


int main(int argc, char *argv[])
{
  parse_arguments(argc, argv);

  sparse_matrix A = load_sparse_matrix(params.matrix_file,
                                       params.num_blocks, params.mode);

  double *b = malloc(A.N*sizeof(double));
  double *x = malloc(A.N*sizeof(double));
  double *r = malloc(A.N*sizeof(double));
  double *p = malloc(A.N*sizeof(double));
  double *w = malloc(A.N*sizeof(double));

  // Initialize vectors b and x
  for (unsigned y = 0; y < A.N; y++)
  {
    b[y] = rand() / (double)RAND_MAX;
    x[y] = 0.0;
  }

  printf("\n");
  int block_size = A.N/params.num_blocks;
  printf("matrix size           = %u x %u\n", A.N, A.N);
  printf("matrix block size     = %u x %u\n", block_size, block_size);
  printf("number of non-zeros   = %u (%.4f%%)\n",
         A.nnz, A.nnz/((double)A.N*(double)A.N)*100);
  printf("maximum iterations    = %u\n", params.max_itrs);
  printf("convergence threshold = %g\n", params.conv_threshold);
  // TODO: Print ABFT mode
  printf("\n");

  if (params.num_bit_flips)
  {
    // Flip a set of random (consecutive) bits in a random matrix element
    srand(time(NULL));
    int index = rand() % A.nnz;
    int region_size = (params.bitflip_region_end - params.bitflip_region_start);
    region_size -= params.num_bit_flips;
    int start_bit   = (rand() % region_size) + params.bitflip_region_start;
    for (int bit = start_bit; bit < start_bit + params.num_bit_flips; bit++)
    {
#if COO
      flip_bit(A.elements+index, bit);
#elif CSR
      csr_colval colval;
      colval.value = A.values[index];
      colval.column = A.cols[index];
      flip_bit(&colval, bit);
      A.values[index] = colval.value;
      A.cols[index] = colval.column;
#endif
      printf("*** flipping bit %d at index %d ***\n", bit, index);
    }
  }

  double start = get_timestamp();

  // r = b - Ax;
  // p = r
  spmv(A, x, r);
  for (unsigned i = 0; i < A.N; i++)
  {
    p[i] = r[i] = b[i] - r[i];
  }

  // rr = rT * r
  double rr = 0.0;
  for (unsigned i = 0; i < A.N; i++)
  {
    rr += r[i] * r[i];
  }

  unsigned itr = 0;
  for (; itr < params.max_itrs && rr > params.conv_threshold; itr++)
  {
    // w = A*p
    spmv(A, p, w);

    // pw = pT * A*p
    double pw = 0.0;
    for (unsigned i = 0; i < A.N; i++)
    {
      pw += p[i] * w[i];
    }

    double alpha = rr / pw;

    // x = x + alpha * p
    // r = r - alpha * A*p
    // rr_new = rT * r
    double rr_new = 0.0;
    for (unsigned i = 0; i < A.N; i++)
    {
      x[i] += alpha * p[i];
      r[i] -= alpha * w[i];

      rr_new += r[i] * r[i];
    }

    double beta = rr_new / rr;

    // p = r + beta * p
    for (unsigned  i = 0; i < A.N; i++)
    {
      p[i] = r[i] + beta*p[i];
    }

    rr = rr_new;

    if (itr % 1 == 0)
      printf("iteration %5u :  rr = %12.4lf\n", itr, rr);
  }

  double end = get_timestamp();

  printf("\n");
  printf("ran for %u iterations\n", itr);

  printf("\ntime taken = %7.2lf ms\n\n", (end-start)*1e-3);

  // Compute Ax
  double *Ax = malloc(A.N*sizeof(double));
  spmv(A, x, Ax);

  // Compare Ax to b
  double err_sq = 0.0;
  double max_err = 0.0;
  for (unsigned i = 0; i < A.N; i++)
  {
    double err = fabs(b[i] - Ax[i]);
    err_sq += err*err;
    max_err = err > max_err ? err : max_err;
  }
  printf("total error = %lf\n", sqrt(err_sq));
  printf("max error   = %lf\n", max_err);
  printf("\n");

  // TODO
  //free(A.elements);
  free(b);
  free(x);
  free(r);
  free(p);
  free(w);
  free(Ax);

  return 0;
}

double get_timestamp()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_usec + tv.tv_sec*1e6;
}

double parse_double(const char *str)
{
  char *next;
  double value = strtod(str, &next);
  return strlen(next) ? -1 : value;
}

int parse_int(const char *str)
{
  char *next;
  int value = strtoul(str, &next, 10);
  return strlen(next) ? -1 : value;
}

void parse_arguments(int argc, char *argv[])
{
  // Set defaults
  params.max_itrs       = 1000;
  params.conv_threshold = 0.001;
  params.num_bit_flips = 0;
  params.bitflip_region_start = 0;
#if COO
  params.bitflip_region_end   = 128;
#elif CSR
  params.bitflip_region_end   = 96;
#endif

  params.num_blocks = 25;
  params.matrix_file = "matrices/shallow_water1/shallow_water1.mtx";
  params.mode = NONE;

  for (int i = 1; i < argc; i++)
  {
    if (!strcmp(argv[i], "--convergence") || !strcmp(argv[i], "-c"))
    {
      if (++i >= argc || (params.conv_threshold = parse_double(argv[i])) < 0)
      {
        printf("Invalid convergence threshold\n");
        exit(1);
      }
    }
    else if (!strcmp(argv[i], "--iterations") || !strcmp(argv[i], "-i"))
    {
      if (++i >= argc || (params.max_itrs = parse_int(argv[i])) < 0)
      {
        printf("Invalid number of iterations\n");
        exit(1);
      }
    }
    else if (!strcmp(argv[i], "--num-blocks") || !strcmp(argv[i], "-b"))
    {
      if (++i >= argc || (params.num_blocks = parse_int(argv[i])) < 1)
      {
        printf("Invalid number of blocks\n");
        exit(1);
      }
    }
    else if (!strcmp(argv[i], "--matrix-file") || !strcmp(argv[i], "-f"))
    {
      if (++i >= argc)
      {
        printf("Matrix filename required\n");
        exit(1);
      }
      params.matrix_file = argv[i];
    }
    else if (!strcmp(argv[i], "--mode") || !strcmp(argv[i], "-m"))
    {
      if (++i >= argc)
      {
        printf("ABFT mode required\n");
        exit(1);
      }

      if (!strcmp(argv[i], "NONE"))
        params.mode = NONE;
      else if (!strcmp(argv[i], "CONSTRAINTS"))
        params.mode = CONSTRAINTS;
      else if (!strcmp(argv[i], "SED"))
        params.mode = SED;
      else if (!strcmp(argv[i], "SEC7"))
        params.mode = SEC7;
      else if (!strcmp(argv[i], "SEC8"))
        params.mode = SEC8;
      else if (!strcmp(argv[i], "SECDED"))
        params.mode = SECDED;
      else
      {
        printf("Invalid ABFT mode\n");
        exit(1);
      }
    }
    else if (!strcmp(argv[i], "--inject-bitflip") || !strcmp(argv[i], "-x"))
    {
      params.num_bit_flips = 1;
      while ((i+1) < argc && argv[i+1][0] != '-')
      {
        i++;
        if (!strcmp(argv[i], "INDEX"))
        {
          params.bitflip_region_start = 0;
#if COO
          params.bitflip_region_end   = 64;
#elif CSR
          params.bitflip_region_end   = 32;
#endif
        }
        else if (!strcmp(argv[i], "VALUE"))
        {
#if COO
          params.bitflip_region_start = 64;
          params.bitflip_region_end   = 128;
#elif CSR
          params.bitflip_region_start = 32;
          params.bitflip_region_end   = 96;
#endif
        }
        else if ((params.num_bit_flips = parse_int(argv[i])) < 1)
        {
          printf("Invalid bit-flip parameter\n");
          exit(1);
        }
      }
    }
    else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h"))
    {
      printf("\n");
      const char *exe = strrchr(argv[0], '/');
      printf("Usage: %s [OPTIONS]\n\n", exe ? exe+1 : argv[0]);
      printf("Options:\n");
      printf(
        "  -h  --help                  Print this message\n"
        "  -b  --num-blocks      B     Number of times to block input matrix\n"
        "  -c  --convergence     C     Convergence threshold\n"
        "  -f  --matrix-file     M     Path to matrix-market format file\n"
        "  -i  --iterations      I     Maximum number of iterations\n"
        "  -m  --mode            MODE  ABFT mode\n"
        "  -x  --inject-bitflip        Inject a random bit-flip into A\n"
        "\n"
        "  The -m|--mode argument controls which scheme to use for protecting\n"
        "  the sparse matrix data. The available options are:\n"
        "    - NONE (default)\n"
        "    - CONSTRAINTS\n"
        "    - SED\n"
        "    - SEC7\n"
        "    - SEC8\n"
        "    - SECDED\n"
        "\n"
        "  The -x|--inject-bitflip argument optionally takes a number to \n"
        "  control how many bits to flip, and either INDEX or VALUE to \n"
        "  restrict the region of bits in the matrix element to target.\n"
      );
      printf("\n");
      exit(0);
    }
    else
    {
      printf("Unrecognized argument '%s' (try '--help')\n", argv[i]);
      exit(1);
    }
  }
}
