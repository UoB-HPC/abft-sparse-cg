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

#include "common.h"

struct
{
  int    num_blocks;
  int    max_itrs;       // max iterations to run
  double conv_threshold; // convergence threshold to stop CG
  const char *matrix_file;

  int    num_bit_flips;        // number of bits to flip in a matrix element
  int    bitflip_region_start; // start of region in matrix element to bit-flip
  int    bitflip_region_end;   // end of region in matrix element to bit-flip
} params;

double        get_timestamp();
sparse_matrix load_sparse_matrix(int *N);
void          parse_arguments(int argc, char *argv[]);


int main(int argc, char *argv[])
{
  parse_arguments(argc, argv);

  int N;
  sparse_matrix A = load_sparse_matrix(&N);

  init_matrix_ecc(A);

  double *b = malloc(N*sizeof(double));
  double *x = malloc(N*sizeof(double));
  double *r = malloc(N*sizeof(double));
  double *p = malloc(N*sizeof(double));
  double *w = malloc(N*sizeof(double));

  // Initialize vectors b and x
  for (unsigned y = 0; y < N; y++)
  {
    b[y] = rand() / (double)RAND_MAX;
    x[y] = 0.0;
  }

  printf("\n");
  int block_size = N/params.num_blocks;
  printf("matrix size           = %u x %u\n", N, N);
  printf("matrix block size     = %u x %u\n", block_size, block_size);
  printf("number of non-zeros   = %u (%.4f%%)\n",
         A.nnz, A.nnz/((double)N*(double)N)*100);
  printf("maximum iterations    = %u\n", params.max_itrs);
  printf("convergence threshold = %g\n", params.conv_threshold);
  printf("\n");

  if (params.num_bit_flips)
  {
    // Flip a set of random (consecutive) bits in a random matrix element
    srand(time(NULL));
    int index = rand() % A.nnz;
    int col   = A.elements[index].col & 0x00FFFFFF;
    int row   = A.elements[index].row & 0x00FFFFFF;
    int region_size = (params.bitflip_region_end - params.bitflip_region_start);
    region_size -= params.num_bit_flips;
    int start_bit   = (rand() % region_size) + params.bitflip_region_start;
    for (int bit = start_bit; bit < start_bit + params.num_bit_flips; bit++)
    {
      int word  = bit / 32;
      ((uint32_t*)(A.elements+index))[word] ^= 1<<(bit%32);
      printf("*** flipping bit %d of element (%d,%d) ***\n", bit, col, row);
    }
  }

  double start = get_timestamp();

  // r = b - Ax;
  // p = r
  spmv(A, x, r, N);
  for (unsigned i = 0; i < N; i++)
  {
    p[i] = r[i] = b[i] - r[i];
  }

  // rr = rT * r
  double rr = 0.0;
  for (unsigned i = 0; i < N; i++)
  {
    rr += r[i] * r[i];
  }

  unsigned itr = 0;
  for (; itr < params.max_itrs && rr > params.conv_threshold; itr++)
  {
    // w = A*p
    spmv(A, p, w, N);

    // pw = pT * A*p
    double pw = 0.0;
    for (unsigned i = 0; i < N; i++)
    {
      pw += p[i] * w[i];
    }

    double alpha = rr / pw;

    // x = x + alpha * p
    // r = r - alpha * A*p
    // rr_new = rT * r
    double rr_new = 0.0;
    for (unsigned i = 0; i < N; i++)
    {
      x[i] += alpha * p[i];
      r[i] -= alpha * w[i];

      rr_new += r[i] * r[i];
    }

    double beta = rr_new / rr;

    // p = r + beta * p
    for (unsigned  i = 0; i < N; i++)
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
  double *Ax = malloc(N*sizeof(double));
  spmv(A, x, Ax, N);

  // Compare Ax to b
  double err_sq = 0.0;
  double max_err = 0.0;
  for (unsigned i = 0; i < N; i++)
  {
    double err = fabs(b[i] - Ax[i]);
    err_sq += err*err;
    max_err = err > max_err ? err : max_err;
  }
  printf("total error = %lf\n", sqrt(err_sq));
  printf("max error   = %lf\n", max_err);
  printf("\n");

  free(A.elements);
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
  params.bitflip_region_end   = 128;

  params.num_blocks = 25;
  params.matrix_file = "matrices/shallow_water1/shallow_water1.mtx";

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
    else if (!strcmp(argv[i], "--matrix-file") || !strcmp(argv[i], "-m"))
    {
      if (++i >= argc)
      {
        printf("Matrix filename required\n");
        exit(1);
      }
      params.matrix_file = argv[i];
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
          params.bitflip_region_end   = 64;
        }
        else if (!strcmp(argv[i], "VALUE"))
        {
          params.bitflip_region_start = 64;
          params.bitflip_region_end   = 128;
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
      printf("Usage: ./cg [OPTIONS]\n\n");
      printf("Options:\n");
      printf(
        "  -h  --help                 Print this message\n"
        "  -b  --num-blocks      B    Number of times to block input matrix\n"
        "  -c  --convergence     C    Convergence threshold\n"
        "  -i  --iterations      I    Maximum number of iterations\n"
        "  -m  --matrix-file     M    Path to matrix-market format file\n"
        "  -p  --percent-nzero   P    Percentage of A to be non-zero (approx)\n"
        "  -x  --inject-bitflip       Inject a random bit-flip into A\n"
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

int compare_matrix_elements(const void *a, const void *b)
{
  matrix_entry _a = *(matrix_entry*)a;
  matrix_entry _b = *(matrix_entry*)b;

  if (_a.row < _b.row)
  {
    return -1;
  }
  else if (_a.row > _b.row)
  {
    return 1;
  }
  else
  {
    return _a.col - _b.col;
  }
}

// Load a sparse matrix from a matrix-market format file
sparse_matrix load_sparse_matrix(int *N)
{
  sparse_matrix M = {0, NULL};

  FILE *file = fopen(params.matrix_file, "r");
  if (file == NULL)
  {
    printf("Failed to open '%s'\n", params.matrix_file);
    exit(1);
  }

  int width, height, nnz;
  mm_read_mtx_crd_size(file, &width, &height, &nnz);
  if (width != height)
  {
    printf("Matrix is not square\n");
    exit(1);
  }

  M.nnz = 0;
  M.elements = malloc(params.num_blocks*2*nnz*sizeof(matrix_entry));
  for (int i = 0; i < nnz; i++)
  {
    matrix_entry element;

    int col, row;

    if (fscanf(file, "%d %d %lg\n", &col, &row, &element.value) != 3)
    {
      printf("Failed to read matrix data\n");
      exit(1);
    }
    col--; /* adjust from 1-based to 0-based */
    row--;

    element.col = col;
    element.row = row;
    M.elements[M.nnz] = element;
    M.nnz++;

    if (element.col == element.row)
      continue;

    element.row = col;
    element.col = row;
    M.elements[M.nnz] = element;
    M.nnz++;
  }

  qsort(M.elements, M.nnz, sizeof(matrix_entry), compare_matrix_elements);

  nnz = M.nnz;
  for (int j = 1; j < params.num_blocks; j++)
  {
    for (int i = 0; i < nnz; i++)
    {
      matrix_entry element = M.elements[i];
      element.col = element.col + j*width;
      element.row = element.row + j*height;
      M.elements[M.nnz] = element;
      M.nnz++;
    }
  }

  *N = width*params.num_blocks;

  return M;
}
