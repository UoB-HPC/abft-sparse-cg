//
// Simple conjugate gradient solver
//

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <sys/time.h>

#include "CGContext.h"

extern "C"
{
  #include "mmio.h"
}

struct
{
  int    num_blocks;
  int    max_itrs;       // max iterations to run
  double conv_threshold; // convergence threshold to stop CG
  const char *matrix_file;

  const char *target;
  const char *mode;

  int    num_bit_flips;  // number of bits to flip in a matrix element
  CGContext::BitFlipKind bitflip_kind;
} params;

double            get_timestamp();
static cg_matrix* load_sparse_matrix(CGContext *context, const char *filename,
                                     int num_blocks, int *N, int *nnz);
void              parse_arguments(int argc, char *argv[]);

int main(int argc, char *argv[])
{
  parse_arguments(argc, argv);

  CGContext *context = CGContext::create(params.target, params.mode);

  int N, nnz;
  cg_matrix *A = load_sparse_matrix(context, params.matrix_file,
                                    params.num_blocks, &N, &nnz);

  cg_vector *b = context->create_vector(N);
  cg_vector *x = context->create_vector(N);
  cg_vector *r = context->create_vector(N);
  cg_vector *p = context->create_vector(N);
  cg_vector *w = context->create_vector(N);

  // Initialize vectors b and x
  double *h_b = context->map_vector(b);
  double *h_x = context->map_vector(x);
  for (unsigned y = 0; y < N; y++)
  {
    h_b[y] = rand() / (double)RAND_MAX;
    h_x[y] = 0.0;
  }
  context->unmap_vector(b, h_b);
  context->unmap_vector(x, h_x);

  printf("\n");
  int block_size = N/params.num_blocks;
  printf("matrix size           = %u x %u\n", N, N);
  printf("matrix block size     = %u x %u\n", block_size, block_size);
  printf("number of non-zeros   = %u (%.4f%%)\n",
         nnz, nnz/((double)N*(double)N)*100);
  printf("maximum iterations    = %u\n", params.max_itrs);
  printf("convergence threshold = %g\n", params.conv_threshold);
  // TODO: Print ABFT mode
  printf("\n");

  // Inject bitflip if required
  if (params.num_bit_flips)
  {
    srand(time(NULL));
    context->inject_bitflip(A, params.bitflip_kind, params.num_bit_flips);
  }

  double start = get_timestamp();

  // r = b - Ax
  // p = r
  context->copy_vector(r, b); // Ax is all zero, if x is all zero
  context->copy_vector(p, r);

  // rr = rT * r
  double rr = context->dot(r, r);

  unsigned itr = 0;
  for (; itr < params.max_itrs && rr > params.conv_threshold; itr++)
  {
    // w = A*p
    context->spmv(A, p, w);

    // pw = pT * A*p
    double pw = context->dot(p, w);

    double alpha = rr / pw;

    // x = x + alpha * p
    // r = r - alpha * A*p
    // rr_new = rT * r
    double rr_new = context->calc_xr(x, r, p, w, alpha);

    double beta = rr_new / rr;

    // p = r + beta * p
    context->calc_p(p, r, beta);

    rr = rr_new;

    if (itr % 1 == 0)
      printf("iteration %5u :  rr = %12.4lf\n", itr, rr);
  }

  double end = get_timestamp();

  printf("\n");
  printf("ran for %u iterations\n", itr);

  printf("\ntime taken = %7.2lf ms\n\n", (end-start)*1e-3);

  // Compute r = Ax
  context->spmv(A, x, r);

  // Compare Ax to b
  double err_sq = 0.0;
  double max_err = 0.0;
  double *h_r = context->map_vector(r);
  h_b = context->map_vector(b);
  for (unsigned i = 0; i < N; i++)
  {
    double err = fabs(h_b[i] - h_r[i]);
    err_sq += err*err;
    max_err = err > max_err ? err : max_err;
  }
  context->unmap_vector(b, h_b);
  context->unmap_vector(r, h_r);
  printf("total error = %lf\n", sqrt(err_sq));
  printf("max error   = %lf\n", max_err);
  printf("\n");

  context->destroy_matrix(A);
  context->destroy_vector(b);
  context->destroy_vector(x);
  context->destroy_vector(r);
  context->destroy_vector(p);
  context->destroy_vector(w);

  delete context;

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
  params.bitflip_kind  = CGContext::ANY;

  params.num_blocks = 25;
  params.matrix_file = "matrices/shallow_water1/shallow_water1.mtx";

  params.target = "cpu";
  params.mode   = "none";

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
    else if (!strcmp(argv[i], "--list") || !strcmp(argv[i], "-l"))
    {
      CGContext::list_contexts();
      exit(0);
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

      params.mode = argv[i];
    }
    else if (!strcmp(argv[i], "--target") || !strcmp(argv[i], "-t"))
    {
      if (++i >= argc)
      {
        printf("Implementation target required\n");
        exit(1);
      }

      params.target = argv[i];
    }
    else if (!strcmp(argv[i], "--inject-bitflip") || !strcmp(argv[i], "-x"))
    {
      params.num_bit_flips = 1;
      while ((i+1) < argc && argv[i+1][0] != '-')
      {
        i++;
        if (!strcmp(argv[i], "INDEX"))
        {
          params.bitflip_kind = CGContext::INDEX;
        }
        else if (!strcmp(argv[i], "VALUE"))
        {
          params.bitflip_kind = CGContext::VALUE;
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

struct matrix_entry
{
  uint32_t col;
  uint32_t row;
  double value;
};

struct matrix_block
{
  matrix_entry *elements;
};

static int compare_matrix_elements(const void *a, const void *b)
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

cg_matrix* load_sparse_matrix(CGContext *context, const char *filename,
                              int num_blocks, int *N, int *nnz)
{
  matrix_block *M = new matrix_block;

  FILE *file = fopen(filename, "r");
  if (file == NULL)
  {
    printf("Failed to open '%s'\n", filename);
    exit(1);
  }

  int width, height, input_nnz;
  mm_read_mtx_crd_size(file, &width, &height, &input_nnz);
  if (width != height)
  {
    printf("Matrix is not square\n");
    exit(1);
  }

  // Load the matrix block from the input file
  int block_nnz = 0;
  M->elements = new matrix_entry[2*input_nnz];
  for (int i = 0; i < input_nnz; i++)
  {
    matrix_entry element;

    int col, row;

    if (fscanf(file, "%d %d %lg\n", &col, &row, &element.value) != 3)
    {
      printf("Failed to read matrix data\n");
      exit(1);
    }
    // adjust from 1-based to 0-based
    col--;
    row--;

    element.col = col;
    element.row = row;
    M->elements[block_nnz] = element;
    block_nnz++;

    if (element.col == element.row)
      continue;

    element.row = col;
    element.col = row;
    M->elements[block_nnz] = element;
    block_nnz++;
  }

  qsort(M->elements, block_nnz, sizeof(matrix_entry), compare_matrix_elements);

  uint32_t *columns = new uint32_t[block_nnz * num_blocks];
  uint32_t *rows    = new uint32_t[block_nnz * num_blocks];
  double   *values  = new double[block_nnz * num_blocks];

  // Duplicate block across diagonal of full matrix
  *nnz = 0;
  for (int j = 0; j < num_blocks; j++)
  {
    for (int i = 0; i < block_nnz; i++)
    {
      matrix_entry element = M->elements[i];

      columns[*nnz] = element.col + j*width;
      rows[*nnz]    = element.row + j*height;
      values[*nnz]  = element.value;

      (*nnz)++;
    }
  }

  *N = width*num_blocks;

  cg_matrix *result = context->create_matrix(columns, rows, values, *N, *nnz);

  delete[] columns;
  delete[] rows;
  delete[] values;

  return result;
}
