#include "CGContext.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct cg_vector
{
  int N;
  double *data;
};

// 128-bit matrix element
// Bits  0 to  31 are the colum index
// Bits 32 to  63 are the row index
// Bits 64 to 127 are the floating point value
struct matrix_entry
{
  uint32_t col;
  uint32_t row;
  double value;
};

struct cg_matrix
{
  unsigned N;
  unsigned nnz;
  //abft_mode mode;
  matrix_entry *elements;
};

class CPUContext : public CGContext
{
  using CGContext::CGContext;

  cg_matrix* create_matrix(const uint32_t *columns,
                           const uint32_t *rows,
                           const double *values,
                           int N, int nnz)
  {
    cg_matrix *M = new cg_matrix;

    M->N         = N;
    M->nnz       = nnz;
    M->elements  = new matrix_entry[nnz];

    for (int i = 0; i < nnz; i++)
    {
      M->elements[i].col   = columns[i];
      M->elements[i].row   = rows[i];
      M->elements[i].value = values[i];
    }

    // Initialize ECC bits
    // for (int i = 0; i < M->nnz; i++)
    // {
    //   matrix_entry element = M->elements[i];
    //
    //   switch (mode)
    //   {
    //   case NONE:
    //   case CONSTRAINTS:
    //     break;
    //   case SED:
    //     element.col |= ecc_compute_overall_parity(element) << 31;
    //     break;
    //   case SEC7:
    //     element.col |= ecc_compute_col8(element);
    //     break;
    //   case SEC8:
    //   case SECDED:
    //     element.col |= ecc_compute_col8(element);
    //     element.col |= ecc_compute_overall_parity(element) << 24;
    //     break;
    //   }
    //
    //   M->elements[i] = element;
    // }

    return M;
  }

  void destroy_matrix(cg_matrix *mat)
  {
    delete[] mat->elements;
    delete mat;
  }

  cg_vector* create_vector(int N)
  {
    cg_vector *result = new cg_vector;
    result->N    = N;
    result->data = new double[N];
    return result;
  }

  void destroy_vector(cg_vector *vec)
  {
    delete[] vec->data;
    delete vec;
  }

  double* map_vector(cg_vector *v)
  {
    return v->data;
  }

  void unmap_vector(cg_vector *v, double *h)
  {
  }

  void copy_vector(cg_vector *dst, cg_vector *src)
  {
    memcpy(dst->data, src->data, dst->N*sizeof(double));
  }

  double dot(cg_vector *a, cg_vector *b)
  {
    double ret = 0.0;
    for (int i = 0; i < a->N; i++)
    {
      ret += a->data[i] * b->data[i];
    }
    return ret;
  }

  double calc_xr(cg_vector *x, cg_vector *r,
                 cg_vector *p, cg_vector *w,
                 double alpha)
  {
    double ret = 0.0;
    for (int i = 0; i < x->N; i++)
    {
      x->data[i] += alpha * p->data[i];
      r->data[i] -= alpha * w->data[i];

      ret += r->data[i] * r->data[i];
    }
    return ret;
  }

  void calc_p(cg_vector *p, cg_vector *r, double beta)
  {
    for (int i = 0; i < p->N; i++)
    {
      p->data[i] = r->data[i] + beta*p->data[i];
    }
  }

  void spmv(cg_matrix *mat, cg_vector *vec, cg_vector *result)
  {
    // Initialize result vector to zero
    for (unsigned i = 0; i < mat->N; i++)
      result->data[i] = 0.0;

    // Loop over non-zeros in matrix
    for (unsigned i = 0; i < mat->nnz; i++)
    {
      // Load non-zero element
      matrix_entry element = mat->elements[i];

      // Multiply element value by the corresponding vector value
      // and accumulate into result vector
      result->data[element.col] += element.value * vec->data[element.row];
    }
  }
};

namespace
{
  static CGContext::RegisterContext<CPUContext> X("cpu", "none");
}
