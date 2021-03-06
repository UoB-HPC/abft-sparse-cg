#include "CPUContext.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

void CPUContext::generate_ecc_bits(coo_element& element)
{
}

cg_matrix* CPUContext::create_matrix(const uint32_t *columns,
                                     const uint32_t *rows,
                                     const double *values,
                                     int N, int nnz)
{
  cg_matrix *M = new cg_matrix;

  M->N         = N;
  M->nnz       = nnz;
  M->elements  = new coo_element[nnz];

  for (int i = 0; i < nnz; i++)
  {
    coo_element element;
    element.col   = columns[i];
    element.row   = rows[i];
    element.value = values[i];

    generate_ecc_bits(element);

    M->elements[i] = element;
  }

  return M;
}

void CPUContext::destroy_matrix(cg_matrix *mat)
{
  delete[] mat->elements;
  delete mat;
}

cg_vector* CPUContext::create_vector(int N)
{
  cg_vector *result = new cg_vector;
  result->N    = N;
  result->data = new double[N];
  return result;
}

void CPUContext::destroy_vector(cg_vector *vec)
{
  delete[] vec->data;
  delete vec;
}

double* CPUContext::map_vector(cg_vector *v)
{
  return v->data;
}

void CPUContext::unmap_vector(cg_vector *v, double *h)
{
}

void CPUContext::copy_vector(cg_vector *dst, const cg_vector *src)
{
  memcpy(dst->data, src->data, dst->N*sizeof(double));
}

double CPUContext::dot(const cg_vector *a, const cg_vector *b)
{
  double ret = 0.0;
  for (int i = 0; i < a->N; i++)
  {
    ret += a->data[i] * b->data[i];
  }
  return ret;
}

double CPUContext::calc_xr(cg_vector *x, cg_vector *r,
                           const cg_vector *p, const cg_vector *w,
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

void CPUContext::calc_p(cg_vector *p, const cg_vector *r, double beta)
{
  for (int i = 0; i < p->N; i++)
  {
    p->data[i] = r->data[i] + beta*p->data[i];
  }
}

void CPUContext::spmv(const cg_matrix *mat, const cg_vector *vec,
                      cg_vector *result)
{
  // Initialize result vector to zero
  for (unsigned i = 0; i < mat->N; i++)
    result->data[i] = 0.0;

  // Loop over non-zeros in matrix
  for (unsigned i = 0; i < mat->nnz; i++)
  {
    // Load non-zero element
    coo_element element = mat->elements[i];

    // Multiply element value by the corresponding vector value
    // and accumulate into result vector
    result->data[element.col] += element.value * vec->data[element.row];
  }
}

void CPUContext::inject_bitflip(cg_matrix *mat, BitFlipKind kind, int num_flips)
{
  int index = rand() % mat->nnz;

  int start = 0;
  int end   = 128;
  if (kind == VALUE)
    start = 64;
  else if (kind == INDEX)
    end = 64;

  for (int i = 0; i < num_flips; i++)
  {
    int bit   = (rand() % (end-start)) + start;
    printf("*** flipping bit %d at index %d ***\n", bit, index);
    ((uint32_t*)(mat->elements+index))[bit/32] ^= 0x1U << (bit % 32);
  }
}

void CPUContext_Constraints::spmv(const cg_matrix *mat, const cg_vector *vec,
                                  cg_vector *result)
{
  // Initialize result vector to zero
  for (unsigned i = 0; i < mat->N; i++)
    result->data[i] = 0.0;

  // Loop over non-zeros in matrix
  for (unsigned i = 0; i < mat->nnz; i++)
  {
    // Load non-zero element
    coo_element element = mat->elements[i];

    // Check index size constraints
    if (element.row >= mat->N)
    {
      printf("row size constraint violated for index %d\n", i);
      exit(1);
    }
    if (element.col >= mat->N)
    {
      printf("column size constraint violated for index %d\n", i);
      exit(1);
    }

    // Check index order constraints
    // Skip last row
    if (i < mat->nnz - 1)
    {
      // Compare this row index to the next row index
      uint32_t next_row = mat->elements[i+1].row;
      if (element.row > next_row)
      {
        printf("row index order violated at index %d\n", i);
        exit(1);
      }
      else if (element.row == next_row)
      {
        // Compare this column index to the next column index (if same row)
        uint32_t next_col = mat->elements[i+1].col;
        if (element.col >= next_col)
        {
          printf("column index order violated at index %d\n", i);
          exit(1);
        }
      }
    }

    // Multiply element value by the corresponding vector value
    // and accumulate into result vector
    result->data[element.col] += element.value * vec->data[element.row];
  }
}

void CPUContext_SED::generate_ecc_bits(coo_element& element)
{
  element.col |= ecc_compute_overall_parity(element) << 31;
}

void CPUContext_SED::spmv(const cg_matrix *mat, const cg_vector *vec,
                          cg_vector *result)
{
  // Initialize result vector to zero
  for (unsigned i = 0; i < mat->N; i++)
    result->data[i] = 0.0;

  // Loop over non-zeros in matrix
  for (unsigned i = 0; i < mat->nnz; i++)
  {
    // Load non-zero element
    coo_element element = mat->elements[i];

    // Check overall parity bit
    if (ecc_compute_overall_parity(element))
    {
      printf("[ECC] error detected at index %d\n", i);
      exit(1);
    }

    // Mask out ECC from high order column bits
    element.col &= 0x00FFFFFF;

    // Multiply element value by the corresponding vector value
    // and accumulate into result vector
    result->data[element.col] += element.value * vec->data[element.row];
  }
}

void CPUContext_SEC7::generate_ecc_bits(coo_element& element)
{
  element.col |= ecc_compute_col8(element);
}

void CPUContext_SEC7::spmv(const cg_matrix *mat, const cg_vector *vec,
                           cg_vector *result)
{
  // Initialize result vector to zero
  for (unsigned i = 0; i < mat->N; i++)
    result->data[i] = 0.0;

  // Loop over non-zeros in matrix
  for (unsigned i = 0; i < mat->nnz; i++)
  {
    // Load non-zero element
    coo_element element = mat->elements[i];

    // Check ECC
    uint32_t syndrome = ecc_compute_col8(element);
    if (syndrome)
    {
      // Unflip bit
      uint32_t bit = ecc_get_flipped_bit_col8(syndrome);
      ((uint32_t*)(&element))[bit/32] ^= 0x1U << (bit % 32);
      mat->elements[i] = element;

      printf("[ECC] corrected bit %u at index %d\n", bit, i);
    }

    // Mask out ECC from high order column bits
    element.col &= 0x00FFFFFF;

    // Multiply element value by the corresponding vector value
    // and accumulate into result vector
    result->data[element.col] += element.value * vec->data[element.row];
  }
}

void CPUContext_SEC8::generate_ecc_bits(coo_element& element)
{
  element.col |= ecc_compute_col8(element);
  element.col |= ecc_compute_overall_parity(element) << 24;
}

void CPUContext_SEC8::spmv(const cg_matrix *mat, const cg_vector *vec,
                           cg_vector *result)
{
  // Initialize result vector to zero
  for (unsigned i = 0; i < mat->N; i++)
    result->data[i] = 0.0;

  // Loop over non-zeros in matrix
  for (unsigned i = 0; i < mat->nnz; i++)
  {
    // Load non-zero element
    coo_element element = mat->elements[i];

    // Check overall parity bit
    if (ecc_compute_overall_parity(element))
    {
      // Compute error syndrome from hamming bits
      uint32_t syndrome = ecc_compute_col8(element);
      if (syndrome)
      {
        // Unflip bit
        uint32_t bit = ecc_get_flipped_bit_col8(syndrome);
        ((uint32_t*)(&element))[bit/32] ^= 0x1U << (bit % 32);

        printf("[ECC] corrected bit %u at index %d\n", bit, i);
      }
      else
      {
        // Correct overall parity bit
        element.col ^= 0x1U << 24;

        printf("[ECC] corrected overall parity bit at index %d\n", i);
      }
      mat->elements[i] = element;
    }

    // Mask out ECC from high order column bits
    element.col &= 0x00FFFFFF;

    // Multiply element value by the corresponding vector value
    // and accumulate into result vector
    result->data[element.col] += element.value * vec->data[element.row];
  }
}

void CPUContext_SECDED::generate_ecc_bits(coo_element& element)
{
  element.col |= ecc_compute_col8(element);
  element.col |= ecc_compute_overall_parity(element) << 24;
}

void CPUContext_SECDED::spmv(const cg_matrix *mat, const cg_vector *vec,
                             cg_vector *result)
{
  // Initialize result vector to zero
  for (unsigned i = 0; i < mat->N; i++)
    result->data[i] = 0.0;

  // Loop over non-zeros in matrix
  for (unsigned i = 0; i < mat->nnz; i++)
  {
    // Load non-zero element
    coo_element element = mat->elements[i];

    // Check parity bits
    uint32_t overall_parity = ecc_compute_overall_parity(element);
    uint32_t syndrome = ecc_compute_col8(element);
    if (overall_parity)
    {
      if (syndrome)
      {
        // Unflip bit
        uint32_t bit = ecc_get_flipped_bit_col8(syndrome);
        ((uint32_t*)(&element))[bit/32] ^= 0x1U << (bit % 32);

        printf("[ECC] corrected bit %u at index %d\n", bit, i);
      }
      else
      {
        // Correct overall parity bit
        element.col ^= 0x1U << 24;

        printf("[ECC] corrected overall parity bit at index %d\n", i);
      }
      mat->elements[i] = element;
    }
    else
    {
      if (syndrome)
      {
        // Overall parity fine but error in syndrom
        // Must be double-bit error - cannot correct this
        printf("[ECC] double-bit error detected\n");
        exit(1);
      }
    }

    // Mask out ECC from high order column bits
    element.col &= 0x00FFFFFF;

    // Multiply element value by the corresponding vector value
    // and accumulate into result vector
    result->data[element.col] += element.value * vec->data[element.row];
  }
}

namespace
{
  static CGContext::Register<CPUContext> A("cpu", "none");
  static CGContext::Register<CPUContext_Constraints> B("cpu", "constraints");
  static CGContext::Register<CPUContext_SED> C("cpu", "sed");
  static CGContext::Register<CPUContext_SEC7> D("cpu", "sec7");
  static CGContext::Register<CPUContext_SEC8> E("cpu", "sec8");
  static CGContext::Register<CPUContext_SECDED> F("cpu", "secded");
}
