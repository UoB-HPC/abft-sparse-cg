#include "OCLContext.h"

OCLContext::OCLContext()
{
  // TODO: initialize OpenCL context
}

OCLContext::~OCLContext()
{
  // TODO: release OpenCL resources
}

cg_matrix* OCLContext::create_matrix(const uint32_t *columns,
                                     const uint32_t *rows,
                                     const double *values,
                                     int N, int nnz)
{
  // TODO: implement
  return NULL;
}

void OCLContext::destroy_matrix(cg_matrix *mat)
{
  // TODO: implement
}

cg_vector* OCLContext::create_vector(int N)
{
  // TODO: implement
  return NULL;
}

void OCLContext::destroy_vector(cg_vector *vec)
{
  // TODO: implement
}

double* OCLContext::map_vector(cg_vector *v)
{
  // TODO: implement
  return NULL;
}

void OCLContext::unmap_vector(cg_vector *v, double *h)
{
  // TODO: implement
}

void OCLContext::copy_vector(cg_vector *dst, cg_vector *src)
{
  // TODO: implement
}

double OCLContext::dot(cg_vector *a, cg_vector *b)
{
  // TODO: implement
  return 0.0;
}

double OCLContext::calc_xr(cg_vector *x, cg_vector *r,
               cg_vector *p, cg_vector *w,
               double alpha)
{
  // TODO: implement
  return 0.0;
}

void OCLContext::calc_p(cg_vector *p, cg_vector *r, double beta)
{
  // TODO: implement
}

void OCLContext::spmv(cg_matrix *mat, cg_vector *vec, cg_vector *result)
{
  // TODO: implement
}

void OCLContext::inject_bitflip(cg_matrix *mat, BitFlipKind kind, int num_flips)
{
  // TODO: implement
}

namespace
{
  static CGContext::Register<OCLContext> A("ocl", "none");
}
