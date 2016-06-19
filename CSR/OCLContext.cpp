#include "CGContext.h"

#ifdef __APPLE__
  #include <OpenCL/cl.h>
#else
  #include <CL/cl.h>
#endif

struct cg_vector
{
  int N;
  cl_mem data;
};

struct cg_matrix
{
  unsigned N;
  unsigned nnz;
  cl_mem cols;
  cl_mem rows;
  cl_mem values;
};

class OCLContext : public CGContext
{
public:
  OCLContext()
  {
    // TODO: initialize OpenCL context
  }

  virtual ~OCLContext()
  {
    // TODO: release OpenCL resources
  }

  cg_matrix* create_matrix(const uint32_t *columns,
                           const uint32_t *rows,
                           const double *values,
                           int N, int nnz)
  {
    // TODO: implement
    return NULL;
  }

  void destroy_matrix(cg_matrix *mat)
  {
    // TODO: implement
  }

  cg_vector* create_vector(int N)
  {
    // TODO: implement
    return NULL;
  }

  void destroy_vector(cg_vector *vec)
  {
    // TODO: implement
  }

  double* map_vector(cg_vector *v)
  {
    // TODO: implement
    return NULL;
  }

  void unmap_vector(cg_vector *v, double *h)
  {
    // TODO: implement
  }

  void copy_vector(cg_vector *dst, cg_vector *src)
  {
    // TODO: implement
  }

  double dot(cg_vector *a, cg_vector *b)
  {
    // TODO: implement
    return 0.0;
  }

  double calc_xr(cg_vector *x, cg_vector *r,
                 cg_vector *p, cg_vector *w,
                 double alpha)
  {
    // TODO: implement
    return 0.0;
  }

  void calc_p(cg_vector *p, cg_vector *r, double beta)
  {
    // TODO: implement
  }

  void spmv(cg_matrix *mat, cg_vector *vec, cg_vector *result)
  {
    // TODO: implement
  }

  void inject_bitflip(cg_matrix *mat, BitFlipKind kind, int num_flips)
  {
    // TODO: implement
  }
};

namespace
{
  static CGContext::Register<OCLContext> A("ocl", "none");
}
