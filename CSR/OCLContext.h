#include "CGContext.h"

#ifdef __APPLE__
  #include <OpenCL/cl.h>
#else
  #include <CL/cl.h>
#endif

class OCLContext : public CGContext
{
public:
  OCLContext();
  virtual ~OCLContext();

  cg_matrix* create_matrix(const uint32_t *columns,
                           const uint32_t *rows,
                           const double *values,
                           int N, int nnz);
  void destroy_matrix(cg_matrix *mat);

  cg_vector* create_vector(int N);
  void destroy_vector(cg_vector *vec);
  double* map_vector(cg_vector *v);
  void unmap_vector(cg_vector *v, double *h);
  void copy_vector(cg_vector *dst, cg_vector *src);

  double dot(cg_vector *a, cg_vector *b);
  double calc_xr(cg_vector *x, cg_vector *r,
                 cg_vector *p, cg_vector *w,
                 double alpha);
  void calc_p(cg_vector *p, cg_vector *r, double beta);

  void spmv(cg_matrix *mat, cg_vector *vec, cg_vector *result);

  void inject_bitflip(cg_matrix *mat, BitFlipKind kind, int num_flips);
};
