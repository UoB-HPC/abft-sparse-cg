# abft-sparse-cg

This project implements a simple sparse matrix CG solver for
experimenting with ABFT techniques. In particular, it implements
several software ECC schemes that can detect and (in some cases)
correct single and double bit errors. Both COO and CSR implementations
are provided.

# Building

Running `make` in the top-level directory will build all the
implementation and also download a test matrix to use as input data.

The executables built will be named in the form:

    cg-[coo|csr]-[impl]

Where `impl` is	one of

- `c`
- `arm32`

Running `make test` will perform some quick sanity check on each
executable that is produced.

# Running

    Usage: cg-coo-c [OPTIONS]

    Options:
      -h  --help                  Print this message
      -b  --num-blocks      B     Number of times to block input matrix
      -c  --convergence     C     Convergence threshold
      -f  --matrix-file     M     Path to matrix-market format file
      -i  --iterations      I     Maximum number of iterations
      -m  --mode            MODE  ABFT mode
      -x  --inject-bitflip        Inject a random bit-flip into A

      The -m|--mode argument controls which scheme to use for protecting
      the sparse matrix data. The available options are:
        - NONE (default)
        - CONSTRAINTS
        - SED
        - SEC7
        - SEC8
        - SECDED

      The -x|--inject-bitflip argument optionally takes a number to
      control how many bits to flip, and either INDEX or VALUE to
      restrict the region of bits in the matrix element to target.
