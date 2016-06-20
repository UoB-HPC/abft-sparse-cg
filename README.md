# abft-sparse-cg

This project implements a simple sparse matrix CG solver for
experimenting with ABFT techniques. In particular, it implements
several software ECC schemes that can detect and (in some cases)
correct single and double bit errors. Both COO and CSR implementations
are provided.

# Building

Running `make` in the top-level directory will build all the
implementations and also download a test matrix to use as input data.

The executables built will cg-coo and cg-csr.

Running `make test` will perform some quick sanity check on each
executable that is produced.

# Running

    Usage: cg-csr [OPTIONS]

    Options:
      -h  --help                  Print this message
      -b  --num-blocks      B     Number of times to block input matrix
      -c  --convergence     C     Convergence threshold
      -f  --matrix-file     M     Path to matrix-market format file
      -i  --iterations      I     Maximum number of iterations
      -l  --list                  List available implementations
      -m  --mode            MODE  ABFT mode
      -t  --target          TARG  Implementation target
      -x  --inject-bitflip        Inject a random bit-flip into A

      The -l|--list argument will provide a list of tuples that describe
      which implementations are available to be passed to the
      -t|--target and -m|--mode arguments.

      The -x|--inject-bitflip argument optionally takes a number to
      control how many bits to flip, and either INDEX or VALUE to
      restrict the region of bits in the matrix element to target.
