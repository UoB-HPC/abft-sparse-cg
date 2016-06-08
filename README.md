# abft-sparse-cg

This project implements a simple sparse matrix CG solver for
experimenting with ABFT techniques. In particular, it implements
several software ECC schemes that can detect and (in some cases)
correct single and double bit errors. Both COO and CSR implementations
are provided.

# Building

Running `make` in the top-level directory will build all the
implementation variants and also download a test matrix to use as
input data.

The executables built will be named in the form:

    cg-[coo|csr]-[impl]-[scheme]

Where `impl` is	one of	

- `c`
- `arm32`

and `scheme` is	one of

- `baseline`
- `constraints`
- `sed`
- `sec7`
- `sec8`
- `secded`

Running	`make test` will perform a quick sanity	check on each
executable that is produced.

# Running

    Usage: cg-coo-c-baseline [OPTIONS]

    Options:
      -h  --help                 Print this message
      -b  --num-blocks      B    Number of times to block input matrix
      -c  --convergence     C    Convergence threshold
      -i  --iterations      I    Maximum number of iterations
      -m  --matrix-file     M    Path to matrix-market format file
      -p  --percent-nzero   P    Percentage of A to be non-zero (approx)
      -x  --inject-bitflip       Inject a random bit-flip into A

      The -x|--inject-bitflip argument optionally takes a number to
      control how many bits to flip, and either INDEX or VALUE to
      restrict the region of bits in the matrix element to target.
