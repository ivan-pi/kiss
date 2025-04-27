# kiss

krylov-based iterative sparse solvers

## Solvers

Systems of linear equations:
* [x] bicg
* [x] bicgstab
* [x] cg
* [x] cgs
* [ ] gmres
* [ ] lgmres
* [ ] minres
* [ ] gcrotmk
* [x] qmr
* [x] tfqmr

## Constraints

- vectors are contiguous Fortran arrays
- routines must be interchangeable
- feature parity with `scipy.sparse.linalg`

## Future ideas

- OpenMP for GPU parallelism
- verbose mode with formatted output