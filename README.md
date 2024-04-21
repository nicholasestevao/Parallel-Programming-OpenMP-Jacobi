# Parallel-Programming-OpenMP-Jacobi
This project involves developing a parallel application with shared address space, in C/OpenMP, to solve a linear system (Ax=b), according to the Iterative Method of Jacobi-Richardson (also known as Gauss-Jacobi).

## To Compile
### Only Sequential:
```bash
$ make seq
```

### Only Parallel:
```bash
$ make par
```

### Both:
```bash
$ make all
```

## To Run
**Notes:** If you do not use Windows, append .out to the end of the file name.

### Sequential:
``` bash
$ ./jacobiseq <matrix_order> <seed> <option_debug>
```
### Parallel:
``` bash
$ ./jacobiseq <matrix_order> <num_threads> <seed> <option_debug>
```

