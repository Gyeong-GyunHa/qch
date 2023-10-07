# Quad clustering coefficient in a hypergraph
The quad clustering coefficient, one of the generalisations of the clustering coefficient, is the fraction of quads around a node in a hypergraph.

This `Fortran` code implements the quad clustering coefficient in a hypergraph as presented in this paper.

## Disclaimer
The code has only been tested and compiled with Apple clang version 14.0.0 (GNU Fortran).

## The quad clustering coefficient in a nondirected hypergraph

### Usage
```zsh
> git clone https://github.com/Gyeong-GyunHa/qch
> cd qch
> gfortran nondir_quadclustering.f90
```

If you want to compute it with parallel computation using openmp environment, you can compile it using following command:
```zsh
> gfortran -fopenmp nondir_quadclustering.f90
```

Once the application has been compiled you can run it using the following command:
```zsh
> ./a.out
```

In this code, all options are hardcoded.

### Input
The incidence list of a hypergraph in whitespace-separated values format, with one connection information per line
The meaning of the columns in sample_hypergraph.txt are:
First column: Index of node 
Second column: Index of hyperedge

Note that the index must start from 1 (NOT ZERO)

### Output
The local quad clustering coefficient vector

### Example



## The quad clustering coefficient in a directed hypergraph

### Usage
```zsh
> gfortran -
```

This code also supports parallel computation using openmp environment.
```zsh
> gfortran -
```

Once the application has been compiled you can run it using the following command:
```zsh
> gfortran -
```

### Input
The incidence list of a directed hypergraph

### Output
The local directed quad clustering coefficient vector


## Reference


Gyeong-Gyun Ha, KCL, 2023




