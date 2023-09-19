# Quad clustering coefficient in a hypergraph
The quad clustering coefficient, one of the generalisations of the clustering coefficient, is the fraction of quads around a node in a hypergraph.

This Fortran code implements the quad clustering coefficient in a hypergraph as presented in this paper.

## Disclaimer
The code has only been tested and compiled with Apple clang version 14.0.0 (GNU Fortran).

## The quad clustering coefficient in a nondirected hypergraph

### Usage
```zsh
> git clone https://github.com/Gyeong-GyunHa/qch
> cd qch
> gfortran -
```

If you want to compute it with parallel computation using openmp environment, you can compile it using following command:
```zsh
> gfortran -
```

Once the application has been compiled you can run it using the following command:
```zsh
> gfortran -
```

In this code, all options are hardcoded.

### Input
The incidence list of a hypergraph

### Output
The local quad clustering coefficient vector

### Example



## The quad clustering coefficient in a directed hypergraph

### Usage
```zsh
> gfortran -
```

This code also supports parallel computation using openmp environment.
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




