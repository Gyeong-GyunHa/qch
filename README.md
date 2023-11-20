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

If you want to compute it with parallel computation using openmp environment, you can compile it using the following command:
```zsh
> gfortran -fopenmp nondir_quadclustering.f90
```

Once the application has been compiled you can run it using the following command:
```zsh
> ./a.out
```

In this code, all options are hardcoded.

### Input
The incidence list of a hypergraph with one connection information per line

The meaning of the columns in *sample_hypergraph.txt* are:

- First column: Index of node 

- Second column: Index of hyperedge

Note that the index must start from 1 (NOT ZERO)

### Output
The local quad clustering coefficient vector

### Example
The sample data (*sample_hypergraph.txt*) is a hypergraph with five nodes and four hyperedges as
![nondirected_example](https://github.com/Gyeong-GyunHa/qch/assets/25861047/a79881c5-c4ca-4d00-8ae3-e174860652f6)

The quad clustering coefficient for a few nodes can be calculated  as :

Node 1: Node 1 is connected to hyperedges 1 and 2, and they have a common connection to nodes 3 and 4, making a total of two quads (numerator is 2.) And the maximum possible number of quads with the hyperedges 1 and 2 is 2 (denominator.) Therefore, the quad clustering coefficient is $`C^{\rm q}_{1}=\frac{2}{2}=1`$.

Node 3: Node 3 is connected to hyperedges 1, 2 and 3. There are 2 common neighbors 1 and 4 for the hyperedges 1 and 2, and 0 common neighbor for the hyperedges 2 and 3, and 1 common neighbor 2 for the hyperedges 1 and 3. Then we can compute the total number of quads centered on the node 3 is 3 (numerator of the quad clustering coefficient.) And the maximum possible number of quads centered on the node 3 is six (denominator): two each for hyperedges 1 and 2, hyperedges 2 and 3, and hyperedges 1 and 3. Therefore, the quad clustering coefficient is $`C^{\rm q}_{3}=\frac{1}{2}`$.

We can get the same result using this code. Here is the result file *quadclustering_dist.txt*

>   1.00000000    
>  0.500000000    
>  0.500000000    
>  0.500000000    
>   0.00000000    
>

This indicates that the quad clustering coefficients for nodes 1, 2, 3, 4, and 5 are 1, 1/2, 1/2, 1/2, and 0, respectively.


## The quad clustering coefficient in a directed hypergraph

### Usage
```zsh
> git clone https://github.com/Gyeong-GyunHa/qch
> cd qch
> gfortran dir_quadclustering.f90
```

This code doesn't support parallel computation for calculating the directed quad clustering coefficient but has parallel computation in the preprocessing step before calculating the coefficient.
To use parallel computation, you can compile it with the following command:
```zsh
> gfortran -fopenmp dir_quadclustering.f90
```

Once the application has been compiled you can run it using the following command:
```zsh
> ./a.out
```

### Input
The incidence list of a directed hypergraph with one connection information per line

The meaning of the columns in *dir_hyperg_sample.txt* are:

- First column: Index of node 

- Second column: Index of hyperedge

- Third column: Direction information of the connection (1 denotes the direction from node to hyperedge, 2 indicates the direction from hyperedge to node. To represent a bidirectional connection, use two lines. - *See sample data*)


### Output
The local directed quad clustering coefficient vector


## Reference

> G.-G. Ha, I. Neri, & A. Annibale (2023). Clustering coefficients for networks with higher order interactions. arXiv preprint arXiv:2311.08563.

Contact for inquiries about the code: gyeonggyun.ha_at_gmail.com
Gyeong-Gyun Ha, KCL, 2023




