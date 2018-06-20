## SimplForman3D - discrete Morse theory on simplicial complexes

SimplForman3D provides a set of methods for computing a discrete Morse complex over a scalar field defined on a simplicial complex.

The main features of SimplForman3D are:
- the capability of handling manifold and non-manifold simplicial complexes up to dimension 3
- a compact representation for encoding the underlying simplicial complex [2]
- a compact representation for encoding the Forman gradient [1]


## Requirements and Installing instructions

SimplForman3D has been tested under both Ubuntu Unix and MacOS systems.

The required libraries are:
- [boost](https://www.boost.org)
- [PHAT](https://github.com/blazs/phat)

Once all the required libraries have been installed SimplForman3D can be compiled by using Cmake simply typing the following commands

```
mkdir build
cd build
cmake ..
```

Before running the `./SimplForman3D` executable remember to copy the two files `tables/table2D.txt` and `tables/table3D.txt` in same folder of `./SimplForman3D`.

## Input Format

The library assumes to receive in input an OFF file modified for handling non-manifold simplicial complexes. The accepted format indicate the number of vertices and top simplices in the simplicial complex.

- Each vertex is represented by a list of coordinates. The last coordinate is used as field value.

- Each top simplex is represented by the number of vertices and the indices of such vertices in the previous list.

For visualization purposes, the first three coordinates of each vertex are used as x,y,z coordinates.

A simple example for a simplicial complex formed by 1 top edge and 1 triangle is

```
OFF
4 2 0
0.0 0.0 1.0
-1.0 0.0 2.0
0.0 1.0 3.0
0.3 0.5 4.0
1.0 0.0 5.0
2 0 1
3 0 2 3
```

## Quick start

The library provides three main functions for producing results for data analysis or visualization. The following examples assume to receive the scalar field in input by running the library as

```
./SimplForman3D simplComplex.off
```

#### Efficiently computing persistent homology

```c++
#include "forman/formangradient.h"

using namespace std;
int main(int argc, char* argv[])
{
    //read the input
    FormanGradient grad = FormanGradient(argc,argv);

    //compute the Forman gradient
    grad.computeFormanGradient(true);

    //compute persistent homology
    grad.computePersistentHomology();

    return 0;
}
```
The result is the list of persistence pairs written in a txt file (persistence_pairs.txt)

#### Visualizing critical points and the discrete Morse complex

```c++
#include "forman/formangradient.h"

using namespace std;
int main(int argc, char* argv[])
{
    //read the input
    FormanGradient grad = FormanGradient(argc,argv);

    //compute the Forman gradient
    grad.computeFormanGradient(true);

    //write the critical cells
    grad.outputCriticalPoints();

    //write the cells of the discrete Morse complex
    grad.outputDescendingMorse();

    return 0;
}
```
The results are written in different .vtk files that can be visualized with Paraview.

<img src="https://imgur.com/ydMZ7Kl.png" alt="Descending2cells" width="300px"/>
<img src="https://imgur.com/CY60cSe.png" alt="Descending1cells" width="285px"/>

## Attribution

If you use our library please consider referring our work

[1] K. Weiss, F. Iuricich, R. Fellegara, and L. De Floriani, “A primal/dual representation for discrete Morse complexes on tetrahedral meshes,” Comput. Graph. Forum, vol. 32, no. 3, pp. 361–370, 2013.

[2] D. Canino, L. De Floriani, and K. Weiss, “IA*: An adjacency-based representation for non-manifold simplicial shapes in arbitrary dimensions,” Comput. Graph., vol. 35, no. 3, pp. 747–753, 2011.
