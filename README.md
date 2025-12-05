# Min-Knapsack-Voronoi-Diagram-Algorithm
C++ implementation of an original algorithm that computes the Min Knapsack variant of the Voronoi diagram.
Given a set of points, their individual capacity and a total threshold to reach, it outputs the corresponding Min Knapsack Voronoi Diagram (MKVD).

For a more detailed definition of the structure read the article written by professor Plastria in 2016 **"Up- and downgrading the euclidean 1-median problem and knapsack Voronoi diagrams"**.

> Note: in this implementation the distance used is the Euclidean distance.

This project uses the following external libraries/components: 
- **MyGal** library by user pvigier to compute the Voronoi diagram of first order (https://github.com/pvigier/MyGAL)
- **SFML** library to visualize the diagram (https://www.sfml-dev.org/index.php)

## How to use
Compile command:
```bash
g++ -std=c++17 -IMyGAL/include Main.cpp MinKnapsack.cpp Point2D.cpp -o MyProgram -lsfml-graphics -lsfml-window -lsfml-system
```

Run command:
```bash
./MyProgram --file $(FILE) --visualize 1 --minKnapsack 1
```
The **flags** do the following:

- **file**: name of the *.txt* file that contains the capacity to reach, the positions and weights of each point in the following format:
    ```
    <capacity>

    <total number of points>

    <coordinate x> <coordinate y> <weight>

    ...

    <coordinate x> <coordinate y> <weight>

- **visualize**: 1/0 value to choose if to visualize the graph or not
- **minKnapsack**: 1/0 value to choose if to apply the minKnapsack procedure or to just visualize the order-1 Voronoi diagram of the input points

This work is still in progress, many tests are available in the Data section.

© 2025. This work is openly licensed via CC BY-NC