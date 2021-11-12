# Harmonic-Coordinates-for-Character-Animation

This project was a final project for INF 585 - Computer Animation (2020-2021) course at Ecole Polytechnique de Paris. This project is implementation of the paper by Pushkar Joshi - [Harmonic-Coordinates-for-Character-Animation](https://www.cs.jhu.edu/~misha/Fall07/Papers/Joshi07.pdf). 

The details of our implementation are described in the [report](https://github.com/iradab/Harmonic-Coordinates-for-Character-Animation/blob/main/Harmonic_coordinates_report.pdf). 

## Description of the project

In computer animation, mesh deformation using
a grid allows for the manipulation of the space
containing the mesh, which renders certain
operations much easier. To deform the space
within the grid when grid points are moved,
many different functions have been proposed.
Here, we show an implementation of a method
called Harmonic coordinates, which allows for
smoother and more intuitive mesh point
manipulation using a cage.

## User Interface
Our first task was to implement a simple
interface for the user to interact with the 2D
board. We wrote code so that the user could
change the mouse click behavior using the
number row from 1 to 6. Here is what we
implemented and what number it is mapped to:
1. Place the mesh vertices, edges are
created between subsequent vertices,
and between the first and last one

2. Place the cage vertices, edges are
created between subsequent vertices,
and between the first and last one

3. Allows you to drag cage vertices. If one
of the cage vertices is selected by mouse
click, then the harmonic coordinates
related to that vertex are displayed. The
inside and boundary cells are displayed
by coloring them in blue and green
respectively. The intensity of the color is
dependent on the value of harmonic
coordinate.

4. Toggles the visibility of mesh and cage
edges

5. Runs all the grid calculations required to
compute harmonic coordinates.

6. Shows the grid points, with different
colors for grid points inside, outside and
at the boundary of the cage.


## Results 

Here is an example of using the program: 


https://user-images.githubusercontent.com/48281612/141531015-d81c5204-e742-4f34-b48c-30932bddf6ac.mp4



## Authors

* **Irada Bunyatova**     [iradab](https://github.com/iradab)
* **Gaël Van der Lee**       [Gvanderl](https://github.com/Gvanderl)
