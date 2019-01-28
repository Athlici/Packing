Overview
--------

![Packing](https://github.com/Athlici/Packing/blob/master/Packing.png)

Ever needed to arrange irregular objects with rotations and translations into as small a container as possible?

Worry no more for I have written the code to do just that so you don't have to.

Given a set of irregular 2D objects this program finds an arrangement such that the size of their container is minimized. This can be useful to optimize material or space usage in cutting and packing applications.
The objects can be unions of convex polygons as well as convex and concave circle segments which are subject to continuous translation and rotations.
The container can currently be a disk subject to scaling.

This program is at a Proof of Concept stage, but given interest in it I'm willing to assist and extend it. For example the infrastructure to extend this to an intersection of disks and half-spaces for the container or other scenarios and support for different objective functions is either already partially present or can be added easily.

So, how is this different from these great projects doing packing/nesting:

* [SVGnest](https://github.com/Jack000/SVGnest)
* [libnest2d](https://github.com/tamasmeszaros/libnest2d)
* [2D-Bin-Packing](https://github.com/mses-bly/2D-Bin-Packing)

All of them can only handle a discrete (and usually small) number of possible rotations. In return (and through heuristics) they can solve bigger problem instances.

Installation
------------
This project uses [Ipopt](https://projects.coin-or.org/Ipopt) for the actual numerical calculations. As such it and it's dependencies, in particular a sparse matrix library are required (and if that isn't ma57 but mumps change the line in main.cpp). Follow it's guide or pray your distro has a package for it.

Furthermore it uses the Eigen matrix library, so get that one as well.

After that building should be as simple as running

    make

In theory there should be no reason this couldn't run under windows, although getting ipopt to compile there already proves your skills superior to mine.

Usage
-----
Data In- and Export is currently implemented via XML files. The test.xml example contains all functions currently implemented and is given as first and only parameter when running the program. The output will be a copy of the input with the parameters of the best solution found added in out.xml. This solution can be visualized with the Mathematica notebook in this repo.

The defining points of the basic shapes (convex polygon, circle segment, hat) are expected in mathematically positive (counter-clockwise) order with the numbering for the latter two starting at the beginning of the line cutting of a circle.

TODO: Infographics

How it works
------------

Given the objects the program builds their distance functions based on and as described in this [paper](https://pdfs.semanticscholar.org/8c7b/92a7379af4c87b4cd5900745d48d62fbd954.pdf). This results in a lot of linear functions in a Min-Max-Tree with the constraint that the result of the tree is non-negative (no overlap).

Choosing a branch of this tree on which Ipopt can work is done by starting with a heuristic solution which is then successively refined as long as it reaches a new branching of the tree.

The heuristic I developed for this starts by finding a bounding circle for every object and then stringing these along a chain, inserting a circle at the first position where a circle can fit while touching both its predecessor and successor.

![Circles](https://github.com/Athlici/Packing/blob/master/Circles.png)

Known Bugs
----------

Somehow Ipopt doesn't want to work sometimes (regularly) so no solution is found. Restarting the program usually fixes this.
