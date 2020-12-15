# Interpolated corrected curvature measures for polygonal surfaces

This repository is a demo of the paper:

**Interpolated corrected curvature measures for polygonal surfaces**,
Jacques-Olivier Lachaud, Pascal Romon, Boris Thibert, David
Coeurjolly,
Symposium on Geometry Processing, Computer Graphics Forum, 2020.

``` bibtex
@article{cnc2020,
  author  =    {Jacques-Olivier Lachaud, Pascal Romon, Boris Thibert and David
Coeurjolly},
  title   =    {Interpolated corrected curvature measures for polygonal surfaces},
  journal =    {Computer Graphics Forum (Proceedings of Symposium on Geometry Processing 2020)},
  year    =    {2020},
  volume  =    {39},
  number  =    {5},
}
```

In this project, we compare our stable curvature measure with existing
approaches (Normal Cycles, Rusinkiewicz's formula, polynomial fitting
via Jet Fitting). Our closed form formulas are given in the
[CorrectedNormalCurrentFormula.h](https://github.com/dcoeurjo/CorrectedNormalCurrent/blob/master/CorrectedNormalCurrentFormula.h)
file. We rely on [CGAL](https://cgal.org) for the Monge form via Jet Fitting
approach. This project also relies on the [DGtal](dgtal.org) library
for basic linear algebra. The project heavily uses [polyscope](http://polyscope.run) for the visualization and UI.

If you would like to include our CNC estimators in your project, you
would just need to copy the
[CorrectedNormalCurrentFormula.h](https://github.com/dcoeurjo/CorrectedNormalCurrent/blob/master/CorrectedNormalCurrentFormula.h)
header and provide your own `RealVector` and 3x3 matrix
operation implementations.

As an example, we also provide a stand-alone implementation (BSD License 2.0) of some formulas using [eigen](https://eigen.tuxfamily.org) in [CorrectedNormalCurrentFormulaEigen.h](https://github.com/dcoeurjo/CorrectedNormalCurrent/blob/master/CorrectedNormalCurrentFormulaEigen.h).


To use this project, just clone it (with submodule):

```
git clone --recursive https://github.com/dcoeurjo/CorrectedNormalCurrent.git
```

Then you can compile the example file using:

```
mkdir build
cd build
cmake ..  -DCMAKE_BUILD_TYPE=Release
make
```

*Note*: to compile DGtal, you would need [boost](boost.org) (only headers) and  [zlib](https://www.zlib.net).


Once the code has been built, just run *simpleTest* on a triangulated
OBJ file:

```
simpleTest ../spot.obj
```

![](spot.gif)

## Authors

* [Jacques-Olivier Lachaud](http://www.lama.univ-savoie.fr/pagesmembres/lachaud/People/LACHAUD-JO/person.html)
* [Boris Thibert](https://ljk.imag.fr/membres/Boris.Thibert/)
* [Pascal Romon](https://perso.math.u-pem.fr/romon.pascal/)
* [David Coeurjolly](http://perso.liris.cnrs.fr/david.coeurjolly)

## License

GNU LGPL (see header files)
