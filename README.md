# Interpolated corrected curvature measures for polygonal surfaces

This repository is a demo of the paper:

**Interpolated corrected curvature measures for polygonal surfaces**,
Jacques-Olivier Lachaud, Pascal Romon, Boris Thibert, David
Coeurjolly,
Symposium on Geometry Processing, Computer Graphics Forum, 2020.

In this project, we compare our stable curvature measure with existing
approaches (Normal Cycles, Rusinkiewicz's formula, polynomial fitting
via Jet Fitting). Our closed form forumulas are given in the
[CorrectedNormalCurrentFormula.h](https://github.com/dcoeurjo/CorrectedNormalCurrent/blob/master/CorrectedNormalCurrentFormula.h)
file. We rely on [CGAL](https://cgal.org) for the Monge form via Jet Fitting
approach. This project also relies on the [DGtal](dgtal.org) library
for basic linear algebra.

If you would like to include our CNC forumlas in your project, you
would just need to copy the
[CorrectedNormalCurrentFormula.h](https://github.com/dcoeurjo/CorrectedNormalCurrent/blob/master/CorrectedNormalCurrentFormula.h)
header and provide your own `RealVector` and 3x3 matrix
operation implementations.


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

![](screenshot.png)

## Authors

* [Jacques-Olivier Lachaud](http://www.lama.univ-savoie.fr/pagesmembres/lachaud/People/LACHAUD-JO/person.html)
* [Boris Thibert](https://ljk.imag.fr/membres/Boris.Thibert/)
* [Pascal Romon](https://perso.math.u-pem.fr/romon.pascal/)
* [David Coeurjolly](http://perso.liris.cnrs.fr/david.coeurjolly)
