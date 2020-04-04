# polyscope-dgtal

This repository is an integration example between [polyscope](https://polyscope.run) and [DGtal](https://dgtal.org).

To use this project, just clone it (with submodule):

```
git clone --recursive https://github.com/dcoeurjo/PolyscopeCNC.git
```

Then you can compile the example file using:

```
mkdir build
cd build
cmake ..
make
```
Consider using  `cmake -DCMAKE_BUILD_TYPE=Release` for faster computation.

*Note*: to compile DGtal, you would need [boost](boost.org) (only headers) and  [zlib](https://www.zlib.net).


The example code extracts the digital surface of an implicit shape and computes some differential quantities.


## Authors

* [David Coeurjolly](http://perso.liris.cnrs.fr/david.coeurjolly) (CNRS)


![](http://polyscope.run/media/teaser.svg) ![](https://dgtal.org/img/logoDGtal.svg)
