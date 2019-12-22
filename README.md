# polyscope-dgtal

This repository is an integration example between [polyscope](https://polyscope.run) and [DGtal](https://dgtal.org).

To use this project, just clone it (with submodule):

```
git clone --recursive https://github.com/dcoeurjo/polyscope-dgtal.git
```

Then you can compile the example file using:

```
mkdir buikd
cd build
cmake ..
make
```


The example code extracts the digital surface of an implicit shape and computes some differential quantities.


## Authors

* [David Coeurjolly](http://perso.liris.cnrs.fr/david.coeurjolly) (CNRS)
