# **M**agnetic **C**oordinate **T**ransformer

Given poloidal magnetic flux $\psi_{p}$ (analytic or numeric), calculate its Boozer coordinate representation.

## Quick Start

```bash
# also clone submodules
git clone https://github.com/12ff54e/mct --recurse-submodules
# configure & build
cmake -B ./build -DCMAKE_BUILD_TYPE=Release
cmake --build ./build
# convert specified gfile to spdata.dat, all options are optional
./build/mct path/to/gfile --output-path path/to/spdata.dat --radial 129 --poloidal 255 --use-si
```

## Known Issue

* Can not deal with negative poloidal current
