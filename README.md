# **M**agnetic **C**oordinate **T**ransformer

Given poloidal magnetic flux $\psi_{p}$ (analytic or numeric), calculate its Boozer coordinate representation.

## Quick Start

```bash
# configure & build
cmake -B ./build -DCMAKE_BUILD_TYPE:STRING=Release
cmake --build ./build
# convert specified gfile to spdata.dat
./build/mct file-path-to-gfile
```
