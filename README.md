# ntmesh_nearest
ntmesh distance queries for pointclouds and polymeshes

`src/ntmesh_nearpc.h` builds a flat tree over a point array. Each node splits the longer axis at the halfway point, and the two halves build on separate cores. `nearest` returns the point index. Equal distances keep the lower index. Call `clear()` when finished. The build uses threads.

## Roots

`NTV_ROOT` is the ntv repository, the directory that contains `src/nvec3.h`. `NTMESH_ROOT` is the ntmesh repository, the directory that contains `src/ntmesh.h`. The scene programs read `samplegeo/light_bunny.geo` in this repo.

A `-D` argument wins, then the environment variable of the same name, then a sibling checkout (`../ntv/ntv` and `../ntmesh`). Configuration stops and names the variable when a root is still missing.

```
cmake -S . -B build
cmake --build build
./build/ntmesh_nearest_tests
```

`build/nearpc_scene` writes `build/test_nearpc.geo` when run from the repository root. `build/nearpc_scene_big` is the large grid.
