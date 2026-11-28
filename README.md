# ntmesh_nearest

Nearest point in a point cloud, and closest point on mesh primitives. `src/ntmesh_nearpc.h` builds a flat tree over an array of `ntv::nvec3` and returns the index of the closest sample. `src/ntmesh_primscan.h` walks triangles, quads, n-gons, and polylines. `src/ntmesh_nearprim.h` builds a median tree over those same primitives.

`samples/nearpc` places bunnies, queries the point cloud, and writes a `.geo`. `samples/primscan` uses the same grids and query points, and walks every primitive. Read the docs, then run a sample.

## Documentation

- [Overview](docs/overview.md) — the tree, the nearest contract, and what is safe to call from several threads
- [ntmesh_nearpc](docs/ntmesh_nearpc.md) — `NearPC`, `build`, `nearest`, `clear`
- [ntmesh_primscan](docs/ntmesh_primscan.md) — closest point on triangles, quads, n-gons, and polylines
- [nearpc sample](docs/nearpc_sample.md) — `s`, `m`, `l`, `bench`, and `--threads`
- [primscan sample](docs/primscan_sample.md) — the same scenes, walking every primitive
- [nearprim sample](docs/nearprim_sample.md) — the same scenes, answered by the median tree
- [quadscan sample](docs/quadscan_sample.md) — the same scenes, four closest-point methods on quads

## Sample

`build/nearpc_scene` reads `samplegeo/light_bunny.geo`. No arguments runs the small scene on one thread.

| Arg | Scene | Queries | Output |
|-----|-------|---------|--------|
| `s`, `small` | flat 4×4 | 1,000 | `samples/nearpc/small.geo` |
| `m`, `medium` | 4×4×4, plus the box edges | 10,000 | `samples/nearpc/medium.geo` |
| `l`, `large` | 16×16×16 | 100,000 | `samples/nearpc/large.geo` |
| `b`, `bench` | time all three, serial then every core | the three counts above | a report, and no `.geo` |

`--threads N` sets the query threads for `s`, `m`, and `l`. `1` is serial. `-1` uses every core. Bench ignores the flag. The [sample doc](docs/nearpc_sample.md) has the geometry, the report, and the exit codes.

`build/primscan_scene` takes the same scene words and writes `samples/primscan/<name>.geo`. Bench times `s` and `m` only. `l` is still a mode, and each query walks every triangle, quad, n-gon, and polyline, so that grid is much slower than the point-cloud scene. The line in the `.geo` runs from the query to the closest point on the mesh. The [primscan sample doc](docs/primscan_sample.md) has the scaling.

`build/quadscan_scene` places copies of `samplegeo/twistedCube.geo` on the same grids and query points, and writes `samples/quadscan/<name>.geo` from the triangle fan. Bench times that fan, Newton, the degree-5 polynomial, and the planar inverse on the small scene, on one thread. The [quadscan sample doc](docs/quadscan_sample.md) has the report.

`build/nearprim_scene` uses the same grids and query points and writes `samples/nearprim/<name>.geo`. There is no `--threads` flag. The batch query stays on one core below 1,024 points and uses every core from there up, so `small` is serial and `medium` and `large` are not. Bench times all three: one query at a time, then the batch. The [nearprim sample doc](docs/nearprim_sample.md) has the report.

## Build

`NTV_ROOT` is the ntv repository, the directory that contains `src/nvec3.h`. `NTMESH_ROOT` is the ntmesh repository, the directory that contains `src/ntmesh.h`.

A `-D` argument wins, then the environment variable of the same name, then a sibling checkout (`../ntv/ntv` and `../ntmesh`). Configuration stops and names the variable when a root is still missing.

```
cmake -S . -B build
cmake --build build
./build/ntmesh_nearest_tests
./build/nearpc_scene
./build/nearpc_scene bench
./build/primscan_scene
./build/nearprim_scene
./build/quadscan_scene
./build/quadscan_scene bench
```

Tests live in `test/test_near.C`, `test/test_primscan.C`, and `test/test_nearprim.C`, driven by `test/main.C`.
