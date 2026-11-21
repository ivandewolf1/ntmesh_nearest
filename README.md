# ntmesh_nearest

Nearest point in a point cloud. One header, `src/ntmesh_nearpc.h`, builds a flat tree over an array of `ntv::nvec3` and returns the index of the closest sample.

The header is the library. The sample in `samples/nearpc` is the program that places bunnies, queries them, and writes a `.geo`. Read the docs, then run the sample.

## Documentation

- [Overview](docs/overview.md) — the tree, the nearest contract, and what is safe to call from several threads
- [ntmesh_nearpc](docs/ntmesh_nearpc.md) — `NearPC`, `build`, `nearest`, `clear`
- [nearpc sample](docs/nearpc_sample.md) — `s`, `m`, `l`, `bench`, and `--threads`

## Sample

`build/nearpc_scene` reads `samplegeo/light_bunny.geo`. No arguments runs the small scene on one thread.

| Arg | Scene | Queries | Output |
|-----|-------|---------|--------|
| `s`, `small` | flat 4×4 | 1,000 | `samples/nearpc/small.geo` |
| `m`, `medium` | 4×4×4, plus the box edges | 10,000 | `samples/nearpc/medium.geo` |
| `l`, `large` | 16×16×16 | 100,000 | `samples/nearpc/large.geo` |
| `b`, `bench` | time all three, serial then every core | the three counts above | a report, and no `.geo` |

`--threads N` sets the query threads for `s`, `m`, and `l`. `1` is serial. `-1` uses every core. Bench ignores the flag. The [sample doc](docs/nearpc_sample.md) has the geometry, the report, and the exit codes.

## Build

`NTV_ROOT` is the ntv repository, the directory that contains `src/nvec3.h`. `NTMESH_ROOT` is the ntmesh repository, the directory that contains `src/ntmesh.h`.

A `-D` argument wins, then the environment variable of the same name, then a sibling checkout (`../ntv/ntv` and `../ntmesh`). Configuration stops and names the variable when a root is still missing.

```
cmake -S . -B build
cmake --build build
./build/ntmesh_nearest_tests
./build/nearpc_scene
./build/nearpc_scene bench
```

Tests live in `test/test_near.C`, driven by `test/main.C`.
