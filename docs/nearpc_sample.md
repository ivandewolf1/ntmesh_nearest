# nearpc sample

`samples/nearpc` is the program that exercises [`ntmesh_nearpc`](ntmesh_nearpc.md). It reads `samplegeo/light_bunny.geo`, lays down rotated copies, finds the nearest bunny vertex for each random point, and joins that point to the vertex with a polyline.

Source: [`samples/nearpc/src/main.C`](../samples/nearpc/src/main.C). The executable is `build/nearpc_scene`.

The bunny is about a unit across. Centers are 2 apart, so neighboring copies have a gap between them. Each copy gets a random axis and a random angle. The random sequence starts at the same seed for every preset, so `s`, `m`, and `l` rebuild the same scene every run.

## Build and run

From the repository root, after the [README](../README.md) build:

```
./build/nearpc_scene
./build/nearpc_scene m
./build/nearpc_scene large --threads -1
./build/nearpc_scene bench
```

The `.geo` path is compiled in. The working directory does not matter. The files are gitignored.

A bad argument prints the usage line and exits 1. `--threads 0`, any negative count other than `-1`, and a second mode word are bad arguments.

## Arguments

One mode word, then an optional `--threads` count. The mode may follow the flag. No arguments means `s` on one thread.

| Word | Also | Grid | Queries | Query box | Written file |
|------|------|------|---------|-----------|--------------|
| `s` | `small` | 4×4 on z = 0 | 1,000 | [-4.5, 4.5]³ | `samples/nearpc/small.geo` |
| `m` | `medium` | 4×4×4 | 10,000 | [-4.5, 4.5]³ | `samples/nearpc/medium.geo` |
| `l` | `large` | 16×16×16 | 100,000 | [-16, 16]³ | `samples/nearpc/large.geo` |
| `b` | `bench` | all three, one after another | | | nothing |

Centers for a grid of size 4 are -3, -1, 1, 3 on each axis that the preset uses. The flat grid keeps z at 0. The 16³ grid runs from -15 to 15.

`medium` draws the twelve edges of its query box after the search. Those corners are not samples, so a nearest query cannot land on them. `small` and `large` leave the box undrawn. The lines that come from the bunny file itself stay in the mesh; the new polylines are added on top.

`large` prints a line as each z slice is merged. That grid is 4096 bunnies and about ten million points. The `.geo` is about a gigabyte.

## `--threads`

The flag sets how many threads run `nearest`. It does not change the tree build. `build` already splits large ranges across cores.

| Value | Queries |
|-------|---------|
| omitted, or `1` | one thread |
| `-1` | `hardware_concurrency`, or 1 when that count is 0 |
| `N` > 1 | `N` threads, and no more threads than there are queries |

Each thread writes only its own slice of the index and distance arrays. The mesh is updated after the queries return, on the calling thread.

Bench ignores the flag. If you pass it anyway, the first line of the report is `bench ignores --threads`, and the six runs still use one thread and then every core.

## Bench

`bench` builds each preset once, times every query on one thread, then times the same queries on every core. It compares the index and the squared distance of every query. The tree is released before the next preset, so the 16³ cloud is not held twice.

A line looks like this:

```
l  1 thread  query 68.086s  680.86 us/q
l  32 threads  query 4.606s  46.06 us/q  ratio 14.78  match
```

`us/q` is the pass time divided by the query count, in microseconds. `ratio` is the one-thread per-query time divided by the many-thread per-query time. `match` means every index and every squared distance agreed. The first disagreement prints `mismatch`, names that query on stderr, and the process exits 1 after the rest of the sizes finish.

The last line is the sum of the six query passes. Build time and file output are outside that sum. Bench writes no `.geo`.

One run on a 32-core machine:

```
s  1 thread  query 0.059s  58.71 us/q
s  32 threads  query 0.004s  3.63 us/q  ratio 16.18  match
m  1 thread  query 0.957s  95.68 us/q
m  32 threads  query 0.051s  5.12 us/q  ratio 18.70  match
l  1 thread  query 68.086s  680.86 us/q
l  32 threads  query 4.606s  46.06 us/q  ratio 14.78  match
6 runs        query 73.763s
```

The `l` lines above omit the sixteen `large z` progress lines that print while that grid is merged.

## A single scene

`s`, `m`, and `l` print the resolved scene, the query time, and the file:

```
medium: 4x4x4, 10000 queries, 2 threads
query 0.475s
points 173336 lines 10076 -> /path/to/samples/nearpc/medium.geo
```

The point count is the bunny copies, the query points, and, for `medium`, the eight box corners. The line count is the polylines that came in with the bunny, one segment per query, and the twelve box edges on `medium`.

## See also

- [overview](overview.md) · [ntmesh_nearpc](ntmesh_nearpc.md)
- [README](../README.md)
