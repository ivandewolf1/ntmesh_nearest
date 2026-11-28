# primscan sample

`samples/primscan` exercises [`ntmesh_primscan`](ntmesh_primscan.md). It reads `samplegeo/light_bunny.geo`, lays down the same rotated copies as the [nearpc sample](nearpc_sample.md), and walks every triangle, quad, n-gon, and polyline for each query. The `.geo` draws a segment from the query to that closest point.

Source: [`samples/primscan/src/main.C`](../samples/primscan/src/main.C). The executable is `build/primscan_scene`.

The grids, the query boxes, and the seed match `nearpc_scene`. Centers are 2 apart. `medium` draws the twelve box edges after the search, so those edges are not targets. Bench times `s` and `m` only. `l` is still a mode, and it is the long one.

## Build and run

From the repository root, after the [README](../README.md) build:

```
./build/primscan_scene
./build/primscan_scene m
./build/primscan_scene large --threads -1
./build/primscan_scene bench
```

The `.geo` path is compiled in. The working directory does not matter. The files are gitignored.

A bad argument prints the usage line and exits 1. `--threads 0`, any negative count other than `-1`, and a second mode word are bad arguments.

## Arguments

One mode word, then an optional `--threads` count. The mode may follow the flag. No arguments means `s` on one thread.

| Word | Also | Grid | Queries | Query box | Written file |
|------|------|------|---------|-----------|--------------|
| `s` | `small` | 4×4 on z = 0 | 1,000 | [-4.5, 4.5]³ | `samples/primscan/small.geo` |
| `m` | `medium` | 4×4×4 | 10,000 | [-4.5, 4.5]³ | `samples/primscan/medium.geo` |
| `l` | `large` | 16×16×16 | 100,000 | [-16, 16]³ | `samples/primscan/large.geo` |
| `b` | `bench` | `s` and `m`, serial then every core | | | nothing |

`--threads` is the query thread count. `1` is serial. `-1` is every core. Bench ignores the flag and says so when the flag was passed. Each thread writes only its own slice of the hit array. The mesh is updated after the queries return.

## What a query costs

There is no tree. Every query runs the closest-point test on every triangle, quad half, n-gon fan triangle, and line segment. Cost tracks that count times the query count.

One bunny is 2,761 polygons. The small grid is 16 copies, medium is 64, and large is 4,096. Against small, medium is 40 times the primitive-visits (4 times the copies, 10 times the queries). Large is 25,600 times small, and 640 times medium.

One run on this machine, one thread:

```
small: 4x4, 1000 queries, 1 thread
prims 44176  tris 12176  quads 30896  ngons 1088  lines 16
query 33.928s
points 42832 lines 1016 -> /path/to/samples/primscan/small.geo
```

Large is 25,600 times those primitive-visits. At this per-visit cost that is about 10 days on one thread (11,309,056 prims, 100,000 queries). Medium is 40 times small, about 23 minutes for the one-thread pass.

The point count is the bunny copies plus two points per query (the query and the hit), and eight box corners on `medium`. The line count adds one segment per query, and the twelve box edges on `medium`.

## Bench

`bench` builds small, then medium. Each one is timed on one thread, then on every core. It compares kind, primitive index, sub-index, squared distance, the hit point, and the weights. The mesh is released before the next preset. Large is left out.

A one-thread line from the small run, and the shape of the threaded line that follows it:

```
s  1 thread  query 33.928s  33928.00 us/q
s  N threads  query <time>s  <us/q> us/q  ratio <serial/threaded>  match
```

`ratio` is the one-thread per-query time divided by the many-thread per-query time. `match` means every compared field agreed. The first disagreement prints `mismatch`, names that query on stderr as `kind:prim`, and the process exits 1 after the rest of the sizes finish. The last line is the sum of the four query passes. Bench writes no `.geo`.

`l` is not part of that sum. On one thread it is the 10 day estimate above.

## See also

- [overview](overview.md) · [ntmesh_primscan](ntmesh_primscan.md) · [nearpc sample](nearpc_sample.md)
- [README](../README.md)
