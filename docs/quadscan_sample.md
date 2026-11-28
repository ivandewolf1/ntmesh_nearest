# quadscan sample

`samples/quadscan` exercises `src/ntmesh_quadscan.h`. It reads `samplegeo/twistedCube.geo`, lays down rotated copies on the same grids as the [nearpc sample](nearpc_sample.md), and finds the closest point on a quad. The cube is six twisted quads. The `.geo` draws a segment from the query to the fan hit.

Source: [`samples/quadscan/src/main.C`](../samples/quadscan/src/main.C). The executable is `build/quadscan_scene`.

The grids, the query boxes, and the seed match `nearpc_scene`. Centers are 2 apart. `medium` draws the twelve box edges after the search, so those edges are not targets. `s`, `m`, and `l` use the triangle fan, method 0. Bench times all four methods on `small` only, on one thread.

## Build and run

From the repository root, after the [README](../README.md) build:

```
./build/quadscan_scene
./build/quadscan_scene m
./build/quadscan_scene bench
```

The `.geo` path is compiled in. The working directory does not matter. The files are gitignored.

A bad argument prints the usage line and exits 1. `--threads 0`, any negative count other than `-1`, and a second mode word are bad arguments.

## Arguments

One mode word, then an optional `--threads` count. The mode may follow the flag. No arguments means `s` on one thread, using the fan.

| Word | Also | Grid | Queries | Query box | Written file |
|------|------|------|---------|-----------|--------------|
| `s` | `small` | 4×4 on z = 0 | 1,000 | [-4.5, 4.5]³ | `samples/quadscan/small.geo` |
| `m` | `medium` | 4×4×4 | 10,000 | [-4.5, 4.5]³ | `samples/quadscan/medium.geo` |
| `l` | `large` | 16×16×16 | 100,000 | [-16, 16]³ | `samples/quadscan/large.geo` |
| `b` | `bench` | `small`, four methods, one thread | | | nothing |

`--threads` is the query thread count for `s`, `m`, and `l`. `1` is serial. `-1` is every core. Bench ignores the flag and says so when the flag was passed. Bench stays on one thread so the four times compare the solvers.

## Methods

| Method | What it solves |
|--------|----------------|
| `fan` | Two triangles, `(0,1,2)` and `(0,2,3)` |
| `newton` | Bilinear patch, eight Newton steps, then the four edges |
| `algebraic` | The same patch, by a degree-5 polynomial in `v` |
| `planar` | Project onto the average plane and invert the bilinear map, then the four edges |

`x fan` on a later line is that method's time divided by the fan time. `max d2` is the largest absolute gap in squared distance against the fan, over the 1,000 queries. `prim differ` counts queries whose winning quad index is not the fan's quad.

## Bench

`bench` builds the small grid once. It times the fan, then Newton, then the polynomial, then the planar inverse, on the same 1,000 queries. One run on this machine:

```
s: 96 quads, 1000 queries, 1 thread
s  fan         query 0.085s  85.20 us/q
s  newton      query 0.293s  293.11 us/q  3.44x fan  max d2 0.105  prim differ 67
s  algebraic   query 0.217s  217.31 us/q  2.55x fan  max d2 0.105  prim differ 67
s  planar      query 0.095s  95.00 us/q  1.12x fan  max d2 0.113  prim differ 73
4 runs        query 0.691s
```

The last line is the sum of the four query passes. Building the grid is outside that sum. Bench writes no `.geo`. A method that misses a query prints that index on stderr and the process exits 1.

`m` and `l` are not part of the bench. Each query still walks every quad.

## See also

- [overview](overview.md) · [primscan sample](primscan_sample.md) · [nearpc sample](nearpc_sample.md)
- [README](../README.md)
