# nearprim sample

`samples/nearprim` exercises `src/ntmesh_nearprim.h`. It reads `samplegeo/light_bunny.geo`, lays down the same rotated copies as the [nearpc sample](nearpc_sample.md), builds one median tree, and finds the closest point on a triangle, quad, n-gon, or polyline. The `.geo` draws a segment from the query to that point.

Source: [`samples/nearprim/src/main.C`](../samples/nearprim/src/main.C). The executable is `build/nearprim_scene`.

The grids, the query boxes, and the seed match `nearpc_scene`. Centers are 2 apart. `medium` draws the twelve box edges after the search, so those edges are not targets. Each leaf stores its own corners. The mesh can be released after `build`.

## Build and run

From the repository root, after the [README](../README.md) build:

```
./build/nearprim_scene
./build/nearprim_scene m
./build/nearprim_scene large
./build/nearprim_scene bench
```

The `.geo` path is compiled in. The working directory does not matter. The files are gitignored.

A bad argument prints the usage line and exits 1. A second word is a bad argument. There is no `--threads` flag. The batch query chooses the thread count.

## Arguments

One mode word. No arguments means `s`.

| Word | Also | Grid | Queries | Query box | Written file |
|------|------|------|---------|-----------|--------------|
| `s` | `small` | 4×4 on z = 0 | 1,000 | [-4.5, 4.5]³ | `samples/nearprim/small.geo` |
| `m` | `medium` | 4×4×4 | 10,000 | [-4.5, 4.5]³ | `samples/nearprim/medium.geo` |
| `l` | `large` | 16×16×16 | 100,000 | [-16, 16]³ | `samples/nearprim/large.geo` |
| `b` | `bench` | all three, one after another | | | nothing |

A list shorter than 1,024 queries stays on the calling thread. `small` is that list. `medium` and `large` use every core. The header comment on `query(queries, n, hits)` is the rule the sample prints.

The scene mode calls that batch query, prints the build and the query, then writes the file. The tree is released before the lines are added.

## Bench

`bench` builds each preset once. It times every query on one call at a time, then times `query(queries, n, hits)` on the same tree. It compares kind, primitive index, sub-index, squared distance, the hit point, and the weights. The tree and the mesh are released before the next preset.

`small` prints `1 thread` on both lines, because 1,000 queries never leave the calling thread. `medium` and `large` print every core on the second line.

One run of `s` on this machine:

```
small: 4x4, 1000 queries, 1 thread
prims 44176  tris 12176  quads 30896  ngons 1088  lines 16
build 0.037s  leaves 78544
query 0.025s
points 42832 lines 1016 -> /path/to/samples/nearprim/small.geo
```

A bench block is a build line, the one-at-a-time pass, then the batch pass. `small` stays on one thread for both passes:

```
s build <time>s  leaves <count>
s  1 thread  query <time>s  <us/q> us/q
s  1 thread  query <time>s  <us/q> us/q  ratio <serial/batch>  match
```

`ratio` is the one-at-a-time per-query time divided by the batch per-query time. `match` means every compared field agreed. The first disagreement prints `mismatch`, names that query on stderr as `kind:prim`, and the process exits 1 after the rest of the sizes finish. The last line is the sum of the six query passes. Build time is outside that sum. Bench writes no `.geo`.

## See also

- [overview](overview.md) · [primscan sample](primscan_sample.md) · [nearpc sample](nearpc_sample.md)
- [README](../README.md)
