# Overview

**ntmesh_nearest** answers two questions. Given a cloud of points, which point is closest to this query? And given triangles, quads, n-gons, and polylines, which primitive is closest, and where on it? Positions are `ntv::nvec3<T>` from [ntv](https://github.com/ivandewolf1/ntv). The sample that draws the point-cloud result uses [ntmesh](https://github.com/ivandewolf1/ntmesh) to read a bunny and write a `.geo`.

[`src/ntmesh_nearpc.h`](../src/ntmesh_nearpc.h) is the point cloud. Include it, build a tree, call `nearest`. Call `clear()` when you are done. [`src/ntmesh_primscan.h`](../src/ntmesh_primscan.h) walks triangles, quads, n-gons, and polylines. There is no tree.

## The tree

`build` lays the points out in one array of nodes. There is no pointer chasing and no per-node allocation. A cloud of `n` points becomes `2n - 1` nodes. Each node stores an axis-aligned box.

An internal node cuts along its longer side. The points on that axis are split by count, half on the left and the rest on the right. A leaf stores the index of one input point, and its box is that point.

The two halves of a split share nothing, so a large node builds its right half on another core while the caller builds the left. The split into threads starts once a side holds at least 1024 points. Smaller ranges stay on the calling thread.

## `nearest`

`nearest(q, dist2)` returns the index of the closest sample. The squared distance is written to `dist2`.

A tie keeps the lower index. An empty tree returns `-1` and leaves `dist2` alone.

The walk prunes any node whose box is already farther than the best point found so far. `nearBoxDist2` is that box test.

`nearest` only reads the finished node array. Its stack and its best distance live in the call. After `build` has returned, any number of threads may call `nearest` on the same tree. `build` and `clear` are the calls that change the array. Keep those on one thread, and finish them before the next query.

A plain copy of a `NearPC` copies the pointer. `clear()` frees that array. Call it once for a set of copies.

## What stays on the host

`build` and `clear` allocate, sort, and spawn threads. They are host-only.

`nearest` and `nearBoxDist2` are marked `NTL_HD`. Define `NTL_HD` as `__host__ __device__` before the include when a CUDA translation unit needs the query. Left undefined, the macro is empty.

## Docs

- [ntmesh_nearpc](ntmesh_nearpc.md) — types and calls
- [ntmesh_primscan](ntmesh_primscan.md) — closest point on triangles, quads, n-gons, and polylines
- [nearpc sample](nearpc_sample.md) — the bunny scenes, `bench`, and `--threads`
- [primscan sample](primscan_sample.md) — the same scenes, walking every primitive
- [nearprim sample](nearprim_sample.md) — the same scenes, answered by the median tree
- [quadscan sample](quadscan_sample.md) — the same scenes, four closest-point methods on quads

See the [README](../README.md) for the roots and the build.
