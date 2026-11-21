# ntmesh_nearpc

Nearest point in a point cloud. Header: [`src/ntmesh_nearpc.h`](../src/ntmesh_nearpc.h).

```cpp
#include "ntmesh_nearpc.h"

ntv::nvec3<float> pts[3] = {
    {0, 0, 0},
    {2, 0, 0},
    {0, 3, 0},
};
ntm::NearPC<float> cloud;
cloud.build(pts, 3);
float d2 = 0.f;
int hit = cloud.nearest(ntv::nvec3<float>(1, 0, 0), d2);  // 0, d2 == 1
cloud.clear();
```

`build` reads `pts` only while it runs. The finished tree keeps its own copy of each point in the leaf box, so the input array can go away after `build` returns. See [overview](overview.md) for the split and the threading rule.

---

## `NearNode<T>`

One node in the flat array.

| Field | Meaning |
|-------|---------|
| `mn`, `mx` | Axis-aligned box. On a leaf, both are the point |
| `left`, `right` | Child slots. `-1` on a leaf |
| `point` | Index into the array passed to `build`. `-1` on an internal node |

Children are stored in the same array. The left child of slot `s` is at `s + 1`. The right child follows the whole left subtree.

---

## `NearPC<T>`

| Field | Meaning |
|-------|---------|
| `nodes` | The node array, or null |
| `count` | Number of nodes. `2n - 1` after a successful build of `n` points, otherwise `0` |

A plain copy shares `nodes`. `clear()` frees that block. Call it once.

---

## `clear()`

Frees `nodes`, sets the pointer to null, and sets `count` to 0. Safe on an empty cloud.

---

## `build(pts, n)`

Replaces any tree already held. A null `pts` or an `n` less than 1 leaves the cloud empty.

On success the node count is `2n - 1`. The build picks the longer axis of the current box and partitions the points so the lower half by count is the left child. The two halves of a range of 1024 points or more build on separate cores when more than one hardware thread is available.

If the allocation fails, the cloud is left empty.

---

## `nearest(q, dist2)`

Returns the index of the nearest sample and writes the squared distance to `dist2`.

| Case | Result |
|------|--------|
| empty cloud | `-1`, and `dist2` is left as the caller passed it |
| one nearest point | that point's index, and its squared distance |
| two points at the same distance | the lower index |

The walk uses a fixed stack of 64 nodes. The median split keeps the tree shallow; a cloud that fits in memory stays inside that stack.

`nearest` is const. Concurrent calls on one finished tree are safe. Finish `build` and `clear` before any thread enters `nearest`.

```cpp
float keep = 5.f;
int miss = cloud.nearest(ntv::nvec3<float>(0, 0, 0), keep);
// empty cloud: miss == -1 and keep is still 5
```

---

## `nearBoxDist2(node, q)`

Squared distance from `q` to the node's box. Zero when `q` is inside the box. `nearest` uses it to skip a subtree that cannot improve the best hit.

---

## Host and CUDA

`build` and `clear` are host-only. They use allocation, `std::nth_element`, and `std::thread`.

`nearest` and `nearBoxDist2` carry `NTL_HD`. For a device compile, define it before the include:

```cpp
#define NTL_HD __host__ __device__
#include "ntmesh_nearpc.h"
```

The default is empty, which is the host build.

---

## See also

- [overview](overview.md) · [nearpc sample](nearpc_sample.md)
- [README](../README.md)
