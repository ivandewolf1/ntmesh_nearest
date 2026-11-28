# ntmesh_primscan

Closest point on triangles, quads, n-gons, and polylines. Header: [`src/ntmesh_primscan.h`](../src/ntmesh_primscan.h).

`primScan` walks every list and keeps the closest. Nothing is built. An equal distance keeps the earlier primitive. Across lists the order is triangles, quads, n-gons, then lines. An empty mesh returns `-1` and leaves the hit unchanged.

```cpp
#include "ntmesh_primscan.h"

ntv::nvec3<float> pts[3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
int tris[] = {0, 1, 2};
ntm::PrimHit<float> hit;
int kind = ntm::primScan(pts, tris, 1, nullptr, 0, nullptr, nullptr, 0, nullptr, nullptr, 0,
                         ntv::nvec3<float>(0.2f, 0.2f, 1), hit);
// kind == ntm::ScanTri, hit.p is (0.2, 0.2, 0), hit.dist2 == 1
```

Pass `triangleIndices`, `quadIndices`, and the n-gon and line index and offset arrays from an `ntmesh`. A list you do not have is a null pointer and a count of 0.

---

## `PrimHit<T>`

| Field | Meaning |
|-------|---------|
| `p` | Closest point, the blend of `w` |
| `w` | Barycentric weights. An unused weight is an exact `0` |
| `dist2` | Squared distance from the query to `p` |
| `kind` | `ScanTri`, `ScanQuad`, `ScanNgon`, or `ScanLine`. `-1` is a miss |
| `prim` | Index in that list |
| `sub` | Quad half, n-gon fan triangle, or line segment. `0` on a triangle |

A quad fans from corner 0. Half `0` is corners `(0,1,2)`, half `1` is `(0,2,3)`, and `w` is on that half. An n-gon fans the same way: fan `f` is corners `(0, f+1, f+2)`. A face with fewer than three corners is skipped. A polyline is an open chain. Segment weights are `(1-s, s, 0)`. A zero-length segment returns its first vertex with weights `(1, 0, 0)`.

---

## The triangle

The query is projected onto the plane. That point is kept when every weight is non-negative. Otherwise the closest of the three edges is kept. A zero-area triangle has no plane, so the edges return the nearest vertex.

This is the direct form of the closest point in Christer Ericson, *Real-Time Collision Detection* (Morgan Kaufmann, 2005), section 5.1.5.

---

## See also

- [overview](overview.md) · [ntmesh_nearpc](ntmesh_nearpc.md) · [primscan sample](primscan_sample.md)
- [README](../README.md)
