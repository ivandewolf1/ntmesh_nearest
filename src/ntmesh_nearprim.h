// Closest point on triangles, quads, n-gons, and polylines.
// 2026
// Ivan DeWolf
//
// A flat median-split tree. Each leaf stores its own corners, so the mesh
// can be cleared after build. A triangle leaf holds three positions. A
// segment holds two. Internal nodes hold the box and the children.
// Equal distance keeps the earlier primitive. Order is triangles, quads,
// n-gons, then lines. Quads and n-gons are fans from corner 0. A polyline
// does not wrap. An empty tree returns -1 and leaves the hit unchanged.
//
// The walk visits the closer child first. A node whose nearest corner is
// strictly farther than the best hit is skipped. The search radius also
// shrinks to a box's farthest corner, since some point of that box is at
// most that far. The leaf test is considerTri or considerEdge.
//
// query(q, hit) is one point. query(queries, n, hits) Morton-sorts the
// list, then walks it on every core. Call clear() when finished. A plain
// copy shares the arrays.
//
// Later: a device build can fill these same nodes. A 4-wide tree for CPU
// queries would be a separate header.


#pragma once
#include "ntmesh_primscan.h"
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <thread>
#include <vector>
#ifndef NTL_HD
#  define NTL_HD
#endif

namespace ntm {

template <typename T>
struct NearPrimNode {
    ntv::nvec3<T> mn;
    ntv::nvec3<T> mx;
    int left = -1;
    int right = -1;
    int leaf = -1;
};

template <typename T>
struct NearPrimLeaf {
    ntv::nvec3<T> a;
    ntv::nvec3<T> b;
    ntv::nvec3<T> c;
    int kind = -1;
    int prim = -1;
    int sub = 0;
    int ncorn = 0;
};

template <typename T>
struct NearPrim {
    NearPrimNode<T>* nodes = nullptr;
    NearPrimLeaf<T>* leaves = nullptr;
    int count = 0;
    int nleaves = 0;

    void clear();
    void build(const ntv::nvec3<T>* pts,
               const int* tris, int ntris,
               const int* quads, int nquads,
               const int* ngonIdx, const int* ngonOff, int nngons,
               const int* lineIdx, const int* lineOff, int nlines);
    // Kind of the closest primitive, or -1. A miss leaves hit unchanged.
    NTL_HD int query(ntv::nvec3<T> q, PrimHit<T>& hit) const;
    // Host list. Hits keep the caller's order. A miss leaves that hit unchanged.
    void query(const ntv::nvec3<T>* queries, int n, PrimHit<T>* hits) const;
};

template <typename T>
NTL_HD void nearPrimBoxDist2(const NearPrimNode<T>& nd, ntv::nvec3<T> q, T& dMin, T& dMax) {
    T lo = T(0);
    T hi = T(0);
    for (int a = 0; a < 3; ++a) {
        T u = nd.mn[a] - q[a];
        T v = q[a] - nd.mx[a];
        T near = u > v ? u : v;
        if (near < T(0)) near = T(0);
        T far = u < v ? u : v;
        lo += near * near;
        hi += far * far;
    }
    dMin = lo;
    dMax = hi;
}

template <typename T>
NTL_HD void keepPrim(PrimHit<T>& best, const PrimHit<T>& cand) {
    if (cand.kind < 0) return;
    if (best.kind < 0 || cand.dist2 < best.dist2) {
        best = cand;
        return;
    }
    if (cand.dist2 > best.dist2) return;
    if (cand.kind < best.kind
        || (cand.kind == best.kind && cand.prim < best.prim)
        || (cand.kind == best.kind && cand.prim == best.prim && cand.sub < best.sub))
        best = cand;
}

template <typename T>
NTL_HD void applyLeaf(const NearPrimLeaf<T>& lf, ntv::nvec3<T> q, PrimHit<T>& best) {
    PrimHit<T> local;
    if (lf.ncorn <= 2)
        considerEdge(q, lf.a, lf.b, 0, 1, local, lf.kind, lf.prim, lf.sub);
    else
        considerTri(q, lf.a, lf.b, lf.c, local, lf.kind, lf.prim, lf.sub);
    keepPrim(best, local);
}

template <typename T>
NTL_HD int NearPrim<T>::query(ntv::nvec3<T> q, PrimHit<T>& hit) const {
    if (!nodes || count <= 0 || !leaves) return -1;
    PrimHit<T> best;
    T bound = T(0);
    int bounded = 0;
    int stack[64];
    int sp = 0;
    stack[sp++] = 0;
    while (sp > 0) {
        const NearPrimNode<T>& nd = nodes[stack[--sp]];
        T dMin = T(0);
        T dMax = T(0);
        nearPrimBoxDist2(nd, q, dMin, dMax);
        if (bounded && dMin > bound) continue;
        if (!bounded || dMax < bound) {
            bound = dMax;
            bounded = 1;
        }
        if (nd.leaf >= 0) {
            applyLeaf(leaves[nd.leaf], q, best);
            if (best.kind >= 0 && best.dist2 < bound) bound = best.dist2;
            continue;
        }
        T lMin = T(0), lMax = T(0), rMin = T(0), rMax = T(0);
        nearPrimBoxDist2(nodes[nd.left], q, lMin, lMax);
        nearPrimBoxDist2(nodes[nd.right], q, rMin, rMax);
        int lok = !bounded || lMin <= bound;
        int rok = !bounded || rMin <= bound;
        if (lok && lMax < bound) bound = lMax;
        if (rok && rMax < bound) bound = rMax;
        lok = lMin <= bound;
        rok = rMin <= bound;
        if (lok && rok) {
            int leftCloser = lMin < rMin || (lMin == rMin && lMax <= rMax);
            int first = leftCloser ? nd.left : nd.right;
            int second = leftCloser ? nd.right : nd.left;
            if (sp < 63) stack[sp++] = second;
            if (sp < 64) stack[sp++] = first;
        } else if (lok) {
            if (sp < 64) stack[sp++] = nd.left;
        } else if (rok) {
            if (sp < 64) stack[sp++] = nd.right;
        }
    }
    if (best.kind < 0) return -1;
    hit = best;
    return best.kind;
}

template <typename T>
void growPrim(ntv::nvec3<T>& mn, ntv::nvec3<T>& mx, ntv::nvec3<T> p) {
    for (int a = 0; a < 3; ++a) {
        if (p[a] < mn[a]) mn[a] = p[a];
        if (p[a] > mx[a]) mx[a] = p[a];
    }
}

template <typename T>
void leafBounds(const NearPrimLeaf<T>& lf, ntv::nvec3<T>& mn, ntv::nvec3<T>& mx) {
    mn = lf.a;
    mx = lf.a;
    growPrim(mn, mx, lf.b);
    if (lf.ncorn >= 3) growPrim(mn, mx, lf.c);
}

template <typename T>
ntv::nvec3<T> leafCentroid(const NearPrimLeaf<T>& lf) {
    ntv::nvec3<T> c = lf.a + lf.b;
    T n = T(2);
    if (lf.ncorn >= 3) {
        c = c + lf.c;
        n = T(3);
    }
    return c * (T(1) / n);
}

template <typename T>
void rangeBounds(const NearPrimLeaf<T>* leaves, const int* order, int begin, int end,
                 ntv::nvec3<T>& mn, ntv::nvec3<T>& mx) {
    leafBounds(leaves[order[begin]], mn, mx);
    for (int i = begin + 1; i < end; ++i) {
        ntv::nvec3<T> a, b;
        leafBounds(leaves[order[i]], a, b);
        growPrim(mn, mx, a);
        growPrim(mn, mx, b);
    }
}

template <typename T>
void buildPrimRange(NearPrimNode<T>* nodes, NearPrimLeaf<T>* leaves, ntv::nvec3<T>* cents,
                    int* order, int slot, int begin, int end, int threads) {
    int count = end - begin;
    NearPrimNode<T>& nd = nodes[slot];
    rangeBounds(leaves, order, begin, end, nd.mn, nd.mx);
    if (count == 1) {
        nd.left = -1;
        nd.right = -1;
        nd.leaf = order[begin];
        return;
    }
    int axis = 0;
    T span = nd.mx[0] - nd.mn[0];
    for (int a = 1; a < 3; ++a) {
        T s = nd.mx[a] - nd.mn[a];
        if (s > span) {
            span = s;
            axis = a;
        }
    }
    int mid = begin + count / 2;
    std::nth_element(order + begin, order + mid, order + end, [&](int ia, int ib) {
        return cents[ia][axis] < cents[ib][axis];
    });
    int leftCount = mid - begin;
    if (leftCount < 1) leftCount = 1;
    if (leftCount >= count) leftCount = count - 1;
    mid = begin + leftCount;
    int leftSlot = slot + 1;
    int rightSlot = slot + 2 * leftCount;
    nd.left = leftSlot;
    nd.right = rightSlot;
    nd.leaf = -1;
    if (threads > 1 && count >= 1024) {
        int rightThreads = threads / 2;
        if (rightThreads < 1) rightThreads = 1;
        std::thread worker([&]() {
            buildPrimRange(nodes, leaves, cents, order, rightSlot, mid, end, rightThreads);
        });
        buildPrimRange(nodes, leaves, cents, order, leftSlot, begin, mid, threads - rightThreads);
        worker.join();
    } else {
        buildPrimRange(nodes, leaves, cents, order, leftSlot, begin, mid, 1);
        buildPrimRange(nodes, leaves, cents, order, rightSlot, mid, end, 1);
    }
}

template <typename T>
int countPrimLeaves(const int* tris, int ntris, const int* quads, int nquads,
                    const int* ngonIdx, const int* ngonOff, int nngons,
                    const int* lineIdx, const int* lineOff, int nlines) {
    int n = 0;
    if (tris && ntris > 0) n += ntris;
    if (quads && nquads > 0) n += nquads * 2;
    if (ngonIdx && ngonOff && nngons > 0) {
        for (int i = 0; i < nngons; ++i) {
            int k = ngonOff[i + 1] - ngonOff[i];
            if (k > 2) n += k - 2;
        }
    }
    if (lineIdx && lineOff && nlines > 0) {
        for (int i = 0; i < nlines; ++i) {
            int k = lineOff[i + 1] - lineOff[i];
            if (k > 1) n += k - 1;
        }
    }
    return n;
}

template <typename T>
void fillPrimLeaves(NearPrimLeaf<T>* leaves, const ntv::nvec3<T>* pts,
                    const int* tris, int ntris, const int* quads, int nquads,
                    const int* ngonIdx, const int* ngonOff, int nngons,
                    const int* lineIdx, const int* lineOff, int nlines) {
    int w = 0;
    if (pts && tris) {
        for (int i = 0; i < ntris; ++i) {
            const int* id = tris + i * 3;
            NearPrimLeaf<T>& lf = leaves[w++];
            lf.a = pts[id[0]];
            lf.b = pts[id[1]];
            lf.c = pts[id[2]];
            lf.kind = ScanTri;
            lf.prim = i;
            lf.sub = 0;
            lf.ncorn = 3;
        }
    }
    if (pts && quads) {
        for (int i = 0; i < nquads; ++i) {
            const int* id = quads + i * 4;
            NearPrimLeaf<T>& h0 = leaves[w++];
            h0.a = pts[id[0]];
            h0.b = pts[id[1]];
            h0.c = pts[id[2]];
            h0.kind = ScanQuad;
            h0.prim = i;
            h0.sub = 0;
            h0.ncorn = 3;
            NearPrimLeaf<T>& h1 = leaves[w++];
            h1.a = pts[id[0]];
            h1.b = pts[id[2]];
            h1.c = pts[id[3]];
            h1.kind = ScanQuad;
            h1.prim = i;
            h1.sub = 1;
            h1.ncorn = 3;
        }
    }
    if (pts && ngonIdx && ngonOff) {
        for (int i = 0; i < nngons; ++i) {
            int begin = ngonOff[i];
            int k = ngonOff[i + 1] - begin;
            const int* id = ngonIdx + begin;
            for (int f = 0; f < k - 2; ++f) {
                NearPrimLeaf<T>& lf = leaves[w++];
                lf.a = pts[id[0]];
                lf.b = pts[id[f + 1]];
                lf.c = pts[id[f + 2]];
                lf.kind = ScanNgon;
                lf.prim = i;
                lf.sub = f;
                lf.ncorn = 3;
            }
        }
    }
    if (pts && lineIdx && lineOff) {
        for (int i = 0; i < nlines; ++i) {
            int begin = lineOff[i];
            int k = lineOff[i + 1] - begin;
            for (int s = 0; s < k - 1; ++s) {
                NearPrimLeaf<T>& lf = leaves[w++];
                lf.a = pts[lineIdx[begin + s]];
                lf.b = pts[lineIdx[begin + s + 1]];
                lf.c = lf.a;
                lf.kind = ScanLine;
                lf.prim = i;
                lf.sub = s;
                lf.ncorn = 2;
            }
        }
    }
}

template <typename T>
void NearPrim<T>::clear() {
    std::free(nodes);
    std::free(leaves);
    nodes = nullptr;
    leaves = nullptr;
    count = 0;
    nleaves = 0;
}

template <typename T>
void NearPrim<T>::build(const ntv::nvec3<T>* pts,
                        const int* tris, int ntris,
                        const int* quads, int nquads,
                        const int* ngonIdx, const int* ngonOff, int nngons,
                        const int* lineIdx, const int* lineOff, int nlines) {
    clear();
    int n = countPrimLeaves<T>(tris, ntris, quads, nquads, ngonIdx, ngonOff, nngons, lineIdx, lineOff, nlines);
    if (!pts || n <= 0) return;
    leaves = static_cast<NearPrimLeaf<T>*>(std::malloc(sizeof(NearPrimLeaf<T>) * (size_t)n));
    nodes = static_cast<NearPrimNode<T>*>(std::malloc(sizeof(NearPrimNode<T>) * (size_t)(2 * n - 1)));
    int* order = static_cast<int*>(std::malloc(sizeof(int) * (size_t)n));
    ntv::nvec3<T>* cents = static_cast<ntv::nvec3<T>*>(std::malloc(sizeof(ntv::nvec3<T>) * (size_t)n));
    if (!leaves || !nodes || !order || !cents) {
        std::free(order);
        std::free(cents);
        clear();
        return;
    }
    fillPrimLeaves(leaves, pts, tris, ntris, quads, nquads, ngonIdx, ngonOff, nngons, lineIdx, lineOff, nlines);
    for (int i = 0; i < n; ++i) {
        order[i] = i;
        cents[i] = leafCentroid(leaves[i]);
    }
    int threads = (int)std::thread::hardware_concurrency();
    if (threads < 1) threads = 1;
    buildPrimRange(nodes, leaves, cents, order, 0, 0, n, threads);
    count = 2 * n - 1;
    nleaves = n;
    std::free(order);
    std::free(cents);
}

template <typename T>
std::uint32_t spreadPrim10(std::uint32_t v) {
    v &= 1023u;
    v = (v * 0x00010001u) & 0xFF0000FFu;
    v = (v * 0x00000101u) & 0x0F00F00Fu;
    v = (v * 0x00000011u) & 0xC30C30C3u;
    v = (v * 0x00000005u) & 0x49249249u;
    return v;
}

template <typename T>
std::uint32_t queryMorton(ntv::nvec3<T> q, ntv::nvec3<T> mn, ntv::nvec3<T> mx) {
    std::uint32_t code = 0;
    for (int a = 0; a < 3; ++a) {
        T span = mx[a] - mn[a];
        std::uint32_t u = 0;
        if (span > T(0)) {
            T t = (q[a] - mn[a]) / span * T(1023);
            if (t < T(0)) t = T(0);
            if (t > T(1023)) t = T(1023);
            u = (std::uint32_t)t;
        }
        code |= spreadPrim10<T>(u) << a;
    }
    return code;
}

template <typename T>
struct QueryOrd {
    std::uint32_t code;
    int index;
};

template <typename T>
void queryPrimSlice(const NearPrim<T>* tree, const ntv::nvec3<T>* queries, PrimHit<T>* hits,
                    const QueryOrd<T>* ord, int begin, int end) {
    for (int i = begin; i < end; ++i)
        tree->query(queries[ord[i].index], hits[ord[i].index]);
}

template <typename T>
void NearPrim<T>::query(const ntv::nvec3<T>* queries, int n, PrimHit<T>* hits) const {
    if (!queries || !hits || n <= 0) return;
    if (!nodes || count <= 0) return;
    QueryOrd<T>* ord = static_cast<QueryOrd<T>*>(std::malloc(sizeof(QueryOrd<T>) * (size_t)n));
    if (!ord) {
        for (int i = 0; i < n; ++i) query(queries[i], hits[i]);
        return;
    }
    ntv::nvec3<T> mn = queries[0];
    ntv::nvec3<T> mx = queries[0];
    for (int i = 1; i < n; ++i) growPrim(mn, mx, queries[i]);
    for (int i = 0; i < n; ++i) {
        ord[i].code = queryMorton<T>(queries[i], mn, mx);
        ord[i].index = i;
    }
    std::sort(ord, ord + n, [](const QueryOrd<T>& a, const QueryOrd<T>& b) {
        if (a.code != b.code) return a.code < b.code;
        return a.index < b.index;
    });
    int threads = (int)std::thread::hardware_concurrency();
    if (threads < 1) threads = 1;
    if (threads > n) threads = n;
    if (n < 1024) threads = 1;
    if (threads == 1) {
        queryPrimSlice(this, queries, hits, ord, 0, n);
    } else {
        int chunk = (n + threads - 1) / threads;
        std::vector<std::thread> pool;
        for (int t = 0; t < threads; ++t) {
            int begin = t * chunk;
            if (begin >= n) break;
            int end = begin + chunk;
            if (end > n) end = n;
            pool.emplace_back(queryPrimSlice<T>, this, queries, hits, ord, begin, end);
        }
        for (std::thread& th : pool) th.join();
    }
    std::free(ord);
}

} // namespace ntm
