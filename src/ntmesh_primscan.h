// Closest point on triangles, quads, n-gons, and polylines.
// 2026
// Ivan DeWolf
//
// a simplest-possible PAINFULLY slow algorithm, 
// tests every triangle, quad, n-gon, and polyline.
// One walk, nothing built. Equal distance keeps the earlier primitive.
// Order is triangles, quads, n-gons, then lines.
// Quads and n-gons are fans from corner 0. A polyline does not wrap.
// An empty mesh returns -1 and leaves the hit unchanged.
//
// A triangle projects the query onto its plane and keeps that point when
// every barycentric weight is non-negative. Otherwise it keeps the closest
// of the three edges. A zero-area triangle has no plane, so the edges
// return the nearest vertex. This is the direct form of the closest point
// in Christer Ericson, Real-Time Collision Detection (Morgan Kaufmann,
// 2005), section 5.1.5.

#pragma once
#include "nvec3.h"
#ifndef NTL_HD
#  define NTL_HD
#endif

namespace ntm {

enum ScanList { ScanTri = 0, ScanQuad = 1, ScanNgon = 2, ScanLine = 3 };

template <typename T>
struct PrimHit {
    ntv::nvec3<T> p;
    ntv::nvec3<T> w;
    T dist2 = T(0);
    int kind = -1;
    int prim = -1;
    int sub = 0;
};

template <typename T>
NTL_HD void takeHit(PrimHit<T>& best, T d, ntv::nvec3<T> p, ntv::nvec3<T> w, int kind, int prim, int sub) {
    if (best.kind >= 0 && !(d < best.dist2)) return;
    best.p = p;
    best.w = w;
    best.dist2 = d;
    best.kind = kind;
    best.prim = prim;
    best.sub = sub;
}

template <typename T>
NTL_HD void considerEdge(ntv::nvec3<T> q, ntv::nvec3<T> a, ntv::nvec3<T> b, int ia, int ib,
                         PrimHit<T>& best, int kind, int prim, int sub) {
    ntv::nvec3<T> ab = b - a;
    T len2 = ab.length2();
    T s = T(0);
    if (len2 != T(0)) {
        s = ab.dot(q - a) / len2;
        if (s < T(0)) s = T(0);
        else if (s > T(1)) s = T(1);
    }
    ntv::nvec3<T> p = a + ab * s;
    ntv::nvec3<T> w(T(0), T(0), T(0));
    w[ia] = T(1) - s;
    w[ib] = s;
    takeHit(best, (q - p).length2(), p, w, kind, prim, sub);
}

template <typename T>
NTL_HD void considerTri(ntv::nvec3<T> q, ntv::nvec3<T> a, ntv::nvec3<T> b, ntv::nvec3<T> c,
                        PrimHit<T>& best, int kind, int prim, int sub) {
    ntv::nvec3<T> ab = b - a;
    ntv::nvec3<T> ac = c - a;
    ntv::nvec3<T> n = ab.cross(ac);
    T n2 = n.length2();
    if (n2 != T(0)) {
        ntv::nvec3<T> proj = q - n * (n.dot(q - a) / n2);
        T d00 = ab.dot(ab);
        T d01 = ab.dot(ac);
        T d11 = ac.dot(ac);
        ntv::nvec3<T> ap = proj - a;
        T denom = d00 * d11 - d01 * d01;
        if (denom != T(0)) {
            T v = (d11 * ap.dot(ab) - d01 * ap.dot(ac)) / denom;
            T wv = (d00 * ap.dot(ac) - d01 * ap.dot(ab)) / denom;
            T u = T(1) - v - wv;
            if (u >= T(0) && v >= T(0) && wv >= T(0)) {
                ntv::nvec3<T> p = a * u + b * v + c * wv;
                takeHit(best, (q - p).length2(), p, ntv::nvec3<T>(u, v, wv), kind, prim, sub);
                return;
            }
        }
    }
    considerEdge(q, a, b, 0, 1, best, kind, prim, sub);
    considerEdge(q, a, c, 0, 2, best, kind, prim, sub);
    considerEdge(q, b, c, 1, 2, best, kind, prim, sub);
}

template <typename T>
NTL_HD int primScan(const ntv::nvec3<T>* pts,
                    const int* tris, int ntris,
                    const int* quads, int nquads,
                    const int* ngonIdx, const int* ngonOff, int nngons,
                    const int* lineIdx, const int* lineOff, int nlines,
                    ntv::nvec3<T> q, PrimHit<T>& hit) {
    PrimHit<T> best;
    if (pts && tris) {
        for (int i = 0; i < ntris; ++i) {
            const int* id = tris + i * 3;
            considerTri(q, pts[id[0]], pts[id[1]], pts[id[2]], best, ScanTri, i, 0);
        }
    }
    if (pts && quads) {
        for (int i = 0; i < nquads; ++i) {
            const int* id = quads + i * 4;
            considerTri(q, pts[id[0]], pts[id[1]], pts[id[2]], best, ScanQuad, i, 0);
            considerTri(q, pts[id[0]], pts[id[2]], pts[id[3]], best, ScanQuad, i, 1);
        }
    }
    if (pts && ngonIdx && ngonOff) {
        for (int i = 0; i < nngons; ++i) {
            int begin = ngonOff[i];
            int n = ngonOff[i + 1] - begin;
            const int* id = ngonIdx + begin;
            for (int f = 0; f < n - 2; ++f)
                considerTri(q, pts[id[0]], pts[id[f + 1]], pts[id[f + 2]], best, ScanNgon, i, f);
        }
    }
    if (pts && lineIdx && lineOff) {
        for (int i = 0; i < nlines; ++i) {
            int begin = lineOff[i];
            int n = lineOff[i + 1] - begin;
            for (int s = 0; s < n - 1; ++s)
                considerEdge(q, pts[lineIdx[begin + s]], pts[lineIdx[begin + s + 1]], 0, 1, best, ScanLine, i, s);
        }
    }
    if (best.kind < 0) return -1;
    hit = best;
    return best.kind;
}

} // namespace ntm
