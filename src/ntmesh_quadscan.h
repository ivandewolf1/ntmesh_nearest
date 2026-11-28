// Closest point on quads.
// 2026
// Ivan DeWolf
//
// method 0 splits a quad into triangles (0,1,2) and (0,2,3).
// methods 1, 2, and 3 use the bilinear patch
//   S(u,v) = (1-u)(1-v) p0 + u(1-v) p1 + uv p2 + (1-u)v p3.
// Newton starts at the center and at the planar (u,v). Eight steps.
// Algebraic clears u from the footpoint equations and solves a degree-5
// polynomial in v, then polishes each root. Planar projects onto the
// average plane and inverts the bilinear map. All three also test the
// four straight edges. w stores (u, v, 0).

#pragma once
#include "nvec3.h"
#include <cmath>
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
NTL_HD ntv::nvec3<T> bilinear(ntv::nvec3<T> a, ntv::nvec3<T> e, ntv::nvec3<T> f, ntv::nvec3<T> g, T u, T v) {
    return a + e * u + f * v + g * (u * v);
}

template <typename T>
NTL_HD void recordQuad(ntv::nvec3<T> q, ntv::nvec3<T> p, T u, T v,
                       PrimHit<T>& best, int kind, int prim, int sub) {
    takeHit(best, (q - p).length2(), p, ntv::nvec3<T>(u, v, T(0)), kind, prim, sub);
}

template <typename T>
NTL_HD void considerQuadEdge(ntv::nvec3<T> q, ntv::nvec3<T> a, ntv::nvec3<T> b,
                             T u0, T v0, T u1, T v1,
                             PrimHit<T>& best, int kind, int prim, int sub) {
    ntv::nvec3<T> ab = b - a;
    T len2 = ab.length2();
    T s = T(0);
    if (len2 != T(0)) {
        s = ab.dot(q - a) / len2;
        if (s < T(0)) s = T(0);
        else if (s > T(1)) s = T(1);
    }
    recordQuad(q, a + ab * s, u0 + (u1 - u0) * s, v0 + (v1 - v0) * s, best, kind, prim, sub);
}

template <typename T>
NTL_HD void considerQuadEdges(ntv::nvec3<T> q, ntv::nvec3<T> a, ntv::nvec3<T> b, ntv::nvec3<T> c, ntv::nvec3<T> d,
                              PrimHit<T>& best, int kind, int prim, int sub) {
    considerQuadEdge(q, a, b, T(0), T(0), T(1), T(0), best, kind, prim, sub);
    considerQuadEdge(q, b, c, T(1), T(0), T(1), T(1), best, kind, prim, sub);
    considerQuadEdge(q, c, d, T(1), T(1), T(0), T(1), best, kind, prim, sub);
    considerQuadEdge(q, d, a, T(0), T(1), T(0), T(0), best, kind, prim, sub);
}

template <typename T>
NTL_HD T cross2(T ax, T ay, T bx, T by) {
    return ax * by - ay * bx;
}

template <typename T>
NTL_HD void axisOf(ntv::nvec3<T> p, ntv::nvec3<T> o, ntv::nvec3<T> ex, ntv::nvec3<T> ey, T& x, T& y) {
    ntv::nvec3<T> d = p - o;
    x = d.dot(ex);
    y = d.dot(ey);
}

template <typename T>
NTL_HD ntv::nvec3<T> quadNormal(ntv::nvec3<T> a, ntv::nvec3<T> b, ntv::nvec3<T> c, ntv::nvec3<T> d) {
    ntv::nvec3<T> p[4] = {a, b, c, d};
    ntv::nvec3<T> n(T(0), T(0), T(0));
    for (int i = 0; i < 4; ++i) {
        ntv::nvec3<T> p0 = p[i];
        ntv::nvec3<T> p1 = p[(i + 1) % 4];
        n[0] += (p0[1] - p1[1]) * (p0[2] + p1[2]);
        n[1] += (p0[2] - p1[2]) * (p0[0] + p1[0]);
        n[2] += (p0[0] - p1[0]) * (p0[1] + p1[1]);
    }
    return n;
}

// Up to two (u, v) pairs from the bilinear inverse in the average plane.
template <typename T>
NTL_HD int planarUV(ntv::nvec3<T> q, ntv::nvec3<T> a, ntv::nvec3<T> b, ntv::nvec3<T> c, ntv::nvec3<T> d,
                    T& u0, T& v0, T& u1, T& v1) {
    ntv::nvec3<T> nrm = quadNormal(a, b, c, d);
    if (nrm.length2() < T(1e-20)) return 0;
    ntv::nvec3<T> mid = (a + b + c + d) * T(0.25);
    ntv::nvec3<T> on = q - nrm * (nrm.dot(q - mid) / nrm.length2());
    ntv::nvec3<T> ex = b - a;
    if (ex.length2() < T(1e-20)) ex = d - a;
    if (ex.length2() < T(1e-20)) ex = c - a;
    if (ex.length2() < T(1e-20)) return 0;
    ex.normalize();
    ntv::nvec3<T> ey = nrm.cross(ex);
    if (ey.length2() < T(1e-20)) return 0;
    ey.normalize();

    T b1x, b1y, c2x, c2y, d3x, d3y, hx, hy;
    axisOf(b, a, ex, ey, b1x, b1y);
    axisOf(c, a, ex, ey, c2x, c2y);
    axisOf(d, a, ex, ey, d3x, d3y);
    axisOf(on, a, ex, ey, hx, hy);
    T exx = b1x;
    T exy = b1y;
    T fx = d3x;
    T fy = d3y;
    T gx = c2x - b1x - d3x;
    T gy = c2y - b1y - d3y;

    T k2 = cross2(fx, fy, gx, gy);
    T k1 = cross2(fx, fy, exx, exy) - cross2(hx, hy, gx, gy);
    T k0 = -cross2(hx, hy, exx, exy);
    T span = exx * exx + exy * exy + fx * fx + fy * fy + gx * gx + gy * gy;
    T vs[2];
    int nv = 0;
    if (std::fabs(k2) <= span * T(1e-8)) {
        if (std::fabs(k1) > span * T(1e-10)) vs[nv++] = -k0 / k1;
    } else {
        T disc = k1 * k1 - T(4) * k2 * k0;
        if (disc >= -span * T(1e-8)) {
            if (disc < T(0)) disc = T(0);
            T sign = k1 < T(0) ? T(-1) : T(1);
            T qq = -T(0.5) * (k1 + sign * std::sqrt(disc));
            if (qq != T(0)) {
                vs[nv++] = qq / k2;
                T other = k0 / qq;
                if (std::fabs(other - vs[0]) > T(1e-6)) vs[nv++] = other;
            } else {
                vs[nv++] = T(0);
                vs[nv++] = -k1 / k2;
            }
        }
    }

    int n = 0;
    for (int i = 0; i < nv; ++i) {
        T v = vs[i];
        T sux = exx + v * gx;
        T suy = exy + v * gy;
        T den = sux * sux + suy * suy;
        T u = T(0);
        if (den > span * T(1e-12)) {
            T px = hx - v * fx;
            T py = hy - v * fy;
            u = (px * sux + py * suy) / den;
        }
        if (n == 0) { u0 = u; v0 = v; }
        else { u1 = u; v1 = v; }
        ++n;
    }
    return n;
}

template <typename T>
NTL_HD bool inSquare(T u, T v) {
    return u >= T(0) && u <= T(1) && v >= T(0) && v <= T(1);
}

template <typename T>
NTL_HD void clampSquare(T& u, T& v) {
    if (u < T(0)) u = T(0);
    else if (u > T(1)) u = T(1);
    if (v < T(0)) v = T(0);
    else if (v > T(1)) v = T(1);
}

// Eight Newton steps from one seed. Records the landing when it stays near the square.
template <typename T>
NTL_HD void newtonFrom(ntv::nvec3<T> q, ntv::nvec3<T> a, ntv::nvec3<T> e, ntv::nvec3<T> f, ntv::nvec3<T> g,
                       T u, T v, PrimHit<T>& best, int kind, int prim, int sub) {
    T scale = e.length2() + f.length2() + g.length2();
    for (int iter = 0; iter < 8; ++iter) {
        ntv::nvec3<T> su = e + g * v;
        ntv::nvec3<T> sv = f + g * u;
        ntv::nvec3<T> r = bilinear(a, e, f, g, u, v) - q;
        T fu = r.dot(su);
        T fv = r.dot(sv);
        T j00 = su.length2();
        T j11 = sv.length2();
        T j01 = su.dot(sv) + r.dot(g);
        T det = j00 * j11 - j01 * j01;
        if (std::fabs(det) <= scale * scale * T(1e-12)) break;
        T du = (-fu * j11 + fv * j01) / det;
        T dv = (-fv * j00 + fu * j01) / det;
        T step = T(1);
        T nu = u;
        T nv = v;
        for (int h = 0; h < 6; ++h) {
            nu = u + du * step;
            nv = v + dv * step;
            if (nu >= T(-1) && nu <= T(2) && nv >= T(-1) && nv <= T(2)) break;
            step *= T(0.5);
        }
        T adu = nu - u;
        T adv = nv - v;
        u = nu;
        v = nv;
        if (adu * adu + adv * adv < T(1e-14)) break;
        if (u < T(-2) || u > T(3) || v < T(-2) || v > T(3)) break;
    }
    if (u >= T(-1e-5) && u <= T(1) + T(1e-5) && v >= T(-1e-5) && v <= T(1) + T(1e-5)) {
        clampSquare(u, v);
        recordQuad(q, bilinear(a, e, f, g, u, v), u, v, best, kind, prim, sub);
    }
}

template <typename T>
NTL_HD T polyAt(const T* c, int deg, T x) {
    T y = c[deg];
    for (int i = deg - 1; i >= 0; --i) y = y * x + c[i];
    return y;
}

template <typename T>
NTL_HD T polyMag(const T* c, int deg) {
    T m = T(0);
    for (int i = 0; i <= deg; ++i) {
        T a = std::fabs(c[i]);
        if (a > m) m = a;
    }
    if (m < T(1)) m = T(1);
    return m;
}

template <typename T>
NTL_HD int sortUnique(T* r, int n) {
    for (int i = 1; i < n; ++i) {
        T v = r[i];
        int j = i;
        while (j > 0 && r[j - 1] > v) {
            r[j] = r[j - 1];
            --j;
        }
        r[j] = v;
    }
    int w = 0;
    for (int i = 0; i < n; ++i) {
        if (w == 0 || r[i] - r[w - 1] > T(1e-5)) r[w++] = r[i];
    }
    return w;
}

template <typename T>
NTL_HD int quadRoots01(T c0, T c1, T c2, T* roots) {
    T mag = std::fabs(c0) + std::fabs(c1) + std::fabs(c2);
    if (mag < T(1)) mag = T(1);
    if (std::fabs(c2) <= mag * T(1e-10)) {
        if (std::fabs(c1) <= mag * T(1e-10)) return 0;
        T x = -c0 / c1;
        if (x >= T(0) && x <= T(1)) {
            roots[0] = x;
            return 1;
        }
        return 0;
    }
    T disc = c1 * c1 - T(4) * c2 * c0;
    if (disc < T(0)) return 0;
    T sign = c1 < T(0) ? T(-1) : T(1);
    T qq = -T(0.5) * (c1 + sign * std::sqrt(disc));
    int n = 0;
    if (qq != T(0)) {
        T x1 = qq / c2;
        T x2 = c0 / qq;
        if (x1 >= T(0) && x1 <= T(1)) roots[n++] = x1;
        if (x2 >= T(0) && x2 <= T(1) && (n == 0 || std::fabs(x2 - roots[0]) > T(1e-6))) roots[n++] = x2;
    } else if (c0 == T(0)) {
        roots[n++] = T(0);
        T x = -c1 / c2;
        if (x > T(0) && x <= T(1)) roots[n++] = x;
    }
    return n;
}

// Real roots of c[0] + c[1] x + ... + c[deg] x^deg on [0, 1].
template <typename T>
NTL_HD int polyRoots01(const T* c, int deg, T* roots) {
    T mag = polyMag(c, deg);
    while (deg > 0 && std::fabs(c[deg]) <= mag * T(1e-10)) --deg;
    if (deg <= 0) return 0;
    if (deg == 1) {
        if (c[1] == T(0)) return 0;
        T x = -c[0] / c[1];
        if (x >= T(0) && x <= T(1)) {
            roots[0] = x;
            return 1;
        }
        return 0;
    }
    if (deg == 2) return quadRoots01(c[0], c[1], c[2], roots);

    T dc[6];
    for (int i = 1; i <= deg; ++i) dc[i - 1] = c[i] * T(i);
    T crit[6];
    int nc = polyRoots01(dc, deg - 1, crit);
    T pts[8];
    int np = 0;
    pts[np++] = T(0);
    for (int i = 0; i < nc && np < 7; ++i) pts[np++] = crit[i];
    pts[np++] = T(1);

    T found[12];
    int nf = 0;
    T tol = mag * T(1e-6);
    for (int i = 0; i < np; ++i) {
        if (std::fabs(polyAt(c, deg, pts[i])) <= tol && nf < 12) found[nf++] = pts[i];
    }
    for (int i = 0; i + 1 < np; ++i) {
        T lo = pts[i];
        T hi = pts[i + 1];
        if (hi - lo < T(1e-8)) continue;
        T yLo = polyAt(c, deg, lo);
        T yHi = polyAt(c, deg, hi);
        if (yLo == T(0) || yHi == T(0)) continue;
        if ((yLo > T(0) && yHi > T(0)) || (yLo < T(0) && yHi < T(0))) continue;
        for (int it = 0; it < 50; ++it) {
            T mid = (lo + hi) * T(0.5);
            T yMid = polyAt(c, deg, mid);
            if (yMid == T(0) || hi - lo < T(1e-12)) {
                lo = mid;
                break;
            }
            if ((yLo <= T(0) && yMid <= T(0)) || (yLo > T(0) && yMid > T(0))) {
                lo = mid;
                yLo = yMid;
            } else {
                hi = mid;
            }
        }
        if (nf < 12) found[nf++] = (lo + hi) * T(0.5);
    }
    int n = sortUnique(found, nf);
    if (n > 8) n = 8;
    for (int i = 0; i < n; ++i) roots[i] = found[i];
    return n;
}

template <typename T>
NTL_HD void considerQuadNewton(ntv::nvec3<T> q, ntv::nvec3<T> a, ntv::nvec3<T> b, ntv::nvec3<T> c, ntv::nvec3<T> d,
                               PrimHit<T>& best, int kind, int prim, int sub) {
    ntv::nvec3<T> e = b - a;
    ntv::nvec3<T> f = d - a;
    ntv::nvec3<T> g = a - b + c - d;
    newtonFrom(q, a, e, f, g, T(0.5), T(0.5), best, kind, prim, sub);
    T u0, v0, u1, v1;
    int n = planarUV(q, a, b, c, d, u0, v0, u1, v1);
    if (n >= 1) newtonFrom(q, a, e, f, g, u0, v0, best, kind, prim, sub);
    if (n >= 2) newtonFrom(q, a, e, f, g, u1, v1, best, kind, prim, sub);
    considerQuadEdges(q, a, b, c, d, best, kind, prim, sub);
}

template <typename T>
NTL_HD void considerQuadPlanar(ntv::nvec3<T> q, ntv::nvec3<T> a, ntv::nvec3<T> b, ntv::nvec3<T> c, ntv::nvec3<T> d,
                               PrimHit<T>& best, int kind, int prim, int sub) {
    ntv::nvec3<T> e = b - a;
    ntv::nvec3<T> f = d - a;
    ntv::nvec3<T> g = a - b + c - d;
    T u0, v0, u1, v1;
    int n = planarUV(q, a, b, c, d, u0, v0, u1, v1);
    if (n >= 1 && inSquare(u0, v0))
        recordQuad(q, bilinear(a, e, f, g, u0, v0), u0, v0, best, kind, prim, sub);
    if (n >= 2 && inSquare(u1, v1))
        recordQuad(q, bilinear(a, e, f, g, u1, v1), u1, v1, best, kind, prim, sub);
    considerQuadEdges(q, a, b, c, d, best, kind, prim, sub);
}

template <typename T>
NTL_HD void considerQuadAlgebraic(ntv::nvec3<T> q, ntv::nvec3<T> a, ntv::nvec3<T> b, ntv::nvec3<T> c, ntv::nvec3<T> d,
                                  PrimHit<T>& best, int kind, int prim, int sub) {
    ntv::nvec3<double> qd(q[0], q[1], q[2]);
    ntv::nvec3<double> ad(a[0], a[1], a[2]);
    ntv::nvec3<double> bd(b[0], b[1], b[2]);
    ntv::nvec3<double> cd(c[0], c[1], c[2]);
    ntv::nvec3<double> dd(d[0], d[1], d[2]);
    ntv::nvec3<double> e = bd - ad;
    ntv::nvec3<double> f = dd - ad;
    ntv::nvec3<double> g = ad - bd + cd - dd;
    ntv::nvec3<double> w = ad - qd;
    double d0 = e.dot(e);
    double d1 = 2.0 * e.dot(g);
    double d2 = g.dot(g);
    double n0 = -w.dot(e);
    double n1 = -(w.dot(g) + f.dot(e));
    double n2 = -f.dot(g);
    ntv::nvec3<double> l0 = w * d0 + e * n0;
    ntv::nvec3<double> l1 = f * d0 + w * d1 + g * n0 + e * n1;
    ntv::nvec3<double> l2 = f * d1 + w * d2 + g * n1 + e * n2;
    ntv::nvec3<double> l3 = f * d2 + g * n2;
    ntv::nvec3<double> r0 = f * d0 + g * n0;
    ntv::nvec3<double> r1 = f * d1 + g * n1;
    ntv::nvec3<double> r2 = f * d2 + g * n2;
    double coeff[6];
    coeff[0] = l0.dot(r0);
    coeff[1] = l0.dot(r1) + l1.dot(r0);
    coeff[2] = l0.dot(r2) + l1.dot(r1) + l2.dot(r0);
    coeff[3] = l1.dot(r2) + l2.dot(r1) + l3.dot(r0);
    coeff[4] = l2.dot(r2) + l3.dot(r1);
    coeff[5] = l3.dot(r2);

    double roots[8];
    int nr = polyRoots01(coeff, 5, roots);
    ntv::nvec3<T> eT = b - a;
    ntv::nvec3<T> fT = d - a;
    ntv::nvec3<T> gT = a - b + c - d;
    double scale = w.length2() + e.length2() + f.length2() + g.length2();
    for (int i = 0; i < nr; ++i) {
        double v = roots[i];
        ntv::nvec3<double> su = e + g * v;
        double den = su.length2();
        ntv::nvec3<double> aa = w + f * v;
        double u;
        if (den <= scale * 1e-14) {
            double ag = aa.dot(g);
            if (std::fabs(ag) <= scale * 1e-8) continue;
            u = -aa.dot(f) / ag;
        } else {
            u = -aa.dot(su) / den;
        }
        if (u < -1e-6 || u > 1.0 + 1e-6 || v < -1e-6 || v > 1.0 + 1e-6) continue;
        if (u < 0.0) u = 0.0;
        else if (u > 1.0) u = 1.0;
        if (v < 0.0) v = 0.0;
        else if (v > 1.0) v = 1.0;
        recordQuad(q, bilinear(a, eT, fT, gT, T(u), T(v)), T(u), T(v), best, kind, prim, sub);
        newtonFrom(q, a, eT, fT, gT, T(u), T(v), best, kind, prim, sub);
    }
    considerQuadEdges(q, a, b, c, d, best, kind, prim, sub);
}

template <typename T>
NTL_HD int quadScan(const ntv::nvec3<T>* pts,
                    const int* quads, int nquads,
                    int method,
                    ntv::nvec3<T> q, PrimHit<T>& hit) {
    PrimHit<T> best;

    if (pts && quads) {
        for (int i = 0; i < nquads; ++i) {
            const int* id = quads + i * 4;
            if (method == 0) {
                considerTri(q, pts[id[0]], pts[id[1]], pts[id[2]], best, ScanQuad, i, 0);
                considerTri(q, pts[id[0]], pts[id[2]], pts[id[3]], best, ScanQuad, i, 1);
            }
            if (method == 1) {
                considerQuadNewton(q, pts[id[0]], pts[id[1]], pts[id[2]], pts[id[3]], best, ScanQuad, i, 0);
            }
            if (method == 2) {
                considerQuadAlgebraic(q, pts[id[0]], pts[id[1]], pts[id[2]], pts[id[3]], best, ScanQuad, i, 0);
            }
            if (method == 3) {
                considerQuadPlanar(q, pts[id[0]], pts[id[1]], pts[id[2]], pts[id[3]], best, ScanQuad, i, 0);
            }
        }
    }
    if (best.kind < 0) return -1;
    hit = best;
    return best.kind;
}

} // namespace ntm
