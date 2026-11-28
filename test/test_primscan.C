#include <iostream>
#include "tests.h"
#include "ntmesh_primscan.h"

using ntv::nvec3;

static int tri(const nvec3<float>* pts, const int* idx, int n, nvec3<float> q, ntm::PrimHit<float>& hit) {
  return ntm::primScan(pts, idx, n, nullptr, 0, nullptr, nullptr, 0, nullptr, nullptr, 0, q, hit);
}
static int quad(const nvec3<float>* pts, const int* idx, int n, nvec3<float> q, ntm::PrimHit<float>& hit) {
  return ntm::primScan(pts, nullptr, 0, idx, n, nullptr, nullptr, 0, nullptr, nullptr, 0, q, hit);
}
static int ngon(const nvec3<float>* pts, const int* idx, const int* off, int n, nvec3<float> q, ntm::PrimHit<float>& hit) {
  return ntm::primScan(pts, nullptr, 0, nullptr, 0, idx, off, n, nullptr, nullptr, 0, q, hit);
}
static int line(const nvec3<float>* pts, const int* idx, const int* off, int n, nvec3<float> q, ntm::PrimHit<float>& hit) {
  return ntm::primScan(pts, nullptr, 0, nullptr, 0, nullptr, nullptr, 0, idx, off, n, q, hit);
}

void primScanTest()
{
  const int fails_before = g_fails;
  std::cout << "*************** test ntmesh_primscan *********************\n";

  nvec3<float> pts[8] = {
    nvec3<float>(0, 0, 0), nvec3<float>(1, 0, 0), nvec3<float>(0, 1, 0),
    nvec3<float>(10, 0, 0), nvec3<float>(11, 0, 0), nvec3<float>(10, 1, 0),
    nvec3<float>(1, 1, 0), nvec3<float>(3, 0, 0)
  };
  int one[] = {0, 1, 2};
  int two[] = {0, 1, 2, 3, 4, 5};
  int twin[] = {0, 1, 2, 0, 1, 2};
  ntm::PrimHit<float> hit;
  hit.dist2 = 5.f;
  hit.kind = 8;
  hit.prim = 7;
  expect_true(tri(pts, one, 0, nvec3<float>(0, 0, 0), hit) == -1
              && hit.kind == 8 && hit.prim == 7 && hit.dist2 == 5.f, "an empty triangle list leaves the hit");

  expect_true(tri(pts, one, 1, pts[0], hit) == ntm::ScanTri && hit.prim == 0 && hit.dist2 == 0.f,
              "a query on a vertex has distance zero");
  expect_vec(hit.p, pts[0], 1e-6f, "the hit point is that vertex");
  expect_vec(hit.w, nvec3<float>(1, 0, 0), 1e-6f, "a vertex hit is weight one");

  expect_true(tri(pts, one, 1, nvec3<float>(0.5f, -1, 0), hit) == ntm::ScanTri, "a query beside an edge hits that edge");
  expect_vec(hit.p, nvec3<float>(0.5f, 0, 0), 1e-5f, "the edge point is the projection");
  expect_vec(hit.w, nvec3<float>(0.5f, 0.5f, 0), 1e-5f, "an edge leaves the opposite weight at zero");
  expect_true(hit.dist2 > 0.99f && hit.dist2 < 1.01f, "the edge distance is one");

  expect_true(tri(pts, one, 1, nvec3<float>(0.2f, 0.2f, 1), hit) == ntm::ScanTri, "a query over the face hits the interior");
  expect_vec(hit.p, nvec3<float>(0.2f, 0.2f, 0), 1e-5f, "the interior point is the projection");
  expect_vec(hit.w, nvec3<float>(0.6f, 0.2f, 0.2f), 1e-5f, "interior weights match the projection");

  expect_true(tri(pts, one, 1, nvec3<float>(2, 2, 0), hit) == ntm::ScanTri, "a query past a corner lands on the near edge");
  expect_vec(hit.p, nvec3<float>(0.5f, 0.5f, 0), 1e-5f, "the near edge point is the midpoint");
  expect_vec(hit.w, nvec3<float>(0, 0.5f, 0.5f), 1e-5f, "the far vertex weight is zero");

  expect_true(tri(pts, two, 2, nvec3<float>(10, 0, 1), hit) == ntm::ScanTri && hit.prim == 1,
              "the second triangle wins when it is closer");
  expect_vec(hit.p, nvec3<float>(10, 0, 0), 1e-5f, "the second hit is its nearest vertex");

  expect_true(tri(pts, twin, 2, nvec3<float>(0, 0, 1), hit) == ntm::ScanTri && hit.prim == 0,
              "a tie keeps the lower triangle");

  nvec3<float> flat[3] = {nvec3<float>(3, 4, 5), nvec3<float>(3, 4, 5), nvec3<float>(3, 4, 5)};
  int collapsed[] = {0, 1, 2};
  expect_true(tri(flat, collapsed, 1, nvec3<float>(3, 4, 6), hit) == ntm::ScanTri,
              "a zero-area triangle still returns a hit");
  expect_vec(hit.p, flat[0], 1e-6f, "a collapsed triangle returns its vertex");
  expect_true(hit.dist2 > 0.99f && hit.dist2 < 1.01f, "a collapsed triangle keeps a finite distance");

  int quadIdx[] = {0, 1, 6, 2};
  hit.kind = 4;
  expect_true(quad(pts, quadIdx, 0, nvec3<float>(0, 0, 0), hit) == -1 && hit.kind == 4, "an empty quad list leaves the hit");
  expect_true(quad(pts, quadIdx, 1, nvec3<float>(0.9f, 0.1f, 2), hit) == ntm::ScanQuad && hit.sub == 0,
              "a quad query hits the near half");
  expect_vec(hit.p, nvec3<float>(0.9f, 0.1f, 0), 1e-5f, "the quad point is the projection");
  expect_vec(hit.w, nvec3<float>(0.1f, 0.8f, 0.1f), 1e-4f, "the quad weights are on that half");
  expect_true(hit.dist2 > 3.99f && hit.dist2 < 4.01f, "the quad is two units away");

  expect_true(quad(pts, quadIdx, 1, nvec3<float>(0.5f, 0.5f, 1), hit) == ntm::ScanQuad && hit.sub == 0,
              "a diagonal tie keeps half zero");
  expect_vec(hit.w, nvec3<float>(0.5f, 0, 0.5f), 1e-4f, "a diagonal hit has a zero middle weight");

  int ngonIdx[] = {0, 1, 2};
  int ngonOff[] = {0, 3};
  hit.kind = 9;
  expect_true(ngon(pts, ngonIdx, ngonOff, 0, nvec3<float>(0, 0, 1), hit) == -1 && hit.kind == 9,
              "an empty n-gon list leaves the hit");
  expect_true(ngon(pts, ngonIdx, ngonOff, 1, nvec3<float>(0, 0, 1), hit) == ntm::ScanNgon && hit.sub == 0,
              "a three-corner fan matches the triangle");
  expect_vec(hit.p, nvec3<float>(0, 0, 0), 1e-5f, "the fan hit is the vertex");
  expect_vec(hit.w, nvec3<float>(1, 0, 0), 1e-5f, "the fan weight is on corner zero");

  int lineIdx[] = {0, 1, 7};
  int lineOff[] = {0, 3};
  hit.kind = 3;
  hit.prim = 3;
  expect_true(line(pts, lineIdx, lineOff, 0, nvec3<float>(0.5f, -1, 0), hit) == -1 && hit.prim == 3,
              "an empty line list leaves the hit");
  expect_true(line(pts, lineIdx, lineOff, 1, nvec3<float>(0.5f, -1, 0), hit) == ntm::ScanLine && hit.sub == 0,
              "a polyline query hits the near segment");
  expect_vec(hit.p, nvec3<float>(0.5f, 0, 0), 1e-5f, "the segment point is the projection");
  expect_vec(hit.w, nvec3<float>(0.5f, 0.5f, 0), 1e-5f, "a segment stores two weights and a zero");

  expect_true(line(pts, lineIdx, lineOff, 1, nvec3<float>(3, 1, 0), hit) == ntm::ScanLine && hit.sub == 1,
              "the second segment wins");
  expect_vec(hit.p, nvec3<float>(3, 0, 0), 1e-5f, "the second segment hit is its end");
  expect_vec(hit.w, nvec3<float>(0, 1, 0), 1e-5f, "the end of a segment is weight on B");

  hit.kind = 8;
  hit.prim = 8;
  hit.dist2 = 5.f;
  expect_true(ntm::primScan(pts, one, 0, quadIdx, 0, ngonIdx, ngonOff, 0, lineIdx, lineOff, 0,
                            nvec3<float>(0, 0, 1), hit) == -1
              && hit.kind == 8 && hit.dist2 == 5.f, "an empty mesh leaves the hit");

  expect_true(ntm::primScan(pts, one, 1, nullptr, 0, nullptr, nullptr, 0, lineIdx, lineOff, 1,
                            nvec3<float>(3, 1, 0), hit) == ntm::ScanLine,
              "a closer line beats the triangle");
  expect_true(hit.prim == 0 && hit.sub == 1 && hit.dist2 > 0.99f && hit.dist2 < 1.01f,
              "the combined line is one unit away");
  expect_vec(hit.p, nvec3<float>(3, 0, 0), 1e-5f, "the line hit keeps the segment point");
  expect_vec(hit.w, nvec3<float>(0, 1, 0), 1e-5f, "the winning line keeps its weights");

  expect_true(ntm::primScan(pts, one, 1, nullptr, 0, nullptr, nullptr, 0, lineIdx, lineOff, 1,
                            nvec3<float>(0, 0, 1), hit) == ntm::ScanTri && hit.prim == 0,
              "an equal distance keeps the triangle");
  expect_vec(hit.w, nvec3<float>(1, 0, 0), 1e-5f, "the triangle hit holds its weights");

  const int nfail = g_fails - fails_before;
  if (nfail == 0)
    std::cout << "ntmesh_primscan: ok\n";
  else
    std::cout << "ntmesh_primscan: " << nfail << " failed\n";
}
