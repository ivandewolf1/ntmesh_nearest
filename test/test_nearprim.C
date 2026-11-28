#include <iostream>
#include "tests.h"
#include "ntmesh_nearprim.h"

using ntv::nvec3;

static int sameHit(const ntm::PrimHit<float>& a, const ntm::PrimHit<float>& b) {
  return a.kind == b.kind && a.prim == b.prim && a.sub == b.sub
      && a.p.almostEqual(b.p, 1e-5f) && a.w.almostEqual(b.w, 1e-4f)
      && (a.dist2 > b.dist2 ? a.dist2 - b.dist2 : b.dist2 - a.dist2) < 1e-4f;
}

static void checkQuery(ntm::NearPrim<float>& tree, const nvec3<float>* pts,
                       const int* tris, int ntris, const int* quads, int nquads,
                       const int* ngonIdx, const int* ngonOff, int nngons,
                       const int* lineIdx, const int* lineOff, int nlines,
                       nvec3<float> q, const char* what) {
  ntm::PrimHit<float> scan;
  ntm::PrimHit<float> treeHit;
  scan.kind = 3;
  scan.prim = 9;
  scan.dist2 = 4.f;
  treeHit = scan;
  int ks = ntm::primScan(pts, tris, ntris, quads, nquads, ngonIdx, ngonOff, nngons,
                         lineIdx, lineOff, nlines, q, scan);
  int kt = tree.query(q, treeHit);
  expect_true(ks == kt && sameHit(scan, treeHit), what);
}

void nearPrimTest()
{
  const int fails_before = g_fails;
  std::cout << "*************** test ntmesh_nearprim *********************\n";

  nvec3<float> pts[8] = {
    nvec3<float>(0, 0, 0), nvec3<float>(1, 0, 0), nvec3<float>(0, 1, 0),
    nvec3<float>(10, 0, 0), nvec3<float>(11, 0, 0), nvec3<float>(10, 1, 0),
    nvec3<float>(1, 1, 0), nvec3<float>(3, 0, 0)
  };
  int tris[] = {0, 1, 2, 0, 1, 2};
  int quads[] = {0, 1, 6, 2};
  int ngonIdx[] = {3, 4, 5};
  int ngonOff[] = {0, 3};
  int lineIdx[] = {0, 1, 7};
  int lineOff[] = {0, 3};

  ntm::NearPrim<float> tree;
  ntm::PrimHit<float> keep;
  keep.kind = 8;
  keep.prim = 7;
  keep.dist2 = 5.f;
  expect_true(tree.query(nvec3<float>(0, 0, 1), keep) == -1
              && keep.kind == 8 && keep.prim == 7 && keep.dist2 == 5.f,
              "an empty tree leaves the hit");

  tree.build(pts, tris, 2, quads, 1, ngonIdx, ngonOff, 1, lineIdx, lineOff, 1);
  expect_true(tree.nleaves == 2 + 2 + 1 + 2 && tree.count == 2 * tree.nleaves - 1,
              "fans expand and the node count is 2n-1");

  nvec3<float> saved = pts[0];
  ntm::PrimHit<float> before;
  nvec3<float> qCopy(0, 0, 1);
  tree.query(qCopy, before);
  pts[0] = nvec3<float>(50, 50, 50);
  ntm::PrimHit<float> after;
  tree.query(qCopy, after);
  expect_true(sameHit(before, after), "a copied corner survives a change to the mesh");
  pts[0] = saved;

  nvec3<float> queries[8] = {
    nvec3<float>(0, 0, 0), nvec3<float>(0.2f, 0.2f, 1), nvec3<float>(0.5f, -1, 0),
    nvec3<float>(2, 2, 0), nvec3<float>(0.9f, 0.1f, 2), nvec3<float>(0.5f, 0.5f, 1),
    nvec3<float>(10, 0, 1), nvec3<float>(3, 1, 0)
  };
  const char* names[8] = {
    "vertex", "interior", "edge", "far corner", "quad half", "diagonal", "n-gon", "line end"
  };
  for (int i = 0; i < 8; ++i)
    checkQuery(tree, pts, tris, 2, quads, 1, ngonIdx, ngonOff, 1, lineIdx, lineOff, 1,
               queries[i], names[i]);

  ntm::PrimHit<float> batch[8];
  for (int i = 0; i < 8; ++i) {
    batch[i].kind = 4;
    batch[i].dist2 = 9.f;
  }
  tree.query(queries, 8, batch);
  for (int i = 0; i < 8; ++i) {
    ntm::PrimHit<float> one;
    tree.query(queries[i], one);
    expect_true(sameHit(one, batch[i]), "a short list matches the single query");
  }
  tree.clear();

  const int nq = 1200;
  nvec3<float>* manyPts = new nvec3<float>[64];
  int* manyTris = new int[96];
  unsigned seed = 1;
  for (int i = 0; i < 64; ++i) {
    seed = seed * 1664525u + 1013904223u;
    float x = (seed >> 8) * (1.f / 16777216.f);
    seed = seed * 1664525u + 1013904223u;
    float y = (seed >> 8) * (1.f / 16777216.f);
    seed = seed * 1664525u + 1013904223u;
    float z = (seed >> 8) * (1.f / 16777216.f);
    manyPts[i] = nvec3<float>(x * 8.f - 4.f, y * 8.f - 4.f, z * 2.f - 1.f);
  }
  for (int i = 0; i < 32; ++i) {
    manyTris[i * 3] = i % 64;
    manyTris[i * 3 + 1] = (i * 3 + 1) % 64;
    manyTris[i * 3 + 2] = (i * 5 + 2) % 64;
  }
  tree.build(manyPts, manyTris, 32, nullptr, 0, nullptr, nullptr, 0, nullptr, nullptr, 0);
  nvec3<float>* qs = new nvec3<float>[nq];
  ntm::PrimHit<float>* hits = new ntm::PrimHit<float>[nq];
  for (int i = 0; i < nq; ++i) {
    seed = seed * 1664525u + 1013904223u;
    float x = (seed >> 8) * (1.f / 16777216.f);
    seed = seed * 1664525u + 1013904223u;
    float y = (seed >> 8) * (1.f / 16777216.f);
    seed = seed * 1664525u + 1013904223u;
    float z = (seed >> 8) * (1.f / 16777216.f);
    qs[i] = nvec3<float>(x * 10.f - 5.f, y * 10.f - 5.f, z * 6.f - 3.f);
    hits[i].kind = -2;
  }
  tree.query(qs, nq, hits);
  int mismatch = 0;
  for (int i = 0; i < nq; ++i) {
    ntm::PrimHit<float> scan;
    int ks = ntm::primScan(manyPts, manyTris, 32, nullptr, 0, nullptr, nullptr, 0,
                           nullptr, nullptr, 0, qs[i], scan);
    if (ks != hits[i].kind || !sameHit(scan, hits[i])) ++mismatch;
  }
  expect_true(mismatch == 0, "a threaded Morton list matches the scan");
  tree.clear();
  delete[] manyPts;
  delete[] manyTris;
  delete[] qs;
  delete[] hits;

  const int nfail = g_fails - fails_before;
  if (nfail == 0)
    std::cout << "ntmesh_nearprim: ok\n";
  else
    std::cout << "ntmesh_nearprim: " << nfail << " failed\n";
}
