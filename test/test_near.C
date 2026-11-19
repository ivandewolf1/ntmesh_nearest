#include <iostream>
#include "tests.h"
#include "ntmesh_nearpc.h"

using ntv::nvec3;

void nearestTest()
{
  const int fails_before = g_fails;
  std::cout << "*************** test ntmesh_nearest *********************\n";

  {
    ntm::NearPC<float> cloud;
    float keep = 5.f;
    expect_true(cloud.nearest(nvec3<float>(0, 0, 0), keep) == -1 && keep == 5.f, "empty cloud leaves the distance");
    nvec3<float> one(2, 3, 4);
    cloud.build(&one, 1);
    float d2 = -1.f;
    expect_true(cloud.nearest(one, d2) == 0 && d2 == 0.f, "the only point is nearest");
    cloud.clear();

    nvec3<float> same[3] = {nvec3<float>(1, 1, 1), nvec3<float>(1, 1, 1), nvec3<float>(1, 1, 1)};
    cloud.build(same, 3);
    expect_true(cloud.count == 5, "three points make five nodes");
    d2 = -1.f;
    expect_true(cloud.nearest(nvec3<float>(1, 1, 1), d2) == 0 && d2 == 0.f, "a tie keeps the lower index");
    cloud.clear();

    const int n = 20000;
    nvec3<float>* pts = new nvec3<float>[n];
    unsigned seed = 1;
    for (int i = 0; i < n; ++i) {
      float p[3];
      for (int a = 0; a < 3; ++a) {
        seed = seed * 1664525u + 1013904223u;
        p[a] = (seed >> 8) * (1.f / 16777216.f) * 20.f - 10.f;
      }
      pts[i] = nvec3<float>(p[0], p[1], p[2]);
    }
    pts[100] = pts[50];
    cloud.build(pts, n);
    expect_true(cloud.count == 2 * n - 1, "node count is two per point minus one");
    nvec3<float> queries[4] = {pts[50], pts[15000], nvec3<float>(0, 0, 0), nvec3<float>(40, -30, 12)};
    for (int k = 0; k < 4; ++k) {
      int scan = -1;
      float scanD = 0.f;
      for (int i = 0; i < n; ++i) {
        float d = (pts[i] - queries[k]).length2();
        if (scan < 0 || d < scanD || (d == scanD && i < scan)) {
          scan = i;
          scanD = d;
        }
      }
      float gotD = -1.f;
      int got = cloud.nearest(queries[k], gotD);
      expect_true(got == scan && gotD == scanD, "tree matches a full scan");
    }
    cloud.clear();
    expect_true(cloud.nearest(queries[0], keep) == -1, "cleared cloud is empty");
    delete[] pts;
  }

  const int nfail = g_fails - fails_before;
  if (nfail == 0)
    std::cout << "ntmesh_nearest: ok\n";
  else
    std::cout << "ntmesh_nearest: " << nfail << " failed\n";
}
