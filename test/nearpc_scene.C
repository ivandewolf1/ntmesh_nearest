// Sixteen rotated copies of the bunny on a 4x4 grid, then 1000 nearest-point lines.
// Grid centers are 2 units apart: x and y at -3, -1, 1, 3, z = 0.
// Sample points are uniform in the cube [-4.5, 4.5]^3, which is 9 units on a side
// and centered on that grid.
// this is the first real code to exercise the ntmesh library, 
// and it was generated with AI from this prompt:
//
// light_bunny.geo is roughly the size of a unit cube. 
// Could you write me a test app that reads this bunny, 
// creates 16 randomly rotated bunnies in a 4  by 4 by 4 cube, 
// offset by 2 units between buunnies, 
// and merges them all  into a single ntmesh. 
// It should then generate 1000 random bpoints scattered in this 9 by 9 by 9 cube, 
// finds the closest point in the pointcloud,
// and appends a poline connecting the new point to the closest-in-poiintcloud. 
// Then dump this geometry out to build/test_nearpc.geo

#include <cstdio>
#include <filesystem>
#include "ntmesh.h"
#include "ntmesh_file.h"
#include "ntmesh_nearpc.h"

using ntm::ntmesh;
using ntv::nmat4;
using ntv::nvec3;

static unsigned g_seed = 1;

static float rnd() {
    g_seed = g_seed * 1664525u + 1013904223u;
    return (g_seed >> 8) * (1.f / 16777216.f);
}

int main() {
    ntmesh<float> bunny;
    if (!ntm::readFile(bunny, LIGHT_BUNNY)) {
        std::fprintf(stderr, "could not read %s\n", LIGHT_BUNNY);
        return 1;
    }
    nvec3<float> mn, mx;
    bunny.bounds(mn, mx);
    nvec3<float> center = (mn + mx) * 0.5f;

    ntmesh<float> scene;
    for (int iy = 0; iy < 4; ++iy) {
        for (int ix = 0; ix < 4; ++ix) {
            nvec3<float> axis(rnd() * 2.f - 1.f, rnd() * 2.f - 1.f, rnd() * 2.f - 1.f);
            if (axis.length2() < 1e-8f) axis = nvec3<float>(0, 1, 0);
            else axis.normalize();
            nmat4<float> m;
            m.setRotation(axis, rnd() * 6.2831853f);
            nvec3<float> place((ix - 1.5f) * 2.f, (iy - 1.5f) * 2.f, 0.f);
            m.translate(place - (m * center));
            ntmesh<float> inst;
            ntm::copyMesh(inst, bunny);
            inst.transform(m);
            scene.merge(inst);
            inst.clear();
        }
    }
    bunny.clear();

    ntm::NearPC<float> cloud;
    cloud.build(scene.points.data, scene.pointCount());

    for (int i = 0; i < 1000; ++i) {
        nvec3<float> q(rnd() * 9.f - 4.5f, rnd() * 9.f - 4.5f, rnd() * 9.f - 4.5f);
        float d2 = 0.f;
        int near = cloud.nearest(q, d2);
        int qid = scene.addPoint(q);
        int line[] = {qid, near};
        scene.addLine(line, 2);
    }
    cloud.clear();

    std::filesystem::create_directories("build");
    if (!ntm::writeFile(scene, "build/test_nearpc.geo")) {
        std::fprintf(stderr, "could not write build/test_nearpc.geo\n");
        scene.clear();
        return 1;
    }
    std::printf("points %d lines %d -> build/test_nearpc.geo\n", scene.pointCount(), scene.lineCount());
    scene.clear();
    return 0;
}
