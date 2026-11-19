// Stress variant of nearpc_scene.C.
// A 40x40x40 bunny grid is about 163 million points. The .geo writer builds the
// whole document in memory before it writes a byte, and that document does not
// fit. This run is a 16x16x16 grid: 4096 bunnies, spacing 2, random rotations.
// Sample points fill the cube [-16, 16]^3, which is 32 units on a side and
// centered on the grid. The original 4x4 scene is left in nearpc_scene.C.

#include <chrono>
#include <cstdio>
#include <filesystem>
#include "ntmesh.h"
#include "ntmesh_file.h"
#include "ntmesh_nearpc.h"

using ntm::ntmesh;
using ntv::nmat4;
using ntv::nvec3;

static const int kGrid = 16;
static const int kQueries = 100000;
static unsigned g_seed = 1;

static float rnd() {
    g_seed = g_seed * 1664525u + 1013904223u;
    return (g_seed >> 8) * (1.f / 16777216.f);
}

static double seconds(std::chrono::steady_clock::time_point a, std::chrono::steady_clock::time_point b) {
    return std::chrono::duration<double>(b - a).count();
}

int main() {
    auto tRun = std::chrono::steady_clock::now();
    ntmesh<float> bunny;
    auto tRead0 = std::chrono::steady_clock::now();
    if (!ntm::readFile(bunny, LIGHT_BUNNY)) {
        std::fprintf(stderr, "could not read %s\n", LIGHT_BUNNY);
        return 1;
    }
    auto tRead1 = std::chrono::steady_clock::now();
    nvec3<float> mn, mx;
    bunny.bounds(mn, mx);
    nvec3<float> center = (mn + mx) * 0.5f;

    const int copies = kGrid * kGrid * kGrid;
    ntmesh<float> scene;
    scene.reserve(bunny.pointCount() * copies + kQueries,
                  bunny.vertexCount() * copies + kQueries * 2,
                  bunny.primCount() * copies + kQueries);

    float origin = -(float)(kGrid - 1);
    for (int iz = 0; iz < kGrid; ++iz) {
        for (int iy = 0; iy < kGrid; ++iy) {
            for (int ix = 0; ix < kGrid; ++ix) {
                nvec3<float> axis(rnd() * 2.f - 1.f, rnd() * 2.f - 1.f, rnd() * 2.f - 1.f);
                if (axis.length2() < 1e-8f) axis = nvec3<float>(0, 1, 0);
                else axis.normalize();
                nmat4<float> m;
                m.setRotation(axis, rnd() * 6.2831853f);
                nvec3<float> place(ix * 2.f + origin, iy * 2.f + origin, iz * 2.f + origin);
                m.translate(place - (m * center));
                ntmesh<float> inst;
                ntm::copyMesh(inst, bunny);
                inst.transform(m);
                scene.merge(inst);
                inst.clear();
            }
        }
        std::printf("placed z %d / %d, points %d\n", iz + 1, kGrid, scene.pointCount());
        std::fflush(stdout);
    }
    bunny.clear();

    ntm::NearPC<float> cloud;
    auto t0 = std::chrono::steady_clock::now();
    cloud.build(scene.points.data, scene.pointCount());
    auto t1 = std::chrono::steady_clock::now();

    nvec3<float> queries[kQueries];
    int nearest[kQueries];
    float half = (float)kGrid;
    for (int i = 0; i < kQueries; ++i)
        queries[i] = nvec3<float>(rnd() * 2.f * half - half, rnd() * 2.f * half - half, rnd() * 2.f * half - half);
    auto t2 = std::chrono::steady_clock::now();
    for (int i = 0; i < kQueries; ++i) {
        float d2 = 0.f;
        nearest[i] = cloud.nearest(queries[i], d2);
    }
    auto t3 = std::chrono::steady_clock::now();
    cloud.clear();

    for (int i = 0; i < kQueries; ++i) {
        int qid = scene.addPoint(queries[i]);
        int line[] = {qid, nearest[i]};
        scene.addLine(line, 2);
    }

    std::printf("points %d\n", scene.pointCount());
    std::printf("bvh_build_s %.3f\n", seconds(t0, t1));
    std::printf("queries %d query_s %.3f\n", kQueries, seconds(t2, t3));
    std::fflush(stdout);

    std::filesystem::create_directories("build");
    auto tWrite0 = std::chrono::steady_clock::now();
    if (!ntm::writeFile(scene, "build/BIG_test_nearpoints.geo")) {
        std::fprintf(stderr, "could not write build/BIG_test_nearpoints.geo\n");
        scene.clear();
        return 1;
    }
    auto tWrite1 = std::chrono::steady_clock::now();
    std::printf("read_s %.3f\n", seconds(tRead0, tRead1));
    std::printf("write_s %.3f\n", seconds(tWrite0, tWrite1));
    std::printf("wrote build/BIG_test_nearpoints.geo lines %d\n", scene.lineCount());
    scene.clear();
    std::printf("total_s %.3f\n", seconds(tRun, std::chrono::steady_clock::now()));
    return 0;
}
