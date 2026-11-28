// One closest-point scene over quads.
// s / small   4x4 twisted cubes on z = 0, 1000 queries in [-4.5, 4.5]^3, fan
// m / medium  4x4x4 twisted cubes, 10000 queries in [-4.5, 4.5]^3, fan, then the box edges
// l / large   16x16x16 twisted cubes, 100000 queries in [-16, 16]^3, fan
// b / bench   small only, one thread, methods fan, newton, algebraic, planar. No .geo.
// Centers are 2 apart. Each copy gets a random rotation.
// --threads N is the query thread count for s, m, l (1 serial, -1 every core).
// Bench ignores it. No arguments runs small on 1 thread.
// s, m, and l write samples/quadscan/<name>.geo

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>
#include <vector>
#include "ntmesh.h"
#include "ntmesh_file.h"
#include "ntmesh_quadscan.h"

using ntm::ntmesh;
using ntv::nmat4;
using ntv::nvec3;

struct Preset {
    const char* tag;
    const char* name;
    const char* file;
    int n;
    int nz;
    bool flat;
    int queries;
    float half;
    bool box;
};

static const Preset kSmall = {"s", "small", "small.geo", 4, 1, true, 1000, 4.5f, false};
static const Preset kMedium = {"m", "medium", "medium.geo", 4, 4, false, 10000, 4.5f, true};
static const Preset kLarge = {"l", "large", "large.geo", 16, 16, false, 100000, 16.f, false};

static const char* kMethodName[4] = {"fan", "newton", "algebraic", "planar"};

enum Mode { ModeSmall, ModeMedium, ModeLarge, ModeBench };

struct Args {
    Mode mode = ModeSmall;
    int threads = 1;
    bool threadsGiven = false;
};

static unsigned g_seed = 1;

static float rnd() {
    g_seed = g_seed * 1664525u + 1013904223u;
    return (g_seed >> 8) * (1.f / 16777216.f);
}

static double seconds(std::chrono::steady_clock::time_point a, std::chrono::steady_clock::time_point b) {
    return std::chrono::duration<double>(b - a).count();
}

static void usage() {
    std::fprintf(stderr,
        "usage: quadscan_scene [s|m|l|b] [--threads N]\n"
        "  s, small    4x4, 1000 queries, fan\n"
        "  m, medium   4x4x4, 10000 queries, fan\n"
        "  l, large    16x16x16, 100000 queries, fan\n"
        "  b, bench    time fan, newton, algebraic, and planar on small\n"
        "  --threads N query threads for s, m, l (1 serial, -1 every core)\n"
        "bench ignores --threads and stays on one thread.\n"
        "no arguments runs s on 1 thread.\n");
}

static int hardwareThreads() {
    int n = (int)std::thread::hardware_concurrency();
    if (n < 1) n = 1;
    return n;
}

static int resolveThreads(int requested) {
    if (requested == -1) return hardwareThreads();
    return requested;
}

static bool parseInt(const char* s, long& out) {
    if (!s || !s[0]) return false;
    char* end = nullptr;
    out = std::strtol(s, &end, 10);
    return end && *end == 0;
}

static int modeOf(const char* s) {
    if (std::strcmp(s, "s") == 0 || std::strcmp(s, "small") == 0) return ModeSmall;
    if (std::strcmp(s, "m") == 0 || std::strcmp(s, "medium") == 0) return ModeMedium;
    if (std::strcmp(s, "l") == 0 || std::strcmp(s, "large") == 0) return ModeLarge;
    if (std::strcmp(s, "b") == 0 || std::strcmp(s, "bench") == 0) return ModeBench;
    return -1;
}

static bool parseArgs(int argc, char** argv, Args& args) {
    bool modeSet = false;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--threads") == 0) {
            if (i + 1 >= argc) {
                usage();
                return false;
            }
            long v = 0;
            if (!parseInt(argv[++i], v) || v == 0 || v < -1 || v > 100000) {
                usage();
                return false;
            }
            args.threads = (int)v;
            args.threadsGiven = true;
            continue;
        }
        int mode = modeOf(argv[i]);
        if (mode < 0 || modeSet) {
            usage();
            return false;
        }
        args.mode = (Mode)mode;
        modeSet = true;
    }
    return true;
}

static void placeCopy(ntmesh<float>& scene, const ntmesh<float>& cube, const nvec3<float>& pivot, nvec3<float> at) {
    nvec3<float> axis(rnd() * 2.f - 1.f, rnd() * 2.f - 1.f, rnd() * 2.f - 1.f);
    if (axis.length2() < 1e-8f) axis = nvec3<float>(0, 1, 0);
    else axis.normalize();
    nmat4<float> m;
    m.setRotation(axis, rnd() * 6.2831853f);
    m.translate(at - (m * pivot));
    ntmesh<float> inst;
    ntm::copyMesh(inst, cube);
    inst.transform(m);
    scene.merge(inst);
    inst.clear();
}

static void buildGrid(ntmesh<float>& scene, const ntmesh<float>& cube, const nvec3<float>& pivot, const Preset& preset) {
    g_seed = 1;
    int copies = preset.n * preset.n * preset.nz;
    int boxPts = preset.box ? 8 : 0;
    int boxLines = preset.box ? 12 : 0;
    scene.reserve(cube.pointCount() * copies + preset.queries * 2 + boxPts,
                  cube.vertexCount() * copies + preset.queries * 2 + boxLines * 2,
                  cube.primCount() * copies + preset.queries + boxLines);
    float origin = -(float)(preset.n - 1);
    for (int iz = 0; iz < preset.nz; ++iz) {
        for (int iy = 0; iy < preset.n; ++iy) {
            for (int ix = 0; ix < preset.n; ++ix) {
                float z = preset.flat ? 0.f : iz * 2.f + origin;
                placeCopy(scene, cube, pivot, nvec3<float>(ix * 2.f + origin, iy * 2.f + origin, z));
            }
        }
        if (preset.n >= 16) {
            std::printf("%s z %d / %d, points %d\n", preset.name, iz + 1, preset.nz, scene.pointCount());
            std::fflush(stdout);
        }
    }
}

static void makeQueries(std::vector<nvec3<float>>& queries, const Preset& preset) {
    queries.resize(preset.queries);
    float half = preset.half;
    for (int i = 0; i < preset.queries; ++i)
        queries[i] = nvec3<float>(rnd() * 2.f * half - half, rnd() * 2.f * half - half, rnd() * 2.f * half - half);
}

static void querySlice(const nvec3<float>* pts, const int* quads, int nquads, int method,
                       const nvec3<float>* q, ntm::PrimHit<float>* hits, int begin, int end) {
    for (int i = begin; i < end; ++i)
        ntm::quadScan(pts, quads, nquads, method, q[i], hits[i]);
}

static void queryAll(const ntmesh<float>& scene, const std::vector<nvec3<float>>& queries,
                     std::vector<ntm::PrimHit<float>>& hits, int method, int threads) {
    int n = (int)queries.size();
    if (threads > n) threads = n;
    if (threads < 1) threads = 1;
    hits.assign(n, ntm::PrimHit<float>());
    const nvec3<float>* pts = scene.points.data;
    const int* quads = scene.quadIndices.data;
    int nquads = scene.quadCount();
    if (threads == 1) {
        querySlice(pts, quads, nquads, method, queries.data(), hits.data(), 0, n);
        return;
    }
    int chunk = (n + threads - 1) / threads;
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; ++t) {
        int begin = t * chunk;
        if (begin >= n) break;
        int end = begin + chunk;
        if (end > n) end = n;
        pool.emplace_back(querySlice, pts, quads, nquads, method, queries.data(), hits.data(), begin, end);
    }
    for (std::thread& th : pool) th.join();
}

static void addBox(ntmesh<float>& scene, float half) {
    int corner[8];
    for (int i = 0; i < 8; ++i) {
        float x = (i & 1) ? half : -half;
        float y = (i & 2) ? half : -half;
        float z = (i & 4) ? half : -half;
        corner[i] = scene.addPoint(nvec3<float>(x, y, z));
    }
    const int edges[12][2] = {
        {0, 1}, {2, 3}, {4, 5}, {6, 7},
        {0, 2}, {1, 3}, {4, 6}, {5, 7},
        {0, 4}, {1, 5}, {2, 6}, {3, 7},
    };
    for (int e = 0; e < 12; ++e) {
        int line[] = {corner[edges[e][0]], corner[edges[e][1]]};
        scene.addLine(line, 2);
    }
}

static void gridText(const Preset& preset, char* buf, int n) {
    if (preset.flat) std::snprintf(buf, n, "%dx%d", preset.n, preset.n);
    else std::snprintf(buf, n, "%dx%dx%d", preset.n, preset.n, preset.nz);
}

static int writeScene(ntmesh<float>& scene, const Preset& preset,
                      const std::vector<nvec3<float>>& queries, const std::vector<ntm::PrimHit<float>>& hits) {
    for (int i = 0; i < preset.queries; ++i) {
        int qid = scene.addPoint(queries[i]);
        int hid = scene.addPoint(hits[i].p);
        int line[] = {qid, hid};
        scene.addLine(line, 2);
    }
    if (preset.box) addBox(scene, preset.half);
    std::string path = std::string(OUT_DIR) + "/" + preset.file;
    if (!ntm::writeFile(scene, path.c_str())) {
        std::fprintf(stderr, "could not write %s\n", path.c_str());
        return 1;
    }
    std::printf("points %d lines %d -> %s\n", scene.pointCount(), scene.lineCount(), path.c_str());
    return 0;
}

static int runScene(const ntmesh<float>& cube, const nvec3<float>& pivot, const Preset& preset, int threads) {
    char grid[64];
    gridText(preset, grid, (int)sizeof(grid));
    std::printf("%s: %s, %d queries, %d thread%s, fan\n",
                preset.name, grid, preset.queries, threads, threads == 1 ? "" : "s");
    std::fflush(stdout);

    ntmesh<float> scene;
    buildGrid(scene, cube, pivot, preset);
    std::printf("prims %d  tris %d  quads %d  ngons %d  lines %d\n",
                scene.primCount(), scene.triangleCount(), scene.quadCount(),
                scene.ngonCount(), scene.lineCount());
    std::fflush(stdout);
    std::vector<nvec3<float>> queries;
    makeQueries(queries, preset);
    std::vector<ntm::PrimHit<float>> hits;
    auto t0 = std::chrono::steady_clock::now();
    queryAll(scene, queries, hits, 0, threads);
    auto t1 = std::chrono::steady_clock::now();
    std::printf("query %.3fs\n", seconds(t0, t1));
    int rc = writeScene(scene, preset, queries, hits);
    scene.clear();
    return rc;
}

static int runBench(const ntmesh<float>& cube, const nvec3<float>& pivot, bool threadsGiven) {
    if (threadsGiven) std::printf("bench ignores --threads\n");
    const Preset& preset = kSmall;
    ntmesh<float> scene;
    buildGrid(scene, cube, pivot, preset);
    std::vector<nvec3<float>> queries;
    makeQueries(queries, preset);
    std::printf("s: %d quads, %d queries, 1 thread\n", scene.quadCount(), preset.queries);
    std::fflush(stdout);

    std::vector<ntm::PrimHit<float>> fanHits;
    std::vector<ntm::PrimHit<float>> hits;
    double total = 0;
    double fanSec = 0;
    int misses = 0;
    for (int method = 0; method < 4; ++method) {
        auto t0 = std::chrono::steady_clock::now();
        queryAll(scene, queries, hits, method, 1);
        auto t1 = std::chrono::steady_clock::now();
        double sec = seconds(t0, t1);
        total += sec;
        double us = sec * 1e6 / (double)preset.queries;
        if (method == 0) {
            fanHits.swap(hits);
            fanSec = sec;
            std::printf("s  %-10s  query %.3fs  %.2f us/q\n", kMethodName[method], sec, us);
        } else {
            float maxD2 = 0.f;
            int primDiffer = 0;
            for (int i = 0; i < preset.queries; ++i) {
                if (hits[i].kind < 0) {
                    ++misses;
                    std::fprintf(stderr, "%s query %d missed\n", kMethodName[method], i);
                    continue;
                }
                float gap = hits[i].dist2 - fanHits[i].dist2;
                if (gap < 0.f) gap = -gap;
                if (gap > maxD2) maxD2 = gap;
                if (hits[i].prim != fanHits[i].prim) ++primDiffer;
            }
            double ratio = sec > 0.0 && fanSec > 0.0 ? sec / fanSec : 0.0;
            std::printf("s  %-10s  query %.3fs  %.2f us/q  %.2fx fan  max d2 %.3g  prim differ %d\n",
                        kMethodName[method], sec, us, ratio, maxD2, primDiffer);
        }
        std::fflush(stdout);
    }
    std::printf("4 runs        query %.3fs\n", total);
    scene.clear();
    return misses ? 1 : 0;
}

int main(int argc, char** argv) {
    Args args;
    if (!parseArgs(argc, argv, args)) return 1;

    ntmesh<float> cube;
    if (!ntm::readFile(cube, TWISTED_CUBE)) {
        std::fprintf(stderr, "could not read %s\n", TWISTED_CUBE);
        return 1;
    }
    nvec3<float> mn, mx;
    cube.bounds(mn, mx);
    nvec3<float> pivot = (mn + mx) * 0.5f;

    int rc = 0;
    if (args.mode == ModeBench) rc = runBench(cube, pivot, args.threadsGiven);
    else if (args.mode == ModeMedium) rc = runScene(cube, pivot, kMedium, resolveThreads(args.threads));
    else if (args.mode == ModeLarge) rc = runScene(cube, pivot, kLarge, resolveThreads(args.threads));
    else rc = runScene(cube, pivot, kSmall, resolveThreads(args.threads));
    cube.clear();
    return rc;
}
