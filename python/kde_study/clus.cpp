// Standalone 1-D (optionally 2-D) track-time clustering for the offline workbench.
// Build: clang++ -O3 -std=c++17 -shared -fPIC -o libclus.dylib clus.cpp
// All algorithms work per event on tracks [off[e], off[e+1]) and write a cluster label per track.
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstdint>

namespace {
struct C { double t, s, z, sz; std::vector<int> m; bool merged = false; };

inline double dist(const C& a, const C& b, bool useZ, double floor2) {
  double dt = (a.t - b.t) / std::sqrt(a.s * a.s + b.s * b.s + floor2);
  double d2 = dt * dt;
  if (useZ) { double dz = (a.z - b.z) / std::sqrt(a.sz * a.sz + b.sz * b.sz); d2 += dz * dz; }
  return std::sqrt(d2);
}
inline C merge(const C& a, const C& b) {
  C r;
  double wa = 1 / (a.s * a.s), wb = 1 / (b.s * b.s);
  r.t = (a.t * wa + b.t * wb) / (wa + wb); r.s = 1 / std::sqrt(wa + wb);
  double za = 1 / (a.sz * a.sz), zb = 1 / (b.sz * b.sz);
  r.z = (a.z * za + b.z * zb) / (za + zb); r.sz = 1 / std::sqrt(za + zb);
  r.m = a.m; r.m.insert(r.m.end(), b.m.begin(), b.m.end());
  return r;
}
}  // namespace

extern "C" {
// method: 0 iterative, 1 simultaneous, 2 cone, 3 mean-shift modes, 4 iterative with fixed-window (cut in ps around weighted centroid)
// seedw: seed-ordering weight (pT for the production algorithm). kw: kernel / centroid weight for methods 3-4.
// par: method-specific (3: kernel width ps; 4: unused). floor: extra ps added in quadrature to the distance denominator.
void cluster(int nev, const int64_t* off, const double* t, const double* s, const double* z, const double* sz,
             const double* seedw, const double* kw, int method, double cut, int useZ, double floor, double par,
             int32_t* label) {
  const double floor2 = floor * floor;
  for (int e = 0; e < nev; ++e) {
    const int64_t b = off[e]; const int n = (int)(off[e + 1] - b);
    if (n == 0) continue;
    std::vector<C> col(n);
    for (int i = 0; i < n; ++i) { col[i].t = t[b + i]; col[i].s = s[b + i]; col[i].z = z[b + i]; col[i].sz = sz[b + i]; col[i].m = {i}; }
    std::vector<C> res;
    if (method == 0) {
      std::vector<char> used(n, 0);
      while (true) {
        int seed = -1; double mx = -1;
        for (int i = 0; i < n; ++i) if (!used[i] && seedw[b + i] > mx) { mx = seedw[b + i]; seed = i; }
        if (seed < 0) break;
        used[seed] = 1; C cur = col[seed];
        while (true) {
          int best = -1; double md = cut;
          for (int i = 0; i < n; ++i) { if (used[i]) continue; double d = dist(cur, col[i], useZ, floor2); if (d < md) { md = d; best = i; } }
          if (best < 0) break;
          used[best] = 1; cur = merge(cur, col[best]);
        }
        res.push_back(std::move(cur));
      }
    } else if (method == 1) {
      res = col;
      while (res.size() > 1) {
        int i0 = 0, j0 = 0; double d0 = 1e30;
        for (size_t i = 0; i < res.size(); ++i) for (size_t j = i + 1; j < res.size(); ++j) {
          double d = dist(res[i], res[j], useZ, floor2); if (d < d0) { d0 = d; i0 = (int)i; j0 = (int)j; }
        }
        if (d0 >= cut) break;
        C nc = merge(res[i0], res[j0]);
        res.erase(res.begin() + j0); res.erase(res.begin() + i0); res.push_back(std::move(nc));
      }
    } else if (method == 2) {
      std::vector<char> used(n, 0);
      while (true) {
        int seed = -1; double mx = -1;
        for (int i = 0; i < n; ++i) if (!used[i] && seedw[b + i] > mx) { mx = seedw[b + i]; seed = i; }
        if (seed < 0) break;
        used[seed] = 1; C cur = col[seed];
        for (int i = 0; i < n; ++i) { if (used[i]) continue; if (dist(col[seed], col[i], useZ, floor2) < cut) { used[i] = 1; cur = merge(cur, col[i]); } }
        res.push_back(std::move(cur));
      }
    } else if (method == 3) {
      // every track ascends the kernel density (weights kw, Gaussian width par); modes within `cut` ps are one cluster
      std::vector<double> mode(n);
      for (int i = 0; i < n; ++i) {
        double m = t[b + i];
        for (int it = 0; it < 30; ++it) {
          double num = 0, den = 0;
          for (int j = 0; j < n; ++j) { double u = (t[b + j] - m) / par; double g = kw[b + j] * std::exp(-0.5 * u * u); num += g * t[b + j]; den += g; }
          double m2 = den > 0 ? num / den : m; if (std::abs(m2 - m) < 0.01) { m = m2; break; } m = m2;
        }
        mode[i] = m;
      }
      std::vector<int> ord(n); for (int i = 0; i < n; ++i) ord[i] = i;
      std::sort(ord.begin(), ord.end(), [&](int a, int c) { return mode[a] < mode[c]; });
      for (int k = 0; k < n; ++k) {
        int i = ord[k];
        if (k == 0 || mode[i] - mode[ord[k - 1]] > cut) res.push_back(col[i]);
        else res.back() = merge(res.back(), col[i]);
      }
    } else if (method == 4) {
      // seed by seedw; absorb the nearest track while |t - centroid| < cut ps; centroid = kw-weighted mean
      std::vector<char> used(n, 0);
      while (true) {
        int seed = -1; double mx = -1;
        for (int i = 0; i < n; ++i) if (!used[i] && seedw[b + i] > mx) { mx = seedw[b + i]; seed = i; }
        if (seed < 0) break;
        used[seed] = 1; C cur = col[seed]; double sw = kw[b + seed], swt = kw[b + seed] * t[b + seed];
        while (true) {
          int best = -1; double md = cut; const double c = swt / sw;
          for (int i = 0; i < n; ++i) { if (used[i]) continue; double d = std::abs(t[b + i] - c); if (d < md) { md = d; best = i; } }
          if (best < 0) break;
          used[best] = 1; cur = merge(cur, col[best]); sw += kw[b + best]; swt += kw[b + best] * t[b + best];
        }
        res.push_back(std::move(cur));
      }
    }
    for (size_t k = 0; k < res.size(); ++k) for (int i : res[k].m) label[b + i] = (int32_t)k;
  }
}
}
