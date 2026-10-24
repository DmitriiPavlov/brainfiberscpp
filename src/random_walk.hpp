#pragma once
#include "graph.hpp"

struct Lcg {
    uint64_t s;
    explicit Lcg(uint64_t seed = 1) : s(seed) {}
    uint32_t operator()() {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        return uint32_t(s >> 32);   
    }
};

inline Lcg rand_particle(1234); //chosen completely randomly

class Particle{
public:
    double s;
    int a;
    int b;

    Particle() = default;

    void step(Graph& g, double dist) {
        while (dist > 0) {
            const double len = (g[b] - g[a]).norm();   // cache per edge if this shows up in profiles
            const double remaining = len - s;

            if (dist < remaining) {
                s += dist;
                return;
            }

            dist -= remaining;
            a = b;
            choose_rand_nb(g);
            s = 0;
        }
    }

    void step(PlanarGraphView& g, double dist) {
        while (dist > 0) {
            const double len = (g[b] - g[a]).norm();   // cache per edge if this shows up in profiles
            const double remaining = len - s;

            if (dist < remaining) {
                s += dist;
                return;
            }

            dist -= std::max(remaining, 1e-12);
            a = b;
            choose_rand_nb(g);
            s = 0;
        }
    }

    void choose_rand_nb(Graph& g){
        const auto& nb = g.neighbors(a);
        b = nb[rand_particle() % nb.size()];
    }

    void choose_rand_nb(PlanarGraphView& g){
        const auto& nb = g.neighbors(a);
        b = nb[rand_particle() % nb.size()];
    }

    Vec2 position(const Graph& g) const {
        const double len = (g[b] - g[a]).norm();
        return g[a] + (g[b] - g[a]) * (len > 0 ? s / len : 0.0);
    }
    
    Vec2 position(const PlanarGraphView& g) const {
        const double len = (g[b] - g[a]).norm();
        return g[a] + (g[b] - g[a]) * (len > 0 ? s / len : 0.0);
    }

};

class ParticleAggregate {
public:
    ParticleAggregate(PlanarGraphView& g, int i, int n)
        : graph(g), total(n), concentration(g.h * g.w), view(concentration.data(), g.h, g.w) {
        point_mass(i, n);
    }
    ParticleAggregate(const ParticleAggregate&) = delete;

    //reading triggers a refresh if we mutated with step_all
    double operator[](int x, int y) const {
        refresh();
        return view[y, x];
    }

    const std::vector<double>& field() const {   
        refresh();
        return concentration;
    }

    void point_mass(int i, int n) {
        total = n;
        particles.clear();
        particles.reserve(n);
        for (int k = 0; k < n; k++) {
            Particle p;
            p.a = p.b = i;
            p.s = 0;
            p.choose_rand_nb(graph);
            particles.push_back(p);
        }
        dirty = true;
    }

    void step_all(double v) {
        for (auto& p : particles) p.step(graph, v);
        dirty = true;
    }

private:
    void refresh() const {
        if (!dirty) return;
        std::ranges::fill(concentration, 0.0);
        const double w = 1.0 / total;
        for (const auto& p : particles) {
            auto [x, y] = graph.to_cell(p.position(graph));
            view[y, x] += w;
        }
        dirty = false;
    }

    PlanarGraphView& graph;
    int total;
    mutable std::vector<double> concentration;
    std::vector<Particle> particles;
    mutable std::mdspan<double, std::dextents<size_t, 2>> view;   // declared after concentration
    mutable bool dirty = true;
};



















