#include "integrators/Integrator.hpp"
#include "data/ModelData.hpp"
#include "data/Topology.hpp"
#include "data/SharedData.hpp"
#include "forces/Force.hpp"
#include "reporters/Reporter.hpp"
#include "data/Results.hpp"
#include "utils/RNG.hpp"
using namespace std;

namespace cg::reference {
    class Simulation {
    public:
        toolkit::RNG &rng;
        ModelData modelData;
        Topology topology;
        Integrator *integrator = 0;
        Results results;
        SharedData sharedData;

        const Real delta = 0.005 * toolkit::nanosecond; // dt

        Simulation(toolkit::Model const &, toolkit::RNG &);

        void attachForce(Force *);
        void attachReporter(Reporter *, int period = 1);

        void run(int max_steps);

    private:
        Reals3 forces;
        vector<Force*> forceObjects;
        vector<pair<Reporter*, int>> reporters;

        void calcForces();
        void dumpResults();
    };
}
