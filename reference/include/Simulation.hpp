#include "integrators/Integrator.hpp"
#include "data/ModelData.hpp"
#include "data/Topology.hpp"
#include "forces/Force.hpp"
#include "reporters/Reporter.hpp"
#include "data/Results.hpp"
using namespace std;

namespace cg::reference {
    class Simulation {
    public:
        ModelData modelData;
        Topology topology;
        Integrator *integrator = 0;
        Results results;
        const Real delta = 0.001; // dt

        Simulation(toolkit::Model const&);

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
