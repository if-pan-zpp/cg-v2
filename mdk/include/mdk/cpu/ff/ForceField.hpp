#pragma once
#include <mdk/cpu/core/Types.hpp>

namespace mdk {
    template<typename T>
    class ForceField {
    public:
        void evalFF(double& V, Vectors& gradV) const {
            static_cast<T const&>(*this).evalFF(V, dV);
        }
    };

    template<typename T1, typename T2>
    class SumFF: public ForceField<SumFF<T1, T2>> {
    private:
        T1 const& ff1;
        T2 const& ff2;

    public:
        SumFF(T1 const& ff1, T2 const& ff2):
            ff1 {ff1}, ff2 {ff2} {};

        void evalFF(double& V, Vectors& gradV) const {
            ff1.evalFF(V, gradV);
            ff2.evalFF(V, gradV);
        }
    };

    template<typename T1, typename T2>
    SumFF<T1, T2> operator+(ForceField<T1> const& ff1, ForceField<T2> const& ff2) {
        return SumFF(static_cast<T1 const&>(ff1), static_cast<T2 const&>(ff2));
    }
}
