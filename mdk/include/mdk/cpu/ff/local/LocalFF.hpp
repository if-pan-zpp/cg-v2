#pragma once
#include <mdk/cpu/ff/ForceField.hpp>
#include <mdk/cpu/core/Particles.hpp>

namespace mdk {
    template<typename T>
    class LocalFF: public ForceField<LocalFF<T>> {
    private:
        Particles const& pts;

    public:
        LocalFF(Particles const& pts):
            pts {pts} {};

        void evalFF(double &V, Vectors& gradV) const {
            for (auto const& [beg, end]: pts->chains) {
                for (int idx = beg; idx < end; ++idx) {
                    localTerm(V, dV, idx);
                }
            }
        }

        void localTerm(double &V, Vectors& gradV, int idx, int beg, int end) const {
            static_cast<T const&>(*this).kernel(V, gradV, idx, beg, end);
        }
    };

    template<typename T1, typename T2>
    class SumLocalFF: public LocalFF<SumLocalFF<T1, T2>> {
    private:
        T1 const& ff1;
        T2 const& ff2;

    public:
        SumLocalFF(T1 const& ff1, T2 const& ff2):
            ff1 {ff1}, ff2 {ff2} {};

        void localTerm(double &V, Vectors& gradV, int idx, int beg, int end) const {
            ff1.localTerm(V, gradV, idx);
            ff2.localTerm(V, gradV, idx);
        }
    };

    template<typename T1, typename T2>
    SumLocalFF<T1, T2> operator+(LocalFF<T1> const& ff1, LocalFF<T2> const& ff2) {
        return SumLocalFF(static_cast<T1 const&>(ff1), static_cast<T2 const&>(ff2));
    }
}