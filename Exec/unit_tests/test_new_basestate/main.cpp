
#include <BaseState.H>
using namespace amrex;

std::string inputs_name = "";

int main(int argc, char* argv[]) {
    // in AMReX.cpp
    Initialize(argc, argv);

    // timer for profiling
    BL_PROFILE_VAR("main()", main);

    // wallclock time
    const Real strt_total = ParallelDescriptor::second();

    {
        BaseState<Real> base_state;

        const int nlevs = 3;
        const int len = 20;
        const int ncomp = 5;

        Print() << "Defining base state with " << nlevs << " levels, length "
                << len << " and " << ncomp << " components." << std::endl;

        base_state.define(nlevs, len, ncomp);

        Print() << "setting value to 1.0" << std::endl;

        base_state.setVal(1.0);

        Print() << "print an element: " << base_state.array()(0, 0)
                << std::endl;

        Print() << "defining another base state" << std::endl;

        BaseState<Real> other_base_state(nlevs, len, ncomp);

        for (auto comp = 0; comp < ncomp; ++comp) {
            other_base_state.setVal(comp, 1.0);
        }

        Print() << "are the two base states equal? "
                << (base_state == other_base_state) << std::endl;

        Print() << "making a deep copy" << std::endl;

        BaseState<Real> deep_copy(base_state);

        Print() << "are the two base states equal? "
                << (base_state == deep_copy) << std::endl;

        Print() << "Change the copy" << std::endl;

        deep_copy.setVal(1, 0.0);

        Print() << "print an element: " << deep_copy.array()(0, 0, 1)
                << std::endl;

        Print() << "are the two base states equal? "
                << (base_state == deep_copy) << std::endl;

        Print() << "Initialize from a Fortran-style 1d base state" << std::endl;

        Gpu::ManagedVector<Real> f_base(nlevs * len * ncomp);

        BaseState<Real> c_base(f_base, nlevs, len, ncomp);

        Print() << "let's iterate over the base state using a BaseStateArray"
                << std::endl;
        BaseStateArray<Real> base_arr = base_state.array();
        for (auto l = 0; l < nlevs; ++l) {
            AMREX_PARALLEL_FOR_1D(len, n, { base_arr(l, n) = 0.5; });
        }
        Gpu::synchronize();

        Print() << "base_state = " << base_arr(0, 0) << std::endl;

        Print() << "we can also iterate over the components" << std::endl;
        for (auto l = 0; l < nlevs; ++l) {
            AMREX_PARALLEL_FOR_1D(len, n, {
                for (auto comp = 0; comp < ncomp; ++comp) {
                    base_arr(l, n, comp) = Real(comp);
                }
            });
        }
        Gpu::synchronize();

        Print() << "base_state = " << base_arr(0, 0, 0) << std::endl;

        Print() << "multiply by scalar" << std::endl;

        base_state *= 2.0;
        base_state += 2.0;

        Print() << "add two base states" << std::endl;

        BaseState<Real> summed_state = base_state + other_base_state;

        Print() << "subtract in place" << std::endl;

        base_state -= other_base_state;
    }

    // destroy timer for profiling
    BL_PROFILE_VAR_STOP(main);

    // in AMReX.cpp
    Finalize();
}
