
#include <BaseState.H>
using namespace amrex;

std::string inputs_name = "";

int main(int argc, char* argv[])
{

    // in AMReX.cpp
    Initialize(argc, argv);

    // timer for profiling
    BL_PROFILE_VAR("main()", main);

    // wallclock time
    const Real strt_total = ParallelDescriptor::second();

    {

        // we need to make a dummy mf for the MFIter
        RealBox real_box;
        for (int n = 0; n < BL_SPACEDIM; n++)
        {
            real_box.setLo(n, 0.0);
            real_box.setHi(n, 1.0);
        }

        IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
        IntVect domain_hi(AMREX_D_DECL(1,1,1));
        const Box domain(domain_lo, domain_hi);

        int coord = 0;
        int is_per[BL_SPACEDIM];
        for (int i = 0; i < BL_SPACEDIM; i++)
            is_per[i] = 0;
        Geometry geom(domain, &real_box, coord, is_per);
        
        BoxArray ba(domain);
        ba.maxSize(16);
        DistributionMapping dm(ba);        

        MultiFab dummy_mf(ba, dm, 1, 0);
        // const Box& tile_box  = mfi.tilebox();

        // BaseState<Real> base_state;

        const int nlevs = 3;
        const int len = 20;
        const int ncomp = 5;

        Print() << "Defining base state with " << nlevs << " levels, length " << len << " and " << ncomp << " components." << std::endl;

        // base_state.define(nlevs, len, ncomp);
        // BaseState<Real> base_state(nlevs, len, ncomp);

    // for(MFIter mfi = MFIter(dummy_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    // {
        BaseState<Real> base_state(nlevs,len,ncomp);

        Print() << "setting value to 1.0" << std:: endl;

        // amrex::Gpu::ManagedVector<Real> base_data(ncomp*len*nlevs);


        // Real * AMREX_RESTRICT ptr = base_state.base_data.dataPtr();
        // AMREX_PARALLEL_FOR_1D(ncomp*len*nlevs, i, {
        //     ptr[i] = 1.0;
        // });
        base_state.setVal(1.0);

        Print() << "print an element: " << base_state(0,0) << std::endl;

        Print() << "defining another base state" << std::endl;
        BaseState<Real> other_base_state(nlevs, len, ncomp);

        amrex::Gpu::synchronize();

        for (auto comp = 0; comp < ncomp; ++comp) {
            other_base_state.setVal(comp, 1.0);


            amrex::Gpu::synchronize();
        }

        // Print() << "are the two base states equal? " << (base_state == other_base_state) << std::endl;

        // Print() << "making a deep copy" << std::endl;

        // BaseState<Real> deep_copy(base_state);

        // Print() << "are the two base states equal? " << (base_state == deep_copy) << std::endl;

        // Print() << "Change the copy" << std::endl;

        // deep_copy.setVal(1, 0.0);

        // Print() << "print an element: " << deep_copy(0,0,1) << std::endl;

        // Print() << "are the two base states equal? " << (base_state == deep_copy) << std::endl;

        // Print() << "let's iterate over the base state" << std::endl;
        // for (auto l = 0; l < nlevs; ++l) {
        //     AMREX_PARALLEL_FOR_1D(base_state.length(), n, {
        //         base_state(l,n) = 0.5;
        //     });
        // }

        // Print() << "we can also iterate over the components" << std::endl;
        // for (auto l = 0; l < nlevs; ++l) {
        //     AMREX_PARALLEL_FOR_1D(base_state.length(), n, {
        //         for (auto comp = 0; comp < base_state.nComp(); ++comp) {
        //             base_state(l,n,comp) = Real(comp);
        //         }
        //     });
        // }

        // Print() << "multiply by scalar" << std::endl;

        // base_state *= 2.0;

        // Print() << "add two base states" << std::endl;

        // BaseState<Real> summed_state = base_state + other_base_state;

        // Print() << "subtract in place" << std::endl;

        // base_state -= other_base_state;
    // }
    }

    // destroy timer for profiling
    BL_PROFILE_VAR_STOP(main);

    // in AMReX.cpp
    Finalize();
}
