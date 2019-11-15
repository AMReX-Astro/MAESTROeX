void Put1dSpecArrayOnCart (const RealVector& s0,
                           amrex::Vector<amrex::MultiFab>& s0_cart,
                           int is_input_edge_centered,
                           const amrex::Vector<amrex::BCRec>& bcs = amrex::Vector<amrex::BCRec>(),
                           int sbccomp = 0, int variable_type=0);

void Put1dSpecArrayOnCart (int level, const RealVector& s0,
                        amrex::Vector<amrex::MultiFab>& s0_cart,
                        int is_input_edge_centered,
                        const amrex::Vector<amrex::BCRec>& bcs = amrex::Vector<amrex::BCRec>(),
                        int sbccomp = 0);