#include "unit_test.hpp"
#include "test_helpers.hpp"
#include "plato_helpers.hpp"

#include "Physics.hpp"
#include "GasState.hpp"
#include "BasicOperations.hpp"
#include "LteEOS.hpp"

using real_t = Prandtl::real_t;

using namespace Prandtl;

// -----------------------------------------------------------------------------
// EOS test : Thermodynamic properties of a gas-mixture in LTE
// The lte table is populated using PLATO
// -----------------------------------------------------------------------------
TEST(PLATO_Tablelookup_test)
{
    Mpi::Init();

    // Range and resolution of the table (Evenly spaced table in this test)
    int nx = 101, ny = 101;
    real_t rho_min  = 0.01  , rho_max  = 0.11 ;
    real_t T_min = 250.0, T_max = 350.0;

    int num_properties = 9;

    mfem::Vector lte_table( (num_properties) * (nx*ny) ), inv_table( nx*ny );
    mfem::Vector rho_grid, T_grid, e_grid;

    const int dim = 3, ndofs = 1;
    StateLayout L(dim, ndofs, nx, ny);
    const int num_eq = L.eq_energy + 1;
    const real_t u1[3] = {10.0, -3.0, 5.0};
    std::vector<real_t> U(num_eq * ndofs);

    // Generation of LTE table
        std::string solver  = "LTE_table_rhoT_(air5)";
        std::string mixture = "air5";
        std::string path    = "/home/cherith2/Workspace/SOURCE_CODES/database";

        uniform_grid(nx, rho_min, rho_max, rho_grid);
        uniform_grid(ny, T_min, T_max, T_grid);
        real_t e_min, e_max;

        fill_lte_table(L, rho_grid.GetData(), T_grid.GetData(), 
                        solver.c_str(), mixture.c_str(), path.c_str(),
                        lte_table.GetData(), e_min, e_max, MPI_COMM_SELF);

        uniform_grid(ny, e_min, e_max, e_grid);
        fill_inv_table(L, rho_grid.GetData(), e_grid.GetData(), inv_table.GetData(), MPI_COMM_SELF);
        std::shared_ptr<PhysicsConstants> phys =std::make_shared<PhysicsConstants>(lte_table.HostRead(), inv_table.HostRead(),
                                                rho_grid.HostRead(), T_grid.HostRead(), e_grid.HostRead());
        LTEGasEOS eos;

    // ------------------------------------------------------------------------------------
    
    // PLATO setup
        std::string empty_str = "empty";
        plato_initialize(solver.c_str(), mixture.c_str(), empty_str.c_str(), empty_str.c_str(), path.c_str());
        PlatoMixture mix;

    // Test 1 : (Table values at the corner of the table)
        real_t rho  = rho_grid[3];
        real_t T = T_grid[7];
        plato_set_state(rho, T, mix);
        real_t e = mix.e;

        real_t rhoe = rho * e;
        fill_single_dof_state(L, U, dim, rho, u1, rhoe);
        DofStateView S1(U.data(), 0);

        real_t P_table  = eos.pressure(*phys, L, S1);
        real_t P_corner = lte_table[L.lte_property_index(L.P_idx, 3, 7)]; 
        real_t rel_err = std::abs(mix.P - P_table)/std::abs(mix.P);

        EXPECT_SMALL(rel_err, 1e-15);
        EXPECT_CLOSE(mix.P, P_corner, 0.0);

    // TEST 2 : (Table values at a mid-point of a cell in the table)
        rho = 0.5*(rho_grid[3] + rho_grid[4]);
        T   = 0.5*(T_grid[7] + T_grid[8]);
        plato_set_state(rho, T, mix);
        e = mix.e;

        rhoe = rho * e;
        fill_single_dof_state(L, U, dim, rho, u1, rhoe);
        DofStateView S2(U.data(), 0);

        P_table  = eos.pressure(*phys, L, S2);
        rel_err = std::abs(mix.P - P_table)/std::abs(mix.P);

        real_t P_expected = 0;
        int l_x = hunt(rho_grid.Read(), nx, rho, 0), l_y = hunt(T_grid.Read(), ny, T, 0);
        for(int i=0; i < 2; i++)
        {
            for(int j=0; j < 2; j++)
            {
                P_expected += lte_table[L.lte_property_index(L.P_idx, l_x + i, l_y + j)];
            }
        }
        P_expected /= 4.0;

        EXPECT_SMALL(rel_err, 1e-7);
        EXPECT_CLOSE(P_table/P_expected, 1.0, 1e-7);

    // TEST 3 : (Table values at an arbitrary point in the table)
        rho = 0.233525*rho_grid[3] + 0.766475*rho_grid[4];
        T   = 0.768256*T_grid[7] + 0.231744*T_grid[8];
        plato_set_state(rho, T, mix);
        e = mix.e;
        rhoe = rho * e;
        fill_single_dof_state(L, U, dim, rho, u1, rhoe);
        DofStateView S3(U.data(), 0);
        P_table  = eos.pressure(*phys, L, S3);
        rel_err = std::abs(mix.P - P_table)/std::abs(mix.P);
        EXPECT_SMALL(rel_err, 1e-7);

    plato_finalize();
    return 0;
}