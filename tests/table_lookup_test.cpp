// tests/table_lookup_test.cpp
#include "unit_test.hpp"
#include "test_helpers.hpp"

#include "Physics.hpp"
#include "GasState.hpp"
#include "BasicOperations.hpp"
#include "LteEOS.hpp"

#include <mutation++/mutation++.h>
#include "../libs/mfem/general/forall.hpp"

using real_t = Prandtl::real_t;

using namespace Prandtl;
using namespace Mutation;

TEST(hunt_cpu_test)
{
    int n = 100;
    Prandtl::Vector arr(n);
    for (int i = 0; i < n; i++)
    {
        arr[i] = i + 1;
    }

    real_t x = 34.5;

    // Hunt Right (guess near left)
    int ind_lo = 2;
    ind_lo = hunt(arr.Read(), n, x, ind_lo);
    EXPECT_CLOSE(ind_lo, 33, 1e-14);

    // Hunt Left (guess near right)
    ind_lo = 70;
    ind_lo = hunt(arr.Read(), n, x, ind_lo);
    EXPECT_CLOSE(ind_lo, 33, 1e-14);

    // Left Boundary
    x = 1;
    ind_lo = 50;
    ind_lo = hunt(arr.Read(), n, x, ind_lo);
    EXPECT_CLOSE(ind_lo, 0, 1e-14);

    // Right Boundary
    x = n;
    ind_lo = 20;
    ind_lo = hunt(arr.Read(), n, x, ind_lo);
    EXPECT_CLOSE(ind_lo, n-2, 1e-14);

    return 0;
}

TEST(hunt_gpu_test)
{
    const int n = 100;
    mfem::Vector arr(n);
    for (int i = 0; i < n; i++) { arr[i] = i + 1; }

    const real_t *a = arr.Read(); // device-safe read pointer
    arr.UseDevice();

    mfem::Vector outv(4);
    outv.UseDevice();
    real_t *out_d = outv.Write(); // device-safe write pointer

    mfem::forall(4, [=] MFEM_HOST_DEVICE (int i)
    {
        real_t x;
        int guess;

        if (i == 0)      { x = 34.5;  guess = 2;  }   // hunt right
        else if (i == 1) { x = 34.5;  guess = 70; }   // hunt left
        else if (i == 2) { x = 1.0;   guess = 50; }   // left boundary
        else             { x = 100.0; guess = 20; }   // right boundary

        int idx = hunt(a, n, x, guess);
        out_d[i] = (real_t) idx;
    });

    const real_t *out_h = outv.HostRead();

    EXPECT_CLOSE(out_h[0], 33.0, 1e-14);
    EXPECT_CLOSE(out_h[1], 33.0, 1e-14);
    EXPECT_CLOSE(out_h[2],  0.0, 1e-14);
    EXPECT_CLOSE(out_h[3], 98.0, 1e-14); // n-2

    return 0;
}

// -----------------------------------------------------------------------------
// EOS test : Thermodynamic properties of a gas-mixture in LTE
// The lte table is populated using Mutation++ property = property(rho x rhoe)
// -----------------------------------------------------------------------------
TEST(LTEGasEOS_Tablelookup_test)
{
    // Range and resolution of the table (Evenly spaced table in this test)
    int nx = 11, ny = 11;
    real_t rho_min  = 0.5  , rho_max  = 1.5  , rho_step  = (rho_max-rho_min)/(nx-1);
    real_t rhoe_min = 2.0e5, rhoe_max = 1.1e6, rhoe_step = (rhoe_max-rhoe_min)/(ny-1);

    // Mutation++ object for Equilibrium state
    MixtureOptions opts("air_5");          // Gas-Mixture (N, O, NO, N2, O2)
    opts.setStateModel("Equil");           // Equilibrium state model
    opts.setThermodynamicDatabase("RRHO"); // Default thermodynamic data base
    opts.setViscosityAlgorithm("Chapmann-Enskog_LDLT"); // Viscosity algorithm
    Mixture mix(opts);                     // Initializing mixture object
        mix.addComposition("N:0.8, O:0.2", true); // composition
    int num_properties = 9;

    const int dim = 3, ndofs = 1;

    // Data arrays for table lookup for LTE properties
    mfem::Vector lte_table( (num_properties) * (nx*ny) ); // LTE-Tables
    mfem::Vector rho_grid(nx), rhoe_grid(ny); // 1-D grids of rho and rhoE

    StateLayout L(dim, ndofs, nx, ny);
    const int num_eq = L.eq_energy + 1;
    LTEGasEOS eos;

    // Populating the 1D grid of density and rho*internal energy
    for(int ind_x=0; ind_x < nx; ind_x++) rho_grid[ind_x] = rho_min + ind_x * rho_step;
    for(int ind_y=0; ind_y < ny; ind_y++) rhoe_grid[ind_y] = rhoe_min + ind_y * rhoe_step;

    // Generation of LTE table
    for(int ind_y=0; ind_y < ny; ind_y++)
    {
        for(int ind_x=0; ind_x < nx; ind_x++)
        {
            mix.setState(&rho_grid[ind_x], &rhoe_grid[ind_y], 0);
            lte_table[L.lte_property_index(L.P_idx, ind_x, ind_y)] = mix.P();
            lte_table[L.lte_property_index(L.T_idx, ind_x, ind_y)] = mix.T();
        }
    }

    // Setup the LTE EOS
    std::shared_ptr<PhysicsConstants> phys =
      std::make_shared<PhysicsConstants>(lte_table.Read(), rho_grid.Read(), rhoe_grid.Read());

    // ------------------------------------------------------------------------------------


    // Populate the 1-DOF state (Ensure the values are in ranges of the table)
    real_t rho  = 0.75; // 0.753
    real_t rhoe = 580000.0; // 600000.0
    const real_t u1[3] = {10.0, -3.0, 5.0};

    std::vector<real_t> U(num_eq * ndofs);
    fill_single_dof_state(L, U, dim, rho, u1, rhoe);
    DofStateView S1(U.data(), 0);

    // Populate the hunt array with initial guesses for the hunt algorithm
    int l_x = hunt(rho_grid.Read(), nx, rho, 0);
    int l_y = hunt(rhoe_grid.Read(), ny, rhoe, 0);
    int u_x = l_x + 1      ; int u_y = l_y + 1;
    EXPECT_CLOSE(l_x, 2, 0);
    EXPECT_CLOSE(l_y, 4, 0);

    int P_idx = L.P_idx, T_idx = L.T_idx;

    // Test to check the values of the corners of the 2D cell in P and T tables
    real_t P_true, T_true, P_table, T_table;
    for(int i=0; i < 2; i++)
    {
        for(int j=0; j < 2; j++)
        {
            mix.setState(&rho_grid[2 + i], &rhoe_grid[4 + j], 0);
            P_true = mix.P();   P_table = phys->lte_table[L.lte_property_index(P_idx, l_x + i, l_y + j)];
            T_true = mix.T();   T_table = phys->lte_table[L.lte_property_index(T_idx, l_x + i, l_y + j)];
            EXPECT_CLOSE(P_true, P_table, 1e-6);
            EXPECT_CLOSE(T_true, T_table, 1e-6);
        }
    }

    return 0;
}

TEST(LTEGasEOS_BilinearInterpolation_test)
{
    // Range and resolution of the table (Evenly spaced table in this test)
    int nx = 11, ny = 1001;
    real_t rho_min  = 0.5  , rho_max  = 1.5  , rho_step  = (rho_max-rho_min)/(nx-1);
    real_t e_min = 2.0e5, e_max = 1.2e6, e_step = (e_max-e_min)/(ny-1);

    // Mutation++ object for Equilibrium state
    MixtureOptions opts("air_5");          // Gas-Mixture (N, O, NO, N2, O2)
    opts.setStateModel("Equil");           // Equilibrium state model
    opts.setThermodynamicDatabase("RRHO"); // Default thermodynamic data base
    opts.setViscosityAlgorithm("Chapmann-Enskog_LDLT"); // Viscosity algorithm
    Mixture mix(opts);                     // Initializing mixture object
        mix.addComposition("N:0.8, O:0.2", true); // composition
    int num_properties = 9;

    const int dim = 3, ndofs = 1;

    // Data arrays for table lookup for LTE properties
    mfem::Vector lte_table( (num_properties) * (nx*ny) ); // LTE-Tables
    mfem::Vector rho_grid(nx), e_grid(ny); // 1-D grids of rho and rhoE

    StateLayout L(dim, ndofs, nx, ny);
    const int num_eq = L.eq_energy + 1;
    LTEGasEOS eos;

    // Populating the 1D grid of density and rho*internal energy
    for(int ind_x=0; ind_x < nx; ind_x++) rho_grid[ind_x] = rho_min + ind_x * rho_step;
    for(int ind_y=0; ind_y < ny; ind_y++) e_grid[ind_y] = e_min + ind_y * e_step + mix.mixtureHMass(0.0000001);

    // Generation of LTE table
    for(int ind_y=0; ind_y < ny; ind_y++)
    {
        for(int ind_x=0; ind_x < nx; ind_x++)
        {
            real_t rhoe = rho_grid[ind_x] * e_grid[ind_y];
            mix.setState(&rho_grid[ind_x], &rhoe, 0);
            lte_table[L.lte_property_index(L.P_idx, ind_x, ind_y)] = mix.P();
            lte_table[L.lte_property_index(L.T_idx, ind_x, ind_y)] = mix.T();
        }
    }

    // Setup the LTE EOS
    std::shared_ptr<PhysicsConstants> phys =
      std::make_shared<PhysicsConstants>(lte_table.Read(), rho_grid.Read(), e_grid.Read());

    // Testing Bi-linear Interpolation 
    // A (Choosing one corner)
    real_t rho  = rho_grid[5];
    real_t rhoe = rho*(e_grid[5]);
    const real_t u1[3] = {10.0, -3.0, 5.0};

    std::vector<real_t> U(num_eq * ndofs);
    fill_single_dof_state(L, U, dim, rho, u1, rhoe);
    DofStateView S1(U.data(), 0);
    
    real_t P_interpolated = eos.pressure(*phys, L, S1);
    real_t T_interpolated = eos.temperature(*phys, L, S1);

    mix.setState(&rho, &rhoe, 0);
    real_t P_true = mix.P();
    real_t T_true = mix.T();

    // std::cout << "P_true - P_table: " << (P_true - P_interpolated)/P_true << std::endl;
    // std::cout << "T_true - T_table: " << (T_true - T_interpolated)/T_true << std::endl;

    EXPECT_CLOSE(P_true/P_true, P_interpolated/P_true, 1e-12);
    EXPECT_CLOSE(T_true/T_true, T_interpolated/T_true, 1e-12);


    // Testing Bi-linear Interpolation 
    // B (Mid point of a cell)
    rho  = 0.75;
    real_t e = (e_grid[5] + e_grid[6])/2.0;
    rhoe = rho*e;
    
    fill_single_dof_state(L, U, dim, rho, u1, rhoe);
    DofStateView S2(U.data(), 0);

    P_interpolated = eos.pressure(*phys, L, S2);
    T_interpolated = eos.temperature(*phys, L, S2);

    int l_x = hunt(rho_grid.Read(), nx, rho, 0), l_y = hunt(e_grid.Read(), ny, e, 0);

    P_true = 0;
    T_true = 0;
    for(int i=0; i < 2; i++)
    {
        for(int j=0; j < 2; j++)
        {
            P_true += phys->lte_table[L.lte_property_index(L.P_idx, l_x + i, l_y + j)];
            T_true += phys->lte_table[L.lte_property_index(L.T_idx, l_x + i, l_y + j)];
        }
    }

    P_true /= 4.0;
    T_true /= 4.0;

    EXPECT_CLOSE(P_true, P_interpolated, 0.0);
    EXPECT_CLOSE(T_true, T_interpolated, 0.0);

    // Testing the accuracy of this approach at an arbitary point
    rho = 0.2*rho_grid[5] + 0.8*rho_grid[6];
    e = 0.3*e_grid[5] + 0.7*e_grid[6];
    rhoe = rho*e;

    fill_single_dof_state(L, U, dim, rho, u1, rhoe);
    DofStateView S3(U.data(), 0);
    P_interpolated = eos.pressure(*phys, L, S3);
    T_interpolated = eos.temperature(*phys, L, S3);

    mix.setState(&rho, &rhoe, 0);
    P_true = mix.P();
    T_true = mix.T();

    // std::cout << "P_true - P_table: " << (P_true - P_interpolated)/P_true << std::endl;
    // std::cout << "T_true - T_table: " << (T_true - T_interpolated)/T_true << std::endl;

    EXPECT_CLOSE(P_true/P_true, P_interpolated/P_true, 1e-7);
    EXPECT_CLOSE(T_true/T_true, T_interpolated/T_true, 1e-7);

    return 0;
}

TEST(InverseLTETable_test)
{
    // -------------------------------- Mixtures Definitions -----------------------
    // Mutation++ object for Equilibrium state
    MixtureOptions opts("air_5");
    opts.setStateModel("Equil");
    opts.setThermodynamicDatabase("RRHO");
    opts.setViscosityAlgorithm("Chapmann-Enskog_LDLT");
    Mixture mix1(opts);
    mix1.addComposition("N:0.8, O:0.2", true);

    MixtureOptions opts2("air_5");
    opts2.setStateModel("EquilTP");
    opts2.setThermodynamicDatabase("RRHO");
    Mixture mix2(opts2);
    mix2.addComposition("N:0.8, O:0.2", true);

    // -------------------------------- LTE Table Generation -----------------------

    // Range and resolution of the table (Evenly spaced table in this test)
    int nx = 100, ny = 100;
    real_t rho_min  = 0.5  , rho_max  = 1.5  , rho_step  = (rho_max-rho_min)/(nx-1);
    real_t e_min = 2.0e5, e_max = 2.0e6, e_step = (e_max-e_min)/(ny-1);

    int num_properties = 9;
    const int dim = 3, ndofs = 1;

    // Data arrays for table lookup for LTE properties
    mfem::Vector lte_table( (num_properties) * (nx*ny) ); // LTE-Tables
    mfem::Vector rho_grid(nx), e_grid(ny); // 1-D grids of rho and e

    StateLayout L(dim, ndofs, nx, ny);
    const int num_eq = L.eq_energy + 1;
    LTEGasEOS eos;

    // Populating the 1D grid of density and rho*internal energy
    for(int ind_x=0; ind_x < nx; ind_x++) rho_grid[ind_x] = rho_min + ind_x * rho_step;
    for(int ind_y=0; ind_y < ny; ind_y++) e_grid[ind_y] = e_min + ind_y * e_step + mix1.mixtureHMass(0.0000001);

    // Generation of LTE table
    for(int ind_y=0; ind_y < ny; ind_y++)
    {
        for(int ind_x=0; ind_x < nx; ind_x++)
        {
            mix1.setState(&rho_grid[ind_x], &e_grid[ind_y], 0);
            lte_table[L.lte_property_index(L.P_idx, ind_x, ind_y)] = mix1.P();
            lte_table[L.lte_property_index(L.T_idx, ind_x, ind_y)] = mix1.T();
        }
    }

    std::shared_ptr<PhysicsConstants> phys =
      std::make_shared<PhysicsConstants>(lte_table.Read(), rho_grid.Read(), e_grid.Read());

    // -------------------------------- LTE TABLE Look-up --------------------------

    real_t rho_true  = 0.753;
    real_t rhoe_true = rho_true * (580000.0 + mix1.mixtureHMass(0.0000001)); 
    const real_t u1[3] = {10.0, -3.0, 5.0};

    std::vector<real_t> U(num_eq * ndofs);
    fill_single_dof_state(L, U, dim, rho_true, u1, rhoe_true);
    DofStateView S1(U.data(), 0);

    mix1.setState(&rho_true, &rhoe_true, 0);
    real_t P_true = mix1.P();
    real_t T_true = mix1.T();

    real_t P_interpolated = eos.pressure(*phys, L, S1);
    real_t T_interpolated = eos.temperature(*phys, L, S1);


    // ------------------------------------ TESTS ----------------------------------

    // TEST 1 - Mutation++ Inverse
    mix2.setState(&T_true, &P_true, 0);
    real_t rho_inverse = mix2.density();
    EXPECT_CLOSE(rho_true, rho_inverse, 1e-12);

    // TEST 2 - Inverse using Interpolated Pressure and Temperature
    mix2.setState(&T_interpolated, &P_interpolated, 0);
    real_t rho_inverse_interp = mix2.density();
    EXPECT_CLOSE(rho_true, rho_inverse_interp, 1e-5);

    // TEST 3 - Obtaining internal energy from the from pressure (inverse table lookup)
    real_t rhoe_new = rho_true*( 700000.0 + mix1.mixtureHMass(0.0000001) ); // Initial guess for rho*e
    fill_single_dof_state(L, U, dim, rho_true, u1, rhoe_new);
    PointStateView S(U.data());
    real_t rhoe_inverse = eos.internal_energy_from_pressure(*phys, L, S, P_interpolated);
    EXPECT_CLOSE(rhoe_true, rhoe_inverse, 1e-8);

    return 0;
}

TEST(LTETable_diagnostics_test)
{
    MixtureOptions opts("air_5");
    opts.setStateModel("Equil");
    opts.setThermodynamicDatabase("RRHO");
    opts.setViscosityAlgorithm("Chapmann-Enskog_LDLT");
    Mixture mix1(opts);
    mix1.addComposition("N:0.79, O:0.21", true);

    MixtureOptions opts2("air_5");
    opts2.setStateModel("EquilTP");
    opts2.setThermodynamicDatabase("RRHO");
    Mixture mix2(opts2);
    mix2.addComposition("N:0.79, O:0.21", true);

    double T = 8500.0;
    double P = 1e6;

    for(double T = 5000.0; T <= 8000.0; T += 1000.0){
        for(double P = 1e5; P <= 1e6; P *= 10){
            mix2.setState(&T, &P);
            double rho = mix2.density();
            double rhoe = rho * mix2.mixtureEnergyMass();

            // std::cout << "(T,P) = " << "(" << T << ", " << P << ")" << std::endl;
            // CL ALERT : For T>= 8500 K and P >= 1e6 no solution is obtained from Mutation++ using (rho, rhoe)
            mix1.setState(&rho, &rhoe, 0);
            double P_recovered = mix1.P();
            double T_recovered = mix1.T();

            EXPECT_CLOSE(P, P_recovered, 1e-5);
            EXPECT_CLOSE(T, T_recovered, 1e-8);
        }
    }

    return 0;
}

TEST(TableSelfConsistency_test)
{
    // -------------------------------- Mixtures Definitions -----------------------
    // Mutation++ object for Equilibrium state
    MixtureOptions opts("air_5");
    opts.setStateModel("Equil");
    opts.setThermodynamicDatabase("RRHO");
    opts.setViscosityAlgorithm("Chapmann-Enskog_LDLT");
    Mixture mix1(opts);
    mix1.addComposition("N:0.8, O:0.2", true);

    // -------------------------------- LTE Table Generation -----------------------

    // Range and resolution of the table (Evenly spaced table in this test)
    int nx = 100, ny = 100;
    real_t rho_min  = 0.5  , rho_max  = 1.5  , rho_step  = (rho_max-rho_min)/(nx-1);
    real_t e_min = 2.0e5, e_max = 2.0e6, e_step = (e_max-e_min)/(ny-1);

    int num_properties = 9;
    const int dim = 3, ndofs = 1;

    // Data arrays for table lookup for LTE properties
    mfem::Vector lte_table( (num_properties) * (nx*ny) ); // LTE-Tables
    mfem::Vector rho_grid(nx), e_grid(ny); // 1-D grids of rho and e
    mfem::Vector P_table(nx*ny); // 2D tables for P

    StateLayout L(dim, ndofs, nx, ny);
    const int num_eq = L.eq_energy + 1;
    LTEGasEOS eos;

    // Populating the 1D grid of density and internal energy
    for(int ind_x=0; ind_x < nx; ind_x++) rho_grid[ind_x] = rho_min + ind_x * rho_step;
    for(int ind_y=0; ind_y < ny; ind_y++) e_grid[ind_y] = e_min + ind_y * e_step + mix1.mixtureHMass(0.0000001);

    // Generation of LTE table
    for(int ind_y=0; ind_y < ny; ind_y++)
    {
        for(int ind_x=0; ind_x < nx; ind_x++)
        {
            mix1.setState(&rho_grid[ind_x], &e_grid[ind_y], 0);
            lte_table[L.lte_property_index(L.P_idx, ind_x, ind_y)] = mix1.P();
            lte_table[L.lte_property_index(L.T_idx, ind_x, ind_y)] = mix1.T();

            P_table[ind_y*nx + ind_x] = mix1.P();
        }
    }

    std::shared_ptr<PhysicsConstants> phys =
      std::make_shared<PhysicsConstants>(lte_table.Read(), rho_grid.Read(), e_grid.Read());

    const real_t u1[3] = {10.0, -3.0, 5.0};
    std::vector<real_t> U(num_eq * ndofs);

    // TEST 1 : Set (rho,rhoe) Obtain P and then see if we can obtain same rhoe from P using the inverse table lookup
    for(real_t rho_true = 0.7532; rho_true <= 0.8982 ; rho_true += 0.0001)
    {
        // real_t rho_true  = 0.753;
        real_t rhoe_true = rho_true * (580000.0 + mix1.mixtureHMass(0.0000001));

        fill_single_dof_state(L, U, dim, rho_true, u1, rhoe_true);
        DofStateView S1(U.data(), 0);

        real_t P_interpolated = eos.pressure(*phys, L, S1);
        real_t T_interpolated = eos.temperature(*phys, L, S1);


        real_t rhoe_new = rho_true*( 700000.0 + mix1.mixtureHMass(0.0000001) ); // Initial guess for rho*e
        fill_single_dof_state(L, U, dim, rho_true, u1, rhoe_new);
        PointStateView S(U.data());
        real_t rhoe_inverse = eos.internal_energy_from_pressure(*phys, L, S, P_interpolated);
        EXPECT_CLOSE(rhoe_true, rhoe_inverse, 1e-9);
    }

    // TEST 2 : Set (rho, P) and obtain rhoe from inverse table lookup and then see if we can obtain same P from (rho, rhoe)
    for(real_t fac=0.235; fac <= 0.7; fac += 0.001)
    {
        real_t P_target = P_table.Min() + fac*(P_table.Max() - P_table.Min());

        real_t rho_true = 1.0;
        fill_single_dof_state(L, U, dim, rho_true, u1, rho_true*e_grid[50]);
        DofStateView S1(U.data(), 0);
        real_t rhoe_inverse = eos.internal_energy_from_pressure(*phys, L, S1, P_target);

        fill_single_dof_state(L, U, dim, rho_true, u1, rhoe_inverse);
        DofStateView S2(U.data(), 0);
        real_t P_interpolated = eos.pressure(*phys, L, S2);
        // std::cout << "P_target - P_interpolated = " << (P_target - P_interpolated)/P_target << std::endl;
        EXPECT_CLOSE(P_target/P_target, P_interpolated/P_target, 1e-14);
    }

    return 0;
}