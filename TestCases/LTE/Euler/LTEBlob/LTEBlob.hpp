#pragma once

#include "ConditionFactory.hpp"

namespace Prandtl
{

// LTE Blob initial condition
std::function<void(const Vector&, Vector&)> LTEBlobIC(real_t radius,
                                                        real_t T_inf,
                                                        real_t T_blob,
                                                        real_t P_inf)
{
    return [=](const Vector &x, Vector &y)
    {
        MFEM_ASSERT(x.Size() == 2, "");

        const real_t xc = 0.0;
        const real_t yc = 0.0;

        real_t dx = x(0) - xc;
        real_t dy = x(1) - yc;

        real_t r2rad = (dx*dx + dy*dy) / (radius * radius);

        const real_t exp_full = std::exp(-r2rad);

        // Quiescent velocity field
        const real_t velX = 0.0;
        const real_t velY = 0.0;
        const real_t vel2 = velX * velX + velY * velY;

        // Gaussian temperature perturbation
        const real_t T = T_inf + (T_blob - T_inf) * exp_full;

        // LTE mixture using Mutation++
        Mutation::MixtureOptions opts("air_5");
        opts.setStateModel("EquilTP");
        opts.setThermodynamicDatabase("RRHO");
        Mutation::Mixture mix(opts);
        mix.addComposition("N:0.79, O:0.21", true);

        mix.setState(&T, &P_inf);

        const real_t den  = mix.density();
        const real_t rhoe = den * mix.mixtureEnergyMass();
        const real_t rhoE = rhoe + 0.5 * den * vel2;

        y(0) = den;
        y(1) = den * velX;
        y(2) = den * velY;
        y(3) = rhoE;
    };
}

// Registration helper that automatically registers these functions
struct RegisterLTEBlob
{
    RegisterLTEBlob()
    {
        // Register initial condition.
        ConditionFactory::Instance().RegisterInitialCondition4("LTEBlobIC", LTEBlobIC);
    }
};
// Global static instance to ensure registration happens at startup.
static RegisterLTEBlob regLTEBlob;

}