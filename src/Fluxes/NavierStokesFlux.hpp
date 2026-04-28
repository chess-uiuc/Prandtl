#pragma once

#include "mfem.hpp"
#include "GasModel.hpp"

namespace Prandtl
{

using namespace mfem;

class NavierStokesFlux : public FluxFunction
{
private:
  const ActiveGasModel gasModel;
public:
  explicit NavierStokesFlux(const ActiveGasModel &gasModel_)
    : FluxFunction(gasModel_.num_equations(), gasModel_.dim()), gasModel(gasModel_){};
  void ComputeViscousFlux(const Vector &state, const Vector &dqdx, const Vector &dqdy, const Vector &dqdz, DenseMatrix &flux) const;
  void ComputeViscousFlux(const Vector &state, const Vector &dqdx, const Vector &dqdy, DenseMatrix &flux) const;
  void ComputeViscousFlux(const Vector &state, const Vector &dqdx, DenseMatrix &flux) const;
  MFEM_HOST_DEVICE inline real_t pressure(const real_t *state) const
  {
    Prandtl::PointStateView S{state};
    return gasModel.pressure(S);
  }

  /**
   * @brief Compute inviscid flux from conserved state
   *
   * @param state conserved state at current integration point
   * @param Tr current element transformation with the integration point
   * @param flux inviscid flux (ex, ideal single gas: F(ρ, ρu, E) = [ρuᵀ; ρuuᵀ + pI; uᵀ(E + p)])
   * @return real_t maximum characteristic speed, c + |u| (c = speed of sound)
   */
  real_t ComputeFlux(const Vector &state, ElementTransformation &Tr,
                     DenseMatrix &flux) const override;
  
  /**
   * @brief Compute inviscid flux along normal
   *
   * @param x conserved state at current integration point
   * @param normal normal vector, usually not a unit vector
   * @param Tr current element transformation with the integration point
   * @param fluxN inviscid flux dotted with normal
   * @return real_t maximum characteristic speed, c + |u.n|
   */
  real_t ComputeFluxDotN(const Vector &x, const Vector &normal,
                         FaceElementTransformations &Tr,
                         Vector &fluxN) const override;

    // Inviscid / Euler Flux
    template<typename GasT>
    MFEM_HOST_DEVICE inline static void
    ComputeInviscidFluxKernel(const GasT &gas,
                              const real_t *state,
                              real_t inv_flux[Prandtl::MAXEQ][Prandtl::MAXDIM])
    { 
      PointStateView S{state};
      
      // 1. Get states
      const int dim = gas.dim();
      const real_t density = gas.density(S);
      const real_t spec_vol = 1.0/density;
      real_t momentum[Prandtl::MAXDIM] = {0.,0.,0.};
      for(int idim = 0;idim < dim;idim++){
        momentum[idim] = gas.momentum(S, idim);
      }

      const real_t energy = gas.energy(S);
      const real_t pressure = gas.pressure(S);
      const real_t ke = gas.kinetic_energy_density(S);
      const int eq_mass = gas.L.eq_mass;
      const int eq_mom0 = gas.L.eq_mom0;
      const int eq_ener = gas.L.eq_energy;
      const int eq_spec = gas.L.eq_scalar0;

      const real_t H = (energy + pressure)*spec_vol;
      // 2. Compute Flux
      for (int d = 0; d < dim; d++)
        {
          inv_flux[eq_mass][d] = momentum[d];
          for (int i = 0; i < dim; i++)
            {
              // ρuuᵀ
              inv_flux[eq_mom0+i][d] = momentum[i]*momentum[d]*spec_vol;
            }
          // (ρuuᵀ) + p
          inv_flux[eq_mom0+d][d] += pressure;
          inv_flux[eq_ener][d] = momentum[d]*H;
          // for(int s = 0;s < gas.L.num_scalars;s++){
          //   inv_flux[eq_spec+s][d] = gas.scalar(S, s) * momentum[d] * spec_vol;
          // }
        }
      // 3. Compute maximum characteristic speed 
      // const real_t sound = gas.sound_speed(S);
      // fluid speed |u|
      // const real_t speed = Prandtl::Kernels::rsqrt(2.0 * ke / density);
      // max characteristic speed = fluid speed + sound speed
      // return speed + sound;
    }                        
};
  
}
