/*
  Copyright (C) 2024 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/material_model/rheology/compositional_viscosity_prefactors.h>
#include <aspect/utilities.h>
#include <aspect/global.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/adiabatic_conditions/interface.h>


#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>
#include <aspect/simulator_signals.h>
#include <aspect/parameters.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      template <int dim>
      CompositionalViscosityPrefactors<dim>::CompositionalViscosityPrefactors ()
        = default;

      template <int dim>
      double
      CompositionalViscosityPrefactors<dim>::compute_viscosity (const MaterialModel::MaterialModelInputs<dim> &in,
                                                                const double base_viscosity,
                                                                const unsigned int composition_index,
                                                                const unsigned int q,
                                                                const ModifiedFlowLaws &modified_flow_laws) const
      {

        double factored_viscosities = base_viscosity;
        switch (viscosity_prefactor_scheme)
          {
            case none:
            {
              factored_viscosities = base_viscosity;
              break;
            }
            case hk04_olivine_hydration:
            {
              // We calculate the atomic H/Si ppm (C_OH) at each point to compute the water fugacity of
              // olivine assuming a composition of 90 mol% Forsterite and 10 mol% Fayalite from Hirth
              // and Kohlstedt 2004 10.1029/138GM06.
              const double temperature_for_fugacity = (this->simulator_is_past_initialization())
                                                      ?
                                                      in.temperature[q]
                                                      :
                                                      this->get_adiabatic_conditions().temperature(in.position[q]);
              AssertThrow(temperature_for_fugacity != 0, ExcMessage(
                            "The temperature used in the calculation for determining the water fugacity is zero. "
                            "This is not allowed, because this value is used to divide through. It is probably "
                            "being caused by the temperature being zero somewhere in the model. The relevant "
                            "values for debugging are: temperature (" + Utilities::to_string(in.temperature[q]) +
                            "), and adiabatic temperature ("
                            + Utilities::to_string(this->get_adiabatic_conditions().temperature(in.position[q])) +
                            "). If the adiabatic temperature is 0, double check that you are correctly defining an "
                            " 'Adiabatic conditions model'."));

              const unsigned int bound_fluid_idx = this->introspection().compositional_index_for_name("bound_fluid");
              const double mass_fraction_H2O = std::max(minimum_mass_fraction_water_for_dry_creep[composition_index], in.composition[q][bound_fluid_idx]); // mass fraction of bound water
              const double mass_fraction_olivine = 1 - mass_fraction_H2O; // mass fraction of olivine
              const double COH = (mass_fraction_H2O/molar_mass_H2O) / (mass_fraction_olivine/molar_mass_olivine) * 1e6; // COH in H / Si ppm
              const double point_water_fugacity = COH / A_H2O *
                                                  std::exp((activation_energy_H2O + in.pressure[q]*activation_volume_H2O)/
                                                           (constants::gas_constant * temperature_for_fugacity));
              const double r = modified_flow_laws == diffusion
                               ?
                               -diffusion_water_fugacity_exponents[composition_index]
                               :
                               -dislocation_water_fugacity_exponents[composition_index];
              factored_viscosities = base_viscosity*std::pow(point_water_fugacity, r);
              break;
            }

            case peng_robinson85_fugacity:
            {
              // Use the reference adiabatic pressure rather than the solved
              // pressure, which can contain a dynamic component that is not
              // part of the thermodynamic reference state used by this EOS.
              const double adiabatic_pressure = this->get_adiabatic_conditions().pressure(in.position[q]);

              // The parameters contain r/n, so the viscosity dependence is
              // f^(-r/n).
              const double viscosity_fugacity_exponent = modified_flow_laws == diffusion
                                                         ?
                                                         -diffusion_water_fugacity_exponents[composition_index]
                                                         :
                                                         -dislocation_water_fugacity_exponents[composition_index];

              const double point_water_fugacity =
                compute_fugacity(in.temperature[q], adiabatic_pressure);

              // Apply the raw fugacity in Pa without reference-fugacity
              // normalization. Consequently, creep-law prefactors must be
              // calibrated for fugacity expressed in Pa.
              factored_viscosities =
                base_viscosity
                * std::pow(point_water_fugacity, viscosity_fugacity_exponent);

              break;
            }

          }
        return factored_viscosities;
      }

      template <int dim>
      bool
      CompositionalViscosityPrefactors<dim>::uses_peng_robinson_fugacity () const
      {
        return viscosity_prefactor_scheme == peng_robinson85_fugacity;
      }



      template <int dim>
      double
      CompositionalViscosityPrefactors<dim>::compute_fugacity(
        const double temperature, const double pressure) const
      {
        AssertThrow(uses_peng_robinson_fugacity(),
                    ExcMessage("Water fugacity can only be evaluated when the "
                               "viscosity prefactor scheme is "
                               "'peng_robinson85_fugacity'."));

        // Use short names below to keep the Peng-Robinson equations close to
        // their conventional notation.
        const double T_c = critical_temperature;
        const double P_c = critical_pressure;
        const double R = constants::gas_constant;

        // Dimensionless constants in the original Peng-Robinson formulation.
        const double a_coefficient = 0.45724;
        const double b_coefficient = 0.07780;

        AssertThrow(temperature > 0.0 && pressure >= 0.0,
                    ExcMessage("Temperature must be positive and absolute "
                               "pressure must be non-negative."));

        // The zero-pressure limit is an ideal, infinitely dilute vapor.
        if (pressure == 0.0)
          return 0.0;

        // Cap the EOS pressure at 2.5 GPa; the Peng-Robinson model is not
        // calibrated for higher pressures.
        const double capped_pressure = std::min(pressure, 2.5e9);

        // With R in J/(mol K) and P_c in Pa, b is in m^3/mol and a is
        // in Pa (m^3/mol)^2.
        const double b = b_coefficient * R * T_c / P_c;
        const double alpha = std::pow(1.0 + kappa
                                      * (1.0 - std::sqrt(temperature/T_c)),
                                      2.0);
        const double a = a_coefficient * alpha * R*R*T_c*T_c/P_c;


        // Form the dimensionless Peng-Robinson coefficients used in the
        // compressibility-factor cubic.
        const double A = a * capped_pressure / (R*R*temperature*temperature);
        const double B = b * capped_pressure / (R*temperature);
        const double a0 = (-1)*A*B + B*B + B*B*B;
        const double a1 = A - 3*B*B - 2*B;
        const double a2 = B - 1;

        // Transform Z^3 + a2*Z^2 + a1*Z + a0 = 0 to the depressed
        // cubic y^3 + p*y + q = 0 and calculate all real roots.
        const double p =   a1 - a2*a2/3.0;
        const double q = 2.0*a2*a2*a2/27.0 - a2*a1/3.0 + a0;
        const double discriminant = q*q/4.0 + p*p*p/27.0;

        std::array<double,3> roots = {{0.0, 0.0, 0.0}};
        unsigned int n_roots = 0;

        if (discriminant > 0.0)
          {
            roots[0] = std::cbrt(-q/2.0 + std::sqrt(discriminant))
                       + std::cbrt(-q/2.0 - std::sqrt(discriminant))
                       - a2/3.0;
            n_roots = 1;
          }
        else
          {
            const double amplitude = 2.0*std::sqrt(std::max(0.0, -p/3.0));

            if (amplitude == 0.0)
              {
                roots[0] = -a2/3.0;
                n_roots = 1;
              }
            else
              {
                const double cosine_argument =
                  std::clamp(3.0*q/(p*amplitude), -1.0, 1.0);
                const double theta = std::acos(cosine_argument)/3.0;

                for (unsigned int i=0; i<3; ++i)
                  roots[i] = amplitude
                             * std::cos(theta + 2.0*numbers::PI*i/3.0)
                             - a2/3.0;
                n_roots = 3;
              }
          }


        // Retain physical roots and order them from liquid-like (small Z) to
        // vapor-like (large Z).
        std::vector<double> physical_roots;
        for (unsigned int i=0; i<n_roots; ++i)
          if (roots[i] > B)
            physical_roots.push_back(roots[i]);

        std::sort(physical_roots.begin(), physical_roots.end());
        physical_roots.erase(std::unique(physical_roots.begin(),
                                         physical_roots.end(),
                                         [](const double left, const double right)
        {
          return std::abs(left-right)
                 < 1e-12*std::max(1.0, std::abs(left));
        }),
        physical_roots.end());

        AssertThrow(!physical_roots.empty(),
                    ExcMessage("The Peng-Robinson equation did not produce "
                               "a physical compressibility-factor root."));

        const auto fugacity_for_root = [A, B, capped_pressure](const double Z)
        {
          const double sqrt_two = std::sqrt(2.0);
          const double logarithm_argument =
            (Z + (sqrt_two + 1.0)*B)
            / (Z - (sqrt_two - 1.0)*B);

          AssertThrow(Z > B && logarithm_argument > 0.0,
                      ExcMessage("A compressibility-factor root is not valid "
                                 "for calculating fugacity."));

          const double ln_phi =
            (Z - 1.0) - std::log(Z - B)
            - A/(2.0*sqrt_two*B)*std::log(logarithm_argument);
          return capped_pressure*std::exp(ln_phi);
        };

        // When more than one physical root exists, compare the liquid-like
        // and vapor-like roots and retain the phase with the lower fugacity.
        const double liquid_Z = physical_roots.front();
        const double vapor_Z = physical_roots.back();
        const double liquid_fugacity = fugacity_for_root(liquid_Z);
        const double vapor_fugacity = fugacity_for_root(vapor_Z);
        const double equilibrium_fugacity =
          std::min(liquid_fugacity, vapor_fugacity);

        AssertThrow(std::isfinite(equilibrium_fugacity),
                    ExcMessage("The Peng-Robinson equation produced a "
                               "non-finite fugacity."));

        return equilibrium_fugacity;
      }


      template <int dim>
      void
      CompositionalViscosityPrefactors<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Minimum mass fraction bound water content for fugacity", "6.15e-6",
                           Patterns::List(Patterns::Double(0)),
                           "The minimum water content for the HK04 olivine hydration viscosity "
                           "prefactor scheme. This acts as the cutoff between 'dry' creep and 'wet' creep "
                           "for olivine, and the default value is chosen based on the value reported by "
                           "Hirth & Kohlstedt 2004. For a mass fraction of bound water beneath this value, "
                           "this value is used instead to compute the water fugacity. Units: \\si{\\kg} / \\si{\\kg} %.");

        prm.declare_entry ("Water fugacity exponents for diffusion creep", "0.0",
                           Patterns::List(Patterns::Double(0)),
                           "List of water fugacity exponents for diffusion creep for "
                           "background material and compositional fields, for a total of N+1 "
                           "where N is the number of all compositional fields or only those "
                           "corresponding to chemical compositions. This is only applied when using the "
                           "Viscosity prefactor scheme 'HK04 olivine hydration' or "
                           "'peng_robinson85_fugacity'. Note, the water fugacity exponent "
                           "required by ASPECT for diffusion creep is r/n, where n is the stress exponent "
                           "for diffusion creep, which typically is 1. Units: none.");

        prm.declare_entry ("Water fugacity exponents for dislocation creep", "0.0",
                           Patterns::List(Patterns::Double(0)),
                           "List of water fugacity exponents for dislocation creep for "
                           "background material and compositional fields, for a total of N+1 "
                           "where N is the number of all compositional fields or only those "
                           "corresponding to chemical compositions. This is only applied when using the "
                           "Viscosity prefactor scheme 'HK04 olivine hydration' or "
                           "'peng_robinson85_fugacity'. Note, the water fugacity exponent "
                           "required by ASPECT for dislocation creep is r/n, where n is the stress exponent "
                           "for dislocation creep, which typically is 3.5. Units: none.");

        prm.declare_entry ("Viscosity prefactor scheme", "none",
                           Patterns::Selection("none|HK04 olivine hydration|peng_robinson85_fugacity"),
                           "Select what type of viscosity multiplicative prefactor scheme to apply. "
                           "Allowed entries are 'none', 'HK04 olivine hydration', and "
                           "'peng_robinson85_fugacity'. 'HK04 olivine hydration' calculates "
                           "the viscosity change due to hydrogen incorporation into olivine "
                           "following Hirth & Kohlstedt 2004 (10.1029/138GM06). "
                           "'peng_robinson85_fugacity' computes pure-water fugacity from "
                           "temperature and adiabatic pressure with the Peng-Robinson equation "
                           "of state and applies the resulting fugacity prefactor to dislocation "
                           "creep. 'none' does not modify the viscosity. Units: none.");
        prm.declare_entry ("Critical temperature", "647.3",
                           Patterns::Double (0.0),
                           "Critical temperature of the fluid used by the Peng-Robinson "
                           "equation of state. This parameter is only used when "
                           "'Viscosity prefactor scheme' is 'peng_robinson85_fugacity'. "
                           "Units: K.");
        prm.declare_entry ("Critical pressure", "22.12e6",
                           Patterns::Double (0.0),
                           "Critical pressure of the fluid used by the Peng-Robinson "
                           "equation of state. This parameter is only used when "
                           "'Viscosity prefactor scheme' is 'peng_robinson85_fugacity'. "
                           "Units: Pa.");
        prm.declare_entry ("Acentric factor", "0.344",
                           Patterns::Double (),
                           "Acentric factor of the fluid. The value is retained with the "
                           "fluid parameters, while the current equation-of-state calculation "
                           "uses the separately supplied 'Kappa' value in its temperature "
                           "correction. This parameter is only read for the "
                           "'peng_robinson85_fugacity' scheme. Units: none.");
        prm.declare_entry ("Kappa", "0.873236",
                           Patterns::Double (),
                           "Peng-Robinson kappa coefficient in the temperature-dependent "
                           "attraction parameter alpha. The default is the value corresponding "
                           "to water with an acentric factor of 0.344. This parameter is only "
                           "used by the 'peng_robinson85_fugacity' scheme. Units: none.");
      }



      template <int dim>
      void
      CompositionalViscosityPrefactors<dim>::parse_parameters (ParameterHandler &prm)
      {
        if (prm.get ("Viscosity prefactor scheme") == "none")
          viscosity_prefactor_scheme = none;
        if (prm.get ("Viscosity prefactor scheme") == "HK04 olivine hydration")
          {
            // Retrieve the list of compositional names
            std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();
            AssertThrow(this->introspection().compositional_name_exists("bound_fluid"),
                        ExcMessage("The HK04 olivine hydration pre-exponential factor only works if "
                                   "there is a compositional field called bound_fluid."));
            viscosity_prefactor_scheme = hk04_olivine_hydration;

            // Retrieve the list of chemical names
            std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();

            // Establish that a background field is required here
            compositional_field_names.insert(compositional_field_names.begin(), "background");
            chemical_field_names.insert(chemical_field_names.begin(),"background");

            Utilities::MapParsing::Options options(chemical_field_names, "Minimum mass fraction bound water content for fugacity");
            options.list_of_allowed_keys = compositional_field_names;
            minimum_mass_fraction_water_for_dry_creep = Utilities::MapParsing::parse_map_to_double_array (prm.get("Minimum mass fraction bound water content for fugacity"),
                                                        options);
          }
        if (prm.get ("Viscosity prefactor scheme") == "peng_robinson85_fugacity")
          {
            viscosity_prefactor_scheme = peng_robinson85_fugacity;

            // Read the thermodynamic constants only when this scheme is
            // selected. Other viscosity-prefactor schemes continue to use
            // their existing parameters and defaults.
            critical_temperature = prm.get_double ("Critical temperature");
            critical_pressure = prm.get_double ("Critical pressure");
            acentric_factor = prm.get_double ("Acentric factor");
            kappa = prm.get_double ("Kappa");

          }

        if (viscosity_prefactor_scheme == hk04_olivine_hydration
            || viscosity_prefactor_scheme == peng_robinson85_fugacity)
          {
            std::vector<std::string> compositional_field_names =
              this->introspection().get_composition_names();
            std::vector<std::string> chemical_field_names =
              this->introspection().chemical_composition_field_names();

            compositional_field_names.insert(compositional_field_names.begin(), "background");
            chemical_field_names.insert(chemical_field_names.begin(), "background");

            Utilities::MapParsing::Options options(
              chemical_field_names,
              "Water fugacity exponents for diffusion creep");
            options.list_of_allowed_keys = compositional_field_names;
            diffusion_water_fugacity_exponents =
              Utilities::MapParsing::parse_map_to_double_array(
                prm.get("Water fugacity exponents for diffusion creep"), options);

            options.property_name = "Water fugacity exponents for dislocation creep";
            dislocation_water_fugacity_exponents =
              Utilities::MapParsing::parse_map_to_double_array(
                prm.get("Water fugacity exponents for dislocation creep"), options);
          }
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
#define INSTANTIATE(dim) \
  template class CompositionalViscosityPrefactors<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }
}
