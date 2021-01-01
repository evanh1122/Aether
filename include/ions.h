// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef AETHER_INCLUDE_IONS_H_
#define AETHER_INCLUDE_IONS_H_

#include <armadillo>

#include "inputs.h"
#include "report.h"
#include "grid.h"

using namespace arma;

class Ions {

 public:

  struct species_chars {

    std::string cName;
    float mass;
    int charge;
    
    int DoAdvect;
    
    // Sources and Losses:

    fcube density_scgc;
    fcube par_velocity_vcgc;
    fcube perp_velocity_vcgc;

    fcube temperature_scgc;
    
    // Sources and Losses:

    fcube ionization_scgc;

    fcube sources_scgc;
    fcube losses_scgc;
    
  };

  // bulk quantities (states):
  fcube density_scgc;

  fcube ion_temperature_scgc;
  fcube electron_temperature_scgc;
  
  std::vector<species_chars> species;
  
  // ------------------------------
  // Functions:
  
  Ions(Grid grid, Inputs input, Report report);
  species_chars create_species(Grid grid);
  int read_planet_file(Inputs input, Report report);
  void fill_electrons(Report &report);

};
#endif // AETHER_INCLUDE_NEUTRALS_H_
