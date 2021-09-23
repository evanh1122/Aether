// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_CONSTANTS_H_
#define INCLUDE_CONSTANTS_H_

/*! \file constants.h
    \brief Physical consants.
    
    Contains constants used around the code.
*/

#include <vector>

// -------------------------------------------------------------------------
// Physical Constants
//   - Naming standards:
//     - Names start with a "c" to indicate they are constants
//     - Names should be all CAPS (except for "c" at start)
//     - Should use common abbreviations for names and not whole names
// -------------------------------------------------------------------------

// avogadros_number (per mole)
const float cNA = 6.02214086e23;

// boltzmanns_constant (m2 kg /s2 /K)
const float cKB = 1.38064852e-23;

// Universal Gas Constant (J/K/mol)
const float cR = cNA * cKB;

// Mass of a proton (kg)
const float cMP = 1.6726219e-27;

// Atomic Mass Unit (kk)
const float cAMU = cMP;

// Mass of an electron (kg)
const float cME = 9.1094e-31;

// Plancks Constant (Js)
const float cH = 6.6261e-34;

// Elemental Charge (C (J/eV))
const float cE = 1.6022e-19;

// Speed of Light (m/s)
const float cC = 2.9979e8;

// -------------------------------------------------------------------------
// Stellar constants:
// -------------------------------------------------------------------------

// gravitational_constant (m3/kg/s2)
const float cG = 6.67408e-11;

// This is for our sun (solar mass in kg):
const float cMSOL = 1.989e30;

// -------------------------------------------------------------------------
// Some constants for time
// -------------------------------------------------------------------------

const std::vector<int> cDAYS {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
const int cREFYEAR = 1965;
const double cREFJULIAN = 2438762.0;
const double cJULIAN2000 = 2451545.0;

// -------------------------------------------------------------------------
// Some Constants for Angles
// -------------------------------------------------------------------------

const float cPI = 3.141592653589793;
const float cTWOPI = 2*cPI;

// -------------------------------------------------------------------------
// Conversion Constants:
// - Naming standards:
//   - Names start with a "c" to indicate they are constants
//   - Names have a lower-case "to" in them to indicate conversion
//   - Names are all UPPER CASE otherwise
// -------------------------------------------------------------------------

const float cDtoR = cPI/180.0;
const float cRtoD = 180.0/cPI;

// -------------------------------------------------------------------------
// converting time between seconds and other units of time:
// -------------------------------------------------------------------------

// Years <-> Seconds:
const double cYtoS = 31536000.0;
const double cStoY = 1.0 / cYtoS;

// Days <-> Seconds:
const double cDtoS = 86400.0;
const double cStoD = 1.0 / cDtoS;

// Hours <-> Seconds:
const double cHtoS = 3600.0;
const double cStoH = 1.0 / cHtoS;

// Minutes <-> Seconds:
const double cMtoS = 60.0;
const double cStoM = 1.0 / cMtoS;

// MilliSeconds <-> Seconds:
const double cMStoS = 1.0/1000.0;
const double cStoMS = 1000.0;

// -------------------------------------------------------------------------
// converting lengths between meters and other units of measurement:
// -------------------------------------------------------------------------

// kilometer <-> meter:
const float cKMtoM = 1000.0;
const float cMtoKM = 1.0 / cKMtoM;

// Angstrom <-> Meter
const float cAtoM = 1.0e-10;

// AU <-> Meter
const float cAUtoM = 1.495978707e11;

const float pcm3topm3 = 1e6;  // /cm3 to /m3
const float pcm2topm2 = 1e4;  // /cm2 to /m2
const float pcmtopm = 100.0;  // /cm to /m

#endif  // INCLUDE_CONSTANTS_H_
