#ifndef FPC_H
#define FPC_H

/**********************************************************
FUNDAMENTAL PHYSICAL CONSTANTS
**********************************************************/

// Numbers
//=============================================================
const double PI = 3.141592653589793238;

// Fondamental constants
//=============================================================
const double AMU = 1.660531e-24; // atomic mass unit in grams
                                 // amu is 1/12 of mass of carbon atom

const double CCGS = 2.99792458e10; // light speed in cm/s
const double CSI = 2.99792458e8;   // light speed in m/s

const double GCGS = 6.67259e-8; // gravitational constant in cm3/g/s2
const double GSI = 6.67259e-11; // gravitational constant in m3/kg/s2

const double QCGS = 4.8032068e-10; // electron charge in ESU
const double QSI = 1.60217733e-19; // electron charge in SI

const double R0SI = 2.817940325e-15; // electron classical radius in SI

const double SIGMA = 5.67051e-5; // Stefan-Boltzmann constant in erg/cm3/K4

const double XKB = 1.380658e-16; // kB in erg/K
const double XKBSI = 1.380658e-23; // kB in J/K

const double XH = 6.6260755e-27; // h in erg/s

const double XME = 9.1093897e-28; // electron mass in g
const double XMP = 1.6726231e-24; // proton mass in g
const double XMH = 1.6605402e-24; // atomic mass unit in g
const double XMU = 1.6605402e-24; // atomic mass unit in g ( p214 Lang )
const double MP_EV = 0.9383e9; // proton mass in eV
const double MP_TEV = 0.9383e-3; // proton mass in eV

const double ThomsonCS_CM = 6.652e-25; // Thomson Cross section in cm2
const double PLANCK_CGS = 6.6260688e-27; // Planck constant in erg s 

// Conversion
//=============================================================
const double ARCMIN_TO_DEG = 1.0 / 60.0; // arcminute to degree
const double ARCSECOND_TO_RADIAN = 4.8481368e-6; // arcsecond to radian

const double CM_TO_KM = 1.e-5; // cm to km
const double CM_TO_M = 1.e-2; // cm to m
const double CM_TO_PC = 1.0 / 3.086e18; // cm to pc

const double DEG_TO_RAD = 1.7453e-2; // degree to radian

const double ERG_TO_EV  = 6.242e11;  // erg to eV
const double ERG_TO_KEV = 6.242e8;  // erg to keV
const double ERG_TO_MEV = 6.242e5;  // erg to MeV
const double ERG_TO_GEV = 6.242e2;  // erg to GeV
const double ERG_TO_TEV = 6.242e-1; // erg to TeV

const double ERG_TO_HZ  = 1.509e26;  // erg to Hz
const double ERG_TO_GHZ = 1.509e17; // erg to GHz

const double ESU_TO_C = 2.9979e9; // ESU to coulomb

const double GAUSS_TO_TESLA = 1e-4; // Gauss to Tesla (Tesla = SI unit)
const double GAUSS_TO_mGAUSS = 1.e3; // Gauss to microGauss
const double GAUSS_TO_uGAUSS = 1.e6; // Gauss to microGauss

const double GHZ_TO_ERG = 6.627e-18; // GHz to erg

const double EV_TO_ERG  = 1.602e-12; // eV to erg
const double KEV_TO_ERG = 1.602e-9; // keV to erg
const double MEV_TO_ERG = 1.602e-6; // MeV to erg
const double GEV_TO_ERG = 1.602e-3; // GeV to erg
const double TEV_TO_ERG = 1.602e0;  // TeV to erg

const double KEV_TO_HZ = 2.41774040e17; // keV to Hz

const double KG_TO_G = 1.e3; // kg to g

const double KM_TO_CM = 1.e5; // km to cm
const double KPC_TO_CM = 3.086e21; // kpc to cm

const double MBARN_TO_CM2 = 1.e-27; // mbarn to cm^2

const double M_TO_CM = 1.e2; // m to cm
const double M_TO_PC = 1.0 / 3.086e16; // m to pc

const double PC_TO_CM = 3.086e18; // pc to cm
const double PC_TO_KM = 3.086e13; // pc to km

const double S_TO_YR = 1.0 / 3.15576e7; // s to year

const double YR_TO_S = 3.15576e7; // year to s

// Units
//=============================================================
const double XJANSKY = 1.e-23; // in erg/s/cm2/Hz

// Main constants in astrophysics
//=============================================================
const double SUN_MASS = 1.989e33; // sun mass in g
const double SUN_RADIUS = 6.955e11; // sun radius in cm
const double SUN_LUMINOSITY = 3.8515e33; // sun luminosity in erg/s

const double EARTH_MASS = 5.976e27; // earth mass in g
const double EARTH_RADIUS = 6.378e8; // earth radius in cm


#endif

/*****************************************************************
COPYRIGHT:
   Gamil CASSAM-CHENAI ( gcc@cea.fr )
   Service d'Astrophysique, Bat 709
   CEA Saclay - Orme des Merisiers
   F-91191 Gif / Yvette - France
******************************************************************/
