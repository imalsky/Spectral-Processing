/*----------------------- input.h ------------------------------------

Defines input values and files for 3-D emission spectra

---------------------------------------------------------------------- */

#ifndef __INPUT_H__
#define __INPUT_H__

/* I/O SETTINGS. */

/* File names */
#define OUTPUT_PREFIX <<output_file>>      /* output name */
#define T_P_3D_FILE <<input_file>>         /* input file */

/* Output settings */
#define N_PHASE 1                          /* Number of phases [96 max; lon grid in increments of 3.75] */
#define DOPPLER <<doppler>>                /* 0:Off; 1:On */
#define CLOUDS <<CLOUDS>>                           /* 0:Off; 1:On */

/* Grid settings */
#define NTAU <<NTAU>>                            /* Number of altitude points in grid      */
#define NLAT  <<NLAT>>                           /* Number of latitude points in 3-D  grid */
#define NLON  <<NLON>>                           /* Number of longitude points in 3-D grid */

#define NTEMP 30                           /* Number of temperature points in grid   */
#define NLAMBDA 2598                       /* Number of wavelength points in grid [4616/2598]   */

// This is the Npressure for low res
#define NPRESSURE <<num_pressure_points>>    /* Number of pressure points in grid   [13/17]   */

#define W0_VAL <<W0_VAL>>
#define G0_VAL <<G0_VAL>>


/* Planet parameters */
#define INPUT_INCLINATION <<inclination>>  /* Planet inclination in radians            */
#define INPUT_PHASE <<phase>>              /* Planet inclination in degrees           */
#define G 12.932                             /* Planet surface gravity                 */
#define R_PLANET 1.287e+08                 /* Planet radius at base of atmosphere      */

#define ORB_SEP 8.901e+09                  // This is some distance
#define STELLAR_TEMP 6213                // Stellar Blackbody temperature
#define R_STAR 1.029e+09                    /* Stellar radius                         */
#define P_ROT  4.6171                        /* Rotation period in days (= P_ORB for tidally locked planet)    */
#define R_VEL 0.0                          /* Radial Velocity                        */
#define MU 2.36                            /* Mean molecular weight                  */
#define FORMAT 2                           /* FORMAT=1 -> small opacity table        */
                                           /* FORMAT=2 -> large opacity table        */

/* Aerosol properties (calculated by the Mischenko Mie code) */


#define PI0_MgSiO3 0.99
#define G0_MgSiO3 0.14
#define QE_MgSiO3 0.07

#define PI0_Fe 0.71
#define G0_Fe -0.13
#define QE_Fe 1.25

#define PI0_Al2O3 0.74
#define G0_Al2O3 0.15
#define QE_Al2O3 0.12

#define PI0_MnS 0.99
#define G0_MnS 0.35
#define QE_MnS 0.56

/* Opacities for spectra */
#define CHEM_FILE   <<CHEM_FILE>>
#define CH4_FILE    <<CH4_FILE>>
#define CO2_FILE    <<CO2_FILE>>
#define CO_FILE     <<CO_FILE>>
#define H2O_FILE    <<H2O_FILE>>
#define NH3_FILE    <<NH3_FILE>>
#define O2_FILE     <<O2_FILE >>
#define O3_FILE     <<O3_FILE>>

#endif /* !__INPUT_H__ */

/* ------- end ---------------------------- input.h  ----------------- */