//header of MTCLM43.cpp

/* parameters for the Tair algorithm */
#define TDAYCOEF     0.45  /* (dim) daylight air temperature coefficient (dim) */

/* parameters for the snowpack algorithm */
#define SNOW_TCRIT   -6.0  /* (deg C) critical temperature for snowmelt */   
#define SNOW_TRATE  0.042  /* (cm/degC/day) snowmelt rate */

/* parameters for the radiation algorithm */
#define TBASE       0.870  /* (dim) max inst. trans., 0m, nadir, dry atm */
#define ABASE     -6.1e-5  /* (1/Pa) vapor pressure effect on transmittance */
#define C             1.5  /* (dim) radiation parameter */
#define B0          0.013  /* (dim) radiation parameter */
#define B1          0.201  /* (dim) radiation parameter */
#define B2          0.185  /* (dim) radiation parameter */
#define RAIN_SCALAR  0.75  /* (dim) correction to trans. for rain day 0.75 */
#define DIF_ALB       0.6  /* (dim) diffuse albedo for horizon correction */
#define SC_INT       1.32  /* (MJ/m2/day) snow correction intercept */
#define SC_SLOPE    0.096  /* (MJ/m2/day/cm) snow correction slope */

/* output file extension */
//#define POSTFIX  ".mtc43"  /* extension added to output filename prefix */

// constants
#define SECPERRAD 13750.9871     /* seconds per radian of hour angle */
#define RADPERDAY 0.017214       /* radians of Earth orbit per julian day */
#define RADPERDEG 0.01745329     /* radians per degree */
#define MINDECL -0.4092797       /* minimum declination (radians) */
#define DAYSOFF 11.25            /* julian day offset of winter solstice */
#define SRADDT 600.0             /* timestep for radiation routine (seconds) */

#define MA       28.9644e-3      /* (kg mol-1) molecular weight of air */
#define MW       18.0148e-3      /* (kg mol-1) molecular weight of water */
#define RGAS        8.3143          /* (m3 Pa mol-1 K-1) gas law constant */
#define G_STD    9.80665         /* (m s-2) standard gravitational accel. */ 
#define P_STD    101325.0        /* (Pa) standard pressure at 0.0 m elevation */
#define T_STD    288.15          /* (K) standard temp at 0.0 m elevation  */ 
#define CP       1010.0          /* (J kg-1 K-1) specific heat of air */
#define LR_STD   0.0065          /* (-K m-1) standard temperature lapse rate */
#define EPSILON      0.62196351      /* (MW/MA) unitless ratio of molec weights */
//#define PI       3.14159265      /* pi */



