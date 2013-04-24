//----------------------------------------------------------------------
//
//  FILE:   lpjio.cpp
//  AUTHOR: I. Ross
//  DATE:   26-MAY-2005
//  Modifed by Sarah 10-April-06 for updated version of LPJ-SPITFIRE
//  New version contains fire model in spin up & lightening input file
//
//  Modifed again by Doug, April/MAy-09 for including agriculture inputs
//  (crops and pasture) for SPITFIRE.
//----------------------------------------------------------------------
//
//  I/O file for standalone LPJ build using monthly climatology of
//  temperature, pressure and cloudiness from MOTIF format netCDF
//  files.
//
//----------------------------------------------------------------------

//======================================================================
//
//  HEADER FILES
//
//======================================================================

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <set>
#include <vector>
#include <cmath>
#include <limits>

using namespace std;

#include <netcdf.hh>
#include "Array.hh"


//======================================================================
//
//  CONSTANT DEFINITIONS
//
//======================================================================

// Array sizes etc.

const int NMON = 12, NDAYS = 365, NTRACE = 6, NTRANSIENT = 1;
const int NPFT = 9;
const int CLIMBUF = 20;
const int NCO2 = 3;                             // Number of C variables : C02, C13, C14
const int co2_dim_values[NCO2] = { 1, 2, 3 };  // Values for the C dimension

// Missing values to use for output files.

const float FMISSING = 1.0E35;
const int IMISSING = -999;


//======================================================================
//
//  TYPE DEFINITIONS
//
//======================================================================

//----------------------------------------------------------------------
//
//  Run parameters: read from configuration file.
//

struct RunParameters {

 RunParameters() :
    temp_var("tas"), prec_var("pr"), sun_var("clt"), wetdays_var("wetdays"),
    tmin_var("tasmin"), tmax_var("tasmax"), windsp_var("uvas"),
    mask_var("mask"), sun_is_cloud(true), mask_type(BOOLEAN),
    wetdays_fixed(false),
    soil_type_fixed(true), fixed_soil_type(2), soil_type_var("soiltype"),
    ann_popdens_var("popdens"),
    crop_var("crop"), crop_fixed(false),   //Doug 04/09
    pas_var("pas"), pas_fixed(false),      //Doug 04/09
    lightn_var("lightn"),    
    a_nd_fixed(false), fixed_a_nd(0.0), a_nd_var("a_nd"),
    spinup_years(1000), spinup_data_years(50),
    spinup_file("data_after_step1a.txt"), rampup_file("data_after_step1b.txt"),    //Doug 07/09
    rampup_years(0), rampup_data_years(0),  //Doug 07/09
    run_years(30),
    spinup_out_freq(0),
    rampup_out_freq(0),                        //Doug 07/09
    run_out_freq(1),
    run_out_years(0)                    //Doug 06/09: Parameter describing the year from when outputs start.
    ,co2_var("CO2"),c13_var("C13"),c14_var("C14")
    { }

  enum MASK_TYPE { BOOLEAN, PERCENT_LAND };


  // Driver data filenames and variable names (variable names default
  // to PMIP2 values).  The sun_is_cloud flag is used to indicate that
  // the input sun values are actually cloud.

  string temp_file, prec_file, sun_file, wetdays_file, mask_file;
  string tmin_file, tmax_file, windsp_file;
  string temp_var, prec_var, sun_var, wetdays_var, mask_var;
  string tmin_var, tmax_var, windsp_var;
  bool sun_is_cloud;
  MASK_TYPE mask_type;
  bool wetdays_fixed;
  float fixed_wetdays;


  // Soil type data - either a single fixed value or a filename for a
  // file on the same grid as the climate data.

  bool soil_type_fixed;
  int fixed_soil_type;
  string soil_type_file;
  string soil_type_var;

  // Population density and human fire ignition data.

  string ann_popdens_file;
  string ann_popdens_var;
  string crop_file,crop_var;    //Doug 04/09
  bool crop_fixed;              //Doug 04/09
  float fixed_crop;             //Doug 04/09
  string pas_file,pas_var;      //Doug 04/09
  bool pas_fixed;               //Doug 04/09
  float fixed_pas;              //Doug 04/09
  string lightn_file;
  string lightn_var;
  bool a_nd_fixed;
  float fixed_a_nd;
  string a_nd_file;
  string a_nd_var;

  // Number of years of spin-up to do (default to 1000), the number of
  // years of the climate data to use for the spinup (zero means use
  // all years), and the number of years of transient run to do (zero
  // means use all the available data once).

  int spinup_years, spinup_data_years, run_years;


  // Doug 07/09: Number of yearts of "ramping up" from end of the (often
  // detrended) spin up, to the point where you want to start doing
  // experiments (i.e, spin up in pre-industrial, ramp up to present day,
  // then main run for future experiments after present day). Zero ramp_up
  // litrally means no ramp up (i.e, for experiments that start at teh end
  // of spin up). zero rampup_data_years means use all available years in
  // input (provided rampup_years is not equal to zero). Zero
  // rampup_out_freq means no ouput during ramp_up.

  int rampup_years, rampup_data_years,  rampup_out_freq;


  // The interval in years between output records written during
  // spin-up (default is no output during spin-up, indicated by a zero
  // value) and the interval in years between output records written
  // during the main part of the run (output values are averaged over
  // these intervals, and the default, indicated by a zero value, is
  // to write a single averaged record at the end of the run).
  // Doug 06/09: run_out_years added. Year after-which ouputs occur.
  // initalized at 1 to allow all possible ouputs.

  int spinup_out_freq, run_out_freq, run_out_years;


  // Doug 07/09: Files and file location where spinup and rampup is I/Oed

  string spinup_file, rampup_file;


  // CO2 concentration.
  
  float co2_conc[NCO2];
  string co2_file, c13_file, c14_north_file, c14_equator_file, c14_south_file;
  string co2_var, c13_var, c14_var;

  // Output filename.
  string output_file;
  // Output directory.
  string output_dir;


  // Output variables.

  struct OutVariable {
    OutVariable() :
      name(""), mean(false), sdev(false), max(false), min(false) { }
    OutVariable(string n) :
      name(n), mean(true), sdev(false), max(false), min(false) { }
    string name;
    bool mean, sdev, max, min;
  };

  vector<OutVariable> output_vars;
};


//----------------------------------------------------------------------
//
//  Variable grid definition: derived from driver files.
//

struct Grid {
  int nlat, nlon;
  vector<float> lats, lons;

  Grid() : nlat(0), nlon(0) { }
  Grid(string file);
  Grid(const Grid &other);
  Grid &operator=(const Grid &other);

  bool operator==(const Grid &other);
  bool operator!=(const Grid &other) { return !(*this == other); }
};


//----------------------------------------------------------------------
//
//  Types for handling LPJ output variables.
//

enum LPJ_VAR_TYPE {
  SIMPLE,
  SIMPLE_MONTHLY,
  SIMPLE_DAILY,
  PER_PFT,
  PER_PFT_MONTHLY,
  PER_PFT_DAILY,
  TRACE,
  TRACE_MONTHLY,
  PER_CO2,
  PER_CO2_MONTHLY,
  //PER_CO2_DAILY,
  PER_CO2_PFT,
};

struct LPJVariable {
  const char *name;
  LPJ_VAR_TYPE type;
};

struct LPJOutputVariable {
  LPJOutputVariable(LPJVariable *in_def, RunParameters::OutVariable &collect);
  ~LPJOutputVariable();
  void set_nc(int idx, NcFile &nc);
  void reset(void) { clear_buffs(); out_index = 0; }
  void clear_buffs(void);
  bool average_and_output(void);

  void transpose_buffers(int time_axis_len, int other_axis_len);
  
  LPJVariable *def;
  NcVar **annual_nc, **mean_nc, **sdev_nc, **min_nc, **max_nc;
  int buff_size;
  float *annual_buff, *avg_buff, *sdev_buff, *min_buff, *max_buff;
  bool annual,mean, sdev, min, max;
  int accum_years;
  int out_index;
};


//======================================================================
//
//  LPJ VARIABLE DEFINITIONS
//
//======================================================================

// All the variable names we recognise from LPJ.
LPJVariable LPJ_VARIABLES[] = {
  { "dhuman_ign",    SIMPLE_DAILY    },   // Yan: 24/10/07
  { "dlightn",       SIMPLE_DAILY    },   // Yan: 24/10/07
  { "nind",          PER_PFT         },
  { "lm_ind",        PER_CO2_PFT     },
  { "rm_ind",        PER_CO2_PFT     },
  { "sm_ind",        PER_CO2_PFT     },
  { "hm_ind",        PER_CO2_PFT     },
  { "fpc_grid",      PER_PFT         },
  { "anpp_sum",      PER_CO2         },
  { "anpp",          PER_CO2_PFT     },   // Yan: 14/10/07
  { "anpp_pft_sum",  PER_CO2         },   // Yan: 14/10/07
  { "acflux_estab",  PER_CO2         },
  { "litter_ag",     PER_CO2_PFT     },
  { "litter_bg",     PER_CO2_PFT     },
  { "cpool_fast",    PER_CO2         },
  { "cpool_slow",    PER_CO2         },
  { "arh",           PER_CO2         },
  { "afire_frac",    SIMPLE          },
  { "afire_frac_human", SIMPLE       },  // Yan: 22/11/07
  { "afire_frac_lightn", SIMPLE      },  // Yan: 22/11/07
  { "acflux_fire",   PER_CO2         },
  { "mcflux_fire",   PER_CO2_MONTHLY },
  { "arunoff",       SIMPLE          },
  { "sla",           PER_PFT         },
  { "mpar",          SIMPLE_MONTHLY  },
  { "mapar",         SIMPLE_MONTHLY  },
  { "mphen",         PER_PFT_MONTHLY },
  { "anpp_add",      PER_CO2         },
  { "mnpp",          PER_PFT_MONTHLY },
  { "mnpp_grid",     SIMPLE_MONTHLY  }, //Doug 06/09
  { "mrunoff",       SIMPLE_MONTHLY  },
  { "aaet",          SIMPLE          },
  { "mrh",           PER_CO2_MONTHLY },
  { "mnpp_add",      PER_CO2_MONTHLY },
  { "maet",          SIMPLE_MONTHLY  },
  { "mgpp",          PER_PFT_MONTHLY },
  { "mtemp_soil",    SIMPLE_MONTHLY  },
  { "mcica",         PER_PFT_MONTHLY },
  { "lresp",         PER_PFT_MONTHLY },
  { "sresp",         PER_PFT_MONTHLY },
  { "rresp",         PER_PFT_MONTHLY },
  { "gresp",         PER_PFT_MONTHLY },
  { "aresp",         PER_PFT_MONTHLY },
  { "mw1",           SIMPLE_MONTHLY  },
  { "mw2",           SIMPLE_MONTHLY  },
  { "maep",          SIMPLE_MONTHLY  },
  { "aaep",          SIMPLE          },
  { "mpet_grid",     SIMPLE_MONTHLY  },
  { "apet_grid",     SIMPLE          },
  { "mintc",         SIMPLE_MONTHLY  },
  { "aintc",         SIMPLE          },
  { "meangc",        PER_PFT_MONTHLY },
  { "mgp",           PER_PFT_MONTHLY },
  { "num_fire",      SIMPLE_MONTHLY  },
  { "num_fire_human",      SIMPLE_MONTHLY  },  // Yan 22/11/07
  { "num_fire_lightn",     SIMPLE_MONTHLY  },  // Yan 22/11/07
  { "annum_fire",    SIMPLE          },
  { "annum_fire_human",    SIMPLE    },        // Yan 22/11/07
  { "annum_fire_lightn",   SIMPLE    },        // Yan 22/11/07
  { "area_burnt",    SIMPLE_MONTHLY  },        
  { "area_burnt_human",    SIMPLE_MONTHLY  },  // Yan 22/11/07
  { "area_burnt_lightn",   SIMPLE_MONTHLY  },  // Yan 22/11/07
  { "an_areafires",  SIMPLE          },
  { "an_areafires_human",  SIMPLE    },  // Yan 22/11/07
  { "an_areafires_lightn", SIMPLE    },  // Yan 22/11/07
  { "char_net_fuel_0", SIMPLE          },  
  { "livegrass_0", SIMPLE          },  
  { "dead_fuel_0", SIMPLE          },  
  { "dead_fuel_all_0", SIMPLE          },  
  { "fuel_all_0", SIMPLE          },  
  { "fuel_1hr_total_0", SIMPLE          },  
  { "fuel_10hr_total_0", SIMPLE          },  
  { "fuel_100hr_total_0", SIMPLE          },
  { "fuel_1000hr_total_0", SIMPLE          },
  { "mfuel_1hr_total",     SIMPLE_MONTHLY  }, //Doug 03/09
  { "mfuel_10hr_total",    SIMPLE_MONTHLY  }, //Doug 03/09
  { "mfuel_100hr_total",   SIMPLE_MONTHLY  }, //Doug 03/09
  { "mfuel_1000hr_total",  SIMPLE_MONTHLY  }, //Doug 03/09
  { "mlivegrass",  SIMPLE_MONTHLY  }, //Doug 05/09
//  { "fuel_1000hr_total_0", SIMPLE          },  
  { "deltaa", SIMPLE          },  
  { "deltaa_fpc", SIMPLE          },  
  { "delt_c13_fpc", SIMPLE          },  
  { "fbdep", SIMPLE          },  
  { "mfdi",          SIMPLE_MONTHLY  },
  { "mfire_frac",          SIMPLE_MONTHLY  }, 
  { "an_fdi",        SIMPLE          },
  { "litter_decom_ave",        SIMPLE          },
  { "turnover_ind",        SIMPLE          },
  { "an_fseason",    SIMPLE          },
  { "acflux_trace",  TRACE           },
  { "mcflux_trace",  TRACE_MONTHLY   },
  { "m_fc_crown",    SIMPLE_MONTHLY  },
  { "an_fc_crown",   SIMPLE          },
  { "m_i_surface",   SIMPLE_MONTHLY  },
  { "an_i_surface",  SIMPLE          },
  { "gdd",           PER_PFT         },
  { "height",        PER_PFT         },
//  { "mgpp",          PER_PFT_MONTHLY },
  { "wu",            SIMPLE          },
  { "wl",            SIMPLE          },
  { "dphen",         PER_PFT_DAILY   },
  { "dphen_change",  PER_PFT_MONTHLY },    //Doug 05/09
  { "fbare",         SIMPLE          },
  //sarah
//  { "soilc",           PER_CO2         },
  { "crop",          SIMPLE          },    //Doug 05/09: proportion of cell with crops
  { "pas",           SIMPLE          },    //Doug 05/09: proportion of cell with pasture
  { "anpp_grid",     PER_PFT         },    //Doug 07/09: cheat annual npp without the deltaC's
  { "arh_grid",      SIMPLE          },    //Doug 07/09: cheat again, for hetro. resp.
  { "acflux_fire_grid", SIMPLE       },    //Doug 07/09: another cherat for carbon burnt
  { "gdd_grid",      SIMPLE          },    //Doug 07/09: bioclimatic varible growing degree days base 5
  { "alpha_ws",      SIMPLE          },    //Doug 07/09: Bioclimatic varible for water stress
  { "cgf",  SIMPLE_MONTHLY },               //Doug 12/10: cload-to-ground fraction
  { "fdry",  SIMPLE_MONTHLY },               //Doug 12/10: cload-to-ground fraction
  { "lt_days",  SIMPLE_MONTHLY }               //Doug 12/10: cload-to-ground fraction
};



//======================================================================
//
//  GLOBAL VARIABLE DEFINITIONS
//
//======================================================================

RunParameters params;           // Run parameters.

Grid grid;                      // Common variable grid definition.

NcFile *clim_temp_nc;           // Climate data files.
NcFile *clim_precip_nc;
NcFile *clim_sun_nc;
NcFile *clim_wetdays_nc;
NcFile *clim_tmin_nc;
NcFile *clim_tmax_nc;
NcFile *clim_windsp_nc;
NcFile *clim_lightn_nc;
NcFile *clim_ann_popdens_nc;
NcFile *clim_crop_nc;            // Doug 04/09
NcFile *clim_pas_nc;             // Doug 04/09

NcVar *clim_temp_var;           // Climate data variables.
NcVar *clim_precip_var;
NcVar *clim_sun_var;
NcVar *clim_wetdays_var;
NcVar *clim_tmin_var;
NcVar *clim_tmax_var;
NcVar *clim_windsp_var;
NcVar *clim_lightn_var;
NcVar *clim_ann_popdens_var;
NcVar *clim_crop_var;            // Doug 04/09
NcVar *clim_pas_var;             // Doug 04/09

BoolArray *land_mask;           // Land mask.

IntArray *soil_type;            // Soil type data.

Array *a_nd_vals;               // Human fire ignition data.

int no_of_months;               // Climate data size (in months).

float *clim_temp;               // Climate data for current point.
float *clim_precip;
float *clim_sun;
float *clim_wetdays;
float *clim_tmin;
float *clim_tmax;
float *clim_windsp;
float *clim_lightn;
float *clim_ann_popdens;
float *clim_crop;                // Doug 04/09
float *clim_pas;                 // Doug 04/09

int land_points, points_done;   // Point counters.
int latidx, lonidx;

NcFile **out_nc;                // Output files.
int out_year_count;
int out_idx;
int *out_years;

NcDim *lat_dim;                 // Output dimensions.
NcDim *lon_dim;
NcDim *pft_dim;
NcDim *month_dim;               // Month index.
bool month_dim_needed = false;
NcDim *day_dim;                 // Day-of-year index.
bool day_dim_needed = false;
NcDim *trace_dim;               // Trace gas index.
bool trace_dim_needed = false;
NcDim *co2_dim;
bool co2_dim_needed = false;

LPJOutputVariable **out_vars = 0; // Output variables.
int out_var_count = 0;

float * co2_vals;
float * c13_vals;
float * c14_north_vals;
float * c14_equator_vals;
float * c14_south_vals;


//======================================================================
//
//  LOCAL FUNCTION PROTOTYPES
//
//======================================================================

static void read_config(RunParameters &params);
static bool boolean_string(string value, int line, bool &error_flag);
static void read_mask(BoolArray *mask, string mask_file, Grid &grid,
                      string var, RunParameters::MASK_TYPE mask_type);
static void read_soil_type(IntArray *soil_type, string soil_type_file,
                           Grid &grid, string soil_type_var);
static void read_a_nd(Array *a_nd, string a_nd_file,
                      Grid &grid, string a_nd_var);
static void check_var_grid(NcVar *var, NcFile *nc,
                           string var_name, string file_name);
static void open_output_files(string file_name, string directory_name);
static void handle_output_record(string name, float *data);
static int read_co2(string var_name, string file_name, float **data);

//======================================================================
//
//  LPJ INTERFACE FUNCTION DEFINITIONS
//
//======================================================================

//--------------------------------------------------------------------
//
//  initio
//
//  Reads input climate and land mask files, initialises for looping
//  over the grid cells and time.
//   

// In fortran : real dummy_array(1:2,1:3,1:4)
typedef float (*DummyFortranArray)[3][2];

extern "C" void initio_(DummyFortranArray dummy_array)
{

  // Check array order
  float dummy = 1.0;

  for( int i = 0; i < 2; ++i ) {
    for( int j = 0; j < 3; ++j ) {
      for( int k = 0; k < 4; ++k ) {
        if(  dummy_array[k][j][i] != dummy ) {
          cerr << "FATAL ERROR : layout of fortran arrays in C++ is not what was expected\n";
          exit( -1 );
        }
        // As a float[]
        float * d = (float *)dummy_array;
        if( d[(( k * 3 ) + j ) * 2 + i] != dummy ) {
          cerr << "FATAL ERROR : layout of fortran arrays in C++ is not what was expected\n";
          exit( -1 );
        }
        dummy += 1;
      }
    }

  }

  // Read the configuration file.
  read_config(params);

  // Extract grid information from input files.
  Grid mask_grid(params.mask_file);
  Grid temp_grid(params.temp_file);
  Grid prec_grid(params.prec_file);
  Grid sun_grid(params.sun_file);
  Grid tmin_grid(params.tmin_file);
  Grid tmax_grid(params.tmax_file);
  Grid windsp_grid(params.windsp_file);
  Grid popdens_grid(params.ann_popdens_file);
  Grid lightn_grid(params.lightn_file);
  if (mask_grid != temp_grid || prec_grid != sun_grid ||
      mask_grid != prec_grid || mask_grid != tmin_grid ||
      mask_grid != tmax_grid || mask_grid != windsp_grid ||
      mask_grid != popdens_grid || mask_grid != lightn_grid )
  {
    cerr << "ERROR:c Incompatible grid definitions in driver files" << endl;
    exit(1);
  }

  grid = mask_grid;
  if (!params.wetdays_fixed) {
    Grid wetdays_grid(params.wetdays_file);
    if (grid != wetdays_grid) {
      cerr << "ERROR:wd Incompatible grid definitions in driver file wetdays" << endl;
      exit(1);
    }
  }

  //Doug 04/09: crops and pasture annual fraction
  if (!params.crop_fixed) {
    Grid crop_grid(params.crop_file);
    if (grid != crop_grid) {
      cerr << "ERROR:wd Incompatible grid definitions in driver file crops" << endl;
      exit(1);
    }
  }

  if (!params.pas_fixed) {
    Grid pas_grid(params.pas_file);
    if (grid != pas_grid) {
      cerr << "ERROR:wd Incompatible grid definitions in driver file pasture" << endl;
      exit(1);
    }
  }
  if (!params.soil_type_fixed) {
    Grid soil_type_grid(params.soil_type_file);
    if (soil_type_grid != grid) {
      cerr << "ERROR:st Incompatible grid definitions in driver file soil type" << endl;
      exit(1);
    }
  }

  if (!params.a_nd_fixed) {
    Grid a_nd_grid(params.a_nd_file);
    if (a_nd_grid != grid) {
      cerr << "ERROR: Incompatible grid definitions in driver files a(Nd)" << endl;
      exit(1);
    }
  }

  // Read the land/sea mask and convert to a boolean array.
  land_mask = new BoolArray(grid.nlon, grid.nlat);
  read_mask(land_mask, params.mask_file, grid,
            params.mask_var, params.mask_type);

  // Read the soil type data, if necessary.
  if (!params.soil_type_fixed) {
    soil_type = new IntArray(grid.nlon, grid.nlat);
    read_soil_type(soil_type, params.soil_type_file,
                   grid, params.soil_type_var);
  }

  // Read the population density data,human ignition risk and lightening
  // 
  if (!params.a_nd_fixed) {
    a_nd_vals = new Array(grid.nlon, grid.nlat);
    read_a_nd(a_nd_vals, params.a_nd_file, grid, params.a_nd_var);
  }
  // Open all climate input data files.

  clim_temp_nc = new NcFile(params.temp_file.c_str());
  if (!clim_temp_nc->is_valid()) {
    cerr << "ERROR: could not open climate file " << params.temp_file << endl;
    exit(1);
  }
  if (params.prec_file == params.temp_file)
    clim_precip_nc = clim_temp_nc;
  else {
    clim_precip_nc = new NcFile(params.prec_file.c_str());
    if (!clim_precip_nc->is_valid()) {
      cerr << "ERROR: could not open climate file "
           << params.prec_file << endl;
      exit(1);
    }
  }
  if (params.sun_file == params.temp_file)
    clim_sun_nc = clim_temp_nc;
  else if (params.sun_file == params.prec_file)
    clim_sun_nc = clim_precip_nc;
  else {
    clim_sun_nc = new NcFile(params.sun_file.c_str());
    if (!clim_sun_nc->is_valid()) {
      cerr << "ERROR: could not open climate file "
           << params.sun_file << endl;
      exit(1);
    }
  }

  if (!params.wetdays_fixed) {
    if (params.wetdays_file == params.temp_file)
      clim_wetdays_nc = clim_temp_nc;
    else if (params.wetdays_file == params.prec_file)
      clim_wetdays_nc = clim_precip_nc;
    else if (params.wetdays_file == params.sun_file)
      clim_wetdays_nc = clim_sun_nc;
    else {
      clim_wetdays_nc = new NcFile(params.wetdays_file.c_str());
      if (!clim_wetdays_nc->is_valid()) {
        cerr << "ERROR: could not open climate file "
             << params.wetdays_file << endl;
        exit(1);
      }
    }
  }

  if (params.tmin_file == params.temp_file)
    clim_tmin_nc = clim_temp_nc;
  else if (params.tmin_file == params.prec_file)
    clim_tmin_nc = clim_precip_nc;
  else if (params.tmin_file == params.sun_file)
    clim_tmin_nc = clim_sun_nc;
  else if (!params.wetdays_fixed && params.tmin_file == params.wetdays_file)
    clim_tmin_nc = clim_wetdays_nc;
  else {
    clim_tmin_nc = new NcFile(params.tmin_file.c_str());
    if (!clim_tmin_nc->is_valid()) {
      cerr << "ERROR: could not open climate file "
           << params.tmin_file << endl;
      exit(1);
    }
  }

  if (params.tmax_file == params.temp_file)
    clim_tmax_nc = clim_temp_nc;
  else if (params.tmax_file == params.prec_file)
    clim_tmax_nc = clim_precip_nc;
  else if (params.tmax_file == params.sun_file)
    clim_tmax_nc = clim_sun_nc;
  else if (!params.wetdays_fixed && params.tmax_file == params.wetdays_file)
    clim_tmax_nc = clim_wetdays_nc;
  else if (params.tmax_file == params.tmin_file)
    clim_tmax_nc = clim_tmin_nc;
  else {
    clim_tmax_nc = new NcFile(params.tmax_file.c_str());
    if (!clim_tmax_nc->is_valid()) {
      cerr << "ERROR: could not open climate file "
           << params.tmax_file << endl;
      exit(1);
    }
  }

  if (params.windsp_file == params.temp_file)
    clim_windsp_nc = clim_temp_nc;
  else if (params.windsp_file == params.prec_file)
    clim_windsp_nc = clim_precip_nc;
  else if (params.windsp_file == params.sun_file)
    clim_windsp_nc = clim_sun_nc;
  else if (!params.wetdays_fixed && params.windsp_file == params.wetdays_file)
    clim_windsp_nc = clim_wetdays_nc;
  else if (params.windsp_file == params.tmin_file)
    clim_windsp_nc = clim_tmin_nc;
  else if (params.windsp_file == params.tmax_file)
    clim_windsp_nc = clim_tmax_nc;
  else {
    clim_windsp_nc = new NcFile(params.windsp_file.c_str());
    if (!clim_windsp_nc->is_valid()) {
      cerr << "ERROR: could not open climate file "
           << params.windsp_file << endl;
      exit(1);
    }
  }

  //Kirsten: 12-month lightning climatology
  if (params.lightn_file == params.temp_file)
    clim_lightn_nc = clim_temp_nc;
  else if (params.lightn_file == params.prec_file)
    clim_lightn_nc = clim_precip_nc;
  else if (params.lightn_file == params.sun_file)
    clim_lightn_nc = clim_sun_nc;
  else if (!params.wetdays_fixed && params.lightn_file == params.wetdays_file)
    clim_lightn_nc = clim_wetdays_nc;
  else if (params.lightn_file == params.tmin_file)
    clim_lightn_nc = clim_tmin_nc;
  else if (params.lightn_file == params.tmax_file)
    clim_lightn_nc = clim_tmax_nc;
  else if (params.lightn_file == params.windsp_file)
    clim_lightn_nc = clim_windsp_nc;
  else {
    clim_lightn_nc = new NcFile(params.lightn_file.c_str());
    if (!clim_lightn_nc->is_valid()) {
      cerr << "ERROR: could not open climate file "
           << params.lightn_file << endl;
      exit(1);
    }
  }

  //Kirsten: annual population density
  if (params.ann_popdens_file == params.temp_file)
    clim_ann_popdens_nc = clim_temp_nc;
  else if (params.ann_popdens_file == params.prec_file)
    clim_ann_popdens_nc = clim_precip_nc;
  else if (params.ann_popdens_file == params.sun_file)
    clim_ann_popdens_nc = clim_sun_nc;
  else if (!params.wetdays_fixed && params.ann_popdens_file == params.wetdays_file)
    clim_ann_popdens_nc = clim_wetdays_nc;
  else if (params.ann_popdens_file == params.tmin_file)
    clim_ann_popdens_nc = clim_tmin_nc;
  else if (params.ann_popdens_file == params.tmax_file)
    clim_ann_popdens_nc = clim_tmax_nc;
  else if (params.ann_popdens_file == params.windsp_file)
    clim_ann_popdens_nc = clim_windsp_nc;
  else {
    clim_ann_popdens_nc = new NcFile(params.ann_popdens_file.c_str());
    if (!clim_ann_popdens_nc->is_valid()) {
      cerr << "ERROR: could not open climate file "
           << params.ann_popdens_file << endl;
      exit(1);
    }
  }

  //Doug 04/09: crops and pasture fraction
  if (!params.crop_fixed) {
    if (params.crop_file == params.temp_file)
      clim_crop_nc = clim_temp_nc;
    else if (params.crop_file == params.prec_file)
      clim_crop_nc = clim_precip_nc;
    else if (params.crop_file == params.sun_file)
      clim_crop_nc = clim_sun_nc;
    else if (!params.wetdays_fixed && params.crop_file == params.wetdays_file)
      clim_crop_nc = clim_wetdays_nc;
    else if (params.crop_file == params.tmin_file)
      clim_crop_nc = clim_tmin_nc;
    else if (params.crop_file == params.tmax_file)
      clim_crop_nc = clim_tmax_nc;
    else if (params.crop_file == params.windsp_file)
      clim_crop_nc = clim_windsp_nc;
    else if (params.crop_file == params.ann_popdens_file)
      clim_crop_nc = clim_ann_popdens_nc; 
    else {
      clim_crop_nc = new NcFile(params.crop_file.c_str());
      if (!clim_crop_nc->is_valid()) {
        cerr << "ERROR: could not open climate file "
             << params.crop_file << endl;
        exit(1);
      }
    }
  }

  if (!params.pas_fixed) {
    if (params.pas_file == params.temp_file)
      clim_pas_nc = clim_temp_nc;
    else if (params.pas_file == params.prec_file)
      clim_pas_nc = clim_precip_nc;
    else if (params.pas_file == params.sun_file)
      clim_pas_nc = clim_sun_nc;
    else if (!params.wetdays_fixed && params.pas_file == params.wetdays_file)
      clim_pas_nc = clim_wetdays_nc;
    else if (params.pas_file == params.tmin_file)
      clim_pas_nc = clim_tmin_nc;
    else if (params.pas_file == params.tmax_file)
      clim_pas_nc = clim_tmax_nc;
    else if (params.pas_file == params.windsp_file)
      clim_pas_nc = clim_windsp_nc;
    else if (params.pas_file == params.ann_popdens_file)
      clim_pas_nc = clim_ann_popdens_nc;
    else if (!params.crop_fixed && params.pas_file == params.crop_file) 
      clim_pas_nc = clim_wetdays_nc;
    else {
      clim_pas_nc = new NcFile(params.pas_file.c_str());
      if (!clim_pas_nc->is_valid()) {
        cerr << "ERROR: could not open climate file "
             << params.pas_file << endl;
        exit(1);
      }
    }
  }

   
  // Access all climate variables.

  clim_temp_var = clim_temp_nc->get_var(params.temp_var.c_str());
  if (!clim_temp_var) {
    cerr << "ERROR: could not access climate variable '"
         << params.temp_var << "'" << endl;
    exit(1);
  }

  check_var_grid(clim_temp_var, clim_temp_nc,
                 params.temp_var, params.temp_file);
  clim_precip_var = clim_precip_nc->get_var(params.prec_var.c_str());
  if (!clim_precip_var) {
    cerr << "ERROR: could not access climate variable '"
         << params.prec_var << "'" << endl;
    exit(1);
  }

  check_var_grid(clim_precip_var, clim_precip_nc,
                 params.prec_var, params.prec_file);
  clim_sun_var = clim_sun_nc->get_var(params.sun_var.c_str());
  if (!clim_sun_var) {
    cerr << "ERROR: could not access climate variable '"
         << params.sun_var << "'" << endl;
    exit(1);
  }

  check_var_grid(clim_sun_var, clim_sun_nc,
                 params.sun_var, params.sun_file);
  if (!params.wetdays_fixed) {
    clim_wetdays_var = clim_wetdays_nc->get_var(params.wetdays_var.c_str());
    if (!clim_wetdays_var) {
      cerr << "ERROR: could not access climate variable '"
           << params.wetdays_var << "'" << endl;
      exit(1);
    }
    check_var_grid(clim_wetdays_var, clim_wetdays_nc,
                   params.wetdays_var, params.wetdays_file);
  }

  clim_tmin_var = clim_tmin_nc->get_var(params.tmin_var.c_str());
  if (!clim_tmin_var) {
    cerr << "ERROR: could not access climate variable '"
         << params.tmin_var << "'" << endl;
    exit(1);
  }
  check_var_grid(clim_tmin_var, clim_tmin_nc,
                 params.tmin_var, params.tmin_file);

  clim_tmax_var = clim_tmax_nc->get_var(params.tmax_var.c_str());
  if (!clim_tmax_var) {
    cerr << "ERROR: could not access climate variable '"
         << params.tmax_var << "'" << endl;
    exit(1);
  }
  check_var_grid(clim_tmax_var, clim_tmax_nc,
                 params.tmax_var, params.tmax_file);

  clim_windsp_var = clim_windsp_nc->get_var(params.windsp_var.c_str());
  if (!clim_windsp_var) {
    cerr << "ERROR: could not access climate variable '"
         << params.windsp_var << "'" << endl;
    exit(1);
  }
  check_var_grid(clim_windsp_var, clim_windsp_nc,
                 params.windsp_var, params.windsp_file);

  clim_lightn_var = clim_lightn_nc->get_var(params.lightn_var.c_str());
  if (!clim_lightn_var) {
    cerr << "ERROR: could not access climate variable '"
         << params.lightn_var << "'" << endl;
    exit(1);
  }
  check_var_grid(clim_lightn_var, clim_lightn_nc,
                 params.lightn_var, params.lightn_file);

  clim_ann_popdens_var = clim_ann_popdens_nc->get_var(params.ann_popdens_var.c_str());
  if (!clim_ann_popdens_var) {
    cerr << "ERROR: could not access climate variable '"
         << params.ann_popdens_var << "'" << endl;
    exit(1);
  }

  check_var_grid(clim_ann_popdens_var, clim_ann_popdens_nc,
                 params.ann_popdens_var, params.ann_popdens_file);

  //Doug 04/09: Annual crops and pasture
  if (!params.crop_fixed) {
    clim_crop_var = clim_crop_nc->get_var(params.crop_var.c_str());
    if (!clim_crop_var) {
      cerr << "ERROR: could not access climate variable '"
           << params.crop_var << "'" << endl;
      exit(1);
    }
  }

  check_var_grid(clim_crop_var, clim_crop_nc,
                 params.crop_var, params.crop_file);

  if (!params.pas_fixed) {
    clim_pas_var = clim_pas_nc->get_var(params.pas_var.c_str());
    if (!clim_pas_var) {
      cerr << "ERROR: could not access climate variable '"
           << params.pas_var << "'" << endl;
      exit(1);
    }
  }

  check_var_grid(clim_pas_var, clim_pas_nc,
                 params.pas_var, params.pas_file);


  // Check the time bounds on the climate variables.  There should be
  // a whole number of years of data, and there should be the same
  // amount of data for each variable.

  NcDim *temp_time_dim = clim_temp_var->get_dim(0);
  NcDim *precip_time_dim = clim_precip_var->get_dim(0);
  NcDim *sun_time_dim = clim_sun_var->get_dim(0);
  NcDim *tmin_time_dim = clim_tmin_var->get_dim(0);
  NcDim *tmax_time_dim = clim_tmax_var->get_dim(0);
  NcDim *windsp_time_dim = clim_windsp_var->get_dim(0);
  NcDim *wetdays_time_dim = clim_wetdays_var->get_dim(0);

// Yan: put windsp the same as other climate varaibles
// if (windsp_time_dim->size() != NMON) {
//   cerr << "ERROR: bad time dimension in windspeed file" << endl;
//   exit(1);
// }

  NcDim *lightn_time_dim = clim_lightn_var->get_dim(0);
  if (lightn_time_dim->size() != NMON) {
    cerr << "ERROR: bad time dimension in lightn file" << endl;
    exit(1);
  }

  NcDim *ann_popdens_time_dim = clim_ann_popdens_var->get_dim(0);
  NcDim *crop_time_dim = clim_crop_var->get_dim(0);  //Doug 04/09
  NcDim *pas_time_dim = clim_pas_var->get_dim(0);    //Doug 04/09

  if (temp_time_dim->size() != precip_time_dim->size() ||
      temp_time_dim->size() != sun_time_dim->size() ||
      temp_time_dim->size() != windsp_time_dim->size() ||
      temp_time_dim->size() != wetdays_time_dim->size() ||
      temp_time_dim->size() != tmin_time_dim->size() ||
      temp_time_dim->size() != tmax_time_dim->size()
      || temp_time_dim->size() != ann_popdens_time_dim ->size() * NMON
      || temp_time_dim->size() != crop_time_dim->size()*NMON    //Doug 04/09
      || temp_time_dim->size() != pas_time_dim->size()*NMON    //Doug 04/09
      ) {
    cerr << "ERROR: incompatible time dimensions in popdens climate files" << endl;
    exit(1);
  }

//  if (!params.wetdays_fixed) {
//    NcDim *wetdays_time_dim = clim_wetdays_var->get_dim(0);
//    if (wetdays_time_dim->size() != NMON) {
//      cerr << "ERROR: bad time dimension in wet days file" << endl;
//      exit(1);
//    }
//  }

  no_of_months = temp_time_dim->size();
  if (no_of_months % NMON != 0) {
    cerr << "ERROR: non whole number of years of climate data!" << endl;
    exit(1);
  }

  clim_temp = new float[no_of_months];
  clim_precip = new float[no_of_months];
  clim_sun = new float[no_of_months];
  clim_wetdays = new float[no_of_months];
  clim_tmin = new float[no_of_months];
  clim_tmax = new float[no_of_months];
  clim_windsp = new float[no_of_months];
  clim_lightn = new float[NMON];
  clim_ann_popdens = new float[no_of_months / NMON];
  clim_crop = new float[no_of_months / NMON];   //Doug 04/09
  clim_pas = new float[no_of_months / NMON];   //Doug 04/09

  if( no_of_months / NMON != read_co2( params.co2_var, params.co2_file, &co2_vals )) {
      cerr << "ERROR: incompatible time dimensions in co2 files" << endl;
      exit(1);
  }
  if( no_of_months / NMON != read_co2( params.c13_var, params.c13_file, &c13_vals )) {
      cerr << "ERROR: incompatible time dimensions in co2 files" << endl;
      exit(1);
  }
  if( no_of_months / NMON != read_co2( params.c14_var, params.c14_north_file, &c14_north_vals )) {
      cerr << "ERROR: incompatible time dimensions in co2 files" << endl;
      exit(1);
  }
  if( no_of_months / NMON != read_co2( params.c14_var, params.c14_equator_file, &c14_equator_vals )) {
      cerr << "ERROR: incompatible time dimensions in co2 files" << endl;
      exit(1);
  }
  if( no_of_months / NMON != read_co2( params.c14_var, params.c14_south_file, &c14_south_vals )) {
      cerr << "ERROR: incompatible time dimensions in co2 files" << endl;
      exit(1);
  }

  // Count land points.

  land_points = 0;
  for (int lat = 0; lat < grid.nlat; ++lat)
    for (int lon = 0; lon < grid.nlon; ++lon)
      if ((*land_mask)(lon, lat)) ++land_points;

  // Determine number of output files.
  out_year_count = 0;

#ifdef LPJ_STEP_1A
  if (params.spinup_out_freq != 0)
    out_year_count += params.spinup_years / params.spinup_out_freq;

#elif defined(LPJ_STEP_1B)
  if (params.rampup_years !=0) {                                         //Doug 07/09
    if (params.rampup_out_freq != 0)                                     //Doug 07/09
      out_year_count += params.rampup_years / params.rampup_out_freq;    //Doug 07/09
  }
#elif defined(LPJ_STEP_2)
  if (params.run_out_freq == 0)
    out_year_count += 1;
  else {
    if (params.run_years == 0)
      out_year_count += ((no_of_months / NMON)-params.run_out_years) / params.run_out_freq;
    else
      out_year_count += (params.run_years-params.run_out_years) / params.run_out_freq;
                        //Doug 06/09: add -parasm.run_out_years to reduce number of
                        //outputs by years before ouput start
  }

#else
  if (params.spinup_out_freq != 0)
    out_year_count += params.spinup_years / params.spinup_out_freq;

  if (params.rampup_years !=0) {                                         //Doug 07/09
    if (params.rampup_out_freq != 0)                                     //Doug 07/09
      out_year_count += params.rampup_years / params.rampup_out_freq;    //Doug 07/09
  }

  if (params.run_out_freq == 0)
    out_year_count += 1;
  else {
    if (params.run_years == 0)
      out_year_count += ((no_of_months / NMON)-params.run_out_years) / params.run_out_freq;
    else
      out_year_count += (params.run_years-params.run_out_years) / params.run_out_freq;
                        //Doug 06/09: add -parasm.run_out_years to reduce number of
                        //outputs by years before ouput start
  }
#endif

    cerr << out_year_count;
  out_years = new int[out_year_count];
  int out_year_idx = 0;

#ifdef LPJ_STEP_1A
  if (params.spinup_out_freq != 0) {
    for (int idx = 1;
         idx <= params.spinup_years / params.spinup_out_freq; ++idx)
      out_years[out_year_idx++] = idx * params.spinup_out_freq;
  }

#elif defined(LPJ_STEP_1B)
  if (params.rampup_years !=0 ) {
    if (params.rampup_out_freq == 0) {
      out_years[out_year_idx] = params.spinup_years + params.run_years;
    } else {
      for (int out_year = params.spinup_years + params.rampup_out_freq;
           out_year_idx < out_year_count; out_year += params.rampup_out_freq)
        out_years[out_year_idx++] = out_year;
    }
  }

#elif defined(LPJ_STEP_2)
  if (params.run_out_freq == 0) {    //Doug 07/09: added instances of rampup_years
    if (params.run_years == 0)
      out_years[out_year_idx] = params.spinup_years + params.rampup_years + no_of_months / NMON;
    else                                      
      out_years[out_year_idx] = params.spinup_years + params.rampup_years + params.run_years;
  } else {
    for (int out_year = params.spinup_years + params.rampup_years + params.run_out_freq + params.run_out_years;
         out_year_idx < out_year_count; out_year += params.run_out_freq)
      out_years[out_year_idx++] = out_year;
  }

#else
  if (params.spinup_out_freq != 0) {
    for (int idx = 1;
         idx <= params.spinup_years / params.spinup_out_freq; ++idx)
      out_years[out_year_idx++] = idx * params.spinup_out_freq;
  }

  if (params.rampup_years !=0 ) {
    if (params.rampup_out_freq == 0) {
      out_years[out_year_idx] = params.spinup_years + params.run_years;
    } else {
      for (int out_year = params.spinup_years + params.rampup_out_freq;
           out_year_idx < out_year_count; out_year += params.rampup_out_freq)
        out_years[out_year_idx++] = out_year;
    }
  }

  if (params.run_out_freq == 0) {    //Doug 07/09: added instances of rampup_years
    if (params.run_years == 0)
      out_years[out_year_idx] = params.spinup_years + params.rampup_years + no_of_months / NMON;
    else                                      
      out_years[out_year_idx] = params.spinup_years + params.rampup_years + params.run_years;
  } else {
    for (int out_year = params.spinup_years + params.rampup_years + params.run_out_freq + params.run_out_years;
         out_year_idx < out_year_count; out_year += params.run_out_freq)
      out_years[out_year_idx++] = out_year;
  }
#endif


  // Diagnostic output.

// Doug 07/09: if rampup years havent been set, the stop rpogram if in step_1b mode
#if LPJ_STEP_1B
  if ( params.rampup_years==0) {
    cerr << "ERROR: number of rampup years undefined or set to zero" << endl;
    //stop;
  }
#endif

  cerr << "Total land points:  " << land_points << endl;

#ifdef LPJ_STEP_1A
  cerr << "Spinup written too: " << params.spinup_file << endl;
#endif

#ifdef LPJ_STEP_1B
  cerr << "Spinup read from:   " << params.spinup_file << endl;
#endif

  cerr << "Spinup years:       " << params.spinup_years << endl;
  cerr << "Spinup output freq: ";
  if (params.spinup_out_freq != 0)
    cerr << params.spinup_out_freq << endl;
  else
    cerr << "n/a" << endl;
  if (params.rampup_years != 0) {
#ifdef LPJ_STEP_1B
    cerr << "Rampup written too: " << params.rampup_file << endl;
#endif

#ifdef LPJ_STEP_2
    cerr << "Rampup read from:   " << params.rampup_file << endl;
#endif

    cerr << "Rampup years:       " << params.rampup_years << endl;
    cerr << "Rampup output freq: ";
    if (params.rampup_out_freq != 0)
      cerr << params.rampup_out_freq << endl;
    else
      cerr << "n/a" << endl;
  }
  cerr << "Run years:          ";
  if (params.run_years != 0)
    cerr << params.run_years << endl;
  else
    cerr << "all available" << endl;
  cerr << "Run output freq:    ";
  if (params.run_out_freq != 0) {
    cerr << params.run_out_freq << endl;
    cerr << "  out_years:  " << out_year_idx << endl;
    if (params.run_out_years > 0) {                                // Doug 06/09
      cerr << "  from year:  " << params.run_out_years+1 << endl;  // Doug 06/09
    }
  }
  else
    cerr << "single output file" << endl;
  if( params.output_dir.empty() ) {
      cerr << "output files in working directory" << endl;
  }
  else {
      cerr << "output files in " << params.output_dir << endl;
  }

  // Initialise the grid point counter.

  points_done = 0;


  // Find first valid point.

  latidx = lonidx = 0;
  while (latidx < grid.nlat && !(*land_mask)(lonidx, latidx)) {
    ++lonidx;
    if (lonidx >= grid.nlon) { lonidx = 0; ++latidx; }
  }


  // Open the output file (also initialises variable definitions and
  // buffers).

  open_output_files(params.output_file, params.output_dir);
}

#if defined(LPJ_STEP_1A) || defined(LPJ_STEP_1B) || defined(LPJ_STEP_2)

//--------------------------------------------------------------------
//
//  get_spinup_years
//
//  When doing a 2 step simulation, gets the length of the first step.
//  Doug 07/09: or, now, 1b step as well

extern "C" int get_spinup_years_( int * spinup_years )
{
    *spinup_years = params.spinup_years;
}

//--------------------------------------------------------------------
//
//  get_rampup_years
//
//  Doug 07/09: When doing a 2 step simulation, gets the length of the rampup step.
//

extern "C" int get_rampup_years_( int * rampup_years )
{
    *rampup_years = params.rampup_years;
}

//--------------------------------------------------------------------
//
//  get_run_years
//
//  When doing a 2 step simulation, gets the length of the second step.
//   

extern "C" int get_run_years_( int * run_years )
{
  *run_years = params.run_years;
}

#endif

//----------------------------------------------------------------------
//
//  getgrid
//
//  Find the location of the next land point to process and pass the
//  latitude and soil type to the main LPJ module.  Also deal with
//  end condition if there are no more points to process.

extern "C" int getgrid_(float *latitude, float *longitude, int *soiltype,
                    float *lightn, float *a_nd, int *dogridcell)
{
  // Set up return values.

  if (latidx >= grid.nlat)
    *dogridcell = 0;
  else {
    *dogridcell = 1;
    if (params.soil_type_fixed)
      *soiltype = params.fixed_soil_type;
    else
      *soiltype = (*soil_type)(lonidx, latidx);

    if (params.a_nd_fixed)
      *a_nd = params.fixed_a_nd;
    else
      *a_nd = (*a_nd_vals)(lonidx, latidx);
    *latitude = grid.lats[latidx];
    *longitude = grid.lons[lonidx];

    // Read climate data for this point (note that the windspeed data
    // we're using here is only a seasonal cycle, not a full monthly
    // time series like the other data).

    if (!clim_temp_var->set_cur(0, latidx, lonidx) ||
        !clim_temp_var->get(clim_temp, no_of_months, 1, 1)) {
      cerr << "Failed to read temperature data" << endl;
      exit(1);
    }
    if (!clim_precip_var->set_cur(0, latidx, lonidx) ||
        !clim_precip_var->get(clim_precip, no_of_months, 1, 1)) {
      cerr << "Failed to read precipitation data" << endl;
      exit(1);
    }
    if (!clim_sun_var->set_cur(0, latidx, lonidx) ||
        !clim_sun_var->get(clim_sun, no_of_months, 1, 1)) {
      cerr << "Failed to read cloudiness data" << endl;
      exit(1);
    }
    if (!params.wetdays_fixed) {
      if (!clim_wetdays_var->set_cur(0, latidx, lonidx) ||
          !clim_wetdays_var->get(clim_wetdays, NMON, 1, 1)) {
        cerr << "Failed to read wet days data" << endl;
        exit(1);
      }
    } else {
      for (int idx = 0; idx < NMON; ++idx)
        clim_wetdays[idx] = params.fixed_wetdays;
    }
    if (!clim_tmin_var->set_cur(0, latidx, lonidx) ||
        !clim_tmin_var->get(clim_tmin, no_of_months, 1, 1)) {
      cerr << "Failed to read minimum temperature data" << endl;
      exit(1);
    }
    if (!clim_tmax_var->set_cur(0, latidx, lonidx) ||
        !clim_tmax_var->get(clim_tmax, no_of_months, 1, 1)) {
      cerr << "Failed to read maximum temperature data" << endl;
      exit(1);
    }
    if (!clim_windsp_var->set_cur(0, latidx, lonidx) ||
        !clim_windsp_var->get(clim_windsp, NMON, 1, 1)) {
      cerr << "Failed to read windspeed data" << endl;
      exit(1);
    }
    if (!clim_lightn_var->set_cur(0, latidx, lonidx) ||
        !clim_lightn_var->get(clim_lightn, NMON, 1, 1)) {
      cerr << "Failed to read lightning data" << endl;
      exit(1);
    }
    // DM : lightning is independant of year, no need to read it in getclimate
    for (int i = 0; i < NMON; ++i) {
        lightn[i] = clim_lightn[i];
    }
    if (!clim_ann_popdens_var->set_cur(0, latidx, lonidx) ||
        !clim_ann_popdens_var->get(clim_ann_popdens, no_of_months / NMON, 1, 1)) {
      //check alternatively !clim_ann_popden_var->get(clim_ann_popden, 1, 1, 1))
      cerr << "Failed to read annual population density data" << endl;
      exit(1);
    }
    // Doug 04/09: crops and pasture
    if (!params.crop_fixed) {
      if (!clim_crop_var->set_cur(0, latidx, lonidx) ||
          !clim_crop_var->get(clim_crop, no_of_months / NMON, 1, 1)) {
        //check alternatively !clim_crop_var->get(clim_crop, 1, 1, 1))
        cerr << "Failed to read annual population density data" << endl;
        exit(1);
      }
    }
    if (!params.pas_fixed) {
      if (!clim_pas_var->set_cur(0, latidx, lonidx) ||
          !clim_pas_var->get(clim_pas, no_of_months / NMON, 1, 1)) {
        //check alternatively !clim_pas_var->get(clim_pas, 1, 1, 1))
        cerr << "Failed to read annual population density data" << endl;
        exit(1);
      }
    }

    // Convert units: temperature should be in deg. C (input is K);
    // precipitation should be in mm/month (input is in kg m-2 s-1);
    // sunshine should be in percent (input is percentage cloud).

    for (int mon = 0; mon < no_of_months; ++mon) {
      clim_temp[mon] -= 273.16;
      clim_precip[mon] *= 3600 * 24 * 30;
      
      if (params.sun_is_cloud) clim_sun[mon] = 100.0 - clim_sun[mon];
      //Kirsten: correct cloudiness if >100 contemporary!
      if ((clim_sun[mon]>100.0) && (clim_sun[mon]<1099)) {
	//	cerr << "sunshine > 100 " << clim_sun[mon] << endl;
        clim_sun[mon]=100.0;
        //cerr << "now corrected: " << clim_sun[mon] << endl;
      }
      //Kirsten: correct cloudiness if <0 contemporary!
      if ((clim_sun[mon]<0.0) && (clim_sun[mon]<1099)) {
	//cerr << "sunshine < 0" << clim_sun[mon] << endl;
        clim_sun[mon]=0.0;
        //cerr << "now corrected: " << clim_sun[mon] << endl;
      }

      clim_tmin[mon] -= 273.16;
      clim_tmax[mon] -= 273.16;
    }


    // Reset all output variables.

    for (int idx = 0; idx < out_var_count; ++idx) out_vars[idx]->reset();
    out_idx = 0;
    //printf("out_va_count, %d\n", out_var_count);

  }
  return 0;
}


//----------------------------------------------------------------------
//
//  getclimate
//
//  Sends this year's climatology to the main program.  Also deals
//  with end condition for time.
//  Doug 07/09: add instances of rampup_years

extern "C" int getclimate_(int *year,
                           float *mtemp,
                           float *mtemp_dmin, float *mtemp_dmax,
                           float *mprec, float *mwet, float *msun,
                           float *mwindsp,
                           float *popden,
                           float *crop,	//Doug 04/09
                           float *pas,	//Doug 04/09
                           float *co2, int *doyear)
{
  // The total length of the simulation is params.spinup_years of
  // spin-up, followed by as many years as we have data for.

  int run_years =
    params.run_years == 0 ? no_of_months / NMON : params.run_years;
  if (*year > params.spinup_years + params.rampup_years + run_years) {
    *doyear = 0;
    return 0;
  } else
    *doyear = 1;
  int yearidx;
  if (*year > params.rampup_years + params.spinup_years)
    yearidx = *year - params.spinup_years - params.rampup_years - 1;
  else if (*year > params.spinup_years) 
    yearidx = *year - params.spinup_years - 1;
  else
    yearidx = (*year - 1) % (params.spinup_data_years ?
                             params.spinup_data_years : (no_of_months / NMON));
  if (yearidx >= no_of_months / NMON) yearidx %= no_of_months / NMON; 
  int timeidx = yearidx * NMON;
 

  // Copy climate data for the current month.

  for (int mon = 0; mon < NMON; ++mon) {
    mtemp[mon] = clim_temp[timeidx + mon];
    mprec[mon] = clim_precip[timeidx + mon];
    msun[mon] = clim_sun[timeidx + mon];
    mwet[mon] = clim_wetdays[mon];
    mtemp_dmin[mon] = clim_tmin[timeidx + mon];
    mtemp_dmax[mon] = clim_tmax[timeidx + mon];
    mwindsp[mon] = clim_windsp[mon];
  }

  // Set up the (constant) CO2 concentration.

  co2[0] = co2_vals[yearidx];
  co2[1] = c13_vals[yearidx];
  if( grid.lats[latidx] > 30 ) {
    co2[2] = c14_north_vals[yearidx];
  }
  else if( grid.lats[latidx] > -30 ) {
      co2[2] = c14_equator_vals[yearidx];
  }
  else {
      co2[2] = c14_south_vals[yearidx];
  }
  
  *popden = clim_ann_popdens[yearidx];
  *crop = clim_crop[yearidx];    //Doug 04/09
  *pas = clim_pas[yearidx];      //Doug 04/09

      return 0;
}

#if defined(LPJ_STEP_1A) || defined(LPJ_STEP_1B)

template< typename OUT, typename T >
void output_to_text_file_0( OUT & out, const char * name, T * value )
{
    out << name << " " << *value << "\n";
}

template< typename OUT >
void output_to_text_file_0( OUT & out, const char * name, float * value )
{
    bool ok = isfinite( *value );
#ifdef WARN_ABNORMAL_VALUES
    if( !ok ) cerr << "WARNING : " << name << " has an abnormal value" << endl;
#endif
    float val = ok ? *value : 0.0;
    out << name << " " << scientific << setprecision(12) << val << "\n";
}

template< typename OUT, typename T >
void output_to_text_file_1( OUT & out, const char * name, T * values, int size )
{
    out << name;
    for( int i = 0; i < size; ++i ) {
        out << " " << values[i];
    }
    out << "\n";
}

template< typename OUT >
void output_to_text_file_1( OUT & out, const char * name, float * values, int size )
{
    out << name;
    for( int i = 0; i < size; ++i ) {
        bool ok = isfinite( values[i] );
#ifdef WARN_ABNORMAL_VALUES
        if( !ok ) cerr << "WARNING : " << name << " has an abnormal value at " << i << endl;
#endif
        float val = ok ? values[i] : 0.0;
        out << " " << scientific << setprecision(12) << val;
    }
    out << "\n";
}

template< typename OUT, typename T >
void output_to_text_file_2( OUT & out, const char * name, T * values, int size_1, int size_2 )
{
    out << name;
    for( int i = 0; i < size_1; ++i ) {
        for( int j = 0; j < size_2; ++j ) {
            out << " " << values[j * size_1 + i];
            //out << " " << values[j][i];
        }
    }
    out << "\n";
}

template< typename OUT >
void output_to_text_file_2( OUT & out, const char * name, float * values, int size_1, int size_2 )
{
    out << name;
    for( int i = 0; i < size_1; ++i ) {
        for( int j = 0; j < size_2; ++j ) {
            bool ok = isfinite( values[j * size_1 + i] );
#ifdef WARN_ABNORMAL_VALUES
            if( !ok ) cerr << "WARNING : " << name << " has an abnormal value at (" << i << ", " << j << ")" << endl;
#endif
            float val = ok ? values[j * size_1 + i] : 0;
            out << " " << scientific << setprecision(12) << val;
            //out << " " << values[j][i];
        }
    }
    out << "\n";
}

template< typename OUT, typename T >
void output_to_text_file_3( OUT & out, const char * name, T * values, int size_1, int size_2, int size_3 )
{
    out << name;
    for( int i = 0; i < size_1; ++i ) {
        for( int j = 0; j < size_2; ++j ) {
            for( int k = 0; k < size_3; ++k ) {
                out << " " << values[(( k * size_2 ) + j ) * size_1 + i];
                //out << " " << values[k][j][i];
            }
        }
    }
    out << "\n";
}

template< typename OUT >
void output_to_text_file_3( OUT & out, const char * name, float * values, int size_1, int size_2, int size_3 )
{
    out << name;
    for( int i = 0; i < size_1; ++i ) {
        for( int j = 0; j < size_2; ++j ) {
            for( int k = 0; k < size_3; ++k ) {
                bool ok = isfinite( values[(( k * size_2 ) + j ) * size_1 + i] );
#ifdef WARN_ABNORMAL_VALUES
                if( !ok ) cerr << "WARNING : " << name << " has an abnormal value at (" << i << ", " << j << ", " << k << ")" << endl;
#endif
                float val = ok ? values[(( k * size_2 ) + j ) * size_1 + i] : 0;
                out << " " << scientific << setprecision(12) << val;
                //out << " " << values[k][j][i];
            }
        }
    }
    out << "\n";
}

ofstream year1000out;

//----------------------------------------------------------------------
//
//  init_saved_dataa
//

extern "C" int init_saved_dataa_()
{
    /*
    string out_file_name;
    if( params.output_dir.empty() ) {
        out_file_name  = "Input1/spin_up/europe_data_after_step1_2.txt";
    }
    else if( params.output_dir[params.output_dir.size() - 1] == '/' ) {
        out_file_name  = params.output_dir + "Input1/spin_up/europe_data_after_step1_2.txt";
    }
    else {
        out_file_name  = params.output_dir + "/" + "Input1/spin_up/europe_data_after_step1_2.txt";
    }
    year1000.open( out_file_name.c_str() );
    */

#ifdef LPJ_STEP_1A
    year1000out.open( params.spinup_file.c_str() );
#else
    year1000out.open( params.rampup_file.c_str() );
#endif

    return 0;
}

//----------------------------------------------------------------------
//
//  end_saved_dataa
//

extern "C" int end_saved_dataa_()
{
    year1000out.close();
    return 0;
}
                               

//----------------------------------------------------------------------
//
//  put_saved_data
//

extern "C" int put_saved_data_(
             int *year,
             // initgrid variables
             // tree,  Constant array
             float *k_fast_ave,
             float *k_slow_ave,
             float *litter_decom_ave,
             int   *present,
             float *litter_ag,
             float *fuel_1hr,
             float *fuel_10hr,
             float *fuel_100hr,
             float *fuel_1000hr,
             float *litter_bg,
             float *crownarea,
             // float *mprec,    used by initgrid
             float *w,
             float *w_t,
             float *dwscal365,
             float *lm_ind,
             float *sm_ind,
             float *hm_ind,
             float *rm_ind,
             float *fpc_grid,
             float *mnpp,
             float *anpp,
             int   *leafondays,
             int   *leafoffdays,
             int   *leafon,
             float *snowpack,
             float *mtemp_old,
             // float *mtemp,    used only
             float *maxthaw_old,
             // float *delta_thawday,  unused in lpjmain, removed
             float *mw1,
             float *mw2,
             float *mw1_t,
             float *mw2_t,
             float *uw1,
             float *uw2,
             float *fw,
             // float *soilpar,   used only
             float *mcica,
             float *mgpp,
             float *lresp,
             float *sresp,
             float *rresp,
             float *gresp,
             float *aresp,
             float *dbh,
             float *tau_c,
             float *cl_t,
             float *height_class,
             float *agpp,
             // Others
             float *cpool_fast,
             float *cpool_slow,
             float *bm_inc,
             float *nind,
             float *gdd,
             float *lai_ind,
             float *height,
             float *w_ep,
             float *gdd_buf,
             float *mtemp_min_buf,
             float *mtemp_max_buf,
             float *lm_sapl,
             float *sm_sapl,
             float *hm_sapl,
             float *rm_sapl,
             float *meangc
        )
{
    output_to_text_file_0( year1000out, "latidx", &latidx );
    output_to_text_file_0( year1000out, "lonidx", &lonidx );
    
    output_to_text_file_0( year1000out, "year", year );
    
    // real k_fast_ave   !running average k_fast for subroutine littersom
    output_to_text_file_0( year1000out, "k_fast_ave", k_fast_ave );

    // real k_slow_ave        !running average k_slow for subroutine littersom
    output_to_text_file_0( year1000out, "k_slow_ave", k_slow_ave );

    // real litter_decom_ave(1:nco2)  !running average litter_decom for subroutine littersom
    output_to_text_file_1( year1000out, "litter_decom_ave", litter_decom_ave, NCO2 );

    // logical present(1:npft)         !whether PFT present in gridcell
    output_to_text_file_1( year1000out, "present", present, NPFT );
    
    // real litter_ag(1:npft,1:nco2)          !gridcell above-ground litter (gC/m2)
    output_to_text_file_2( year1000out, "litter_ag", litter_ag, NPFT, NCO2 );
    
    // real fuel_1hr(1:npft,1:nco2)    !1hr dead fuel: dead grass,"cured" grass,shed tree leaves, small twigs
    output_to_text_file_2( year1000out, "fuel_1hr", fuel_1hr, NPFT, NCO2 );
    
    // real fuel_10hr(1:npft,1:nco2)    !1hr dead fuel: dead grass,"cured" grass,shed tree leaves, small twigs
    output_to_text_file_2( year1000out, "fuel_10hr", fuel_10hr, NPFT, NCO2 );
    
    // real fuel_100hr(1:npft,1:nco2)    !1hr dead fuel: dead grass,"cured" grass,shed tree leaves, small twigs
    output_to_text_file_2( year1000out, "fuel_100hr", fuel_100hr, NPFT, NCO2 );
    
    // real fuel_1000hr(1:npft,1:nco2)    !1hr dead fuel: dead grass,"cured" grass,shed tree leaves, small twigs
    output_to_text_file_2( year1000out, "fuel_1000hr", fuel_1000hr, NPFT, NCO2 );
    
    // real litter_bg(1:npft,1:nco2)          !gridcell below-ground litter (gC/m2)
    output_to_text_file_2( year1000out, "litter_bg", litter_bg, NPFT, NCO2 );
    
    // real crownarea(1:npft)          !crown area (m2)
    output_to_text_file_1( year1000out, "crownarea", crownarea, NPFT );

    // real w(1:2)   !soil layer 1 and 2 water content (fraction of available water holding capacity)
    output_to_text_file_1( year1000out, "w", w, 2 );
    
    // real w_t(1:2)
    output_to_text_file_1( year1000out, "w_t", w_t, 2 );
    
    // real dwscal365(1:npft)          !daily water scalar for day 365
    output_to_text_file_1( year1000out, "dwscal365", dwscal365, NPFT );
    
    
    // real lm_ind(1:npft,1:nco2)             !individual leaf mass (gC)
    output_to_text_file_2( year1000out, "lm_ind", lm_ind, NPFT, NCO2 );
    
    // real sm_ind(1:npft,1:nco2)             !individual sapwood mass (gC)
    output_to_text_file_2( year1000out, "sm_ind", sm_ind, NPFT, NCO2 );
    
    // real hm_ind(1:npft,1:nco2)             !individual heartwood mass (gC)
    output_to_text_file_2( year1000out, "hm_ind", hm_ind, NPFT, NCO2 );
    
    // real rm_ind(1:npft,1:nco2)             !individual fine root mass (gC)
    output_to_text_file_2( year1000out, "rm_ind", rm_ind, NPFT, NCO2 );
    
    // real fpc_grid(1:npft)           !gridcell foliar projective cover (FPC)
    output_to_text_file_1( year1000out, "fpc_grid", fpc_grid, NPFT );
    
    // real mnpp(1:12,1:npft,1:nco2)   !monthly gridcell NPP (gC/m2)
    output_to_text_file_3( year1000out, "mnpp", mnpp, 12, NPFT, NCO2 );
    
    // real anpp(1:npft,1:nco2)               !annual gridcell NPP (gC/m2)
    output_to_text_file_2( year1000out, "anpp", anpp, NPFT, NCO2 );
    
    // integer leafondays(1:npft)
    output_to_text_file_1( year1000out, "leafondays", leafondays, NPFT );

    // integer leafoffdays(1:npft)
    output_to_text_file_1( year1000out, "leafoffdays", leafoffdays, NPFT );
    
    // logical leafon(1:npft)
    output_to_text_file_1( year1000out, "leafon", leafon, NPFT );
    
    // real snowpack                   !storage of precipitation as snow (mm)
    output_to_text_file_0( year1000out, "snowpack", snowpack );

    // real mtemp_old(1:12)            !last year's monthly temperatures (deg C)
    output_to_text_file_1( year1000out, "mtemp_old", mtemp_old, 12 );
    
    // real maxthaw_old
    output_to_text_file_0( year1000out, "maxthaw_old", maxthaw_old );
    
    // real delta_thawday,
    //output_to_text_file_0( year1000out, "delta_thawday", delta_thawday );
    
    // real mw1(1:12)       ! monthly values w(1), fraction avail. water
    output_to_text_file_1( year1000out, "mw1", mw1, 12 );
    
    // real mw2(1:12)       ! monthly values w(2), fraction avail. water
    output_to_text_file_1( year1000out, "mw2", mw2, 12 );
    
    // real mw1_t(1:12)     ! monthly values w(1), total water (liquid+ice)
    output_to_text_file_1( year1000out, "mw1_t", mw1_t, 12 );
    
    // real mw2_t(1:12)     ! monthly values w(2)
    output_to_text_file_1( year1000out, "mw2_t", mw2_t, 12 );
    
    // real uw1
    output_to_text_file_0( year1000out, "uw1", uw1 );
    
    // real uw2
    output_to_text_file_0( year1000out, "uw2", uw2 );
    
    // real fw(1:2)
    output_to_text_file_1( year1000out, "fw", fw, 2 );
    
    // real mcica(1:12,1:npft)
    output_to_text_file_2( year1000out, "mcica", mcica, 12, NPFT );
    
    // real mgpp(1:12,1:npft,1:nco2)   !monthly grid cell GPP (gC/m2)
    output_to_text_file_3( year1000out, "mgpp", mgpp, 12, NPFT, NCO2 );
    
    // real lresp(1:12,1:npft) !monthly leaf respiration
    output_to_text_file_2( year1000out, "lresp", lresp, 12, NPFT );
    
    // real sresp(1:12,1:npft)  ! monthly sapwood respiration
    output_to_text_file_2( year1000out, "sresp", sresp, 12, NPFT );
    
    // real rresp(1:12,1:npft)  ! monthly root respiration
    output_to_text_file_2( year1000out, "rresp", rresp, 12, NPFT );
    
    // real gresp(1:12,1:npft)  ! monthly growth respiration
    output_to_text_file_2( year1000out, "gresp", gresp, 12, NPFT );
    
    // real aresp(1:12,1:npft) !monthly autotrophic respiration
    output_to_text_file_2( year1000out, "aresp", aresp, 12, NPFT );
    
    // real dbh(1:npft)                ! stem diameter per PFT, Kirsten
    output_to_text_file_1( year1000out, "dbh",dbh, NPFT );
    
    // real tau_c(0:4,1:npft)    ! critical time to cambial kill, Kirsten
    output_to_text_file_2( year1000out, "tau_c", tau_c, 4, NPFT );
    
    // real cl_t(0:4,1:npft)           !crown length per height class, Kirsten
    output_to_text_file_2( year1000out, "cl_t", cl_t, 4, NPFT );
    
    // real height_class(0:4,1:npft)
    output_to_text_file_2( year1000out, "height_class", height_class, 4, NPFT );
    
    // real agpp(1:npft,1:nco2)       !annual gridcell GPP (gC/m2)
    output_to_text_file_2( year1000out, "agpp", agpp, NPFT, NCO2 );
    
    // Others
    // real cpool_fast(1:nco2)                 !fast-decomposing soil C pool (gC/m2)
    output_to_text_file_1( year1000out, "cpool_fast", cpool_fast, NCO2 );
    
    // real cpool_slow(1:nco2)                !slow-decomposing soil C pool (gC/m2)
    output_to_text_file_1( year1000out, "cpool_slow", cpool_slow, NCO2 );
    
    // real bm_inc(1:npft,1:nco2)      !annual biomass increment (gC/m2)
    output_to_text_file_2( year1000out, "bm_inc", bm_inc, NPFT, NCO2 );

    // real nind(1:npft)               !gridcell individual density (indiv/m2)
    output_to_text_file_1( year1000out, "nind", nind, NPFT );

    // real gdd(1:npft)                !current-year growing degree days
    output_to_text_file_1( year1000out, "gdd", gdd, NPFT );
    
    // real lai_ind(1:npft)            !individual leaf area index
    output_to_text_file_1( year1000out, "lai_ind", lai_ind, NPFT );
    
    // real height(1:npft)             !tree height (m)
    output_to_text_file_1( year1000out, "height", height, NPFT );
    
    // real w_ep                       !fraction of available water holding capiacity evap-layer
    output_to_text_file_0( year1000out, "w_ep", w_ep );
    
    // real gdd_buf(1:npft,1:climbuf)  !buffer to store 'climbuf' years of GDD totals
    output_to_text_file_2( year1000out, "gdd_buf", gdd_buf, NPFT, CLIMBUF );
    
    // real mtemp_min_buf(1:climbuf)   !buffer to store 'climbuf' years of coldest month temperatures
    output_to_text_file_1( year1000out, "mtemp_min_buf", mtemp_min_buf, CLIMBUF );
    
    // real mtemp_max_buf(1:climbuf)   !buffer to store 'climbuf' years of coldest month temperatures
    output_to_text_file_1( year1000out, "mtemp_max_buf", mtemp_max_buf, CLIMBUF );
    
    // real lm_sapl(1:npft,1:nco2)      ! initial (sapling) leaf mass (gC/m2)
    output_to_text_file_2( year1000out, "lm_sapl", lm_sapl, NPFT, NCO2 );
    // real sm_sapl(1:npft,1:nco2)      ! initial (sapling) sapwood mass (gC/m2)
    output_to_text_file_2( year1000out, "sm_sapl", sm_sapl, NPFT, NCO2 );
    // real hm_sapl(1:npft,1:nco2)      ! initial (sapling) heartwood mass (gC/m2)
    output_to_text_file_2( year1000out, "hm_sapl", hm_sapl, NPFT, NCO2 );
    // real rm_sapl(1:npft,1:nco2)      ! initial (sapling) fine root mass (gC/m2)
    output_to_text_file_2( year1000out, "rm_sapl", rm_sapl, NPFT, NCO2 );
    
    // real meangc(1:12,1:npft)
    output_to_text_file_2( year1000out, "meangc", meangc, NMON, NPFT );
    
    return 0;
}

#endif

#if defined(LPJ_STEP_1B) || defined(LPJ_STEP_2)

template< typename IN, typename T >
void input_from_text_file_0( IN & in, const char * name, T value )
{
    string key;
    in >> key >> *value;
    if( key != name ) {
        cerr << name << " expected, got " << key << "\n";
    }
}

template< typename IN, typename T >
void input_from_text_file_1( IN & in, const char * name, T * values, int size )
{
    string key;
    in >> key;
    if( key != name ) {
        cerr << name << " expected, got " << key << "\n";
    }
    for( int i = 0; i < size; ++i ) {
        in >> values[i];
    }
}

template< typename IN, typename T >
void input_from_text_file_2( IN & in, const char * name, T * values, int size_1, int size_2 )
{
    string key;
    in >> key;
    if( key != name ) {
        cerr << name << " expected, got " << key << "\n";
    }
    for( int i = 0; i < size_1; ++i ) {
        for( int j = 0; j < size_2; ++j ) {
            in >> values[j * size_1 + i];
            //in >> values[j][i];
        }
    }
}

template< typename IN, typename T >
void input_from_text_file_3( IN & in, const char * name, T * values, int size_1, int size_2, int size_3 )
{
    string key;
    in >> key;
    if( key != name ) {
        cerr << name << " expected, got " << key << "\n";
    }
    for( int i = 0; i < size_1; ++i ) {
        for( int j = 0; j < size_2; ++j ) {
            for( int k = 0; k < size_3; ++k ) {
                in >> values[(( k * size_2 ) + j ) * size_1 + i];
                //in >> values[k][j][i];
            }
        }
    }
}

ifstream year1000in;

//----------------------------------------------------------------------
//
//  init_saved_datab
//

extern "C" int init_saved_datab_()
{
    /*
    string out_file_name;
    if( params.output_dir.empty() ) {
        out_file_name  = "Input1/spin_up/europe_data_after_step1_2.txt";
    }
    else if( params.output_dir[params.output_dir.size() - 1] == '/' ) {
        out_file_name  = params.output_dir + "Input1/spin_up/europe_data_after_step1_2.txt";
    }
    else {
        out_file_name  = params.output_dir + "/" + "Input1/spin_up/europe_data_after_step1_2.txt";
    }
    year1000.open( out_file_name.c_str() );
    */

#ifdef LPJ_STEP_2
    if (params.rampup_years!=0)
       year1000in.open( params.rampup_file.c_str() );
    else
       year1000in.open( params.spinup_file.c_str() );
    
#else
    year1000in.open( params.spinup_file.c_str() );
#endif
    

    return 0;
}
//----------------------------------------------------------------------
//
//  end_saved_datab
//

extern "C" int end_saved_datab_()
{
    year1000in.close();
    return 0;
}


//----------------------------------------------------------------------
//
//  get_saved_data
//

extern "C" int get_saved_data_(
             int *year,
             // initgrid variables
             // tree,  Constant array
             float *k_fast_ave,
             float *k_slow_ave,
             float *litter_decom_ave,
             int   *present,
             float *litter_ag,
             float *fuel_1hr,
             float *fuel_10hr,
             float *fuel_100hr,
             float *fuel_1000hr,
             float *litter_bg,
             float *crownarea,
             // float *mprec,    used by initgrid
             float *w,
             float *w_t,
             float *dwscal365,
             float *lm_ind,
             float *sm_ind,
             float *hm_ind,
             float *rm_ind,
             float *fpc_grid,
             float *mnpp,
             float *anpp,
             int   *leafondays,
             int   *leafoffdays,
             int   *leafon,
             float *snowpack,
             float *mtemp_old,
             // float *mtemp,    used only
             float *maxthaw_old,
             // float *delta_thawday,  unused in lpjmain, removed
             float *mw1,
             float *mw2,
             float *mw1_t,
             float *mw2_t,
             float *uw1,
             float *uw2,
             float *fw,
             // float *soilpar,   used only
             float *mcica,
             float *mgpp,
             float *lresp,
             float *sresp,
             float *rresp,
             float *gresp,
             float *aresp,
             float *dbh,
             float *tau_c,
             float *cl_t,
             float *height_class,
             float *agpp,
             // Others
             float *cpool_fast,
             float *cpool_slow,
             float *bm_inc,
             float *nind,
             float *gdd,
             float *lai_ind,
             float *height,
             float *w_ep,
             float *gdd_buf,
             float *mtemp_min_buf,
             float *mtemp_max_buf,
             float *lm_sapl,
             float *sm_sapl,
             float *hm_sapl,
             float *rm_sapl,
             float *meangc
        )
{
    int latidx_read, lonidx_read;
    input_from_text_file_0( year1000in, "latidx", &latidx_read );
    input_from_text_file_0( year1000in, "lonidx", &lonidx_read );
    if(( latidx_read != latidx ) || ( lonidx_read != lonidx )) {
        cerr << "get_saved_data_ : the expected point (" << latidx
             << ", " << lonidx << ") is not the one read ("
             << latidx_read << ", " << lonidx_read << ")\n";
        return -1;
    }
    
    int year_read;
    input_from_text_file_0( year1000in, "year", &year_read );
    if( year_read != *year ) {
        cerr << "get_saved_data_ : the expected year ("
             << *year << ") is not the one read ("
             << year_read << ")\n";
        return -1;
    }
    
    // real k_fast_ave   !running average k_fast for subroutine littersom
    input_from_text_file_0( year1000in, "k_fast_ave", k_fast_ave );

    // real k_slow_ave        !running average k_slow for subroutine littersom
    input_from_text_file_0( year1000in, "k_slow_ave", k_slow_ave );

    // real litter_decom_ave(1:nco2)  !running average litter_decom for subroutine littersom
    input_from_text_file_1( year1000in, "litter_decom_ave", litter_decom_ave, NCO2 );

    // logical present(1:npft)         !whether PFT present in gridcell
    input_from_text_file_1( year1000in, "present", present, NPFT );
    
    // real litter_ag(1:npft,1:nco2)          !gridcell above-ground litter (gC/m2)
    input_from_text_file_2( year1000in, "litter_ag", litter_ag, NPFT, NCO2 );
    
    // real fuel_1hr(1:npft,1:nco2)    !1hr dead fuel: dead grass,"cured" grass,shed tree leaves, small twigs
    input_from_text_file_2( year1000in, "fuel_1hr", fuel_1hr, NPFT, NCO2 );
    
    // real fuel_10hr(1:npft,1:nco2)    !1hr dead fuel: dead grass,"cured" grass,shed tree leaves, small twigs
    input_from_text_file_2( year1000in, "fuel_10hr", fuel_10hr, NPFT, NCO2 );
    
    // real fuel_100hr(1:npft,1:nco2)    !1hr dead fuel: dead grass,"cured" grass,shed tree leaves, small twigs
    input_from_text_file_2( year1000in, "fuel_100hr", fuel_100hr, NPFT, NCO2 );
    
    // real fuel_1000hr(1:npft,1:nco2)    !1hr dead fuel: dead grass,"cured" grass,shed tree leaves, small twigs
    input_from_text_file_2( year1000in, "fuel_1000hr", fuel_1000hr, NPFT, NCO2 );
    
    // real litter_bg(1:npft,1:nco2)          !gridcell below-ground litter (gC/m2)
    input_from_text_file_2( year1000in, "litter_bg", litter_bg, NPFT, NCO2 );
    
    // real crownarea(1:npft)          !crown area (m2)
    input_from_text_file_1( year1000in, "crownarea", crownarea, NPFT );

    // real w(1:2)   !soil layer 1 and 2 water content (fraction of available water holding capacity)
    input_from_text_file_1( year1000in, "w", w, 2 );
    
    // real w_t(1:2)
    input_from_text_file_1( year1000in, "w_t", w_t, 2 );
    
    // real dwscal365(1:npft)          !daily water scalar for day 365
    input_from_text_file_1( year1000in, "dwscal365", dwscal365, NPFT );
    
    // real lm_ind(1:npft,1:nco2)             !individual leaf mass (gC)
    input_from_text_file_2( year1000in, "lm_ind", lm_ind, NPFT, NCO2 );
    
    // real sm_ind(1:npft,1:nco2)             !individual sapwood mass (gC)
    input_from_text_file_2( year1000in, "sm_ind", sm_ind, NPFT, NCO2 );
    
    // real hm_ind(1:npft,1:nco2)             !individual heartwood mass (gC)
    input_from_text_file_2( year1000in, "hm_ind", hm_ind, NPFT, NCO2 );
    
    // real rm_ind(1:npft,1:nco2)             !individual fine root mass (gC)
    input_from_text_file_2( year1000in, "rm_ind", rm_ind, NPFT, NCO2 );
    
    // real fpc_grid(1:npft)           !gridcell foliar projective cover (FPC)
    input_from_text_file_1( year1000in, "fpc_grid", fpc_grid, NPFT );
    
    // real mnpp(1:12,1:npft,1:nco2)   !monthly gridcell NPP (gC/m2)
    input_from_text_file_3( year1000in, "mnpp", mnpp, 12, NPFT, NCO2 );
    
    // real anpp(1:npft,1:nco2)               !annual gridcell NPP (gC/m2)
    input_from_text_file_2( year1000in, "anpp", anpp, NPFT, NCO2 );
    
    // integer leafondays(1:npft)
    input_from_text_file_1( year1000in, "leafondays", leafondays, NPFT );

    // integer leafoffdays(1:npft)
    input_from_text_file_1( year1000in, "leafoffdays", leafoffdays, NPFT );
    
    // logical leafon(1:npft)
    input_from_text_file_1( year1000in, "leafon", leafon, NPFT );
    
    // real snowpack                   !storage of precipitation as snow (mm)
    input_from_text_file_0( year1000in, "snowpack", snowpack );

    // real mtemp_old(1:12)            !last year's monthly temperatures (deg C)
    input_from_text_file_1( year1000in, "mtemp_old", mtemp_old, 12 );
    
    // real maxthaw_old
    input_from_text_file_0( year1000in, "maxthaw_old", maxthaw_old );
    
    // real delta_thawday,
    // input_from_text_file_0( year1000in, "delta_thawday", delta_thawday );
    
    // real mw1(1:12)       ! monthly values w(1), fraction avail. water
    input_from_text_file_1( year1000in, "mw1", mw1, 12 );
    
    // real mw2(1:12)       ! monthly values w(2), fraction avail. water
    input_from_text_file_1( year1000in, "mw2", mw2, 12 );
    
    // real mw1_t(1:12)     ! monthly values w(1), total water (liquid+ice)
    input_from_text_file_1( year1000in, "mw1_t", mw1_t, 12 );
    
    // real mw2_t(1:12)     ! monthly values w(2)
    input_from_text_file_1( year1000in, "mw2_t", mw2_t, 12 );
    
    // real uw1
    input_from_text_file_0( year1000in, "uw1", uw1 );
    
    // real uw2
    input_from_text_file_0( year1000in, "uw2", uw2 );
    
    // real fw(1:2)
    input_from_text_file_1( year1000in, "fw", fw, 2 );
    
    // real mcica(1:12,1:npft)
    input_from_text_file_2( year1000in, "mcica", mcica, 12, NPFT );
    
    // real mgpp(1:12,1:npft,1:nco2)   !monthly grid cell GPP (gC/m2)
    input_from_text_file_3( year1000in, "mgpp", mgpp, 12, NPFT, NCO2 );
    
    // real lresp(1:12,1:npft) !monthly leaf respiration
    input_from_text_file_2( year1000in, "lresp", lresp, 12, NPFT );
    
    // real sresp(1:12,1:npft)  ! monthly sapwood respiration
    input_from_text_file_2( year1000in, "sresp", sresp, 12, NPFT );
    
    // real rresp(1:12,1:npft)  ! monthly root respiration
    input_from_text_file_2( year1000in, "rresp", rresp, 12, NPFT );
    
    // real gresp(1:12,1:npft)  ! monthly growth respiration
    input_from_text_file_2( year1000in, "gresp", gresp, 12, NPFT );
    
    // real aresp(1:12,1:npft) !monthly autotrophic respiration
    input_from_text_file_2( year1000in, "aresp", aresp, 12, NPFT );
    
    // real dbh(1:npft)                ! stem diameter per PFT, Kirsten
    input_from_text_file_1( year1000in, "dbh", dbh, NPFT );
    
    // real tau_c(0:4,1:npft)    ! critical time to cambial kill, Kirsten
    input_from_text_file_2( year1000in, "tau_c", tau_c, 4, NPFT );
    
    // real cl_t(0:4,1:npft)           !crown length per height class, Kirsten
    input_from_text_file_2( year1000in, "cl_t", cl_t, 4, NPFT );
    
    // real height_class(0:4,1:npft)
    input_from_text_file_2( year1000in, "height_class", height_class, 4, NPFT );
    
    // real agpp(1:npft,1:nco2)       !annual gridcell GPP (gC/m2)
    input_from_text_file_2( year1000in, "agpp", agpp, NPFT, NCO2 );
    
    // Others
    // real cpool_fast(1:nco2)                 !fast-decomposing soil C pool (gC/m2)
    input_from_text_file_1( year1000in, "cpool_fast", cpool_fast, NCO2 );
    
    // real cpool_slow(1:nco2)                !slow-decomposing soil C pool (gC/m2)
    input_from_text_file_1( year1000in, "cpool_slow", cpool_slow, NCO2 );
    
    
    // real bm_inc(1:npft,1:nco2)      !annual biomass increment (gC/m2)
    input_from_text_file_2( year1000in, "bm_inc", bm_inc, NPFT, NCO2 );
    
    // real nind(1:npft)               !gridcell individual density (indiv/m2)
    input_from_text_file_1( year1000in, "nind", nind, NPFT );
    
    // real gdd(1:npft)                !current-year growing degree days
    input_from_text_file_1( year1000in, "gdd", gdd, NPFT );
    
    // real lai_ind(1:npft)            !individual leaf area index
    input_from_text_file_1( year1000in, "lai_ind", lai_ind, NPFT );
    
    // real height(1:npft)             !tree height (m)
    input_from_text_file_1( year1000in, "height", height, NPFT );
    
    // real w_ep                       !fraction of available water holding capiacity evap-layer
    input_from_text_file_0( year1000in, "w_ep", w_ep );
    
    // real gdd_buf(1:npft,1:climbuf)  !buffer to store 'climbuf' years of GDD totals
    input_from_text_file_2( year1000in, "gdd_buf", gdd_buf, NPFT, CLIMBUF );
    
    // real mtemp_min_buf(1:climbuf)   !buffer to store 'climbuf' years of coldest month temperatures
    input_from_text_file_1( year1000in, "mtemp_min_buf", mtemp_min_buf, CLIMBUF );
    
    // real mtemp_max_buf(1:climbuf)   !buffer to store 'climbuf' years of coldest month temperatures
    input_from_text_file_1( year1000in, "mtemp_max_buf", mtemp_max_buf, CLIMBUF );
    
    // real lm_sapl(1:npft,1:nco2)      ! initial (sapling) leaf mass (gC/m2)
    input_from_text_file_2( year1000in, "lm_sapl", lm_sapl, NPFT, NCO2 );
    // real sm_sapl(1:npft,1:nco2)      ! initial (sapling) sapwood mass (gC/m2)
    input_from_text_file_2( year1000in, "sm_sapl", sm_sapl, NPFT, NCO2 );
    // real hm_sapl(1:npft,1:nco2)      ! initial (sapling) heartwood mass (gC/m2)
    input_from_text_file_2( year1000in, "hm_sapl", hm_sapl, NPFT, NCO2 );
    // real rm_sapl(1:npft,1:nco2)      ! initial (sapling) fine root mass (gC/m2)
    input_from_text_file_2( year1000in, "rm_sapl", rm_sapl, NPFT, NCO2 );
    
    // real meangc(1:12,1:npft)
    input_from_text_file_2( year1000in, "meangc", meangc, NMON, NPFT );
    
    return 0;
}

#endif

//----------------------------------------------------------------------
//
//  outannual
//
//  Collect a year's worth of output for one grid cell.

extern "C" int outannual_(int *year, int *present,
                          float *nind,
                          float *lm_ind, float *rm_ind,
                          float *sm_ind, float *hm_ind,
                          float *fpc_grid,
                          float *anpp, float *acflux_estab,
                          float *litter_ag, float *litter_bg,
                          float *cpool_fast, float *cpool_slow,
                          float *arh, float *afire_frac, float *acflux_fire,
                          float *mcflux_fire,
                          float *arunoff, float *sla, float *mpar,
                          float *mapar, float *mphen, float *anpp_add,
                          float *mnpp,
                          float *mnpp_grid,    //Doug 06/09
                          float *mrunoff, float *aaet,
                          float *mrh, float *mnpp_add, float *maet,
                          float *mtemp_soil, float *mcica,
                          float *lresp, float *sresp, float *rresp,
                          float *gresp, float *aresp,
                          float *mauw1, float *mauw2,
                          float *maep, float *aaep,
                          float *mpet_grid, float *apet_grid,
                          float *mintc, float *aintc,
                          float *meangc, float *mgp,
                          float *num_fire, float *annum_fire,
                          float *area_burnt, float *an_areafires, 
                          float *mfdi, float *an_fdi,
                          float *an_fseason, 
                          float *acflux_trace, float *mcflux_trace, 
                          float *m_fc_crown, float *an_fc_crown, 
                          float *m_i_surface, float *an_i_surface, 
                          float *gdd, float *height, float *mgpp,
                          float *wu, float *wl, float *dphen,
                          float *dphen_change,                    // Doug 05/09
                          float *dhuman_ign, float *dlightn,      // Yan: 24/10/07
                          float *mw1, float *mw2,      // Yan: 22/11/07
                          float *num_fire_human, float *num_fire_lightn,             // Yan: 22/11/07
                          float *annum_fire_human, float *annum_fire_lightn,         // Yan: 22/11/07
                          float *area_burnt_human, float *area_burnt_lightn,         // Yan: 22/11/07
                          float *an_areafires_human, float *an_areafires_lightn,     // Yan: 22/11/07
                          float *afire_frac_human, float *afire_frac_lightn,         // Yan: 22/11/07
			  float *char_net_fuel_0,		
			  float *livegrass_0, float *dead_fuel_0,             
                          float *dead_fuel_all_0, float *fuel_all_0,	
                          float *fuel_1hr_total_0, float *fuel_10hr_total_0,		
                          float *fuel_100hr_total_0, float *fuel_1000hr_total_0,
                          float *mfuel_1hr_total, float *mfuel_10hr_total,           //Doug 03/09: monthly fuel
                          float *mfuel_100hr_total, float *mfuel_1000hr_total,       //Doug 03/09: monthly fuel
                          float *mlivegrass,                                         //Doug 05/09: monthly fuel
                          float *deltaa, float *deltaa_fpc, float *delt_c13_fpc, 	
                          float *mfire_frac, 
                          float *fbdep, float * litter_decom_ave, float * turnover_ind,
                          float *crop, float *pas,                                    //Doug 05/09: crop and pasture proportions
                          float *anpp_grid, float *arh_grid, float *acflux_fire_grid,//Doug 07/09: the cheats
                          float *gdd_grid, float *alpha_ws, float *cgf, float *fdry,
                          float *lt_days)                          //Doug 07/09: vioclimatic varibles         
{
  // Skip output during spin-up if required.
  // Doug 06/09: Skip output if before fist year of ouput
  // Doug 07/09: Add params.rampup_years and move all into one logic test to save runtime
  if (*year <= params.spinup_years) {
    if (params.spinup_out_freq == 0)
      return 0;
  } else if (*year <= params.rampup_years+params.spinup_years) {
    if (params.rampup_out_freq == 0)
      return 0;
  } else if (*year < params.spinup_years + params.rampup_years + params.run_out_years + 1)
    return 0;



  // Calculate bare ground fraction.
  float fbare = 0.0;
  for (int idx = 0; idx < NPFT; ++idx) fbare += fpc_grid[idx];
  fbare = 1.0 - fbare;
  if (fbare > 1.0) fbare = 1.0;
  if (fbare < 0.0) fbare = 0.0;

  //Sarah soilc=litter_bg(1:9) + cpool_slow + cpool_fast
  float litter_bg_sum[NCO2] = { 0.0, 0.0, 0.0 };
  float soilc[NCO2] = { 0.0, 0.0, 0.0 };
  float anpp_pft_sum[NCO2] = { 0.0, 0.0, 0.0 };
  float anpp_sum[NCO2] = { 0.0, 0.0, 0.0 };

  
  for (int idx = 0; idx < NPFT; ++idx) {
    litter_bg_sum[0] += litter_bg[0 * NPFT + idx];
  }
    
  for (int jdx = 1; jdx < NCO2; ++jdx) {
    for (int idx = 0; idx < NPFT; ++idx) {
      litter_bg_sum[jdx] += litter_bg[jdx * NPFT + idx] * litter_bg[0 * NPFT + idx];
      anpp_pft_sum[jdx] += anpp[jdx * NPFT + idx] * anpp[0 * NPFT + idx];
    }
    litter_bg_sum[jdx] /= litter_bg_sum[0];
    anpp_pft_sum[jdx] /= anpp_pft_sum[0];
  }
   
  //for (int jdx = 0; jdx < NCO2; ++jdx) {
  //  soilc[jdx] = litter_bg_sum[jdx] + cpool_slow[jdx] + cpool_fast[jdx]; 
  //  anpp_sum[jdx] = anpp_pft_sum[jdx] + anpp_add[jdx] + acflux_estab[jdx];
  //}
  // Marko's mail
  soilc[0]    = litter_bg_sum[0] + cpool_slow[0] + cpool_fast[0]; 
  anpp_sum[0] = anpp_pft_sum[0] + anpp_add[0] + acflux_estab[0];
  for (int jdx = 1; jdx < NCO2; ++jdx) {
    soilc[jdx] = litter_bg_sum[jdx] * litter_bg_sum[0]
               + cpool_slow[jdx]    * cpool_slow[0]
               + cpool_fast[jdx]    * cpool_fast[0]; 
    anpp_sum[jdx] = anpp_pft_sum[jdx] * anpp_pft_sum[0]
                  + anpp_add[jdx]     * anpp_add[0]
                  + acflux_estab[jdx] * acflux_estab[0];
  }

  for (int jdx = 1; jdx < NCO2; ++jdx) {
    soilc[jdx]    /= soilc[0];
    anpp_sum[jdx] /= anpp_sum[0];
  }


  // Accumulate values for averaging.
  handle_output_record("dhuman_ign", dhuman_ign);   // Yan: 24/10/07
  handle_output_record("dlightn", dlightn);         // Yan: 24/10/07
  handle_output_record("nind", nind);
  handle_output_record("lm_ind", lm_ind);
  handle_output_record("sm_ind", sm_ind);
  handle_output_record("hm_ind", hm_ind);
  handle_output_record("rm_ind", rm_ind);
  handle_output_record("fpc_grid", fpc_grid);
  handle_output_record("anpp_sum", anpp_sum);
  handle_output_record("acflux_estab", acflux_estab);
  handle_output_record("litter_ag", litter_ag);
  handle_output_record("litter_bg", litter_bg);
  handle_output_record("cpool_fast", cpool_fast);
  handle_output_record("cpool_slow", cpool_slow);
  handle_output_record("arh", arh);
  handle_output_record("afire_frac", afire_frac);
  handle_output_record("acflux_fire", acflux_fire);
  handle_output_record("mcflux_fire", mcflux_fire);
  handle_output_record("arunoff", arunoff);
  handle_output_record("sla", sla);
  handle_output_record("mpar", mpar);
  handle_output_record("mapar", mapar);
  handle_output_record("mphen", mphen);
  handle_output_record("anpp_add", anpp_add);
  handle_output_record("mnpp", mnpp);
  handle_output_record("mnpp_grid", mnpp_grid);  //Doug 06/09
  handle_output_record("anpp", anpp);        // Yan 20/11/07
  handle_output_record("mrunoff", mrunoff);
  handle_output_record("aaet", aaet);
  handle_output_record("mrh", mrh);
  handle_output_record("mnpp_add", mnpp_add);
  handle_output_record("maet", maet);
  handle_output_record("mgpp", mgpp);
  handle_output_record("mtemp_soil", mtemp_soil);
  handle_output_record("mcica", mcica);
  handle_output_record("lresp", lresp);
  handle_output_record("sresp", sresp);
  handle_output_record("rresp", rresp);
  handle_output_record("gresp", gresp);
  handle_output_record("aresp", aresp);
  handle_output_record("mauw1", mauw1);
  handle_output_record("mauw2", mauw2);
  handle_output_record("maep", maep);
  handle_output_record("aaep", aaep);
  handle_output_record("mpet_grid", mpet_grid);
  handle_output_record("apet_grid", apet_grid);
  handle_output_record("mintc", mintc);
  handle_output_record("aintc", aintc);
  handle_output_record("meangc", meangc);
  handle_output_record("mgp", mgp);
  handle_output_record("num_fire", num_fire);
  handle_output_record("annum_fire", annum_fire);
  handle_output_record("area_burnt", area_burnt);
  handle_output_record("an_areafires", an_areafires);
  handle_output_record("mfdi", mfdi);
  handle_output_record("an_fdi", an_fdi);
  handle_output_record("an_fseason", an_fseason);
  handle_output_record("acflux_trace", acflux_trace);
  handle_output_record("mcflux_trace", mcflux_trace);
  handle_output_record("m_fc_crown", m_fc_crown);
  handle_output_record("an_fc_crown", an_fc_crown);
  handle_output_record("m_i_surface", m_i_surface);
  handle_output_record("an_i_surface", an_i_surface);
  handle_output_record("gdd", gdd);
  handle_output_record("height", height);
//  handle_output_record("mgpp", mgpp);
  handle_output_record("wu", wu);
  handle_output_record("wl", wl);
  handle_output_record("dphen", dphen);
  handle_output_record("dphen_change", dphen_change);
  handle_output_record("fbare", &fbare);
  handle_output_record("mw1", mw1);    // Yan 22/11/07
  handle_output_record("mw2", mw2);    // Yan 22/11/07
  handle_output_record("num_fire_human", num_fire_human);            // Yan 22/11/07
  handle_output_record("num_fire_lightn", num_fire_lightn);          // Yan 22/11/07
  handle_output_record("annum_fire_human", annum_fire_human);        // Yan 22/11/07
  handle_output_record("annum_fire_lightn", annum_fire_lightn);      // Yan 22/11/07
  handle_output_record("area_burnt_human", area_burnt_human);        // Yan 22/11/07
  handle_output_record("area_burnt_lightn", area_burnt_lightn);      // Yan 22/11/07
  handle_output_record("an_areafires_human", an_areafires_human);    // Yan 22/11/07
  handle_output_record("an_areafires_lightn", an_areafires_lightn);  // Yan 22/11/07
  handle_output_record("afire_frac_human", afire_frac_human);        // Yan 22/11/07
  handle_output_record("afire_frac_lightn", afire_frac_lightn);      // Yan 22/11/07
  handle_output_record("char_net_fuel_0",char_net_fuel_0);       
  handle_output_record("dead_fuel_0",dead_fuel_0);       
  handle_output_record("dead_fuel_all_0",dead_fuel_all_0);       
  handle_output_record("livegrass_0",livegrass_0);       
  handle_output_record("fuel_all_0",fuel_all_0);       
  handle_output_record("fuel_1hr_total_0",fuel_1hr_total_0);       
  handle_output_record("fuel_10hr_total_0",fuel_10hr_total_0);      
  handle_output_record("fuel_100hr_total_0",fuel_100hr_total_0);     
  handle_output_record("fuel_1000hr_total_0",fuel_1000hr_total_0);
  handle_output_record("mfuel_1hr_total",mfuel_1hr_total);           // Doug 03/09
  handle_output_record("mfuel_10hr_total",mfuel_10hr_total);         // Doug 03/09
  handle_output_record("mfuel_100hr_total",mfuel_100hr_total);       // Doug 03/09
  handle_output_record("mfuel_1000hr_total",mfuel_1000hr_total);     // Doug 03/09
  handle_output_record("mlivegrass",mlivegrass);                     // Doug 05/09
  handle_output_record("deltaa",deltaa);       
  handle_output_record("deltaa_fpc",deltaa_fpc);       
  handle_output_record("delt_c13_fpc",delt_c13_fpc);  
  handle_output_record("mfire_frac",mfire_frac);     
  handle_output_record("fbdep",fbdep);       
  handle_output_record("litter_decom_ave",litter_decom_ave);    
  handle_output_record("turnover_ind",turnover_ind);       
  // Sarah
  //handle_output_record("soilc", soilc);
  handle_output_record("crop", crop);                                //Doug 05/09: proportion of cropland
  handle_output_record("pas", pas);                                  //Doug 05/09: proportion of pasture
  handle_output_record("anpp_grid",anpp_grid);                       //Doug 07/09: a cheat
  handle_output_record("arh_grid",arh_grid);                         //Doug 07/09: another cheat
  handle_output_record("acflux_fire_grid",acflux_fire_grid);         //Doug 07/09: another cheat
  handle_output_record("gdd_grid",gdd_grid);                         //Doug 07/09: growing degree days base 5
  handle_output_record("alpha_ws",  alpha_ws);                       //Doug 07/09: bioclimatic alpha
  handle_output_record("cgf",  cgf); 
  handle_output_record("fdry",  fdry); 
  handle_output_record("lt_days",  lt_days);   //Doug 07/09: bioclimatic alpha

  // Output is written every params.spinup_out_freq years during
  // spin-up (averaged over that period).  During the transient part
  // of the run, data is either collected to write an average value
  // for each relevant variable at the end of the run, or average
  // values are written out every params.run_out_freq years.  A zero
  // value for params.spinup_out_freq means to write no output during
  // spin-up.


  if (*year <= params.spinup_years) {
    if (*year % params.spinup_out_freq != 0)
      return 0;
  } else if (*year <= params.rampup_years + params.spinup_years) {           //Doug 07/09
    if ((*year-params.spinup_years) % params.rampup_out_freq != 0)    //Doug 07/09
      return 0;                                                       //Doug 07/09
  } else {
    int run_years =
      params.run_years == 0 ? no_of_months / NMON : params.run_years;
    if (params.run_out_freq == 0) {
      if (*year < params.spinup_years + params.rampup_years + run_years)
        return 0;
    } else {
  //    if (*year < params.spinup_years + params.run_out_years+1)  //Doug 06/09
  //      return 0;                        //Doug 06/09
  //    else {
        if((*year-params.spinup_years-params.rampup_years)%params.run_out_freq!=0)
          return 0;
      
    }
  }


  // Average and output variables and reinitialise averaging buffers.
  for (int idx = 0; idx < out_var_count; ++idx) {
    if (!out_vars[idx]->average_and_output()) {
      cerr << "ERROR: failed writing variable: "
           << out_vars[idx]->def->name << endl;
      exit(1);
    }
  }

  out_nc[out_idx++]->sync();


  return 0;

}


//----------------------------------------------------------------------
//
//  outgrid
//
//  Write progress message and find next valid point.
//

extern "C" int outgrid_(void)
{
  // Progress message.

  if (++points_done % 10 == 0)
    cerr << points_done << " / " << land_points << " points processed / total land points ("
         << points_done * 100 / land_points << "%)" << endl;
  
  // Find next valid point.

  ++lonidx;
  if (lonidx >= grid.nlon) { lonidx = 0; ++latidx; }
  while (latidx < grid.nlat && !(*land_mask)(lonidx, latidx)) {
    ++lonidx;
    if (lonidx >= grid.nlon) { lonidx = 0; ++latidx; }
  }

  return 0;
}


//----------------------------------------------------------------------
//
//  termio
//
//  Close output file.
//

extern "C" int termio_(void)
{
  for (int idx = 0; idx < out_year_count; ++idx) delete out_nc[idx];
  delete [] out_nc;
  return 0;
}


//======================================================================
//
//  LOCAL FUNCTION DEFINITIONS
//
//======================================================================

//----------------------------------------------------------------------
//
//  LPJOutputVariable::LPJOutputVariable
//
//  Constructor to create a single variable entry - contains
//  information about the variable definition, the output netCDF
//  variable, and an averaging buffer.
//

LPJOutputVariable::LPJOutputVariable(LPJVariable *in_def,
                                     RunParameters::OutVariable &collect)
{
  def = in_def;
  switch (def->type) {
    case SIMPLE:          buff_size = 1;             break;
    case SIMPLE_MONTHLY:  buff_size = NMON;          break;
    case SIMPLE_DAILY:    buff_size = NDAYS;         break;
    case PER_PFT:         buff_size = NPFT;          break;
    case PER_PFT_MONTHLY: buff_size = NMON * NPFT;   break;
    case PER_PFT_DAILY:   buff_size = NDAYS * NPFT;  break;
    case TRACE:           buff_size = NTRACE;        break;
    case TRACE_MONTHLY:   buff_size = NTRACE * NMON; break;
    case PER_CO2:         buff_size = NCO2;          break;
    case PER_CO2_MONTHLY: buff_size = NCO2 * NMON;   break;
    case PER_CO2_PFT:     buff_size = NPFT * NCO2;   break;
  }
  if (def->type == SIMPLE_MONTHLY  ||
      def->type == PER_PFT_MONTHLY ||
      def->type == TRACE_MONTHLY   
      || def->type == PER_CO2_MONTHLY
      )
    month_dim_needed = true;
  if (def->type == SIMPLE_DAILY || def->type == PER_PFT_DAILY)
    day_dim_needed = true;
  if (def->type == TRACE || def->type == TRACE_MONTHLY)
    trace_dim_needed = true;
  if (def->type == PER_CO2  ||
      def->type == PER_CO2_MONTHLY ||
      def->type == PER_CO2_PFT)
    co2_dim_needed = true;
  mean = collect.mean;
  sdev = collect.sdev;
  min = collect.min;
  max = collect.max;
  avg_buff = sdev_buff = min_buff = max_buff = 0;
  if (mean || sdev) avg_buff = new float[buff_size];
  if (sdev) sdev_buff = new float[buff_size];
  if (min) min_buff = new float[buff_size];
  if (max) max_buff = new float[buff_size];
  if (mean) mean_nc = new NcVar*[out_year_count];
  if (sdev) sdev_nc = new NcVar*[out_year_count];
  if (min) min_nc = new NcVar*[out_year_count];
  if (max) max_nc = new NcVar*[out_year_count];
  clear_buffs();
  out_index = 0;
}



LPJOutputVariable::~LPJOutputVariable()
{
  delete [] avg_buff;
  delete [] sdev_buff;
  delete [] min_buff;
  delete [] max_buff;
}


//----------------------------------------------------------------------
//
//  LPJOutputVariable::set_nc
//
//  Set up a netCDF representation of a variable in a single output
//  file.
//

void LPJOutputVariable::set_nc(int idx, NcFile &n)
{
  string sdev_nm = string(def->name) + "_SDEV";
  string min_nm = string(def->name) + "_MIN";
  string max_nm = string(def->name) + "_MAX";
  switch (def->type) {
  case SIMPLE:
    if (mean)
      mean_nc[idx] = n.add_var(def->name, ncFloat, lat_dim, lon_dim);
    if (sdev)
      sdev_nc[idx] = n.add_var(sdev_nm.c_str(), ncFloat,
                               lat_dim, lon_dim);
    if (min)
      min_nc[idx] = n.add_var(min_nm.c_str(), ncFloat,
                              lat_dim, lon_dim);
    if (max)
      max_nc[idx] = n.add_var(max_nm.c_str(), ncFloat,
                              lat_dim, lon_dim);
    break;
  case SIMPLE_MONTHLY:
    if (mean)
      mean_nc[idx] = n.add_var(def->name, ncFloat,
                               month_dim, lat_dim, lon_dim);
    if (sdev)
      sdev_nc[idx] = n.add_var(sdev_nm.c_str(), ncFloat,
                               month_dim, lat_dim, lon_dim);
    if (min)
      min_nc[idx] = n.add_var(min_nm.c_str(), ncFloat,
                              month_dim, lat_dim, lon_dim);
    if (max)
      max_nc[idx] = n.add_var(max_nm.c_str(), ncFloat,
                              month_dim, lat_dim, lon_dim);
    break;
  case SIMPLE_DAILY:
    if (mean)
      mean_nc[idx] = n.add_var(def->name, ncFloat,
                               day_dim, lat_dim, lon_dim);
    if (sdev)
      sdev_nc[idx] = n.add_var(sdev_nm.c_str(), ncFloat,
                               day_dim, lat_dim, lon_dim);
    if (min)
      min_nc[idx] = n.add_var(min_nm.c_str(), ncFloat,
                              day_dim, lat_dim, lon_dim);
    if (max)
      max_nc[idx] = n.add_var(max_nm.c_str(), ncFloat,
                              day_dim, lat_dim, lon_dim);
    break;
  case PER_PFT:
    if (mean)
      mean_nc[idx] = n.add_var(def->name, ncFloat,
                               pft_dim, lat_dim, lon_dim);
    if (sdev)
      sdev_nc[idx] = n.add_var(sdev_nm.c_str(), ncFloat,
                               pft_dim, lat_dim, lon_dim);
    if (min)
      min_nc[idx] = n.add_var(min_nm.c_str(), ncFloat,
                              pft_dim, lat_dim, lon_dim);
    if (max)
      max_nc[idx] = n.add_var(max_nm.c_str(), ncFloat,
                              pft_dim, lat_dim, lon_dim);
    break;
  case PER_PFT_MONTHLY:
    if (mean)
      mean_nc[idx] = n.add_var(def->name, ncFloat,
                               month_dim, pft_dim, lat_dim, lon_dim);
    if (sdev)
      sdev_nc[idx] = n.add_var(sdev_nm.c_str(), ncFloat,
                               month_dim, pft_dim, lat_dim, lon_dim);
    if (min)
      min_nc[idx] = n.add_var(min_nm.c_str(), ncFloat,
                              month_dim, pft_dim, lat_dim, lon_dim);
    if (max)
      max_nc[idx] = n.add_var(max_nm.c_str(), ncFloat,
                              month_dim, pft_dim, lat_dim, lon_dim);
    break;
  
  case PER_PFT_DAILY:
    if (mean)
      mean_nc[idx] = n.add_var(def->name, ncFloat,
                               day_dim, pft_dim, lat_dim, lon_dim);
    if (sdev)
      sdev_nc[idx] = n.add_var(sdev_nm.c_str(), ncFloat,
                               day_dim, pft_dim, lat_dim, lon_dim);
    if (min)
      min_nc[idx] = n.add_var(min_nm.c_str(), ncFloat,
                              day_dim, pft_dim, lat_dim, lon_dim);
    if (max)
      max_nc[idx] = n.add_var(max_nm.c_str(), ncFloat,
                              day_dim, pft_dim, lat_dim, lon_dim);
    break;
  case TRACE:
    if (mean)
      mean_nc[idx] = n.add_var(def->name, ncFloat,
                               trace_dim, lat_dim, lon_dim);
    if (sdev)
      sdev_nc[idx] = n.add_var(sdev_nm.c_str(), ncFloat,
                               trace_dim, lat_dim, lon_dim);
    if (min)
      min_nc[idx] = n.add_var(min_nm.c_str(), ncFloat,
                              trace_dim, lat_dim, lon_dim);
    if (max)
      max_nc[idx] = n.add_var(max_nm.c_str(), ncFloat,
                              trace_dim, lat_dim, lon_dim);
    break;
  case TRACE_MONTHLY:
    if (mean)
      mean_nc[idx] = n.add_var(def->name, ncFloat,
                               month_dim, trace_dim, lat_dim, lon_dim);
    if (sdev)
      sdev_nc[idx] = n.add_var(sdev_nm.c_str(), ncFloat,
                               month_dim, trace_dim, lat_dim, lon_dim);
    if (min)
      min_nc[idx] = n.add_var(min_nm.c_str(), ncFloat,
                              month_dim, trace_dim, lat_dim, lon_dim);
    if (max)
      max_nc[idx] = n.add_var(max_nm.c_str(), ncFloat,
                              month_dim, trace_dim, lat_dim, lon_dim);
    break;
  case PER_CO2:
    if (mean)
      mean_nc[idx] = n.add_var(def->name, ncFloat,
                               co2_dim, lat_dim, lon_dim);
    if (sdev)
      sdev_nc[idx] = n.add_var(sdev_nm.c_str(), ncFloat,
                               co2_dim, lat_dim, lon_dim);
    if (min)
      min_nc[idx] = n.add_var(min_nm.c_str(), ncFloat,
                              co2_dim, lat_dim, lon_dim);
    if (max)
      max_nc[idx] = n.add_var(max_nm.c_str(), ncFloat,
                              co2_dim, lat_dim, lon_dim);
    break;
  case PER_CO2_MONTHLY:
    if (mean)
      mean_nc[idx] = n.add_var(def->name, ncFloat,
                               month_dim, co2_dim, lat_dim, lon_dim);
    if (sdev)
      sdev_nc[idx] = n.add_var(sdev_nm.c_str(), ncFloat,
                               month_dim, co2_dim, lat_dim, lon_dim);
    if (min)
      min_nc[idx] = n.add_var(min_nm.c_str(), ncFloat,
                              month_dim, co2_dim, lat_dim, lon_dim);
    if (max)
      max_nc[idx] = n.add_var(max_nm.c_str(), ncFloat,
                              month_dim, co2_dim, lat_dim, lon_dim);
    break;
  case PER_CO2_PFT:
    if (mean)
      mean_nc[idx] = n.add_var(def->name, ncFloat,
                               pft_dim, co2_dim, lat_dim, lon_dim);
    if (sdev)
      sdev_nc[idx] = n.add_var(sdev_nm.c_str(), ncFloat,
                               pft_dim, co2_dim, lat_dim, lon_dim);
    if (min)
      min_nc[idx] = n.add_var(min_nm.c_str(), ncFloat,
                              pft_dim, co2_dim, lat_dim, lon_dim);
    if (max)
      max_nc[idx] = n.add_var(max_nm.c_str(), ncFloat,
                              pft_dim, co2_dim, lat_dim, lon_dim);
    break;
  }
  if (mean && !mean_nc[idx]) {
    cerr << "ERROR: Failed to add variable " << def->name << endl;
    exit(1);
  }
  if (sdev && !sdev_nc[idx]) {
    cerr << "ERROR: Failed to add variable " << sdev_nm.c_str() << endl;
    exit(1);
  }
  if (min && !min_nc[idx]) {
    cerr << "ERROR: Failed to add variable " << min_nm.c_str() << endl;
    exit(1);
  }
  if (max && !max_nc[idx]) {
    cerr << "ERROR: Failed to add variable " << max_nm.c_str() << endl;
    exit(1);
  }

  if (mean) {
    mean_nc[idx]->add_att("missing_value", 1.0E20f);
    mean_nc[idx]->add_att("_FillValue", 1.0E20f);
  }
  if (sdev) {
    sdev_nc[idx]->add_att("missing_value", 1.0E20f);
    sdev_nc[idx]->add_att("_FillValue", 1.0E20f);
  }
  if (min) {
    min_nc[idx]->add_att("missing_value", 1.0E20f);
    min_nc[idx]->add_att("_FillValue", 1.0E20f);
  }
  if (max) {
    max_nc[idx]->add_att("missing_value", 1.0E20f);
    max_nc[idx]->add_att("_FillValue", 1.0E20f);
  }
}


//----------------------------------------------------------------------
//
//  Clear the averaging buffer for a variable.
//

void LPJOutputVariable::clear_buffs(void)
{
  for (int idx = 0; idx < buff_size; ++idx) {
    if (avg_buff) avg_buff[idx] = 0.0f;
    if (sdev_buff) sdev_buff[idx] = 0.0f;
    if (min_buff) min_buff[idx] = 1.0E20f;
    if (max_buff) max_buff[idx] = -1.0E20f;
  }

  accum_years = 0;
}


//----------------------------------------------------------------------
//
//  Average and output values for a variable.
//

bool LPJOutputVariable::average_and_output(void)
{
  bool retval = true;
  // Calculate mean and standard deviations as required.
  if (accum_years > 0) {
    for (int idx = 0; idx < buff_size; ++idx) {
      // At this point, sdev_buff[idx] contains the sum of squares of
      // the data and avg_buff[idx] the sum of the data.
      if (sdev) {
        // rounding may lead to negative value
        float v = sdev_buff[idx] * accum_years - avg_buff[idx] * avg_buff[idx];
          if( v < 0 ) {
            //if( v < -0.001 ) cerr << "Unexpected negative value (" << v << ") in average_and_output for " << def->name << "\n";
            v = 0;
          }
        sdev_buff[idx] = sqrt(v / (accum_years * accum_years));
      }
      if (mean) avg_buff[idx] /= accum_years;
    }
  }
  // Output the required data.


  switch (def->type) {
  case SIMPLE:
    if (mean)
      if (!mean_nc[out_index]->set_cur(latidx, lonidx) ||
          !mean_nc[out_index]->put(avg_buff, 1, 1))
        retval = false;
    if (sdev)
      if (!sdev_nc[out_index]->set_cur(latidx, lonidx) ||
          !sdev_nc[out_index]->put(sdev_buff, 1, 1))
        retval = false;
    if (min)
      if (!min_nc[out_index]->set_cur(latidx, lonidx) ||
          !min_nc[out_index]->put(min_buff, 1, 1))
        retval = false;
    if (max)
      if (!max_nc[out_index]->set_cur(latidx, lonidx) ||
          !max_nc[out_index]->put(max_buff, 1, 1))
        retval = false;
    break;
  case SIMPLE_MONTHLY:
    if (mean)
      if (!mean_nc[out_index]->set_cur(0, latidx, lonidx) ||
          !mean_nc[out_index]->put(avg_buff, NMON, 1, 1))
        retval = false;
    if (sdev)
      if (!sdev_nc[out_index]->set_cur(0, latidx, lonidx) ||
          !sdev_nc[out_index]->put(sdev_buff, NMON, 1, 1))
        retval = false;
    if (min)
      if (!min_nc[out_index]->set_cur(0, latidx, lonidx) ||
          !min_nc[out_index]->put(min_buff, NMON, 1, 1))
        retval = false;
    if (max)
      if (!max_nc[out_index]->set_cur(0, latidx, lonidx) ||
          !max_nc[out_index]->put(max_buff, NMON, 1, 1))
        retval = false;
    break;
  case SIMPLE_DAILY:
    if (mean)
      if (!mean_nc[out_index]->set_cur(0, latidx, lonidx) ||
          !mean_nc[out_index]->put(avg_buff, NDAYS, 1, 1))
        retval = false;
    if (sdev)
      if (!sdev_nc[out_index]->set_cur(0, latidx, lonidx) ||
          !sdev_nc[out_index]->put(sdev_buff, NDAYS, 1, 1))
        retval = false;
    if (min)
      if (!min_nc[out_index]->set_cur(0, latidx, lonidx) ||
          !min_nc[out_index]->put(min_buff, NDAYS, 1, 1))
        retval = false;
    if (max)
      if (!max_nc[out_index]->set_cur(0, latidx, lonidx) ||
          !max_nc[out_index]->put(max_buff, NDAYS, 1, 1))
        retval = false;
    break;
  case PER_PFT:
    if (mean)
      if (!mean_nc[out_index]->set_cur(0, latidx, lonidx) ||
          !mean_nc[out_index]->put(avg_buff, NPFT, 1, 1))
        retval = false;
    if (sdev)
      if (!sdev_nc[out_index]->set_cur(0, latidx, lonidx) ||
          !sdev_nc[out_index]->put(sdev_buff, NPFT, 1, 1))
        retval = false;
    if (min)
      if (!min_nc[out_index]->set_cur(0, latidx, lonidx) ||
          !min_nc[out_index]->put(min_buff, NPFT, 1, 1))
        retval = false;
    if (max)
      if (!max_nc[out_index]->set_cur(0, latidx, lonidx) ||
          !max_nc[out_index]->put(max_buff, NPFT, 1, 1))
        retval = false;
    break;
  case PER_PFT_MONTHLY:
    transpose_buffers(NMON, NPFT);
    if (mean)
      if (!mean_nc[out_index]->set_cur(0, 0, latidx, lonidx) ||
          !mean_nc[out_index]->put(avg_buff, NMON, NPFT, 1, 1))
        retval = false;
    if (sdev)
      if (!sdev_nc[out_index]->set_cur(0, 0, latidx, lonidx) ||
          !sdev_nc[out_index]->put(sdev_buff, NMON, NPFT, 1, 1))
        retval = false;
    if (min)
      if (!min_nc[out_index]->set_cur(0, 0, latidx, lonidx) ||
          !min_nc[out_index]->put(min_buff, NMON, NPFT, 1, 1))
        retval = false;
    if (max)
      if (!max_nc[out_index]->set_cur(0, 0, latidx, lonidx) ||
          !max_nc[out_index]->put(max_buff, NMON, NPFT, 1, 1))
        retval = false;
    break;
  case PER_PFT_DAILY:
    transpose_buffers(NDAYS, NPFT);
    if (mean)
      if (!mean_nc[out_index]->set_cur(0, 0, latidx, lonidx) ||
          !mean_nc[out_index]->put(avg_buff, NDAYS, NPFT, 1, 1))
        retval = false;
    if (sdev)
      if (!sdev_nc[out_index]->set_cur(0, 0, latidx, lonidx) ||
          !sdev_nc[out_index]->put(sdev_buff, NDAYS, NPFT, 1, 1))
        retval = false;
    if (min)
      if (!min_nc[out_index]->set_cur(0, 0, latidx, lonidx) ||
          !min_nc[out_index]->put(min_buff, NDAYS, NPFT, 1, 1))
        retval = false;
    if (max)
      if (!max_nc[out_index]->set_cur(0, 0, latidx, lonidx) ||
          !max_nc[out_index]->put(max_buff, NDAYS, NPFT, 1, 1))
        retval = false;
    break;
  case TRACE:
    if (mean)
      if (!mean_nc[out_index]->set_cur(0, latidx, lonidx) ||
          !mean_nc[out_index]->put(avg_buff, NTRACE, 1, 1))
        retval = false;
    if (sdev)
      if (!sdev_nc[out_index]->set_cur(0, latidx, lonidx) ||
          !sdev_nc[out_index]->put(sdev_buff, NTRACE, 1, 1))
        retval = false;
    if (min)
      if (!min_nc[out_index]->set_cur(0, latidx, lonidx) ||
          !min_nc[out_index]->put(min_buff, NTRACE, 1, 1))
        retval = false;
    if (max)
      if (!max_nc[out_index]->set_cur(0, latidx, lonidx) ||
          !max_nc[out_index]->put(avg_buff, NTRACE, 1, 1))
        retval = false;
    break;
  case TRACE_MONTHLY:
    transpose_buffers(NMON, NTRACE);
    if (mean)
      if (!mean_nc[out_index]->set_cur(0, 0, latidx, lonidx) ||
          !mean_nc[out_index]->put(avg_buff, NMON, NTRACE, 1, 1))
        retval = false;
    if (sdev)
      if (!sdev_nc[out_index]->set_cur(0, 0, latidx, lonidx) ||
          !sdev_nc[out_index]->put(sdev_buff, NMON, NTRACE, 1, 1))
        retval = false;
    if (min)
      if (!min_nc[out_index]->set_cur(0, 0, latidx, lonidx) ||
          !min_nc[out_index]->put(min_buff, NMON, NTRACE, 1, 1))
        retval = false;
    if (max)
      if (!max_nc[out_index]->set_cur(0, 0, latidx, lonidx) ||
          !max_nc[out_index]->put(max_buff, NMON, NTRACE, 1, 1))
        retval = false;
    break;
  case PER_CO2:
    if (mean)
      if (!mean_nc[out_index]->set_cur(0, latidx, lonidx) ||
          !mean_nc[out_index]->put(avg_buff, NCO2, 1, 1))
        retval = false;
    if (sdev)
      if (!sdev_nc[out_index]->set_cur(0, latidx, lonidx) ||
          !sdev_nc[out_index]->put(sdev_buff, NCO2, 1, 1))
        retval = false;
    if (min)
      if (!min_nc[out_index]->set_cur(0, latidx, lonidx) ||
          !min_nc[out_index]->put(min_buff, NCO2, 1, 1))
        retval = false;
    if (max)
      if (!max_nc[out_index]->set_cur(0, latidx, lonidx) ||
          !max_nc[out_index]->put(max_buff, NCO2, 1, 1))
        retval = false;
    break;
  case PER_CO2_MONTHLY:
    transpose_buffers(NMON, NCO2);
    if (mean)
      if (!mean_nc[out_index]->set_cur(0, 0, latidx, lonidx) ||
          !mean_nc[out_index]->put(avg_buff, NMON, NCO2, 1, 1))
        retval = false;
    if (sdev)
      if (!sdev_nc[out_index]->set_cur(0, 0, latidx, lonidx) ||
          !sdev_nc[out_index]->put(sdev_buff, NMON, NCO2, 1, 1))
        retval = false;
    if (min)
      if (!min_nc[out_index]->set_cur(0, 0, latidx, lonidx) ||
          !min_nc[out_index]->put(min_buff, NMON, NCO2, 1, 1))
        retval = false;
    if (max)
      if (!max_nc[out_index]->set_cur(0, 0, latidx, lonidx) ||
          !max_nc[out_index]->put(max_buff, NMON, NCO2, 1, 1))
        retval = false;
    break;
  case PER_CO2_PFT:
    transpose_buffers(NPFT, NCO2);
    if (mean)
      if (!mean_nc[out_index]->set_cur(0, 0, latidx, lonidx) ||
          !mean_nc[out_index]->put(avg_buff, NPFT, NCO2, 1, 1))
        retval = false;
    if (sdev)
      if (!sdev_nc[out_index]->set_cur(0, 0, latidx, lonidx) ||
          !sdev_nc[out_index]->put(sdev_buff, NPFT, NCO2, 1, 1))
        retval = false;
    if (min)
      if (!min_nc[out_index]->set_cur(0, 0, latidx, lonidx) ||
          !min_nc[out_index]->put(min_buff, NPFT, NCO2, 1, 1))
        retval = false;
    if (max)
      if (!max_nc[out_index]->set_cur(0, 0, latidx, lonidx) ||
          !max_nc[out_index]->put(max_buff, NPFT, NCO2, 1, 1))
        retval = false;
    break;
  }

  ++out_index;
  clear_buffs();
  return retval;
}

void LPJOutputVariable::transpose_buffers(int time_axis_len,
                                          int other_axis_len)
{
  if (mean) {
    float *tmp = new float[buff_size];
    for (int pft = 0; pft < other_axis_len; ++pft)
      for (int tim = 0; tim < time_axis_len; ++tim)
        tmp[tim * other_axis_len + pft] = avg_buff[pft * time_axis_len + tim];
    delete [] avg_buff;
    avg_buff = tmp;
  }
  if (sdev) {
    float *tmp = new float[buff_size];
    for (int pft = 0; pft < other_axis_len; ++pft)
      for (int tim = 0; tim < time_axis_len; ++tim)
        tmp[tim * other_axis_len + pft] =
          sdev_buff[pft * time_axis_len + tim];
    delete [] sdev_buff;
    sdev_buff = tmp;
  }
  if (min) {
    float *tmp = new float[buff_size];
    for (int pft = 0; pft < other_axis_len; ++pft)
      for (int tim = 0; tim < time_axis_len; ++tim)
        tmp[tim * other_axis_len + pft] = min_buff[pft * time_axis_len + tim];
    delete [] min_buff;
    min_buff = tmp;
  }
  if (max) {
    float *tmp = new float[buff_size];
    for (int pft = 0; pft < other_axis_len; ++pft)
      for (int tim = 0; tim < time_axis_len; ++tim)
        tmp[tim * other_axis_len + pft] = max_buff[pft * time_axis_len + tim];
    delete [] max_buff;
    max_buff = tmp;
  }
}

//----------------------------------------------------------------------
//
//  read_config
//
//  Read lpj.cfg configuration file and set up internal parameter
//  structure.
//

static void read_config(RunParameters &params)
{
  // Keep track of the keywords used in the file.

  set<string> keywords;

  // Open the configuration file and process a line at a time.

  FILE *fp = fopen("lpj.cfg", "r");
  if (!fp) {
    cerr << "Configuration file lpj.cfg not found!" << endl;
    exit(1);
  }
  int line = 0;
  bool error = false;
  for (;;) {
    // Read lines, removing leading and trailing blanks and skipping
    // blank lines and comment lines (introduced by a # character).

    char buff[132], *pt = buff;
    fgets(buff, 132, fp);
    if (feof(fp)) break;
    ++line;
    if (strchr(buff, '\n')) *strchr(buff, '\n') = '\0';
    string in = pt;
    while (isspace(in[0])) in.erase(0, 1);
    while (isspace(in[in.size() - 1])) in.erase(in.size() - 1);
    if (in.size() == 0 || in[0] == '#') continue;


    // Split the line into the keyword and value (separated by a colon).

    string::size_type colon = in.find(':');
    if (colon == string::npos) {
      cerr << "lpj.cfg:" << line
           << ": no colon separator between keyword and value" << endl;
      error = true;
      continue;
    }
    string kw = in.substr(0, colon);
    string value = in.substr(colon + 1);
    while (isspace(value[0])) value.erase(0, 1);
    while (isspace(kw[kw.size() - 1])) kw.erase(kw.size() - 1);
    for (string::size_type idx = 0; idx < kw.size(); ++idx)
      kw[idx] = toupper(kw[idx]);


    // Now split on the keyword.

    if (kw != "OUTPUT_VAR" && keywords.find(kw) != keywords.end()) {
      cerr << "lpj.cfg:" << line
           << ": duplicate keyword '" << kw << "'" << endl;
      error = true;
      continue;
    }
    keywords.insert(kw);
    if (kw == "TEMP")             params.temp_file = value;
    else if (kw == "PREC")        params.prec_file = value;
    else if (kw == "SUN")         params.sun_file = value;
    else if (kw == "WETDAYS") {
      bool valid_number = true;
      for (int idx = 0; idx < value.size(); ++idx)
        if (!isspace(value[idx]) && !isdigit(value[idx]) &&
            value[idx] != '.') {
          valid_number = false;
          break;
        }
      if (!valid_number)
        params.wetdays_file = value;
      else {
        params.wetdays_fixed = true;
        params.fixed_wetdays = atof(value.c_str());
      }
    }
    else if (kw == "TMIN")        params.tmin_file = value;
    else if (kw == "TMAX")        params.tmax_file = value;
    else if (kw == "WINDSPEED")   params.windsp_file = value;
    else if (kw == "MASK")        params.mask_file = value;
    else if (kw == "TEMP_VAR")    params.temp_var = value;
    else if (kw == "PREC_VAR")    params.prec_var = value;
    else if (kw == "SUN_VAR")     params.sun_var = value;
    else if (kw == "WETDAYS_VAR") params.wetdays_var = value;
    else if (kw == "TMIN_VAR")      params.tmin_var = value;
    else if (kw == "TMAX_VAR")      params.tmax_var = value;
    else if (kw == "WINDSPEED_VAR") params.windsp_var = value;
    else if (kw == "MASK_VAR")    params.mask_var = value;
    else if (kw == "SUN_IS_CLOUD")
      params.sun_is_cloud = boolean_string(value, line, error);
    else if (kw == "MASK_TYPE") {
      if (value == "boolean" || value == "BOOLEAN")
        params.mask_type = RunParameters::BOOLEAN;
      else if (value == "percent_land" || value == "PERCENT_LAND")
        params.mask_type = RunParameters::PERCENT_LAND;
      else {
        cerr << "lpj.cfg:" << line << ": invalid mask type" << endl;
        error = true;
      }
    } else if (kw == "SOIL_TYPE") {
      params.soil_type_fixed = value.size() == 1 && isdigit(value[0]);
      if (params.soil_type_fixed)
        params.fixed_soil_type = value[0] - '0';
      else
        params.soil_type_file = value;
    } else if (kw == "SOIL_TYPE_VAR")  params.soil_type_var = value;
    else if (kw == "POPDENS")          params.ann_popdens_file = value;
    else if (kw == "POPDENS_VAR")      params.ann_popdens_var = value;
    else if (kw == "CROP") {    // Doug 04/09
      bool valid_number = true;
      for (int idx = 0; idx < value.size(); ++idx)
        if (!isspace(value[idx]) && !isdigit(value[idx]) &&
            value[idx] != '.') {
          valid_number = false;
          break;
        }
      if (!valid_number)
        params.crop_file = value;
      else {
        params.crop_fixed = true;
        params.fixed_crop = atof(value.c_str());
      }
    }
    else if (kw == "CROP_VAR")          params.crop_var = value;     // Doug 04/09
    else if (kw == "PAS") {    // Doug 04/09
      bool valid_number = true;
      for (int idx = 0; idx < value.size(); ++idx)
        if (!isspace(value[idx]) && !isdigit(value[idx]) &&
            value[idx] != '.') {
          valid_number = false;
          break;
        }
      if (!valid_number)
        params.pas_file = value;
      else {
        params.pas_fixed = true;
        params.fixed_pas = atof(value.c_str());
      }
    }
    else if (kw == "PAS_VAR")          params.pas_var = value;     // Doug 04/09
    else if (kw == "LIGHTN")           params.lightn_file = value;
    else if (kw == "LIGHTN_VAR")       params.lightn_var = value;
    else if (kw == "A_ND") {
      params.a_nd_fixed = true;
      for (int idx = 0; idx < value.size(); ++idx)
        if (!isdigit(value[idx]) && value[idx] != '.') {
          params.a_nd_fixed = false;
          break;
        }
      if (params.a_nd_fixed)
        params.fixed_a_nd = atof(value.c_str());
      else
        params.a_nd_file = value;
    }
    else if (kw == "A_ND_VAR")  params.a_nd_var = value;
    else if (kw == "SPINUP_YEARS")
      params.spinup_years = atoi(value.c_str());
    else if (kw == "RAMPUP_YEARS")                         //Doug 07/09
      params.rampup_years = atoi(value.c_str());           //Doug 07/09
    else if (kw == "RUN_YEARS")
      params.run_years = atoi(value.c_str());
    else if (kw == "SPINUP_DATA_YEARS")
      params.spinup_data_years = atoi(value.c_str());
    else if (kw == "RAMPUP_DATA_YEARS")                    //Doug 07/09
      params.rampup_years = atoi(value.c_str());           //Doug 07/09
    else if (kw == "SPINUP_OUT_FREQ")
      params.spinup_out_freq = atoi(value.c_str());
    else if (kw == "RAMPUP_OUT_FREQ")                      //Doug 07/09
      params.rampup_out_freq = atoi(value.c_str());        //Doug 07/09
    else if (kw == "SPINUP_FILE")                          //Doug 07/09
      params.spinup_file = value;            //Doug 07/09
    else if (kw == "RAMPUP_FILE")                          //Doug 07/09
      params.rampup_file = value;           //Doug 07/09
    else if (kw == "RUN_OUT_FREQ")
      params.run_out_freq = atoi(value.c_str());
    else if (kw == "RUN_OUT_YEARS")                //Doug 06/09
      params.run_out_years = atoi(value.c_str())-1;  //Doug 06/09
    else if (kw == "CO2_FILE")           params.co2_file         = value;
    else if (kw == "C13_FILE")           params.c13_file         = value;
    else if (kw == "C14_NORTH_FILE")     params.c14_north_file   = value;
    else if (kw == "C14_EQUATOR_FILE")   params.c14_equator_file = value;
    else if (kw == "C14_SOUTH_FILE")     params.c14_south_file   = value;
    else if (kw == "CO2_VAR")            params.co2_var          = value;
    else if (kw == "C13_VAR")            params.c13_var          = value;
    else if (kw == "C14_VAR")            params.c14_var          = value;
    else if (kw == "OUTPUT_FILE")  params.output_file = value;
    else if (kw == "OUTPUT_DIR")  params.output_dir = value;
    else if (kw == "OUTPUT_VAR") {
      RunParameters::OutVariable out_var;
      string::size_type space = value.find(' ');
      if (space == string::npos) {
        // No specifiers - just collect the mean value.
        out_var.name = value;
        out_var.mean = true;
      } else {
        out_var.name = value.substr(0, space);
        value.erase(0, space);
        while (value[0] == ' ') value.erase(0, 1);
        if (value.size() == 0 ||
            value[0] != '[' || value[value.size() - 1] != ']') {
          cerr << "lpj.cfg:" << line << ": invalid variable specifier" << endl;
          error = true;
        } else {
          value.erase(0, 1);
          value.erase(value.size() - 1, 1);
          vector<string> specifiers;
          while (value.find(',') != string::npos) {
            specifiers.push_back(value.substr(0, value.find(',')));
            value.erase(0, value.find(',') + 1);
          }
          specifiers.push_back(value);
          for (vector<string>::iterator it = specifiers.begin();
               it != specifiers.end(); ++it)
            if (*it == "MEAN") out_var.mean = true;
            else if (*it == "SDEV") out_var.sdev = true;
            else if (*it == "MAX") out_var.max = true;
            else if (*it == "MIN") out_var.min = true;
            else {
              cerr << "lpj.cfg:" << line
                   << ": invalid variable specifier" << endl;
              error = true;
            }
        }
      }
      if (!error) params.output_vars.push_back(out_var);
    }
  }
  // Doug 08/09: If ouput years in cfg file are set funny, then just ouput last year
  if(params.run_out_freq>params.run_years) {
    cerr << "ERROR: RUN_OUT_FREQ in cfg file greater then RUN_YEARS \n";
    cerr << "       RUN_OUT_FREQ set to average over entire run. \n";
    cerr << "       Press enter to continue \n";
    cin.get ();
    params.run_out_freq=params.run_years;
    params.run_out_years=0;
  }
  if(params.run_out_years>params.run_years) {
    cerr << "ERROR: RUN_OUT_YEARS in cfg file greater then RUN_YEARS \n";
    cerr << "       RUN_OUT_YEARS set to last year of run: \n";
    cerr << "                                             ";
    cerr << params.run_years << "\n";
    cerr << "       Press enter to continue \n";
    cin.get ();
    params.run_out_years=params.run_years-1;
  }

  // These entries are obligatory in the configuration file.

  char *mandatory[] = { "TEMP", "TMIN", "TMAX", "PREC", "WETDAYS", "SUN",
                        "WINDSPEED", "MASK",
                        "POPDENS",
                        "CROP",                                  // Doug 04/09
                        "PAS",                                   // Doug 04/09
                        "LIGHTN", "A_ND",                        // delect CO2_CONC
                        "CO2_FILE", "C13_FILE",
                        "C14_NORTH_FILE", "C14_EQUATOR_FILE", "C14_SOUTH_FILE",
                        "OUTPUT_FILE" };
  for (int check = 0; check < sizeof(mandatory) / sizeof(char *); ++check) {
    char *check_word = mandatory[check];
    if (keywords.find(check_word) == keywords.end()) {
      cerr << "lpj.cfg: mandatory keywords '"<< check_word
           << "' missing" << endl;
      error = true;
    }
  }


  if (params.output_vars.size() != 0) {
    // Check that the specified output variables exist.

    for (vector<RunParameters::OutVariable>::iterator it =
           params.output_vars.begin();
         it != params.output_vars.end(); ++it) {
      string check_str = it->name;
      bool found = false;
      for (int chk = 0;
           chk < sizeof(LPJ_VARIABLES) / sizeof(LPJVariable); ++chk)
        if (check_str == LPJ_VARIABLES[chk].name) {
          found = true;
          break;
        }
      if (!found) {
        cerr << "lpj.cfg: unrecognised variable '"<< check_str << "'" << endl;
        error = true;
      }
    }
  } else {
    // Add a default minimum set of output variables if none have been
    // specified.

    params.output_vars.push_back(RunParameters::OutVariable("fpc_grid"));
    params.output_vars.push_back(RunParameters::OutVariable("anpp"));
    params.output_vars.push_back(RunParameters::OutVariable("nind"));
    params.output_vars.push_back(RunParameters::OutVariable("gdd"));
    params.output_vars.push_back(RunParameters::OutVariable("height"));
    params.output_vars.push_back(RunParameters::OutVariable("fbare"));
  }

  // Quit on errors.

  if (error) {
    cerr << "Error in configuration file: exiting" << endl;
    exit(1);
  }
}


//----------------------------------------------------------------------
//
//  boolean_string
//
//  Convert a configuration file string value to a boolean, detecting
//  invalid values.
//

static bool boolean_string(string value, int line, bool &error_flag)
{
  if (value == "true" || value == "TRUE" ||
      value == "yes" || value == "YES" ||
      value == "t" || value == "T" ||
      value == "y" || value == "Y" ||
      value == "1")
    return true;
  else if (value == "false" || value == "FALSE" ||
      value == "no" || value == "NO" ||
      value == "f" || value == "F" ||
      value == "n" || value == "N" ||
      value == "0")
    return false;
  else {
    cerr << "lpj.cfg:" << line << ": invalid boolean value" << endl;
    error_flag = true;
    return false;
  }
}


//----------------------------------------------------------------------
//
//  Grid member functions
//

Grid::Grid(string file)
{
  // Open the netCDF file.

  NcFile nc(file.c_str());
  if (!nc.is_valid()) {
    cerr << "ERROR: could not open netCDF file " << file << endl;
    exit(1);
  }


  // Access latitude and longitude dimensions and variables.

  NcDim *lat_dim = nc.get_dim("lat");
  if (!lat_dim) {
    cerr << "ERROR: no 'lat' dimension in netCDF file " << file << endl;
    exit(1);
  }
  NcDim *lon_dim = nc.get_dim("lon");
  if (!lon_dim) {
    cerr << "ERROR: no 'lon' dimension in netCDF file " << file << endl;
    exit(1);
  }
  NcVar *lat_var = nc.get_var("lat");
  if (!lat_var) {
    cerr << "ERROR: no 'lat' variable in netCDF file " << file << endl;
    exit(1);
  }
  NcVar *lon_var = nc.get_var("lon");
  if (!lon_var) {
    cerr << "ERROR: no 'lon' variable in netCDF file " << file << endl;
    exit(1);
  }


  // Check that the dimensions of the coordinate variables are
  // correct.

  if (lat_var->num_dims() != 1 || lat_var->get_dim(0) != lat_dim) {
    cerr << "ERROR: bad dimensions for 'lat' variable in netCDF file "
         << file << endl;
    exit(1);
  }
  if (lon_var->num_dims() != 1 || lon_var->get_dim(0) != lon_dim) {
    cerr << "ERROR: bad dimensions for 'lon' variable in netCDF file "
         << file << endl;
    exit(1);
  }


  // Set up the latitude and longitude arrays.

  lats.resize(nlat = lat_dim->size());
  lons.resize(nlon = lon_dim->size());
  for (int idx = 0; idx < nlat; ++idx) lats[idx] = lat_var->as_float(idx);
  for (int idx = 0; idx < nlon; ++idx) lons[idx] = lon_var->as_float(idx);
}


Grid::Grid(const Grid &other) :
  nlat(other.nlat), nlon(other.nlon), lats(other.lats), lons(other.lons)
{}


Grid &Grid::operator=(const Grid &other)
{
  if (this != &other) {
    nlat = other.nlat;  nlon = other.nlon;
    lats = other.lats;  lons = other.lons;
  }
  return *this;
}


bool Grid::operator==(const Grid &other)
{
  return (nlat == other.nlat && nlon == other.nlon &&
          lats == other.lats && lons == other.lons);
}


//----------------------------------------------------------------------
//
//  read_mask
//
//  Read land/sea mask from netCDF file.
//

static void read_mask(BoolArray *mask, string mask_file, Grid &grid,
                      string mask_var, RunParameters::MASK_TYPE mask_type)
{
  // Open the mask file.

  NcFile nc(mask_file.c_str());
  if (!nc.is_valid()) {
    cerr << "ERROR: could not open mask file " << mask_file << endl;
    exit(1);
  }


  // Access the mask variable and check the dimensions are correct.

  NcVar *var = nc.get_var(mask_var.c_str());
  NcDim *lat_dim = nc.get_dim("lat");
  NcDim *lon_dim = nc.get_dim("lon");
  if (!var) {
    cerr << "ERROR: mask variable '" << mask_var
         << "' not found in mask file " << mask_file << endl;
    exit(1);
  }
  if (var->num_dims() != 2 ||
      var->get_dim(0) != lat_dim || var->get_dim(1) != lon_dim) {
    cerr << "ERROR: mask variable '" << mask_var
         << "' has invalid dimensions" << endl;
    exit(1);
  }


  // Get the mask data and convert to a boolean mask.

  if (mask_type == RunParameters::BOOLEAN) {
    int tmp[grid.nlon];
    for (int latidx = 0; latidx < grid.nlat; ++latidx) {
      if (!var->set_cur(latidx, 0) || !var->get(tmp, 1, grid.nlon)) {
        cerr << "ERROR: failed reading land mask" << endl;
        exit(1);
      }
      for (int lonidx = 0; lonidx < grid.nlon; ++lonidx)
        (*mask)(lonidx, latidx) = (tmp[lonidx] == 1);
    }
  } else {
    float tmp[grid.nlon];
    for (int latidx = 0; latidx < grid.nlat; ++latidx) {
      if (!var->set_cur(latidx, 0) || !var->get(tmp, 1, grid.nlon)) {
        cerr << "ERROR: failed reading land mask" << endl;
        exit(1);
      }
      for (int lonidx = 0; lonidx < grid.nlon; ++lonidx)
        (*mask)(lonidx, latidx) = (tmp[lonidx] > 25.0);
    }
  }
}


//----------------------------------------------------------------------
//
//  read_soil_type
//
//  Read soil type field from netCDF file.
//

static void read_soil_type(IntArray *soil_type, string soil_type_file,
                           Grid &grid, string soil_type_var)
{
  // Open the soil type file.

  NcFile nc(soil_type_file.c_str());
  if (!nc.is_valid()) {
    cerr << "ERROR: could not open soil type file " << soil_type_file << endl;
    exit(1);
  }


  // Access the soil type variable and check the dimensions are correct.

  NcVar *var = nc.get_var(soil_type_var.c_str());
  NcDim *lat_dim = nc.get_dim("lat");
  NcDim *lon_dim = nc.get_dim("lon");
  if (!var) {
    cerr << "ERROR: soil type variable '" << soil_type_var
         << "' not found in soil type file " << soil_type_file << endl;
    exit(1);
  }
  if (var->num_dims() != 2 ||
      var->get_dim(0) != lat_dim || var->get_dim(1) != lon_dim) {
    cerr << "ERROR: soil type variable '" << soil_type_var
         << "' has invalid dimensions" << endl;
    exit(1);
  }


  // Get the soil type data.

  for (int latidx = 0; latidx < grid.nlat; ++latidx) {
    if (!var->set_cur(latidx, 0) ||
        !var->get((*soil_type)(latidx), 1, grid.nlon)) {
      cerr << "ERROR: failed reading soil type" << endl;
      exit(1);
    }
  }
}


//----------------------------------------------------------------------
// 
//  read_lightn
//
//  Read lightening data from netCDF file.
//

static void read_lightn(Array *lightn, string lightn_file,
                         Grid &grid, string lightn_var)
{
  // Open the lightn file.

  NcFile nc(lightn_file.c_str());
  if (!nc.is_valid()) {
    cerr << "ERROR: could not open lightning file "
         << lightn_file << endl;
    exit(1);
  }


  // Access the lightening variable and check the dimensions
  // are correct.

  NcVar *var = nc.get_var(lightn_var.c_str());
  NcDim *lat_dim = nc.get_dim("lat");
  NcDim *lon_dim = nc.get_dim("lon");
  if (!var) {
    cerr << "ERROR: lightening  variable '" << lightn_var
         << "' not found in file " << lightn_file << endl;
    exit(1);
  }
  if (var->num_dims() != 2 ||
      var->get_dim(0) != lat_dim || var->get_dim(1) != lon_dim) {
    cerr << "ERROR: lightening variable '" << lightn_var
         << "' has invalid dimensions" << endl;
    exit(1);
  }


  // Get the lightening data.

  for (int latidx = 0; latidx < grid.nlat; ++latidx) {
    if (!var->set_cur(latidx, 0) ||
        !var->get((*lightn)(latidx), 1, grid.nlon)) {
      cerr << "ERROR: failed reading lightning" << endl;
      exit(1);
    }
  }
}


//----------------------------------------------------------------------
//
//  read_a_nd
//
//  Read human ignition data from netCDF file.
//

static void read_a_nd(Array *a_nd, string a_nd_file,
                      Grid &grid, string a_nd_var)
{
  // Open the human ignition data file.

  NcFile nc(a_nd_file.c_str());
  if (!nc.is_valid()) {
    cerr << "ERROR: could not open human ignition data file "
         << a_nd_file << endl;
    exit(1);
  }


  // Access the human ignition data variable and check the dimensions
  // are correct.

  NcVar *var = nc.get_var(a_nd_var.c_str());
  NcDim *lat_dim = nc.get_dim("lat");
  NcDim *lon_dim = nc.get_dim("lon");
  if (!var) {
    cerr << "ERROR: human ignition data variable '" << a_nd_var
         << "' not found in file " << a_nd_file << endl;
    exit(1);
  }
  if (var->num_dims() != 2 ||
      var->get_dim(0) != lat_dim || var->get_dim(1) != lon_dim) {
    cerr << "ERROR: human ignition data variable '" << a_nd_var
         << "' has invalid dimensions" << endl;
    exit(1);
  }


  // Get the human ignition data data.

  for (int latidx = 0; latidx < grid.nlat; ++latidx) {
    if (!var->set_cur(latidx, 0) ||
        !var->get((*a_nd)(latidx), 1, grid.nlon)) {
      cerr << "ERROR: failed reading human ignition data" << endl;
      exit(1);
    }
  }
}


//----------------------------------------------------------------------
//
//  read_co2
//
//  Read co2 concentration data from netCDF file.
//

static int read_co2(string var_name, string file_name, float **data)
{
  // Open the human ignition data file.

  NcFile nc(file_name.c_str());
  if (!nc.is_valid()) {
    cerr << "ERROR: could not open co2 concentration data file "
         << file_name << endl;
    exit(1);
  }


  // Access the co2 concentration data variable and check the dimensions
  // are correct.

  NcVar *var = nc.get_var(var_name.c_str());
  NcDim *time_dim = nc.get_dim("time");
  if (!var) {
    cerr << "ERROR: co2 concentration data variable '" << var_name
         << "' not found in file " << file_name << endl;
    exit(1);
  }
  if (var->num_dims() != 1 || var->get_dim(0) != time_dim) {
    cerr << "ERROR: co2 concentration data variable '" << var_name
         << "' has invalid dimensions" << endl;
    exit(1);
  }


  // Get the co2 concentration data data.
  *data = new float[time_dim->size()];
  if (!var->set_cur(long(0)) ||
      !var->get(*data, time_dim->size())) {
      cerr << "ERROR: failed reading co2 concentration data" << endl;
      exit(1);
  }
  return time_dim->size();
}


//----------------------------------------------------------------------
//
//  check_var_grid
//
//  Check grid definition for a climate input variable.  Each variable
//  should have three dimensions, time, lat and lon, and the lat and
//  lon dimensions should match the overall grid.
//

static void check_var_grid(NcVar *var, NcFile *nc,
                           string var_name, string file_name)
{
  if (var->num_dims() != 3) {
    cerr << "ERROR: bad dimensions for '" << var_name
         << "' variable in netCDF file " << file_name << endl;
    exit(1);
  }
  NcDim *time_dim = nc->get_dim("time");
  NcDim *lat_dim = nc->get_dim("lat");
  NcDim *lon_dim = nc->get_dim("lon");
  if (!time_dim || !lat_dim || !lon_dim) {
    cerr << "ERROR: bad dimensions in netCDF file " << file_name << endl;
    exit(1);
  }

  Grid tst_grid(file_name);
  if (var->get_dim(0) != time_dim || var->get_dim(1) != lat_dim ||
      var->get_dim(2) != lon_dim || tst_grid != grid) {
    cerr << "ERROR: bad dimensions for '" << var_name
         << "' variable in netCDF file " << file_name << endl;
    exit(1);
  }
}


//----------------------------------------------------------------------
//
//  open_output_file
//
//  Open and set up the output files.
//

static void open_output_files(string file_name, string directory_name)
{
  // Create data variable data structures.

  out_var_count = params.output_vars.size();
  out_vars = new LPJOutputVariable* [out_var_count];
  int var_idx = 0;
  for (vector<RunParameters::OutVariable>::iterator it =
         params.output_vars.begin();
       it != params.output_vars.end(); ++it) {
    string name = it->name;
    LPJVariable *def = 0;
    for (int idx = 0; idx < sizeof(LPJ_VARIABLES) / sizeof(LPJVariable);
         ++idx)
      if (name == LPJ_VARIABLES[idx].name) {
        def = &LPJ_VARIABLES[idx];
        break;
      }
    if (!def) { cerr << "INTERNAL ERROR #1" << endl;  exit(1); }
    out_vars[var_idx++] = new LPJOutputVariable(def, *it);
  }


  // Allocate space for output files.

  out_nc = new NcFile*[out_year_count];


  // Open netCDF output files.

  for (int idx = 0; idx < out_year_count; ++idx) {
    // Sort out filenames: the input string should just be a name
    // without any file extension.  Here we create an output filename
    // incorporating the index of the file.

    char index_str[32];
    sprintf(index_str, "%d", out_years[idx]);
    string out_file_name;
    if( directory_name.empty() ) {
      out_file_name  = file_name + "-" + index_str + ".nc";
    }
    else if( directory_name[directory_name.size() - 1] == '/' ) {
      out_file_name  = directory_name + file_name + "-" + index_str + ".nc";
    }
    else {
      out_file_name  = directory_name + "/" + file_name + "-" + index_str + ".nc";
    }

    NcFile *nc = new NcFile(out_file_name.c_str(), NcFile::Replace);
    out_nc[idx] = nc;
    if (!out_nc[idx]->is_valid()) {
      cerr << "ERROR: failed to open netCDF file " << out_file_name << endl;
      exit(1);
    }


    // Create dimensions.

    lat_dim = nc->add_dim("lat", grid.nlat);
    lon_dim = nc->add_dim("lon", grid.nlon);
    pft_dim = nc->add_dim("pft", NPFT);
    if (month_dim_needed) month_dim = nc->add_dim("month", NMON);
    if (day_dim_needed) day_dim = nc->add_dim("day", NDAYS);
    if (trace_dim_needed) trace_dim = nc->add_dim("trace", NTRACE);
    if (co2_dim_needed) co2_dim = nc->add_dim("co2", NCO2);

    // Create coordinate variables.

    NcVar *lat_var = nc->add_var("lat", ncFloat, lat_dim);
    lat_var->add_att("long_name", "Latitude");
    lat_var->add_att("units", "degrees_north");
    NcVar *lon_var = nc->add_var("lon", ncFloat, lon_dim);
    lon_var->add_att("long_name", "Longitude");
    lon_var->add_att("units", "degrees_east");
    NcVar *pft_var = nc->add_var("pft", ncInt, pft_dim);
    pft_var->add_att("long_name", "Plant Functional Type");
    NcVar *time_var = nc->add_var("time", ncInt);
    time_var->add_att("long_name", "Time");
    time_var->add_att("units", "years");
    NcVar *month_var, *day_var, *trace_var;
    NcVar *co2_var;
    if (month_dim_needed) {
      month_var = nc->add_var("month", ncInt, month_dim);
      month_var->add_att("long_name", "Month of year");
      month_var->add_att("units", "months");
    }
    if (day_dim_needed) {
      day_var = nc->add_var("day", ncInt, day_dim);
      day_var->add_att("long_name", "Day of year");
      day_var->add_att("units", "days");
    }
    if (trace_dim_needed) {
      trace_var = nc->add_var("trace", ncInt, trace_dim);
      trace_var->add_att("long_name", "Trace gas index");
    }
    if (co2_dim_needed) {
      co2_var = nc->add_var("co2", ncInt, co2_dim);
      co2_var->add_att("long_name", "CO2 concentration");
    }

    // Create data variables.

    for (int var_idx = 0; var_idx < out_var_count; ++var_idx)
      out_vars[var_idx]->set_nc(idx, *nc);


    // Write values for coordinate variables.
    
    lat_var->put(&grid.lats[0], grid.nlat);
    lon_var->put(&grid.lons[0], grid.nlon);
    int *pfts = new int[NPFT];
    for (int pft = 0; pft < NPFT; ++pft) pfts[pft] = pft + 1;
    pft_var->put(&pfts[0], NPFT);
    delete [] pfts;
    if (month_dim_needed) {
      int months[NMON] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
      month_var->put(&months[0], NMON);
    }
    if (day_dim_needed) {
      int days[NDAYS];
      for (int day = 0; day < NDAYS; ++day) days[day] = day + 1;
      day_var->put(&days[0], NDAYS);
    }
    if (trace_dim_needed) {
      int traces[NTRACE];
      for (int tr = 0; tr < NTRACE; ++tr) traces[tr] = tr + 1;
      trace_var->put(&traces[0], NTRACE);
    }
    if (co2_dim_needed) {
      co2_var->put(&co2_dim_values[0], NCO2);
    }
    int years[params.run_years];
    for (int yr = 0; yr < params.run_years; ++yr) years[yr] = yr + 1;
    time_var->put(&years[0],params.run_years);

    //    time_var->put(&out_years[idx]);

    nc->sync();
  }


}


//----------------------------------------------------------------------
//
//  handle_output_record
//
//  Handle averaging of data passed from a single call to outannual.
//

static void handle_output_record(string name, float *data)
{
  for (int idx = 0; idx < out_var_count; ++idx) {
    if (name == out_vars[idx]->def->name) {
      // Found the entry, so accumulate the data.
      for (int buf = 0; buf < out_vars[idx]->buff_size; ++buf) {
        if (out_vars[idx]->avg_buff)
          out_vars[idx]->avg_buff[buf] += data[buf];
        if (out_vars[idx]->sdev_buff)
          out_vars[idx]->sdev_buff[buf] += data[buf] * data[buf];
        if (out_vars[idx]->min_buff)
          if (data[buf] < out_vars[idx]->min_buff[buf])
            out_vars[idx]->min_buff[buf] = data[buf];
        if (out_vars[idx]->max_buff)
          if (data[buf] > out_vars[idx]->max_buff[buf])
            out_vars[idx]->max_buff[buf] = data[buf];
      }
      ++out_vars[idx]->accum_years;
    }
  }
}


//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
