
 

c
c CODE VERSION OF JULY 12, 2000.
c
c History of change:
c     12. 7.2000    W. Lucht baselined this version (10 pft, permafrost)
c     24.11.2000    W. Lucht and S. Venevsky fixed a bug that caused
c                   a rarely occurring variable overflow in the permafrost
c                   routine with the consequence of overflows in most other
c                   parameters. This situation occured when the
c                   condition maxthaw_depth-dthaw_depth(d).eq.d1 occurred.
c                   Now when this occurs, dthaw_depth(d) is slightly decreased
c                   to avoid the sitation.
c
c     ADD FUTURE CHANGES HERE
c     06.12.2000    Output additional monthly variables for the O18 analyses
c                   (GPP, ci/ca, individual respiration terms, soil temp)
c     23.05.2001    Output additional monthly variables (mauw1(1:12),mauw2(1:12)))
c     28.10.2002    Output additional annual and monthly variables (intc,aep,pet);
c                   !! completely altered waterbalance subroutine:
c                   interception and soil evaporation included; supply area-weighted;
c                   revised consideration of permafrost in water routine;
c                   demand function now hyperbolic, gm and alpham changed;
c                   percolation now occurs only on days with precipitation !!
c
c                   new subroutine "prdaily" (weather generator)
c
c                   changes in PFT parameters 1,3,6,9,10,12;
c                   new PFT parameters 35,36 (Emax,intc);
c
c                   alphaa=0.5;
c
c                   tsnow changed from -2?°C to 0?°C;
c                   snow melt: km changed to 1.5+0.007*dprec(d)
c starting 28.4.2004 KIRSTEN: implementation of SPITFIRE plus additional output
c                            and information from other routines
c
c DM starting june 2007 D. Marcadet for Y. Zhao
c DM                use of conditional compilation to merge several versions
c DM                intermediate save/restart for 2-steps simulation
c DM                several bugs corrected
c DM                clean some parts of the code
c
c
c
c 12/08-present    Doug Kelley:
c     As well as generral tidying, bug fixing (major ones listed below) and commenting:
c
c     12/08-01/09:    REDISTRIBUTION OF LIGHTENING ONTO WET/DRY DAYS
c           Lightening is redistributed amongst wet and dry days, with the ratio of
c         lightening falling on dry compared to wet days=(Wet days/dry days)**lc with
c         lc being a calibrated parameter.
c           Lightning occuring on wetdays is considered to be ineffective and is set to
c         zero.
c           Redistribution is conducted throiugh the new
c         daily2 subroutine, called in the main program.
c
c     01/09-06/09:    DAILY FUEL LOAD UPDATE
c           Switch from annual fuel load update to daily fuel load update. Changes to
c         fuel load, for each pft, is recorded annual from the start of establishment
c         through to fire in fuel_xhr_inc (where x is the fuel class) in each subroutine
c         where changes in fuel load occure, on the timestep that subroutine operates (so,
c         most changes are recorded annually, but decay in littersom, for example, is
c         currently recorded monthly). The change in then divided between the days in that
c         timestep (so for annual subroutines, divided equally amongst all 365 days, on
c         monthly, divided amongst the 31, 30 or 28 days in that month).
c           1hour fuel also records if the change is posatve or negative (through
c         fuel_1hr_inc_pos or ""_neg).
c              ______________________________________________
c             | subroutine    | fuel effected     | timestep |
c             |_______________|___________________|__________|
c             | reproduction  | 1hr (posative and |          |
c             |               |  negative)        |          |
c             |---------------|-------------------|----------|
c             | turnover      | All               |          |
c             |---------------|-------------------|----------|
c             | littersom     | All               | Monthly  |
c             |---------------|-------------------|----------|
c             | kill_pft      | All               |          |
c             |---------------|-------------------|----------|
c             | allocation    | All               |          |
c             |---------------|-------------------|----------|
c             | light         | All               |          |
c             |---------------|-------------------|----------|
c             | mortality     | All               |          |
c             |---------------|-------------------|----------|
c             | fire          | All - effects     | Daily    |
c             |               |   fire dynmaics   |          |
c             |----------------------------------------------|
c             | fuel_xhr_inc(pos/neg) set to zero ready to   |
c             | record the next years change                 |
c             |----------------------------------------------|
c             | establishment | All               |          |
c             |_______________|___________________|__________|
c
c           In SPITFIRE, the recorded changes are taken off the fuel load (as the actual
c         fuel load passed into the fire subroutine contains all these incraments) for fuels
c         from 10 to 1000hr, and negative 1hr are added daily, at the beginning of the daily
c         timestep. Posative changes of 1hr fuels are only added on days with a negative
c         change in phenology (in dephen)) for that pft. This is calculted in the
c         fuel_1hr_redist before the daily timestep within SPITFIRE.
c           Fuel is also reduced each day by the a factor (1-previous days burning), and
c         has added to it the daily fuel update from plant fire mortaility
c        (dfuel_update_xhr, which is calculated as the daily change in the annual varibale,
c        fuel_update_xhr)
c           These changes are done on fuel_xhr_left, which now stroes the daily fuel load
c         per, and replaces fuel_xhr_0 (the annual fuel load) in:
c               @ fuel bed depth ration calcualtions (ratio_fbd);
c               @ fuel consumption (fuel_consumption subroutine).
c                                              All other fuel effects on fire are still
c         carried out through the now modified fuel_xhr_total, which is the sum of
c         fuel_xhr_left over all tree pfts.
C           Before any fuel changes are calculted, fuel_xhr_inc(_pos/neg) are /0.45, and
c         fuel_left_xhr is * or / 0.45 at various points to conserve units.
c         
c
c     03/09:          ERROR FILE
c           If a minor error occurs that does not shut down the model (such as invalid soil
c         type), the error is now written in an error file, fort.10, in same directory as
c         the model operates.
c
c     03/09-06/09:    MASKING FIRE IN FUEL LIMITED AREAS
c           In SPITFIRE in fire subroutine, fuel limited areas are not allowd to burn.
c         'Masking' occurs in daily timestep, after the days fuel load has been calculated.
c         If the fuel load is less then then a threshold of
c         fuel_1hr+fuel_10hr+livegrass*365<2000, then no fires occur and that day, and the 
c         loop daily loop is cycled round.
c
c     03/09-05/09:    MONTHLY FUEL I/0
C           Added monthly fuel array (mfuel_xhr_total & mlivegrass) for output,
c         calculted in daily timestep in SPITFIRE section of fire subroutine.
c
c     04/09-05/09:    AGRICULTURE I/O
C           Switched crop and pasture from being read in through main code to
c         driver (currently lpjio.cpp), and made trasnient. Main now gets crop and pas value
c         through the drivers getclimate subroutine.
c
c     05/09:    DPHEN_CHANGE I/O
C           Included date of change of phen into output (in dphen_change). For each month,
c         records day of the year dphen 1st changes. A posative number reprents a posative
c         change in dephen, negative reprents negative. dphen_change is calculted in gpp.
c
c     05/09-06/05: (Stephen & Doug)    IMPLEMENT LAND USE CHANGE 
c
c	  AGRICULTURAL VERSION (CCMLP Grand Slam expt.)  
c  	Extra input variables for years 1901-1998 (=98 years) '96-98' no landuse conversion
c  	grid cells cleared updated annually
c  	grid cells abandoned updated annually
c  	rap (relative agricultural productivity) for agricultural areas
c  	updated annually              
c
c  	The model runs as normal assuming natural vegetation. When the
c  	gridcell changes state then the productivity of the natural 
c  	vegetation is taken and the RAP value to receive the agricultural
c  	productivity.
c
c
c
c  	main program calls the ioprogram to find out whether in this
c  	year the gridcell is abandoned or cleared and its state
c  	agricultural or not. In addition this years RAP is supplied.
c  	A new subroutine is constructed (agric) which updates the grid cell
c  	pools if clearance or abandonment in this year. 
c  	The decomposition routine is altered to account for agricultural
c  	above and below ground litter. 
c  	A second subroutine is constructed (agriprod) which calculates
c  	agricultural production and product decompostion (prod1, prod10,
c  	prod100).
c
c  	Send following information back to the io program:
c  	conversion flux
c  	product flux
c  	pools: agprod1, prod10, prod100
c   
c  	Fire routine is kept in based on the natural conditions. Natural litter
c  	is also calculated and used in the fire routine, but not included
c  	as output in agricultural gridcell C-balance.
c
c     06/09:          REDUCED ESTABILMENT FOR TREE PFTS ON PASTURE LAND
C           In estabilment subroutine, establishment of tree pft is reduced by the
c         proportion of pastural land, s.t.:
c              estab_grid=estab_rate*(1.0-fpc_tree_total)
c           -->estab_grid=estab_rate*(1.0-fpc_tree_total)*(1-pas)
c  
c     06/09:          FIRE INENSITY THRESHOLD FIX
c           Fix of fire intensity threshold so that the thrshold works. Intensity below
c         threshold will extinguish fires. Currently set to zero, but that may change...
*******************************************************************************
c *                                                                           *
c *                                  L P J                                    *
c *             Lund-Potsdam-Jena Dynamic Global Vegetation Model             *
c *                                                                           *
c *                  Principal Author: Stephen Sitch (1,2)                    *
c *       Associate Author and Artistic Director: Benjamin Smith (1,3)        *
c *                                                                           *
c * (1) Climate Impacts Group, Department of Ecology, Plant Ecology,          *
c *     Lund University, Ekologihuset, Lund, S-22362 Sweden                   *
c * (2) Potsdam Institute for Climate Impact Research, P.O. Box 60 12 03,     *
c *     D-14412 Potsdam, Germany                                              *
c * (3) Max Planck Institute for Biogeochemistry, P.O. Box 100164, D-07701,   *
c *     Jena, Germany                                                         *
c *                                                                           *
c * Correspondence should be directed to the principal author,                *
c * Potsdam Institute for Climate Impact Research,                            *
c * Fax +49 331 288 2600, E-mail Stephen.Sitch@pik-potsdam.de                 *
c *****************************************************************************

c -----------------------------------------------------------------------------
c                  TECHNICAL INFORMATION ABOUT THIS FILE
c -----------------------------------------------------------------------------

c File name:     lpjmain.f
c Language:      FORTRAN 77 with extensions (do ... enddo, do while ... enddo,
c                long variable names)
c Compilation:   f77 -o lpj lpjmain.f lpjio.f
c                (assumes input/output module called 'lpjio.f', generates
c                executable file 'lpj')
c Version dated: 23/02/99


c *****************************************************************************
c *                                                                           *
c *                                INPUT/OUTPUT                               *
c *                                                                           *
c * The model source code comprises two modules: the main program (this file) *
c * and an input/output module. The two modules should be compiled and linked *
c * to produce a single executable file.  Normally, only the main program     *
c * will be distributed by the authors: IT IS EACH USER'S OWN RESPONSIBILITY  *
c * TO PROVIDE SUITABLE INPUT/OUTPUT CODE. The main program and input/output  *
c * module are linked through calls, in the main program, to six subroutines  *
c * in the input/output module.  The structure of the main program, including *
c * calls to input/output subroutines, is as follows:                         *
c *                                                                           *
c *       call initio                                                         *
c *       call getgrid                                                        *
c *       do while (dogridcell)  !loop through grid cells                     *
c *                      :                                                    *
c *         [ grid cell initialisation ]                                      *
c *                      :                                                    *
c *         call getclimate                                                   *
c *         do while (doyear)    !loop through simulation years               *
c *                     :                                                     *
c *           [ annual processes ]                                            *
c *                     :                                                     *
c *           call outannual                                                  *
c *           call getclimate                                                 *
c *         enddo                                                             *
c *         call outgrid                                                      *
c *         call getgrid                                                      *
c *       enddo                                                               *
c *       call termio                                                         *
c *                                                                           *
c * Specifications of input/output subroutines are given below. Declarations  *
c * and descriptions of argument variables are to be found in the             *
c * variable declarations section of the main program.                        *
c *                                                                           *
c * subroutine initio [no arguments]                                          *
c *   Initialisation of input/output, e.g. opening of files, run-time         *
c *   interaction with user.                                                  *
c *                                                                           *
c * subroutine getgrid(lat,soilcode,dogridcell)                               *
c *   Should return latitude and LPJ soil code (see subroutine                *
c *   soilparameters) for next grid cell. Logical 'dogridcell' should be set  *
c *   to .false. if no further grid cells remain to be processed, and .true.  *
c *   otherwise.                                                              *
c *                                                                           *
c * subroutine getclimate(year,mtemp,mprec,msun,co2,doyear)                   *
c *   Should return monthly temperature, precipitation and sunshine, and      *
c *   atmospheric CO2 concentration, for current grid cell and year of        *
c *   simulation. Year counter 'year' is updated in the main program and      *
c *   should not be modified within the subroutine. Logical 'doyear' should   *
c *   be set to .false. following the last simulation year, and .true.        *
c *   otherwise.                                                              *
c *                                                                           *
c * subroutine outannual(year,present,nind,lm_ind,rm_ind,sm_ind,hm_ind,       *
c *   fpc_grid,anpp,acflux_estab,litter_ag,litter_bg,cpool_fast,cpool_slow,   *
c *   arh,afire_frac,acflux_fire)                                             *
c *   Output subroutine called at the end of each simulation year. The        *
c *   suggested argument list here may be supplemented by other variables     *
c *   if desired (by modifying argument lists both in the subroutine header,  *
c *   and the call in the main program). None of the argument variables       *
c *   should normally be modified within the subroutine.                      *
c *                                                                           *
c * subroutine outgrid(year,present,nind,lm_ind,rm_ind,sm_ind,hm_ind,         *
c *   fpc_grid,anpp,acflux_estab,litter_ag,litter_bg,cpool_fast,cpool_slow,   *
c *   arh,afire_frac,acflux_fire)                                             *
c *   Output subroutine called at the end of the simulation for each grid     *
c *   cell. The suggested argument list here may be supplemented by other     *
c *   variables if desired (by modifying argument lists both in the           *
c *   subroutine header, and the call in the main program). None of the       *
c *   argument variables should normally be modified within the subroutine.   *
c *                                                                           *
c * subroutine termio [no arguments]                                          *
c *   Closing of files etc. at end of model run.                              *
c *                                                                           *
c * NOTES                                                                     *
c *                                                                           *
c *   1. All six subroutines must appear in the input/output module.          *
c *      However, subroutines that are not required (potentially any of the   *
c *      above, except getgrid and getclimate) may be empty, i.e. comprise    *
c *      only a subroutine header, declarations section, and the statements   *
c *      'return' and 'end'.                                                  *
c *                                                                           *
c *   2. It is intended that input/output should be handled independently of  *
c *      the main program code. This will normally require incorporation of   *
c *      one or more common blocks to provide data links between input/       *
c *      output subroutines.                                                  *
c *                                                                           *
c *   3. PFT state variables are defined only for PFTs present at the end of  *
c *      a given simulation year. Among-PFT totals may therefore be           *
c *      calculated within output subroutines using code similar to the       *
c *      following (in this case, for calculation of annual ecosystem NPP):   *
c *                                                                           *
c *      integer pft,npft                                                     *
c *      parameter (npft=13)                                                   *
c *      real anpp_total                                                      *
c *                                                                           *
c *      npp_total=0.0                                                        *
c *      do pft=1,npft                                                        *
c *        if (present(pft)) anpp_total=anpp_total+anpp(pft)                  *
c *      enddo                                                                *
c *                                                                           *
c *   4. Grid cell total vegetation carbon (gC/m2) may be calculated as       *
c *      follows:                                                             *
c *                                                                           *
c *      real vegc                                                            *
c *                                                                           *
c *      vegc=0.0                                                             *
c *      do pft=1,npft                                                        *
c *        if (present(pft)) vegc=vegc+                                       *
c *          (lm_ind(pft)+rm_ind(pft)+sm_ind(pft)+hm_ind(pft))*nind(pft)      *
c *      enddo                                                                *
c *                                                                           *
c *   5. Annual net ecosystem production (NEP, gC/m2) may be calculated as    *
c *      follows:                                                             *
c *                                                                           *
c *      anep = anpp_total + acflux_estab - arh - acflux_fire                 *
c *      where 'anpp_total' is ecosystem annual NPP, as calculated            *
c *      according to Note 3 (above).                                         *
c *                                                                           *
c *****************************************************************************


c -----------------------------------------------------------------------------
c                            SCIENTIFIC SUMMARY
c -----------------------------------------------------------------------------

c The model simulates vegetation dynamics, hydrology and soil organic matter
c dynamics on an area-averaged grid cell basis using a one-year time step.
c Input parameters are monthly mean air temperature, total precipitation and
c percentage of full sunshine, annual atmospheric CO2 concentration and soil
c texture class. The simulation for each grid cell begins from "bare ground",
c requiring a "spin up" (under non-transient climate) of c. 1000 years to
c develop equilibrium vegetation and soil structure.

c Vegetation is represented by a combination of the following plant functional
c types (PFTs)

c       1. tropical broadleaved evergreen tree
c       2. tropical broadleaved raingreen tree
c       3. temperate needleleaved evergreen tree
c       4. temperate broadleaved evergreen tree
c       5. temperate broadleaved summergreen tree
c       6. boreal needleleaved evergreen tree
c       7. boreal needleleaved summergreen tree
c       8. boreal broadleaved summergreen tree
c       9. C3 perennial grass
c      10. C4 perennial grass

c A minimum set of bioclimatic limits (minimum coldest month temperature for
c survival, maximum coldest month temperature for establishment, minimum
c growing degree days [GDD, on 5 degree base except for larches (2 deggrees)]
c for survival) constrain PFT
c occurrence and regeneration within grid cells. Woody PFTs (1-8 above) present
c within a given grid cell are represented by area-averages for the following
c state variables

c       nind        individual density (/m2)
c       lm_ind      individual leaf carbon mass (gC)
c       rm_ind      individual fine root carbon mass (gC)
c       sm_ind      individual sapwood carbon mass (gC)
c       hm_ind      individual heartwood carbon mass (gC)
c       crownarea   individual crown area (m2)

c Individual structure is ignored for grass PFTs (9, 10). Individual density is
c arbitrarily set to the value 1/m2, so that lm_ind and rm_ind represent grid
c cell area-averages (gC/m2). The variables sm_ind and hm_ind are undefined for
c grasses, while crownarea is set to the proportion of the grid cell not
c occupied by trees, and thereby available for colonisation by grasses (see
c below).

c Individual allometry for woody PFTs is defined by the following relations,
c which constrain allocation of annual PFT carbon increments to the three
c living biomass compartments (leaves, fine roots and sapwood):

c      (A) (leaf area) = latosa * (sapwood xs area)
c             (Pipe Model, Shinozaki et al. 1964a,b; Waring et al 1982)
c      (B) lm_ind = lmtorm * rm_ind
c      (C) height = allom2 * (stem diameter)**allom3
c      (D) crownarea = min (allom1 * (stem diameter)**reinickerp,
c                           crownarea_max)

c      where

c      (sapwood xs area) = lm_ind * sla / latosa
c      (stem diameter) = ( 4 * (sm_ind + hm_ind ) / wooddens / pi /
c                        allom2 )**( 1 / (2 + allom3) )
c      latosa, allom1, allom2, allom3, reinickerp, crownarea_max, sla and
c        wooddens are constants.
c      lmtorm is calculated annually based on root-zone soil water
c        availability.

c Grass PFT structure is based solely on equation (B) above.

c PFT proportional cover is represented by foliar projective cover (FPC), which
c in turn is a function of leaf area index (LAI). The relations are:

c      lai_ind  = lm_ind * sla / crownarea      (individual LAI)
c      fpc_ind  = 1 - exp (-0.5 * lai_ind)      (individual FPC)
c      fpc_grid = (crownarea * nind) * fpc_ind  (grid cell FPC)

c Competition for space among woody PFTs is assumed to constrain the sum of
c their grid cell FPCs to marginally less than 1. Establishment of new
c individuals (saplings) is confined to the area of the grid cell not already
c occupied by woody PFTs. The establishment rate per unit area is reduced by
c shading as total woody FPC approaches 1. Competition for light among trees
c and grasses causes grass PFT cover (represented by the variable crownarea) to
c be constrained to the proportion of the grid cell not occupied by trees.

c Carbon uptake through photosynthesis, plant evapotranspiration and soil water
c content (of an upper and lower soil layer) are calculated by a coupled
c photosynthesis and water balance model, closely based on that described by
c Haxeltine & Prentice (1996).  PFTs capture incoming radiation in proportion
c to their grid cell FPCs. The annual biomass increment is based on gross
c primary production (GPP), maintenance and growth respiration, reproduction
c costs, turnover of sapwood to heartwood, and production of leaf and root
c litter.

c Disturbance by fire may cause additional annual biomass losses.  The
c proportion of a grid cell burnt each year is based on litter density and
c moisture content. PFTs may differ in the flammability of their litter.

c Decomposition of litter and soil organic matter is based on temperature
c (Lloyd & Taylor 1994) and soil layer 1 moisture content (Foley 1995).
c Differential decomposition rates apply to the litter, 'fast' and 'slow'
c organic matter (carbon) pools.

c Additional documentation is provided by commenting throughout the source
c code.  References are cited where appropriate and listed at the end of this
c file.

c -----------------------------------------------------------------------------




#define CO2_THREE_VALUES

#define POPULATION_DENSITY_IN_FILE


      program lpj

      implicit none

c     PARAMETERS
      integer npft,npftpar,nsoilpar
        parameter (npft=13)             ! number of PFTs

        parameter (npftpar=59)         ! number of PFT parameters
      
        parameter (nsoilpar=7)         ! number of soil parameters
      integer nco2
        parameter (nco2=3)             ! number of C variables: C,13C,14C
      integer climbuf                  ! number of years over which climate is
        parameter (climbuf=20)         ! averaged to implement bioclimatic limits
      real allom1,allom2,allom3,allom4
        parameter (allom1=100.0)       ! parameters
        parameter (allom2=40.0)        ! in
        parameter (allom3=0.67)        ! allometric
        parameter (allom4=0.3)         ! equations
      real reinickerp
        parameter (reinickerp=1.6)     ! parameter in allometric equation
      real latosa                      ! ratio of leaf area to sapwood
        parameter (latosa=6.0e3)       ! cross-sectional area (Shinozaki et al
                                       ! 1964a,b)
      real wooddens
        parameter (wooddens=2.0e5)     ! wood density (gC/m3)

      real d1,d2                       ! depth of upper and lower buckets (mm)
        parameter (d1=500.0,d2=1000.0)
      real d_evap                      ! depth of evaporation soil layer (mm)
        parameter (d_evap=200.0)

c     INPUT PARAMETERS
      integer soilcode                 ! LPJ soil code (see subr soilparameters)
      real popden                      ! human population density
      real a_nd                        ! human-caused, potential ignitions
      real mlightn(1:12)               ! needed by this version of getclimate
      real dlightn(1:365), dlightn_control(1:365)              ! filled by interpolation
      REAL cgf(1:12)				   ! Doug 09/12: fraction of cloud-to-ground-lightning
      REAL fdry(1:12)                  ! Doug 04/13: fraction of lighting striking on a dry day 
      REAL lt_days(1:12)               ! Doug 04/13: fraction of dry days with lightning 
      real co2(1:nco2)                 ! atmospheric CO2 concentration (ppmv)
      real lat                         ! latitude (degrees +=N, -=S)
      real lon                         ! longditude (degrees +=E, -=W)
                                       ! Doug 04/09: agriculture variables changed from (720,360)
                                       ! array to scalar due to changeing the I/O method to match
                                       ! standard I/O for other variables
      real crop                        ! Proportion of cell containing cropland agriculture. 
      real pas                         ! Proportion of cell containing pasture agriculture.
      integer l,ll,day
      real ii,jj
      real mprec(1:12)                 ! monthly precipitation (mm)
      real aprec                       ! Doug 06/13: annual precipitation (mm)
      real mtemp_dmin(1:12)            ! monthly minimum temperature
      real mtemp_dmax(1:12)            ! monthly maximum temperature
      real mwet(1:12)                  ! monthly number of wet days  !DIETER
      real msun(1:12)                  ! monthly sunshine (% full sunshine)
      real mtemp(1:12)                 ! monthly temperature (deg C)
      real mwindsp(1:12)               ! monthly average wind speed (m/s)

c     LOCAL VARIABLES
      integer year                     ! simulation year
      logical dogridcell               ! whether more gridcell(s)
      logical doyear                   ! whether more year(s)
      logical estab(1:npft)            ! whether PFT within bioclimatic limits
                                       ! for establishment
      logical evergreen(1:npft)        ! whether PFT is evergreen
      logical needle(1:npft)           ! whether PFT is needleleaved
                                       ! (alternative: broadleaved)
      logical boreal(1:npft)           ! whether PFT is boreal
      logical present(1:npft)          ! whether PFT present in gridcell
      logical raingreen(1:npft)        ! whether PFT is raingreen
      logical summergreen(1:npft)      ! whether PFT is summergreen
      logical survive(1:npft)          ! whether PFT within bioclimatic limits
                                       ! for survival
      logical tree(1:npft)             ! whether PFT is a tree (alternative:
                                       ! grass)
      real agpp(1:npft,1:nco2)         ! annual gridcell GPP (gC/m2)
      real alresp(1:npft,1:nco2)       ! annual gridcell leaf respiration (gC/m2)
      real anpp(1:npft,1:nco2)         ! annual gridcell NPP (gC/m2)
      REAL anpp_grid(1:npft)           ! Doug 08/09: annual npp per pft (same as above,  but cutting out the nc dimension, smaller for ouputs)
                                       !     name of array probabily not right. Change _grid bit.
      real arh(1:nco2)                 ! annual heterotrophic respiration (gC/m2)
      REAL arh_grid                    ! Doug 08/09: as above, cutting out nc dimension, smaller for ouputs
      real arunoff                     ! total annual runoff (mm)
      real arunoff_drain               ! annual drainage (layer 2 runoff) (mm)
      real arunoff_surf                ! annual surface (layer 1) runoff (mm)
      real bm_inc(1:npft,1:nco2)       ! annual biomass increment (gC/m2)
      real acflux_estab(1:nco2)        ! annual biomass increment due to
                                       ! establishment (gC/m2)
      real cpool_fast(1:nco2)          ! fast-decomposing soil C pool (gC/m2)
      real cpool_slow(1:nco2)          ! slow-decomposing soil C pool (gC/m2)
      real crownarea(1:npft)           ! crown area (m2)
      real dmelt(1:365)                ! daily snowmelt (mm)
      real dpet(1:365)                 ! daily potential evapotranspiration (mm)
      real dphen(1:365,1:npft)         ! net daily leaf-on fraction
      real dphen_t(1:365,1:npft)       ! daily leaf-on fraction due to
                                       ! temperature phenology
      real dphen_w(1:365,1:npft)       ! daily leaf-on fraction due to
                                       ! drought phenology
      REAL dphen_change(1:12,1:npft)   ! Doug 05/09: record day number in the month
                                         ! and for each pft when dphen (1st) changes.
                                         ! Stored as 0 if no change, +day of year for
                                         ! postaive change, -day for negative change
      REAL dprec(1:365),dprec_out(1:365)! daily precipitation (mm)
      real dsun(1:365)                 ! daily sunshine (% full sunshine)
      real dtemp(1:365)                ! daily temperature (deg C)
      real dtemp_soil(1:365)           ! daily soil temperature (deg C)
c     additions by Kirsten
      real dtemp_min(1:365)            ! daily temperature minimum (deg C)
      real dtemp_max(1:365)            ! daily temperature maximum (deg C)
      real dwindsp(1:365)              ! daily windspeed (m/s)

      real dayl(1:365)
      real dw1(1:365)                  ! daily soil layer 1 water content
      real dwscal365(1:npft)           ! daily water scalar for day 365
      real fpc_grid(1:npft)            ! gridcell foliar projective cover (FPC)
      real fpc_inc(1:npft)             ! increment (if +ve) in FPC since last
                                       ! year
      real gdd(1:npft)                 ! current-year growing degree days
      real gdd_buf(1:npft,1:climbuf)   ! buffer to store 'climbuf' years of GDD
                                       ! totals
      real gdd20(1:npft)               ! 20-year average growing degree days
      REAL gdd_grid                    ! Doug 07/09: gdd with base defined below for each grid (i.e, not just for pfts). Used for output only
      REAL gddbase                     ! Doug 08/09
          pARAMETER(gddbase=5)         ! Doug 08/09
      REAL alpha_ws                    ! Doug 07/09: Bioclimatic indexfor moiture availability indicating water stress=sum(actual)/sum(equilibrium evapotranspiration) over the year
      REAL deet(1:365)                 ! Doug 08/09: Equilibrium evapotrainspiration for calculating alpha_ws
      real height(1:npft)              ! tree height (m)
      real height_class(0:4,1:npft)
      real dbh(1:npft)                 ! stem diameter per PFT, Kirsten.
      REAL dbh_class(0:4,1:npft)       ! Doug 02/13: dbh for each height class
      REAL BTparam1(1:npft,1:3)        ! Doug 02/13: lower, medium and upper bound describing 1-50-99% quantile p1 in Bark thickness equation BT=p1+p2*DBH
      REAL BTparam2(1:npft,1:3)        ! Doug 02/13: "" for p2
      REAL BTmode0(1:npft,1:2)         ! Doug 02/13: 50% quatile starting point incase of establishment or tree death.
      real tau_c(0:4,1:npft)           ! critical time to cambial kill, Kirsten.
      real cl_t(0:4,1:npft)            ! crown length per height class, Kirsten.
      real hm_ind(1:npft,1:nco2)       ! individual heartwood mass (gC)
      real hm_sapl(1:npft,1:nco2)      ! initial (sapling) heartwood mass (gC/m2)
      real k_fast_ave                  ! running average k_fast for subroutine
                                       ! littersom
      real k_slow_ave                  ! running average k_slow for subroutine
                                       ! littersom
      real lai_ind(1:npft)             ! individual leaf area index
      REAL litter_ag(1:npft,1:nco2)	   ! gridcell above-ground litter (gC/m2)
									   ! Doug 11/12: split above-ground litter into grass and wood
									   ! components to allow for different decomosition rates
      REAL litter_ag_leaf(1:npft,1:nco2) !Doug 11/12: grass and leaf litter
      REAL litter_ag_wood(1:npft,1:nco2) !Doug 11/12: wood litter
      real litter_bg(1:npft,1:nco2)    ! gridcell below-ground litter (gC/m2)
      real litter_decom_ave(1:nco2)    ! running average litter_decom for
                                       ! subroutine littersom
      real lm_ind(1:npft,1:nco2)       ! individual leaf mass (gC)
      REAL lm_inc(1:npft)              ! Doug 06/13: +ive incremment of leaf matter
      real lm_sapl(1:npft,1:nco2)      ! initial (sapling) leaf mass (gC/m2)
      real mdayl(1:12)                 ! mid-month daylength
      real mgpp(1:12,1:npft,1:nco2)    ! monthly grid cell GPP (gC/m2)
      real mlresp(1:12,1:npft,1:nco2)  ! monthly grid cell leaf respiration
                                       ! (gC/m2)
      real mnpp(1:12,1:npft,1:nco2)    ! monthly gridcell NPP (gC/m2)
      REAL mnpp_grid(1:12)             ! Doug 06/09: A bit of a cheat that really should be sorted out.
                                         ! Bits of mnpp (i.e. just actaul carbon, and for all cell, not
                                         ! per pft) that we actually want for ouput.
      real mpar_day(1:12)              ! mid-monthly PAR flux (J/m2/day)
      real mrh(1:12,1:nco2)            ! monthly heterotrophic respiration
                                       ! (gC/m2)
      real mrunoff(1:12)               ! total monthly runoff (mm)
      real mtemp_min_buf(1:climbuf)    ! buffer to store 'climbuf' years of
                                       ! coldest month temperatures
      real mtemp_min20                 ! 20-year average minimum monthly
                                       ! temperature (deg C)
      real mtemp_old(1:12)             ! last year's monthly temperatures (deg C)
      real mtemp_soil(1:12)            ! monthly soil temperature (deg C)
      real mauw1(1:12)                 ! monthly average values uw1, soil layer 1 water content (mm)
      real mauw2(1:12)                 ! monthly average values uw2, soil layer 2 water content (mm)
      real nind(1:npft)                ! gridcell individual density (indiv/m2)
      real pftpar(1:npft,1:npftpar)    ! PFT parameters
      real rm_ind(1:npft,1:nco2)       ! individual fine root mass (gC)
      real rm_sapl(1:npft,1:nco2)      ! initial (sapling) fine root mass (gC/m2)
      real sla(1:npft)                 ! PFT specific leaf area (m2/gC)
      real sm_ind(1:npft,1:nco2)       ! individual sapwood mass (gC)
      real sm_sapl(1:npft,1:nco2)      ! initial (sapling) sapwood mass (gC/m2)
      real snowpack                    ! storage of precipitation as snow (mm)
      real soilpar(1:nsoilpar)         ! soil parameters
      real turnover_ind(1:npft)        ! total turnover of living biomass per
                                       ! individual (gC)
      real w(1:2)                      ! soil layer 1 and 2 water content
                                       ! (fraction of available water holding
                                       ! capacity)
      real w_t(1:2)                    ! dto.,frozen+liquid water
      real w_ep                        ! fraction of available water holding capiacity evap-layer
      real wscal(1:npft)               ! mean daily water scalar (among leaf-on
                                       ! days) (0-1 scale)
      real aaet                        ! annual actual evapotranspiration (mm/year)
      real aintc                       ! annual interception evaporation (mm)
      real aaep                        ! annual actual evaporation (mm/year)
      real apet                        ! annual potential evaporation (mm/year) /*DIETER*/
      real apet_grid                   ! annual PET as grid cell means, incl. demand (mm)
      real mtemp_max                   ! warmest-month temperature (deg C)
c DM      real heatstress(1:npft)          ! reduction in individual density (& establishment)
                                       ! due to heat induced mortality  (indiv/m2)

      real u(1:12),v(1:12)
      real sinelat,cosinelat
      real sinehh(1:12),hh(1:12),qo(1:12)
      integer pft,m,i,d
      integer nc                       ! nc is added for nco2
      integer leafondays(1:npft),leafoffdays(1:npft)
      logical leafon(1:npft)
c     mpar=monthly grid par (mpar_day(m)*ndays(month))
c     mapar=monthly grid apar (sum of mpar_day(m)*fpar(m)*alphaa over pfts)
c     mphen=monthly average pheno-state (sum_d(m) of dphen(d,pft)/ndays(m))
      real mphen(1:12,1:npft)
      real mpar(1:12)
      real mapar(1:12)
c DM      real xmphen(1:12*npft)           ! work buffer for driver interface
c DM      real xmnpp(1:12*npft)            ! work buffer for driver interface
      real anpp_add(1:nco2)            ! annual gridc. NPP of killed PFTs(gC/m2)
      real mnpp_add(1:12,1:nco2)       ! monthly gridc. NPP of killed PFTs(gC/m2)
      real maet(1:12)                  ! monthly actual transpiration (mm)
      real mintc(1:12)                 ! monthly interception loss (mm)
      real maep(1:12)                  ! monhtly actual evaporation (mm)
      real mpet(1:12)                  ! monthly potential evaporation (mm)
      real mpet2(1:12)                 ! monthly potential evaporation (mm) for output file
      real mmelt(1:12)                 ! monthly snowmelt (mm)
      real mpet_grid(1:12)             ! monthly PET as grid cell means, incl. demand (mm)
      real mtemp_max20                 ! 20-year average maximum monthly
      real mtemp_max_buf(1:climbuf)    ! buffer to store 'climbuf' years of
                                       ! coldest month temperatures
                                       ! temperature (deg C)
c DM      real gminp(1:npft)

c     variables for fire simulation Kirsten
      real acflux_fire(1:nco2)         ! C flux to atmosphere due to fire (gC/m2)
      REAL acflux_fire_grid            ! Doug 08/09: as above, but with no nc dimension. Smaller for outputs
      real mcflux_fire(1:12,1:nco2)    ! C flux to atmosphere due to fire (gC/m?²)
      real an_fseason
c DM      real fire_length
      real afire_frac                  ! fraction of gridcell burnt this year
      REAL afire_frac_afap_old		   ! Doug 09/12: Correccting error. Lastyears fire_frac fed through each year
      REAL fuel_1hr_leaf(1:npft,1:nco2)! 1hr dead fuel: dead grass,"cured" grass,shed tree leaves, small twigs
      REAL fuel_1hr_wood(1:npft,1:nco2)
      real fuel_10hr(1:npft,1:nco2)    ! 10hr dead fuel: large twigs
      real fuel_100hr(1:npft,1:nco2)   ! 100hr dead fuel: small branches
      real fuel_1000hr(1:npft,1:nco2)  ! 1000hr dead fuel: logs, bole, large branches
      REAL pfuel_limit				   ! Doug 12/12: proportion of days where fire is limited by fuel
	  
      REAL fuel_1hr_leaf_inc_pos(1:npft,1:nco2,1:12)	
      REAL fuel_1hr_wood_inc_pos(1:npft,1:nco2,1:12)	
      REAL fuel_1hr_leaf_inc_neg(1:npft,1:nco2,1:12)	!Doug 01/09 MFL: to record changes in each
      REAL fuel_1hr_wood_inc_neg(1:npft,1:nco2,1:12)	!Doug 01/09 MFL: to record changes in each
      REAL fuel_10hr_inc(1:npft,1:nco2,1:12)	!hr fuel load for each month in the
      REAL fuel_100hr_inc(1:npft,1:nco2,1:12)	!yearly to allow monthly difference to
      REAL fuel_1000hr_inc(1:npft,1:nco2,1:12)	!considered in fire subroutine

      REAL mfuel_1hr_total(1:12)                !Doug 03/09: monthly fuel load (output only)
      REAL mfuel_1hr_leaf_total(1:12)           !Doug 11/12: monthly grass and leaf fuel load (output only)
      REAL mfuel_1hr_wood_total(1:12)           !Doug 11/12: monthly wood fuel load (output only)
      REAL mfuel_10hr_total(1:12)
      REAL mfuel_100hr_total(1:12)
      REAL mfuel_1000hr_total(1:12)
      REAL mlivegrass(1:12)                     !Doug 05/09: monthly livegrass (for ouput only)

      REAL fuel_1hr_del(1:npft,1:nco2)		!Doug 01/09 MFL1: record annual change
      REAL fuel_10hr_del(1:npft,1:nco2)
      REAL fuel_100hr_del(1:npft,1:nco2)
      REAL fuel_1000hr_del(1:npft,1:nco2)

      real num_fire(1:12)
      real annum_fire
      real area_burnt(1:12)
      real an_areafires
      real mfdi(1:12)
      real an_fdi
      real m_fc_crown(1:12),an_fc_crown
      real m_i_surface(1:12),an_i_surface
      real acflux_trace(1:6)           ! flux of trace gas species x depending on acflux_fire (g/m?²)
      real mcflux_trace(1:12,1:6)      ! flux of trace gas species x, monthly resolution (g/m?²)
c DM      real xmcflux_trace(1:12*6)       ! for passage to c driver
c DM      real nep,sum_anpp,sum_vegc,sum_litc,sum_soilc

c     variables for permafrost simulation
      real littercpool_ag
      real msnowpack(1:12)             ! monthly aver. snowpack depth in m
      real mw1(1:12)                   ! monthly values w(1), fraction avail. water
      real mw2(1:12)                   ! monthly values w(2), fraction avail. water
      real mw1_t(1:12)                 ! monthly values w(1), total water (liquid+ice)
      real mw2_t(1:12)                 ! monthly values w(2)
      real dthaw_depth(1:365)
      real maxthaw_depth               ! maximal permafrost thaw depth
      logical perm_attr
      real uw1,uw2,fw(1:2)
c DM      real delta_thawday
      real maxthaw_old

c     variables for the o18 analyses
c DM      real xmcica(1:12*npft)           ! monthly average ci/ca

      real mcica(1:12,1:npft)
      real meangc(1:12,1:npft)
      real mgp(1:12,1:npft)
c DM      real xmgpp(1:12*npft)            ! monthly gpp
      real gresp(1:12,1:npft,1:nco2)   ! monthly growth respiration
      real lresp(1:12,1:npft,1:nco2)   ! monthly leaf respiration
      real sresp(1:12,1:npft,1:nco2)   ! monthly sapwood respiration
      real rresp(1:12,1:npft,1:nco2)   ! monthly root respiration
      real aresp(1:12,1:npft,1:nco2)   ! monthly autotrophic respiration

c DM      real xmlresp(1:12*npft)          ! monthly leaf respiration
c DM      real xmsresp(1:12*npft)          ! monthly sapwood respiration
c DM      real xmrresp(1:12*npft)          ! monthly root respiration
c DM      real xmgresp(1:12*npft)          ! monthly growth respiration
c DM      real xmaresp(1:12*npft)          ! monthly autotrophic respiration
c DM      real xmeangc(1:12*npft)          ! monthly actual canopy conductance
c DM      real xmgp(1:12*npft)             ! monthly potential canopy conductance

#if defined(LPJ_STEP_1A) || defined(LPJ_STEP_1B) || defined(LPJ_STEP_2)
       integer spinup_years,rampup_years,run_years    !Doug 07/09: include rampup
#endif
       integer start_year                             !Doug 07/09: number of years in spinup and rampup
       
      real dhuman_ign(1:365)              !number of human-set fires per mln.ha
C Yan
       real num_fire_human(12),num_fire_lightn(12)
       real m_i_surface_human(12),m_i_surface_lightn(12)
       real area_burnt_human(12),area_burnt_lightn(12)
       real annum_fire_human,annum_fire_lightn
       real an_i_surface_human,an_i_surface_lightn
       real afire_frac_human,afire_frac_lightn
       real an_areafires_human,an_areafires_lightn
	
		real char_net_fuel
          real char_net_fuel_0

          real livegrass_0,dead_fuel_0,dead_fuel_all_0
         real fuel_all_0,fuel_1hr_total_0,fuel_10hr_total_0
         real fuel_100hr_total_0,fuel_1000hr_total_0	   


!           real deltaa
           real deltaa(1:npft),deltaa_fpc,delt_c13_fpc

           real mfire_frac(1:12)
           real  fbdep, ni_acc
      REAL dlm_1hr_old,dlm_10hr_old,dlm_100hr_old,dlm_1000hr_old
      REAL mlm(1:12)    !Doug 08/13: monthly litter moisture
c    Doug & Stephen 05/09: landuse change (for things like deforestation and furniture)
      LOGICAL landuse_on            !Doug 06/09: says if landuse is on or off. I'm taking this out soon to make it a little more elagent
          PARAMETER(landuse_on=.FALSE.)
      LOGICAL clear,abandon,agri    !Doug 06/09: clear,if the cells been cleared for agrculture; abandon,
                                    !  if the cells been gone from agrculture to natural (or noval...);
                                    !  agri, if the cell currently has over 50% crops or pasture
      REAL agri_burn,convflux
      REAL agri_litter_ag,agri_litter_bg
      REAL prod10(0:10),prod100(0:100),prod1
      REAL cflux_prod_total,prod10_total,prod100_total
      REAL anpp_agri,rap,x10(0:10),x100(0:100)
          PARAMETER(rap=1)
      REAl agri_t

c -----------------------------------------------------------------------------
c                                MAIN PROGRAM
c -----------------------------------------------------------------------------


c     Initialise input/output

c DM : to avoid converting from Fortran arrays to C/C++ arrays, we assume
c DM   a certain layout, but check its correctness at the beginning

      real dummy_array(1:2,1:3,1:4)
      real dummy1
      integer j,k
      dummy1=1.0

      do i=1,2
        do j=1,3
          do k=1,4
            dummy_array(i,j,k)=dummy1
            dummy1=dummy1+1
          enddo
        enddo
      enddo

      call initio(dummy_array)

c#if defined(LPJ_STEP_1A) || defined(LPJ_STEP_1B) || defined(LPJ_STEP_2)
c      call init_saved_data()
c      call get_spinup_years(spinup_years)
c      CALL get_rampup_years(rampup_years)    !Doug 07/09
c      call get_run_years(run_years)
c#endif


#ifdef LPJ_STEP_1A
      call init_saved_dataa()
      call get_spinup_years(spinup_years)
      CALL get_rampup_years(rampup_years)    !Doug 07/09
      call get_run_years(run_years)
#endif


#ifdef LPJ_STEP_1B
      call init_saved_dataa()
      call init_saved_datab()
      call get_spinup_years(spinup_years)
      CALL get_rampup_years(rampup_years)    !Doug 07/09
      call get_run_years(run_years)
#endif

#ifdef LPJ_STEP_2
      call init_saved_datab()
      call get_spinup_years(spinup_years)
      CALL get_rampup_years(rampup_years)    !Doug 07/09
      call get_run_years(run_years)
#endif

c    Doug 07/09: set starting year. For spinup (step 1a) or full run, starting year is zero. For
c      rampup (step 1b), spinup years (i.e, starts after spinup). For main run (step 2)
c      starts at spinup+rampup years.

c    if step1a or full run
      start_year=0;

#ifdef LPJ_STEP_1B
      start_year=spinup_years
#endif

#ifdef LPJ_STEP_2
      start_year=spinup_years+rampup_years
#endif

c     Obtain latitude and soil type for next gridcell from
c     input/output module
      call getgrid(lat,lon,soilcode,mlightn,a_nd,dogridcell)
      do while (dogridcell)


c       -------------------------------
c       Start of LOOP THROUGH GRIDCELLS
c       -------------------------------

c       Obtain soil parameters

        call soilparameters(soilcode,soilpar)

c       Obtain plant functional type parameters, define sapling and initial
c       grass mass structure and calculate maximum crown area
        call pftparameters(pftpar,sla,tree,evergreen,summergreen,
     *    raingreen,needle,boreal,lm_sapl,sm_sapl,hm_sapl,rm_sapl,
     *    latosa,allom1,allom2,allom3,allom4,wooddens,reinickerp
     *    ,co2,BTparam1,BTparam2,BTmode0)

c       Initialise year counter
#if defined(LPJ_STEP_2) || defined (LPJ_STEP_1B)

#ifdef LPJ_STEP_2 
c       Doug 07/09: Set to spinup_years  + rampup_years for reading back
        year=spinup_years+rampup_years
#else
c DM    Set to spinup_years for reading back
        year=spinup_years
#endif
c DM   mtemp has been saved, but is put back in mtemp_old, because
c DM   year 1000 is now the previous year
c DM   TODO : gdd has been added, but I am not sure it is needed
        call get_saved_data(year,k_fast_ave,k_slow_ave,
     *    litter_decom_ave,present,litter_ag_leaf,
     *    litter_ag_wood,
     *    fuel_1hr_leaf,fuel_1hr_wood,fuel_10hr,fuel_100hr,fuel_1000hr,
     *    litter_bg,crownarea,w,w_t,dwscal365,
     *    lm_ind,sm_ind,hm_ind,rm_ind,fpc_grid,mnpp,anpp,
     *    leafondays,leafoffdays,leafon,snowpack,mtemp_old,
     *    maxthaw_old,mw1,mw2,mw1_t,mw2_t,uw1,uw2,fw,
     *    mcica,mgpp,lresp,sresp,rresp,gresp,aresp,dbh,
     *    tau_c,cl_t,height_class,agpp,
     *    cpool_fast,cpool_slow,bm_inc,nind,gdd,lai_ind,
     *    height,w_ep,gdd_buf,mtemp_min_buf,mtemp_max_buf,
     *    lm_sapl,sm_sapl,hm_sapl,rm_sapl,meangc)

c DM    Year spinup_years as been done
        year=year+1
        call getclimate(year,mtemp,mtemp_dmin,mtemp_dmax,mprec,mwet,
     *    msun,mwindsp,
     *    popden,
     *    crop,pas,	!Doug 04/09
     *    co2,doyear)   !DIETER,kirsten
        !mtemp=mtemp+273.15
        !mtemp_dmin=mtemp_dmin+273.15
        !mtemp_dmax=mtemp_dmax+273.15

#else
        year=1

c       Obtain annual climate input parameters and atmospheric CO2
c       concentration from input/output module

c DM    TODO : should be done in initgrid
c       Initialise cpool:
        cpool_slow(1) = 0.
        cpool_slow(2) = 0.
        cpool_slow(3) = 0.

        cpool_fast(1) = 0.
        cpool_fast(2) = 0.
        cpool_fast(3) = 0.
c DM   was commented out
        w_ep       = 0.

        do pft=1,npft
          nind(pft)=0.0
c DM   Added : initialisation of lai_ind
          lai_ind(pft)=0.0
c DM   To avoid random values when tree(pft) is false
          height(pft)=0.0
          bm_inc(pft,1)=0.
          bm_inc(pft,2)=0.
          bm_inc(pft,3)=0.
          do m=1,12
            meangc(m,pft)=0.0
          enddo

        enddo

        call getclimate(year,mtemp,mtemp_dmin,mtemp_dmax,mprec,
     *    mwet,msun,mwindsp,
     *    popden,
     *    crop,pas,	!Doug 04/09
     *    co2,doyear)   !DIETER,kirsten

        !mtemp=mtemp+273.15
        !mtemp_dmin=mtemp_dmin+273.15
        !mtemp_dmax=mtemp_dmax+273.15


c       Gridcell initialisations

       ni_acc=0.0
       dlm_1hr_old=0.0
       dlm_10hr_old=0.0
       dlm_100hr_old=0.0
       dlm_1000hr_old=0.0
       aprec=0.0
       
        call initgrid(tree,k_fast_ave,k_slow_ave,litter_decom_ave,
     *    present,litter_ag,litter_ag_leaf,litter_ag_wood, !Doug 11/12: seperate out grass and wood litter
     *    fuel_1hr_leaf,fuel_1hr_wood,fuel_10hr,fuel_100hr,
     *    fuel_1000hr,litter_bg,crownarea,mprec,w,w_t,
     *    dwscal365,lm_ind,sm_ind,hm_ind,rm_ind,fpc_grid,mnpp,anpp,
     *    leafondays,leafoffdays,leafon,snowpack,mtemp_old,mtemp,
     *    maxthaw_old,mw1,mw2,mw1_t,mw2_t,uw1,uw2,fw,
     *    soilpar,mcica,mgpp,lresp,sresp,rresp,gresp,aresp,
     *    dbh,dbh_class, ! Doug 02/13: dbh_class added
     *    tau_c,cl_t,height_class,agpp)

#endif

c    Doug 05/09: land use, deciding if pixel has been recently abandoned, deforested ect
        IF (crop>0.5.AND.landuse_on) THEN    !Doug 06/08: Maybe try crop+pas?
          agri=.TRUE.
        ELSE
          agri=.FALSE.
        END IF

c    Doug 06/09: initalise land use logicals
        clear=.FALSE.
        abandon=.FALSE.

        IF (landuse_on) THEN   !Doug 06/09 if land use is off, then agriculture for land use is regared as zero
          agri_t=crop    !Doug 06/08: Maybe try crop+pas?
        ELSE
          agri_t=0
        END IF



c         Interpolate monthly climate data to daily values
c         for variables that do not depend on year
c	  Doug 01/09: moved down so daily lightening strikes s dependant on wet days
c
c          call daily1(mlightn,dlightn_control)
          fuel_1hr_del(:,:)=0.0
          fuel_10hr_del(:,:)=0.0
          fuel_100hr_del(:,:)=0.0
          fuel_1000hr_del(:,:)=0.0
		  
          fuel_1hr_del=fuel_1hr_leaf+fuel_1hr_wood		!Doug MFL1
          fuel_10hr_del=fuel_10hr

          fuel_100hr_del=fuel_100hr

          fuel_1000hr_del=fuel_1000hr


          afire_frac_afap_old=0.0
		  
c         Doug 03/09: MFL; set fuel incraments to zero
          fuel_1hr_leaf_inc_pos(:,:,:)=0.0
          fuel_1hr_wood_inc_pos(:,:,:)=0.0
          fuel_1hr_leaf_inc_neg(:,:,:)=0.0
          fuel_1hr_wood_inc_neg(:,:,:)=0.0
		  
        do while (doyear)
c         --------------------------------------
c         Start of LOOP THROUGH SIMULATION YEARS
c         --------------------------------------
           
          fuel_10hr_inc(:,:,:)=0.0
          fuel_100hr_inc(:,:,:)=0.0
          fuel_1000hr_inc(:,:,:)=0.0
		  
c         Interpolate monthly climate data to daily values
          call daily(mtemp,dtemp)
          call daily(mtemp_dmin,dtemp_min)
          call daily(mtemp_dmax,dtemp_max)
          call prdaily(mprec,dprec,mwet,year)
          aprec=SUM(mprec)
c Doug 07/09: Calculate a GDD for each grid cell. Used for ouput only.
          gdd_grid=0
          DO day=1,365
            IF (dtemp(day)>gddbase) THEN
               gdd_grid=gdd_grid+dtemp(day)-gddbase
            END IF
          END DO
           
          call daily_lightning(lat,lon,mlightn,dprec,dlightn,
     *      cgf,fdry,lt_days)

						!Doug 01/09: functions
							!distributed lighting
							!differently on days with & 
							!without precipitation
						!Doug 09/12: Name of function changed to
							!make more sense(itself function
							!changed as well, see Doug 09/12
							! comments below)

          call daily(msun,dsun)
          call daily(mwindsp,dwindsp)

c         Calculate mid-month daily photosynthetically active radiation flux,
c         daylength and mid-month daily potential evapotranspiration

c#ifdef LPJ_STEP_2
c
c          call petpar(year,dtemp,dsun,lat,mdayl,mpar_day,mpet,
c     *      u,v,sinelat,cosinelat,sinehh,hh,qo,spinup_years)
c#else
c
c          call petpar(year,dtemp,dsun,lat,mdayl,mpar_day,mpet,
c     *      u,v,sinelat,cosinelat,sinehh,hh,qo)
c#endif

          call petpar(year,dtemp,dsun,lat,mdayl,mpar_day,mpet,
     *      deet,    !Doug 08/09
     *      u,v,sinelat,cosinelat,sinehh,hh,qo,start_year)

          call daily(mpet,dpet)
          call daily(mdayl,dayl)
           
c         Adjust daily precipitation by snowmelt and accumulation in snowpack
          call snow(dtemp,dprec,snowpack,dmelt,msnowpack)

c         Calculate permafrost
          call permafrost(soilcode,mtemp,mw1,mw2,mw1_t,mw2_t,
     *      msnowpack,dtemp,dthaw_depth,maxthaw_depth,perm_attr,
     *      year,maxthaw_old)


c         Calculate summergreen phenology
          call summerphenology(pftpar,mtemp,dtemp,gdd,dphen_t,
     *      summergreen,tree)

c         Calculation of 20-year average climate variables
          call climate20(mtemp,dtemp,gdd,mtemp_min_buf,mtemp_max_buf,
     *      gdd_buf,year,mtemp_min20,mtemp_max20,gdd20,mtemp_max,
     *      pftpar)

c         Implement PFT bioclimatic limits
          call bioclim(pftpar,mtemp_min20,mtemp_max20,gdd,mtemp_max,
     *      survive,estab,year,present)
            

c        Doug 05/09: agriculture routine, updates pools due to abandonment
c        or clearance of grid cell 
	  
C          call agriculture(clear,mtemp_min20,tree,agri_burn, 
C     *      lm_ind,sm_ind,hm_ind,rm_ind,nind,convflux,agri_litter_ag,
C     *      agri_litter_bg,
C     *      litter_ag_leaf,litter_ag_wood, ! Doug 11/12: seperate out wood and grass litter
C     *      litter_bg,abandon,prod10,prod100,
C     *      fpc_grid,present)


c         Calculation of GPP and soil water balance
          call gpp(present,co2,soilpar,pftpar,lai_ind,fpc_grid,
     *      mdayl,mtemp,mpar_day,dphen_t,w,w_t,dpet,dprec,dmelt,sla,
     *      agpp,alresp,arunoff_surf,arunoff_drain,arunoff,mrunoff,
     *      dwscal365,dphen_w,dphen,
     *      dphen_change,wscal,mgpp,mlresp,mw1,dw1,aaet, !Doug 05/09: inc. dphen_change
     *      leafondays,leafoffdays,leafon,tree,raingreen,mpar,
     *      mapar,mphen,year,maet,littercpool_ag,mw2,maxthaw_depth,
     *      dthaw_depth,uw1,uw2,fw,
     *      litter_ag_leaf,litter_ag_wood,	!Seperate out grass and wood litter
     *      mcica,mauw1,mauw2,
     *      lat,lon,w_ep,d_evap,maep,aaep,mintc,aintc,needle,boreal,
     *      dayl,mpet_grid,apet_grid,mpet2,apet,mpet,mmelt,
     *      mw1_t,mw2_t,meangc,mgp,deltaa,deltaa_fpc,
     *      deet,alpha_ws)    !Doug 07-08/09

          
c         Calculation of mid-month soil temperatures
          call soiltemp(soilpar,mtemp,mtemp_old,mtemp_soil,mw1,mw1_t)

c         Interpolate monthly soil temperature to daily values
          call daily(mtemp_soil,dtemp_soil)
          

          call npp(pftpar,dtemp,dtemp_soil,tree,dphen,nind,
     *      lm_ind,sm_ind,rm_ind,mgpp,anpp,mnpp,bm_inc,present,
     * lresp,sresp,rresp,gresp,aresp,year,agpp,delt_c13_fpc,fpc_grid)

           
c         Allocation to reproduction
          call reproduction(bm_inc,lm_sapl,sm_sapl,hm_sapl,rm_sapl,
     *        litter_ag_leaf,    !Doug 11/12: seperate out grass and wood litter
     *        fuel_1hr_leaf,
     *        fuel_1hr_leaf_inc_pos,fuel_1hr_leaf_inc_neg,    !Doug 01/09: fuel_1hr_inc_i
     *        present,tree,co2)
	 
c         Calculation of leaf, sapwood, and fine-root turnover


          call turnover(pftpar,present,tree,lm_ind,sm_ind,hm_ind,
     *      rm_ind,
     *      litter_ag_leaf, !Doug 11/12: seperate out grass and wood litter
     *      litter_bg,fuel_1hr_leaf,
     *      fuel_10hr,fuel_100hr,
     *      fuel_1hr_leaf_inc_pos,fuel_1hr_leaf_inc_neg,
     *      fuel_10hr_inc, fuel_100hr_inc,	!Doug 01/09: fuel_xhr_inc
     *      nind,turnover_ind)

            
c         Litter and soil decomposition calculations
c         This is done before fire, so that fire probability is calculated
c         on litter remaining after year's decomposition
          !PRINT*, "******************"
          !PRINT*, fuel_1hr_leaf(:,1)

          call littersom(pftpar,litter_ag_leaf,litter_ag_wood, !Doug 11/12: seperate out grass and wood litter
     *      litter_bg,fuel_1hr_leaf,fuel_1hr_wood,fuel_10hr,
     *      fuel_100hr,fuel_1000hr,
     *      fuel_1hr_leaf_inc_pos, fuel_1hr_leaf_inc_neg,
     *      fuel_1hr_wood_inc_pos, fuel_1hr_wood_inc_neg,
     *      fuel_10hr_inc,fuel_100hr_inc,fuel_1000hr_inc,       !Doug 01/09: fuel_xhr_inc
     *      mw1,mw1_t, mtemp_soil,cpool_fast,cpool_slow,
     *      arh,mrh,year,k_fast_ave,k_slow_ave,litter_decom_ave,
     *      agri_litter_ag,agri_litter_bg,agri)		!Doug 05/09: agriculture varibles added for land use change stuff

c         Removal of PFTs with negative C increment this year
          call kill_pft(bm_inc,present,tree,lm_ind,rm_ind,hm_ind,
     *      sm_ind,nind,litter_ag_leaf,litter_ag_wood, !Doug 11/12: seperate out grass and wood litter
     *      litter_bg,year,fuel_1hr_leaf,fuel_1hr_wood,fuel_10hr,
     *      fuel_100hr,fuel_1000hr,
     *      fuel_1hr_leaf_inc_pos,fuel_1hr_leaf_inc_neg,
     *      fuel_1hr_wood_inc_pos,fuel_1hr_wood_inc_neg,
     *      fuel_10hr_inc, fuel_100hr_inc,fuel_1000hr_inc)              !Doug 01/09: fuel_xhr_inc

     
            
c         Allocation of annual carbon increment to leaf, stem and fine root
c         compartments


C            IF (lat>-25) THEN 
C              PRINT*,"kill_pft"
C              PRINT*, fuel_1hr_leaf(:,1)
C            END IF

          call allocation(pftpar,allom1,allom2,allom3,latosa,
     *      wooddens,reinickerp,tree,sla,wscal,nind,
     *      bm_inc,lm_ind,lm_inc,sm_ind,hm_ind,rm_ind,
     *      crownarea,fpc_grid,lai_ind,height,
     *      height_class,dbh,dbh_class,tau_c,cl_t,!Doug 02/13: dbh_class added
     *      litter_ag_leaf,litter_ag_wood,        !Doug 11/12: seperate out grass and wood litter
     *      litter_bg,fuel_1hr_leaf,fuel_1hr_wood,
     *      fuel_10hr,fuel_100hr,fuel_1000hr,
     *      fuel_1hr_leaf_inc_pos, fuel_1hr_leaf_inc_neg,
     *      fuel_1hr_wood_inc_pos, fuel_1hr_wood_inc_neg,
     *      fuel_10hr_inc,fuel_100hr_inc,fuel_1000hr_inc,!Doug 01/09: fuel_xhr_inc
     *      fpc_inc,present,year, evergreen)

   
c         Implement light competition between trees and grasses
          call light(present,tree,lm_ind,sm_ind,hm_ind,rm_ind,
     *      crownarea,fpc_grid,fpc_inc,nind,
     *      litter_ag_leaf,litter_ag_wood, !Doug 11/12: seperate out grass and wood litter
     *      litter_bg,
     *      fuel_1hr_leaf,fuel_1hr_wood,
     *      fuel_10hr,fuel_100hr,fuel_1000hr,
     *      fuel_1hr_leaf_inc_pos, fuel_1hr_leaf_inc_neg, 
     *      fuel_1hr_wood_inc_pos, fuel_1hr_wood_inc_neg, 
     *      fuel_10hr_inc,fuel_100hr_inc,fuel_1000hr_inc,!Doug 01/09: fuel_xhr_inc
     *      sla,year)

            
c         Implement light competition and background mortality among tree PFTs
c         (including heat damage and due to lower limit of npp for boreal trees)
          call mortality(pftpar,present,tree,boreal,bm_inc,
     *      turnover_ind,sla,lm_ind,sm_ind,hm_ind,rm_ind,nind,
     *      litter_ag_leaf,litter_ag_wood, !Doug 11/12: seperate out grass and wood litter
     *      litter_bg,fuel_1hr_leaf,fuel_1hr_wood,fuel_10hr,fuel_100hr,
     *      fuel_1000hr,
     *      fuel_1hr_leaf_inc_pos, fuel_1hr_leaf_inc_neg,
     *      fuel_1hr_wood_inc_pos, fuel_1hr_wood_inc_neg,
     *      fuel_10hr_inc,fuel_100hr_inc,fuel_1000hr_inc,!Doug 01/09: fuel_xhr_inc
     *      dtemp,anpp,mtemp_max,year)

            
c         Calculation of biomass destruction by fire disturbance
        fuel_1hr_del=fuel_1hr_del-fuel_1hr_leaf-fuel_1hr_wood		!Doug MFL1
        fuel_10hr_del=fuel_10hr_del-fuel_10hr
        fuel_100hr_del=fuel_100hr_del-fuel_100hr
        fuel_1000hr_del=fuel_1000hr_del-fuel_1000hr

           
         call fire(year,start_year,                                  !Doug 07/09: added start year
     *      pftpar,dtemp,dtemp_min,dtemp_max,dprec,
     *      dwindsp,dlightn,dphen,dphen_change,                      !Doug 06/09: dphen_change added for fire paradox experiments
     *      litter_ag_leaf,litter_ag_wood,							 !Doug 11/12: Doug 11/12: seperate out grass and wood litter
     *      litter_bg,fuel_1hr_leaf,fuel_1hr_wood,
     *      fuel_10hr,fuel_100hr,fuel_1000hr,
     *      pfuel_limit,
     *      fuel_1hr_leaf_inc_pos, fuel_1hr_leaf_inc_neg,
     *      fuel_1hr_wood_inc_pos, fuel_1hr_wood_inc_neg,
     *      fuel_10hr_inc,fuel_100hr_inc, fuel_1000hr_inc,	                     !Doug 01/09 MFL: fuel_xhr_inc
     *      fuel_1hr_del, fuel_10hr_del, fuel_100hr_del, 
     *      fuel_1000hr_del,                                         !Doug 01/09 MFL1
     *      mfuel_1hr_total,
     *      mfuel_1hr_leaf_total,mfuel_1hr_wood_total,
     *      mfuel_10hr_total,                       !Doug 03/09: outputs for monthly fuel loads
     *      mfuel_100hr_total, mfuel_1000hr_total, 
     *      mlivegrass,
     *      acflux_fire,mcflux_fire,
     *      afire_frac,lm_ind,rm_ind,sm_ind,hm_ind,nind,dw1,present,
     *      tree,lat,mw1,fpc_grid,popden,a_nd,height,height_class,
     *      dbh,dbh_class,tau_c,cl_t,  !Doug 02/13: dbh_class added
     *      BTparam1,BTparam2,BTmode0, !Doug 02/13
     *      num_fire,annum_fire,area_burnt,
     *      an_areafires,mfdi,an_fdi,an_fseason,mcflux_trace,
     *      acflux_trace,m_fc_crown,an_fc_crown,m_i_surface,
     *      an_i_surface,
     *      dhuman_ign,
     *      num_fire_human,num_fire_lightn,
     *      annum_fire_human,annum_fire_lightn,
     *      area_burnt_human,area_burnt_lightn,
     *      an_areafires_human,an_areafires_lightn,
     *      afire_frac_human,afire_frac_lightn,char_net_fuel_0,
     * livegrass_0,dead_fuel_0,dead_fuel_all_0,
     * fuel_all_0,fuel_1hr_total_0,fuel_10hr_total_0,
     * fuel_100hr_total_0,fuel_1000hr_total_0,
     * mfire_frac,lon,crop,pas,fbdep,ni_acc,afire_frac_afap_old,
     * dlm_1hr_old,dlm_10hr_old,dlm_100hr_old,dlm_1000hr_old,mlm,
     * dpet,aprec)

          fuel_1hr_del=fuel_1hr_leaf+fuel_1hr_wood		!Doug MFL1
          fuel_10hr_del=fuel_10hr

          fuel_100hr_del=fuel_100hr

          fuel_1000hr_del=fuel_1000hr

c         Doug 03/09: MFL; set fuel incraments to zero
          fuel_1hr_leaf_inc_pos=0.0
          fuel_1hr_wood_inc_pos=0.0
          fuel_1hr_leaf_inc_neg=0.0
          fuel_1hr_wood_inc_neg=0.0
          fuel_10hr_inc=0.0
          fuel_100hr_inc=0.0
          fuel_1000hr_inc=0.0
 
c         Establishment of new individuals (saplings) of woody PFTs,
c         grass establishment, removal of PFTs not adapted to current climate,
c         update of individual structure and FPC.
          call establishment(pftpar,present,survive,estab,nind,	
     *      lm_ind,sm_ind,rm_ind,hm_ind,lm_sapl,sm_sapl,rm_sapl,
     *      hm_sapl,crownarea,fpc_grid,lai_ind,height,dbh,dbh_class, ! Doug 03/13: dbh_class added
     *      tau_c,cl_t,BTparam1,BTparam2,BTmode0,		!Doug 02/13: described BT spread
     *      sla,wooddens,latosa,mprec,reinickerp,
     *      litter_ag_leaf,litter_ag_wood, ! Doug 11/12: seperate out grass and wood litter
     *      litter_bg,fuel_1hr_leaf,fuel_1hr_wood,
     *      fuel_10hr,fuel_100hr,fuel_1000hr,
     *      fuel_1hr_leaf_inc_pos,fuel_1hr_leaf_inc_neg,               !Doug 03/09
     *      fuel_1hr_wood_inc_pos,fuel_1hr_wood_inc_neg,               !Doug 03/09
     *      fuel_10hr_inc,fuel_100hr_inc,fuel_1000hr_inc,    !Doug 03/09
     *      pas,crop,                                             !Doug 06/09
     *      tree,allom1,allom2,allom3,acflux_estab,leafondays,
     *      leafoffdays,leafon,mnpp,anpp,mnpp_add,anpp_add,year)

            
     
c       Doug 05/09: agricultural production and product decomposition     
C          call agriprod(present,agri,anpp,cflux_prod_total,
C     *   prod10_total,prod100_total,prod1,anpp_agri,rap,x10,
C     *   x100,agri_litter_ag,agri_litter_bg,prod10,prod100)




c         radioactive decay of 14C
          
          call decay(present,lm_ind,sm_ind,hm_ind,rm_ind,
     *      litter_ag_leaf,litter_ag_wood, !Doug 11/12: seperate out grass and wood litter
     *      litter_bg,cpool_fast,cpool_slow,
     *      fuel_1hr_leaf,fuel_1hr_wood,
     *      fuel_10hr,fuel_100hr,fuel_1000hr)


c KIRSTEN 25 July 2006: set pft array to zero, where pft not present before passing on to cpp driver.
c Not very familiar with cpp, normally it should be done in the driver.

          do pft=1,npft
            if (.not.present(pft)) then
              nind(pft)=0.0
              do nc=1,nco2
                lm_ind(pft,nc)=0.0
                rm_ind(pft,nc)=0.0
                sm_ind(pft,nc)=0.0
                hm_ind(pft,nc)=0.0
                litter_ag(pft,nc)=0.0
                litter_ag_leaf(pft,nc)=0.0
                litter_ag_wood(pft,nc)=0.0
                litter_bg(pft,nc)=0.0
                anpp(pft,nc)=0.0
              enddo
              fpc_grid(pft)=0.0
c DM  sla is only set in pftparameters : if present(pft) becomes true
c DM  again, its value is lost !!
c DM              sla(pft)=0.0
              gdd(pft)=0.0
              height(pft)=0.0
              do m=1,12
                do nc=1,nco2
                  mgpp(m,pft,nc)=0.0
                enddo
              enddo
              do d=1,365
                dphen(d,pft)=0.0
              enddo
            endif
          enddo

c  Doug 06/09: Find monthly mnpp for entire cell for ouput only.
c  Doug 07/09: For future output, with no deltaC needed, cheat
         DO m=1,12
           mnpp_grid(m)=SUM(mnpp(m,:,1))
         END DO

c  anpp(pft)
         anpp_grid(:)=anpp(:,1)
c  arf(1)
         arh_grid=arh(1)
c  acflux_fire(1)

C	Doug 11/12: Model changed to seperate out grass/leaf and wood above ground litter.
C				So litter_ag is summed again for output
         acflux_fire_grid=acflux_fire(1)
          DO pft=1,npft
          litter_ag(pft,1)=litter_ag_leaf(pft,1)+litter_ag_wood(pft,1)
            DO nc=2,nco2
              litter_ag(pft,nc)=((litter_ag_leaf(pft,nc)*
     *          litter_ag_leaf(pft,1))+(litter_ag_wood(pft,nc)*
     *          litter_ag_wood(pft,1)))/(litter_ag_wood(pft,1)+
     *          litter_ag_leaf(pft,1))
	        END DO
          END DO
          dprec_out=dprec

          call outannual(year,present,nind,lm_ind,lm_inc,rm_ind,sm_ind,
     *      hm_ind,fpc_grid,anpp,acflux_estab,
     *      litter_ag,litter_ag_leaf,litter_ag_wood, !Doug 11/12: Output seperate ggrass and wood litter
     *      litter_bg,
     *      cpool_fast,cpool_slow,arh,afire_frac,acflux_fire,
     *      mcflux_fire,arunoff,sla,mpar,mapar,mphen,anpp_add,mnpp,
     *      mnpp_grid,    !Doug 06/09: monthly npp for entire cell
     *      mrunoff,aaet,mrh,mnpp_add,maet,mtemp_soil, 
     *      mcica,lresp,
     *      sresp,rresp,gresp,aresp,mauw1,mauw2,maep,aaep,mpet_grid,
     *      apet_grid,mintc,aintc,meangc,mgp,num_fire,annum_fire,
     *      area_burnt,an_areafires,mfdi,an_fdi,an_fseason,
     *      acflux_trace,mcflux_trace,m_fc_crown,an_fc_crown,
     *      m_i_surface,an_i_surface,gdd,height,mgpp,
     *      w(1),w(2),dphen,
     *      dphen_change, !Doug 05/09: inc. dphen_change
     *      dhuman_ign,dlightn,mw1,mw2,num_fire_human,num_fire_lightn,
     *      annum_fire_human,annum_fire_lightn,area_burnt_human,
     *      area_burnt_lightn,an_areafires_human,an_areafires_lightn,
     *      afire_frac_human,afire_frac_lightn,char_net_fuel_0,
     *      livegrass_0,dead_fuel_0,dead_fuel_all_0,
     *      fuel_all_0,fuel_1hr_total_0,fuel_10hr_total_0,
     *      fuel_100hr_total_0,fuel_1000hr_total_0,
     *      mfuel_1hr_total,
     *      mfuel_1hr_leaf_total,mfuel_1hr_wood_total,
     *      mfuel_10hr_total, mfuel_100hr_total, !Doug 04/09
     *      mfuel_1000hr_total,mlivegrass,                        !Doug 04/09
     *      deltaa,deltaa_fpc,delt_c13_fpc,
     *      mfire_frac,
     *      fbdep,litter_decom_ave,turnover_ind,
     *      crop,pas,	!Doug 05/09: just checking crops and pasture are implimented properly
     *      anpp_grid,arh_grid,acflux_fire_grid,                  !Doug 07/09: cheating future run ouputs without deltaC's
     *      gdd_grid,alpha_ws,pfuel_limit,dprec_out,!Doug 07/09: biocliamtic variables for heat and water stress
     *      BTparam1,BTparam2,cgf,fdry,lt_days,dlm_1hr_old,mlm)                                    

cccc note: you can alternatively use mpet2 and apet to get PET*1.32
cccc       instead of mpet_grid and apet_grid, respectively

           
#ifdef LPJ_STEP_1A
          if (year.eq.spinup_years) then
c DM      Saved variables are the ones initialized in initgrid plus
c DM      thoses that have been found necessary to get right results
            call put_saved_data(year,k_fast_ave,k_slow_ave,
     *        litter_decom_ave,present,
     *        litter_ag_leaf,litter_ag_wood,   !Doug 11/12: seperate out grass and wood litter
     *        fuel_1hr_leaf,fuel_1hr_wood,
     *        fuel_10hr,fuel_100hr,fuel_1000hr,
     *        litter_bg,crownarea,w,w_t,dwscal365,
     *        lm_ind,sm_ind,hm_ind,rm_ind,fpc_grid,mnpp,anpp,
     *        leafondays,leafoffdays,leafon,snowpack,mtemp,
     *        maxthaw_old,mw1,mw2,mw1_t,mw2_t,uw1,uw2,fw,
     *        mcica,mgpp,lresp,sresp,rresp,gresp,aresp,dbh,
     *        tau_c,cl_t,height_class,agpp,
     *        cpool_fast,cpool_slow,bm_inc,nind,gdd,lai_ind,
     *        height,w_ep,gdd_buf,mtemp_min_buf,mtemp_max_buf,
     *        lm_sapl,sm_sapl,hm_sapl,rm_sapl,meangc)
          endif
#endif

#ifdef LPJ_STEP_1B
          IF (year.eq.spinup_years+rampup_years) THEN
c DM      Saved variables are the ones initialized in initgrid plus
c DM      thoses that have been found necessary to get right results
            CALL put_saved_data(year,k_fast_ave,k_slow_ave,
     *        litter_decom_ave,present,
     *        litter_ag_leaf,litter_ag_wood, !Doug 11/12: Fill:
     *        fuel_1hr_leaf,fuel_1hr_wood,
     *        fuel_10hr,fuel_100hr,fuel_1000hr,
     *        litter_bg,crownarea,w,w_t,dwscal365,
     *        lm_ind,sm_ind,hm_ind,rm_ind,fpc_grid,mnpp,anpp,
     *        leafondays,leafoffdays,leafon,snowpack,mtemp,
     *        maxthaw_old,mw1,mw2,mw1_t,mw2_t,uw1,uw2,fw,
     *        mcica,mgpp,lresp,sresp,rresp,gresp,aresp,dbh,
     *        tau_c,cl_t,height_class,agpp,
     *        cpool_fast,cpool_slow,bm_inc,nind,gdd,lai_ind,
     *        height,w_ep,gdd_buf,mtemp_min_buf,mtemp_max_buf,
     *        lm_sapl,sm_sapl,hm_sapl,rm_sapl,meangc)
          endif
#endif

c         Increment year counter
          year=year+1

c         Obtain annual climate input parameters and atmospheric CO2
c         concentration from input/output module
          call getclimate(year,mtemp,mtemp_dmin,mtemp_dmax,mprec,mwet,
     *      msun,mwindsp,
     *      popden,
     *      crop, pas,    !Doug 05/09
     *      co2,doyear)   !DIETER,kirsten

        !mtemp=mtemp+273.15
        !mtemp_dmin=mtemp_dmin+273.15
        !mtemp_dmax=mtemp_dmax+273.15

c         Doug 05/09: land use, deciding if pixel has been recently abandoned, deforested ect
        IF (crop>0.5.AND.landuse_on) THEN  !Doug 06/09: if cell is greater then 50% agrculture
                            !Doug 06/08: Maybe try crop+pas?
          agri=.TRUE.
        ELSE
          agri=.FALSE.
        END IF

        clear=.FALSE.
        abandon=.FALSE.
        IF (crop>0.5.AND.agri_t<0.5.AND.landuse_on) THEN      !Doug 06/09: Checks to see if cell has gone from none
                                                              !Doug 06/08: Maybe try crop+pas?
                                                                !agriculture (agri_t, last timestep) to agriculture
                                                                !(crop+pas, this timestep)
          clear=.TRUE.
        ELSE IF(crop<0.5.AND.agri_t>0.5.AND.landuse_on) THEN  !Doug 06/09: Checks if cell has gone from agriculture
                                                                !none agriculture

          abandon=.TRUE.
        END IF
        agri_t=crop    !Doug 06/08: Maybe try crop+pas?

#ifdef LPJ_STEP_1A
c DM    If simulation in 2 steps, don't believe getclimate about when it's done
          if (year.eq.(spinup_years+1)) doyear=.false.
#endif

#ifdef LPJ_STEP_1B
c DM    If simulation in 2 steps, don't believe getclimate about when it's done
          IF (year==(spinup_years+rampup_years+1)) doyear=.false.
#endif

c          pause

c    Kirsten: Why are leafondays and leafoffdays not set to zero annuallY? 
c             they accumulate over years, what does that mean?
c          do pft=1,npft
c             leafondays(pft)=0
c             leafoffdays(pft)=0
c          enddo

        enddo   !loop through simulation years

c       Final output to input/output module

        call outgrid()

c       Obtain latitude and soil type for next gridcell from
c       input/output module

!        call getgrid(lat,soilcode,mlightn,a_nd,dogridcell)
        call getgrid(lat,lon,soilcode,mlightn,a_nd,dogridcell)


      enddo     !loop through gridcells


c     Terminate input/output


!!        close(1050)
c#if defined(LPJ_STEP_1A) || defined(LPJ_STEP_1B) || defined(LPJ_STEP_2)
c      call end_saved_data()
c#endif

#ifdef LPJ_STEP_1A
      call end_saved_dataa()
#endif

#ifdef LPJ_STEP_1B
      call end_saved_dataa()
      call end_saved_datab()
#endif

#ifdef LPJ_STEP_2
      call end_saved_datab()
#endif

      call termio
      end



c -----------------------------------------------------------------------------
c                                SUBROUTINES
c -----------------------------------------------------------------------------
c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE DAILY2
c     Doug 01/09: For calculating lightning strikes. Based on linear
c	interpolation, but with grouping around rain events
c
c	Completed in two steps. Step 1 calculates the daily spread of lighning if 
c	lighning was constant, with more ligtning on days with rainfall then days 
c	without. Step2 to multiplys by linearly interpolated lighting strkes, and 
c	scales so total lighting remains the same.  

	SUBROUTINE daily2(lat,lon,mval,dval1,dval2)

	IMPLICIT NONE

	INTEGER	nmonths			! number of months in a year
		PARAMETER (nmonths=12)
	INTEGER	ndayyear		! number of days in a year
		PARAMETER (ndayyear=365)
	INTEGER	ndaymonth(1:nmonths)	!number of days in each month
	INTEGER	m
        	DATA (ndaymonth(m),m=1,12)
     *    / 31,28,31,30,31,30,31,31,30,31,30,31 /


C	LOCAL VARAIBLES:
	INTEGER	month,mday, day	!month and day of month and day of year

c     I/O:
      	REAL mval(nmonths),dval1(ndayyear), dval2(ndayyear)

c	For step 1:
	REAL	fdal(365)		!the value of the function that 
					!redistributes lighting for each day
	REAL	nwd, pwd		!number and proption of wet days in the month
	REAL	gceo			!coeffiant used to calculate proportion of
					!lighting strikes on wet days
		PARAMETER (gceo=0.00001)

c	For step 2:
	REAL	lmval(2)		!for calculating scaling of rediustributed
					!lightning

	REAL	lat,lon




c       --------------------------------------------------------------------------
c       Step 1

        day=0

        DO month=1,nmonths	!Month of year
          nwd=0



          DO mday=1,ndaymonth(month)	!day of month
            day=day+1
            IF(dval1(day)>0) nwd=nwd+1
          END DO	!day of month

          pwd=nwd/ndaymonth(month) !pwd= fraction of wet days in month 
            day=day-ndaymonth(month)
            DO mday=1,ndaymonth(month)	!day of month
              day=day+1    !Seaperates days and calcualates lighting for days with and 
                           !without rain

              IF(dval1(day)>0) THEN
                fdal(day)=(pwd**gceo)/nwd
              ELSE
                fdal(day)=(1-pwd**gceo)/(ndaymonth(month)-nwd)
              END IF

            END DO
          END DO	!Month of year
	
c-------------------------------------------------------------------------------
c Step 2

	
	CALL	daily1(mval,dval2)		!dval2=lighting for a simple
                                                !linear interpolation between
                                                !montrhly lighting values


	day=0
        lmval(:)=0.0

	DO month=1,nmonths !month of year

		
          DO mday=1,ndaymonth(month)	!day of month
            day=day+1
            lmval(1)=lmval(1)+dval2(day)
            dval2(day)=dval2(day)*fdal(day)
            lmval(2)=lmval(2)+dval2(day)

          END DO 	!day of month
	
          DO mday=day+1-ndaymonth(month),day
            IF (lmval(1)==0.AND.lmval(2)==0) THEN
              dval2(mday)=0
            ELSE IF(lmval(2)==0.AND.lmval(1)/=0)	THEN
              PRINT*, "error:redistrbuting lighting around"
              PRINT*, "wet days in daily2 subroutine"
              STOP
            END IF

            dval2(mday)=dval2(mday)*lmval(1)/lmval(2)
              IF (dval1(mday)>0) dval2(mday)=0	!removes all lightin
                                                !on wet days

         END DO


	END DO !month of year




	return
	END !SUBROUTINE	!Daily2
c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE DAILY_LIGHTNING
c     Doug 01/09: For calculating lightning strikes. Based on linear
c	interpolation, but with grouping around rain events
c
c	Completed in two steps. Step 1 calculates the daily spread of lighning if 
c	lighning was constant, with more ligtning on days with rainfall then days 
c	without. Step2 to multiplys by linearly interpolated lighting strkes, and 
c	scales so total lighting remains the same.  
c
c		Doug 09/12: function changed to remove iter-cloud
c	lightning. See Doug 09/12 comments below for deatils.

      SUBROUTINE daily_lightning(lat,lon,mval,dval1,dval2
     *  ,cgf,fdry,lt_days)

      IMPLICIT NONE

C   Obvious parameters
      INTEGER	nmonths			! number of months in a year
        PARAMETER (nmonths=12)
      INTEGER	ndayyear		! number of days in a year
        PARAMETER (ndayyear=365)
      INTEGER	ndaymonth(1:nmonths)	!number of days in each month
      INTEGER	m
        	DATA (ndaymonth(m),m=1,12)
     *    / 31,28,31,30,31,30,31,31,30,31,30,31 /


C	LOCAL VARAIBLES:
      INTEGER	month,mday, day	!month and day of month and day of year
      REAL cgf(1:nmonths) ! Doug 09/12: fraction of lighting that is cloud-to-ground.
      REAL daily_stikes !Doug 04/13: number of strikes per km per day

c     I/O:
      REAL mval(nmonths),dval1(ndayyear)
      REAL dval2(ndayyear)

c	For step 1:
      REAL	fdal(365)		!the value of the function that 
					!redistributes lighting for each day
      REAL	nwd, pwd	!number and proption of wet days in the month, and number of dry days
      INTEGER  ndd, nld !Doug 11/12 number of dry days and number of lighting days
      INTEGER  r !conatains random interger
      REAL lt_days(1:nmonths) !Doug 11/12 proportion of dry days where lightning will strikes
      INTEGER ltk(1:31) !Doug 11/12 stores the day of the month when there is no rain
		
c	For step 2:
      REAL	lmval(2)		!for calculating scaling of rediustributed
					!lightning

      REAL	lat,lon
      REAL flt, fdry(1:nmonths),fwet
      
C   Doug 04/13: Parameters for calculating scaling factors
C       Doug 04/13: CG-IC ration
c       Uses fucntion CG=cgp1*LT^cgp2
c		found by comparing US CG and total lighting data from Vaisala/LIS 
      REAL  cgp1,cgp2
        PARAMETER (cgp1=0.0001267)
        PARAMETER (cgp2=-0.4180)
      
C       Doug 04/13: Wet day/Dry day       
      REAL  wdp1,wdp2
        PARAMETER (wdp1=0.85033)
        PARAMETER (wdp2=-2.835)
        
C       Doug 04/13: Wet day/Dry day       
      REAL  ldp1,ldp2
        PARAMETER (ldp1=1.099)
        PARAMETER (ldp2=94678.69)
       

c Doug 11/12 Declarations for the random generator. This is called to randomoize 
c lighting days amoungst the dry days of the month
      integer k10a(4)
       real random
       data k10a /73,24,881,52/


c       --------------------------------------------------------------------------
c       Step 1
        day=0
        fdal(:)=0
        lt_days(:)=0
		
        DO month=1,nmonths	!Month of year
            ltk(:)=0 !Doug 10/12
          
            daily_stikes=mval(month)/
     *          (1000000.0*ndaymonth(month)) !Doug 04/13: daily strikes per m2
     
     
c		Doug 09/12. Removed inter-cload lighting by calculating
c		the fraction of cload ground (CG) from total (LT)
c		lightning.			
            cgf(month)=cgp1*(daily_stikes)**(cgp2)	  

            IF(cgf(month)>1.0) cgf(month)=1.0
            IF(cgf(month)<0.0) cgf(month)=0.0          

                   

C Doug 03/13: claculates the number of wet (nwd) and dry (ndd) days          
            ndd=0.0
            nwd=0.0  
            DO mday=1,ndaymonth(month)	!day of month
                day=day+1
                IF(dval1(day)>0.0) THEN 
                    nwd=nwd+1
                ELSE
                    ndd=ndd+1
                    ltk(ndd)=day
                END IF
            END DO	!day of month
		  
        
            pwd=nwd/ndaymonth(month) !pwd= fraction of wet days in month 
            day=day-ndaymonth(month)
            flt=wdp1*exp(pwd*wdp2)

            IF (flt<0.0) flt=0.0
            IF (flt>1.0) flt=1.0


c       Doug 10/12: Determine number of lightning days
            lt_days(month)=1-1/(ldp1*
     *          (daily_stikes*cgf(month)*flt+1)**ldp2)

        
            IF (lt_days(month)<0.0) lt_days(month)=0.0
            IF (lt_days(month)>1.0) lt_days(month)=1.0

            nld=ANINT(lt_days(month)*ndd)
            IF (nld==0) nld=1

            IF (nwd==0) THEN  
                fdry(month)=1/(nld)
                fwet=0
            ELSE
                fdry(month)=(flt)/(nld)
                fwet=(1-flt)/(nwd)
            END IF		   	
			
        
C Doug 03/13:   1)Randomly clalculted a day that is a) not wet and b) doesn;t have any
c                   strikes		
c               2)Sets lightning strikes on that day
c               3)reduced the amount of available none wet, non lighting days by 1
c               4) doed it again
            DO mday=1,nld         
                r=CEILING(ndd*random(k10a))
                
                IF (r>ndd) r=ndd
                IF (r<1) r=1
                
                fdal(ltk(r))=fdry(month)                
                dval2(ltk(r))=mval(month)*fdry(month)

                ndd=ndd-1
                ltk(r:ndd)=ltk(r:ndd)+1
            END DO
			
            DO mday=1,ndaymonth(month)	!day of month
                day=day+1    !Seperates days and calcualates lighting for days with and 
                           !without rain
            
                IF(dval1(day)>0) THEN
                    fdal(day)=0
              !ELSE
               ! fdal(day)=fdry(month) !Doug 03/13: already calculated, so now comment out
                END IF
            END DO
            
        END DO	!Month of year
          
c      Doug 04/13: Step 2 is now gone! The interpolation it contained didn't make
c      much difference, and is just an added laye of confusion now we have lightning days
          
      return
      
      
      END !SUBROUTINE	!daily_lightning
c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE DAILY1
c     Linear interpolation of monthly to quasi-daily values

      subroutine daily1(mval,dval)

      implicit none

      integer nmonth
        parameter (nmonth=12)        !number of months per year
      integer ndayyear
        parameter (ndayyear=365)     !number of days per year
      integer ndaymonth(1:12)    !number of days in each month
      integer m
        data (ndaymonth(m),m=1,12)
     *    / 31,28,31,30,31,30,31,31,30,31,30,31 /

c     ARGUMENTS:
      real mval(nmonth),dval(ndayyear)

c     LOCAL VARIABLES:
      real dd,todaysval
      integer dayofyear,today,month
      integer dayofmonth,day,nextmonth

      integer midday(nmonth+1)
        data (midday(month),month=1,13)
     *    / 16,44,75,105,136,166,197,228,258,289,319,350,381 /

      dd=0
      todaysval=0
      dayofyear=0
      today=0
      month=0
      dayofmonth=0
      day=0
      nextmonth=0

      do month=1,nmonth
        if (month.lt.nmonth) then
          nextmonth=month+1
        else
          nextmonth=1
        endif
        todaysval=mval(month)
        do dayofyear=midday(month),midday(month+1)-1
          if (dayofyear.le.ndayyear) then
            today=dayofyear
          else
            today=dayofyear-ndayyear
          endif
          dval(today)=todaysval
        enddo
      enddo

      return
      end



c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE DAILY
c     Linear interpolation of monthly to quasi-daily values

      subroutine daily(mval,dval)

      implicit none

      integer nmonth
        parameter (nmonth=12)        !number of months per year
      integer ndayyear
        parameter (ndayyear=365)     !number of days per year
      integer ndaymonth(1:12)    !number of days in each month
      integer m
        data (ndaymonth(m),m=1,12)
     *    / 31,28,31,30,31,30,31,31,30,31,30,31 /

c     ARGUMENTS:
      real mval(nmonth),dval(ndayyear)

c     LOCAL VARIABLES:
      real dd,todaysval
      integer dayofyear,today,month
      integer dayofmonth,day,nextmonth

      integer midday(nmonth+1)
        data (midday(month),month=1,13)
     *    / 16,44,75,105,136,166,197,228,258,289,319,350,381 /

      dd=0
      todaysval=0
      dayofyear=0
      today=0
      month=0
      dayofmonth=0
      day=0
      nextmonth=0

      do month=1,nmonth
        if (month.lt.nmonth) then
          nextmonth=month+1
        else
          nextmonth=1
        endif
        dd=(mval(nextmonth)-mval(month))/
     *    real(midday(month+1)-midday(month))
        todaysval=mval(month)
        do dayofyear=midday(month),midday(month+1)-1
          if (dayofyear.le.ndayyear) then
            today=dayofyear
          else
            today=dayofyear-ndayyear
          endif
          dval(today)=todaysval
          todaysval=todaysval+dd
        enddo
      enddo

      return
      end
c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE PRDAILY
c     Distribution of monthly precipitation totals to quasi-daily values

      subroutine prdaily(mval_prec,dval_prec,mval_wet,year)
      implicit none

      integer year
      integer nmonth
        parameter (nmonth=12)        !number of months per year
      integer ndayyear
        parameter (ndayyear=365)     !number of days per year
      integer m,ndaymonth(1:12)      !number of days in each month
        data (ndaymonth(m),m=1,12)
     *    / 31,28,31,30,31,30,31,31,30,31,30,31 /
      integer dayofmonth,day
      real mval_prec(nmonth),dval_prec(ndayyear),mval_wet(nmonth)
      real prob,prob_rain(nmonth)
      real mprec(nmonth)     !average precipitation on wet days
      real mprecip(nmonth)   !acc. monthly generated precipitation
      integer daysum         !accumulated days at beginning of month

c Declarations for the random generator
      integer k1(4),k2(4),k3(4),k4(4),k5(4),k6(4),k7(4),k8(4),k9(4),
     *         k10(4),k11(4)
       real vv,v1,random,c1,c2
       data k1 /9,98,915,92/, k2 /135,28,203,85/, k3 /43,54,619,33/,
     *      k4 /645,9,948,65/, k5 /885,41,696,62/, k6 /51,78,648,0/,
     *      k7 /227,57,929,37/, k8 /20,90,215,31/, k9 /320,73,631,49/,
     *      k10 /73,24,881,52/, k11 /6,302,597,3/
       c1=1.0  !normalising coefficient for exponential distribution
       c2=1.2   !power for exponential distribution

      day=0
      do m=1,nmonth
      prob_rain(m)=0.0
      mprec(m)=0.0
      mprecip(m)=0.0
      enddo
      prob=0.0
      daysum=0.0

      do m=1,nmonth
         if(mval_wet(m).le.1.0) mval_wet(m)=1.0
         prob_rain(m)=mval_wet(m)/real(ndaymonth(m))
         mprec(m)=mval_prec(m)/mval_wet(m)
           do dayofmonth=1,ndaymonth(m)
              day=day+1

c Transitional probabilities (Geng et al. 1986)
                if(dval_prec(day-1).lt.0.1) then
                  prob=0.75*prob_rain(m)
                else
                  prob=0.25+(0.75*prob_rain(m))
                endif

c Determine wet days randomly and use Krysanova/Cramer estimates of
c parameter values (c1,c2) for an exponential distribution
              vv=random(k1)
              if(vv.gt.prob) then
                dval_prec(day)=0.0
              else
                v1=random(k5)
                dval_prec(day)=((-alog(v1))**c2)*mprec(m)*c1
                if(dval_prec(day).lt.0.1) dval_prec(day)=0.0
              endif
              mprecip(m)=mprecip(m)+dval_prec(day)
           enddo

c normalise generated precipitation by monthly CRU values
           if(m.gt.1) daysum=daysum+ndaymonth(m-1)
           if(mprecip(m).lt.1.0) mprecip(m)=1.0
           do dayofmonth=1,ndaymonth(m)
              day=daysum+dayofmonth
              dval_prec(day)=dval_prec(day)*(mval_prec(m)/mprecip(m))
              if (dval_prec(day).lt.0.1) dval_prec(day)=0.0
c              dval_prec(day)=mval_prec(m)/ndaymonth(m)  !no generator
           enddo


c Alternative: equal distribution of rain for fixed number of wet days
c                prob=prob_rain(m)+prob
c                  if(prob.ge.1.0) then
c                    dval_prec(day)=mprec(m)
c                    prob=prob-1.0
c                  else
c                    dval_prec(day)=0.0
c                    prob=prob
c                  endif
      enddo !month
      return
      end

c##############################
c Random function: provides random numbers {0-1}
       real function random(k)
       implicit none

       integer k(4),i

       k(4)=3*k(4)+k(2)
       k(3)=3*k(3)+k(1)
       k(2)=3*k(2)
       k(1)=3*k(1)
       i=k(1)/1000
       k(1)=k(1)-i*1000
       k(2)=k(2)+i
       i=k(2)/100
       k(2)=k(2)-100*i
       k(3)=k(3)+i
       i=k(3)/1000
       k(3)=k(3)-i*1000
       k(4)=k(4)+i
       i=k(4)/100
       k(4)=k(4)-100*i
       random=(((float(k(1))*0.001+float(k(2)))*0.01+float(k(3)))*0.001
     *        +float(k(4)))*0.01

       end
	   

c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE SOILPARAMETERS
c     Obtain soil parameter values given LPJ soil code

      subroutine soilparameters(soilcode,soilpar)

      implicit none

c     PARAMETERS:
      integer nsoilpar
        parameter (nsoilpar=7)
      real k2
c     k2 changed from 2.0 10/2002
        parameter (k2=2.0)

c     ARGUMENTS:
      integer soilcode
      real soilpar(1:nsoilpar)

c     LOCAL VARIABLES:
      integer i,j
      real whc,store(1:5,1:9)

      i=0
      j=0
      whc=0.0
c      store(:,:)=0.0

c     -----------------------------------------------------------------
c     Conversion table from BIOME3 soil codes to LPJ soil code
c     -----------------------------------------------------------------
c     BIOME3 soil(1)  BIOME3 soil(2)  LPJ soil code  Description
c     --------------  --------------  -------------  ------------------

c           1              any              1        coarse
c           2              any              2        medium
c           3               0               3        fine, non-vertisol
c           4              any              4        medium-coarse
c           5              any              5        fine-coarse
c           6              any              6        fine-medium
c           7              any              7        fine-medium-coarse
c           8              any              8        organic
c           3               1               9        fine, vertisol

c     -----------------------------------------------------------------

c     SOIL PARAMETERS

c     1  k1 in percolation formula: Perc=k1*w1**k2
c     2  k2 in percolation formula:
c     3  soil layer 1 water holding capacity (relative)
c     4  soil layer 2 WHC (relative)
c     5  thermal diffusivity at wilting point (see below)
c     6  thermal diffusivity at 15% WHC
c     7  thermal diffusivity at field capacity

c     -----------------------------------------------------------------

c     Array 'store'

c     1  empirical parameter in percolation equation (K1) (mm/day)
c     2  volumetric water holding capacity at field capacity minus vol water
c        holding capacity at wilting point (Hmax), as fraction of soil layer
c        depth
c     3  thermal diffusivity (mm2/s) at wilting point (0% WHC)
c     4  thermal diffusivity (mm2/s) at 15% WHC
c     5  thermal diffusivity at field capacity (100% WHC)
c        Thermal diffusivities follow van Duin (1963),
c        Jury et al (1991), Fig 5.11.

      data ((store(i,j),i=1,5),j=1,9) /

c    ----------------------------------------------------
c           1       2       3       4       5    soilcode
c    ----------------------------------------------------

     *    5.0,   0.110,   0.2,  0.800,    0.4,        ! 1
     *    4.0,   0.150,   0.2,  0.650,    0.4,        ! 2
     *    3.0,   0.120,   0.2,  0.500,    0.4,        ! 3
     *    4.5,   0.130,   0.2,  0.725,    0.4,        ! 4
     *    4.0,   0.115,   0.2,  0.650,    0.4,        ! 5
     *    3.5,   0.135,   0.2,  0.575,    0.4,        ! 6
     *    4.0,   0.127,   0.2,  0.650,    0.4,        ! 7
     *    9.0,   0.300,   0.1,  0.100,    0.1,        ! 8
     *    0.2,   0.100,   0.2,  0.500,    0.4/        ! 9

c     -----------------------------------------------------------------

c     Define water holding capacity for both soil layers (fraction)

c      if (soilcode.lt.1.or.soilcode.gt.9) stop 'soil type not permitted'
      if (soilcode.lt.0) print*,soilcode
      whc = store(2,soilcode)

c     Define k1 value (texture-dependent)
      soilpar(1) = store(1,soilcode)

c     Define k2 value (not texture-dependant)
      soilpar(2) = k2
      soilpar(3) = whc
      soilpar(4) = whc

      soilpar(5)=store(3,soilcode)
      soilpar(6)=store(4,soilcode)
      soilpar(7)=store(5,soilcode)

      return
      end



c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE PFTPARAMETERS
c     Assignment of PFT-specific parameters and bioclimatic limits
c     Definition of initial sapling and grass mass structure
c     Calculation of maximum crown area for woody PFTs

      subroutine pftparameters(pftpar,sla,tree,evergreen,summergreen,
     *  raingreen,needle,boreal,lm_sapl,sm_sapl,hm_sapl,rm_sapl,
     *  latosa,allom1,allom2,allom3,allom4,wooddens,reinickerp
     *  ,co2, BTparam1, BTparam2, BTmode0) ! Doug 02/13: Add BT params

      implicit none

c     PARAMETERS:
      integer npft,npftpar,nsoilpar
        parameter (npft=13,npftpar=59,nsoilpar=7)
      integer nco2
        parameter (nco2=3)           !number of C variables: C,13C,14C
      real pi
        parameter (pi=3.14159265)

c     ARGUMENTS:
      real pftpar(1:npft,1:npftpar),sla(1:npft)
      logical tree(1:npft),evergreen(1:npft)
      logical summergreen(1:npft),raingreen(1:npft),needle(1:npft)
      logical boreal(1:npft)
      real lm_sapl(1:npft,1:nco2),sm_sapl(1:npft,1:nco2)
      real hm_sapl(1:npft,1:nco2),rm_sapl(1:npft,1:nco2)
      real latosa
      real allom1,allom2,allom3,allom4
      real wooddens,reinickerp
      real co2(1:nco2)
      REAL BTparam1(1:npft,1:3)
      REAL BTparam2(1:npft,1:3)
      REAL BTmode0(1:npft,1:2)

c     LOCAL VARIABLES:
      integer n,pft
      real table(1:npft,1:npftpar)
      real lai_sapl       !sapling or initial grass LAI
      real x
      real lmtorm         !non-waterstressed leafmass to rootmass ratio
      real stemdiam       !sapling stem diameter
      real height_sapl    !sapling height

      n=0
      pft=0
      lai_sapl=0.0
      x=0.0
      lmtorm=0.0
      stemdiam=0.0
      height_sapl=0.0
      BTparam1(:,:)=0.0 !Doug 02/13
      BTparam2(:,:)=0.0 !Doug 02/13
      BTmode0(:,:) =0.0 !Doug 02/13

c-----------------------------------------------------------------------------

c     PFT PARAMETERS

c      1  fraction of roots in upper soil layer
c      2  plants with C4 (1) or C3 (0) photosynthetic pathway
c      3  water scalar value at which leaves shed by drought deciduous PFT
c      4  canopy conductance component (gmin, mm/s) not associated with
c         photosynthesis (Haxeltine & Prentice 1996, Table 4)
c      5  maintenance respiration coefficient
c      6  flammability threshold
c      7  maximum foliar N content (mg/g)
c         (Haxeltine & Prentice 1996a, Fig 4)
c      8  fire resistance nindex
c      9  leaf turnover period (years)
c     10  leaf longevity (years)
c     11  sapwood turnover period (sapwood converted to heartwood) (years)
c     12  root turnover period (years)
c     13  leaf C:N mass ratio
c     14  sapwood C:N mass ratio
c     15  root C:N mass ratio
c     16  leaf type: broadleaved (1), needleleaved (2) or grass (3)
c     17  phenology type: evergreen (1), summergreen (2), raingreen (3),
c         any type (4)
c     18  leaf to root ratio under non-water stressed conditions
c     19  summergreen phenology ramp, GDD requirement to grow full leaf canopy
c     20  tree maximum crown area (m2)
c     21  sapling (or grass on initialisation) LAI
c     22  sapling [(heartwood mass) + (sapwood mass)] / (sapwood mass)
c     23  boreal pft (1), non-boreal pft (0)
c     24  low temperature limit for CO2 uptake
c     25  lower range of temperature optimum for photosynthesis
c     26  upper range of temperature optimum for photosynthesis
c     27  high temperature limit for CO2 unptake

c     BIOCLIMATIC LIMITS

c     28 minimum coldest monthly mean temperature
c     29 maximum coldest monthly mean temperature
c     30 minimum growing degree days (at or above 5 deg C)
c     31 upper limit of temperature of the warmest month
c     32 lower limit of growth efficiency (g/m2)
c
c     PARAMETERS ADDED LATER
c
c     33 GDD base
c     34 20-year average min warmest - coldest month temperature range
c     35 emax (mm/day)
c     36 intc (Interception storage parameter, unitless)
c     37 fuel bulk density  (kg/m?³)
c     38 emission factor CO2
c     39 emission factor CO
c     40 emission factor CH4
c     41 emission factor VOC
c     42 emission factor TPM
c     43 emission factor NOx
c     44 proportion of tree height as crown
c     45 f parameter for flame length
c     46 g parameter for flame length  Kirsten, 4.11.04: not used at the moment
c     47 Doug 02/13: lower param1 for bark thickness
c     48 Doug 02/13: median param1 for bark thickness
c     49 Doug 02/13: upper param1 for bark thickness
c     50 Doug 02/13: lower param2 for bark thickness
c     51 Doug 02/13: median param2 for bark thickness
c     52 Doug 02/13: upper param2 for bark thickness
c	  Doug 02/13: parameters re-numbered
c     53 r(ck) crown damage: postfire mortality
c     54 p   crown damage: postfire mortality
c	  56 Doug 11/12: decomposition rate of leaf at 10 degrees  (kleaf10)
c	  57 Doug 11/12: decomposition rate of wood at 10 degrees  (kwood10)
c     58 Doug 11/12: Q10 for wood decomposition
c     59 Doug 02/03: % of resprouting success post-fire in pft in pft


c NOTE ON THE INTRODUCTION OF THE LARCH PFT:
c assume that Larch canopy conductance is higher than for other boreals
c assume that phenology ramp is quicker than for broadleaved trees
c also that ramp is applied to a GDD base of 2, not 5
      data ((table(pft,n),n=1,8),pft=1,npft) /

c     ---------------------------------------------------------------------
c          1      2      3      4      5      6      7      8          PFT
c     ---------------------------------------------------------------------


     *  0.80,   0.0,  0.00,   0.5,  0.20,  0.25, 100.0,  0.60,        !  1
     *  0.70,   0.0,  0.10,   0.5,  0.20,  0.25, 100.0,  0.70,        !  2
     *  0.85,   0.0,  0.00,   0.3,  1.20,  0.25, 100.0,  0.12,        !  3
     *  0.80,   0.0,  0.10,   0.3,  1.20,  0.25, 100.0,  0.50,        !  4
     *  0.80,   0.0,  0.00,   0.5,  1.20,  0.25, 120.0,  0.12,        !  5
     *  0.85,   0.0,  0.00,   0.3,  1.20,  0.25, 100.0,  0.12,        !  6
c     *  0.90,   0.0,  0.00,   0.3,  1.30,  0.35, 100.0,  0.12,        !  7
     *  0.80,   0.0,  0.00,   0.5,  1.20,  0.3, 100.0,  0.12,        !  8
     *  0.90,   0.0,  0.20,   0.5,  1.30,  0.25, 100.0,  0.01,        !  9
     *  0.85,   1.0,  0.20,   0.5,  0.70,  0.25, 100.0,  0.01,        ! 10
     *  0.80,   0.0,  0.10,   0.5,  0.20,  0.25, 100.0,  0.70,        !  1r     
     *  0.70,   0.0,  0.10,   0.5,  0.20,  0.25, 100.0,  0.70,        !  2r
     *  0.80,   0.0,  0.10,   0.3,  1.20,  0.25, 100.0,  0.50,        !  4r
     *  0.80,   0.0,  0.00,   0.5,  1.20,  0.25, 120.0,  0.12/        !  5r


      data ((table(pft,n),n=9,17),pft=1,npft) /

c     ---------------------------------------------------------------------
c          9     10     11     12     13     14     15     16     17   PFT
c     ---------------------------------------------------------------------

     *   2.0,  2.00,  20.0,   2.0,  29.0, 330.0,  29.0,   1.0,   1.0, !  1
     *   1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   3.0, !  2
     *   4.0,  4.00,  20.0,   4.0,  29.0, 330.0,  29.0,   2.0,   1.0, !  3
     *   1.0,  1.00,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   1.0, !  4
     *   1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   2.0, !  5
     *   4.0,  4.00,  20.0,   4.0,  29.0, 330.0,  29.0,   2.0,   1.0, !  6
c     *   1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   2.0,   2.0, !  7
     *   1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   2.0, !  8
     *   1.0,  0.50,   1.0,   2.0,  29.0,   0.0,  29.0,   3.0,   4.0, !  9
     *   1.0,  0.50,   1.0,   2.0,  29.0,   0.0,  29.0,   3.0,   4.0, ! 10
     *   2.0,  2.00,  20.0,   2.0,  29.0, 330.0,  29.0,   1.0,   1.0, !  1r
     *   1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   3.0, !  2r
     *   1.0,  1.00,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   1.0, !  4r
     *   1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   2.0/ !  5r

      data ((table(pft,n),n=18,23),pft=1,npft) /

c     ---------------------------------------------------
c           18      19     20      21     22     23     pft
c     ---------------------------------------------------

     *     1.0,  1000.0,  15.0,  1.500,  1.2,    0.0,    !  1
     *     1.0,  1000.0,  15.0,  1.500,  1.2,    0.0,    !  2
     *     1.0,  1000.0,  15.0,  1.500,  1.2,    0.0,    !  3
     *     1.0,  1000.0,  15.0,  1.500,  1.2,    0.0,    !  4
     *     1.0,   200.0,  15.0,  1.500,  1.2,    0.0,    !  5
     *     1.0,  1000.0,  15.0,  1.500,  1.2,    1.0,    !  6
c     *     1.0,   100.0,  15.0,  1.500,  1.2,    1.0,    !  7
     *     1.0,   200.0,  15.0,  1.500,  1.2,    1.0,    !  8
     *    0.65,   100.0,   0.0,  0.001,  1.2,    1.0,    !  9
     *    0.75,   100.0,   0.0,  0.001,  1.2,    0.0,    ! 10
     *     1.0,  1000.0,  15.0,  1.500,  1.2,    0.0,    !  1r
     *     1.0,  1000.0,  15.0,  1.500,  1.2,    0.0,    !  2r
     *     1.0,  1000.0,  15.0,  1.500,  1.2,    0.0,    !  4r
     *     1.0,   200.0,  15.0,  1.500,  1.2,    0.0/    !  5r

      data ((table(pft,n),n=24,27),pft=1,npft) /
c     -------------------------------------
c          24     25     26      27    PFT
c     -------------------------------------
     *    2.0,   25.0,  30.0,   55.0, !  1
     *    2.0,   25.0,  30.0,   55.0, !  2
     *   -4.0,   20.0,  30.0,   42.0, !  3
     *   -4.0,   20.0,  30.0,   42.0, !  4
     *   -4.0,   20.0,  25.0,   38.0, !  5
     *   -4.0,   15.0,  25.0,   38.0, !  6
c     *   -4.0,   15.0,  25.0,   38.0, !  7
     *   -4.0,   15.0,  25.0,   38.0, !  8
     *   -4.0,   10.0,  30.0,   45.0, !  9
     *    6.0,   20.0,  45.0,   55.0, ! 10
     *    2.0,   25.0,  30.0,   55.0, !  1r
     *    2.0,   25.0,  30.0,   55.0, !  2r
     *   -4.0,   20.0,  30.0,   42.0, !  4r
     *   -4.0,   20.0,  25.0,   38.0/ !  5r


      data ((table(pft,n),n=28,32),pft=1,npft) /

c     ------------------------------------------------------------
c          28       29       30       31       32       PFT
c     ------------------------------------------------------------
     *    15.5,  1000.0,    0.0,    1000.0,    0.0,     !  1
     *    13.5,  1000.0,    0.0,    1000.0,    0.0,     !  2
     *    -2.0,    15.0,  900.0,    1000.0,    0.0,     !  3
     *     1.0,    15.5, 1200.0,    1000.0,    0.0,     !  4
     *   -17.0,    15.5, 1200.0,    1000.0,    0.0,     !  5
     *   -32.5,    -4.0,  600.0,    1000.0,    0.0,     !  6
c     * -1000.0,    -2.0,  350.0,      23.0,    0.0,     !  7
     * -1000.0,    -4.0,  350.0,    1000.0,    0.0,     !  8
     *   -1000,  1000.0,    0.0,    1000.0,    0.0,     !  9
     *    10.0,  1000.0,    0.0,    1000.0,    0.0,     ! 10
     *    15.5,  1000.0,    0.0,    1000.0,    0.0,     !  1r
     *    13.5,  1000.0,    0.0,    1000.0,    0.0,     !  2r
     *     1.0,    15.5, 1200.0,    1000.0,    0.0,     !  4r
     *   -17.0,    15.5, 1200.0,    1000.0,    0.0/     !  5r


      data ((table(pft,n),n=33,37),pft=1,npft) /

c     ---------------------------------------------
c          33       34        35    36     37     PFT
c     ---------------------------------------------
!     *     5.0,  -1000.0,     7.0,  0.02, 25.0,   !  1
!     *     5.0,  -1000.0,     7.0,  0.02, 25.0,   !  2
!     *     5.0,  -1000.0,     5.0,  0.06, 20.0,   !  3
!     *     5.0,  -1000.0,     5.0,  0.02, 16.0,   !  4
!     *     5.0,  -1000.0,     5.0,  0.02, 22.0,   !  5
!     *     5.0,  -1000.0,     5.0,  0.06, 18.0,   !  6
c     *     2.0,   43.0,       5.0,  0.06, 16.0,   !  7
!     *     5.0,  -1000.0,     6.0,  0.06, 16.0,   !  8
!     *     5.0,  -1000.0,     5.0,  0.01,  2.0,   !  9
!     *     5.0,  -1000.0,     7.0,  0.01,  2.0/   ! 10

     *     5.0,  -1000.0,     7.0,  0.02, 10.0,   !  1
     *     5.0,  -1000.0,     7.0,  0.02, 10.0,   !  2
     *     5.0,  -1000.0,     5.0,  0.06, 16.0,   !  3
     *     5.0,  -1000.0,     5.0,  0.02, 10.0,   !  4
     *     5.0,  -1000.0,     5.0,  0.02, 13.0,   !  5
     *     5.0,  -1000.0,     5.0,  0.06, 18.0,   !  6
c     *     2.0,   43.0,       5.0,  0.06, 16.0,   !  7
     *     5.0,  -1000.0,     6.0,  0.06, 16.0,   !  8
     *     5.0,  -1000.0,     5.0,  0.01,  0.5,   !  9
     *     5.0,  -1000.0,     7.0,  0.01,  2.0,   ! 10
     *     5.0,  -1000.0,     7.0,  0.02, 10.0,   !  1r
     *     5.0,  -1000.0,     7.0,  0.02, 10.0,   !  2r
     *     5.0,  -1000.0,     5.0,  0.02, 10.0,   !  4r
     *     5.0,  -1000.0,     5.0,  0.02, 13.0/   !  5r

      data ((table(pft,n),n=38,45),pft=1,npft) /

c     --------------------------------------------------------------------
c          38       39      40      41      42   43       44      45    PFT
c     --------------------------------------------------------------------
     *    1580.0, 103.0,   6.80,   8.10,   8.50, 1.999, 0.3334, 0.160, !  1
     *    1664.0,  63.0,   2.20,   3.40,   8.50, 2.540, 0.10,   0.351, !  2
     *    1568.0, 106.0,   4.80,   5.70,  17.60, 3.240, 0.3334, 0.094, !  3
     *    1568.0, 106.0,   4.80,   5.70,  17.60, 3.240, 0.3334, 0.070, !  4
     *    1568.0, 106.0,   4.80,   5.70,  17.60, 3.240, 0.3334, 0.094, !  5
     *    1568.0, 106.0,   4.80,   5.70,  17.60, 3.240, 0.6667, 0.094, !  6
c     *     2.0,   43.0,       5.0,  0.06,   !  7
     *    1568.0, 106.0,   4.80,   5.70,  17.60, 3.240, 0.3334, 0.094, !  8
     *    1568.0, 106.0,   4.80,   5.70,  17.60, 3.240, 0.0,    0.0,   !  9
     *    1664.0,  63.0,   2.20,   3.40,   8.50, 2.540, 0.0,    0.0,   ! 10
     *    1580.0, 103.0,   6.80,   8.10,   8.50, 1.999, 0.3334, 0.160, !  1r
     *    1664.0,  63.0,   2.20,   3.40,   8.50, 2.540, 0.10,   0.351, !  2r
     *    1568.0, 106.0,   4.80,   5.70,  17.60, 3.240, 0.3334, 0.070, !  4r
     *    1568.0, 106.0,   4.80,   5.70,  17.60, 3.240, 0.3334, 0.094/ !  5r
c----------------------------------------------------------------------------

      DATA ((table(pft,n),n=46,52),pft=1,npft) /

c     -------------------------------------------------
c         46		47      48    	 49   	50    		 51  		 52    		PFT
c     -------------------------------------------------
     *    0.630,    0.0,	0.00001,0.00002,0.00395,   0.0167,   0.0399,   !  1
     *    0.460,    0.0,	0.00001,0.00002,0.00463,   0.0194,    0.0571,   !  2
     *    0.667,    0.0,	0.00001,0.00002,0.00609,   0.0257,   0.0576,   !  3
     *    0.560,    0.0,	0.00001,0.00002,0.0125,     0.0302,    0.0909,   !  4
     *    0.667,    0.0,	0.00001,0.00002,0.00617,   0.023,   0.0559,   !  5
     *    0.667,    0.0,	0.00001,0.00002,0.0158,    0.0261,  0.0529,   !  6
c     *     2.0,   43.0,       5.0,  0.06,   !  7
     *    0.667,    0.0,	0.00001,0.00002,0.00875,   0.0316,   0.112,   !  8
     *    0.0,      0.0,	0.0,	0.0,    0.0,       0.0,	     0.0,     !  9
     *    0.0,      0.0,	0.0,	0.0,    0.0,       0.0,	     0.0,  ! 10
     *    0.630,    0.0,	0.00001,0.00002,0.0292,   0.0629,    0.182,   !  1r
     *    0.460,    0.0,	0.00001,0.00002,0.0109,   0.0568,    0.188,   !  2r
     *    0.560,    0.0,	0.00001,0.00002,0.0286,   0.0586,    0.156,   !  4r
     *    0.667,    0.0,	0.00001,0.00002,0.0106,   0.0343,    0.106/   !  5r
c----------------------------------------------------------------------------

      DATA ((table(pft,n),n=53,npftpar),pft=1,npft) /

c     -------------------------------------------------
c         53	54	   55	56	    57    58     59     PFT
c     -------------------------------------------------
     *    0.95, 3.00, 0.26, 0.93,  0.039, 2.75,  0.0,  !  1
     *    0.05, 3.00, 0.25, 1.17,  0.039, 2.75,  0.0,  !  2
     *    0.95, 3.75, 0.03, 0.70,  0.041, 1.97,  0.0,  !  3
     *    0.95, 3.00, 0.03, 0.86,  0.104, 1.37,  0.0,  !  4
     *    0.95, 3.00, 0.03, 0.95,  0.104, 1.37,  0.0,  !  5
     *    0.95, 3.00, 0.58, 0.776, 0.041, 1.97,  0.0,  !  6
c     *     2.0,   43.0,       5.0,  0.06,   !  7
     *    0.95, 3.00, 0.18, 0.94,  0.104, 1.37,  0.0,  !  8
     *    0.00, 0.00, 0.05, 1.20,  0.0  , 0.0 ,  0.0,  !  9
     *    0.00, 0.00, 0.06, 0.97,  0.0  , 0.0 ,  0.0,  ! 10
     *    0.95, 3.00, 0.26, 0.93,  0.039, 2.75,  1.0,  !  1r
     *    0.05, 3.00, 0.25, 1.17,  0.039, 2.75,  1.0,  !  2r
     *    0.95, 3.00, 0.03, 0.86,  0.104, 1.37,  1.0,  !  4r
     *    0.95, 3.00, 0.03, 0.95,  0.104, 1.37,  1.0/  !  5r
c----------------------------------------------------------------------------

      do pft=1,npft

c       Transfer parameter values to array pftpar

        do n=1,npftpar
          pftpar(pft,n)=table(pft,n)
        enddo

c       Assign leaf and phenology logicals

        if (pftpar(pft,16).le.2.0) then
          tree(pft)=.true.
          if (pftpar(pft,16).eq.2.0) then
            needle(pft)=.true.
          else
            needle(pft)=.false.
          endif
        else
          tree(pft)=.false.
          needle(pft)=.false.
        endif

        if (pftpar(pft,17).eq.1.0) then
          evergreen(pft)=.true.
          summergreen(pft)=.false.
          raingreen(pft)=.false.
        elseif (pftpar(pft,17).eq.2.0) then
          evergreen(pft)=.false.
          summergreen(pft)=.true.
          raingreen(pft)=.false.
        elseif (pftpar(pft,17).eq.3.0) then
          evergreen(pft)=.false.
          summergreen(pft)=.false.
          raingreen(pft)=.true.
        else
          evergreen(pft)=.true.
          summergreen(pft)=.true.
          raingreen(pft)=.true.
        endif

        if (pftpar(pft,23).eq.1.0) then
          boreal(pft)=.true.
        else
          boreal(pft)=.false.
        endif

c       Calculate specific leaf area (SLA) for each PFT from leaf longevity
c       Include conversion (multiplier of 2.0) from m2/g(dry wt) to m2/gC
c       Equation based on Reich et al 1997, Fig 1f:

c       SLA = 2e-4 * exp(6.15 - 0.46 ln (leaf_longevity * 12))

c       SLA in m2/gC, leaf_longevity in years

        sla(pft)=2.0e-4*exp(6.15-0.46*log(pftpar(pft,10)*12.0))

c       Define initial mass structure

        lai_sapl=pftpar(pft,21)

        if (tree(pft)) then  !woody PFTs

c         Calculate leafmass for a sapling individual
c          (1) lai = leafmass * sla / (crown area)
c          (2) (leaf area) = latosa * (sapwood xs area)
c                 (Pipe Model, Shinozaki et al. 1964a,b; Waring et al 1982)
c          (3) (crown area) = allom1 * (stem diameter) ** reinickerp
c                 (Reinickes theory)
c         From (1),
c          (4) leafmass = lai * (crown area) / sla
c         From (1) & (3),
c          (5) leafmass = lai * allom1 * (stem diameter)**reinickerp / sla
c         From (2),
c          (6) leafmass = latosa * (sapwood xs area) / sla
c          (7) (sapwood xs area) = pi * (sapwood diameter)**2 / 4
c         From (6) and (7),
c          (8) leafmass = latosa * pi * (sapwood diameter)**2 / 4 / sla
c         From (8),
c          (9) (sapwood diameter) = [ 4 * leafmass * sla / pi / latosa ]**0.5
c         (10) (stem diameter) = (sapwood diameter) + (heartwood diameter)
c         Define x,
c         (11) x = [ (sapwood diameter)+(heartwood diameter) ] /
c                  (sapwood diameter)
c         From (10) & (11),
c         (12) (stem diameter) = x * (sapwood diameter)
c         From (5), (9) & (12),
c         (13) leafmass = lai * allom1 * x**reinickerp *
c                       (4*leafmass*sla/pi/latosa)**(reinickerp*0.5) / sla
c         From (13),
c         (14) leafmass = [ lai * allom1 * x**reinickerp *
c              (4*sla/pi/latosa)**(reinickerp*0.5) / sla ]**(2/(2-reinickerp))

          x=pftpar(pft,22)
	 
          lm_sapl(pft,1)=(lai_sapl*allom1*x**reinickerp*
     *      (4.0*sla(pft)/pi/latosa)**(reinickerp*0.5) / sla(pft))**
     *      (2.0/(2.0-reinickerp))  !eqn 14

          lm_sapl(pft,2)=co2(2)-17.8 ! initial 13C value from llyod & farquhar, 1994
          lm_sapl(pft,3)=co2(3)-35.6 ! 2*17.8 see above

c         Calculate sapling stem diameter
c         From (9) & (12),
c         (15) (stem diameter) = x * [ 4 * leafmass * sla / pi / latosa ]**0.5

          stemdiam=x*(4.0*lm_sapl(pft,1)*sla(pft)/pi/latosa)**0.5  !Eqn 15

c         Calculate sapling height
c         (16) height = allom2 * (stem diameter)**allom3 (source?)

          height_sapl=allom2*stemdiam**allom3   !Eqn 16

c         Calculate sapling sapwood mass
c         (17) (sapwood volume) = height * (sapwood xs area)
c         (18) (sapwood xs area) = leafmass * sla / latosa
c         From (17) & (18),
c         (19) (sapwood volume) = height * leafmass * sla / latosa
c         (20) (sapwood mass) = (wood density) * (sapwood volume)
c         From (19) & (20),
c         (21) (sapwood mass) = (wood density) * height * leafmass * sla /
c                latosa

          sm_sapl(pft,1)=wooddens*height_sapl*lm_sapl(pft,1)
     *      *sla(pft)/latosa   !Eqn 21
          sm_sapl(pft,2)=lm_sapl(pft,2) ! 13C value in permille
          sm_sapl(pft,3)=lm_sapl(pft,3)

c         Calculate sapling heartwood mass
c         From (11),
c         (22) (heartwood mass) = (x-1) * (sapwood mass)

          hm_sapl(pft,1)=(x-1.0)*sm_sapl(pft,1)  !Eqn 22
          hm_sapl(pft,2)=sm_sapl(pft,2) ! 13C value in permille
          hm_sapl(pft,3)=sm_sapl(pft,3)

        else ! grass PFT
c     * lm_sapl(pft,1),lai_sapl,sla(pft)
          lm_sapl(pft,1)=lai_sapl/sla(pft)

c         Set initial 13C values for saplings, grass

          if (pftpar(pft,2).eq.1) then   !C4 plants
            lm_sapl(pft,2)=co2(2)-3.6  !from lloyd & farquhar,1994           
            lm_sapl(pft,3)=co2(3)-7.2  ! 2*3.6
          else                           !C3 plpants
            lm_sapl(pft,2)=co2(2)-17.8  !from lloyd & farquhar,1994          
            lm_sapl(pft,3)=co2(3)-35.6  !2*17.8
          endif

          sm_sapl(pft,2)=0.0             ! no sapwood and hartwood
          hm_sapl(pft,2)=0.0             ! for grass PFT
          sm_sapl(pft,3)=0.0             ! no sapwood and hartwood
          hm_sapl(pft,3)=0.0             ! for grass PFT

        endif

c       Calculate sapling or initial grass rootmass
c       (23) lmtorm = (leafmass) / (rootmass)

        lmtorm=pftpar(pft,18)
        rm_sapl(pft,1)=(1.0/lmtorm)*lm_sapl(pft,1)  !From Eqn 23
        rm_sapl(pft,2)=lm_sapl(pft,2) ! 13C value in permille
        rm_sapl(pft,3)=lm_sapl(pft,3) 
		
c        Doug 02/13: Assign lower, median and upper BT paramerter values for BT=p1+p2*DBH
         DO n=1,3
           BTparam1(pft,n)=pftpar(pft,46+n)
           BTparam2(pft,n)=pftpar(pft,49+n)
         END DO
         BTmode0(pft,1)=BTparam1(pft,2)
         BTmode0(pft,2)=BTparam2(pft,2)

      enddo ! pft loop
      
      return
      end


c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE INITGRID
c     Initialise state variables (as required) for current grid cell

      subroutine initgrid(tree,k_fast_ave,k_slow_ave,
     *  litter_decom_ave,present,litter_ag,
     *  litter_ag_leaf,litter_ag_wood,
     *  fuel_1hr_leaf,fuel_1hr_wood,fuel_10hr,
     *  fuel_100hr,fuel_1000hr,litter_bg,crownarea,mprec,w,w_t,
     *  dwscal365,lm_ind,sm_ind,hm_ind,rm_ind,fpc_grid,mnpp,anpp,
     *  leafondays,leafoffdays,leafon,snowpack,mtemp_old,mtemp,
     *  maxthaw_old,mw1,mw2,mw1_t,mw2_t,uw1,uw2,fw,soilpar,
     *  mcica,mgpp,lresp,sresp,rresp,gresp,aresp,dbh,
     *  tau_c,dbh_class, !Doug 02/13
     *  cl_t,height_class,agpp)

      implicit none

c     PARAMETERS
      integer npft,npftpar,nsoilpar
        parameter (npft=13,npftpar=59,nsoilpar=7)
      integer nco2
        parameter (nco2=3)           !number of C variables: C,13C,14C
      real rainscalar
        parameter (rainscalar=2500.)
      real d1,d2
        parameter (d1=500.0,d2=1000.0)


c     ARGUMENTS
      logical tree(1:npft)
      real k_fast_ave,k_slow_ave
      real litter_decom_ave(1:nco2)
      logical present(1:npft)
      REAL litter_ag(1:npft,1:nco2)
      REAL litter_ag_leaf(1:npft,1:nco2) !Doug 11/12:
      REAL litter_ag_wood(1:npft,1:nco2)
      real fuel_1hr_leaf(1:npft,1:nco2)
      REAL fuel_1hr_wood(1:npft,1:nco2)
      REAL fuel_10hr(1:npft,1:nco2)
      real fuel_100hr(1:npft,1:nco2)
      real fuel_1000hr(1:npft,1:nco2)
      real litter_bg(1:npft,1:nco2)
      real crownarea(1:npft)
      real mprec(1:12)
      real w(1:2),w_t(1:2)
      real dwscal365(1:npft)
      real lm_ind(1:npft,1:nco2)
      real sm_ind(1:npft,1:nco2)
      real hm_ind(1:npft,1:nco2)
      real rm_ind(1:npft,1:nco2)
      real fpc_grid(1:npft)
      real mnpp(1:12,1:npft,1:nco2)          !monthly gridcell NPP (gC/m2)
      real anpp(1:npft,1:nco2)
      integer leafondays(1:npft),leafoffdays(1:npft)
      logical leafon(1:npft)
      real snowpack
      real mtemp_old(1:12),mtemp(1:12)
      real maxthaw_old
      real mw1(1:12),mw2(1:12)        ! monthly values of w(2) i.e.fraction of avilable water
      real mw1_t(1:12),mw2_t(1:12)
      real uw1,uw2,fw(1:2)
      
      real soilpar(1:nsoilpar)        !soil parameters
      real mcica(1:12,1:npft)
      real mgpp(1:12,1:npft,1:nco2)
      real lresp(1:12,1:npft,1:nco2)  ! monthly leaf respiration
      real sresp(1:12,1:npft,1:nco2)  ! monthly sapwood respiration
      real rresp(1:12,1:npft,1:nco2)  ! monthly root respiration
      real gresp(1:12,1:npft,1:nco2)  ! monthly growth respiration
      real aresp(1:12,1:npft,1:nco2)  ! monthly autotrophic respiration
      REAL dbh(1:npft), dbh_class(0:4,1:npft) !Doug 02/13 dbh_class added
      REAL tau_c(0:4,1:npft),cl_t(0:4,1:npft)
      real height_class(0:4,1:npft)
      real agpp(1:npft,1:nco2)


c     LOCAL VARIABLES
      integer class
      integer pft,m,i
      integer nc               ! nc is added for nco2
      real aprec,air_temp

      class=0
      pft=0
      m=0
      i=0
      nc=0
      aprec=0
      air_temp=0

      k_fast_ave=0.0
      k_slow_ave=0.0
      litter_decom_ave(1)=0.0
      litter_decom_ave(2)=0.0
      litter_decom_ave(3)=0.0

      snowpack=0.0

      do pft=1,npft

        present(pft)=.false.
        do nc=1,nco2
          litter_ag(pft,nc)=0.0
          litter_ag_leaf(pft,nc)=0.0 ! Doug 11/12
          litter_ag_wood(pft,nc)=0.0 ! Doug 11/12
          litter_bg(pft,nc)=0.0
          fuel_1hr_leaf(pft,nc)=0.0
          fuel_1hr_wood(pft,nc)=0.0
          fuel_10hr(pft,nc)=0.0
          fuel_100hr(pft,nc)=0.0
          fuel_1000hr(pft,nc)=0.0
        enddo
        dbh(pft)=0.0
        do class=0,4
          tau_c(class,pft)=0.0
          dbh_class(class,pft)=0.0 ! Doug 02/13
          cl_t(class,pft)=0.0
          height_class(class,pft)=0.0
        enddo
        crownarea(pft)=0.0
        if (.not.tree(pft)) crownarea(pft)=1.0
        dwscal365(pft)=1.0
        do nc=1,nco2
          lm_ind(pft,nc)=0.0
          sm_ind(pft,nc)=0.0
          hm_ind(pft,nc)=0.0
          rm_ind(pft,nc)=0.0
        enddo
        fpc_grid(pft)=0.0
        leafondays(pft)=0
        leafoffdays(pft)=0
        leafon(pft)=.true.
        do m=1,12
          do nc=1,nco2
            mnpp(m,pft,nc)=0.0
            mgpp(m,pft,nc)=0.0
            lresp(m,pft,nc)=0.0
            sresp(m,pft,nc)=0.0
            rresp(m,pft,nc)=0.0
            gresp(m,pft,nc)=0.0
            aresp(m,pft,nc)=0.0
          enddo  
          mcica(m,pft)=0.0
        enddo
        do nc=1,nco2
          anpp(pft,nc)=0.0
          agpp(pft,nc)=0.0
        enddo
      enddo
      
      do m=1,12
        aprec=aprec+mprec(m)
        mtemp_old(m)=mtemp(m)
      enddo

c     initialize permafrost
      maxthaw_old=0.0
c DM      delta_thawday=0.0


c     initialise soil water buckets
      do i=1,12
        air_temp=air_temp+mtemp(i)/12.0
      enddo

      if ((air_temp.le.0.0.and.mtemp(1).le.0.0)) then
        w(1)=max(aprec/rainscalar,0.1)
        if(w(1).ge.1.0) w(1)=1.0
        w_t(1)=max(aprec/rainscalar,0.1)
        if(w_t(1).ge.1.0) w_t(1)=1.0
        w(2)=1.0
        w_t(2)=1.0
        do m=1,12
          mw1(m)=w(1)
          mw2(m)=w(2)
          mw1_t(m)=w_t(1)
          mw2_t(m)=w_t(2)
        enddo
        uw1=w(1)*(d1*soilpar(3))    ! This is the water content in soil layer 1 (mm).
        uw2=w(2)*(d2*soilpar(4))    ! This is the water content in soil layer 2 (mm).
        fw(1)=0.0
        fw(2)=0.0

      else
        w(1)=max(aprec/rainscalar,0.1)
        if(w(1).ge.1.0) w(1)=1.0
        w_t(1)=max(aprec/rainscalar,0.1)
        if(w_t(1).ge.1.0) w_t(1)=1.0
        w(2)=1.0
        w_t(2)=1.0
c DM   TODO : is it w(2)=1.0 or w(2)=w(1) ??
        w(2)=w(1)
        w_t(2)=w_t(1)
        do m=1,12
          mw1(m)=w(1)
          mw2(m)=w(2)
          mw1_t(m)=w_t(1)
          mw2_t(m)=w_t(2)
        enddo
        uw1=w(1)*(d1*soilpar(3))
        uw2=w(2)*(d2*soilpar(4))
        fw(1)=0.0
        fw(2)=0.0
      endif

      return
      end
c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE PETPAR
c     Calculates mid-month daily photosynthetically active radiation flux,
c     daylength and daily potential evapotranspiration given temperature,
c     sunshine percentage and latitude
c     Doug 08/09: Added equilibrum evapotranspiration calculation.

c#ifdef LPJ_STEP_2
c      subroutine petpar(year,dtemp,dsun,lat,mdayl,mpar_day,mpet,
c     *      u,v,sinelat,cosinelat,sinehh,hh,qo,spinup_years)
c#else
c      subroutine petpar(year,dtemp,dsun,lat,mdayl,mpar_day,mpet,
c     *      u,v,sinelat,cosinelat,sinehh,hh,qo)
c#endif

      subroutine petpar(year,dtemp,dsun,lat,mdayl,mpar_day,mpet,
     *     deet,    !Doug 08/09
     *      u,v,sinelat,cosinelat,sinehh,hh,qo,start_year)

      implicit none

c     PARAMETERS:
      integer npft,npftpar,nsoilpar
        parameter (npft=13,npftpar=59,nsoilpar=7)
      real pi
        parameter (pi=3.14159265)
      real a,b                          !empirical constants (Eqn 18)
        parameter (a=107.0,b=0.2)       !(a in W/m2; b unitless)
      real qoo
        parameter (qoo=1360.0)          !solar constant (1360 W/m2)
      real beta                         !global average short-wave albedo
        parameter (beta=0.17)           !(fraction)
      real c,d                          !empirical constants (Eqn 1)
        parameter (c=0.25,d=0.5)        !(unitless)
      real degtorad                     !conversion factor from degrees to
        parameter (degtorad=0.01745329) !radians (pi/180)
      real k                            !conversion factor from solar angular
        parameter (k=13750.98708)       !units to seconds (12/pi*3600)

c     ARGUMENTS:
      integer year
      real dtemp(1:365),dsun(1:365),lat
      real mdayl(1:12),mpar_day(1:12)
      real mpet(1:12)
      REAL deet(1:365)                  !Doug 08/09: daily equilibrium ET
      real u(1:12),v(1:12)
      real sinelat,cosinelat
      real sinehh(1:12),hh(1:12),qo(1:12)
c#ifdef LPJ_STEP_2
c      integer spinup_years
c#endif
      integer start_year    !Doug 07/09


c     LOCAL VARIABLES:
c      real dpet(1:365)
      integer m,dm,day
      real ni,delta,w,rs_day
      real rl,gamma,lambda,s,uu,vv,hn
      
      integer ndaymonth(1:12)
      data (ndaymonth(m),m=1,12)
     *  / 31,28,31,30,31,30,31,31,30,31,30,31 /
      integer midday(1:12)
      data (midday(m),m=1,12)
     *  / 16,44,75,105,136,166,197,228,258,289,319,350 /

c      dpet(:)=0.0
      m=0
      dm=0
      day=0
      ni=0.0
      delta=0.0
      w=0.0
      rs_day=0.0
      rl=0.0
      gamma=0.0
      lambda=0.0
      s=0.0
      uu=0.0
      vv=0.0
      hn=0.0
      
      
      do m=1,12  !months

        do dm=1,ndaymonth(m)   !days of month

          day=day+1  !increment day of year

          if (day.eq.midday(m)) then
            if (dsun(day).eq.0.0) dsun(day)=10
            ni=dsun(day)/100.0   !proportion of bright sunshine today

c         CALCULATION OF NET DOWNWARD SHORT-WAVE RADIATION FLUX
c         (also known as insolation or incident solar radiation)
c         Refs: Prentice et al 1993, Monteith & Unsworth 1990,
c               Henderson-Sellers & Robinson 1986
c          (1) rs = (c + d*ni) * (1 - beta) * Qo * cos Z * k
c                (Eqn 7, Prentice et al 1993)
c          (2) Qo = Qoo * ( 1 + 2*0.01675 * cos ( 2*pi*i / 365) )
c                (Eqn 8, Prentice et al 1993; angle in radians)
c          (3) cos Z = sin(lat) * sin(delta) + cos(lat) * cos(delta) * cos h
c                (Eqn 9, Prentice et al 1993)
c          (4) delta = -23.4 * pi / 180 * cos ( 2*pi*(i+10) / 365 )
c                (Eqn 10, Prentice et al 1993, angle in radians)
c          (5) h = 2 * pi * t / 24 = pi * t / 12
c              where rs    = instantaneous net downward shortwave radiation
c                            flux (W/m2 = J/m2/s)
c                    c, d  = empirical constants (c+d = clear sky
c                            transmissivity)
c                    ni    = proportion of bright sunshine
c                    beta  = average 'global' value for shortwave albedo
c                           (not associated with any particular vegetation)
c                    i     = julian date, (1-365, 1=1 Jan)
c                    Qoo   = solar constant, 1360 W/m2
c                    Z     = solar zenith angle (angular distance between the
c                            sun's rays and the local vertical)
c                    k     = conversion factor from solar angular units to
c                            seconds, 12 / pi * 3600
c                    lat   = latitude (+=N, -=S, in radians)
c                    delta = solar declination (angle between the orbital
c                            plane and the Earth's equatorial plane) varying
c                            between +23.4 degrees in northern hemisphere
c                            midsummer and -23.4 degrees in N hemisphere
c                            midwinter
c                    h     = hour angle, the fraction of 2*pi (radians) which
c                            the earth has turned since the local solar noon
c                    t     = local time in hours from solar noon
c
c         From (1) and (3), shortwave radiation flux at any hour during the
c         day, any day of the year and any latitude given by
c          (6) rs = (c + d*ni) * (1 - beta) * Qo * ( sin(lat) * sin(delta) +
c                   cos(lat) * cos(delta) * cos h ) * k
c         Solar zenith angle equal to -pi/2 (radians) at sunrise and pi/2 at
c         sunset.  For Z=pi/2 or Z=-pi/2,
c          (7) cos Z = 0
c         From (3) and (7),
c          (8)  cos hh = - sin(lat) * sin(delta) / ( cos(lat) * cos(delta) )
c               where hh = half-day length in angular units
c         Define
c          (9) u = sin(lat) * sin(delta)
c         (10) v = cos(lat) * cos(delta)
c         Thus
c         (11) hh = acos (-u/v)
c         To obtain the net DAILY downward short-wave radiation flux, integrate
c         equation (6) from -hh to hh with respect to h,
c         (12) rs_day = 2 * (c + d*ni) * (1 - beta) * Qo *
c                       ( u*hh + v*sin(hh) ) * k
c         Define
c         (13) w = (c + d*ni) * (1 - beta) * Qo
c         From (12) & (13),
c         (14) rs_day = 2 * w * ( u*hh + v*sin(hh) ) * k
c#ifdef LPJ_STEP_2
c            IF(year==(spinup_years+1)) THEN
c#else
c            if (year.eq.1) then
c#endif

            IF (year==start_year+1) THEN

              qo(m)=qoo*(1.0+2.0*0.01675*cos(2.0*pi*real(day)/365.0))  !Eqn 2
              delta=-23.4*degtorad*cos(degtorad*360.0*                 !Eqn 4
     *               (real(day)+10.0)/365.0)
              if (m.eq.1) then
                sinelat=sin(lat*degtorad)
                cosinelat=cos(lat*degtorad)
              endif
              u(m)=sinelat*sin(delta) !Eqn 9
              v(m)=cosinelat*cos(delta) !Eqn 10

c           Calculate half-day length in angular units, hh
c           In Eqn (11), hh defined for u in range -v to v
c           For u >= v, hh = pi (12 hours, i.e. polar day)
c           For u <= -v, hh = 0 (i.e. polar night)
              if (u(m).ge.v(m)) then       !polar day
                hh(m)=pi
              elseif (u(m).le.-v(m)) then  !polar night
                hh(m)=0.0
              elseif (-u(m)/v(m).gt.1) then     !polar night
                hh(m)=0.0
              elseif (-u(m)/v(m).lt.-1) then    !polar day
                hh(m)=pi
              else                   !normal day and night
                hh(m)=acos(-u(m)/v(m))     !Eqn 11
              endif
              sinehh(m)=sin(hh(m))
            endif

            w=(c+d*ni)*(1.0-beta)*qo(m)                              !Eqn 13

c           Obtain daylength in hours from hh
            mdayl(m)=24.0*hh(m)/pi

c           Calculate total net downward shortwave radiation for this day
            rs_day=2.0*w*(u(m)*hh(m)+v(m)*sinehh(m))*k  !Eqn 14
            if (rs_day.eq.0.0) then
              rs_day=2.0*w*(u(m)*hh(m)+v(m)*sin(0.002))*k  !Eqn 14
            endif
c           Convert to PAR
            mpar_day(m)=rs_day*0.5

c         CALCULATION OF DAILY POTENTIAL EVAPOTRANSPIRATION
c         (PET, equilibrium evapotranspiration, or evaporative demand)
c         Refs: Jarvis & McNaughton 1986, Prentice et al 1993
c         (15) pet = ( s / (s + gamma) ) * rn / lambda
c                (Eqn 5, Prentice et al 1993)
c         (16) s = 2.503E+6 * exp ( 17.269 * temp / (237.3 + temp) ) /
c                  (237.3 + temp)**2
c                (Eqn 6, Prentice et al 1993)
c         (17) rn = rs - rl
c         (18) rl = ( b + (1-b) * ni ) * ( a - temp )
c                (Eqn 11, Prentice et al 1993)
c              where pet    = instantaneous evaporative demand (mm/s)
c                    gamma  = psychrometer constant, c. 65 Pa/K
c                    lambda = latent heat of vapourisation of water,
c                             c. 2.5E+6 J/kg
c                    temp   = temperature (deg C)
c                    rl     = net upward longwave radiation flux
c                             ('terrestrial radiation') (W/m2)
c                    a, b   = empirical constants
c         Note: gamma and lambda are weakly dependent on temperature. Simple
c               linear functions are used to obtain approximate values for a
c               given temperature
c         Define
c         (19) uu = w * u - rl
c         (20) vv = w * v
c         Total daily PET is instantaneous PET integrated over the period
c         during which rn >= 0. Limits for the integration (half-period
c         hn) are obtained by solving for
c         (21) rn=0
c         From (17) & (21),
c         (22) rs - rl = 0
c         From (6), (19), (20) and (22),
c         (23) uu + vv * cos hn = 0
c         From (23),
c         (24) hn = acos ( -uu/vv )
c         Integration of (15) w.r.t. h in the range -hn to hn leads to the
c         following formula for total daily PET (mm)
c         (25) pet=2*(s/(s+gamma)/lambda)*(uu*hn+vv*sin(hn))*k

            rl=(b+(1.0-b)*ni)*(a-dtemp(day))   !Eqn 18
c            rl=32.0*10e-5*(1.0+4.0*ni)*(100.0-dtemp(day))
            gamma=65.05+dtemp(day)*0.064
            lambda=2.495e6-dtemp(day)*2380.0
            s=2.503e6*exp(17.269*dtemp(day)/(237.3+dtemp(day)))/
     *        (237.3+dtemp(day))**2             !Eqn 16
            uu=w*u(m)-rl  !Eqn 19
            vv=w*v(m)     !Eqn 20

c         Calculate half-period with positive net radiation, hn
c         In Eqn (24), hn defined for uu in range -vv to vv
c         For uu >= vv, hn = pi (12 hours, i.e. polar day)
c         For uu <= -vv, hn = 0 (i.e. polar night)

            if (uu.ge.vv) then       !polar day
              hn=pi
            elseif (uu.le.-vv) then  !polar night
              hn=0.0
            else                     !normal day and night
              hn=acos(-uu/vv)        !Eqn 24
            endif

c Calculate total PET for this day
            mpet(m)=2.0*(s/(s+gamma)/lambda)*(uu*hn+vv*sin(hn))*k  !Eqn 25
c Doug 08/09: Calcualte equlibrium ET
            deet(day)=lambda*(rs_day-rl)*(s/(s + gamma))
          endif
        enddo  !ndaymonth loop

      enddo    !month loop

      return
      end

c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE SNOW
c     Adjust daily precipitation by snowmelt and accumulation in snowpack
c     Ref: Haxeltine & Prentice 1996

      subroutine snow(dtemp,dprec,snowpack,dmelt,msnowpack)
      implicit none

c     PARAMETERS:
      real km,tsnow
c     km changed from 0.7 to 3.0  7/00 for non-permafrost as in HBV
c     note below: new algorithm
       parameter (tsnow=0.0, km=3.0)
      integer ndaymonth(1:12)
      integer m
      data (ndaymonth(m),m=1,12)
     *  / 31,28,31,30,31,30,31,31,30,31,30,31 /

c     ARGUMENTS:
      real dtemp(1:365),dprec(1:365)
      real snowpack
      real dmelt(1:365)
      real msnowpack(1:12) ! monthly aver. snowpack depth in m

c     LOCAL VARIABLES:
      logical perm_attr
      integer d,dm
      real snowmelt,newsnow

      perm_attr=.false.
      d=0
      dm=0
      snowmelt=0.0
      newsnow=0.0

      do m=1,12  !months
        msnowpack(m)=0.0
        do dm=1,ndaymonth(m)   !days of month

c       Calculate snow melt and new snow for today
          d=d+1  !increment day of year
          if (dtemp(d).lt.tsnow) then
           newsnow=dprec(d)
           snowmelt=0.0
          else
           newsnow=0.0
c           snowmelt=km*(dtemp(d)-tsnow)
           snowmelt=(1.5+0.007*dprec(d))*(dtemp(d)-tsnow)
          endif

c       Set maximum snowmelt to size of snowpack
          if (snowmelt.gt.snowpack) snowmelt=snowpack

c       Update snowpack store
        snowpack=snowpack+newsnow-snowmelt

c       set maximum snow pack =10m
        snowpack=min(snowpack,10000.0)

c      density water=1000kg/m3, snow=400kg/m3; factor 2.5
          msnowpack(m)=msnowpack(m)+snowpack/ndaymonth(m)/1000.0*2.5

c       Calculate effective water supply (as daily values in mm)
          dprec(d)=dprec(d)-newsnow
          dmelt(d)=snowmelt
        enddo
      enddo

      return
      end


c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c SUBROUTINE PERMAFROST
      subroutine permafrost(soilcode,mtemp,mw1,mw2,mw1_t,mw2_t,
     * msnowpack,dtemp,dthaw_depth,maxthaw_depth,
     * perm_attr,year,maxthaw_old)

      implicit none

c     PARAMETERS

      real snow_heat_cap,snow_dens,adj_ampl_snow,moss_depth,pi,
     *thaw_moss_therdiff,froz_moss_therdiff,period,d1,d2


      parameter (snow_heat_cap=4000.0) ! Snow heat capacity (J/kg*C)
      parameter (snow_dens=300.0)      ! Snow density (kg/m3) (can be correcred by
                                       ! Yoshyda,1955 equation)
*      parameter (adj_ampl_snow=0.0)    ! Temperatute amplidute adjustment of snow
                                       ! (C)(to be corrected)
      parameter (moss_depth=0.1)       ! Depth of moss/litter (m)
      parameter (pi=3.1416)
      parameter (thaw_moss_therdiff=0.0000000556)
               !thermal diffusitivity of 10 cm thawed moss (m2/s)(Kudryavtsev,et al.1974)
      parameter (froz_moss_therdiff=0.00000139)
               !thermal diffusitivity of 10 cm frozen moss (m2/s)
      parameter (period=31536000.)     !year expressed in s
      parameter (d1=500.0,d2=1000.0)

c     ARGUMENTS:
      integer soilcode
      real mtemp(1:12)
      real mw1(1:12)                   !
      real mw2(1:12)                   !monthly fraction of soil moisture in bucket
      real mw1_t(1:12)                 !liquid+frozen water
      real mw2_t(1:12)                 !liquid+frozen water
      real msnowpack(1:12)             !monthly snowpack (m)
      real dtemp(1:365)
      real dthaw_depth(1:365)          !daily thaw depth (mm)
      real maxthaw_depth               !maximum thaw depth (mm)
c DM      real delta_thawday                ! difference in maximum thaw depth (mm/day))
      logical perm_attr                !if true permafrost
      integer year
      real maxthaw_old                   ! maximum thaw depth from the previous (mm) year

c     LOCAL VARIABLES:
      integer i,j,l,m
      real store(1:6,1:9)
      real air_temp,moss_temp,soil_temp,perfrost_temp

      !mean annual air, moss, soil surface and top of the permafrost temp. (C)
      real rel_moist                   !mean annual farction soil moisture
      real snow_depth                  !mean annual snow depth (m)
      real length_posit,length_frost   !lenght of positive and negative
                                       !temperatures (s)
      real ampl_air,ampl_moss,ampl_soil,ampl_perfrost
      !amplitude of air, moss and soil surface temperature (C)

      real adj_temp_snow,adj_ampl_moss,adj_moss_temp
      !adjustments of mean annual temperatures and amplitudes (C)

      real diff_abmoss_temp,diff_blmoss_temp
      !average differences in temperatures above and bellow moss

      real snow_therm_cond,thaw_therm_cond,froz_therm_cond
      ! snow, thawed and frozen soil thermal conductivities (W/m*C)
      
      real lat_heat
      ! latent heat of phase change (J/kg)

      real vol_heat_cap
      !  volumetric heat capacity of soil J/m3*C

      real depth_term     !active layer thickness caused only by heat transfer

      real temp_numerator,liambda                ! auxilary variables
      real mintemp,maxtemp,arcsine
      real alpha1,betha1,gamma1,ksi1,d,b,e

      adj_ampl_snow=0.0
      
      
      i=0
      j=0
      l=0
      m=0
      store(:,:)=0.0
      air_temp=0.0
      moss_temp=0.0
      soil_temp=0.0
      perfrost_temp=0.0
      rel_moist=0.0
      snow_depth=0.0
      length_posit=0.0
      length_frost=0.0
      ampl_air=0.0
      ampl_moss=0.0
      ampl_soil=0.0
      ampl_perfrost=0.0
      adj_temp_snow=0.0
      adj_ampl_moss=0.0
      adj_moss_temp=0.0
      diff_abmoss_temp=0.0
      diff_blmoss_temp=0.0
      snow_therm_cond=0.0
      thaw_therm_cond=0.0
      froz_therm_cond=0.0
      lat_heat=0.0
      vol_heat_cap=0.0
      depth_term=0.0
      temp_numerator=0.0
      liambda=0.0
      mintemp=0.0
      maxtemp=0.0
      arcsine=0.0
      alpha1=0.0
      betha1=0.0
      gamma1=0.0
      ksi1=0.0
      d=0.0
      b=0.0
      e=0.0


c      LPJ soil code  Description        Matter
c     --------------  --------------    ------------------

c            1        coarse             Sand
c            2        medium             Silt
c            3        fine, non-vertisol Clay
c            4        medium-coarse
c            5        fine-coarse
c            6        fine-medium
c            7        fine-medium-coarse
c            8        organic
c            9        fine, vertisol     Ice

c     -----------------------------------------------------------------
c     -----------------------------------------------------------------

c     Array 'store'
c     1  Density (kg/m3)
c     2  Heat capacity(J/kg*C)
c     3  Slope in equation for thermal conductivity (W/m*C) for
c        thawed soil
c     4  Interception in equation for thermal conductivity (W/m*C) for
c        thawed soil
c     5  Slope in equation for thermal conductivity (W/m*C) for
c        frozen soil
c     6  Interception in equation for thermal conductivity (W/m*C) for
c        frozen soil

c        Basic data from Pavlov (1976)

      data ((store(i,j),i=1,6),j=1,9) /

c    ----------------------------------------------------
c           1       2       3       4       5        6
c    ----------------------------------------------------
     *  1300.0,   690.0,   1.10,  1.05,    1.40,     1.25,
     *  1400.0,   730.0,   0.85,  1.05,    1.15,     1.25,
     *  1500.0,   900.0,   0.80,  0.90,    0.85,     1.15,
     *  1350.0,   710.0,   0.975, 1.05,    1.275,    1.25,
     *  1400.0,   795.0,   0.95,  0.975,   1.125,    1.20,
     *  1450.0,   815.0,   0.825, 0.975,   1.00,     1.12,
     *  1400.0,   773.4,   0.917, 1.00,    1.13,     1.22,
     *   150.0,   200.0,   0.25,  0.35,    1.85,     0.80,
     *     0.0,     0.0,   0.0,   0.0,     0.0,      0.0/ ! Ice is ign.

C IF NOT RUNNING PERMAFROST MODEL THEN UNDELETE THE FOLLOWING LINES
      perm_attr=.false.
      maxthaw_depth=d2+d1
      do i=1,365
         dthaw_depth(i)=d2+d1
      enddo
      return

c     ----------------------------------------------------------------
c     Calculate annual mean air temperature

      maxthaw_old=0.0
      air_temp=0.0
      do i=1,12
         air_temp=air_temp+mtemp(i)/12.0
      enddo

      if ((air_temp.le.0.0)) then
c     10.11.12 the condition air_temp.gt.-20 is abandoned
c     continuous or discontinuous permafrost region

        perm_attr=.true.

       if (soilcode.lt.1.or.soilcode.gt.9) stop 'soil type incorrect'

c    Calculate annual mean relative soil moisture and snowpack
        rel_moist=0.0
        do i=1,12
c            rel_moist=rel_moist+(mw2(i)+mw1(i))/24.0
            rel_moist=rel_moist+(mw2_t(i)+mw1_t(i))/24.0
        enddo
        rel_moist=min(rel_moist,1.0)
        snow_depth=0.0
        l=0
        do i=1,12

        if (mtemp(i).le.0.0) then
           l=l+1
           snow_depth=snow_depth+msnowpack(i)
         endif
        enddo
        snow_depth=snow_depth/(real(l))
c     calculate thawed soil thermal conductivity using linear regression
        thaw_therm_cond= store(3,soilcode)*rel_moist +
     *                   store(4,soilcode)
c     calculate frozen soil thermal conductivity using linear regression
        froz_therm_cond= store(5,soilcode)*rel_moist +
     *                   store(6,soilcode)

c    calculate length of frost period and non-frost period
        length_frost=0.0
        do i=1,365
            if (dtemp(i).lt.0.0) length_frost=length_frost+1.0
        enddo
        length_posit=365.0 - length_frost

c     days->s
        length_frost=length_frost*86400.0
        length_posit=length_posit*86400.0

c    calculate air temperature amplitude
*        ampl_air=(abs(mtemp(1))+abs(mtemp(7)))/2.0
       mintemp=100.0
       maxtemp=-100.0
       do i=1,12
         if (mtemp(i).lt.mintemp) mintemp=mtemp(i)
         if (mtemp(i).gt.maxtemp) maxtemp=mtemp(i)
       enddo
c      ampl_air=(abs(mintemp)+abs(maxtemp))/2.0
       ampl_air=(abs(mintemp+maxtemp))/2.0       !!! Yan 22/11/07

c   define snow thermal conductivity (Abels,1892) using snow density
        snow_therm_cond=0.0000029*snow_dens**2.0

c   All mean temerature or amplitude adjustment (damping) calculations
c   based on stationary Stefan's solutions for heat-transfer equations

c   calculate snow temperature adjustment

        adj_temp_snow=ampl_air/2.0*(1.0-exp(((-snow_depth)*
     *                  ((pi*snow_dens*snow_heat_cap)/
     *                   (snow_therm_cond*period))**0.5)))

c   calculate mean understorey (moss) temperature and amplitude
        adj_ampl_snow=adj_temp_snow
        ampl_moss=ampl_air-adj_ampl_snow
        moss_temp=air_temp+adj_temp_snow

c   calculate differences between the average temperatures below and above moss
c   during cold and warm periods

        diff_abmoss_temp=(ampl_moss-moss_temp)*(1.0-exp(((-moss_depth)*
     *            ((pi)/(froz_moss_therdiff*2.0*length_frost))**0.5)))
        diff_blmoss_temp=(ampl_moss+moss_temp)*(1.0-exp(((-moss_depth)*
     *            ((pi)/(thaw_moss_therdiff*2.0*length_posit))**0.5)))

c   calculate amplitude and temperature adjustment due to moss

       adj_ampl_moss=(diff_blmoss_temp*length_posit+
     *                 diff_abmoss_temp*length_frost)/period
       adj_moss_temp=(diff_abmoss_temp*length_frost-
     *                 diff_blmoss_temp*length_posit)/period*
     *                 (2/pi)

c  calculate mean annual surface temperature and amplitude
       soil_temp=air_temp+adj_moss_temp+adj_temp_snow
       ampl_soil=ampl_air-adj_ampl_moss-adj_ampl_snow

c  calculate mean annual temperature at the top of permafrost (Garagulya,1990)
       arcsine=asin(soil_temp/ampl_soil)
       temp_numerator=0.5*soil_temp*(thaw_therm_cond+froz_therm_cond)+
     *                 ampl_soil*(thaw_therm_cond-froz_therm_cond)/pi*
     *                ((soil_temp/ampl_soil)*arcsine+
     *                   (1-(soil_temp)**2.0/(ampl_soil)**2.0)**0.5)

       if (temp_numerator.le.0.0) then
           perfrost_temp=temp_numerator/froz_therm_cond
           liambda=froz_therm_cond
       else
           perfrost_temp=temp_numerator/thaw_therm_cond
           liambda=thaw_therm_cond
       endif

c   calculate latent heat of phase change (J/kg)
      lat_heat=335200000.0*rel_moist

c   calculate volumetric heat capacity of soil J/m3?°C
       if (perfrost_temp.le.0.0) then
          vol_heat_cap=store(2,soilcode)*store(1,soilcode)+
     *                 2025000.0*rel_moist
        else
          vol_heat_cap=store(2,soilcode)*store(1,soilcode)+
     *                 4190000.0*rel_moist
       endif

c   calculate amplitude of temperature at the top of permafrost
      ampl_perfrost=((ampl_soil-abs(perfrost_temp))/
     *              log((ampl_soil+lat_heat/(2.0*vol_heat_cap))/
     *              (abs(perfrost_temp)+lat_heat/(2.0*vol_heat_cap))))
     *             -lat_heat/(2.0*vol_heat_cap)

c calculate active layer thickness caused by thermal transfer
       depth_term=(2.0*(ampl_soil-abs(perfrost_temp))*
     *              ((vol_heat_cap*period*liambda)/(pi))**0.5)/
     *              (2.0*vol_heat_cap*ampl_perfrost+lat_heat)

c   calculate auxilary variables as in the paper
       alpha1=2.0*vol_heat_cap*ampl_perfrost+lat_heat
       betha1=2.0*(ampl_soil-abs(perfrost_temp))*
     *              ((vol_heat_cap*period*liambda)/(pi))**0.5
       gamma1=2.0*vol_heat_cap*ampl_perfrost*depth_term
       ksi1=((period*liambda)/(pi*vol_heat_cap))**0.5
       d=alpha1*lat_heat
       b=gamma1*alpha1+(alpha1**2.0)*ksi1-
     *   betha1*lat_heat-(lat_heat**2.0)*ksi1
       e=gamma1*betha1+alpha1*betha1*ksi1+gamma1*ksi1*lat_heat
c   calculate maximum thaw depth as a positive solution of quadratic equation
      if (year.ne.1) then
         maxthaw_old=maxthaw_depth

c   depth of active layer in mm
         maxthaw_depth=((-1.0*b+(b*b+4.0*d*e)**0.5)/(2.0*d))*1000.0
      else
         maxthaw_depth=(-1.0*b+(b*b+4.0*d*e)**0.5)/(2.0*d)*1000.0
         maxthaw_old=maxthaw_depth
      endif
c DM      delta_thawday=(maxthaw_depth-maxthaw_old)/(length_posit/86400.0)

c this variable is used to estimate additional water lost or gained form
c maximum thaw differences,we assume this water to be equidistributed
c for thaw period
       do i=1,365
*           dthaw_depth(i)=maxthaw_depth
        dthaw_depth(i)=0.1
       enddo

c  find start day of spring
       l=0
       do i=1,183
       if (dtemp(i).lt.0.0) l=l+1
       enddo

c   calculate thaw depth profile
       do i=1,365
       if (dtemp(i).ge.0.0) dthaw_depth(i)=maxthaw_depth*
     *             (sin(real(i-l)/(length_posit/(86400.0))*pi/2.0))
        dthaw_depth(i)=max(0.1,dthaw_depth(i))
       enddo

c5mai depth of freezing
       m = int(length_posit/(86400.0))
       m = l+m

c5mai Now m is a start of fall
       do i=m,365
        dthaw_depth(i)=maxthaw_depth-(maxthaw_depth*
     *             (sin(real(i-m)/(length_frost/(86400.0))*pi/2.0)))
cThis is to fix a metaphysical occurance where the maxthaw_depth*sin()
cis ecaxtly equal to d1 for some reason ...
        if(abs(maxthaw_depth-dthaw_depth(i)-d1).lt.0.0001) then
           dthaw_depth(i)=dthaw_depth(i)-0.0001
        endif
       enddo

       do i=1,l
        dthaw_depth(i)=maxthaw_depth-(maxthaw_depth*
     *        (sin(real(i+364-m)/(length_frost/(86400.0))*pi/2.0)))
cThis is to fix a metaphysical occurance where the maxthaw_depth*sin()
cis ecaxtly equal to d1 for some reason ...
        if(abs(maxthaw_depth-dthaw_depth(i)-d1).lt.0.0001) then
           dthaw_depth(i)=dthaw_depth(i)-0.0001
        endif
       enddo

c4mai temporal fix to get right profile of thaw
        if (abs(soil_temp).gt.abs(ampl_soil)) then
            maxthaw_depth=0.01
            do i=1,365
              dthaw_depth(i)=0.01
            enddo
        endif

      else
c     no permafrost
        perm_attr=.false.
        maxthaw_depth=d2+d1

c    can be used in future for correct thaw modelling in non-permafrost region
*  out  291199
        do i=1,365
           dthaw_depth(i)=d2+d1
        enddo
      endif

      maxthaw_depth=max(0.1,maxthaw_depth)

      return
      end


c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE SUMMERPHENOLOGY
c     Temperature-based phenology for summergreen PFTs

      subroutine summerphenology(pftpar,mtemp,dtemp,gdd,dphen_t,
     *  summergreen,tree)

      implicit none

c     PARAMETERS
      integer npft,npftpar
        parameter (npft=13,npftpar=59)

c     ARGUMENTS
      real pftpar(1:npft,1:npftpar)
      real mtemp(1:12)
      real dtemp(1:365)
      real gdd(1:npft)
      real dphen_t(1:365,1:npft)
      logical summergreen(1:npft),tree(1:npft)

c     LOCAL VARIABLES
      integer ndaymonth(1:12),month
      data (ndaymonth(month),month=1,12)
     *  / 31,28,31,30,31,30,31,31,30,31,30,31 /
      integer midday(1:12)
      data (midday(month),month=1,12)
     *  / 16,44,75,105,136,166,197,228,258,289,319,350 /

      integer warmest,midsummer,firstday,pft,day
      real ramp
      real leafon
      integer cm1,cm2,cm3,coldest
      integer d
      real aphen
      real gddbase          !GDD base for greenup, pft parameter

      warmest=0
      midsummer=0
      firstday=0
      pft=0
      day=0
      ramp=0.0
      leafon=0.0
      cm1=0
      cm2=0
      cm3=0
      coldest=0
      d=0
      aphen=0.0
      gddbase=0.0
      
      dphen_t(:,:)=1.0

c     First find warmest month

      warmest=1
      do month=1,12
        if (mtemp(month).gt.mtemp(warmest)) warmest=month
      enddo

c     find the coldest month
      coldest=1
      do month=1,12
        if (mtemp(month).lt.mtemp(coldest)) coldest=month
      enddo

      midsummer=-15
      do month=1,warmest
        midsummer=midsummer+ndaymonth(month)
      enddo

      do pft=1,npft

c        Now midsummer is middle day (roughly) of warmest month
c        Find day of leaf abscission at end of summer

         gddbase=pftpar(pft,33)
         firstday=midsummer+1
         do while (dtemp(firstday).ge.gddbase.and.
     *        firstday.ne.midsummer)
            firstday=firstday+1
            if (firstday.gt.365) firstday=1
         enddo



        if (summergreen(pft)) then  !summergreen taxa
          ramp=pftpar(pft,19)   !number of GDDs to attain full leaf cover

          if (firstday.eq.midsummer) then  !no leaf abscission

            do day=1,365
              dphen_t(day,pft)=1.0
            enddo

          else

            gdd(pft)=0.0     !accumulated growing degree days
            leafon=0.0  !proportional leaf-on today

            day=firstday+1
            do while (day.ne.firstday)
              if (dtemp(day).gt.gddbase) then  !growing day
                gdd(pft)=gdd(pft)+dtemp(day)-gddbase
                if (ramp.gt.0.0) then
                  leafon=min(gdd(pft)/ramp,1.0)
                else
                  leafon=1.0
                endif
              endif
              dphen_t(day,pft)=leafon
              day=day+1
              if (day.gt.365) day=1
            enddo

          endif

c         constrain woody deciduous phenology to max= 9 months
          if (tree(pft)) then
             aphen=0.
             do day=1,365
                aphen=aphen+dphen_t(day,pft)  ! calculate the total number of days with foliage
             enddo
             if (aphen.gt.210) then   !limit summergreen to 210 days (approx 5 months)
               do d=midday(coldest),midday(coldest)+75   ! set 75 days after coldest with no foliage
                  if (d.le.365) then
                     day=d
                  else
                     day=d-365
                  endif
                  dphen_t(day,pft)=0.0
               enddo
               do d=midday(coldest)-75,midday(coldest)
                  if(d.ge.1) then
                     day=d
                  else
                     day=365+d
                  endif
                  dphen_t(day,pft)=0.0
               enddo
             endif
           endif

        else  !non-summergreen taxa

          do day=1,365
            dphen_t(day,pft)=1.0
          enddo

        endif

      enddo   !pft

      return
      end



c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE CLIMATE20
c     Calculation of 20-year average coldest-month temperatures and
c     growing degree days

      subroutine climate20(mtemp,dtemp,gdd,mtemp_min_buf,mtemp_max_buf,
     *  gdd_buf,year,mtemp_min20,mtemp_max20,gdd20,mtemp_max,pftpar)

      implicit none

c     PARAMETERS
      integer npft,npftpar,climbuf
        parameter (npft=13,npftpar=59,climbuf=20)
                                !NB must be same values as in main program

c     ARGUMENTS
      real mtemp(1:12),dtemp(1:365)
      real gdd(1:npft)
      real mtemp_min_buf(1:climbuf)
      real mtemp_max_buf(1:climbuf)
      real gdd_buf(1:npft,1:climbuf)
      integer year
      real mtemp_min20, mtemp_max20,gdd20(1:npft)
      real mtemp_max
      real pftpar(1:npft,1:npftpar)

c     LOCAL VARIABLES
      integer n,m,bufsize,d,pft
      real mtemp_min
      real gddbase                    !GDD base for greenup, pft parameter


      n=0
      m=0
      bufsize=0
      d=0
      pft=0
      mtemp_min=0.0
      gddbase=0.0


c     Find mean temperature of the coldest month this year

      mtemp_min=mtemp(1)
      do m=2,12
        if (mtemp(m).lt.mtemp_min) mtemp_min=mtemp(m)
      enddo

c     Find mean temperature of the warmest month this year
      mtemp_max=mtemp(1)
      do m=2,12
        if (mtemp(m).gt.mtemp_max) mtemp_max=mtemp(m)
      enddo

c     Find growing degree days

      do pft=1,npft
         gddbase=pftpar(pft,33)
         gdd(pft)=0.0
         do d=1,365
            gdd(pft)=gdd(pft)+max(dtemp(d)-gddbase,0.0)
         enddo
      enddo

      if (year.gt.climbuf) then

c       Shift stored yearly values up, adding new values for this year

        do n=2,climbuf
          mtemp_min_buf(n-1)=mtemp_min_buf(n)
          mtemp_max_buf(n-1)=mtemp_max_buf(n)
          do pft=1,npft
             gdd_buf(pft,n-1)=gdd_buf(pft,n)
          enddo
        enddo

        mtemp_min_buf(climbuf)=mtemp_min
        mtemp_max_buf(climbuf)=mtemp_max
        do pft=1,npft
           gdd_buf(pft,climbuf)=gdd(pft)
        enddo

      else

c       Average over number of simulation years so far if less than
c       climbuf (20)

        mtemp_min_buf(year)=mtemp_min
        mtemp_max_buf(year)=mtemp_max
        do pft=1,npft
           gdd_buf(pft,year)=gdd(pft)
        enddo

      endif

c     Calculate averages

      mtemp_min20=0.0
      mtemp_max20=0.0
      do pft=1,npft
         gdd20(pft)=0.0
      enddo

      bufsize=min(year,climbuf)

      do n=1,bufsize
        mtemp_min20=mtemp_min20+mtemp_min_buf(n)/real(bufsize)
        mtemp_max20=mtemp_max20+mtemp_max_buf(n)/real(bufsize)
        do pft=1,npft
           gdd20(pft)=gdd20(pft)+gdd_buf(pft,n)/real(bufsize)
        enddo
      enddo

      return
      end



c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE BIOCLIM
c     Apply bioclimatic limits on PFT survival and establishment

      subroutine bioclim(pftpar,mtemp_min20,mtemp_max20,gdd,mtemp_max,
     *         survive,estab,year,present)

      implicit none

c     PARAMETERS
      integer npft,npftpar
      parameter (npft=13,npftpar=59)

c     ARGUMENTS
      real pftpar(1:npft,1:npftpar)
      real mtemp_min20
      real mtemp_max20
      real gdd(1:npft)
      real mtemp_max
      logical survive(1:npft),estab(1:npft)
      integer year
      logical present(1:npft)

c     LOCAL VARIABLES
      real twmax
      real tcmin,tcmax,gddmin
      integer pft
      real min_temprange              !minimum 20-year average warmest-coldest
                                      !month temperature range

      twmax=0.0
      tcmin=0.0
      tcmax=0.0
      gddmin=0.0
      pft=0
      min_temprange=0.0


      do pft=1,npft

c       Limits based on 20-year running averages of coldest-month mean
c       temperature and growing degree days (5 degree base, except larches
c       (2 degrees)), and minimal temperature range required (larches).
c       For SURVIVAL, coldest month temperature and GDD should be
c       at least as high as PFT-specific limits.
c       For REGENERATION, PFT must be able to survive AND coldest month
c       temperature should be no higher than a PFT-specific limit.

        tcmin=pftpar(pft,28)   !PFT-specific minimum coldest-month temperature
        tcmax=pftpar(pft,29)   !PFT-specific maximum coldest-month temperature
        gddmin=pftpar(pft,30)  !PFT-specific minimum GDD
        twmax=pftpar(pft,31)  !PFT-specific upper limit of warmest-month temperature
        min_temprange=pftpar(pft,34)!PFT-specific lower limit of
                                    !20-year temperature range

        if (present(pft)) tcmin=tcmin*1.1
        if (mtemp_min20.ge.tcmin) then
          survive(pft)=.true.
c          if (mtemp_min20.le.tcmax.and.gdd(pft).ge.gddmin) then

          if (mtemp_min20.le.tcmax.and.gdd(pft).ge.gddmin.and.
     *              mtemp_max.le.twmax) then
            estab(pft)=.true.
          else
            estab(pft)=.false.
          endif
        else
          survive(pft)=.false.
        endif

        if ((mtemp_max20-mtemp_min20).lt.min_temprange) then
           survive(pft)=.false.
        endif
       enddo

      return
      end


c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE GPP
c     Calculation of GPP, explicitly linking photosynthesis and water
c     balance through canopy conductance feedback
c     timestep: 

      subroutine gpp(present,co2,soilpar,pftpar,lai_ind,fpc_grid,
     *  mdayl,mtemp,mpar_day,dphen_t,w,w_t,dpet,dprec,dmelt,sla,agpp,
     *  alresp,arunoff_surf,arunoff_drain,arunoff,mrunoff,
     *  dwscal365,dphen_w,dphen,dphen_change, !Doug 05/09: dphen_change
     *  wscal,mgpp,mlresp,mw1,dw1,aaet,
     *  leafondays,leafoffdays,leafon,tree,raingreen,mpar,mapar,mphen,
     *  year,maet,littercpool_ag,mw2,maxthaw_depth,dthaw_depth,
     *  uw1,uw2,fw,
     *  litter_ag_leaf,litter_ag_wood,   !Doug 11/12
     *  mcica,mauw1,mauw2,lat,lon,w_ep,
     *  d_evap,maep,aaep,mintc,aintc,needle,boreal,dayl,mpet_grid,
     *  apet_grid,mpet2,apet,mpet,mmelt,mw1_t,mw2_t,meangc,mgp,
     *  deltaa,deltaa_fpc,deet,alpha_ws)    !Doug 07/09, added alpha_ws

      implicit none

c     PARAMETERS
      integer npft,npftpar,nsoilpar
        parameter (npft=13,npftpar=59,nsoilpar=7)
      integer nco2
        parameter (nco2=3)           !number of C variables: C,13C,14C

      real epsilon
        parameter (epsilon=0.1)  !min precision of solution in bisection method
c        parameter (lambdam=0.8)  !optimal Ci/Ca ratio
c commented out, made pft depandant for carbon isotopes, see below

c     ARGUMENTS:
      logical present(1:npft)
      real co2(1:nco2)
      real soilpar(1:7),pftpar(1:npft,1:npftpar)
      real lai_ind(1:npft),fpc_grid(1:npft)
      real mdayl(1:12),mtemp(1:12),mpar_day(1:12)
      real dphen_t(1:365,1:npft)
      real w(1:2),w_t(1:2)
      real dpet(1:365),dprec(1:365),dmelt(1:365)
      real sla(1:npft)
      real agpp(1:npft,1:nco2),alresp(1:npft,1:nco2)
      real arunoff_surf,arunoff_drain
      real arunoff
      real mrunoff(1:12)
      real dwscal365(1:npft)
      real dphen_w(1:365,1:npft)
      real dphen(1:365,1:npft)
      REAL dphen_change(1:12,1:npft) ! Doug 05/09: records timeing of chang in dphen for each month and pft.
      REAL day_dphen
      real wscal(1:npft)
      real mgpp(1:12,1:npft,1:nco2)
      real mlresp(1:12,1:npft,1:nco2)
      real mw1(1:12),dw1(1:365)
      real aaet
      integer leafondays(1:npft),leafoffdays(1:npft)
      logical leafon(1:npft),tree(1:npft),raingreen(1:npft)
c     mpar=monthly grid par (mpar_day(m)*ndays(month))
c     mapar=monthly grid apar (sum of mpar_day(m)*fpar(m)*alphaa over pfts)
c     mphen=monthly average pheno-state (sum_d(m) of dphen(d,pft)/ndays(m))
      real mpar(1:12)
      real mapar(1:12)
      real mphen(1:12,1:npft)
      integer year                    !simulation year
      real maet(1:12)
      real littercpool_ag    !no need for tripling?
      real mw2(1:12)
      real maxthaw_depth,dthaw_depth(1:365)
      real uw1,uw2,fw(1:2)
c DM      real delta_thawday
      real litter_ag_leaf(1:npft,1:nco2)          !gridcell above-ground litter (gC/m2)
      real litter_ag_wood(1:npft,1:nco2)		  ! Doug 11/12
      real mcica(1:12,1:npft)
      real mauw1(1:12)
      real mauw2(1:12)
      REAL lat,lon
      real w_ep
      real d_evap
      real maep(1:12),aaep
      real mintc(1:12)
      real aintc
      logical needle(1:npft), boreal(1:npft)
      real dayl(1:365)
      real mpet_grid(1:12)
      real apet_grid
      real mpet2(1:12),apet
      real mpet(1:12)
      real mmelt(1:12)
      real mw1_t(1:12)
      real mw2_t(1:12)
      real meangc(1:12,1:npft),mgp(1:12,1:npft)
	  real area


c     LOCAL:
      integer leafondays_an(1:npft)
      real muw1(1:12)
      real muw2(1:12)
      integer nc                    ! for nco2
      integer pft,aleafdays(1:npft),m,b,d,dm,i
      real ca,k(1:2),fwhc(1:2)
      real awscal(1:npft),dwscal(1:npft)
      real minwscal(1:npft),gminp(1:npft)
      real rootprop(1:2,1:npft)
      real emax(1:npft)                                       ! Sibyll
      real intc(1:npft)
      real meanfpc(1:12,1:npft)
      real nmax(1:npft),mrunoff_surf(1:12),mrunoff_drain(1:12)
      real meangmin(1:12,1:npft),tsecs(1:12),fpar,gp(1:12)
      real rd,agd,dval(1:365),dgp(1:365,1:npft)
      real drunoff_surf,drunoff_drain,dgc(1:365,1:npft)
      real gpd,x1,x2,xmid,fmid,rtbis,dx,adt1,adt2,adtmm
      logical c4(1:npft)
      real daet,daep,pet_grid,intercp_tot
      REAL deet(1:365)    !Doug 08/09: daily equilibrium evapotranspiration
      REAL alpha_ws       !Doug 07/09: cliamtic moisture availability index (mesureing water stress)
      real longevity(1:npft)
      real inhibx1(1:npft),inhibx2(1:npft)
      real inhibx3(1:npft),inhibx4(1:npft)
      real dw2(1:365)
      real fpc_tree
     
       real deltaa(1:npft),deltaa_fpc
 
      integer ndaymonth(1:12)    !number of days in each month
      data (ndaymonth(m),m=1,12)
     *  / 31,28,31,30,31,30,31,31,30,31,30,31 /
      real lambdam(1:npft)
      data (lambdam(i),i=1,npft)
     *     /0.9,0.9,0.9,0.8,0.8,0.8,0.9,0.65,0.4,0.9,0.9,0.8,0.8/

c     Initialisations of locals
      leafondays_an=0
      muw1(:)=0.0
      muw2(:)=0.0
      nc=0
      pft=0
      aleafdays(:)=0
      m=0
      b=0
      d=0
      dm=0
      i=0
      ca=co2(1)*1.0e-6  !from ppmv to mole fraction
      k(1)=soilpar(1)
      k(2)=soilpar(2)
      fwhc(1)=soilpar(3)
      fwhc(2)=soilpar(4)
      awscal(:)=0.0
      dwscal(:)=0.0
      minwscal(:)=0.0
      gminp(:)=0.0
      rootprop(:,:)=0.0
      emax(:)=0.0
      intc(:)=0.0
      meanfpc(:,:)=0.0
      nmax(:)=0.0
      mrunoff_surf(:)=0.0
      mrunoff_drain(:)=0.0
      meangmin(:,:)=0.0
      tsecs(:)=0.0
      fpar=0.0
      gp(:)=0.0
      rd=0.0
      agd=0.0
      dval(:)=0.0
      dgp(:,:)=0.0
      drunoff_surf=0.0
      drunoff_drain=0.0
      dgc(:,:)=0.0
      gpd=0.0
      x1=0.0
      x2=0.0
      xmid=0.0
      fmid=0.0
      rtbis=0.0
      dx=0.0
      adt1=0.0
      adt2=0.0
      adtmm=0.0
      c4(:)=.false.
      daet=0.0
      daep=0.0
      pet_grid=0.0
      intercp_tot=0.0
      longevity(:)=0.0
      inhibx1(:)=0.0
      inhibx2(:)=0.0
      inhibx3(:)=0.0
      inhibx4(:)=0.0
      dw2(:)=0.0
      fpc_tree=0.0


c     Other Initialisations
      do m=1,12
        mrunoff(m)=0.0
        maet(m)=0.0
        mintc(m)=0.0
        mmelt(m)=0.0
        maep(m)=0.0
        mpet_grid(m)=0.0
        mpet2(m)=0.0
      enddo

      arunoff_surf=0.0
      arunoff_drain=0.0
      aintc=0.0
      aaet=0.0
      aaep=0.0
      apet_grid=0.0
      apet=0.0
      dphen_change(:,:)=0.0 !Doug 05/09
      
c DM  mphen is set to 0 only when present(pft) is true
c DM  Maybe it's OK for the results, but it's not enough if we want
c DM  the exact same results when using the standard version and
c DM  the two steps version.
      mphen(:,:)=0.0
      mcica(:,:)=0.0
      mgp(:,:)=0.0

      alpha_ws=0.0    !Doug 07/09: initialise water stress meaure, alpha

      littercpool_ag=0.0
      do pft=1,npft
         littercpool_ag=littercpool_ag+  !Doug 11/12
     *     litter_ag_leaf(pft,1)+litter_ag_wood(pft,1)
      enddo

      do pft=1,npft
        if (tree(pft).and.present(pft)) fpc_tree=fpc_tree+fpc_grid(pft)
      enddo


      do pft=1,npft
        longevity(pft)=pftpar(pft,10)
c
c      define pft inhibition function parameters
c
        inhibx1(pft)=pftpar(pft,24)
        inhibx2(pft)=pftpar(pft,25)
        inhibx3(pft)=pftpar(pft,26)
        inhibx4(pft)=pftpar(pft,27)

        if (present(pft)) then
          minwscal(pft)=pftpar(pft,3)
          if (pftpar(pft,2).eq.1.0) then
            c4(pft)=.true.
          else
            c4(pft)=.false.
          endif

c         gminp = PFT-specific min canopy conductance scaled by fpc
c         assuming full leaf cover

          gminp(pft)=pftpar(pft,4)*fpc_grid(pft)
          rootprop(1,pft)=pftpar(pft,1)
          rootprop(2,pft)=1.0-pftpar(pft,1)
          emax(pft)=pftpar(pft,35)
          intc(pft)=pftpar(pft,36)
          do nc=1,nco2
            agpp(pft,nc)=0.0
            alresp(pft,nc)=0.0
          enddo
          nmax(pft)=pftpar(pft,7)

        endif
      enddo  !pft

      d=0  !day of year

      do pft=1,npft

        if (present(pft)) then

c         FIND THE POTENTIAL CANOPY CONDUCTANCE REALISABLE UNDER
c         NON-WATER-STRESSED CONDITIONS

          do m=1,12	!month

c           Initialisations

            meanfpc(m,pft)=0.0
            meangc(m,pft)=0.0
            meangmin(m,pft)=0.0
c DM  See above
c            mphen(m,pft)=0.0
c            mcica(m,pft)=0.0
c            mgp(m,pft)=0.0
 
            tsecs(m)=3600.0*mdayl(m)  !number of daylight seconds/day

c           Calculate non-water-stressed net daytime photosynthesis
c           assuming full leaf cover

            fpar = fpc_grid(pft)

            call photosynthesis(ca,mtemp(m),fpar,mpar_day(m),mdayl(m),
     *        c4(pft),sla(pft),nmax(pft),lambdam(pft),rd,agd,adtmm,
     *        inhibx1(pft),inhibx2(pft),inhibx3(pft),inhibx4(pft),pft)

            if (tsecs(m).gt.0.0) then
c             Calculate non-water-stressed canopy conductance
c             (gp) mm/sec basis averaged over entire grid cell
c             Eqn 21 Haxeltine & Prentice 1996

              gp(m)=(((1.6*adtmm)/(ca*(1.0-lambdam(pft))))/tsecs(m))+
     *          gminp(pft)

            else
              gp(m)=0.0
            endif

            mgp(m,pft)=gp(m)
          enddo  !month

c         Linearly interpolate mid-monthly gp to daily values

          call daily(gp,dval)
          do d=1,365	!day
            dgp(d,pft)=dval(d)
          enddo

        endif

      enddo !pft


      d=0 !day of year
c     mpar=monthly grid par (mpar_day(m)*ndays(month))
c     mapar=monthly grid apar (sum of mpar_day(m)*fpar(m)*alphaa over pfts)
      do m=1,12
         mpar(m)=0.0
         mapar(m)=0.0
      enddo

      do m=1,12
c       CALCULATE ACTUAL EVAPOTRANSPIRATION AND SOIL WATER BALANCE TODAY
C       FOR EACH DAY THIS MONTH

        mw1(m)=0.0
        mw2(m)=0.0
        mw1_t(m)=0.0
        mw2_t(m)=0.0
        muw1(m)=0.0
        muw2(m)=0.0
        mauw1(m)=0.0
        mauw2(m)=0.0

        do dm=1,ndaymonth(m)
          d=d+1

          do pft=1,npft
            if (present(pft)) then

c             Use yesterday's potential water scalar to determine
c             today's drought phenology

              if (d.eq.1) then
                dwscal(pft)=dwscal365(pft)
              endif

c             Drought phenology and net phenology for today.
c             Drought deciduous PFTs shed their leaves when their
c             water scalar falls below their PFT specific minimum
c             value (minwscal). Leaves are replaced immediately (i.e. daily)
c             once the minimum water scalar is exceeded.

              if (dwscal(pft).gt.minwscal(pft).and.leafon(pft)) then
                dphen_w(d,pft)=1.
                dphen(d,pft)=dphen_t(d,pft)
                leafondays(pft)=leafondays(pft)+1
              else
                dphen_w(d,pft)=0.0
                dphen(d,pft)=0.0
              endif
             
c
c  stop deciduous vegetation behaving like evergreen when climate
c  permits
c
              if (raingreen(pft).and.tree(pft)) then
              if (real(leafondays(pft)).ge.(365.0*longevity(pft))) then
                 leafoffdays(pft)=leafoffdays(pft)+1
                 leafon(pft)=.false.
              if (real(leafoffdays(pft)).ge.(365.0*longevity(pft))) then
                    leafoffdays(pft)=0
                    leafondays(pft)=0
                    leafon(pft)=.true.
                 endif
              endif
              endif

            endif


c    Doug 05/09: Finds day of (1st) change in dephen in the month, and records it for output
            IF (d>1.AND.dphen(d,pft)/=dphen(d-1,pft)) THEN

              IF(dphen_change(m,pft)==0) THEN !Doug 05/09: see if change in dphen has already occured this month for this pft
                IF(dphen(d,pft)-dphen(d-1,pft)>0) THEN
                  dphen_change(m,pft)=d  !Doug 05/09: if change is posative, day number is recorded as posative
                ELSE
                  dphen_change(m,pft)=-d !Doug 05/09: iof negative, recorded as negative
                END IF
              ELSE                            !Doug 05/09: write to error file

C                WRITE(10,*), "+-+-+-+-+-+-+-+-+-+-+-"
C                WRITE(10,*), "2nd change in dphen in a month for:"
C                WRITE(10,*), "lat:", lat, "lon", lon, "pft:",pft
C                WRITE(10,*), "year:",year,"month:",m,"day:",d
C                WRITE(10,*), "**********************"

              END IF !change in dphen already occured
            END IF !change in pdhen

          enddo  !pft

          call waterbalance(d,present,rootprop,w,w_t,dgp,dpet,dphen,
     *      dgc,dmelt,dprec,k,fwhc,drunoff_drain,drunoff_surf,dwscal,
     *      daet,mpar_day,m,maxthaw_depth,dthaw_depth,uw1,uw2,fw,
     *      year,fpc_grid,daep,d_evap,w_ep,tree,lai_ind,needle,boreal,
     *      intercp_tot,dayl,lat,gminp,emax,intc,pet_grid)

c         Doug 07/09: calculate mositure availability index for waterstress (alpha=sum(daet/peat))
c          alpha_ws=alpha_ws+daet/dpet(d) !Doug 08/09: changed formulea, below and at end of loop
c         Doug 08/09: alpha_ws=sum(daet)/sum(deet). the sum in daet is added in here, and the 2nd
c             part, /sum(deet) is done at the end of the month loop
          alpha_ws=alpha_ws+daet


c         Store today's water content in soil layer 1
          dw1(d)=w(1) 
      

c          dw1(d)=w_ep !upper 20 cm only
          dw2(d)=w(2)

c         Increment monthly totals

          mrunoff_surf(m)=mrunoff_surf(m)+drunoff_surf
          mrunoff_drain(m)=mrunoff_drain(m)+drunoff_drain
          aaet=aaet+daet
          maet(m)=maet(m)+daet
          aaep=aaep+daep
          maep(m)=maep(m)+daep
          aintc=aintc+intercp_tot
          mintc(m)=mintc(m)+intercp_tot
          mpet_grid(m)=mpet_grid(m)+pet_grid !PET grid means
          apet_grid=apet_grid+pet_grid
          mpet2(m)=mpet2(m)+dpet(d)  !"radiative" PET only
          apet=apet+dpet(d)
          mmelt(m)=mmelt(m)+dmelt(d)

c         Increment monthly w total and monthly u total
          mw1(m)=mw1(m)+w(1)
          mw2(m)=mw2(m)+w(2)
          mw1_t(m)=mw1_t(m)+w_t(1)
          mw2_t(m)=mw2_t(m)+w_t(2)
          muw1(m)=muw1(m)+uw1
          muw2(m)=muw2(m)+uw2

          do pft=1,npft
            if (present(pft)) then

c             Accumulate count of days with some leaf cover and
c             pft-specific annual water scalar used in allocation

              if (dphen(d,pft).gt.0.0) then
                aleafdays(pft)=aleafdays(pft)+1
                awscal(pft)=awscal(pft)+dwscal(pft)
              endif

c             Accumulate mean monthly fpc, actual (gc) and minimum
c             (gmin) canopy conductances, incorporating leaf
c             phenology

              meanfpc(m,pft)=meanfpc(m,pft)+
     *          fpc_grid(pft)*dphen(d,pft)/real(ndaymonth(m))
	 
	 

	 
              meangc(m,pft)=meangc(m,pft)+
     *          dgc(d,pft)/real(ndaymonth(m))
              meangmin(m,pft)=meangmin(m,pft)+
     *          gminp(pft)*dphen(d,pft)/real(ndaymonth(m))
              mphen(m,pft)=mphen(m,pft)+dphen(d,pft)/real(ndaymonth(m))


c            Save final daily water scalar for next year
             if (d.eq.365) dwscal365(pft)=dwscal(pft)
            endif
          enddo
        enddo !day of month

c       CALCULATE monthly average soil content for soil layers 1 and 2 (mm)
        mauw1(m)=muw1(m)/real(ndaymonth(m))
        mauw2(m)=muw2(m)/real(ndaymonth(m))

c       Increment annual runoff totals

        mrunoff(m)=mrunoff_surf(m)+mrunoff_drain(m)
        arunoff_surf=arunoff_surf+mrunoff_surf(m)
        arunoff_drain=arunoff_drain+mrunoff_drain(m)

c       CALCULATE GPP FOR EACH PFT
c       Fnind water-limited daily net photosynthesis (And) and
c       ratio of intercellular to ambient partial pressure of CO2
c       (lambda) by solving simultaneously Eqns 2, 18 and 19
c       (Haxeltine & Prentice 1996).

c       Using a tailored implementation of the bisection method
c       with a fixed 10 bisections, assuming root (f(lambda)=0)
c       bracketed by f(0.02)<0 and f(lambdam+0.05)>0.

        do pft=1,npft

          if (present(pft)) then

c           Convert canopy conductance assoc with photosynthesis
c           (actual minus minimum gc) (gpd) from mm/sec to mm/day
            gpd=tsecs(m)*(meangc(m,pft)-meangmin(m,pft))
            fpar=meanfpc(m,pft)  !cover including phenology

            if (gpd.gt.1.0e-5) then

c             Implement numerical solution

              x1=0.02           !minimum bracket of the root
              x2=lambdam(pft)+0.05   !maximum bracket of the root
              rtbis=x1          !root of the bisection
              dx=x2-x1

              b=0  !number of tries towards solution
              fmid=epsilon+1.0
			  
              do while (abs(fmid).gt.epsilon.and.b.le.10)

c               Abort search if >10 iterations needed to fnind solution

                b=b+1
                dx=dx*0.5
                xmid=rtbis+dx

c               Calculate total daytime photosynthesis implied by
c               canopy conductance from water balance routine and
c               current guess for lambda (xmid).  Units are mm/m2/day
c               (mm come from gpd value, mm/day)
c               Eqn 18, Haxeltine & Prentice 1996

                adt1=(gpd/1.6)*(ca*(1.0-xmid))

c               Call photosynthesis to determine alternative total
c               daytime photosynthesis estimate (adt2) implied by
c               Eqns 2 & 19, Haxeltine & Prentice 1996, and current
c               guess for lambda (xmid)

                call photosynthesis(ca,mtemp(m),fpar,mpar_day(m),
     *            mdayl(m),c4(pft),sla(pft),nmax(pft),xmid,rd,agd,adt2,
     *            inhibx1(pft),inhibx2(pft),inhibx3(pft),inhibx4(pft),
     *            pft)   

c               Evaluate fmid at the point lambda=xmid
c               fmid will be an increasing function with xmid, with
c               a solution (fmid=0) between x1 and x2

                fmid=adt2-adt1

                if (fmid.lt.0.0) rtbis=xmid

              enddo !bisection

            else  !infinitesimal canopy conductance

              agd=0.0
              rd=0.0
              xmid=0.0

            endif


c           Estimate monthly gross photosynthesis and monthly leaf
c           respiration from mid-month daily values
c           Agd = And + Rd (Eqn 2 Haxeltine & Prentice 1996)
            mgpp(m,pft,1)=real(ndaymonth(m))*agd
            agpp(pft,1)=agpp(pft,1)+mgpp(m,pft,1)

            mlresp(m,pft,1)=real(ndaymonth(m))*rd
            alresp(pft,1)=alresp(pft,1)+mlresp(m,pft,1)

            mcica(m,pft)=xmid

            mapar(m)=mapar(m)+mpar_day(m)*fpar*real(ndaymonth(m))*0.5
          endif
        enddo  !pft

        mpar(m)=mpar(m)+mpar_day(m)*real(ndaymonth(m))
        mw1(m)=mw1(m)/real(ndaymonth(m))
        mw2(m)=mw2(m)/real(ndaymonth(m))
        mw1_t(m)=mw1_t(m)/real(ndaymonth(m))
        mw2_t(m)=mw2_t(m)/real(ndaymonth(m))
      enddo  !month

c Doug 08/09: part 2 in calculating mositure availability index
      alpha_ws=alpha_ws/SUM(dpet)


c Doug 05/09: For fire paradox experiment (http://www.fireparadox.org/).
c Experiment 1&3: phenology=0 for grass for an extra month.
c 
C     
c      IF (year>=5180.AND.year<5190) THEN
c      IF (year>=5140.AND.year<5190) THEN
c      IF (year>=5160) THEN
c        day_dphen=365
c        d=0
c        DO pft=1,npft
c         IF(tree(pft)==.FALSE.) THEN
c            DO m=1,12
c              DO dm=1,ndaymonth(m)
c                d=d+1
c                IF (d==day_dphen) GOTO 300
c                IF (dphen_change(m,pft)==d) THEN
c                  day_dphen=d
c                  GOTO 300
c                END IF
c              END DO
c            END DO
c          END IF
c300       CONTINUE
c        END DO
c 
c       DO pft=1,npft
c          IF (tree(pft)==.FALSE.) THEN
c            IF (day_dphen+30>365) THEN
c              dphen(day_dphen:365,pft)=0
c            ELSE
c              dphen(day_dphen:day_dphen+30,pft)=0
c            END IF
c          END IF
c        END DO
c      END IF
cc 
c Experiment 2: 
c 
c      IF (year>=5180) THEN
c        DO m=1,12
c          DO pft=1,npft
c            IF(tree(pft)==.FALSE.) THEN
c              dphen(ndaymonth(m):ndaymonth(m)+7,pft)=0
c            END IF
c          END DO
c        END DO
c      END IF
 
 
c     Convert water scalar to average daily basis using phenology
       do pft=1,npft
        if (present(pft)) then
          if (aleafdays(pft).ne.0) then
            wscal(pft)=awscal(pft)/real(aleafdays(pft))
          else
            wscal(pft)=1.0
          endif
        endif
      enddo


c     Calculate the carbon isotope fractionation in plants following Lloyd 
c     & Farquhar 1994
  
      deltaa_fpc=0.0
 
      do pft=1,npft
         if (present(pft)) then          
            if (agpp(pft,1).gt.0.0) then 
               call isotope(pft,mcica,co2,mtemp,mlresp,
     *    c4,mgpp,agpp,fpc_grid,deltaa)


        if(.not.tree(pft))then
          deltaa_fpc=deltaa_fpc+deltaa(pft)*fpc_grid(pft)
!              write(*,*)'deltaa=',deltaa_fpc
        end if


            else
               do m=1,12
                 mgpp(m,pft,2)=0.0
                 mgpp(m,pft,3)=0.0
               enddo
               agpp(pft,2)=0.0
               agpp(pft,3)=0.0
            endif
         endif      
      enddo

c     Calculate annual total runoff

      arunoff=arunoff_surf+arunoff_drain
	  
 
      return
      end



c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE PHOTOSYNTHESIS
c     Adapted from Farquhar (1982) photosynthesis model, as simplified by
c     Collatz et al 1991, Collatz et al 1992 and Haxeltine & Prentice 1996

      subroutine photosynthesis(ca,temp,fpar,par,dayl,c4,sla,nmax,
     *  lambda,rd,agd,adtmm,x1,x2,x3,x4,pfti)

      implicit none

c     PARAMETERS
      real po2,p,bc3,bc4,theta,q10ko,q10kc,q10tau,ko25,kc25,tau25
      real alphaa,alphac3,alphac4,cmass,lambdamc4,lambdamc3,n0,m
      real e0,t0c3,t0c4,cq,tk25,tmc4,tmc3
      integer npft

      parameter (npft=13)
      parameter (po2=20.9e3)!O2 partial pressure in Pa
      parameter (p=1.0e5)   !atmospheric pressure in Pa
      parameter (bc3=0.015) !leaf respiration as fraction of Vmax for C3 plants
      parameter (bc4=0.02)  !leaf respiration as fraction of vmax for C4 plants
      parameter (theta=0.7) !colimitation (shape) parameter
      parameter (q10ko=1.2) !q10 for temperature-sensitive parameter ko
      parameter (q10kc=2.1) !q10 for temperature-sensitive parameter kc
      parameter (q10tau=0.57) !q10 for temperature-sensitive parameter tau
      parameter (ko25=3.0e4)  !value of ko at 25 deg C
      parameter (kc25=30.0)   !value of kc at 25 deg C
      parameter (tau25=2600.0)!value of tau at 25 deg C
      parameter (alphaa=0.5)  !fraction of PAR assimilated at ecosystem level
                              !relative to leaf level
      parameter (alphac3=0.08)  !intrinsic quantum efficiency of CO2 uptake in
                                !C3 plants
      parameter (alphac4=0.053) !C4 intrinsic quantum efficiency
      parameter (cmass=12.0)    !atomic mass of carbon
      parameter (cq=4.6e-6)   !conversion factor for solar radiation at 550 nm
                              !from J/m2 to E/m2 (E=mol quanta)
      parameter (lambdamc4=0.4)!optimal ratio of intercellular to ambient CO2
                               !concentration (lambda) in C4 plants
      parameter (lambdamc3=0.8)!optimal (maximum) lambda in C3 plants
      parameter (n0=7.15)     !leaf N concentration (mg/g) not involved in
                              !photosynthesis
      parameter (m=25.0)      !corresponds to parameter p in Eqn 28, Haxeltine
                              !& Prentice 1996
      parameter (t0c3=250.0)  !base temperature (K) in Arrhenius temperature
                              !response function for C3 plants
      parameter (t0c4=260.0)  !base temperature in Arrhenius func for C4 plants
      parameter (e0=308.56)   !parameter in Arrhenius temp response function
      parameter (tk25=298.15) !25 deg C in Kelvin
      parameter (tmc3=45.0)   !maximum temperature for C3 photosynthesis
      parameter (tmc4=55.0)   !maximum temperature for C4 photosynthesis

c     ARGUMENTS
      real ca,temp,fpar,par,dayl
      logical c4
      real sla,nmax
      real lambda,rd,agd,adtmm
      real x1,x2,x3,x4
      integer pfti

c     LOCAL VARIABLES
      integer i 
      real lambdam(1:npft)
      data (lambdam(i),i=1,npft)
     *     /0.9,0.9,0.9,0.8,0.8,0.8,0.9,0.65,0.4,0.9,0.9,0.8,0.8/

      real apar,pi,ko,kc,tau,gammastar,c1,c2,b,phipi
      real s,sigma,vm,and,pa,je,jc,vmmax,cn,tk,t0,adt
      real tstress,k1,k2,k3,low,high
      
      i=0
      apar=0.0
      pi=0.0
      ko=0.0
      kc=0.0
      tau=0.0
      gammastar=0.0
      c1=0.0
      c2=0.0
      b=0.0
      phipi=0.0
      s=0.0
      sigma=0.0
      vm=0.0
      and=0.0
      pa=0.0
      je=0.0
      jc=0.0
      vmmax=0.0
      cn=0.0
      tk=0.0
      t0=0.0
      adt=0.0
      tstress=0.0
      k1=0.0
      k2=0.0
      k3=0.0
      low=0.0
      high=0.0
      

c     Return without performing calculations if daylength 0 hours

      if (dayl.lt.0.01) then
        agd=0.0
        adtmm=0.0
        rd=0.0
        return
      endif

c     APAR in J/m2/day
c     alphaa = scaling factor for absorbed PAR at ecosystem, versus leaf, scale
c     See Eqn 4, Haxeltine & Prentice 1996

      apar=par*fpar*alphaa

c     calculate temperate inhibition function

      if (temp.lt.x4) then
         k1=2.*ALOG((1./0.99)-1.)/(x1-x2)
         k2=(x1+x2)/2.
         low=1./(1.+ EXP(k1*(k2-temp)))
         k3=ALOG(0.99/0.01)/(x4-x3)
         high=1.-0.01*EXP(k3*(temp-x3))
         tstress=(low*high)
      else
         tstress=0.
      endif
      if (tstress.lt.1e-2) tstress=0.

c     First calculate catalytic capacity of rubisco, Vm, assuming optimal
c     (non-water-stressed) value for lambda, i.e. lambdamc3

      if (.not.c4) then  !C3 photosynthesis

c       Temperature-adjusted values of kinetic parameters, Eqn 22,
c       Haxeltine & Prentice 1996a

        ko=ko25*q10ko**((temp-25.0)/10.0) !Michaelis constant of rubisco for O2
        kc=kc25*q10kc**((temp-25.0)/10.0) !Michaelis constant for CO2
        tau=tau25*q10tau**((temp-25.0)/10.0) !CO2/O2 specificity ratio

c       CO2 compensation point (CO2 partial pressure, Pa)
c       Eqn 8, Haxeltine & Prentice 1996

        gammastar=po2/(2.0*tau)

c       Convert ambient CO2 level, ca, from mole fraction to partial pressure
c       in Pa

        pa=ca*p

c       Non-water-stressed intercellular CO2 partial pressure in Pa
c       Eqn 7, Haxeltine & Prentice 1996

        pi=lambdam(pfti)*pa

c       Calculation of C1C3, Eqn 4, Haxeltine & Prentice 1996
c       Notes: - there is an error in this equation in the above paper (missing
c                2.0* in denominator) which is fixed here (see Eqn A2, Collatz
c                et al 1991)
c              - There is no longer an explicit temperature inhibition function
c                (low-temperature inhibition is now done mechanistically
c                by imposing a temperature-dependent upper limit on Vm, see
c                below)
c              - There is no longer any reduction in maximum photosynthesis due
c                to leaf age (phic)
c              - alphaa, the upscaling parameter accounting for the reduction
c                in PAR utilisation in ecosystems compared with leaf level,
c                appears in the calculation of APAR instead of here
c              - Cmass, the atomic weight of carbon, used in unit conversion
c                from molC to g appears in the calculation of Vm instead of
c                here

        c1=tstress*alphac3*((pi-gammastar)/(pi+2.0*gammastar))

c       High temperature inhibition modelled primarily by suppression of LUE
c       by decreased relative affinity of rubisco for CO2 relative to O2 with
c       increasing temperature, but we also implement a step function to
c       prohibit any C3 photosynthesis above 45 degrees (Table 3.7, Larcher
c       1983)

        if (temp.gt.tmc3) c1=0.0

c       Calculation of C2C3, Eqn 6, Haxeltine & Prentice 1996

        c2=(pi-gammastar)/(pi+kc*(1.0+po2/ko))

        b=bc3   !Choose C3 value of b for Eqn 10, Haxeltine & Prentice 1996
        t0=t0c3 !base temperature for temperature response of rubisco

      else  !C4 photosynthesis

c       Specify C1C4, C2C4
c       Eqns 14,15, Haxeltine & Prentice 1996
c       Notes:
c              - alphaa, the upscaling parameter accounting for the reduction
c                in PAR utilisation in ecosystems compared with leaf level,
c                appears in the calculation of APAR instead of here
c              - Cmass, the atomic weight of carbon, used in unit conversion
c                from molC to g appears in the calculation of Vm instead of
c                here
c              - parameter phipi is not needed for calculation of optimal Vm
c                which assumes optimal intercellular CO2 concentration
c                (lambdamc4)

        c1=tstress*alphac4

c       High-temperature inhibition modelled conservatively as a step function
c       prohibiting photosynthesis above 55 deg C (Table 3.7, Larcher 1983)

        if (temp.gt.tmc4) c1=0.0

        c2=1.0
        b=bc4    !Choose C4 value of b for Eqn 10, Haxeltine & Prentice 1996
        t0=t0c4  !base temperature for temperature response of rubisco

      endif

c     Eqn 13, Haxeltine & Prentice 1996
      s=(24.0/dayl)*b

c     Eqn 12, Haxeltine & Prentice 1996
      sigma=sqrt(max(0.0,1.0-(c2-s)/(c2-theta*s)))

c     Calculation of optimal rubisco capacity, Vm, in gC/m2/day
c     Eqn 11, Haxeltine & Prentice 1996
      vm=(1.0/b)*(c1/c2)*((2.0*theta-1.0)*s-(2.0*theta*s-
     *  c2)*sigma)*apar*cmass*cq

c     Now use this Vm value to calculate actual photosynthesis
      if (.not.c4) then  !C3 photosynthesis

c       Intercellular CO2 partial pressure in Pa
c       Eqn 7, Haxeltine & Prentice 1996
        pi=lambda*pa

c       Recalculation of C1C3, C2C3 with actual pi
        c1=tstress*alphac3*((pi-gammastar)/(pi+2.0*gammastar))
        if (temp.gt.tmc3) c1=0.0  !high-temperature inhibition

        c2=(pi-gammastar)/(pi+kc*(1.0+po2/ko))

      else  !C4 photosynthesis

c       Parameter accounting for effect of reduced intercellular CO2
c       concentration on photosynthesis, Phipi.
c       Eqn 14,16, Haxeltine & Prentice 1996
c       Fig 1b, Collatz et al 1992

        phipi=min(lambda/lambdam(pfti),1.0)
        c1=tstress*phipi*alphac4
        if (temp.gt.tmc4) c1=0.0  !high-temperature inhibition

      endif

c     je is PAR-limited photosynthesis rate molC/m2/h, Eqn 3
c     Convert je from daytime to hourly basis

c     Calculation of PAR-limited photosynthesis rate, JE, molC/m2/h
c     Eqn 3, Haxeltine & Prentice 1996

      je=c1*apar*cmass*cq/dayl

c     Calculation of rubisco-activity-limited photosynthesis rate JC, molC/m2/h
c     Eqn 5, Haxeltine & Prentice 1996

      jc=c2*vm/24.0

c     Calculation of daily gross photosynthesis, Agd, gC/m2/day
c     Eqn 2, Haxeltine & Prentice 1996
c     Note: - there is an error in this equation in the above paper (missing
c             theta in 4*theta*je*jc term) which is fixed here

      agd=(je+jc-sqrt((je+jc)**2.0-4.0*theta*je*jc))/(2.0*theta)*dayl

c     Daily leaf respiration, Rd, gC/m2/day
c     Eqn 10, Haxeltine & Prentice 1996

      rd=b*vm

c     Daily net photosynthesis (at leaf level), And, gC/m2/day

      and=agd-rd

c     Total daytime net photosynthesis, Adt, gC/m2/day
c     Eqn 19, Haxeltine & Prentice 1996

      adt=and+(1.0-dayl/24.0)*rd

c     Convert adt from gC/m2/day to mm/m2/day using
c     ideal gas equation

      adtmm=adt/cmass*8.314*(temp+273.3)/p*1000.0

      return
      end



c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE WATERBALANCE
c     Calculations of soil hydrology, plant water balance, daily canopy
c     conductance and PFT daily water stress factor

      subroutine waterbalance(d,present,rootprop,w,w_t,dgp,dpet,dphen,
     *  dgc,dmelt,dprec,k,fwhc,drunoff_drain,drunoff_surf,dwscal,daet,
     *  mpar_day,m,maxthaw_depth,dthaw_depth,uw1,uw2,fw,
     *  year,fpc_grid,daep,d_evap,w_ep,tree,lai_ind,needle,boreal,
     *  intercp_tot,dayl,lat,gminp,emax,intc,pet_grid)

      implicit none

c     PARAMETERS:
      integer npft,npftpar,nsoilpar
        parameter (npft=13,npftpar=59,nsoilpar=7)
      real alpham,gm,pt          ! ,emax
c        parameter (emax=5.0)    ! maximum daily transpiration rate (mm/day)
        parameter (gm=3.26)      ! empirical parameter in demand function
        parameter (alpham=1.391) ! Priestley-Taylor coefficient (Demand)
        parameter (pt=1.32)      ! Priestley-Taylor coefficient
        real d1,d2
        parameter (d1=500.0,d2=1000.0)

c     ARGUMENTS:
      integer d
      logical present(1:npft)
      real rootprop(1:2,1:npft),w(1:2),w_t(1:2),dgp(1:365,1:npft)
      real dpet(1:365),dphen(1:365,1:npft),dgc(1:365,1:npft)
      real dmelt(1:365),dprec(1:365),k(1:2),fwhc(1:2)
      real drunoff_drain,drunoff_surf,dwscal(1:npft)
      real daet,mpar_day(1:12)
      integer m
      real maxthaw_depth,dthaw_depth(1:365)
      real uw1,uw2
      real fw(1:2)     ! ice content (mm) in soil layers
      integer year     ! simulation year
      real fpc_grid(1:npft)
      real daep        ! daily evaporation (mm)
      real d_evap      ! soil evaporation layer (mm)
      real w_ep        ! fractional soil moisture in evaporation layer
      logical tree(1:npft)
      real lai_ind(1:npft)
      logical needle(1:npft), boreal(1:npft)
      real intercp_tot !total interception loss (mm)
      real dayl(1:365)
      real lat
      real gminp(1:npft)
      real emax(1:npft)
      real intc(1:npft) !LAI parameter (Kergoat 1998, Tab. 1)
      real pet_grid

c     LOCAL
      integer pft
      real aettotal(1:2)
      real wet                   !fraction of day when canopy is wet
      real int_store(1:npft)     !water stored in the canopy
      real intercep(1:npft)
      real aet(1:npft)
      real wr,beta(1:2),perc,perc_evap
      real f1,f2,delta_thaw
      real supply,demand,demandpot,demand_tot,supply_tot
      real maxthaw   ! maximum thawing depth (ltd. to 1500 mm)
      real freez     ! current freezing depth of both layers (mm)
      real dthaw     ! current thawing depth of both layers (mm)
      real thaw(1:2),thaw_old(1:2) ! today's/yesterday's thawing depths (mm)
      real thaw_ep,thawep_old ! today's/yesterday's thawing depth in evap. layer (mm)
      real fw_ep     ! frozen water in evaporation layer (mm)
      real uw_ep     ! water content in evaporation layer (mm)
      real excess    ! excess water in evaporation layer (mm)
      real v(1:2)    ! fraction of thawing depth in each soil layer
      real cover     ! currently plant-covered fraction of soil
      integer case
      real pet_s,pet_p,q,prec


c     Initialisations
      pft=0
      aettotal(:)=0.0
      wet=0.0
      intercep(:)=0.0
      int_store(:)=0.0
      aet(:)=0.0
      wr=0.0
      beta(:)=0.0
      perc=0.0
      perc_evap=0.0
      f1=0.0
      f2=0.0
      delta_thaw=0.0
      supply=0.0
      demand=0.0
      demandpot=0.0
      demand_tot=0.0
      supply_tot=0.0
      maxthaw=0.0
      freez=0.0
      dthaw=0.0
      thaw(:)=0.0
      thaw_old(:)=0.0
      thaw_ep=0.0
      thawep_old=0.0
      fw_ep=0.0
      uw_ep=0.0
      excess=0.0
c DM  TODO reset to 0 ?   Yan: No! We don't take into permafost,v(1)=thaw/d1;v(2)=thaw/d2;==1 
      v(:) =1.0
c Yan
      cover=0.0
      case=0.0
      pet_s=0.0
      pet_p=0.0
      q=0.0
      prec=0.0
      daet=0.0
      intercp_tot=0.0
      pet_grid=0.0
      drunoff_drain=0.0
      drunoff_surf=0.0
      
      if(d.eq.1) then
       dthaw=1500.0
       maxthaw=1500.0
      endif

c We assume that fractional water holding capacity exhibits only minor
c differences between lower and upper bucket

      do pft=1,npft

c   CALCULATE ACTUAL CANOPY CONDUCTANCE, POTENTIAL WATER SCALAR,
c   INTERCEPTION LOSS AND ACTUAL EVAPOTRANSPIRATION (AET) FOR EACH PFT

        if (present(pft)) then

c      Calculate effective supply function in root zone today
c      f2=fraction of active fraction of roots uptaking water
c         from lower soil layer
c      f1=the remainder in the top layer

          if (maxthaw_depth.ge.d1+d2) then
             f2=rootprop(2,pft)
             f1=rootprop(1,pft)
          else
             f2=rootprop(2,pft)*v(2)
             f1=1.0-f2
          endif

          wr=w(1)*f1*v(1)+w(2)*f2*v(2)

          supply=emax(pft)*wr*dphen(d,pft)*fpc_grid(pft)
c          if(wr.le.0.7) then
c            supply=dpet(d)*pt*wr*dphen(d,pft)*fpc_grid(pft)
c          else
c            supply=dpet(d)*pt*dphen(d,pft)*fpc_grid(pft)
c          endif

c Calculate PFT-specific INTERCEPTION LOSS (Kergoat 1998)
          int_store(pft)=min(dprec(d),
     *    (intc(pft)*lai_ind(pft)*dphen(d,pft)*dprec(d)))
          if(dpet(d).lt.0.0001) then
            wet=0.0
          else
             wet=min(int_store(pft)/(dpet(d)*pt),0.99)
          endif
          intercep(pft)=dpet(d)*pt*wet

c  Calculate AET demand function and potential demand
c  assuming full leaf cover (Eqn 23, Haxeltine & Prentice 1996)
c  Note fraction of daylength
          if(dgp(d,pft).eq.0.0) dgp(d,pft)=0.0001
          if (dphen(d,pft).gt.0.0) then
c            demand=(1.0-wet)*dpet(d)*alpham*(1.0-
c     *        exp(-dgp(d,pft)*dphen(d,pft)/gm))
            demand=(1.0-wet)*dpet(d)*alpham/
     *             (1+gm/(dgp(d,pft)*dphen(d,pft)))
          else
            demand=0.0
          endif
c          demandpot=dpet(d)*alpham*(1.0-exp(-dgp(d,pft)/gm))
          demandpot=dpet(d)*alpham/(1+gm/dgp(d,pft))

c         Calculate daily potential water scalar
          if (demand.ne.0.0) then
            dwscal(pft)=min((supply/dphen(d,pft))/demandpot,1.0)
          else
            dwscal(pft)=1.0
          endif
c  Canopy conductance according to balance between supply and demand
          if(supply.ge.demand) then
            dgc(d,pft)=dgp(d,pft)*dphen(d,pft)
          else
            if(dpet(d).gt.0.0) then
c              dgc(d,pft)=max(0.0,-gm*log(1.0-supply/
c     *        ((1.0-wet)*dpet(d)*alpham)))
              dgc(d,pft)=max(0.0,(
     *           gm*supply/((1.0-wet)*dpet(d)*alpham
     *           -supply)))
            else
              dgc(d,pft)=0.0
            endif
          endif
c
          aet(pft)=min(supply,demand)

c         Accumulate total AET and interception loss
          if (wr.eq.0.0) then
            beta(1)=0.0
            beta(2)=0.0
          else
            beta(1)=f1*w(1)*v(1)/wr
            beta(2)=f2*w(2)*v(2)/wr
          endif

          aettotal(1)=aettotal(1)+beta(1)*aet(pft)
          aettotal(2)=aettotal(2)+beta(2)*aet(pft)
          intercp_tot=intercp_tot+intercep(pft)*fpc_grid(pft)
          demand_tot=demand_tot+demand
          supply_tot=supply_tot+supply

c Calculate area covered by PFTs
          cover=cover+(fpc_grid(pft)*dphen(d,pft))
        endif
      enddo     !pft

c Compute total transpiration and effective rainfall
      daet=aettotal(1)+aettotal(2)
      dprec(d)=dprec(d)-intercp_tot

c Permafrost (Northern Hemisphere): thawing or freezing
         if (d.ne.1) then
           delta_thaw=dthaw_depth(d)-dthaw_depth(d-1)
         else
           delta_thaw=-1.0
         endif
      if (maxthaw_depth.eq.1500.and.dthaw_depth(d).eq.1500) then
      delta_thaw=0.0
      endif

c Constrain thawing depths above 1500
        if (maxthaw_depth.ge.1500) then
           maxthaw=1500.0
        else
           maxthaw=maxthaw_depth
        endif
        if (dthaw_depth(d).ge.1500) then
           dthaw=1500.0
        else
           dthaw=dthaw_depth(d)
        endif

c Redefine daily thawing depth and freezing depth
        thawep_old=thaw_ep
        thaw_old(1)=thaw(1)
        thaw_old(2)=thaw(2)
        if (delta_thaw.lt.0.0) then
           freez=maxthaw-dthaw
              if(freez.le.d1) then
                thaw(1)=d1-freez
                thaw(2)=maxthaw-d1
                if(thaw(2).le.0.0) thaw(2)=0.0
                if(thaw(1).eq.d1) then
                 thaw_ep=d_evap
                else
                 thaw_ep=0.0
                endif
                case=1          ! freezing upper layer
              else
                thaw(1)=0.0
                thaw_ep=0.0
                thaw(2)=maxthaw-freez
                case=2          ! freezing lower layer
              endif
        else
              if(dthaw.le.d1) then
                thaw(1)=dthaw
                 if(dthaw.lt.d_evap) then
                  thaw_ep=dthaw
                 else
                  thaw_ep=d_evap
                 endif
                thaw(2)=0.0
                case=3          ! thawing upper layer
              else
                thaw(1)=d1
                thaw_ep=d_evap
                thaw(2)=dthaw-d1
                case=4          ! thawing lower layer
              endif
        endif

        pet_s=dpet(d)*pt*(1.0-cover)
        pet_grid=pet_s+demand_tot+intercp_tot

c CALCULATE SOIL WATER COMPONENTS
      if(case.lt.3) then ! ice near surface,i.e. all prec=runoff
       drunoff_surf=dprec(d)+dmelt(d)
       uw1=uw1-aettotal(1)
       uw2=uw2-aettotal(2)
       fw(1)=uw1*(d1-thaw(1))/d1
       fw(2)=uw2*(d2-thaw(2))/d2
       w(1)=uw1/(fwhc(1)*d1)
       w(2)=uw2/(fwhc(2)*d2)

      else
c       drunoff_surf=(((w(1)+w(2))/2)**2)*(dprec(d)+dmelt(d))
       prec=dprec(d)+dmelt(d)-drunoff_surf
      daep=dpet(d)*pt*(w_ep)*(thaw_ep/d_evap)*(1.0-cover)
       uw1=uw1+prec-aettotal(1)-daep                              !uw1=Pr-ET-ES
       if(thaw(1).le.0.0001) then
          w(1)=1.0
          w_t(1)=1.0
       else
         w(1)=(uw1-fw(1))/(fwhc(1)*thaw(1)) !liquid water only
         w_t(1)=uw1/(fwhc(1)*d1) !both
       endif
       if (w(1).ge.1.0) then
            drunoff_surf=(w(1)-1.0)*fwhc(1)*thaw(1)+drunoff_surf
            uw1=uw1-((w(1)-1.0)*fwhc(1)*thaw(1))
            w(1)=1.0
       endif

c water content in evaporation layer
       uw_ep=uw_ep+prec-(aettotal(1)*(d_evap*1.3/d1))
     *       -daep
       if(uw_ep.le.0.0001) uw_ep=0.0
       uw_ep=min(uw1,uw_ep)
       if(thaw_ep.le.0.0001) then
         w_ep=0.0
       else
         w_ep=(uw_ep-fw_ep)/(fwhc(1)*thaw_ep)
       endif

       if (w_ep.ge.1.0) then
          excess=(w_ep-1.0)*fwhc(1)*thaw_ep
          uw_ep=uw_ep-excess
          w_ep=1.0
       endif

       if(v(2).gt.0.0.and.prec.ge.0.1) then
         perc=k(1)*(w(1)*v(1))**k(2)
         if(perc.gt.prec) perc=prec
         perc_evap=k(1)*2/5*(w_ep*v(1))**k(2)
         if(perc_evap.gt.prec) perc_evap=prec
         q=0.5*k(1)*(w(2)*v(2))**k(2)
         if(q.gt.(dprec(d)-drunoff_surf))
     *   q=dprec(d)-drunoff_surf
       else
         perc=0.0
         perc_evap=0.0
       endif

       uw1=uw1-perc
       uw_ep=uw_ep-perc_evap

       if(thaw_ep.le.0.0001) then
         w_ep=0.0
       else
         w_ep=(uw_ep-fw_ep)/(fwhc(1)*thaw_ep)
       endif
       if(thaw(1).le.0.0001) then
         w(1)=1.0
         w_t(1)=1.0
       else
         w(1)=(uw1-fw(1))/(fwhc(1)*thaw(1))
         w_t(1)=uw1/(fwhc(1)*d1)
       endif

       uw2=uw2+perc-aettotal(2)-q
            if (uw1.le.0.0) uw1=0.0
            if (uw2.le.0.0) uw2=0.0
            if(thaw(2).le.0.0001) then
               w(2)=1.0
            else
               w(2)=(uw2-fw(2))/(fwhc(2)*thaw(2))
            endif

c Calculate runoff from lower layer
        if(w(2).ge.1.0) then
         drunoff_drain=((w(2)-1.0)*fwhc(2)*thaw(2))+q
         uw2=uw2-((w(2)-1.0)*fwhc(2)*thaw(2))
         w(2)=(uw2-fw(2))/(fwhc(2)*thaw(2))
        else
         drunoff_drain=q
         w(2)=(uw2-fw(2))/(fwhc(2)*thaw(2))
        endif

       if((d_evap-thaw_ep).le.0.0001) then
        fw_ep=0.0
       else
        fw_ep=fw_ep-(fw_ep*(thaw_ep-thawep_old)/(d_evap-thaw_ep))
       endif
       if((d1-thaw(1)).le.0.0001) then
        fw(1)=0.0
       else
        fw(1)=fw(1)-(fw(1)*((thaw(1)-thaw_old(1))/(d1-thaw(1))))
       endif
       if((d2-thaw(2)).le.0.0001) then
        fw(2)=0.0
       else
       fw(2)=fw(2)-(fw(2)*((thaw(2)-thaw_old(2))/(d2-thaw(2))))
       endif
      endif  !case

      v(1)=thaw(1)/d1
      v(2)=thaw(2)/d2
      if (w(1).le.0.0001) w(1)=0.0
      if (w(2).le.0.0001) w(2)=0.0
    

      return
      end



c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE ISOTOPE
c     This subroutine is for calculating the total discrimination of 13C
c     as it goes from free air (as 13CO2) to fixed carbon in the leaf.
c     This program is based upon the model used by Lloyd and Farquhar (1994).

      subroutine isotope(pft,cratio,co2,mtemp,rd,c4,mgpp,agpp,
     *                   fpc_grid,deltaa)


      implicit none

c     PARAMETERS
      integer npft,nco2
      parameter (npft=13,nco2=3)
      real a,es,a1,b,b3,f,e
      parameter (a=4.4,es=1.1,a1=0.7,b=27.5,b3=30.0,e=0.0,f=8.0)


   

c     ARGUMENTS
      integer pft
      real cratio(1:12,1:npft),co2(1:nco2),mtemp(1:12)
      real rd(1:12,1:npft,1:nco2)
      logical c4(1:npft)
      real mgpp(1:12,1:npft,1:nco2),agpp(1:npft,1:nco2)
      

c     LOCAL
      integer m     


!      real deltaa,b4
      real deltaa(1:npft),b4

      real leaftemp,gamma
      real catm,k,rdi
      real q,r,s,t
      real phi
      real c3gpp,c4gpp  !gpp for C3 and C4 plants.
      real fpc_grid(1:npft)
      integer i

      m=0

       do i=1,npft
      deltaa(i)=0.0
       end do

      b4=0.0
      leaftemp=0.0
      gamma=0.0
      catm=0.0
      k=0.0
      rdi=0.0
      q=0.0
      r=0.0
      s=0.0
      t=0.0
      phi=0.4                 !as suggested by Lloyd & Farquhar (1994)
      c3gpp=0.0
      c4gpp=0.0
   
      do m=1,12

        if (mgpp(m,pft,1).gt.0.0) then

           if (cratio(m,pft).lt.0.05) cratio(m,pft)=0.05
          
c          define discrimination parameters
     
           if (rd(m,pft,1).le.0.0) rd(m,pft,1)=0.01
           
           leaftemp = 1.05*(mtemp(m)+2.5)
           gamma = 1.54*leaftemp
           rdi = rd(m,pft,1)/(86400.0*12.0)
           catm = co2(1)*1.0e6
           k = rdi/11.0         
           b4=(26.19-(9483/(273.2+leaftemp))) !from Farquhar et al. 1982 p. 126
           
           if (c4(pft)) then
          
c             calculate PHI for C4 pathway 

c              call calcphi(mgpp,pft,phi)

c             calculate discrimination -> C4 pathway           
              

              deltaa(pft)=a*(1-cratio(m,pft)+0.0125)+0.0375*(es+a1)+
     *             (b4+(b3-es-a1)*phi)*(cratio(m,pft)-0.05)
                            
              mgpp(m,pft,2)=co2(2) * (1-deltaa(pft)/1000.)-deltaa(pft)
              agpp(pft,2)=agpp(pft,2)+mgpp(m,pft,2)*mgpp(m,pft,1)
              mgpp(m,pft,3)=co2(3)-2*deltaa(pft)
              agpp(pft,3)=agpp(pft,3)+mgpp(m,pft,3)*mgpp(m,pft,1)
              c4gpp=c4gpp+mgpp(m,pft,1)

           else
                      
c             calculate the discrimination -> C3 pathway
              
              q = a*(1-cratio(m,pft)+0.025)
              r = 0.075*(es+a1)
              s = b*(cratio(m,pft)-0.1)
              t = (e*rdi/k+f*gamma)/catm

              
              deltaa(pft) = q+r+s-t
                       
              mgpp(m,pft,2)=co2(2) * (1 - deltaa(pft)/1000.)
     *                     -deltaa(pft)
              agpp(pft,2)=agpp(pft,2)+mgpp(m,pft,2)*mgpp(m,pft,1)
              mgpp(m,pft,3)=co2(3)-2*deltaa(pft)
              agpp(pft,3)=agpp(pft,3)+mgpp(m,pft,3)*mgpp(m,pft,1)
              c3gpp=c3gpp+mgpp(m,pft,1)

           endif
           
        else       
           mgpp(m,pft,2)=0.0
        endif
      enddo                     !month


      agpp(pft,2)=agpp(pft,2)/(c3gpp+c4gpp)
      agpp(pft,3)=agpp(pft,3)/(c3gpp+c4gpp)
      
      
      return
      end


c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE CALCPHI
c     This subroutine is for calculating the phi variable used in 
c     C4 photosynethsis isotope discrimination calculations

      subroutine calcphi(mgpp,pft,phi)

      implicit none
	  
c	  PARAMETERS
      INTEGER npft
        PARAMETER(npft=13)	  

c     ARGUMENTS
      real mgpp(1:12,1:npft,1:3),phi
      integer pft


C     LOCAL
      real snormavg(1:4),svar(1:4)
      real avar,a,gppm(1:12)
      real totgpp,meangpp,normgpp(1:12)
      integer m,s

      snormavg(:)=0.0
      svar(:)=0.0
      avar=0.0
      a=0.0
      gppm(:)=0.0
      totgpp=0.0
      meangpp=0.0
      normgpp(:)=0.0
      m=0
      s=0

      do m=1,12
         gppm(m)=mgpp(m,pft,1)
      enddo
      totgpp=0.0       !initialize a few variables
      do s=1,4
        svar(s)=0.0
      enddo

c     This first part of the subroutine estimates annual variability of
c     GPP first by normalizing and then summing seasonal variability
c     which compensates for amplitude and seasonal variation in GPP.

      do m=1,12
        totgpp=totgpp+gppm(m)
      enddo

      meangpp=totgpp/12.0

      do m=1,12
        normgpp(m)=gppm(m)/meangpp
      enddo 

      snormavg(1)=(normgpp(12)+normgpp(1)+normgpp(2))/3.0
      snormavg(2)=(normgpp(3)+normgpp(4)+normgpp(5))/3.0
      snormavg(3)=(normgpp(6)+normgpp(7)+normgpp(8))/3.0
      snormavg(4)=(normgpp(9)+normgpp(10)+normgpp(11))/3.0

c     calculate the population variances by season

      do m=1,2
        a=((normgpp(m)-snormavg(1))**2)/3
        svar(1)=svar(1)+a
      enddo
      svar(1)=svar(1)+(((normgpp(12)-snormavg(1))**2)/3)

      do m=3,5
        a=((normgpp(m)-snormavg(2))**2)/3
        svar(2)=svar(2)+a
      enddo

      do m=6,8
        a=((normgpp(m)-snormavg(3))**2)/3
        svar(3)=svar(3)+a
      enddo

      do m=9,11
        a=((normgpp(m)-snormavg(4))**2)/3
        svar(4)=svar(4)+a
      enddo

      avar=svar(1)+svar(2)+svar(3)+svar(4)

c     ------------------------------------------------

c     This part sets the phi value based upon the annual variability.
c     The equation is a simple regresion based upon hypothetical extreme
c     scenarios of phi.

      phi=0.3518717*avar+0.2552359


      return

      end 



c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE SOILTEMP
c     Estimate soil temperatures given air temperatures for current
c     and previous year

      subroutine soiltemp(soilpar,mtemp,mtemp_old,mtemp_soil,mw1,mw1_t)

      implicit none

c     PARAMETERS:
      real pi
        parameter (pi=3.14159265)

c     ARGUMENTS:
      real soilpar(1:7),mtemp(1:12),mtemp_old(1:12),mtemp_soil(1:12)
      real mw1(1:12),mw1_t(1:12)

c     LOCAL VARIABLES:
      real meanw1,diffus,mtempair(1:24)
      real tempthismonth,templastmonth,avetemp
      real alag,amp,lag,lagtemp
      integer m,n,startmonth

      meanw1=0.0
      diffus=0.0
      mtempair(:)=0.0
      tempthismonth=0.0
      templastmonth=0.0
      avetemp=0.0
      alag=0.0
      amp=0.0
      lag=0.0
      lagtemp=0.0
      m=0
      n=0
      startmonth=0

c     Annual cycle of soil temperatures follows surface temperatures with
c     damped oscillation about a common mean, and temporal lag.  Soil
c     temperature at depth z and time t given by (Carslaw & Jaeger 1959;
c     Eqn 5.52, Jury et al 1991):

c     T(z,t) = Tav + A*exp(-z/d)*sin(omega*t - z/d)

c     Tav       = average (base) air/soil temperature (avetemp)
c     A         = amplitude of air temp fluctuation
c     exp(-z/d) = fractional amplitude of temp fluctuation at soil depth z,
c                 relative to surface temp fluctuation (amp)
c     z/d       = oscillation lag in angular units at soil depth z (alag)
c     z         = soil depth = 0.25 m
c     d         = sqrt(2*K/omega), damping depth
c     K         = soil thermal diffusivity, m2/month (diffus)
c     omega     = 2*pi/tau, angular frequency of oscillation
c     tau       = oscillation period = 12 months

c     Assume a sinusoidal cycle, but estimate soil temperatures by
c     linear interpolation between previous monthly air temperatures to
c     implement a lag relative to air temperature, and damping of soil
c     temperature amplitude relative to air temperature amplitude.

c     Soil thermal diffusivities are calculated by linear interpolation
c     between estimates for 0, 15% and 100% soil water content from
c     van Duin (1963) and Jury et al (1991) Fig 5.11.

c     Calculate mean annual water content in soil layer 1

      meanw1=0.0
      do m=1,12
        meanw1=meanw1+mw1(m)/12.0
      enddo

c     In case of zero soil water, return with soil temp = air temp

      if (meanw1.eq.0.0) then
        do m=1,12
          mtemp_soil(m)=mtemp(m)
        enddo
        return
      endif

c     Interpolate thermal diffusivity function against soil water content

      if (meanw1.lt.0.15) then
        diffus=(soilpar(6)-soilpar(5))/0.15*meanw1+soilpar(5)
      else
        diffus=(soilpar(7)-soilpar(6))/0.85*(meanw1-0.15)+soilpar(6)
      endif

c     Convert diffusivity from mm2/s to m2/month
c     multiplication by 1e-6 (-> m2/s) * 2.628e6 (s/month) = 2.628

      diffus=diffus*2.628

c     Record 24 months of air temperatures in mtempair

      do m=1,12
        mtempair(m)=mtemp_old(m)
        mtempair(m+12)=mtemp(m)
      enddo

c     Calculate amplitude fraction and lag at soil depth 0.25 m

      alag=0.25/sqrt(12.0*diffus/pi)
      amp=exp(-alag)
      lag=alag*(6.0/pi)  !convert lag from angular units to months

c     Calculate monthly soil temperatures for this year.  For each month,
c     calculate average air temp for preceding 12 months (including this one)

      do m=1,12

        startmonth=m+1
        avetemp=0.0
        do n=startmonth,startmonth+11
          avetemp=avetemp+mtempair(n)/12.0
        enddo

c       Estimate air temperature "lag" months ago by linear interpolation
c       between air temperatures for this and last month

        tempthismonth=mtempair(m+12)
        templastmonth=mtempair(m+11)
        lagtemp=(tempthismonth-templastmonth)*(1.0-lag)+templastmonth

c       Adjust amplitude of lagged air temp to give estimated soil temp

        mtemp_soil(m)=avetemp+amp*(lagtemp-avetemp)

      enddo

      do m=1,12
         mtemp_old(m)=mtemp(m)
      enddo

      return
      end



c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE NPP
c     Calculation of maintenance and growth respiration and NPP

      subroutine npp(pftpar,dtemp,dtemp_soil,tree,dphen,nind,
     *  lm_ind,sm_ind,rm_ind,mgpp,anpp,mnpp,bm_inc,present,
     *  lresp,sresp,rresp,gresp,aresp,year,agpp,delt_c13_fpc,fpc_grid)

      implicit none

c     PARAMETERS
      integer npft,npftpar,nsoilpar
!        parameter (npft=13,npftpar=50,nsoilpar=7)
        parameter (npft=13,npftpar=59,nsoilpar=7)
      integer nco2
        parameter (nco2=3)           !number of C variables: C,13C,14C
      real k
        parameter (k=0.0548)

c     ARGUMENTS
      real pftpar(1:npft,1:npftpar)
      real dtemp(1:365),dtemp_soil(1:365)
      logical tree(1:npft)
      real dphen(1:365,1:npft)
      real nind(1:npft)
      real lm_ind(1:npft,1:nco2),sm_ind(1:npft,1:nco2)
      real rm_ind(1:npft,1:nco2)
      real bm_inc(1:npft,1:nco2)
      real mgpp(1:12,1:npft,1:nco2),agpp(1:npft,1:nco2)
      real anpp(1:npft,1:nco2),mnpp(1:12,1:npft,1:nco2)
      logical present(1:npft)
      real gresp(1:12,1:npft,1:nco2)   ! monthly growth respiration 
      real lresp(1:12,1:npft,1:nco2)  ! monthly leaf respiration
      real sresp(1:12,1:npft,1:nco2)  ! monthly sapwood respiration
      real rresp(1:12,1:npft,1:nco2)  ! monthly root respiration
      real aresp(1:12,1:npft,1:nco2)  ! monthly autotrophic respiration
      integer year


c     LOCAL VARIABLES
      integer pft
      real l_cton     ! leaf tissue C:N mass ratio
      real s_cton     ! sapwood tissue C:N mass ratio
      real r_cton     ! fineroot tissue C:N mass ratio
      real respcoeff  ! coefficient in respiration equations
      real lresp_m(1:nco2)    ! monthly total leaf respiration
      real sresp_m(1:nco2)    ! monthly total sapwood respiration
      real rresp_m(1:nco2)    ! monthly total root respiration
      real lresp_d    ! daily leaf respiration
      real sresp_d    ! daily leaf respiration
      real rresp_d    ! daily sapwood respiration
      real gtemp_air  ! value of temperature response function
                      !given air temperature
      real gtemp_soil ! value of temperature response function
                      !given soil temperature
      real mresp(1:nco2)      ! total monthly maintenance respiration
      real gresp_m(1:nco2)      ! monthly growth respiration
      real temp
      integer m,d,dm
      integer nc        ! nc is added for nco2

       real delt_c13_fpc,fpc_grid(1:npft)

      integer ndaymonth(1:12)    !number of days in each month
      data (ndaymonth(m),m=1,12)
     *  / 31,28,31,30,31,30,31,31,30,31,30,31 /

c     Calculation of maintenance respiration components for each
c     living tissue.

c     Based on the relations

c     (A) Tissue respiration response to temperature
c         (Sprugel et al, in press, Eqn 7)

c         (A1) Rm = 7.4e-7 * N * f(T)
c         (A2) f(T) = EXP (beta * T)

c           where Rm   = tissue maintenance respiration rate in mol C/sec
c                 N    = tissue nitrogen in mol N
c                 f(T) = temperature response function
c                 beta = ln Q10 / 10
c                 Q10  = change in respiration rate with a 10 K change
c                        in temperature
c                 T    = tissue absolute temperature in K

c     (B) Temperature response of soil respiration across ecosystems
c         incorporating damping of Q10 response due to temperature acclimation
c         (Lloyd & Taylor 1994, Eqn 11)

c         (B1) R = R10 * g(T)
c         (B2) g(T) = EXP [308.56 * (1 / 56.02 - 1 / (T - 227.13))]

c           where R    = respiration rate
c                 R10  = respiration rate at 10 K
c                 g(T) = temperature response function at 10 deg C
c                 T    = soil absolute temperature in K

c     Mathematical derivation:

c     For a tissue with C:N mass ratio cton, and C mass, c_mass, N concentration
c     in mol given by
c       (1) N = c_mass / cton / atomic_mass_N
c     Tissue respiration in gC/day given by
c       (2) R = Rm * atomic_mass_C * seconds_per_day
c     From (A1), (1) and (2),
c       (3) R = 7.4e-7 * c_mass / cton / atomic_mass_N * atomic_mass_C
c             * seconds_per_day * f(T)
c     Let
c       (4) k = 7.4e-7 * atomic_mass_C / atomic_mass_N * seconds_per_day
c             = 0.0548
c     from (3), (4)
c       (5) R = k * c_mass / cton * f(T)
c     substituting ecosystem temperature response function g(T) for f(T)
c     (Eqn B2),
c       (6) R = k * c_mass / cton * EXP [308.56 * (1 / 56.02 - 1 /
c             (T - 227.13))]
c     incorporate PFT-specific respiration coefficient to model acclimation
c     of respiration rates to average (temperature) conditions for PFT (Ryan
c     1991)
c       (7) R_pft = respcoeff_pft * k * c_mass / cton
c             * EXP [308.56 * (1 / 56.02 - 1 / (T - 227.13))]

c     Initializations (see comment at the beginning of allocation)
      pft=0
      l_cton=0.0
      s_cton=0.0
      r_cton=0.0
      respcoeff=0.0
      lresp_m(:)=0.0
      sresp_m(:)=0.0
      rresp_m(:)=0.0
      lresp_d=0.0
      sresp_d=0.0
      rresp_d=0.0
      gtemp_air=0.0
      gtemp_soil=0.0
      mresp(:)=0.0
      gresp_m(:)=0.0
      temp=0.0
      m=0
      d=0
      dm=0
      nc=0

        delt_c13_fpc=0.0  !delt_c13_fpc is the FPC weighted average delta-C13 of C3 and C4 grasses.

      do pft=1,npft
        if (present(pft)) then
          l_cton=pftpar(pft,13)
          s_cton=pftpar(pft,14)
          r_cton=pftpar(pft,15)

          respcoeff=pftpar(pft,5)

          anpp(pft,1)=0.0
          anpp(pft,2)=0.0
          anpp(pft,3)=0.0

          d=0  !day of year

          do m=1,12

           do nc=1,nco2
            lresp_m(nc)=0.0
            sresp_m(nc)=0.0
            rresp_m(nc)=0.0
           enddo



            do dm=1,ndaymonth(m)

              d=d+1  !next day of year

              if (dtemp(d).ge.-40.0) then

                gtemp_air=
     *            exp(308.56*(1.0/56.02-1.0/(dtemp(d)+46.02)))

                gtemp_soil=
     *            exp(308.56*(1.0/56.02-1.0/(dtemp_soil(d)+46.02)))

                    ! Eqn (7)
                    ! 46.02 = -227.13 + 273.15; conversion from deg C to K

              else

                gtemp_air=0.0
                gtemp_soil=0.0

              endif
c             Calculate tissue maintenance respiration values today [Eqn (7)]
c             incorporating daily phenology
              if (tree(pft)) then

                lresp_d=respcoeff*k*lm_ind(pft,1)*nind(pft)/
     *            l_cton*gtemp_air*dphen(d,pft)
                sresp_d=respcoeff*k*sm_ind(pft,1)*nind(pft)/
     *            s_cton*gtemp_air
                rresp_d=respcoeff*k*rm_ind(pft,1)*nind(pft)/
     *            r_cton*gtemp_soil*dphen(d,pft)

                lresp_m(1)=lresp_m(1)+lresp_d
                sresp_m(1)=sresp_m(1)+sresp_d
                rresp_m(1)=rresp_m(1)+rresp_d

              else   !grass

                lresp_d=respcoeff*k*lm_ind(pft,1)*nind(pft)/
     *            l_cton*gtemp_air*dphen(d,pft)
                rresp_d=respcoeff*k*rm_ind(pft,1)*nind(pft)/
     *            r_cton*gtemp_soil*dphen(d,pft)

                lresp_m(1)=lresp_m(1)+lresp_d
                rresp_m(1)=rresp_m(1)+rresp_d

              endif
            enddo  !day of month

c           Incorporate growth respiration = 25% of
c           (GPP - maintenance respiration) in monthly NPP calculation

            if (tree(pft)) then
              mresp(1)=lresp_m(1)+sresp_m(1)+rresp_m(1)
              do nc=2,nco2
                mresp(nc)=mgpp(m,pft,nc) ! assuming that assimilate is directly used for autotr. resp.
                lresp_m(nc)=lm_ind(pft,nc)
                sresp_m(nc)=sm_ind(pft,nc)
                rresp_m(nc)=rm_ind(pft,nc)
              enddo 
            else !grass
              mresp(1)=lresp_m(1)+rresp_m(1)
              do nc=2,nco2
                mresp(nc)=mgpp(m,pft,nc) ! assuming that assimilate is directly used for autotr. resp.
                lresp_m(nc)=lm_ind(pft,nc)
                rresp_m(nc)=rm_ind(pft,nc)  
              enddo
            endif
            gresp_m(1)=max((mgpp(m,pft,1)-mresp(1))*0.25,0.0)
                    !growth resp always non-negativ
            gresp_m(2)=mgpp(m,pft,2)
            gresp_m(3)=mgpp(m,pft,3)

            mnpp(m,pft,1)=mgpp(m,pft,1)-mresp(1)-gresp_m(1) 

            if (mnpp(m,pft,1).gt.1.E-9) then  
               mnpp(m,pft,2)=mgpp(m,pft,2)       
               mnpp(m,pft,3)=mgpp(m,pft,3)
            else
               mnpp(m,pft,2)=0.0   !PRoblem: mnpp can be <0, see also bm_inc
               mnpp(m,pft,3)=0.0
            endif

c           Increment annual NPP
            temp=anpp(pft,1)
            anpp(pft,1)=anpp(pft,1)+mnpp(m,pft,1)
            if (mnpp(m,pft,1).gt.0.0) then
               if(anpp(pft,1).gt.0.0.and.temp.le.0.0) then
                  do nc=2,nco2
                     anpp(pft,nc)=mnpp(m,pft,nc)
                  enddo
               elseif(anpp(pft,1).gt.0.0.and.temp.gt.0.0) then
                  do nc=2,nco2
                     anpp(pft,nc)=(temp*anpp(pft,nc)+mnpp(m,pft,1)
     *                    *mnpp(m,pft,nc))/anpp(pft,1)
                  enddo
               endif
            endif
            
c           Copy NPP to bm_inc, which will store running total of
c           C production remaining for growth            
            do nc=1,nco2
              lresp(m,pft,nc)=lresp_m(nc)
              sresp(m,pft,nc)=sresp_m(nc)
              rresp(m,pft,nc)=rresp_m(nc)
              gresp(m,pft,nc)=gresp_m(nc)
            enddo
            aresp(m,pft,1)=lresp_m(1)+sresp_m(1)+rresp_m(1)+gresp_m(1)
            if (aresp(m,pft,1).gt.0.) then
              do nc=2,nco2
                aresp(m,pft,nc)=(lresp_m(1)*lresp_m(nc)+sresp_m(1)
     *               *sresp_m(nc)+rresp_m(1)*rresp_m(nc)+gresp_m(1)
     *               *gresp_m(nc))/aresp(m,pft,1)
              enddo
            endif

          enddo  !month     

            do nc=1,nco2
              bm_inc(pft,nc)=anpp(pft,nc)
            enddo

c         don't allow negative anpp. set to zero

          if (anpp(pft,1).lt.0.0) then
             bm_inc(pft,2)=0.0
             bm_inc(pft,3)=0.0
             anpp(pft,1)=0.0
             anpp(pft,2)=0.0
             anpp(pft,3)=0.0
             do m=1,12  
               do nc=1,nco2
                mnpp(m,pft,nc)=0.0
               enddo
             enddo  ! month
          endif 

        if(.not.tree(pft))then
!!      delt_C13_pft=delt_C13_pft+anpp(pft,2)*fpc_grid(pft)*anpp(pft,1)
          delt_C13_fpc=delt_C13_fpc+anpp(pft,2)*fpc_grid(pft)
         end if
        endif  !present(pft)

      enddo !pft

C     Doug 02/13: grass pfts hardcoded. This needs redoing
      if(anpp(8,2).ne.0.0.and.anpp(9,2).ne.0.0)then
             delt_c13_fpc=delt_c13_fpc/(fpc_grid(8)+fpc_grid(9))
      else if(anpp(8,2).eq.0.0.and.anpp(9,2).ne.0.0)then
         delt_c13_fpc=delt_c13_fpc/fpc_grid(9)
      else if(anpp(8,2).eq.0.0.and.anpp(9,2).ne.0.0)then
         delt_c13_fpc=delt_c13_fpc/fpc_grid(8)
      else
         delt_c13_fpc=0.0
      end if

      return
      end



c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE REPRODUCTION
c     Deduction of reproduction costs from annual biomass increment
c
c     Doug 01/09: addition of fuel_1hr_inc as I/O

      subroutine reproduction(bm_inc,lm_sapl,sm_sapl,hm_sapl,rm_sapl,
     *     litter_ag_leaf, !Doug 11/12
     *     fuel_1hr_leaf,
     *     fuel_1hr_leaf_inc_pos, fuel_1hr_leaf_inc_neg,
     *     present, tree, co2)

      implicit none

c     PARAMETERS
      integer npft,npftpar,nsoilpar
!      parameter (npft=13,npftpar=50,nsoilpar=7)
      parameter (npft=13,npftpar=59,nsoilpar=7)
      integer nco2
        parameter (nco2=3)           !number of C variables: C,13C,14C
      real reprod_cost              !proportion of NPP lost to reproduction
        parameter (reprod_cost=0.1) !Harper 1977

c     ARGUMENTS
      real bm_inc(1:npft,1:nco2)
      REAL litter_ag_leaf(1:npft,1:nco2) !Doug 11/12
      real fuel_1hr_leaf(1:npft,1:nco2)
      REAL fuel_1hr_leaf_inc_pos(1:npft,1:nco2,1:12)
      REAL fuel_1hr_leaf_inc_neg(1:npft,1:nco2,1:12)	!Doug 01/09: records variations
      REAL temp_fuel(1:nco2)		        !in fuel for each month (for FIRE)
      logical present(1:npft)
      real lm_sapl(1:npft,1:nco2),sm_sapl(1:npft,1:nco2)
      real hm_sapl(1:npft,1:nco2),rm_sapl(1:npft,1:nco2)
      real co2(1:nco2)
      logical tree(1:npft)

c     LOCAL VARIABLES
      integer pft
      integer nc
      real temp
      real reprod        !allocation to reproduction (gC/m2)

      pft=0
      nc=0
      temp=0.0
      reprod=0.0
      
      do pft=1,npft	!pft

        temp_fuel=fuel_1hr_leaf(pft,:)	!Doug 01/09: MFL

        if (present(pft)) then

c         Calculate allocation to reproduction
c         Reproduction costs taken simply as a constant fraction of annual NPP


          reprod=max(bm_inc(pft,1)*reprod_cost,0.0)



c         assume the costs go to reproductive structures which will
c         eventually enter the litter pool

          temp=litter_ag_leaf(pft,1) !Doug 11/12: litter_ag -> litter_ag_leaf. rerpod goes into grass and leaf litter
          litter_ag_leaf(pft,1)=litter_ag_leaf(pft,1)+reprod
          if (reprod.gt.0.0) then
             do nc=2,nco2
                litter_ag_leaf(pft,nc)=(litter_ag_leaf(pft,nc)*temp+
     *               bm_inc(pft,nc)*reprod)/(litter_ag_leaf(pft,1))
             enddo
          endif

 
          fuel_1hr_leaf(pft,1)=fuel_1hr_leaf(pft,1)+reprod

          if (reprod.gt.0.0) then
             do nc=2,nco2
          temp=fuel_1hr_leaf(pft,nc)
                fuel_1hr_leaf(pft,nc)=(fuel_1hr_leaf(pft,nc)*temp+
     *               bm_inc(pft,nc)*reprod)/(fuel_1hr_leaf(pft,1))
             enddo
          endif
     
c         Reduce biomass increment by reproductive cost

          bm_inc(pft,1)=bm_inc(pft,1)-reprod

          if (tree(pft)) then
            if (bm_inc(pft,1).gt.0.0) then
               do nc=2,nco2
                  lm_sapl(pft,nc)=bm_inc(pft,nc)
                  sm_sapl(pft,nc)=bm_inc(pft,nc)
                  hm_sapl(pft,nc)=bm_inc(pft,nc)
                  rm_sapl(pft,nc)=bm_inc(pft,nc)
               enddo
            else              
               lm_sapl(pft,2)=co2(2)-17.8 !C3 value         
               sm_sapl(pft,2)=co2(2)-17.8 !from lloyd & farquhar,1994
               hm_sapl(pft,2)=co2(2)-17.8
               rm_sapl(pft,2)=co2(2)-17.8
               lm_sapl(pft,3)=co2(3)-35.6
               sm_sapl(pft,3)=co2(3)-35.6
               hm_sapl(pft,3)=co2(3)-35.6
               rm_sapl(pft,3)=co2(3)-35.6
            endif
          else
            if (bm_inc(pft,1).gt.0.0) then
               do nc=2,nco2
                  lm_sapl(pft,nc)=bm_inc(pft,nc)
                  rm_sapl(pft,nc)=bm_inc(pft,nc)
               enddo
            else
               if (pft.eq.npft) then
                  lm_sapl(pft,2)=co2(2)-3.6 !C4 value
                  rm_sapl(pft,2)=co2(2)-3.6 !from lloyd & farquhar,1994 
                  lm_sapl(pft,3)=co2(3)-7.2
                  rm_sapl(pft,3)=co2(3)-7.2
               else
                  lm_sapl(pft,2)=co2(2)-17.8 !C3 value
                  rm_sapl(pft,2)=co2(2)-17.8 !from lloyd & farquhar,1994 
                  lm_sapl(pft,3)=co2(3)-35.6
                  rm_sapl(pft,3)=co2(3)-35.6
               endif
            endif
          endif

        endif
        DO nc=1,nco2
          IF (fuel_1hr_leaf(pft,1)-temp_fuel(1)>0) THEN                                    !Doug 01/09: MFL
            fuel_1hr_leaf_inc_pos(pft,nc,:)=
     *        fuel_1hr_leaf_inc_pos(pft,nc,:)+
     *        (fuel_1hr_leaf(pft,nc)-temp_fuel(nc))/12
          ELSE
            fuel_1hr_leaf_inc_neg(pft,nc,:)=
     *        fuel_1hr_leaf_inc_neg(pft,nc,:)+
     *        (fuel_1hr_leaf(pft,nc)-temp_fuel(nc))/12
          END IF
        ENDDO                                                          !Doug 01/09: MFL 

      enddo	!pft

      

      return
      end	!reproduction



c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE TURNOVER
c     Turnover of PFT-specific fraction from each living C pool
c     Leaf and root C transferred to litter, sapwood C to heartwood
c
c     Doug 01/09: addition of fuel_ihr_inc as I/O

      subroutine turnover(pftpar,present,tree,lm_ind,sm_ind, 
     *  hm_ind,rm_ind,
     *  litter_ag_leaf,litter_bg,fuel_1hr_leaf,
     *  fuel_10hr,fuel_100hr,
     *  fuel_1hr_leaf_inc_pos, fuel_1hr_leaf_inc_neg,
     *  fuel_10hr_inc,
     *  fuel_100hr_inc, nind, turnover_ind)

      implicit none

c     PARAMETERS:
      integer npft,npftpar,nsoilpar
!        parameter (npft=13,npftpar=50,nsoilpar=7)
        parameter (npft=13,npftpar=59,nsoilpar=7)
      integer nco2
        parameter (nco2=3)           !number of C variables: C,13C,14C

c     ARGUMENTS:
      real pftpar(1:npft,1:npftpar)
      logical present(1:npft),tree(1:npft)
      real lm_ind(1:npft,1:nco2)
      real sm_ind(1:npft,1:nco2)
      real hm_ind(1:npft,1:nco2)
      real rm_ind(1:npft,1:nco2)
      REAL litter_ag_leaf(1:npft,1:nco2) !Doug 11/12
      real litter_bg(1:npft,1:nco2)
      REAL fuel_1hr_leaf(1:npft,1:nco2)
      REAL fuel_10hr(1:npft,1:nco2)
      real fuel_100hr(1:npft,1:nco2)
      REAL fuel_1hr_leaf_inc_pos(1:npft,1:nco2,1:12)
      REAL fuel_1hr_leaf_inc_neg(1:npft,1:nco2,1:12)	!Doug 01/09: MFL; records variations
      REAL fuel_10hr_inc(1:npft,1:nco2,1:12)	!in fuel for each month (for FIRE)
      REAL fuel_100hr_inc(1:npft,1:nco2,1:12)
      REAL temp_fuel(1:nco2,1:3)
      real nind(1:npft)
      real turnover_ind(1:npft)

c     LOCAL VARIABLES:
      integer pft
      integer nc
      real l_torate
      real s_torate
      real r_torate
      real lm_turn
      real sm_turn
      real rm_turn
      real temp

c     ALLAN
       real small_branch_torate,large_twig_torate,small_twig_torate
       parameter (small_branch_torate=0.02,large_twig_torate=0.0333,
     *    small_twig_torate=0.05)

       real small_branch_drop
       real large_twig_drop
       real small_twig_drop
       real sm_turn_drop
       real hm_turn_drop

      pft=0
      nc=0
      l_torate=0.0
      s_torate=0.0
      r_torate=0.0
      lm_turn=0.0
      sm_turn=0.0
      rm_turn=0.0
      temp=0.0
       small_branch_drop=0.0
       large_twig_drop=0.0
       small_twig_drop=0.0
       sm_turn_drop=0.0
       hm_turn_drop=0.0

      do pft=1,npft

        temp_fuel(:,1)=fuel_1hr_leaf(pft,:)		!Doug 01/09: MFL
        temp_fuel(:,2)=fuel_10hr(pft,:)		!Doug 01/09: MFL
        temp_fuel(:,3)=fuel_100hr(pft,:)	!Doug 01/09: MFL

c       Turnover rates are reciprocals of tissue longevity

        l_torate=1.0/pftpar(pft,9)
        s_torate=1.0/pftpar(pft,11)
        r_torate=1.0/pftpar(pft,12)

        if (present(pft)) then

c         Calculate the biomass turnover in this year

          lm_turn=lm_ind(pft,1)*l_torate
          sm_turn=sm_ind(pft,1)*s_torate
          rm_turn=rm_ind(pft,1)*r_torate

c         Update the pools

          lm_ind(pft,1)=lm_ind(pft,1)-lm_turn
          sm_ind(pft,1)=sm_ind(pft,1)-sm_turn
          rm_ind(pft,1)=rm_ind(pft,1)-rm_turn

c         Convert sapwood to heartwood

          temp=hm_ind(pft,1)
          hm_ind(pft,1)=hm_ind(pft,1)+sm_turn
          if (hm_ind(pft,1).gt.0.0) then
             do nc=2,nco2
                hm_ind(pft,nc)=(hm_ind(pft,nc)*temp+sm_ind(pft,nc)
     *               *sm_turn)/hm_ind(pft,1)
             enddo
          endif

c      Drop small branches, large & small twigs from
c      sapwood, and then from heartwood mass. Update S & H pools.

c      discussion 25/7/05: revise ratios of sm for branch & twig drop. No coarse roots in LPJ      
c       small_branch_drop = (0.1422*sm_ind(pft))*small_branch_torate
c       large_twig_drop = (0.05*sm_ind(pft))*large_twig_torate
c       small_twig_drop = (0.03*sm_ind(pft))*small_twig_torate

c       sm_turn_drop=small_branch_drop +
c     *   large_twig_drop + small_twig_drop

c       sm_ind(pft)=sm_ind(pft)-sm_turn_drop
        
c      discussion 25/7/05: revise ratios of sm for branch & twig drop. No coarse roots in LPJ      

       small_branch_drop = (0.1422*hm_ind(pft,1))*small_branch_torate
       large_twig_drop = (0.05*hm_ind(pft,1))*large_twig_torate
       small_twig_drop = (0.03*hm_ind(pft,1))*small_twig_torate

c       hm_turn_drop=small_branch_drop +
c     *  large_twig_drop + small_twig_drop
  
c       hm_ind(pft)=hm_ind(pft)-hm_turn_drop

c         Transfer to litter pools

       temp=litter_ag_leaf(pft,1)
       litter_ag_leaf(pft,1)=litter_ag_leaf(pft,1)+lm_turn*nind(pft)
       if (litter_ag_leaf(pft,1).gt.0.0) then
          do nc=2,nco2
             litter_ag_leaf(pft,nc)=(litter_ag_leaf(pft,nc)*temp
     *            +lm_ind(pft,nc)*lm_turn*nind(pft))/
     *            litter_ag_leaf(pft,1)
          enddo
       endif
          
       temp=litter_bg(pft,1)
       litter_bg(pft,1)=litter_bg(pft,1)+rm_turn*nind(pft)
       if (litter_bg(pft,1).gt.0.0) then
          do nc=2,nco2
             litter_bg(pft,nc)=(litter_bg(pft,nc)*temp+rm_ind(pft,nc)
     *            *rm_turn*nind(pft))/litter_bg(pft,1)
          enddo
       endif

c       Allan
c       litter_ag(pft)=litter_ag(pft)+sm_turn_drop*nind(pft)
c       litter_ag(pft)=litter_ag(pft)+hm_turn_drop*nind(pft)


c         KIRSTEN: fuel classes

       temp=fuel_1hr_leaf(pft,1)
       fuel_1hr_leaf(pft,1)=fuel_1hr_leaf(pft,1)+lm_turn*nind(pft)
       if (fuel_1hr_leaf(pft,1).gt.0.0) then
          do nc=2,nco2
             fuel_1hr_leaf(pft,nc)=(fuel_1hr_leaf(pft,nc)*temp
     *            +lm_ind(pft,nc)*lm_turn*nind(pft))/
     *            fuel_1hr_leaf(pft,1)
          enddo
       endif

c       if (tree(pft)) then
c       fuel_1hr(pft)=fuel_1hr(pft)+ small_twig_drop*nind(pft) !A
c       fuel_10hr(pft)=fuel_10hr(pft)+large_twig_drop*nind(pft) !A
c       fuel_100hr(pft)=fuel_100hr(pft)+small_branch_drop*nind(pft) !A
c       endif
 
c         Record total turnover

          turnover_ind(pft)=lm_turn+sm_turn+rm_turn
c       Allan
c       turnover_ind(pft)=turnover_ind(pft)+sm_turn_drop+hm_turn_drop

          
          do nc=1,nco2
             if (lm_ind(pft,1).le.0.) lm_ind(pft,nc)=0. 
             if (sm_ind(pft,1).le.0.) sm_ind(pft,nc)=0. 
             if (rm_ind(pft,1).le.0.) rm_ind(pft,nc)=0.
          enddo


        endif

        DO nc=1,nco2
          IF (fuel_1hr_leaf(pft,1)-temp_fuel(1,1)>0) THEN						!Doug 01/09: MFL
            fuel_1hr_leaf_inc_pos(pft,nc,:)=
     *         fuel_1hr_leaf_inc_pos(pft,nc,:)+
     *        (fuel_1hr_leaf(pft,nc)-temp_fuel(nc,1))/12
          ELSE
            fuel_1hr_leaf_inc_neg(pft,nc,:)=
     *         fuel_1hr_leaf_inc_neg(pft,nc,:)+
     *        (fuel_1hr_leaf(pft,nc)-temp_fuel(nc,1))/12
          END IF
          fuel_10hr_inc(pft,nc,:)=fuel_10hr_inc(pft,nc,:)+
     *      (fuel_10hr(pft,nc)-temp_fuel(nc,2))/12
          fuel_100hr_inc(pft,nc,:)=fuel_100hr_inc(pft,nc,:)+
     *      (fuel_100hr(pft,nc)-temp_fuel(nc,3))/12
        END DO							!Doug 01/09: MFL

      enddo  !pft

      return
      end 



c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE LITTERSOM		
c     Litter and soil decomposition
c     Incorporates analytical solution for soil pool sizes once litter inputs
c     are (assumed to be) at equilibrium, reducing spin-up time for carbon
c     fluxes due to soil respiration.
c
c     Doug 01/09: addition of fuel_1hr_inc as I/O

      subroutine littersom(pftpar,litter_ag_leaf,litter_ag_wood, !Doug 11/12
     *  litter_bg,fuel_1hr_leaf,fuel_1hr_wood,fuel_10hr,
     *  fuel_100hr,fuel_1000hr,
     *  fuel_1hr_leaf_inc_pos, fuel_1hr_leaf_inc_neg,
     *  fuel_1hr_wood_inc_pos, fuel_1hr_wood_inc_neg,
     *  fuel_10hr_inc,
     *  fuel_100hr_inc, fuel_1000hr_inc,
     *  mw1,mw1_t,mtemp_soil,cpool_fast,
     *  cpool_slow,arh,mrh,year,k_fast_ave,k_slow_ave,litter_decom_ave,
     *  agri_litter_ag,agri_litter_bg,agri)    !Doug 05/09: landuse variable added

      implicit none

c     PARAMETERS
      integer npft,npftpar,nsoilpar
!      parameter (npft=13,npftpar=50,nsoilpar=7)
      parameter (npft=13,npftpar=59,nsoilpar=7)
      integer nco2
        parameter (nco2=3)           !number of C variables: C,13C,14C
c     DOug 12/09: No longer need the equlibium fix. Infact, the fix was a little bit rubbish.
c      integer soil_equil_year           !number of years until pool sizes for
c        parameter (soil_equil_year=400) !soil decomposition solved analytically
      real k_litter10                   !litter decomposition rate at
c     k_litter10 changed from 0.5 on 6/2000
        parameter (k_litter10=0.35)     !10 deg C (/year)
      real k_soil_fast10                !fast pool decomposition rate at
        parameter (k_soil_fast10=0.03)  !10 deg C (/year)
      real k_soil_slow10                !slow pool decomposition rate at
        parameter (k_soil_slow10=0.001) !10 deg C (/year)
      real fastfrac                     !fraction of litter entering fast soil
        parameter (fastfrac=0.985)      !decomposition pool
      real atmfrac                      !fraction of litter decomposition
        parameter (atmfrac=0.7)         !going directly into the atmosphere
      INTEGER nmonth                    !Doug 06/09: no of months per year
        PARAMETER (nmonth=12)

c     ARGUMENTS
      REAL pftpar(1:npft,1:npftpar)      ! Doug 11/12
      REAL litter_ag_leaf(1:npft,1:nco2) !Doug 11/12
      REAL litter_ag_wood(1:npft,1:nco2)
      real litter_bg(1:npft,1:nco2)
      REAL fuel_1hr_leaf(1:npft,1:nco2)
      REAL fuel_1hr_wood(1:npft,1:nco2)
      REAL fuel_10hr(1:npft,1:nco2)
      real fuel_100hr(1:npft,1:nco2)
      real fuel_1000hr(1:npft,1:nco2)
      REAL fuel_1hr_leaf_inc_pos(1:npft,1:nco2,1:12)
      REAL fuel_1hr_leaf_inc_neg(1:npft,1:nco2,1:12)	!Doug 01/09: records variations
      REAL fuel_1hr_wood_inc_pos(1:npft,1:nco2,1:12)
      REAL fuel_1hr_wood_inc_neg(1:npft,1:nco2,1:12)	!Doug 01/09: records variations
      REAL fuel_10hr_inc(1:npft,1:nco2,1:12)	!in fuel for each month (for FIRE)
      REAL fuel_100hr_inc(1:npft,1:nco2,1:12)
      REAL fuel_1000hr_inc(1:npft,1:nco2,1:12)
      REAL temp_fuel(1:npft,1:nco2,1:5)
      real mw1(1:12),mw1_t(1:12)
      real mtemp_soil(1:12)
      real cpool_fast(1:nco2)
      real cpool_slow(1:nco2)
      real arh(1:nco2)
      real mrh(1:12,1:nco2)
      integer year
      real k_fast_ave
      real k_slow_ave
      real litter_decom_ave(1:nco2)
c    Doug 05/09: landuse change variables
      LOGICAL agri           !Doug 06/09: if land is currently in agrculure
      REAL agri_litter_ag    !Doug 06/09: agricultureal above ground litter
      REAL agri_litter_bg    !Doug 06/09: agricultureal below ground litter

c     LOCAL VARIABLES
      integer pft,m
      integer nc       ! nc is added for nco2
      REAL k_litter_bg         !monthly litter decomposition rate (/month) !Doug 11/12: changed to below ground only
      REAL k_litter_leaf(1:npft)!Doug 11/12: monthly litter decomposition rate (/month) !Doug 11/12: changed to below ground only
      REAL k_litter_wood(1:npft)!Doug 11/12: monthly litter decomposition rate (/month) !Doug 11/12: changed to below ground only
      real k_fast(1:12)        !monthly fast pool decomposition rate (/month)
      real k_slow(1:12)        !monthly slow pool decomposition rate (/month)
      real temp_resp           !monthly temperature response of decomposition
      real moist_resp          !monthly moisture response of decomposition
      REAL leaf_repsonse, wood_response  !Doug 11/12
      real cflux_litter_soil   !litter decomposition flux to soil
      real cflux_litter_atmos  !litter decomposition flux to atmosphere
      real cflux_fast_atmos    !soil fast pool decomposition flux to atmosphere
      real cflux_slow_atmos    !soil slow pool decomposition flux to atmosphere
      real litter_decom(1:12,1:nco2)  !monthly litter decomposition
      REAL litter_decom_ag_leaf     !above-ground component of litter decomposition
      REAL litter_decom_ag_wood     !above-ground component of litter decomposition
      REAL fuel_decom_1hr_leaf(1:npft)   !fuel decomposition, parallel to litter_decom_ag
      REAL fuel_decom_1hr_wood(1:npft)   !fuel decomposition, parallel to litter_decom_ag
      real fuel_decom_10hr(1:npft)
      real fuel_decom_100hr(1:npft)
      real fuel_decom_1000hr(1:npft)
      real fuel_decom_slow     !fuel decomposition, parallel to soil fast decomp.

      real litter_decom_bg     !below-ground component of litter decomposition
      real soilfrac     !fraction of litter decomposition going to soil C pools
      real temp

c    Doug 05/09: landuse ones:
      REAL agri_decom(1:nmonth)    !Doug 06/09: decomposition of land used for agriculture
      REAL nat_decom(1:nmonth,1:nco2)
                                   !Doug 06/09: decomposition of natural (or noval) land.
      REAL agri_decom_ag           !Doug 06/09: decomposition of above ground agricultural litter
      REAL agri_decom_bg           !Doug 06/09: decomposition of below ground agricultural litter

      pft=0
      m=0
      nc=0
      k_litter_bg=0.0
      k_litter_leaf(:)=0.0
      k_litter_wood(:)=0.0
      k_fast(:)=0.0
      k_slow(:)=0.0
      temp_resp=0.0
      moist_resp=0.0
      cflux_litter_soil=0.0
      cflux_litter_atmos=0.0
      cflux_fast_atmos=0.0
      cflux_slow_atmos=0.0
      litter_decom(:,:)=0.0
      litter_decom_ag_leaf=0.0
      litter_decom_ag_wood=0.0
      fuel_decom_1hr_leaf(:)=0.0
      fuel_decom_1hr_wood(:)=0.0
      fuel_decom_10hr(:)=0.0
      fuel_decom_100hr(:)=0.0
      fuel_decom_1000hr(:)=0.0
      fuel_decom_slow=0.0
      litter_decom_bg=0.0
      soilfrac=0.0
      temp=0.0

      soilfrac=1.0-atmfrac

c DM   arh seemed to be never initialised, I think it is the right place
      arh(:)=0.0



      do m=1,12

        temp_fuel(:,:,1)=fuel_1hr_leaf	!Doug 01/09: MFL
        temp_fuel(:,:,2)=fuel_1hr_wood	!Doug 01/09: MFL
        temp_fuel(:,:,3)=fuel_10hr		!Doug 01/09: MFL
        temp_fuel(:,:,4)=fuel_100hr		!Doug 01/09: MFL
        temp_fuel(:,:,5)=fuel_1000hr	!Doug 01/09: MFL




c       Temperature response function is a modified Q10 relationship
c       (Lloyd & Taylor 1994)

        if (mtemp_soil(m).le.-40.0) then !avoid division by zero
           temp_resp=0.0
        else
            ! Doug 06/13: temp respinse for leaf and below ground litter and carbon pools.
           temp_resp=exp(308.56*((1.0/56.02)-
     *       (1.0/(mtemp_soil(m)+273.0-227.13))))  !Lloyd & Taylor 1994
        endif

c       Moisture response based on soil layer 1 moisture content (Foley 1995)
c        moist_resp=(0.25+(0.75*mw1(m)))
c        moist_resp=(0.2+(0.8*mw1(m)))
        moist_resp=(1.0-exp(-1.0*mw1(m)))/(1.0-exp(-1.0))

c       Calculate monthly decomposition rates (k, /month) as a function of
c       temperature and moisture

c       k = k_10 * temp_resp * moist_resp
C		Doug 11/12: seperate out litter and wood responses to temperature using Q10 parameters (pftpar 58)
        leaf_repsonse=temp_resp*moist_resp !Doug 11/12
        IF (mtemp_soil(m)>-40.) THEN
          wood_response=pftpar(pft,58)*(mtemp_soil(m)-10.)/10.
        ELSE
          wood_response=0.0
        END IF
        DO pft=1,npft

           k_litter_leaf(pft)=(pftpar(pft,56)*leaf_repsonse)/12.0
           k_litter_wood(pft)=(pftpar(pft,57)*wood_response)/12.0
        END DO
        k_litter_bg=(k_litter10*temp_resp*moist_resp)/12.0
        k_fast(m)=(k_soil_fast10*temp_resp*moist_resp)/12.0
        k_slow(m)=(k_soil_slow10*temp_resp*moist_resp)/12.0

c       Calculate monthly litter decomposition using equation
c         (1) dc/dt = -kc     where c=pool size, t=time, k=decomposition rate
c       from (1),
c         (2) c = c0*exp(-kt) where c0=initial pool size
c       from (2), decomposition in any month given by
c         (3) delta_c = c0 - c0*exp(-k)
c       from (4)
c         (4) delta_c = c0*(1.0-exp(-k))

C Doug 05/09: initialise landuse variables
        litter_decom(m,1)=0.0
        nat_decom(m,:)=0.0
        agri_decom(m)=0.0
        
        do pft=1,npft
          fuel_decom_1hr_leaf(pft)=0.0
          fuel_decom_1hr_wood(pft)=0.0
          fuel_decom_10hr(pft)=0.0
          fuel_decom_100hr(pft)=0.0
          fuel_decom_1000hr(pft)=0.0
          
          litter_decom_ag_leaf=litter_ag_leaf(pft,1)*
     *       (1.0-exp(-k_litter_leaf(pft)))  !eqn 4
	      litter_decom_ag_wood=litter_ag_wood(pft,1)*
     *       (1.0-exp(-k_litter_wood(pft)))  !eqn 4
          litter_decom_bg=litter_bg(pft,1)*(1.0-exp(-k_litter_bg))

c   Doug 05/09: natural decomposion
          nat_decom(m,1)=nat_decom(m,1)+litter_decom_ag_leaf+
     *      litter_decom_ag_wood+litter_decom_bg

          DO nc=2,nco2
             nat_decom(m,nc)=nat_decom(m,nc)+litter_ag_leaf(pft,nc)
     *         *litter_decom_ag_leaf+litter_ag_wood(pft,nc)
     *         *litter_decom_ag_wood+litter_bg(pft,nc)*litter_decom_bg  
          END DO !StopC

c          litter_decom(m,1)=litter_decom(m,1)+litter_decom_ag+
c     *      litter_decom_bg
c
c          do nc=2,nco2
c             litter_decom(m,nc)=litter_decom(m,nc)+litter_ag(pft,nc)
c     *         *litter_decom_ag+litter_bg(pft,nc)*litter_decom_bg  
c          enddo

c         Update the litter pools

          litter_ag_leaf(pft,1)=litter_ag_leaf(pft,1)-
     *      litter_decom_ag_leaf
          litter_ag_wood(pft,1)=litter_ag_wood(pft,1)-
     *      litter_decom_ag_wood
C          litter_ag(pft,1)=litter_ag(pft,1)-litter_decom_ag
          litter_bg(pft,1)=litter_bg(pft,1)-litter_decom_bg

c       KIRSTEN: calculate decomposition for 1hr fuel class, parallel to litter_decom_ag
c       no influence on carbon balance/fluxes


          fuel_decom_1hr_leaf(pft)=fuel_1hr_leaf(pft,1)*
     *      (1.0-exp(-k_litter_leaf(pft)))
          fuel_decom_1hr_wood(pft)=fuel_1hr_wood(pft,1)*
     *      (1.0-exp(-k_litter_leaf(pft)))			!Doug 01/13: Leaf or wood rate?

          fuel_1hr_leaf(pft,1)=fuel_1hr_leaf(pft,1)
     *      -fuel_decom_1hr_leaf(pft)
	 
          fuel_1hr_wood(pft,1)=fuel_1hr_wood(pft,1)
     *      -fuel_decom_1hr_wood(pft)

          
          fuel_decom_10hr(pft)=fuel_10hr(pft,1)*
     *      (1.0-exp(-k_litter_wood(pft)))
     
          fuel_10hr(pft,1)=fuel_10hr(pft,1)-fuel_decom_10hr(pft)
		  
          
            
c          fuel_decom_1hr(pft)=fuel_1hr(pft,1)*(1.0-exp(-k_litter))

C          fuel_1hr(pft,1)=fuel_1hr(pft,1)-fuel_decom_1hr(pft)

c          fuel_1hr(pft)=fuel_1hr(pft)+0.02*fuel_10hr(pft) !A
c          fuel_10hr(pft)=fuel_10hr(pft)-0.02*fuel_10hr(pft) !A

          fuel_decom_100hr(pft)=fuel_100hr(pft,1)*
     *      (1.0-exp(-k_litter_wood(pft)))
          fuel_100hr(pft,1)=fuel_100hr(pft,1)-fuel_decom_100hr(pft)

c          fuel_1hr(pft)=fuel_1hr(pft)+0.02*fuel_100hr(pft)  !A
c          fuel_100hr(pft)=fuel_100hr(pft)-0.02*fuel_100hr(pft) !A


c      ACHTUNG: Aren't boles decomposed mechanistically first, i.e. enter other
c               fuel classes rather than a decomposed to SOM directly?
          fuel_decom_1000hr(pft)=fuel_1000hr(pft,1)*
     *      (1.0-exp(-k_litter_wood(pft)))
          fuel_1000hr(pft,1)=fuel_1000hr(pft,1)-fuel_decom_1000hr(pft)

c          fuel_1hr(pft)=fuel_1hr(pft) +0.02*fuel_1000hr(pft) !A
c          fuel_1000hr(pft)=fuel_1000hr(pft)-0.02*fuel_1000hr(pft) !A

        enddo

        if (litter_decom(m,1).gt.0.0) then
           litter_decom(m,2)=litter_decom(m,2)/litter_decom(m,1)
           litter_decom(m,3)=litter_decom(m,3)/litter_decom(m,1)
        endif

		IF (1==0) THEN !Doug 11/12: lets turn agriculate off while we're not using it
c       Doug 05/09: agri: calculate the decomposition from the agricultural litter pools 
        agri_decom_ag=agri_litter_ag*(1.0-exp(-k_litter_bg)) 
        agri_decom_bg=agri_litter_bg*(1.0-exp(-k_litter_bg))  
        agri_decom(m)=agri_decom_ag+agri_decom_bg   

c       Doug 06/09 Update the litter pools   
        agri_litter_ag=agri_litter_ag-agri_decom_ag 
        agri_litter_bg=agri_litter_bg-agri_decom_bg 

c       Doug 06/09: pass forward either agricultural decomposition or natural
        if (agri) then   ! agricultural grid cel
           litter_decom(m,1)=agri_decom(m)
           litter_decom(m,2:nco2)=0

        else             ! Doug 06/09:natural grid cell (+remainder of agri litter)
           agri_decom(m)=0
           litter_decom(m,1)=nat_decom(m,1)+agri_decom(m)
           litter_decom(m,2:nco2)=nat_decom(m,2:nco2)

        endif
		
        END IF !(0==1) 

c       Calculate carbon flux to atmosphere and soil
        cflux_litter_atmos=atmfrac*litter_decom(m,1)
        cflux_litter_soil=soilfrac*litter_decom(m,1)

c       Further subdivide soil fraction between fast and slow soil pools

        temp=cpool_fast(1)
        cpool_fast(1)=cpool_fast(1)+fastfrac*cflux_litter_soil
        if (cpool_fast(1).gt.0.0) then
           do nc=2,nco2
              cpool_fast(nc)=(cpool_fast(nc)*temp+fastfrac
     *             *litter_decom(m,nc)*cflux_litter_soil)/cpool_fast(1)
           enddo
        endif
        
        temp=cpool_slow(1)
        cpool_slow(1)=cpool_slow(1)+(1.0-fastfrac)*cflux_litter_soil
        if (cpool_slow(1).gt.0.0) then
           do nc=2,nco2
              cpool_slow(nc) =(cpool_slow(nc)*temp+(1.0-fastfrac)
     *             *cflux_litter_soil*litter_decom(m,nc))/cpool_slow(1)
           enddo
        endif
        
c     Calculate monthly soil decomposition to the atmosphere

        cflux_fast_atmos=cpool_fast(1)*(1.0-exp(-k_fast(m))) !eqn 4
        cflux_slow_atmos=cpool_slow(1)*(1.0-exp(-k_slow(m))) !eqn 4
        
c       Update the soil pools

        cpool_fast(1)=cpool_fast(1)-cflux_fast_atmos
        cpool_slow(1)=cpool_slow(1)-cflux_slow_atmos   

c       Calculate monthly heterotrophic respiration

        mrh(m,1)=cflux_litter_atmos+cflux_fast_atmos+cflux_slow_atmos

        if (mrh(m,1).gt.0.0) then
           do nc=2,nco2
              mrh(m,nc)=(cflux_litter_atmos*litter_decom(m,nc)+
     *             cflux_fast_atmos*cpool_fast(nc)+cflux_slow_atmos*
     *             cpool_slow(nc))/mrh(m,1)
           enddo
        endif

        arh(1)=arh(1)+mrh(m,1)
        arh(2)=arh(2)+mrh(m,2)*mrh(m,1)
        arh(3)=arh(3)+mrh(m,3)*mrh(m,1)

c       Empty soil pools below a minimum threshold

        if (cpool_fast(1).lt.1.0e-5) then
          do nc=1,nco2
           cpool_fast(nc)=0.0
          enddo
        endif
        if (cpool_slow(1).lt.1.0e-5) then
          do nc=1,nco2
           cpool_slow(nc)=0.0
          enddo
        endif


        DO nc=1,nco2						!Doug 01/09: MFL
          DO pft=1,npft

            IF (fuel_1hr_leaf(pft,1)-temp_fuel(pft,1,1)>0) THEN
              fuel_1hr_leaf_inc_pos(pft,nc,m)=
     *          fuel_1hr_leaf_inc_pos(pft,nc,m)+
     *          fuel_1hr_leaf(pft,nc)-temp_fuel(pft,nc,1)
            ELSE
              fuel_1hr_leaf_inc_neg(pft,nc,m)=
     *          fuel_1hr_leaf_inc_neg(pft,nc,m)+
     *          fuel_1hr_leaf(pft,nc)-temp_fuel(pft,nc,1)
            END IF
			
            IF (fuel_1hr_wood(pft,1)-temp_fuel(pft,1,2)>0) THEN
              fuel_1hr_wood_inc_pos(pft,nc,m)=
     *          fuel_1hr_wood_inc_pos(pft,nc,m)+
     *          fuel_1hr_wood(pft,nc)-temp_fuel(pft,nc,2)
            ELSE
              fuel_1hr_wood_inc_neg(pft,nc,m)=
     *          fuel_1hr_wood_inc_neg(pft,nc,m)+
     *          fuel_1hr_wood(pft,nc)-temp_fuel(pft,nc,2)
            END IF
			
c            IF (fuel_1hr(pft,1)-temp_fuel(pft,1,1)>0) THEN
c              fuel_1hr_inc_pos(pft,nc,m)=fuel_1hr_inc_pos(pft,nc,m)+
c     *          fuel_1hr(pft,nc)-temp_fuel(pft,nc,1)
c            ELSE
c              fuel_1hr_inc_neg(pft,nc,m)=fuel_1hr_inc_neg(pft,nc,m)+
c     *          fuel_1hr(pft,nc)-temp_fuel(pft,nc,1)
c            END IF
			
          END DO
	      fuel_10hr_inc(:,nc,m)=fuel_10hr_inc(:,nc,m)+
     *      fuel_10hr(:,nc)-temp_fuel(:,nc,3)
          fuel_100hr_inc(:,nc,m)=fuel_100hr_inc(:,nc,m)+
     *      fuel_100hr(:,nc)-temp_fuel(:,nc,4)
          fuel_1000hr_inc(:,nc,m)=fuel_1000hr_inc(:,nc,m)+
     *      fuel_1000hr(:,nc)-temp_fuel(:,nc,5)
        END DO							!Doug 01/09: MFL

      enddo   !month

      if (arh(1).gt.0.0) then
         arh(2)=arh(2)/arh(1)
         arh(3)=arh(3)/arh(1)
      endif


c     SOIL DECOMPOSITION EQUILIBRIUM CALCULATION
C     Doug 12/09: This next bit has been removed cos its no longer needed
c       and didnt actaul work anyway (although I left in the description,
c       incase you interested in how the fix worked). Please see any veriosn
c       before r180 for the fix
c     Analytical solution of differential flux equations for fast and slow
c     soil carbon pools.  Implemented after (soil_equil_year) simulation
c     years, when annual litter inputs should be close to equilibrium.  Assumes
c     average climate (temperature and soil moisture) from all years up to
c     soil_equil_year.

 

c       Analytically calculate pool sizes this year only

c       Rate of change of soil pool size = litter input - decomposition
c         (5) dc/dt = litter_decom - kc
c       At equilibrium,
c         (6) dc/dt = 0
c       From (5) & (6),
c         (7) c = litter_decom / k

         
         
      return
      end



c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE KILL
c     Removal of PFTs with negative annual C increment
c     Note: Killing of PFTs newly beyond their bioclimatic limits is done in
c     subroutine establishment
c
c     Doug 01/09: addition of fuel_1hr_inc as I/O

      subroutine kill_pft(bm_inc,present,tree,lm_ind,rm_ind,hm_ind,
     *  sm_ind,nind,litter_ag_leaf,litter_ag_wood, !Doug 11/12
     *  litter_bg,year,fuel_1hr_leaf,fuel_1hr_wood,fuel_10hr,
     *  fuel_100hr,fuel_1000hr,
     *  fuel_1hr_leaf_inc_pos, fuel_1hr_leaf_inc_neg, 
     *  fuel_1hr_wood_inc_pos, fuel_1hr_wood_inc_neg, 
     *  fuel_10hr_inc, fuel_100hr_inc, fuel_1000hr_inc)

      implicit none

c     PARAMETERS
      integer npft
        parameter (npft=13)
      integer nco2
        parameter (nco2=3)           !number of C variables: C,13C,14C

c     ARGUMENTS
      real bm_inc(1:npft,1:nco2),lm_ind(1:npft,1:nco2)
      real rm_ind(1:npft,1:nco2),hm_ind(1:npft,1:nco2)
      real sm_ind(1:npft,1:nco2)
      real nind(1:npft)
      REAL litter_ag_leaf(1:npft,1:nco2) !Doug 11/12
      REAL litter_ag_wood(1:npft,1:nco2) !Doug 11/12
      REAL litter_bg(1:npft,1:nco2)
      REAL fuel_1hr_leaf(1:npft,1:nco2)
      REAL fuel_1hr_wood(1:npft,1:nco2)
      REAL fuel_10hr(1:npft,1:nco2)
      real fuel_100hr(1:npft,1:nco2)
      real fuel_1000hr(1:npft,1:nco2)
      REAL fuel_1hr_leaf_inc_pos(1:npft,1:nco2,1:12)
      REAL fuel_1hr_leaf_inc_neg(1:npft,1:nco2,1:12)
      REAL fuel_1hr_wood_inc_pos(1:npft,1:nco2,1:12)
      REAL fuel_1hr_wood_inc_neg(1:npft,1:nco2,1:12)	!Doug 01/09: records variations
      REAL fuel_10hr_inc(1:npft,1:nco2,1:12)	!in fuel for each month (for FIRE)
      REAL fuel_100hr_inc(1:npft,1:nco2,1:12)
      REAL fuel_1000hr_inc(1:npft,1:nco2,1:12)
      REAL temp_fuel(1:nco2,1:5)
      logical present(1:npft),tree(1:npft)
      integer year

c     LOCAL VARIABLES
      integer pft
      integer nc
      INTEGER m
c      real kill_frac
c      real nind_kill
      real temp, litter_inc

      pft=0
      nc=0
      temp=0.0
      litter_inc=0.0
      do pft=1,npft

        temp_fuel(:,1)=fuel_1hr_leaf(pft,:)		!Doug 01/09: MFL
        temp_fuel(:,2)=fuel_1hr_wood(pft,:)		!Doug 01/09: MFL
        temp_fuel(:,3)=fuel_10hr(pft,:)		!Doug 01/09: MFL
        temp_fuel(:,4)=fuel_100hr(pft,:)	!Doug 01/09: MFL
        temp_fuel(:,5)=fuel_1000hr(pft,:)	!Doug 01/09: MFL

        if (present(pft)) then
          if (bm_inc(pft,1).lt.0.0) then  !negative C increment this year

            present(pft)=.false.   !remove PFT

c           Transfer killed biomass to litter

            if (tree(pft)) then
               temp=litter_ag_leaf(pft,1) !Doug 11/12
			   temp=litter_ag_wood(pft,1) !Doug 11/12
			   
			   litter_ag_leaf(pft,1)=litter_ag_leaf(pft,1)+
     *           lm_ind(pft,1)*nind(pft)
	 
               IF (litter_ag_leaf(pft,1) .gt. 0.0) THEN
                 DO nc=2,nco2
                  litter_inc=lm_ind(pft,1)*lm_ind(pft,nc)
                  litter_ag_leaf(pft,nc)=(litter_ag_leaf(pft,nc)*temp
     *                 +litter_inc*nind(pft))/litter_ag_leaf(pft,1)
                 END DO
               END IF
	 
	           litter_ag_wood(pft,1)=litter_ag_wood(pft,1)+
     *           (sm_ind(pft,1)+ hm_ind(pft,1))*nind(pft)
	 
               IF (litter_ag_wood(pft,1) .gt. 0.0) THEN
                 DO nc=2,nco2
                  litter_inc=sm_ind(pft,1)
     *                 *sm_ind(pft,nc)+hm_ind(pft,1)*hm_ind(pft,nc)
                  litter_ag_wood(pft,nc)=(litter_ag_wood(pft,nc)*temp
     *                 +litter_inc*nind(pft))/litter_ag_wood(pft,1)
                 END DO
               END IF
			   
c               temp=litter_ag(pft,1)
c               litter_ag(pft,1)=litter_ag(pft,1)+(lm_ind(pft,1)+
c     *              sm_ind(pft,1)+ hm_ind(pft,1))*nind(pft)
C               if (litter_ag(pft,1) .gt. 0.0) then
C               do nc=2,nco2
C                  litter_inc=lm_ind(pft,1)*lm_ind(pft,nc)+sm_ind(pft,1)
C     *                 *sm_ind(pft,nc)+hm_ind(pft,1)*hm_ind(pft,nc)
C                  litter_ag(pft,nc)=(litter_ag(pft,nc)*temp
C     *                 +litter_inc*nind(pft))/litter_ag(pft,1)
C               enddo
C               endif

c            KIRSTEN: fuel classes

               temp=fuel_1hr_leaf(pft,1)
               fuel_1hr_leaf(pft,1)=fuel_1hr_leaf(pft,1)+lm_ind(pft,1)*
     *            nind(pft)
               IF (fuel_1hr_leaf(pft,1) .gt. 0.0) THEN
                 DO nc=2,nco2
                  litter_inc=lm_ind(pft,1)*lm_ind(pft,nc)
                  fuel_1hr_leaf(pft,nc)=(fuel_1hr_leaf(pft,nc)*temp
     *                 +litter_inc*nind(pft))/fuel_1hr_leaf(pft,1)
                 END DO
               END IF
			   
               temp=fuel_1hr_wood(pft,1)
               fuel_1hr_wood(pft,1)=fuel_1hr_wood(pft,1)+
     *              (0.045*sm_ind(pft,1)+ 0.045*hm_ind(pft,1))*nind(pft)
               IF (fuel_1hr_wood(pft,1) .gt. 0.0) THEN
                 DO nc=2,nco2
                  litter_inc=sm_ind(pft,1)*sm_ind(pft,nc)*0.045
     *                 +hm_ind(pft,1)*hm_ind(pft,nc)+0.045
                  fuel_1hr_wood(pft,nc)=(fuel_1hr_wood(pft,nc)*temp
     *                 +litter_inc*nind(pft))/fuel_1hr_wood(pft,1)
                 END DO
               END IF

c               temp=fuel_1hr(pft,1)
c               fuel_1hr(pft,1)=fuel_1hr(pft,1)+(lm_ind(pft,1)+
c     *              0.045*sm_ind(pft,1)+ 0.045*hm_ind(pft,1))*nind(pft)
c               if (fuel_1hr(pft,1) .gt. 0.0) then
c               do nc=2,nco2
c                  litter_inc=lm_ind(pft,1)*lm_ind(pft,nc)
c     *                 +sm_ind(pft,1)*sm_ind(pft,nc)*0.045
c     *                 +hm_ind(pft,1)*hm_ind(pft,nc)+0.045
c                  fuel_1hr(pft,nc)=(fuel_1hr(pft,nc)*temp
c     *                 +litter_inc*nind(pft))/fuel_1hr(pft,1)
c               enddo
c               endif

               temp=fuel_10hr(pft,1)
               fuel_10hr(pft,1)=fuel_10hr(pft,1)+((0.075*sm_ind(pft,1)
     *              +0.075*hm_ind(pft,1))*nind(pft))
               if (fuel_10hr(pft,1) .gt. 0.0) then
               do nc=2,nco2
                  litter_inc=sm_ind(pft,1)*sm_ind(pft,nc)*0.075
     *                 +hm_ind(pft,1)*hm_ind(pft,nc)+0.075
                  fuel_10hr(pft,nc)=(fuel_10hr(pft,nc)*temp
     *                 +litter_inc*nind(pft))/fuel_10hr(pft,1)
               enddo
               endif

               temp=fuel_100hr(pft,1)
               fuel_100hr(pft,1)=fuel_100hr(pft,1)+((0.21*sm_ind(pft,1)
     *              +0.21*hm_ind(pft,1))*nind(pft))
               if (fuel_100hr(pft,1) .gt. 0.0) then
               do nc=2,nco2
                  litter_inc=sm_ind(pft,1)*sm_ind(pft,nc)*0.21
     *                 +hm_ind(pft,1)*hm_ind(pft,nc)+0.21
                  fuel_100hr(pft,nc)=(fuel_100hr(pft,nc)*temp
     *                 +litter_inc*nind(pft))/fuel_100hr(pft,1)
               enddo
               endif

               temp=fuel_1000hr(pft,1)
               fuel_1000hr(pft,1)=fuel_1000hr(pft,1)+((0.67*
     *              *sm_ind(pft,1) 
     *              +0.67*hm_ind(pft,1))*nind(pft))
               if (fuel_1000hr(pft,1) .gt. 0.0) then
               do nc=2,nco2
                  litter_inc=sm_ind(pft,1)*sm_ind(pft,nc)*0.67
     *                 +hm_ind(pft,1)*hm_ind(pft,nc)+0.67
                  fuel_1000hr(pft,nc)=(fuel_1000hr(pft,nc)*temp
     *                 +litter_inc*nind(pft))/fuel_1000hr(pft,1)
               enddo
               endif

            else  !grasses
               temp=litter_ag_leaf(pft,1) !Doug 11/12
               litter_ag_leaf(pft,1)=litter_ag_leaf(pft,1)+lm_ind(pft,1)
     *              *nind(pft)
               if (litter_ag_leaf(pft,1) .gt. 0.0) then
               do nc=2,nco2
                  litter_inc=lm_ind(pft,1)*lm_ind(pft,nc)
                  litter_ag_leaf(pft,nc)=(litter_ag_leaf(pft,nc)*temp
     *                 +litter_inc*nind(pft))/litter_ag_leaf(pft,1)
               enddo
               endif

c            KIRSTEN: fuel classes
               temp=fuel_1hr_leaf(pft,1)
               fuel_1hr_leaf(pft,1)=fuel_1hr_leaf(pft,1)+lm_ind(pft,1)
     *              *nind(pft)
               if (fuel_1hr_leaf(pft,1) .gt. 0.0) then
               do nc=2,nco2
                  litter_inc=lm_ind(pft,1)*lm_ind(pft,nc)
                  fuel_1hr_leaf(pft,nc)=(fuel_1hr_leaf(pft,nc)*temp
     *                 +litter_inc*nind(pft))/fuel_1hr_leaf(pft,1)
               enddo
               endif

            endif

            temp=litter_bg(pft,1)
            litter_bg(pft,1)=litter_bg(pft,1)+rm_ind(pft,1)*nind(pft)
            if (litter_bg(pft,1) .gt. 0.0) then
            do nc=2,nco2
               litter_bg(pft,nc)=(litter_bg(pft,nc)*temp+rm_ind(pft,nc)
     *              *rm_ind(pft,1)*nind(pft))/litter_bg(pft,1)
            enddo
            endif

          endif

        endif
		
        
        DO nc=1,nco2
          DO m=1,12						!Doug 01/09: MFL !StopC
            IF (fuel_1hr_leaf(pft,1)-temp_fuel(1,1)>0) THEN
              fuel_1hr_leaf_inc_pos(pft,nc,m)=
     *          fuel_1hr_leaf_inc_pos(pft,nc,m)+
     *          (fuel_1hr_leaf(pft,nc)-temp_fuel(nc,1))/12
            ELSE
              fuel_1hr_leaf_inc_neg(pft,nc,m)=
     *          fuel_1hr_leaf_inc_neg(pft,nc,m)+
     *          (fuel_1hr_leaf(pft,nc)-temp_fuel(nc,1))/12
            END IF
		  
            IF (fuel_1hr_wood(pft,1)-temp_fuel(1,2)>0) THEN
              fuel_1hr_wood_inc_pos(pft,nc,m)=
     *          fuel_1hr_wood_inc_pos(pft,nc,m)+
     *          (fuel_1hr_wood(pft,nc)-temp_fuel(nc,2))/12
            ELSE
              fuel_1hr_wood_inc_neg(pft,nc,m)=
     *          fuel_1hr_wood_inc_neg(pft,nc,m)+
     *          (fuel_1hr_wood(pft,nc)-temp_fuel(nc,2))/12
            END IF

            fuel_10hr_inc(pft,nc,:)=fuel_10hr_inc(pft,nc,m)+
     *        (fuel_10hr(pft,nc)-temp_fuel(nc,3))/12
            fuel_100hr_inc(pft,nc,:)=fuel_100hr_inc(pft,nc,m)+
     *        (fuel_100hr(pft,nc)-temp_fuel(nc,4))/12
            fuel_1000hr_inc(pft,nc,:)=fuel_1000hr_inc(pft,nc,m)+
     *        (fuel_1000hr(pft,nc)-temp_fuel(nc,5))/12
          END DO							!Doug 01/09: MFL
        END DO

      enddo	!pft

      return
      end

c      kill_frac=0.1

c      do pft=1,npft

c        if (present(pft)) then

c          if (bm_inc(pft).lt.0.0) then  !negative C increment this year

c            present(pft)=.false.   !remove PFT
c             bm_inc(pft)=0.0

c           Transfer killed biomass to litter


c            if (tree(pft)) then
c              nind_kill=nind(pft)*kill_frac
c              nind(pft)=nind(pft)-nind_kill
c              litter_ag(pft)=litter_ag(pft)+(lm_ind(pft)+sm_ind(pft)+
c     *          hm_ind(pft))*nind_kill
c              litter_bg(pft)=litter_bg(pft)+rm_ind(pft)*nind_kill
c            else  !grasses
c              litter_ag(pft)=litter_ag(pft)+lm_ind(pft)*kill_frac
c              litter_bg(pft)=litter_bg(pft)+rm_ind(pft)*kill_frac
c              lm_ind(pft)=lm_ind(pft)-lm_ind(pft)*kill_frac
c              rm_ind(pft)=rm_ind(pft)-rm_ind(pft)*kill_frac
c            endif


c            if (tree(pft)) then
c              litter_ag(pft)=litter_ag(pft)+(lm_ind(pft)+sm_ind(pft)+
c     *          hm_ind(pft))*nind(pft)
c            else  !grasses
c              litter_ag(pft)=litter_ag(pft)+lm_ind(pft)*nind(pft)
c            endif
c            litter_bg(pft)=litter_bg(pft)+rm_ind(pft)*nind(pft)

c          endif

c        endif

c      enddo

c      return
c      end



c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE ALLOCATION
c     Allocation of annual C increment to leaf, stem and fine root
c     compartments, update of individual structure
c
c     Doug 01/09: addition of fuel_1hr_inc as I/O

      subroutine allocation(pftpar,allom1,allom2,allom3,latosa,wooddens,
     * reinickerp,tree,sla,wscal,nind,
     * bm_inc,lm_ind,lm_inc,sm_ind,hm_ind,rm_ind,
     * crownarea,fpc_grid,lai_ind,height,height_class,dbh,dbh_class, !Doug 02/13: dbh_class added
     * tau_c,cl_t,
     * litter_ag_leaf,litter_ag_wood, ! Doug 11/12
     * litter_bg,fuel_1hr_leaf,fuel_1hr_wood,
     * fuel_10hr,fuel_100hr,fuel_1000hr,
     * fuel_1hr_leaf_inc_pos, fuel_1hr_leaf_inc_neg,
     * fuel_1hr_wood_inc_pos, fuel_1hr_wood_inc_neg,
     * fuel_10hr_inc,
     * fuel_100hr_inc, fuel_1000hr_inc, fpc_inc, present, year,
     * evergreen)

      implicit none

c     PARAMETERS
      integer npft,npftpar
      parameter (npft=13,npftpar=59)
      integer nco2
        parameter (nco2=3)           !number of C variables: C,13C,14C
      real pi,xacc,yacc
      integer nseg
        parameter (pi=3.14159265)
        parameter (xacc=0.1)     !threshold x-axis precision of allocation soln
        parameter (yacc=1.0e-10) !threshold y-axis precision of allocation soln
        parameter (nseg=20)
      integer jmax
        parameter (jmax=40) !max number of iterations in search for alloc soln

c     ARGUMENTS
      real pftpar(1:npft,1:npftpar)
      real allom1,allom2,allom3
      real latosa,wooddens,reinickerp
      logical tree(1:npft)
      real sla(1:npft)
      real wscal(1:npft)
      real nind(1:npft)
      real bm_inc(1:npft,1:nco2)
      REAL lm_ind(1:npft,1:nco2)
      REAL lm_inc(1:npft)   ! Doug 06/13: total fine mass incrememet for outputs
      REAL rm_ind(1:npft,1:nco2)
      real sm_ind(1:npft,1:nco2),hm_ind(1:npft,1:nco2)
      real crownarea(1:npft)
      real lai_ind(1:npft),fpc_grid(1:npft)
      real height(1:npft)
      real height_class(0:4,1:npft)      
      REAL dbh(1:npft), dbh_class(0:4,1:npft) !Doug 02/13:dbh_class added
      REAL tau_c(0:4,1:npft),cl_t(0:4,1:npft)
      REAL litter_ag_leaf(1:npft,1:nco2) ! Doug 11/12
      REAL litter_ag_wood(1:npft,1:nco2) ! Doug 11/12
      REAL litter_bg(1:npft,1:nco2)
      REAL fuel_1hr_leaf(1:npft,1:nco2)
      REAL fuel_1hr_wood(1:npft,1:nco2)
      REAL fuel_10hr(1:npft,1:nco2)
      real fuel_100hr(1:npft,1:nco2)
      real fuel_1000hr(1:npft,1:nco2)
      REAL fuel_1hr_leaf_inc_pos(1:npft,1:nco2,1:12)
      REAL fuel_1hr_leaf_inc_neg(1:npft,1:nco2,1:12)	!Doug 01/09: records variations
      REAL fuel_1hr_wood_inc_pos(1:npft,1:nco2,1:12)
      REAL fuel_1hr_wood_inc_neg(1:npft,1:nco2,1:12)	!Doug 01/09: records variations
      REAL fuel_10hr_inc(1:npft,1:nco2,1:12)	!in fuel for each month (for FIRE)
      REAL fuel_100hr_inc(1:npft,1:nco2,1:12)
      REAL fuel_1000hr_inc(1:npft,1:nco2,1:12)
      REAL temp_fuel(1:nco2,5)
      real fpc_inc(1:npft)
      logical present(1:npft)
      integer year
      logical evergreen(1:npft)


c     LOCAL VARIABLES
      real param1(1:npft),param2(1:npft)
      integer pft,ppft
      integer nc
      real lminc_ind ! individual leafmass increment this year
      real sminc_ind ! individual sapmass increment this year
      real rminc_ind ! individual fineroot mass increment this year
      real bm_inc_ind ! individual total biomass increment this year
      real sap_xsa ! cross sectional area of sapwood
      real stemdiam ! stem diameter
      real x1,x2,rtbis,dx,xmid,sign !working vars in bisection
      double precision fx1,fmid
      real lminc_ind_min  !min leafmass increment to maintain current sapwood
      real rminc_ind_min  !min rootmass increment to support new leafmass
      integer j  !counter in bisection loop
      real lmtorm ! ratio of leafmass to fine rootmass
      real fpc_grid_old  !previous year's FPC
      real fpc_ind  !individual FPC
      real crownarea_max  !maximum crown area (m2)
c  Kirsten
      integer class
      real bt(0:4,1:npft),crown(1:npft) !crown length,bark thickness (cm?²)  
      real temp,lm_temp,rm_temp,sm_temp
      logical normal
      REAL lm_tot(1:npft)

c DM  One very nasty bug was really hard to find because rm_temp was not initialized.
c DM  So, local variables should be initialized everywhere.
c DM  If the initialization value is not correct, you'll have to change it.
c DM  If you are sure the initialization is not needed, comment out
      param1(:)=0.0
      param2(:)=0.0
      pft=0
      nc=0
      lminc_ind=0.0
      sminc_ind=0.0
      rminc_ind=0.0
      bm_inc_ind=0.0
      sap_xsa=0.0
      stemdiam=0.0
      x1=0.0
      x2=0.0
      rtbis=0.0
      dx=0.0
      xmid=0.0
      sign=0.0
      fx1=0.0
      fmid=0.0
      lminc_ind_min=0.0
      rminc_ind_min=0.0
      j=0
      lmtorm=0.0
      fpc_grid_old=0.0
      fpc_ind=0.0
      crownarea_max=0.0
      class=0
      bt(:,:)=0.0
      crown(:)=0.0
      dbh_class(:,:)=0.0 
      lm_inc(:)=0.0      
c#ifdef CO2_THREE_VALUES
c DM   Worse: given how rm_temp is used, it should be initialised
c DM   at the beginning of the main loop (why should the value of
c DM   rm_ind(pft1, nc1) depend on the one of rm_ind(pft2, nc2)?)
c DM   We got infinities for some variables otherwise
c DM   We do the same for the others
c DM   TODO : what other variables should be given the same treatment ?
c      temp=0.0
c      lm_temp=0.0
c      rm_temp=0.0
c      sm_temp=0.0
c#endif

      do pft=1,npft

        temp_fuel(:,1)=fuel_1hr_leaf(pft,:)		!Doug 01/09: MFL
        temp_fuel(:,2)=fuel_1hr_wood(pft,:)		!Doug 01/09: MFL
        temp_fuel(:,3)=fuel_10hr(pft,:)			!Doug 01/09: MFL
        temp_fuel(:,4)=fuel_100hr(pft,:)		!Doug 01/09: MFL
        temp_fuel(:,5)=fuel_1000hr(pft,:)		!Doug 01/09: MFL

        temp=0.0
        lm_temp=0.0
        rm_temp=0.0
        sm_temp=0.0

        if (present(pft)) then

          bm_inc_ind=bm_inc(pft,1)/nind(pft)

c         calculate this year's leaf to fine root mass ratio from mean annual
c         water scalar and pft specific parameter

          lmtorm=pftpar(pft,18)*wscal(pft)

          if (tree(pft)) then

            normal=.false.
            crownarea_max=pftpar(pft,20)

c           PFT parameter crown length & bark thickness, Kirsten            
            crown(pft)=pftpar(pft,44)
            param1(pft)=pftpar(pft,47)
            param2(pft)=pftpar(pft,48)

c           TREE ALLOCATION

c           Allocation of this year's biomass increment (bm_inc_ind) to the
c           three living carbon pools, such that the basic allometric
c           relationships (A-C below) are always satisfied.

c           (A) (leaf area) = latosa * (sapwood xs area)
c                 (Pipe Model, Shinozaki et al. 1964a,b; Waring et al 1982)
c           (B) (leaf mass) = lmtorm * (root mass)
c           (C) height = allom2 * (stem diameter)**allom3
c                 (source?)
c           (D) (crown area) = min (allom1 * (stem diameter)**reinickerp,
c                                   crownarea_max)

c           Mathematical derivation:

c             (1) bm_inc_ind = lminc_ind + sminc_ind + rminc_ind
c             (2) leaf_area_new = latosa * sap_xsa_new   [from (A)]
c             (3) leaf_area_new = (lm_ind + lminc_ind) * sla
c           from (2) & (3),
c             (4) (lm_ind + lminc_ind) * sla = latosa * sap_xsa_new
c           from (4),
c             (5) sap_xsa_new = (lm_ind + lminc_ind) * sla / latosa
c             (6) (lm_ind + lminc_ind) = lmtorm * (rm_ind + rminc_ind)
c                   [from (B)]
c             (7) height_new = allom2 * stemdiam_new**allom3  [from (C)]
c           from (1),
c             (8) sminc_ind = bm_inc_ind - lminc_ind - rminc_ind
c           from (6),
c             (9) rminc_ind=((lm_ind + lminc_ind) / lmtorm) - rm_ind
c           from (8) & (9),
c            (10) sminc_ind = bm_inc_ind - lminc_ind
c                   - ((lm_ind + lminc_ind)  / lmtorm) + rm_ind
c            (11) wooddens = (sm_ind + sminc_ind + hm_ind) / stemvolume_new
c            (12) stemvolume_new = height_new * pi * stemdiam_new**2 / 4
c           from (10), (11) & (12)
c            (13) stemdiam_new = [ ((sm_ind + bm_inc_ind - lminc_ind
c                   - ((lm_ind + lminc_ind) / lmtorm) + rm_ind + hm_ind)
c                   / wooddens) / (height_new * pi / 4) ]**(1/2)
c           combining (7) and (13),
c            (14) height_new = allom2 * [ ((sm_ind + bm_inc_ind - lminc_ind
c                   - ((lm_ind + lminc_ind) / lmtorm) + rm_ind + hm_ind)
c                   / wooddens) / (height_new * pi / 4) ]**(1/2 * allom3)
c           from (14),
c            (15) height_new**(1 + 2 / allom3) = allom2**(2 / allom3)
c                   * ((sm_ind + bm_inc_ind - lminc_ind - ((lm_ind + lminc_ind)
c                   / lmtorm) + rm_ind + hm_ind) / wooddens) / (pi / 4)
c            (16) wooddens = (sm_ind + sminc_ind) / sapvolume_new
c           from (10) and (16),
c            (17) wooddens = (sm_ind + bm_inc_ind - lminc_ind
c                   - ((lm_ind + lminc_ind) / lmtorm) + rm_ind) / sapvolume_new
c            (18) sapvolume_new = height_new * sap_xsa_new
c           from (17) and (18),
c            (19) sap_xsa_new = (sm_ind + bm_inc_ind - lminc_ind
c                   - ((lm_ind + lminc_ind) / lmtorm) + rm_ind)
c                   / (height_new * wooddens)
c           from (19),
c            (20) height_new = (sm_ind + bm_inc_ind - lminc_ind
c                   - ((lm_ind + lminc_ind) / lmtorm) + rm_ind )
c                   / (sap_xsa_new * wooddens)
c           from (5) and (20),
c            (21) height_new**(1 + 2 / allom3) = [ (sm_ind + bm_inc_ind
c                   - lminc_ind - ((lm_ind + lminc_ind) / lmtorm) + rm_ind )
c                   / ((lm_ind + lminc_ind) * sla * wooddens / latosa) ]
c                   **(1 + 2 / allom3)
c           -------------------------------------------------------------------
c            (15) and (21) are two alternative expressions for
c                 height_new**(1 + 2 / allom3). Combining these,

c            (22) allom2**(2 / allom3) * ((sm_ind + bm_inc_ind - lminc_ind
c                   - ((lm_ind + lminc_ind) / lmtorm) + rm_ind + hm_ind)
c                   / wooddens) / (pi / 4) - [ (sm_ind + bm_inc_ind - lminc_ind
c                   - ((lm_ind + lminc_ind) / lmtorm) + rm_ind )
c                   / ((lm_ind + lminc_ind) * sla * wooddens / latosa) ]
c                   **(1 + 2 / allom3)
c                   = 0

c           Equation (22) can be expressed in the form f(lminc_ind)=0.

c           Numerical methods are used to solve the equation for the
c           unknown lminc_ind.
c           -------------------------------------------------------------------

c           Work out minimum leaf production to maintain current sapmass

c            (23) sap_xsa = sm_ind / wooddens / height
c           from (A) and (23),
c            (24) leaf_mass * sla = latosa * sap_mass / wooddens / height
c           from (24),
c            (25) leaf_mass = latosa * sap_mass / (wooddens * height * sla)
c           from (25), assuming sminc_ind=0,
c            (26) lm_ind + lminc_ind_min = latosa * sm_ind
c                   / (wooddens * height * sla)
c           from (26),
c            (27) lminc_ind_min = latosa * sm_ind / (wooddens * height * sla)
c                   - lm_ind

            lminc_ind_min=latosa*sm_ind(pft,1)/
     *        (wooddens*height(pft)*sla(pft))-lm_ind(pft,1)  !eqn (27)

c           Work out minimum root production to support this leaf mass
c           (i.e. lm_ind + lminc_ind_min)
c           May be negative following a reduction in soil water limitation
c           (increase in lmtorm) relative to last year.

c           from (B) and (25),
c            (28) root_mass = latosa * sap_mass / (wooddens * height * sla)
c                   / lmtorm
c           from (28), assuming sminc_ind=0,
c            (29) rm_ind + rminc_ind_min = latosa * sm_ind
c                   / (wooddens * height * sla * lmtorm)
c           from (29),
c            (30) rminc_ind_min = latosa * sm_ind
c                   / (wooddens * height * sla * lmtorm) - rm_ind


            rminc_ind_min=latosa*sm_ind(pft,1)/(wooddens*height(pft)*
     *        sla(pft)*lmtorm)-rm_ind(pft,1)      !eqn (30)

            if (rminc_ind_min.gt.0.0.and.lminc_ind_min.gt.0.0.and.
     *        rminc_ind_min+lminc_ind_min.le.bm_inc_ind.or.
     *        bm_inc_ind.le.0.0) then

               normal=.true.

c             Normal allocation (positive increment to all living C
c             compartments)

c             Calculation of leaf mass increment (lminc_ind) that satisfies
c             Eqn (22) using Bisection Method (Press et al 1986, p 346)

c             Seeking a root for non-negative lminc_ind, rminc_ind and
c             sminc_ind.  There should be exactly one (no proof presented, but
c             Steve has managed one) and it should lie between x1=0 and
c             x2=(bm_inc_ind-(lm_ind/lmtorm-rm_ind))/(1+1/lmtorm).

              x1=0.0
              x2=(bm_inc_ind-(lm_ind(pft,1)/lmtorm-rm_ind(pft,1)))/
     *          (1.0+1.0/lmtorm)
              dx=(x2-x1)/real(nseg)

              if (lm_ind(pft,1).eq.0.0) x1=x1+dx  !to avoid division by zero

c             evaluate f(x1)=LHS of eqn (22) at x1

              fx1=allom2**(2.0/allom3)*((sm_ind(pft,1)+bm_inc_ind-x1
     *         -((lm_ind(pft,1)+x1)/lmtorm)+rm_ind(pft,1)+hm_ind(pft,1))
     *         /wooddens)/(pi/4.0)-((sm_ind(pft,1)+bm_inc_ind-x1
     *         -((lm_ind(pft,1)+x1)/lmtorm)+rm_ind(pft,1))
     *         /((lm_ind(pft,1)+x1)*sla(pft)*wooddens/latosa))
     *         **(1.0+2.0/allom3)

c             Find approximate location of leftmost root on the interval
c             (x1,x2).  Subdivide (x1,x2) into nseg equal segments seeking
c             change in sign of f(xmid) relative to f(x1).

              fmid=fx1
              xmid=x1

              do while (fmid*fx1.gt.0.0.and.xmid.lt.x2)

                xmid=xmid+dx
                fmid=allom2**(2.0/allom3)*((sm_ind(pft,1)+
     *            bm_inc_ind-xmid-((lm_ind(pft,1)+xmid)/lmtorm)
     *            +rm_ind(pft,1)+hm_ind(pft,1))
     *            /wooddens)/(pi/4.0)-((sm_ind(pft,1)+bm_inc_ind-xmid
     *            -((lm_ind(pft,1)+xmid)/lmtorm)+rm_ind(pft,1))
     *            /((lm_ind(pft,1)+xmid)*sla(pft)*wooddens/latosa))
     *            **(1.0+2.0/allom3)

              enddo


              x1=xmid-dx
              x2=xmid

c             Apply bisection method to find root on the new interval (x1,x2)

              fx1=allom2**(2.0/allom3)*((sm_ind(pft,1)+bm_inc_ind-x1
     *          -((lm_ind(pft,1)+x1)/lmtorm)+rm_ind(pft,1)+
     *          hm_ind(pft,1))
     *          /wooddens)/(pi/4.0)-((sm_ind(pft,1)+bm_inc_ind-x1
     *          -((lm_ind(pft,1)+x1)/lmtorm)+rm_ind(pft,1))
     *          /((lm_ind(pft,1)+x1)*sla(pft)*wooddens/latosa))
     *          **(1.0+2.0/allom3)

              if (fx1.ge.0.0) then
                sign=-1.0
              else
                sign=1.0
              endif

              rtbis=x1
              dx=x2-x1

c             Bisection loop
c             Search iterates on value of xmid until xmid lies within
c             xacc of the root, i.e. until |xmid-x|<xacc where f(x)=0

              fmid=1.0  !dummy value to guarantee entry to loop
              j=0.0     !number of iterations so far (maximum tries=jmax)

              do while (dx.ge.xacc.and.abs(fmid).gt.yacc)

                dx=dx*0.5
                xmid=rtbis+dx

c               calculate fmid=f(xmid) [eqn (22)]

                fmid=allom2**(2.0/allom3)*((sm_ind(pft,1)+bm_inc_ind
     *           -xmid-((lm_ind(pft,1)+xmid)/lmtorm)+rm_ind(pft,1)
     *            +hm_ind(pft,1))
     *            /wooddens)/(pi/4.0)-((sm_ind(pft,1)+bm_inc_ind-xmid
     *            -((lm_ind(pft,1)+xmid)/lmtorm)+rm_ind(pft,1))
     *            /((lm_ind(pft,1)+xmid)*sla(pft)*wooddens/latosa))
     *            **(1.0+2.0/allom3)

                if (fmid*sign.le.0.0) rtbis=xmid
                j=j+1

              enddo

c             Now rtbis contains numerical solution for lminc_ind given
c             eqn (22)

              lminc_ind=rtbis

c             Calculate increments in other compartments using allometry
c             relationships

              rminc_ind=(lm_ind(pft,1)+lminc_ind)/lmtorm-rm_ind(pft,1) !eqn (9)
              sminc_ind=bm_inc_ind-lminc_ind-rminc_ind    !eqn (1)

            else

c             Abnormal allocation: reduction in some C compartment(s)
c             to satisfy allometry

c             Attempt to distribute this year's production among leaves and
c             roots only
c              (31) bm_inc_ind = lminc_ind + rminc_ind
c             from (31) and (9),
c              (32) bm_inc_ind = lminc_ind + ((lm_ind + lminc_ind) / lmtorm)
c                     - rm_ind
c             from (32)
c              (33) lminc_ind = (bm_inc_ind - lmind / lmtorm + rm_ind) /
c                     (1 + 1 / lmtorm)

              lminc_ind=(bm_inc_ind-lm_ind(pft,1)/lmtorm+rm_ind(pft,1))
     *          /(1.0+1.0/lmtorm)  !eqn (33)

              if (lminc_ind.ge.0.0) then

c               Positive allocation to leafmass

                rminc_ind=bm_inc_ind-lminc_ind  !eqn (31)

c               Add killed roots (if any) to below-ground litter

                if (rminc_ind.lt.0.0) then
                  lminc_ind=bm_inc_ind
                  rminc_ind=((lm_ind(pft,1) + lminc_ind) / lmtorm)
     *                      - rm_ind(pft,1)
                  if (rminc_ind.lt.0.) then
                     temp=litter_bg(pft,1)
                     litter_bg(pft,1)=litter_bg(pft,1)+(-rminc_ind)
     *                    *nind(pft)                  
                     do nc=2,nco2
                        litter_bg(pft,nc)=(litter_bg(pft,nc)*temp+
     *                       (-rminc_ind)*rm_ind(pft,nc)*nind(pft))/
     *                       litter_bg(pft,1)
                     enddo
                  elseif (rminc_ind.ge.0.) then
                     temp=litter_bg(pft,1)
                     litter_bg(pft,1)=litter_bg(pft,1)+(rminc_ind)
     *                    *nind(pft)                  
                     do nc=2,nco2
                        litter_bg(pft,nc)=(litter_bg(pft,nc)*temp+
     *                       (rminc_ind)*rm_ind(pft,nc)*nind(pft))/
     *                       litter_bg(pft,1)
                     enddo
                  endif        
                endif
              else

c               Negative allocation to leaf mass

                rminc_ind=bm_inc_ind
          lminc_ind=(rm_ind(pft,1)+rminc_ind)*lmtorm-lm_ind(pft,1)
                            !from eqn (9)

c               Add killed leaves to litter

                temp=litter_ag_leaf(pft,1)   !Doug 11/12
                litter_ag_leaf(pft,1)=litter_ag_leaf(pft,1)
     *               +(-lminc_ind)*nind(pft)
                if (litter_ag_leaf(pft,1) .gt. 0.0) then
                do nc=2,nco2
                   litter_ag_leaf(pft,nc)=(litter_ag_leaf(pft,nc)*temp
     *                  +(-lminc_ind)*lm_ind(pft,nc)*nind(pft))
     *                  /litter_ag_leaf(pft,1) 
                enddo
                endif
c            KIRSTEN: fuel classes
                temp=fuel_1hr_leaf(pft,1)
                fuel_1hr_leaf(pft,1)=fuel_1hr_leaf(pft,1)
     *               +(-lminc_ind)*nind(pft)
                if (fuel_1hr_leaf(pft,1) .gt.0.0) then
                do nc=2,nco2
                   fuel_1hr_leaf(pft,nc)=(fuel_1hr_leaf(pft,nc)*temp
     *                  +(-lminc_ind)*lm_ind(pft,nc)*nind(pft))
     *                  /fuel_1hr_leaf(pft,1) 
                enddo
                endif
              endif

c             Calculate sminc_ind (must be negative)

c             from (25),
c              (34) lm_ind + lminc_ind = latosa * (sm_ind + sminc_ind)
c                     / (wooddens * height * sla)
c             from (34),
c              (35) sminc_ind = (lm_ind + lminc_ind) * wooddens * height * sla
c                     / latosa - sm_ind

              sminc_ind=(lm_ind(pft,1)+lminc_ind)*wooddens*height(pft)*
     *          sla(pft)/latosa-sm_ind(pft,1)                         !eqn (35)

c             Convert killed sapwood to heartwood

              temp=hm_ind(pft,1)
              hm_ind(pft,1)=hm_ind(pft,1)+(-sminc_ind)
              if (hm_ind(pft,1) .gt. 0.0) then
              do nc=2,nco2
                 hm_ind(pft,nc)=(hm_ind(pft,nc)*temp
     *                +(-sminc_ind)*sm_ind(pft,nc))/hm_ind(pft,1)
              enddo
              endif
            endif

c           Increment C compartments
            
            lm_temp=lm_ind(pft,1)
            rm_temp=rm_ind(pft,1)
            sm_temp=sm_ind(pft,1)
            
            lm_inc(pft)  =lminc_ind*nind(pft)
            lm_ind(pft,1)=lm_ind(pft,1)+lminc_ind
            rm_ind(pft,1)=rm_ind(pft,1)+rminc_ind
            sm_ind(pft,1)=sm_ind(pft,1)+sminc_ind

            if (normal) then            
               
               do nc=2,nco2
                  lm_ind(pft,nc)=(lm_ind(pft,nc)*lm_temp+lminc_ind
     *                 *bm_inc(pft,nc))/lm_ind(pft,1)
                  rm_ind(pft,nc)=(rm_ind(pft,nc)*rm_temp+rminc_ind
     *                 *bm_inc(pft,nc))/rm_ind(pft,1)
                  sm_ind(pft,nc)=(sm_ind(pft,nc)*sm_temp+sminc_ind
     *                 *bm_inc(pft,nc))/sm_ind(pft,1)
               enddo

            else

              if (lminc_ind.gt.0.0) then
                 do nc=2,nco2
                    lm_ind(pft,nc)=(lm_ind(pft,nc)*lm_temp
     *                   +lminc_ind*bm_inc(pft,nc))/lm_ind(pft,1)
                 enddo
              endif
              
              if (rminc_ind.gt.0.0) then
                 do nc=2,nco2
                    rm_ind(pft,nc)=(rm_ind(pft,nc)*rm_temp
     *                   +rminc_ind*bm_inc(pft,nc))/rm_ind(pft,1)
                 enddo
              endif

            endif
            !PRINT*, "wow"
			!		 print*, lm_ind(pft,1)
			!		 PRINT*, rm_ind(pft,1)					 
			!		 PRINT*, sm_ind(pft,1)
            if (lm_ind(pft,1).lt.0.0.or.rm_ind(pft,1).lt.0.0
     *         .or.sm_ind(pft,1).lt.0.0) then

               temp=litter_ag_leaf(pft,1) !Doug 11/12
               litter_ag_leaf(pft,1)=litter_ag_leaf(pft,1)+lm_ind(pft,1)
     *            *nind(pft)
               if (lm_ind(pft,1).gt.0) then
                  do nc=2,nco2
                     litter_ag_leaf(pft,nc)=(litter_ag_leaf(pft,nc)*temp
     *                    +lm_ind(pft,nc)*lm_ind(pft,1)*nind(pft))
     *                    /litter_ag_leaf(pft,1)
                  enddo
               endif
               lm_ind(pft,1)=0.0
               lm_ind(pft,2)=0.0
               lm_ind(pft,3)=0.0

               temp=litter_ag_wood(pft,1) !Doug 11/12
               litter_ag_wood(pft,1)=litter_ag_wood(pft,1)+sm_ind(pft,1)
     *            *nind(pft)
               if (sm_ind(pft,1).gt.0) then
                  do nc=2,nco2
                     litter_ag_wood(pft,nc)=(litter_ag_wood(pft,nc)*temp
     *                    +sm_ind(pft,nc)*sm_ind(pft,1)*nind(pft))
     *                    /litter_ag_wood(pft,1)
                  enddo
               endif

               temp=litter_ag_wood(pft,1)
               litter_ag_wood(pft,1)=litter_ag_wood(pft,1)+ !Doug 11/12
     *           ((hm_ind(pft,1)+bm_inc_ind)*nind(pft))
               if (litter_ag_wood(pft,1) .gt. 0.0) then
               do nc=2,nco2
                  litter_ag_wood(pft,nc)=(litter_ag_wood(pft,nc)*temp
     *                 +(hm_ind(pft,nc)*hm_ind(pft,1)+bm_inc_ind*
     *                 bm_inc(pft,nc))*nind(pft))/litter_ag_wood(pft,1)
               enddo
               endif

               temp=litter_bg(pft,1)
               litter_bg(pft,1)=litter_bg(pft,1)+(rm_ind(pft,1))
     *              *nind(pft)
               if (rm_ind(pft,1).gt.0) then
                  do nc=2,nco2
                     litter_bg(pft,nc)=(litter_bg(pft,nc)*temp
     *                    +rm_ind(pft,nc)*rm_ind(pft,1)*nind(pft))
     *                    /litter_bg(pft,1)
                  enddo
               endif
               rm_ind(pft,1)=0.0
               rm_ind(pft,2)=0.0
               rm_ind(pft,3)=0.0
               present(pft)=.false. !remove PFT
            endif 

c           Calculate new height, diameter and crown area

            sap_xsa=lm_ind(pft,1)*sla(pft)/latosa  !eqn (5)
            height(pft)=sm_ind(pft,1)/sap_xsa/wooddens
            stemdiam=(height(pft)/allom2)**(1.0/allom3)  !eqn (C)
			                     
            crownarea(pft)=min(allom1*stemdiam**reinickerp,
     *        crownarea_max)  !eqn (D)

c           Kirsten: stemdiam per pft and bark thickness, critical time to cambial kill
            dbh(pft)=stemdiam*100       !in cm
          do class=0,4
              dbh_class(class,pft)=(2.0*dbh(pft)+
     *                  (class*0.25*2.0*dbh(pft))-0.125*2.0*dbh(pft))   
              
              bt(class,pft)=(param1(pft)*dbh_class(class,pft)
     *                    +param2(pft))
              tau_c(class,pft)=2.9*(bt(class,pft)**2.0)

              height_class(class,pft)=2*height(pft)-
     *             (class*0.25*2.0*height(pft))+ 0.125*2.0*height(pft)   
              cl_t(class,pft)=height_class(class,pft)*crown(pft)  
           enddo !pseudo age class

          else

c           GRASS ALLOCATION
c           Distribute this year's production among leaves and fine roots
c           according to leaf to rootmass ratio [eqn (33)]
c           Relocation of C from one compartment to the other not allowed:
c           negative increment in either compartment transferred to litter

            lminc_ind=(bm_inc_ind-lm_ind(pft,1)/lmtorm+rm_ind(pft,1))/
     *        (1.0+1.0/lmtorm)
            rminc_ind=bm_inc_ind-lminc_ind

            if (lminc_ind.ge.0.0) then

c             Add killed roots (if any) to below-ground litter
c
c   CHECK: take out if statement because if rminc is negative than
c   root mass has been translocated to the leaves, therefore mass balance
c   problem since this carbon stays in the vegetation but is in addition
c   added to the litter pool. ALLOW translocation from roots to leaves
c   i.e. assume carbon stores in the roots which can be delivered
c   to the leaves under times of stress.
c
*
*              if (rminc_ind.lt.0.0)
*     *          litter_bg(pft)=litter_bg(pft)+(-rminc_ind)*nind(pft)
*
            else

c             Negative allocation to leaf mass

              rminc_ind=bm_inc_ind
              lminc_ind=(rm_ind(pft,1)+rminc_ind)*lmtorm-lm_ind(pft,1)
                          !from eqn (9)

c             Add killed leaves to litter

              temp=litter_ag_leaf(pft,1) ! Doug 11/12
              litter_ag_leaf(pft,1)=litter_ag_leaf(pft,1)
     *             +(-lminc_ind)*nind(pft)
              if (litter_ag_leaf(pft,1) .gt. 0.0) then
              do nc=2,nco2
                 litter_ag_leaf(pft,nc)=(litter_ag_leaf(pft,nc)*temp
     *                +(-lminc_ind)*lm_ind(pft,nc)*nind(pft))
     *                /litter_ag_leaf(pft,1)
              enddo
              endif

c            KIRSTEN: fuel classes
             
              temp=fuel_1hr_leaf(pft,1)
              fuel_1hr_leaf(pft,1)=fuel_1hr_leaf(pft,1)
     *             +(-lminc_ind)*nind(pft)
              if (fuel_1hr_leaf(pft,1) .gt. 0.) then
              do nc=2,nco2
                 fuel_1hr_leaf(pft,nc)=(fuel_1hr_leaf(pft,nc)*temp
     *                +(-lminc_ind)*lm_ind(pft,nc)*nind(pft))
     *                /fuel_1hr_leaf(pft,1)
              enddo
              endif

            endif

c           Increment C compartments
            lm_inc(pft)  =lminc_ind*nind(pft)
            lm_ind(pft,1)=lm_ind(pft,1)+lminc_ind
            rm_ind(pft,1)=rm_ind(pft,1)+rminc_ind

            if (lminc_ind.gt.0.0) then
               do nc=2,nco2
                  lm_ind(pft,nc)=(lm_ind(pft,nc)*lm_temp
     *                 +lminc_ind*bm_inc(pft,nc))/lm_ind(pft,1)
               enddo
            endif

            if (rminc_ind.gt.0.0) then
               do nc=2,nco2
                  rm_ind(pft,nc)=(rm_ind(pft,nc)*rm_temp
     *                 +rminc_ind*bm_inc(pft,nc))/rm_ind(pft,1)
               enddo
            endif

          endif

c         Update LAI and FPC
C		  Doug 12/12: Using Bern grass "fix" provided by Beni Stocker
C		  			  In old version, indicidual grass lpf grew as if
c                     without compeition from other grass. This fixes
c                     that with scling by total grass cover.
C---------------------------------------------------------------
C                     OLD VERSION
C          if (crownarea(pft).gt.0.0) then
C            lai_ind(pft)=(lm_ind(pft,1)*sla(pft))/crownarea(pft)
C          else
C            lai_ind(pft)=0.0
C          endif
C 
C          fpc_ind=(1.0-exp(-0.5*lai_ind(pft)))
C          fpc_grid_old=fpc_grid(pft)
C          fpc_grid(pft)=crownarea(pft)*nind(pft)*fpc_ind
C---------------------------------------------------------------
C                     NEW VERSION
          lm_tot(pft)=0.0
 
          IF (tree(pft)) THEN
            lm_tot(pft)=lm_ind(pft,1)
          ELSE ! Grass
            DO ppft=1,npft
              IF (.NOT.tree(ppft)) THEN
                lm_tot(pft)=lm_tot(pft)+lm_ind(ppft,1)
              END IF
            END DO
          END IF
 		  
          IF (crownarea(pft).GT.0.0) THEN
            lai_ind(pft)=(lm_tot(pft)*sla(pft))/crownarea(pft)
          ELSE
            lai_ind(pft)=0.0
          END IF
 		  
          fpc_ind=1.0-exp(-0.5*lai_ind(pft))		  
          fpc_grid_old=fpc_grid(pft)
          fpc_grid(pft)=crownarea(pft)*nind(pft)*fpc_ind
 		  
          IF (lm_tot(pft).GT.0.0) THEN
            fpc_grid(pft)=fpc_grid(pft)*lm_ind(pft,1)/lm_tot(pft)
          ELSE
            fpc_grid(pft)=0.0
          END IF
c---------------------------------------------------------------
 
          fpc_inc(pft)=max(0.0,fpc_grid(pft)-fpc_grid_old)
          if(fpc_grid(pft).eq.0.0) present(pft)=.false.
        endif
        DO nc=1,nco2						!Doug 01/09: MFL
          IF (fuel_1hr_leaf(pft,1)-temp_fuel(1,1)>0) THEN
            fuel_1hr_leaf_inc_pos(pft,nc,:)=
     *         fuel_1hr_leaf_inc_pos(pft,nc,:)+
     *        (fuel_1hr_leaf(pft,nc)-temp_fuel(nc,1))/12
          ELSE
            fuel_1hr_leaf_inc_neg(pft,nc,:)=
     *         fuel_1hr_leaf_inc_neg(pft,nc,:)+
     *        (fuel_1hr_leaf(pft,nc)-temp_fuel(nc,1))/12
          END IF
		  
          IF (fuel_1hr_wood(pft,1)-temp_fuel(1,2)>0) THEN
            fuel_1hr_wood_inc_pos(pft,nc,:)=
     *         fuel_1hr_wood_inc_pos(pft,nc,:)+
     *        (fuel_1hr_wood(pft,nc)-temp_fuel(nc,2))/12
          ELSE
            fuel_1hr_wood_inc_neg(pft,nc,:)=
     *         fuel_1hr_wood_inc_neg(pft,nc,:)+
     *        (fuel_1hr_wood(pft,nc)-temp_fuel(nc,2))/12
          END IF

          fuel_10hr_inc(pft,nc,:)=fuel_10hr_inc(pft,nc,:)+
     *      (fuel_10hr(pft,nc)-temp_fuel(nc,3))/12
          fuel_100hr_inc(pft,nc,:)=fuel_100hr_inc(pft,nc,:)+
     *      (fuel_100hr(pft,nc)-temp_fuel(nc,4))/12
          fuel_1000hr_inc(pft,nc,:)=fuel_1000hr_inc(pft,nc,:)+
     *      (fuel_1000hr(pft,nc)-temp_fuel(nc,5))/12
        END DO							!Doug 01/09: MFL


      enddo 	!pft

      return
      end	!allocation



c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE LIGHT
c     Competition for light among PFTs
c
c     Doug 01/09: addition of fuel_1hr_inc as I/O

      subroutine light(present,tree,lm_ind,sm_ind,hm_ind,rm_ind,
     *  crownarea,fpc_grid,fpc_inc,nind,litter_ag_leaf,litter_ag_wood,
     *  litter_bg,fuel_1hr_leaf,fuel_1hr_wood,
     *  fuel_10hr,fuel_100hr,fuel_1000hr,
     *  fuel_1hr_leaf_inc_pos,fuel_1hr_leaf_inc_neg, 
     *  fuel_1hr_wood_inc_pos,fuel_1hr_wood_inc_neg, 
     *  fuel_10hr_inc, fuel_100hr_inc, fuel_1000hr_inc, sla, year)

      implicit none

c     PARAMETERS
      integer npft,npftpar,nsoilpar
        parameter (npft=13,npftpar=59,nsoilpar=7)
      integer nco2
        parameter (nco2=3)           !number of C variables: C,13C,14C
      real fpc_tree_max
        parameter (fpc_tree_max=0.95)  !maximum total tree FPC

c     ARGUMENTS
      logical present(1:npft)
      logical tree(1:npft)
      real lm_ind(1:npft,1:nco2),sm_ind(1:npft,1:nco2)
      real hm_ind(1:npft,1:nco2),rm_ind(1:npft,1:nco2)
      real crownarea(1:npft)
      real fpc_grid(1:npft)
      real fpc_inc(1:npft)
      real nind(1:npft)
      REAL litter_ag_leaf(1:npft,1:nco2)	!Doug 11/12
      REAL litter_ag_wood(1:npft,1:nco2)
      real litter_bg(1:npft,1:nco2)
      REAL fuel_1hr_leaf(1:npft,1:nco2)
      REAL fuel_1hr_wood(1:npft,1:nco2)
      REAL fuel_10hr(1:npft,1:nco2)
      real fuel_100hr(1:npft,1:nco2)
      real fuel_1000hr(1:npft,1:nco2)
      REAL fuel_1hr_leaf_inc_neg(1:npft,1:nco2,1:12)
      REAL fuel_1hr_leaf_inc_pos(1:npft,1:nco2,1:12)	!Doug 01/09: records variations
      REAL fuel_1hr_wood_inc_neg(1:npft,1:nco2,1:12)
      REAL fuel_1hr_wood_inc_pos(1:npft,1:nco2,1:12)	!Doug 01/09: records variations
      REAL fuel_10hr_inc(1:npft,1:nco2,1:12)	!in fuel for each month (for FIRE)
      REAL fuel_100hr_inc(1:npft,1:nco2,1:12)
      REAL fuel_1000hr_inc(1:npft,1:nco2,1:12)
      REAL temp_fuel(1:nco2,1:5)
      real sla(1:npft)
      integer year

c     LOCAL VARIABLES
      integer pft
      integer nc
      integer ntree        !no of tree PFTs currently present
      real fpc_tree_total  !total grid FPC for tree PFTs
      real fpc_grass_total  !total grid FPC for grass PFTs
      real grasscover      !grass PFT proportional cover ("crown area")
      real fpc_inc_tree    !this years total FPC increment for tree PFTs
      real excess          !tree FPC or grass cover to be reduced
      real nind_kill       !reduction in individual density to reduce tree FPC
                           !to permitted maximum (indiv/m2)
      real lm_kill         !reduction in grass PFT leaf mass to reduce grass
                           !cover to permitted maximum (gC)
      real rm_kill         !reduction in grass PFT root mass to reduce grass
                           !cover to permitted maximum (gC)
      real fpc_ind
      real lai_ind
      real lm_old
      real ngrass
      real temp,litter_inc

      pft=0
      nc=0
      ntree=0
      fpc_tree_total=0.0
      fpc_grass_total=0.0
      grasscover=0.0
      fpc_inc_tree=0.0
      excess=0.0
      nind_kill=0.0
      lm_kill=0.0
      rm_kill=0.0
      fpc_ind=0.0
      lai_ind=0.0
      lm_old=0.0
      ngrass=0
      temp=0.0
      litter_inc=0.0

c     Calculate total woody FPC, FPC increment and grass cover (= crown area)

      do pft=1,npft

        if (present(pft)) then
          if (tree(pft)) then
            ntree=ntree+1
            fpc_tree_total=fpc_tree_total+fpc_grid(pft)
            fpc_inc_tree=fpc_inc_tree+fpc_inc(pft)
          else !grasses
            grasscover=crownarea(pft)
            fpc_grass_total=fpc_grass_total+fpc_grid(pft)
            ngrass=ngrass+1
          endif
        endif

      enddo

c     LIGHT COMPETITION
      do pft=1,npft

        temp_fuel(:,1)=fuel_1hr_leaf(pft,:)		!Doug 01/09: MFL
        temp_fuel(:,2)=fuel_1hr_wood(pft,:)		!Doug 01/09: MFL
        temp_fuel(:,3)=fuel_10hr(pft,:)			!Doug 01/09: MFL
        temp_fuel(:,4)=fuel_100hr(pft,:)		!Doug 01/09: MFL
        temp_fuel(:,5)=fuel_1000hr(pft,:)		!Doug 01/09: MFL

        if (present(pft)) then

           if (tree(pft)) then

              if (fpc_tree_total.gt.fpc_tree_max) then    ! case (1)

                   if (fpc_inc_tree.gt.0.0) then
                      excess=(fpc_tree_total-fpc_tree_max)*
     *                  (fpc_inc(pft)/fpc_inc_tree)
                   else
                      excess=(fpc_tree_total-fpc_tree_max)*
     *                  (1.0/real(ntree))
                   endif

c  Reduce individual density (and thereby gridcell-level biomass)
c  so that total tree FPC reduced to 'fpc_tree_max'

                   nind_kill=nind(pft)*(excess/fpc_grid(pft))
                   nind(pft)=nind(pft)-nind_kill

c                 Transfer lost biomass to litter
C 	Doug 11/12
                  temp=litter_ag_leaf(pft,1)
                  litter_ag_leaf(pft,1)=litter_ag_leaf(pft,1)+nind_kill*
     *                 lm_ind(pft,1)
                  IF (litter_ag_leaf(pft,1).gt.0.0) THEN
                     DO nc=2,nco2
                        litter_inc=(lm_ind(pft,1)*lm_ind(pft,nc))
                        litter_ag_leaf(pft,nc)=(litter_ag_leaf(pft,nc)*
     *                       temp+litter_inc*nind_kill)/
     *                       litter_ag_leaf(pft,1)
                     END DO
                  END IF
				  
				  temp=litter_ag_wood(pft,1)
                  litter_ag_wood(pft,1)=litter_ag_wood(pft,1)+nind_kill*
     *                 (sm_ind(pft,1)+hm_ind(pft,1))
                  IF (litter_ag_wood(pft,1).gt.0.0) THEN
                     DO nc=2,nco2
                        litter_inc=(sm_ind(pft,1)*sm_ind(pft,nc)+
     *                       hm_ind(pft,1)*hm_ind(pft,nc))
                        litter_ag_wood(pft,nc)=(litter_ag_wood(pft,nc)*
     *                       temp+litter_inc*nind_kill)/
     *                       litter_ag_wood(pft,1)
                     END DO
                  END IF
				  
c                  temp=litter_ag(pft,1)
c                  litter_ag(pft,1)=litter_ag(pft,1)+nind_kill*
c     *                 (lm_ind(pft,1)+sm_ind(pft,1)+hm_ind(pft,1))
c                  if (litter_ag(pft,1).gt.0.0) then
c                     do nc=2,nco2
c                        litter_inc=(sm_ind(pft,1)*sm_ind(pft,nc)+
c     *                       hm_ind(pft,1)*hm_ind(pft,nc)+
c     *                       lm_ind(pft,1)*lm_ind(pft,nc))
c                        litter_ag(pft,nc)=(litter_ag(pft,nc)*
c     *                       temp+litter_inc*nind_kill)/
c     *                       litter_ag(pft,1)
c                     enddo
c                  endif

c                 KIRSTEN: fuel classes	

                  temp=fuel_1hr_leaf(pft,1)
                  fuel_1hr_leaf(pft,1)=fuel_1hr_leaf(pft,1)+
     *              lm_ind(pft,1)*nind(pft)
                  IF (fuel_1hr_leaf(pft,1) .gt. 0.0) THEN
                    DO nc=2,nco2
                      litter_inc=lm_ind(pft,1)*lm_ind(pft,nc)
                      fuel_1hr_leaf(pft,nc)=(fuel_1hr_leaf(pft,nc)*temp
     *                 +litter_inc*nind(pft))/fuel_1hr_leaf(pft,1)
                    END DO
                  END IF	

                  temp=fuel_1hr_wood(pft,1)
                  fuel_1hr_wood(pft,1)=fuel_1hr_wood(pft,1)+
     *              (0.045*sm_ind(pft,1)+ 0.045*hm_ind(pft,1))*nind(pft)
                  IF (fuel_1hr_wood(pft,1) .gt. 0.0) THEN
                    DO nc=2,nco2
                      litter_inc=sm_ind(pft,1)*sm_ind(pft,nc)*0.045
     *                    +hm_ind(pft,1)*hm_ind(pft,nc)+0.045
                      fuel_1hr_wood(pft,nc)=(fuel_1hr_wood(pft,nc)*temp
     *                 +litter_inc*nind(pft))/fuel_1hr_wood(pft,1)
                    END DO
                  END IF
				  
c                  temp=fuel_1hr(pft,1)
c                  fuel_1hr(pft,1)=fuel_1hr(pft,1)+(lm_ind(pft,1)+
c     *              0.045*sm_ind(pft,1)+ 0.045*hm_ind(pft,1))*nind(pft)
c                  if (fuel_1hr(pft,1) .gt. 0.0) then
c                  do nc=2,nco2
c                     litter_inc=lm_ind(pft,1)*lm_ind(pft,nc)
c     *                    +sm_ind(pft,1)*sm_ind(pft,nc)*0.045
c     *                    +hm_ind(pft,1)*hm_ind(pft,nc)+0.045
c                     fuel_1hr(pft,nc)=(fuel_1hr(pft,nc)*temp
c     *                 +litter_inc*nind(pft))/fuel_1hr(pft,1)
c                  enddo
c                  endif
                  
                  temp=fuel_10hr(pft,1)
                  fuel_10hr(pft,1)=fuel_10hr(pft,1)+((0.075*
     *             sm_ind(pft,1)+0.075*hm_ind(pft,1))*nind(pft))
                  if (fuel_10hr(pft,1) .gt. 0.0) then
                  do nc=2,nco2
                     litter_inc=sm_ind(pft,1)*sm_ind(pft,nc)*0.075
     *                    +hm_ind(pft,1)*hm_ind(pft,nc)+0.075
                     fuel_10hr(pft,nc)=(fuel_10hr(pft,nc)*temp
     *                    +litter_inc*nind(pft))/fuel_10hr(pft,1)
                  enddo
                  endif

                  temp=fuel_100hr(pft,1)
                  fuel_100hr(pft,1)=fuel_100hr(pft,1)+((0.21*
     *              sm_ind(pft,1) +0.21*hm_ind(pft,1))*nind(pft))
                  if (fuel_100hr(pft,1) .gt. 0.0) then
                  do nc=2,nco2
                     litter_inc=sm_ind(pft,1)*sm_ind(pft,nc)*0.21
     *                    +hm_ind(pft,1)*hm_ind(pft,nc)+0.21
                     fuel_100hr(pft,nc)=(fuel_100hr(pft,nc)*temp
     *                    +litter_inc*nind(pft))/fuel_100hr(pft,1)
                  enddo
                  endif

                  temp=fuel_1000hr(pft,1)
                  fuel_1000hr(pft,1)=fuel_1000hr(pft,1)+((0.67
     *                 *sm_ind(pft,1)+0.67*hm_ind(pft,1))*nind(pft))
                  if (fuel_1000hr(pft,1) .gt. 0.0) then
                  do nc=2,nco2
                     litter_inc=sm_ind(pft,1)*sm_ind(pft,nc)*0.67
     *                    +hm_ind(pft,1)*hm_ind(pft,nc)+0.67
                     fuel_1000hr(pft,nc)=(fuel_1000hr(pft,nc)*temp
     *                    +litter_inc*nind(pft))/fuel_1000hr(pft,1)
                  enddo
                  endif

                  temp=litter_bg(pft,1)
                  litter_bg(pft,1)=litter_bg(pft,1)+nind_kill
     *               *rm_ind(pft,1)
                   if (litter_bg(pft,1).gt.0.0) then
                      do nc=2,nco2
                         litter_inc=rm_ind(pft,1)*rm_ind(pft,nc)
                         litter_bg(pft,nc)=(litter_bg(pft,nc)*temp+
     *                        nind_kill*litter_inc)/litter_bg(pft,1)
                      enddo
                   endif    
             endif
         else  ! grass
             if(fpc_grass_total.gt.
     *            (1.0-min(fpc_tree_total,fpc_tree_max))) then

c   grass competes with itself if total fpc exceeds 1
c  NEEDS COMMENTS !!!!!!????!!!!
c
                     excess=(min(fpc_tree_total,fpc_tree_max)
     *                      +fpc_grass_total-1.0)
     *                         *(fpc_grid(pft)/fpc_grass_total)
                     lm_old=lm_ind(pft,1)
                     lm_ind(pft,1)=-2.0*
     *                 alog(1.0-(fpc_grid(pft)-excess))/sla(pft)
                     lm_kill=lm_old-lm_ind(pft,1)

                     rm_kill=rm_ind(pft,1)*(lm_kill/lm_old)
                     rm_ind(pft,1)=rm_ind(pft,1)-rm_kill

c               Transfer lost biomass to litter

                    temp=litter_ag_leaf(pft,1)
                    litter_ag_leaf(pft,1)=litter_ag_leaf(pft,1)+lm_kill
                    IF (litter_ag_leaf(pft,1) .gt. 0.0) THEN
                      DO nc=2,nco2
                       litter_ag_leaf(pft,nc)=(litter_ag_leaf(pft,nc)
     *                     *temp+lm_ind(pft,nc)*lm_kill)/
     *                     litter_ag_leaf(pft,1)
                      END DO
                    END IF

c             KIRSTEN: 1hr fuel

                    temp=fuel_1hr_leaf(pft,1)
                    fuel_1hr_leaf(pft,1)=fuel_1hr_leaf(pft,1)+lm_kill
                    if (fuel_1hr_leaf(pft,1) .gt. 0.0) then
                    do nc=2,nco2
                       fuel_1hr_leaf(pft,nc)=(fuel_1hr_leaf(pft,nc)*temp
     *                     +lm_ind(pft,nc)*lm_kill)/fuel_1hr_leaf(pft,1)
                    enddo
                    endif

                    temp=litter_bg(pft,1)
                    litter_bg(pft,1)=litter_bg(pft,1)+rm_kill
                    if (litter_bg(pft,1) .gt. 0.0) then
                    do nc=2,nco2
                       litter_bg(pft,nc)=(litter_bg(pft,nc)*temp+
     *                      rm_ind(pft,nc)*rm_kill)/litter_bg(pft,1)
                    enddo
                    endif
              endif

          endif

c
c       update fpc (for establishment routine)
c
          if (crownarea(pft).gt.0.0) then
            lai_ind=(lm_ind(pft,1)*sla(pft))/crownarea(pft)
          else
            lai_ind=0.0
          endif
          fpc_ind=(1.0-exp(-0.5*lai_ind))
		  
          fpc_grid(pft)=crownarea(pft)*nind(pft)*fpc_ind
          if(fpc_grid(pft).eq.0.0) present(pft)=.false.
        endif    ! present

        DO nc=1,nco2						!Doug 01/09: MFL
          IF (fuel_1hr_leaf(pft,1)-temp_fuel(1,1)>0) THEN
            fuel_1hr_leaf_inc_pos(pft,nc,:)=
     *         fuel_1hr_leaf_inc_pos(pft,nc,:)+
     *        (fuel_1hr_leaf(pft,nc)-temp_fuel(nc,1))/12
          ELSE
            fuel_1hr_leaf_inc_neg(pft,nc,:)=
     *         fuel_1hr_leaf_inc_neg(pft,nc,:)+
     *        (fuel_1hr_leaf(pft,nc)-temp_fuel(nc,1))/12
          END IF
		  
          IF (fuel_1hr_wood(pft,1)-temp_fuel(1,2)>0) THEN
            fuel_1hr_wood_inc_pos(pft,nc,:)=
     *         fuel_1hr_wood_inc_pos(pft,nc,:)+
     *        (fuel_1hr_wood(pft,nc)-temp_fuel(nc,2))/12
          ELSE
            fuel_1hr_wood_inc_neg(pft,nc,:)=
     *         fuel_1hr_wood_inc_neg(pft,nc,:)+
     *        (fuel_1hr_wood(pft,nc)-temp_fuel(nc,2))/12
          END IF
		  
          fuel_10hr_inc(pft,nc,:)=fuel_10hr_inc(pft,nc,:)+
     *      (fuel_10hr(pft,nc)-temp_fuel(nc,3))/12
          fuel_100hr_inc(pft,nc,:)=fuel_100hr_inc(pft,nc,:)+
     *      (fuel_100hr(pft,nc)-temp_fuel(nc,4))/12
          fuel_1000hr_inc(pft,nc,:)=fuel_1000hr_inc(pft,nc,:)+
     *      (fuel_1000hr(pft,nc)-temp_fuel(nc,5))/12
        END DO							!Doug 01/09: MFL



      enddo  !pft

      return
      end	!light



c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE MORTALITY
c     Tree background and stress mortality
c
c     Doug 01/09: addition of fuel_1hr_inc as I/O

      subroutine mortality(pftpar,present,tree,boreal,bm_inc,
     *  turnover_ind,sla,lm_ind,sm_ind,hm_ind,rm_ind,nind,
     *  litter_ag_leaf,litter_ag_wood,
     *  litter_bg,fuel_1hr_leaf,fuel_1hr_wood,
     *  fuel_10hr,fuel_100hr,fuel_1000hr,
     *  fuel_1hr_leaf_inc_pos, fuel_1hr_leaf_inc_neg,
     *  fuel_1hr_wood_inc_pos, fuel_1hr_wood_inc_neg,
     *  fuel_10hr_inc,
     *  fuel_100hr_inc, fuel_1000hr_inc,
     *  dtemp,anpp,mtemp_max,year)

      implicit none

c     PARAMETERS
      integer npft,npftpar,nsoilpar
      parameter (npft=13,npftpar=59,nsoilpar=7)
      integer nco2
        parameter (nco2=3)           !number of C variables: C,13C,14C
      real mort_max
        parameter (mort_max=0.05) !asymptotic maximum mortality rate (/year)
      real k_mort                 !coefficient of growth efficiency in
        parameter (k_mort=0.5)   !mortality equation
c  ramp for heat damage function. Above 200 growing degree days above the
c  upper limit tw, establishment is zero and mortality 100%
      real ramp_gddtw
        parameter (ramp_gddtw=400.0)

c     ARGUMENTS
      real pftpar(1:npft,1:npftpar)
      logical present(1:npft),tree(1:npft),boreal(1:npft)
      real bm_inc(1:npft,1:nco2)
      real turnover_ind(1:npft)
      real sla(1:npft)
      real lm_ind(1:npft,1:nco2),sm_ind(1:npft,1:nco2)
      real hm_ind(1:npft,1:nco2),rm_ind(1:npft,1:nco2)
      real nind(1:npft)
      REAL litter_ag_leaf(1:npft,1:nco2) !Doug 11/12
      REAL litter_ag_wood(1:npft,1:nco2)
      real litter_bg(1:npft,1:nco2)
      REAL fuel_1hr_leaf(1:npft,1:nco2)
      REAL fuel_1hr_wood(1:npft,1:nco2)
      REAL fuel_10hr(1:npft,1:nco2)
      real fuel_100hr(1:npft,1:nco2)
      real fuel_1000hr(1:npft,1:nco2)
      REAL fuel_1hr_leaf_inc_pos(1:npft,1:nco2,1:12)
      REAL fuel_1hr_leaf_inc_neg(1:npft,1:nco2,1:12)	!Doug 01/09: records variations
      REAL fuel_1hr_wood_inc_pos(1:npft,1:nco2,1:12)
      REAL fuel_1hr_wood_inc_neg(1:npft,1:nco2,1:12)	!Doug 01/09: records variations
      REAL fuel_10hr_inc(1:npft,1:nco2,1:12)	!in fuel for each month (for FIRE)
      REAL fuel_100hr_inc(1:npft,1:nco2,1:12)
      REAL fuel_1000hr_inc(1:npft,1:nco2,1:12)
      REAL temp_fuel(1:nco2,5)
      real dtemp(1:365)
      real anpp(1:npft,1:nco2)
      real mtemp_max
      integer year

c     LOCAL VARIABLES
      real gddtw(1:npft)
      real twmax,greffic(1:npft)
      integer pft
      integer nc
      real bm_delta   !net individual living biomass increment (incorporating
                      !loss through leaf, root and sapwood turnover) (gC)
      real mort       !tree mortality rate
      real nind_kill  !reduction in individual density due to mortality
                      !(indiv/m2)
      integer d
      real heatstress(1:npft)  !reduction in individual density (& establishment) due to heat
                       !induced mortality  (indiv/m2)
      real temp,litter_inc

      gddtw(:)=0.0
      twmax=0.0
      greffic(:)=0.0
      pft=0
      nc=0
      bm_delta=0.0
      mort=0.0
      nind_kill=0.0
      d=0
      heatstress(:)=0.0
      temp=0.0
      litter_inc=0.0

      do pft=1,npft
        temp_fuel(:,1)=fuel_1hr_leaf(pft,:)		!Doug 01/09: MFL
        temp_fuel(:,2)=fuel_1hr_wood(pft,:)		!Doug 01/09: MFL
        temp_fuel(:,3)=fuel_10hr(pft,:)			!Doug 01/09: MFL
        temp_fuel(:,4)=fuel_100hr(pft,:)		!Doug 01/09: MFL
        temp_fuel(:,5)=fuel_1000hr(pft,:)		!Doug 01/09: MFL

c     initialisation
        heatstress(pft)=0.0
        twmax=pftpar(pft,31)  !PFT-specific upper limit of warmest-month temperature

        if (present(pft).and.tree(pft)) then
c         Calculate net individual living biomass increment

          bm_delta=max(0.0,bm_inc(pft,1)/nind(pft)-turnover_ind(pft))

c         Calculate growth efficiency (net biomass increment per unit leaf area)

          greffic(pft)=bm_delta/lm_ind(pft,1)/sla(pft)

c         Mortality rate inversely related to growth efficiency
c         (Prentice et al 1993)
          mort=mort_max/(1.0+k_mort*greffic(pft))
**          mort=0.0

c         heat damage mortality in boreal trees

          if (mtemp_max.gt.twmax) then  ! heat damage
c            calculate growing degree days above twmax
             gddtw(pft)=0.0
             do d=1,365
                gddtw(pft)=gddtw(pft)+max(0.0,dtemp(d)-twmax)
             enddo
             heatstress(pft)=min(1.0,gddtw(pft)/ramp_gddtw)
          else
             heatstress(pft)=0.0
          endif

c         Reduce individual density (and thereby gridcell-level biomass) by
c         mortality rate

          mort=min(1.0,mort+heatstress(pft))
          nind_kill=nind(pft)*mort
          nind(pft)=nind(pft)-nind_kill

c         Transfer lost biomass to litter

          temp=litter_ag_leaf(pft,1)
          litter_ag_leaf(pft,1)=litter_ag_leaf(pft,1)+nind_kill
     *         *lm_ind(pft,1)
          IF (litter_ag_leaf(pft,1) .gt. 0.0) THEN
            DO nc=2,nco2
             litter_inc=(lm_ind(pft,1)*lm_ind(pft,nc))
             litter_ag_leaf(pft,nc)=(litter_ag_leaf(pft,nc)*temp
     *            +litter_inc*nind_kill)/litter_ag_leaf(pft,1)
            END DO
          END IF
		  
          temp=litter_ag_wood(pft,1)
          litter_ag_wood(pft,1)=litter_ag_wood(pft,1)+nind_kill
     *         *(sm_ind(pft,1)+hm_ind(pft,1))
          IF (litter_ag_wood(pft,1) .gt. 0.0) THEN
            DO nc=2,nco2
             litter_inc=(sm_ind(pft,1)*sm_ind(pft,nc)
     *            +hm_ind(pft,1)*hm_ind(pft,nc))
             litter_ag_wood(pft,nc)=(litter_ag_wood(pft,nc)*temp
     *            +litter_inc*nind_kill)/litter_ag_wood(pft,1)
            END DO
          END IF
		  
		  
c          temp=litter_ag(pft,1)
c          litter_ag(pft,1)=litter_ag(pft,1)+nind_kill
c     *         *(sm_ind(pft,1)+hm_ind(pft,1)+lm_ind(pft,1))
c          if (litter_ag(pft,1) .gt. 0.0) then
c          do nc=2,nco2
c             litter_inc=(sm_ind(pft,1)*sm_ind(pft,nc)+hm_ind(pft,1)
c     *            *hm_ind(pft,nc)+lm_ind(pft,1)*lm_ind(pft,nc))
c             litter_ag(pft,nc)=(litter_ag(pft,nc)*temp
c     *            +litter_inc*nind_kill)/litter_ag(pft,1)
c          enddo
c          endif

c         KIRSTEN: fuel classes

c		  Doug 11/12: See commented code below for 2007 changes to isotope code
          temp=fuel_1hr_leaf(pft,1)

          fuel_1hr_leaf(pft,1)=fuel_1hr_leaf(pft,1)+
     *      lm_ind(pft,1)*nind_kill

          IF (fuel_1hr_leaf(pft,1) .gt. 0.0) THEN
            DO nc=2,nco2
              litter_inc=lm_ind(pft,1)*lm_ind(pft,nc)
              fuel_1hr_leaf(pft,nc)=(fuel_1hr_leaf(pft,nc)*temp
     *            +litter_inc*nind_kill)/fuel_1hr_leaf(pft,1)
            END DO
          END IF
		  
          temp=fuel_1hr_wood(pft,1)

          fuel_1hr_wood(pft,1)=fuel_1hr_wood(pft,1)+
     *         (0.045*sm_ind(pft,1)+ 0.045*hm_ind(pft,1))*nind_kill

          IF (fuel_1hr_wood(pft,1) .gt. 0.0) THEN
            DO nc=2,nco2
              litter_inc=
     *            +sm_ind(pft,1)*sm_ind(pft,nc)*0.045
     *            +hm_ind(pft,1)*hm_ind(pft,nc)+0.045
              fuel_1hr_wood(pft,nc)=(fuel_1hr_wood(pft,nc)*temp
     *            +litter_inc*nind_kill)/fuel_1hr_wood(pft,1)
            END DO
          END IF

c          temp=fuel_1hr(pft,1)
c DM  CHANGE to same as CO2_NO_ISOTOPE on 08 JULY 2007
c          fuel_1hr(pft,1)=fuel_1hr(pft,1)+(lm_ind(pft,1)+
c     *         0.045*sm_ind(pft,1)+ 0.045*hm_ind(pft,1))*nind(pft)
c
c          fuel_1hr(pft,1)=fuel_1hr(pft,1)+(lm_ind(pft,1)+
c     *         0.045*sm_ind(pft,1)+ 0.045*hm_ind(pft,1))*nind_kill
c DM  END CHANGE
c
c          if (fuel_1hr(pft,1) .gt. 0.0) then
c          do nc=2,nco2
c DM  CHANGE to same (?) as CO2_NO_ISOTOPE on 08 JULY 2007
c             litter_inc=lm_ind(pft,1)*lm_ind(pft,nc)
c     *            +sm_ind(pft,1)*sm_ind(pft,nc)*0.045
c     *            +hm_ind(pft,1)*hm_ind(pft,nc)+0.045
c             fuel_1hr(pft,nc)=(fuel_1hr(pft,nc)*temp
c     *            +litter_inc*nind(pft))/fuel_1hr(pft,1)
c             litter_inc=lm_ind(pft,1)*lm_ind(pft,nc)
c     *            +sm_ind(pft,1)*sm_ind(pft,nc)*0.045
c     *            +hm_ind(pft,1)*hm_ind(pft,nc)+0.045
c             fuel_1hr(pft,nc)=(fuel_1hr(pft,nc)*temp
c     *            +litter_inc*nind_kill)/fuel_1hr(pft,1)
c DM  END CHANGE
c          enddo
c          endif       

          temp=fuel_10hr(pft,1)
c DM  CHANGE to same as CO2_NO_ISOTOPE on 08 JULY 2007
c          fuel_10hr(pft,1)=fuel_10hr(pft,1)+((0.075*
c     *         sm_ind(pft,1)+0.075*hm_ind(pft,1))*nind(pft))
          fuel_10hr(pft,1)=fuel_10hr(pft,1)+((0.075*
     *         sm_ind(pft,1)+0.075*hm_ind(pft,1))*nind_kill)
c DM  END CHANGE
          if (fuel_10hr(pft,1) .gt. 0.0) then
          do nc=2,nco2
             litter_inc=sm_ind(pft,1)*sm_ind(pft,nc)*0.075
     *            +hm_ind(pft,1)*hm_ind(pft,nc)+0.075
c DM  CHANGE to same (?) as CO2_NO_ISOTOPE on 08 JULY 2007
c             fuel_10hr(pft,nc)=(fuel_10hr(pft,nc)*temp
c     *            +litter_inc*nind(pft))/fuel_10hr(pft,1)
             fuel_10hr(pft,nc)=(fuel_10hr(pft,nc)*temp
     *            +litter_inc*nind_kill)/fuel_10hr(pft,1)
c DM  END CHANGE
          enddo
          endif

          temp=fuel_100hr(pft,1)
c DM  CHANGE to same as CO2_NO_ISOTOPE on 08 JULY 2007
c          fuel_100hr(pft,1)=fuel_100hr(pft,1)+((0.21*
c     *         sm_ind(pft,1) +0.21*hm_ind(pft,1))*nind(pft))
          fuel_100hr(pft,1)=fuel_100hr(pft,1)+((0.21*
     *         sm_ind(pft,1) +0.21*hm_ind(pft,1))*nind_kill)
c DM  END CHANGE
          if (fuel_100hr(pft,1) .gt. 0.0) then
          do nc=2,nco2
             litter_inc=sm_ind(pft,1)*sm_ind(pft,nc)*0.21
     *            +hm_ind(pft,1)*hm_ind(pft,nc)+0.21
c DM  CHANGE to same (?) as CO2_NO_ISOTOPE on 08 JULY 2007
c             fuel_100hr(pft,nc)=(fuel_100hr(pft,nc)*temp
c     *            +litter_inc*nind(pft))/fuel_100hr(pft,1)
             fuel_100hr(pft,nc)=(fuel_100hr(pft,nc)*temp
     *            +litter_inc*nind_kill)/fuel_100hr(pft,1)
c DM  END CHANGE
          enddo
          endif

          temp=fuel_1000hr(pft,1)
c DM  CHANGE to same as CO2_NO_ISOTOPE on 08 JULY 2007
c          fuel_1000hr(pft,1)=fuel_1000hr(pft,1)+((0.67
c     *         *sm_ind(pft,1)+0.67*hm_ind(pft,1))*nind(pft))
          fuel_1000hr(pft,1)=fuel_1000hr(pft,1)+((0.67
     *         *sm_ind(pft,1)+0.67*hm_ind(pft,1))*nind_kill)
c DM  END CHANGE
          if (fuel_1000hr(pft,1) .gt. 0.0) then
          do nc=2,nco2
             litter_inc=sm_ind(pft,1)*sm_ind(pft,nc)*0.67
     *            +hm_ind(pft,1)*hm_ind(pft,nc)+0.67
c DM  CHANGE to same (?) as CO2_NO_ISOTOPE on 08 JULY 2007
c             fuel_1000hr(pft,nc)=(fuel_1000hr(pft,nc)*temp
c     *            +litter_inc*nind(pft))/fuel_1000hr(pft,1)
             fuel_1000hr(pft,nc)=(fuel_1000hr(pft,nc)*temp
     *            +litter_inc*nind_kill)/fuel_1000hr(pft,1)
c DM  END CHANGE
          enddo
          endif


          temp=litter_bg(pft,1)
          litter_bg(pft,1)=litter_bg(pft,1)+nind_kill*rm_ind(pft,1)
          if (litter_bg(pft,1) .gt. 0.0) then
          do nc=2,nco2
             litter_bg(pft,nc)=(litter_bg(pft,nc)*temp+nind_kill
     *            *rm_ind(pft,1)*rm_ind(pft,nc))/litter_bg(pft,1)
          enddo
          endif
        endif

        if (nind(pft).eq.0.) present(pft)=.false.

        DO nc=1,nco2                                            !Doug 01/09: MFL
          IF (fuel_1hr_leaf(pft,1)-temp_fuel(1,1)>0) THEN
            fuel_1hr_leaf_inc_pos(pft,nc,:)=
     *         fuel_1hr_leaf_inc_pos(pft,nc,:)+
     *        (fuel_1hr_leaf(pft,nc)-temp_fuel(nc,1))/12
          ELSE
            fuel_1hr_leaf_inc_neg(pft,nc,:)=
     *         fuel_1hr_leaf_inc_neg(pft,nc,:)+
     *        (fuel_1hr_leaf(pft,nc)-temp_fuel(nc,1))/12
          END IF
          IF (fuel_1hr_wood(pft,1)-temp_fuel(1,2)>0) THEN
            fuel_1hr_wood_inc_pos(pft,nc,:)=
     *         fuel_1hr_wood_inc_pos(pft,nc,:)+
     *        (fuel_1hr_wood(pft,nc)-temp_fuel(nc,2))/12
          ELSE
            fuel_1hr_wood_inc_neg(pft,nc,:)=
     *         fuel_1hr_wood_inc_neg(pft,nc,:)+
     *        (fuel_1hr_wood(pft,nc)-temp_fuel(nc,2))/12
          END IF
	  fuel_10hr_inc(pft,nc,:)=fuel_10hr_inc(pft,nc,:)+
     *      (fuel_10hr(pft,nc)-temp_fuel(nc,3))/12
          fuel_100hr_inc(pft,nc,:)=fuel_100hr_inc(pft,nc,:)+
     *      (fuel_100hr(pft,nc)-temp_fuel(nc,4))/12
          fuel_1000hr_inc(pft,nc,:)=fuel_1000hr_inc(pft,nc,:)+
     *      (fuel_1000hr(pft,nc)-temp_fuel(nc,5))/12
        END DO							!Doug 01/09: MFL


      enddo	!pft

      return
      end	!mortality
	  



c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE FIRE
c     Biomass destruction through disturbance by fire

      subroutine fire(year,start_year,                                       !Doug 07/09: add start_year
     *  pftpar,dtemp,dtemp_min,dtemp_max,dprec,
     *  dwindsp,dlightn,dphen,dphen_change,
     *  litter_ag_leaf,litter_ag_wood,
     *  litter_bg,fuel_1hr_leaf,fuel_1hr_wood,    !Doug 06/09: dphen_change added for fire paradox experiments		
     *  fuel_10hr, fuel_100hr,fuel_1000hr,
     *  pfuel_limit, !Doug 12/12
     *  fuel_1hr_leaf_inc_pos, fuel_1hr_leaf_inc_neg, 
     *  fuel_1hr_wood_inc_pos, fuel_1hr_wood_inc_neg, 
     *  fuel_10hr_inc,fuel_100hr_inc, fuel_1000hr_inc,
     *  fuel_1hr_del,fuel_10hr_del,fuel_100hr_del,fuel_1000hr_del,
     *  mfuel_1hr_total,
     *  mfuel_1hr_leaf_total,mfuel_1hr_wood_total,
     *	mfuel_10hr_total, mfuel_100hr_total,           !Doug 03/09: outputs for monthly fuel loads
     *  mfuel_1000hr_total, mlivegrass,                                  !Doug 05/09: mlivegrass added
     *  acflux_fire,mcflux_fire,afire_frac,lm_ind,rm_ind,
     *  sm_ind,hm_ind,nind,dw1,present,tree,lat,mw1,fpc_grid, popden,
     *  a_nd,height,height_class,dbh,dbh_class,tau_c,cl_t, !Doug 02/13: dbh_class added
     *  BTparam1,BTparam2,BTmode0,num_fire,annum_fire,
     *  area_burnt,an_areafires,mfdi,an_fdi,an_fseason,mcflux_trace,
     *  acflux_trace,m_fc_crown,an_fc_crown,m_i_surface,an_i_surface,
     *  dhuman_ign,         ! human_ignition
     *  num_fire_human,num_fire_lightn,
     *  annum_fire_human,annum_fire_lightn,
     *  area_burnt_human,area_burnt_lightn,
     *  an_areafires_human,an_areafires_lightn,
     *  afire_frac_human,afire_frac_lightn,char_net_fuel_0,
     * livegrass_0,dead_fuel_0,dead_fuel_all_0,
     * fuel_all_0,fuel_1hr_total_0,fuel_10hr_total_0,
     * fuel_100hr_total_0,fuel_1000hr_total_0,
     * mfire_frac,lon,
     * crop,pas,fbdep,ni_acc,afire_frac_afap_old,
     * dlm_1hr_old,dlm_10hr_old,dlm_100hr_old,dlm_1000hr_old,mlm,
     * dpet,aprec)	 
      implicit none

c     PARAMETERS
      integer npft,npftpar,nsoilpar
         parameter (npft=13,npftpar=59,nsoilpar=7)
      real p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12
         parameter (p1=0.17,p2=0.10,p3=0.04,p4=0.01,p5=0.01
     *              ,p6=0.01,p7=0.03,p8=0.04,p9=0.09,p10=0.12,
     *              p11=0.16,p12=0.22)
 
       integer nco2
        parameter (nco2=3)           !number of C variables: C,13C,14C
      real pi
         parameter (pi=3.14159265)
      real minfuel
         parameter (minfuel=100.0)  !fuel threshold to carry a fire (gC/m2)
c     Reg-FIRM
      real wind_speed
*         parameter (wind_speed=76.8)    !wind speed (m/min), =1.28 m/s
c     total mineral content
      real MINER_TOT
         parameter(MINER_TOT=0.055)
c     calorific heat content (kJ/kg)
      real H
         parameter(H=18000.0)
c     surface-area-to-volume ratio (cm/cm?³)
      real sigma_1hr,sigma_10hr,sigma_100hr
         parameter(sigma_1hr=66.0,sigma_10hr=3.58,sigma_100hr=0.98)
c     ALLAN
       real sigma_1000hr
       parameter(sigma_1000hr = 0.5) !highly subjective
       real sigma_livegrass
       parameter (sigma_livegrass=80.0)
c       Implement formula for Me after P&R(1986). Later, weight by dead fuel classes.

c      discussion 25/07/05: physical vs. chemical composition of dead fuel on flammability.
c      With this approach no PFT influence on fire ignition, only physical characteristics
c             of dead fuel
c      Check this approach vs. me as a PFT parameter.
       real moistfactor_livegrass,moistfactor_1hr, moistfactor_10hr
       real moistfactor_100hr,moistfactor_1000hr
c       parameter (moistfactor_livegrass = 0.398)   !0.524-0.066*log(sigma_livegrass)
       parameter (moistfactor_1hr = 0.1)         !0.524-0.066*log(sigma_1hr)
       parameter (moistfactor_10hr = 0.15)        !0.524-0.066*log(sigma_10hr)
       parameter (moistfactor_100hr = 0.25)      !0.524-0.066*log(sigma_100hr)
       parameter (moistfactor_1000hr = 0.3)    !0.524-0.066*log(sigma_1000hr)

       real part_dens
       parameter(part_dens=513.0)
       real pack_ratio,pack_ratio_1hr,pack_ratio_10hr
       real pack_ratio_dead,pack_ratio_100hr
       real pack_ratio_lg

       parameter (pack_ratio_1hr=0.003)
       parameter (pack_ratio_10hr=0.03)
       parameter (pack_ratio_100hr=0.025)
       parameter (pack_ratio_lg=0.002)

c     fire duration (min)
      real fire_durat
c         parameter(fire_durat=420.0)   !8 hrs.
c       parameter (fire_durat=240.0) !4 x 60 mins, Allan 
c     LPJ limited to daily timestep.  This is assumed to be the  
c     maximum fire duration during the 
c     mid-high latitude summer or sub- /equatorial dry season.

       real fbd_a  !scalar used to transform FBD of 10hr fuels to 1hr equivs
       real fbd_b  !scalar used to transform FBD of 100hr fuels to 1hr equivs
       parameter(fbd_a = 1.2, fbd_b = 1.4)

       real fbd_C3_livegrass,fbd_C4_livegrass
       parameter (fbd_C3_livegrass=4.0,fbd_C4_livegrass=4.0) !kg/m3

c     ARGUMENTS
      integer year
      integer start_year !Doug 07/09
      real pftpar(1:npft,1:npftpar)
      real dtemp(1:365),dtemp_min(1:365),dtemp_max(1:365)
      real dprec(1:365),dwindsp(1:365)
      real dhuman_ign(1:365)      ! Yan
      real dlightn(1:365)              !number of lightning fires per mln.ha
      integer dlightn_cons		 !Doug 12/08: coverting number of fires n dlightn to a constant strike count for the year
      real dphen(1:365,1:npft)
      REAL dphen_change(1:12,1:npft)     !Doug 06/09: for fire paradox experiments
      LOGICAL test_growing_start !Doug 06/09: Tests to see if the start of the growing season has happened yet
      REAL litter_ag_leaf(1:npft,1:nco2)
      REAL litter_ag_wood(1:npft,1:nco2)
      REAL litter_bg(1:npft,1:nco2)
      REAL fuel_1hr_leaf(1:npft,1:nco2)
      REAL fuel_1hr_wood(1:npft,1:nco2)
      REAL f1hr_leaf_frac(1:npft)	!Doug 11/12: fraction of 1hr fuel that is grass and leaf
      REAL fuel_10hr(1:npft,1:nco2)
      real fuel_100hr(1:npft,1:nco2)			
      real fuel_1000hr(1:npft,1:nco2)

      REAL fuel_1hr_leaf_inc_pos(1:npft,1:nco2,1:12)
      REAL fuel_1hr_leaf_inc_neg(1:npft,1:nco2,1:12)	!Doug 01/09: records variations
      REAL fuel_1hr_wood_inc_pos(1:npft,1:nco2,1:12)
      REAL fuel_1hr_wood_inc_neg(1:npft,1:nco2,1:12)	!Doug 01/09: records variations
      REAL fuel_10hr_inc(1:npft,1:nco2,1:12)	!in fuel for each month (for FIRE)
      REAL fuel_100hr_inc(1:npft,1:nco2,1:12)
      REAL fuel_1000hr_inc(1:npft,1:nco2,1:12)

      REAL dfuel_leaf(1:365,1:npft)             !Doug 03/09: redistribution of 1 hour
                                                !fuel increments from monthly to daily 
                                                !as a function of pft leaf phenology
      REAL leaf_decayc(1:npft)                  !Doug 03/09: leaf longliverty, for
                                                !calculating daily leaf turnover into
                                                !1 hour fuel

      REAL fuel_1hr_del(1:npft,1:nco2)          !Doug 01/09 MFL1: record annual change
      REAL fuel_10hr_del(1:npft,1:nco2)
      REAL fuel_100hr_del(1:npft,1:nco2)
      REAL fuel_1000hr_del(1:npft,1:nco2)
                                               
      REAL dfuel_update_1hr_leaf(1:npft,1:nco2)             !Doug 03/09: store the daily update of xhr fuel.
      REAL dfuel_update_1hr_wood(1:npft,1:nco2)
      REAL dfuel_update_10hr(1:npft,1:nco2)            !Doug 05/09
      REAL dfuel_update_100hr(1:npft,1:nco2)           !Doug 05/09
      REAL dfuel_update_1000hr(1:npft,1:nco2)          !Doug 05/09

      real acflux_fire(1:nco2)
      real mcflux_fire(1:12,1:nco2)
      real afire_frac
      real lm_ind(1:npft,1:nco2),rm_ind(1:npft,1:nco2)
      real sm_ind(1:npft,1:nco2),hm_ind(1:npft,1:nco2)
      REAL nind_resp(1:npft),prop_fa_rs(1:npft) !Doug 02/13: number of resprouting individuals, and propotion of fire affected plants taht re-sprout per pft
      real nind(1:npft)
      real dw1(1:365)            !dail litter moisture from Nesterov inde
      logical present(1:npft),tree(1:npft)
      real lat,lon
      real mw1(1:12)             !dail litter moisture from Nesterov inde
      real fpc_grid(1:npft)
      real popden,a_nd
      real height(1:npft)       !tree height in m, stemdiam per pft
      real height_class(0:4,1:npft)
      REAL dbh(1:npft), dbh_class(0:4,1:npft)!tree height in m, stemdiam per pft !Doug 02/13: dbh_class added
      real tau_c(0:4,1:npft)                !critical time to cambial kill
      real cl_t(0:4,1:npft)
      REAL BTparam1(1:npft,1:3)        ! Doug 02/13: lower, medium and upper bound describing 1-50-99% quantile p1 in Bark thickness equation BT=p1+p2*DBH
      REAL BTparam2(1:npft,1:3)        ! Doug 02/13: "" for p2
      REAL BTmode0(1:npft,1:2)         ! Doug 02/13: 50% quatile starting point incase of establishment or tree death.
      real num_fire(1:12),annum_fire
      real area_burnt(1:12),an_areafires
      real mfdi(1:12),an_fdi
      real an_fseason
      real mcflux_trace(1:12,1:6)
      real acflux_trace(1:6)
      real m_fc_crown(1:12)
      real an_fc_crown
      real m_i_surface(1:12),an_i_surface  ! surface fire intensity


c     LOCAL VARIABLES
      integer m,month(1:12),month_length(1:12)
      data (month(m),m=1,12) /31,59,90,120,151,181,
     *            212,243,273,304,334,365/
      data (month_length(m),m=1,12)
     *               /31,28,31,30,31,30,31,31,30,31,30,31/

      integer monthno, midday(1:13)					!Doug 01/09: MFL
        data (midday(monthno),monthno=1,13)
     *    / 16,44,75,105,136,166,197,228,258,289,319,350,381 /


      real fuel,fire_prob
      real fire_index
      real fire_length
      real disturb

      real mcflux_fire_pft(1:12,1:npft,1:nco2)
      real mcflux_trace_pft(1:12,1:6,1:npft)
      real dlm(1:365) !dail litter moisture from Nesterov inde
      REAL mlm(1:12)
      real fire_frac(1:365)

c     Reg_FIRM
      real mtemp(1:12),mtemp_dmin(1:12),mtemp_dmax(1:12)
      real d_numfire(1:365),d_numf_old(1:365)
      real d_area_burnt(1:365)
      real d_fdi(1:365)

      integer pft,d,n,i,x,count
      integer count_int,count_fdi,count_yr,n_pft
      real disturb_flux(1:nco2)
      REAL fuel_1hr_total	  
      REAL fuel_1hr_leaf_total
      REAL fuel_1hr_wood_total
      REAL fuel_10hr_total
      real fuel_100hr_total
      real fuel_1000hr_total
      REAL pfuel_limit							!Doug 12/12: proportion of days of fuel limitation
      REAL fuel_left_minus(1:npft,1:5)          !Doug 01/09 MFL


      REAL mfuel_1hr_total(1:12)                !Doug 03/09: stores monthly average
      REAL mfuel_1hr_leaf_total(1:12)
      REAL mfuel_1hr_wood_total(1:12)
      REAL mfuel_10hr_total(1:12)               !values of fuel for output
      REAL mfuel_100hr_total(1:12)              !Doug 03/09
      REAL mfuel_1000hr_total(1:12)             !Doug 03/09
      REAL mlivegrass(1:12)

      real ef_trace(1:npft,1:6)
      real acflux_fire_pft(1:npft,1:nco2)
      real dcflux_fire_pft(1:365,1:npft,1:nco2)
      real resist(1:npft)
      real moistfactor,me,litter_ag_total
      real fire_term,area, area_ha
c     Reg-FIRM
      real m_firelength(1:12),length(0:12)
c      real dtemp_min,dtemp_max
      real U_back,U_front,gamma

      real ros_f(1:365),ros_b(1:365)
      real df(1:365),db(1:365)
      real FDI                         !fire risk function, unitless between 0. and 1.
      real human_ign                   !number of human induced fires per mln.ha (subroutine)
      real fpc_tree_total,fpc_grass_total,lb_grass,lb
      real dens_fuel(1:npft),dens_fuel_tot,dens_fuel_ave,sigma
      real net_fuel,dead_fuel,livegrass
      real wind_forward,wind_backward,back_ws
      real d_i_surface(1:365)  ! surface fire intensity
      real fuel_consum,an_fuel_consum
*      real pot_fc_1hr(1:npft),pot_fc_10hr(1:npft),pot_fc_100hr(1:npft)
      real pot_fc_lg(1:npft,1:nco2)
      real ratio_lg_reduc
      real fc_1hr_total,fc_10hr_total,fc_100hr_total,fc_lg_total
      REAL fc_1hr_leaf(1:npft,1:nco2)
      REAL fc_1hr_wood(1:npft,1:nco2)
      REAL fc_10hr(1:npft,1:nco2)
      real fc_100hr(1:npft,1:nco2)
      real fc_lg(1:npft,1:nco2),fc_1000hr(1:npft,1:nco2)
      real sh(1:npft),f(1:npft),g(1:npft) !scorch height
      real crown(1:npft),ck(1:npft),cl(1:npft)
      real ck_t(1:5,1:npft)  ! crown kill in %  etc.
      real fc_crown(1:365)
      real wc_live(1:npft),postf_mort(1:npft)
      real cf(1:npft),tau_l(1:npft)
      real pm_tau(1:npft),pm_tau_class(0:4,1:npft)
      real pm_ck(1:npft),r_ck(1:npft),p(1:npft)
      real nind_fuel(1:npft),nind_old(1:npft),nind_kill(1:npft)

       real ratio_fbd(npft)   ! proportion of dead fuel in class i PFT j as a proportion of total dead fuel
                              !Doug remove: tun back to scalar
       real ratio_dead_fuel
       real ratio_live_fuel
       real dens_deadgrass_ave
       real dens_livegrass_ave
       real char_dens_fuel_ave
       real char_net_fuel
       real char_net_fuel_0
       real deadgrass
       real char_sigma
       real char_moistfactor
       real ratio_C3_livegrass
       real ratio_C4_livegrass
c       real moistfactor_C3_livegrass
c       real moistfactor_C4_livegrass
c       parameter (moistfactor_C3_livegrass = 0.1, 
c     *       moistfactor_C4_livegrass =0.1)

c       real moistfactor_1hr, moistfactor_10hr
c       real moistfactor_100hr, moistfactor_1000hr
       real   dlm_lg (1:365), dlm_1hr(1:365), dlm_10hr(1:365)  
       real   dlm_100hr (1:365), dlm_1000hr (1:365)
c       real fbd_C3_livegrass
       real nind_fa(1:npft) !# trees affected by fire
       real nind_nkill(1:npft) !# trees not killed by fire
       real ag_non_leaf_ind(1:npft,1:nco2) ! sm_ind + hm_ind
       REAL fuel_update_1hr_leaf(1:npft,1:nco2)
       REAL fuel_update_1hr_wood(1:npft,1:nco2)
       real fuel_update_10hr(1:npft,1:nco2)
       real fuel_update_100hr(1:npft,1:nco2)
       real fuel_update_1000hr(1:npft,1:nco2) ! used to tidy up book-keeping of fuel classes
       real prop_fa_nk(1:npft) ! propn of indivs fire affected but not killed
       real afire_frac_temp,an_areafires_temp
       integer class
       real pot_fc_lg_temp(1:npft)
       real ni_acc
       REAL dlm_1hr_old,dlm_10hr_old,dlm_100hr_old,dlm_1000hr_old
       REAL dpet(1:365)
       REAL aprec
       real temp
       real all_ind,litter_inc
       integer nc               ! nc is added for nco2
C Yan
       real d_numfire_human(1:365),d_numfire_lightn(1:365)
       real d_area_burnt_human(1:365),d_area_burnt_lightn(1:365)
       real d_i_surface_human(1:365),d_i_surface_lightn(1:365)
       real fire_frac_human,fire_frac_lightn 
       real afire_frac_human,afire_frac_lightn 
       real afire_frac_human_temp,afire_frac_lightn_temp 
       real an_areafires_human_temp,an_areafires_lightn_temp 
       real an_areafires_human,an_areafires_lightn 
       real num_fire_human(12),num_fire_lightn(12) 
       real annum_fire_human,annum_fire_lightn
       real m_i_surface_human(12),m_i_surface_lightn(12) 
       real an_i_surface_human,an_i_surface_lightn
       real area_burnt_human(12),area_burnt_lightn(12)
       REAL fuel_1hr_leaf_0(1:npft,1:nco2)
       REAL fuel_1hr_wood_0(1:npft,1:nco2)
       REAL fuel_10hr_0(1:npft,1:nco2)
       real fuel_100hr_0(1:npft,1:nco2),fuel_1000hr_0(1:npft,1:nco2)
       REAL fuel_1hr_leaf_left(1:npft,1:nco2)
       REAL fuel_1hr_wood_left(1:npft,1:nco2)
       REAL fuel_10hr_left(1:npft,1:nco2)
       real fuel_100hr_left(1:npft,1:nco2)
       real fuel_1000hr_left(1:npft,1:nco2)
       real livegrass_left(1:npft,1:nco2)
       real lm_ind_0(1:npft,1:nco2),nind_0(1:npft) ! Initial value at the beginning of the year.
       real pot_fc_lg_0(1:npft,1:nco2) ! This is the initial leaf mass for the live grass at the beginning of the year. 
       real char_dens_fuelbed
       real fbdep 
!       parameter (fbdep=0.18)
       integer co2
       real fbdepth

        real livegrass_0,dead_fuel_0,dead_fuel_all_0,fuel_all_0
        real fuel_1hr_total_0,fuel_10hr_total_0
        real fuel_100hr_total_0,fuel_1000hr_total_0
        real dlightn_sca(1:365)
        real dlightn_ave
        integer dd
        real d_NI,char_alpha_fuel
        real fbdep_dead_fuel,fbdep_live_fuel
         real p_ig,p_lcc,bet,q_ig
         integer an_lightn,dlightn_lcc(1:365)
         integer m_lightn(1:12) 
         integer day
        integer counter_fire
        
        real afire_frac_afap_old,mfire_frac(1:12)
       real crop
       real pas
        integer temp1,j,l,ll
         real ii,jj
         real d_numfire_ud(1:365)
         integer nran
        real ran
         real ran_num
        INTEGER testing
        REAL fuel_load_subtraction(1:npft,1:4)    !Doug 06/09: This is a cheat for the fire
                                                   !paradox experiments. If you see this, and your
                                                   !not me, then get rid of it cos it wont be
                                                   !needed anymore
        INTEGER ndaymonth(1:12)                   !Doug 10/09: number of days in
          DATA(ndaymonth(m),m=1,12)                !each month
     *     / 31,28,31,30,31,30,31,31,30,31,30,31 /

      popden=0.0
      a_nd=0.0

cYan  Initialise
      dhuman_ign(:)=0.0
      d_numfire_human(:)=0.0    ! Yan
      d_numfire_lightn(:)=0.0   ! Yan
	
      d_area_burnt_human(:)=0.0    ! Yan
      d_area_burnt_lightn(:)=0.0   ! Yan

      d_i_surface_human(:)=0.0    ! Yan
      d_i_surface_lightn(:)=0.0   ! Yan

       fire_frac_human=0.0    ! Yan
       fire_frac_lightn=0.0   ! Yan
       afire_frac_human=0.0    ! Yan
       afire_frac_lightn=0.0   ! Yan
       afire_frac_human_temp=0.0    ! Yan
       afire_frac_lightn_temp=0.0   ! Yan
      an_i_surface_human=0.0    ! Yan
      an_i_surface_lightn=0.0   ! Yan
        num_fire_human(:)=0.0   ! Yan
        num_fire_lightn(:)=0.0  ! Yan
        area_burnt_human(:)=0.0   ! Yan
        area_burnt_lightn(:)=0.0  ! Yan
        m_i_surface_human(:)=0.0        ! Yan
        m_i_surface_lightn(:)=0.0       ! Yan
      annum_fire_human=0.0             ! Yan
      annum_fire_lightn=0.0            ! Yan
      an_areafires_human=0.0          ! Yan
      an_areafires_lightn=0.0         ! Yan

c     Initialise

      wind_speed=0.0  !to comment out, when a parameter
      fire_durat=0.0  !to comment out, when a parameter

      fuel=0.0
      fire_prob=0.0
      fire_index=00
      fire_length=0.0
      disturb=0.0
      mcflux_fire_pft(:,:,:)=0.0
      mcflux_trace_pft(:,:,:)=0.0
      dlm(:)=0.0
      mlm(:)=0.0
      fire_frac(:)=0.0
      mtemp(:)=0.0
      mtemp_dmin(:)=0.0
      mtemp_dmax(:)=0.0
      d_numfire(:)=0.0
      d_numf_old(:)=0.0
      d_area_burnt(:)=0.0
      d_fdi(:)=0.0
      pft=0
      d=0
      n=0
      i=0
      x=0
      count=0
      count_int=0
      count_fdi=0
      count_yr=0
      n_pft=0
      disturb_flux(:)=0.0
      fuel_1hr_total=0.0
      fuel_1hr_leaf_total=0.0
      fuel_1hr_wood_total=0.0
      fuel_10hr_total=0.0
      fuel_100hr_total=0.0
      fuel_1000hr_total=0.0
      pfuel_limit=0.0
      ef_trace(:,:)=0.0
      acflux_fire_pft(:,:)=0.0
      dcflux_fire_pft(:,:,:)=0.0
      resist(:)=0.0
      moistfactor=0.0
      me=0.0
      litter_ag_total=0.0
      fire_term=0.0
      area=0.0
      area_ha=0.0
      m_firelength(:)=0.0
      length(:)=0.0
      U_back=0.0
      U_front=0.0
      gamma=0.0
      ros_f(:)=0.0
      ros_b(:)=0.0
      df(:)=0.0
      db(:)=0.0
      FDI=0.0
      fpc_tree_total=0.0
      fpc_grass_total=0.0
      lb_grass=0.0
      lb=0.0
      dens_fuel(:)=0.0
      dens_fuel_tot=0.0
      dens_fuel_ave=0.0
      sigma=0.0
      net_fuel=0.0
      dead_fuel=0.0
      livegrass=0.0
      wind_forward=0.0
      wind_backward=0.0
      back_ws=0.0
      d_i_surface(:)=0.0
      fuel_consum=0.0
      an_fuel_consum=0.0
      pot_fc_lg(:,:)=0.0
      ratio_lg_reduc=0.0
      fc_1hr_total=0.0
      fc_10hr_total=0.0
      fc_100hr_total=0.0
      fc_lg_total=0.0
      fc_1hr_leaf(:,:)=0.0
      fc_1hr_wood(:,:)=0.0
      fc_10hr(:,:)=0.0
      fc_100hr(:,:)=0.0
      fc_lg(:,:)=0.0
      fc_1000hr(:,:)=0.0
      sh(:)=0.0
      f(:)=0.0
      g(:)=0.0
      crown(:)=0.0
      ck(:)=0.0
      cl(:)=0.0
      ck_t(:,:)=0.0
      fc_crown(:)=0.0
      wc_live(:)=0.0
      postf_mort(:)=0.0
      cf(:)=0.0
      tau_l(:)=0.0
      pm_tau(:)=0.0
      pm_tau_class(:,:)=0.0
      pm_ck(:)=0.0
      r_ck(:)=0.0
      p(:)=0.0
      nind_fuel(:)=0.0
      nind_old(:)=0.0
      nind_kill(:)=0.0
       ratio_fbd(:)=0.0	!Doug: remove and turn back to scalar
       ratio_dead_fuel=0.0
       ratio_live_fuel=0.0
       dens_deadgrass_ave=0.0
       dens_livegrass_ave=0.0
       char_dens_fuel_ave=0.0
       char_net_fuel=0.0
       deadgrass=0.0
       char_sigma=0.0
       char_moistfactor=0.0
       ratio_C3_livegrass=0.0
       ratio_C4_livegrass=0.0
         dlm_lg (:)=0.0
        dlm_1hr(:)=0.0
         dlm_10hr(:)  =0.0
         dlm_100hr (:)=0.0
          dlm_1000hr (:)=0.0
       nind_fa(:)=0.0
       nind_nkill(:)=0.0
       nind_resp(:)=0.0 ! Doug 02/13: required for resprouters
       ag_non_leaf_ind(:,:)=0.0
       fuel_update_1hr_leaf(:,:)=0.0
       fuel_update_1hr_wood(:,:)=0.0
       fuel_update_10hr(:,:)=0.0
       fuel_update_100hr(:,:)=0.0
       fuel_update_1000hr(:,:)=0.0
       prop_fa_nk(:)=0.0
       prop_fa_rs(:)=0.0 ! Doug 02/13: required for resprouters
       afire_frac_temp=0.0
       an_areafires_temp=0.0
       class=0
       pot_fc_lg_temp(:)=0.0
c       ni_acc=0.0
       temp=0.0
       all_ind=0.0
       litter_inc=0.0
       nc=0


      an_i_surface=0.0
      an_fc_crown=0.0
        num_fire(:)=0.0
        area_burnt(:)=0.0
        mfire_frac(:)=0.0
        mfdi(:)=0.0
        m_fc_crown(:)=0.0
        m_i_surface(:)=0.0
        mcflux_fire(:,:)=0.0
          mcflux_trace(:,:)=0.0
       acflux_trace(:)=0.0
      acflux_fire(:)=0.0
      annum_fire=0.0
      an_areafires=0.0
      an_fdi=0.0
      an_fseason=0.0
          an_lightn=0
          dlightn_lcc(:)=0
          m_lightn(:)=0
           counter_fire=0
	afire_frac_afap_old=0.0


c     Assign a minimum fire fraction (for presentational purposes)

      afire_frac=0.001

c    ASSIGN PFT PARAMETER REQUIRED IN THE FIRE ROUTINE

c    ASSIGN PFT PARAMETER FOR Glob-FIRM AND Reg-FIRM

c          Calculate total above-ground litter

      DO pft=1,npft
        litter_ag_total=litter_ag_total+litter_ag_leaf(pft,1)
     *    +litter_ag_wood(pft,1)
		  
      END DO
      f1hr_leaf_frac(pft)=fuel_1hr_leaf(pft,1)/
     *  (fuel_1hr_leaf(pft,1)+fuel_1hr_wood(pft,1))

c        Calculate litter moisture weighting factor (moisture of extinction me)

c        Assign emission factors for trace gas emissions resulting from acflux_fire_pft
      do pft=1,npft
         ef_trace(pft,1)=pftpar(pft,38)/1000.0   !CO2
         ef_trace(pft,2)=pftpar(pft,39)/1000.0   !CO
         ef_trace(pft,3)=pftpar(pft,40)/1000.0   !CH4
         ef_trace(pft,4)=pftpar(pft,41)/1000.0   !VOC
         ef_trace(pft,5)=pftpar(pft,42)/1000.0   !TPM
         ef_trace(pft,6)=pftpar(pft,43)/1000.0   !NOx
       enddo
c     ASSIGN PFT PARAMETER FOR Glob-FIRM ONLY !!!
c     Assign PFT resistance to fire

      do pft=1,npft
        resist(pft)=pftpar(pft,8)
      enddo

c     START OF SPITFIRE, when year of observed fire reached
c     otherwise Glob-FIRM

c    !!!!! 8 Feb 2006 NEW: SPITFIRE in SPIN-UP !!!!
c     if (year.lt.0) then    
c     if (year.le.spinup_years) then        !!!! Yan: 26/10/07
      if (year.le.1.0) then        
           
c     Calculate the length of the fire season (units=days)
      i=1
      do d=1,365		!day of year
c       Calculate today's fire probability, fire_prob
c       Assume fire is only possible when temperature is above zero
        if (dtemp(d).gt.0.0.and.moistfactor.gt.0.0) then
          fire_prob=EXP((-pi)*(dw1(d)/moistfactor)**2)
        else
          fire_prob=0.0
        endif

        fire_length=fire_length+fire_prob
      enddo			!day of year
      an_fseason=fire_length

c     Calculate annual fire index

      fire_index=fire_length/365.0

c     Calculate the available fuel (above-ground litter) to carry the fire

      do pft=1,npft
        fuel=fuel+litter_ag_leaf(pft,1)+litter_ag_wood(pft,1)
      enddo
      fire_term=fire_index-1.0
      afire_frac=fire_index*EXP(fire_term/(0.45*fire_term**3
     *    + 2.83*fire_term**2 + 2.96*fire_term + 1.04))

      call pixelarea(lat,area)

        an_areafires=(afire_frac*area)/10000.0

c     Reduce fraction of grid cell affected by fire when fuel
c     becomes limiting (reduced carrying capacity)

        if (fuel.lt.minfuel) then
           an_areafires=0.0
           afire_frac=0.0
         endif

      if (afire_frac.lt.0.001) afire_frac=0.001

c     Implement the effect of the fire on vegetation structure and litter
c     in the disturbed fraction.

c     Each PFT is assigned a resistance to fire, representing the fraction of
c     the PFT which survives a fire. Grasses assumed already to have completed
c     their life cycle and thus are not affected by fire, giving them
c     a competitive advantage against woody PFTs.

c     Calculate trace gas emission resulting from biomass burning

      do pft=1,npft
        if (present(pft).and.tree(pft)) then
c         Calculate the fraction of individuals in grid cell which die
          disturb=(1.0-resist(pft))*afire_frac
c         Calculate carbon flux to atmosphere (gC/m2) due to burnt biomass
          acflux_fire(1)=acflux_fire(1)+disturb*(nind(pft)*
     *      (lm_ind(pft,1)+sm_ind(pft,1)+hm_ind(pft,1)+rm_ind(pft,1)))
          all_ind=lm_ind(pft,1)+sm_ind(pft,1)+hm_ind(pft,1)
     *         +rm_ind(pft,1)
          acflux_fire_pft(pft,1)=disturb*(nind(pft)*all_ind)
          if (all_ind .gt. 0.0) then
          do nc=2,nco2
             acflux_fire(nc)=acflux_fire(nc)+(lm_ind(pft,nc)
     *            * lm_ind(pft,1)
     *            +sm_ind(pft,nc)*sm_ind(pft,1)+hm_ind(pft,nc)*
     *            hm_ind(pft,1)+rm_ind(pft,nc)*rm_ind(pft,1))/all_ind
             acflux_fire_pft(pft,nc)=(lm_ind(pft,nc)*lm_ind(pft,1)
     *            +sm_ind(pft,nc)*sm_ind(pft,1)+hm_ind(pft,nc)*
     *            hm_ind(pft,1)+rm_ind(pft,nc)*rm_ind(pft,1))/all_ind
          enddo
          endif

c         Update the individual density
          nind(pft)=nind(pft)*(1.0-disturb)
        endif

c       Add combusted litter to carbon flux to atmosphere term

        temp=acflux_fire(1)
        acflux_fire(1)=acflux_fire(1)+(afire_frac*
     *    (litter_ag_leaf(pft,1)+litter_ag_leaf(pft,1)))
        if (acflux_fire(1).gt.0.) then
           do nc=2,nco2
              litter_inc=litter_ag_leaf(pft,nc)*litter_ag_leaf(pft,1)+
     *          litter_ag_wood(pft,nc)*litter_ag_wood(pft,1)
              acflux_fire(nc)=(acflux_fire(nc)*temp+
     *             afire_frac*litter_inc)/acflux_fire(1)
           enddo
        else
           do nc=2,nco2
              acflux_fire(nc)=0.0
           enddo
        endif


        temp=acflux_fire_pft(pft,1)
        acflux_fire_pft(pft,1)=acflux_fire_pft(pft,1)+(afire_frac
     *       *(litter_ag_leaf(pft,1)+litter_ag_wood(pft,1)))
        if (acflux_fire_pft(pft,1).gt.0.) then
           do nc=2,nco2
              litter_inc=litter_ag_leaf(pft,nc)*litter_ag_leaf(pft,1)+
     *          litter_ag_wood(pft,nc)*litter_ag_wood(pft,1)
              acflux_fire_pft(pft,nc)=(acflux_fire_pft(pft,nc)*temp+
     *             afire_frac*litter_inc)/acflux_fire_pft(pft,1)
           enddo
        else
           do nc=2,nco2
              acflux_fire_pft(pft,nc)=0.0
           enddo
        endif


c         Calculate trace gas emissions (gSpecies/m?²)

        do x=1,6
           if (acflux_fire(1).gt.0.0) then
               acflux_trace(x)=acflux_trace(x)+
     *                (acflux_fire_pft(pft,1)*ef_trace(pft,x)/0.45)
           else
               acflux_trace(x)=0.0
           endif

        enddo

c       Update the above ground litter term

C        litter_ag(pft,1)=(1.0-afire_frac)*litter_ag(pft,1)
        litter_ag_leaf(pft,1)=(1.0-afire_frac)*litter_ag_leaf(pft,1)
        litter_ag_wood(pft,1)=(1.0-afire_frac)*litter_ag_wood(pft,1)
		
c       fuel classes
		  
          fuel_1hr_leaf(pft,1)=(1.0-afire_frac)*fuel_1hr_leaf(pft,1)
          fuel_1hr_wood(pft,1)=(1.0-afire_frac)*fuel_1hr_wood(pft,1)

          fuel_10hr(pft,1)=(1.0-afire_frac)*fuel_10hr(pft,1)
          fuel_100hr(pft,1)=(1.0-afire_frac)*fuel_100hr(pft,1)
          fuel_1000hr(pft,1)=(1.0-afire_frac)*fuel_1000hr(pft,1)


      enddo

      else !SPITFIRE


ccccccccccccccc  START OF SPITFIRE ccccccccccccccccccccccc
c    ASSIGN PFT PARAMETER VALUES FOR SPITFIRE ONLY!!
     
c    fuel bulk density
*      dens_fuel_ave=0.0
      do pft=1,npft
        if (present(pft)) dens_fuel(pft)=pftpar(pft,37)
      enddo

c    for calculation of fire perimeter calculate tree and grass coverage
        do pft=1,npft
           if (present(pft)) then
            if (tree(pft)) then
              fpc_tree_total=fpc_tree_total+fpc_grid(pft)
            else
              fpc_grass_total=fpc_grass_total+fpc_grid(pft)

            endif
           endif
        enddo


c      f=0.0
      do pft=1,npft

c     crown length as a proportion of tree height
c         crown(pft)=pftpar(pft,44)

c     scorch height parameter for crown fire
c         f=f+(fpc_grid(pft)/fpc_tree_total)*pftpar(pft,45)
          f(pft)=pftpar(pft,45)
          g(pft)=pftpar(pft,46)
c     r(ck) and p parameter for postfire mortality as a result of crown damage
         r_ck(pft)=pftpar(pft,53)
         p(pft)=pftpar(pft,54)

      enddo

       call pixelarea(lat,area)
       area_ha=(area/10000.0)

c---------------------------------------
c  REG-FIRM at daily time step
c---------------------------------------


        m=1

        count=1
        count_int=1
        count_yr=1
        count_fdi=1
        afire_frac=0.0
		
        if(year.eq.1.or.year.eq.1000)afire_frac_afap_old=0.0


       fuel_1000hr_total_0=0.0
       fuel_100hr_total_0=0.0
       fuel_10hr_total_0=0.0
       fuel_1hr_total_0=0.0
       livegrass_0=0.0
       fuel_all_0=0.0
       dead_fuel_all_0=0.0
       dead_fuel_0=0.0
	   net_fuel=0.0
	   
        dead_fuel=0.0
        fuel_1hr_total=0.0
        fuel_1hr_leaf_total=0.0
        fuel_1hr_wood_total=0.0
        fuel_10hr_total=0.0
        fuel_100hr_total=0.0
        fuel_1000hr_total=0.0
        dens_fuel_ave=0.0
        livegrass=0.0


         dfuel_leaf=0.0    !Doug 03/09
        test_growing_start=.TRUE.    !Doug 06/09 fire paradox experiments


        fuel_load_subtraction(:,:)=0.0
      DO pft=1,npft
       IF(fuel_1hr_leaf(pft,1).lt.0.0)fuel_1hr_leaf(pft,1)=0.0
       IF(fuel_1hr_wood(pft,1).lt.0.0)fuel_1hr_wood(pft,1)=0.0
       IF(fuel_10hr(pft,1).lt.0.0)fuel_10hr(pft,1)=0.0
       IF(fuel_100hr(pft,1).lt.0.0)fuel_100hr(pft,1)=0.0
       IF(fuel_1000hr(pft,1).lt.0.0)fuel_1000hr(pft,1)=0.0
      END DO
	  
      DO pft=1,npft
        litter_ag_total=litter_ag_total+litter_ag_leaf(pft,1)
     *    +litter_ag_wood(pft,1)
		  
      END DO
      f1hr_leaf_frac(pft)=fuel_1hr_leaf(pft,1)/
     *  (fuel_1hr_leaf(pft,1)+fuel_1hr_wood(pft,1))
	  
        do pft=1,npft
        do co2=1,nco2
       fuel_1hr_leaf_0(pft,co2)=fuel_1hr_leaf(pft,co2)
       fuel_1hr_wood_0(pft,co2)=fuel_1hr_wood(pft,co2)
       fuel_10hr_0(pft,co2)=fuel_10hr(pft,co2)
       fuel_100hr_0(pft,co2)=fuel_100hr(pft,co2)  
       fuel_1000hr_0(pft,co2)=fuel_1000hr(pft,co2)
       fuel_1hr_leaf_left(pft,co2)=fuel_1hr_leaf(pft,co2)/0.45  	   
       fuel_1hr_wood_left(pft,co2)=fuel_1hr_wood(pft,co2)/0.45  
       fuel_10hr_left(pft,co2)=fuel_10hr(pft,co2)/0.45  
       fuel_100hr_left(pft,co2)=fuel_100hr(pft,co2)/0.45  
       fuel_1000hr_left(pft,co2)=fuel_1000hr(pft,co2)/0.45  
       lm_ind_0(pft,co2)=lm_ind(pft,co2)
       end do
       nind_0(pft)=nind(pft)
       end do

		  
        do pft=1,npft
          if (present(pft)) then

            fuel_1hr_total=fuel_1hr_total+(fuel_1hr_leaf(pft,1)+
     *        fuel_1hr_wood(pft,1))/0.45
	 
            fuel_1hr_leaf_total=fuel_1hr_leaf_total+
     *        fuel_1hr_leaf(pft,1)/0.45
            fuel_1hr_wood_total=fuel_1hr_wood_total+
     *        fuel_1hr_wood(pft,1)/0.45

            if (tree(pft)) then

              fuel_10hr_total=fuel_10hr_total+fuel_10hr(pft,1)/0.45
              fuel_100hr_total=fuel_100hr_total+fuel_100hr(pft,1)/0.45

              fuel_1000hr_total=fuel_1000hr_total+
     *          fuel_1000hr(pft,1)/0.45

            else !grass
c   KIRSTEN:  take proportion of grass leafmass, when green grass leaves are on
c         todays amount of green grass leaves: [gC/m2],influence on ROS only through moist_lg_1hr

           do d=1,365	!day
              livegrass=livegrass+
     *          (lm_ind(pft,1)/0.45*nind(pft))*dphen(d,pft) 
c     *                        (lm_ind(pft)/0.45*nind(pft))  
				
c              used in fire effects section only
c               pot_fc_lg(pft)=(lm_ind(pft)/0.45*nind(pft))   
            end do	!day
			
            endif !if(tree)

          endif ! if(present(pft))
        enddo !do pft

        
        dead_fuel = fuel_1hr_total + fuel_10hr_total
     *     + fuel_100hr_total  !total dead fuel g/m2


         livegrass_0=livegrass*0.45
         dead_fuel_all_0=(dead_fuel+fuel_1000hr_total)*0.45 
         dead_fuel_0=dead_fuel*0.45 ![gC/m2]
         fuel_1hr_total_0=fuel_1hr_total*0.45
         fuel_10hr_total_0=fuel_10hr_total*0.45
         fuel_100hr_total_0=fuel_100hr_total*0.45
         fuel_1000hr_total_0=fuel_1000hr_total*0.45
         fuel_all_0=dead_fuel_all_0+livegrass_0       ! [gC/m2]

        !IF (dead_fuel.le.0.0) THEN
        !  pfuel_limit=1;   !Doug 12/12: propotion of days with fuel limitation
        !                   != 1 as not enough fuel this year for any fire
        !  GOTO 200
        !END IF

c       net fuel load
        !if (dead_fuel.gt.0.0)
     *  !  net_fuel=(1.0-MINER_TOT)*(dead_fuel/1000.0)  ! in kg biomass
        !
        !if (net_fuel.le.0.0) THEN
        !  pfuel_limit=1;   !Doug 12/12: propotion of days with fuel limitation
        !                   != 1 as not enough fuel this year for any fire
        !  GOTO 200
        !END IF
		
c  TODO : remove if not inside if ...

C  Doug 01/09 MFL: change fuel varables to values before this years reproducton, turnover,
C  littersum (dcay), kill_pft, allocaton, light & mortality subroutines did their thing.



       fuel_1hr_leaf_inc_pos=fuel_1hr_leaf_inc_pos/.45
       fuel_1hr_leaf_inc_neg=fuel_1hr_leaf_inc_neg/.45
       fuel_1hr_wood_inc_pos=fuel_1hr_wood_inc_pos/.45
       fuel_1hr_wood_inc_neg=fuel_1hr_wood_inc_neg/.45
       fuel_10hr_inc=fuel_10hr_inc/.45
       fuel_100hr_inc=fuel_100hr_inc/.45
       fuel_1000hr_inc=fuel_1000hr_inc/.45



       dfuel_leaf=0.0

       DO pft=1,npft
         leaf_decayc(pft)=-1/pftpar(pft,10)
       END DO

       CALL fuel_1hr_redist(dfuel_leaf,
     *        fuel_1hr_leaf_inc_pos,fuel_1hr_leaf_inc_neg,
     *        fuel_1hr_wood_inc_pos,fuel_1hr_wood_inc_neg,
     *        leaf_decayc,dphen,present,tree,npft,
     *        nco2)

      DO pft=1,npft
        fuel_1hr_leaf_left(pft,1)=fuel_1hr_leaf_left(pft,1)
     *    -SUM(fuel_1hr_leaf_inc_pos(pft,1,:))
     *    -SUM(fuel_1hr_leaf_inc_neg(pft,1,:))-SUM(dfuel_leaf(:,pft))
	 
        fuel_1hr_wood_left(pft,1)=fuel_1hr_wood_left(pft,1)
     *    -SUM(fuel_1hr_wood_inc_pos(pft,1,:))
     *    -SUM(fuel_1hr_wood_inc_neg(pft,1,:))

        fuel_10hr_left(pft,1)=fuel_10hr_left(pft,1)-
     *    SUM(fuel_10hr_inc(pft,1,:))


        fuel_100hr_left(pft,1)=fuel_100hr_left(pft,1)-
     *    SUM(fuel_100hr_inc(pft,1,:))


        fuel_1000hr_left(pft,1)=fuel_1000hr_left(pft,1)-
     *    SUM(fuel_1000hr_inc(pft,1,:))

      END DO



      fuel_left_minus(:,:)=0.0        !Doug 01/09 MFL

      mfuel_1hr_total(:)=0.0
      mfuel_1hr_leaf_total(:)=0.0
      mfuel_1hr_wood_total(:)=0.0
      mfuel_10hr_total(:)=0.0
      mfuel_100hr_total(:)=0.0
      mfuel_1000hr_total(:)=0.0
      mlivegrass(:)=0.0               !Doug 05/09


      do d=1,365	!day

c        d_numf_old(d)=0.0
        d_i_surface(d)=0.0
        fc_crown(d)=0.0
        wind_speed=0.0
!        fire_frac=0.0
        fire_frac(d)=0.0
        fuel_consum=0.0
!        livegrass=0.0
c        deadgrass=0.0  !not used
        ratio_lg_reduc=0.0
!        ratio_dead_fuel=0.0
!        ratio_live_fuel=0.0

        do pft=1,npft
        do nc=1,nco2 
          fc_lg(pft,nc)=0.0
          fc_1hr_leaf(pft,nc)=0.0
          fc_1hr_wood(pft,nc)=0.0
          fc_10hr(pft,nc)=0.0
          fc_100hr(pft,nc)=0.0
          fc_1000hr(pft,nc)=0.0
          pot_fc_lg(pft,nc)=0.0
         pot_fc_lg_0(pft,nc)=0.0		  
        enddo
        nind_fa(pft)=0.0
        nind_kill(pft)=0.0
        nind_nkill(pft)=0.0
c        do class=1,5
          ck(pft)=0.0
c        enddo
        pm_ck(pft)=0.0
        pm_tau(pft)=0.0
        postf_mort(pft)=0.0
        sh(pft)=0.0
      enddo

       d_numfire(d)=0.0
       d_area_burnt(d)=0.0
     

C    Doug 01/09 MFL: update fuel loads based on the increment of change to the fuel load
C    by this years reproducton, turnover, littersum (dcay), kill_pft, allocaton, light &
C    mortality subroutines.

       dfuel_update_1hr_leaf(:,:)=0.0
       dfuel_update_1hr_wood(:,:)=0.0

       DO pft=1,npft
          IF (present(pft)) THEN

            temp=fuel_1hr_leaf_left(pft,1)
            fuel_1hr_leaf_left(pft,1)=fuel_1hr_leaf_left(pft,1)*
     *        (1-fire_frac(max(1,d-1)))+
     *        fuel_1hr_leaf_inc_pos(pft,1,m)/month_length(m)+
     *        fuel_1hr_leaf_inc_neg(pft,1,m)/month_length(m)+
     *        fuel_left_minus(pft,1)+dfuel_leaf(d,pft)+
     *        dfuel_update_1hr_leaf(pft,1)/0.45	!stopc
	 
            IF (fuel_1hr_leaf_left(pft,1)>0.0) THEN
              DO nc=2,nco2
                fuel_1hr_leaf_left(pft,nc)=(fuel_1hr_leaf_left(pft,nc)
     *            *temp+(fuel_1hr_leaf_inc_pos(pft,1,m)*
     *            fuel_1hr_leaf_inc_pos(pft,nc,m)+
     *            fuel_1hr_leaf_inc_neg(pft,1,m)*
     *            fuel_1hr_leaf_inc_neg(pft,nc,m))/
     *            month_length(m)+
     *            fuel_left_minus(pft,1)*
     *            fuel_1hr_leaf_left(pft,nc)+dfuel_leaf(d,pft)*
     *            fuel_1hr_leaf(pft,nc)+dfuel_update_1hr_leaf(pft,1)*
     *            dfuel_update_1hr_leaf(pft,nc)/0.45)/
     *            fuel_1hr_leaf_left(pft,1)
              END DO
            END IF	 
	 
            temp=fuel_1hr_wood_left(pft,1)

	        fuel_1hr_wood_left(pft,1)=fuel_1hr_wood_left(pft,1)*
     *        (1-fire_frac(max(1,d-1)))+
     *        fuel_1hr_wood_inc_pos(pft,1,m)/month_length(m)+
     *        fuel_1hr_wood_inc_neg(pft,1,m)/month_length(m)+
     *        fuel_left_minus(pft,2)+
     *        dfuel_update_1hr_wood(pft,1)/0.45	!stopc

            IF (fuel_1hr_wood_left(pft,1)>0.0) THEN
              DO nc=2,nco2
                fuel_1hr_wood_left(pft,nc)=(fuel_1hr_wood_left(pft,nc)
     *            *temp+(fuel_1hr_wood_inc_pos(pft,1,m)*
     *            fuel_1hr_wood_inc_pos(pft,nc,m)+
     *            fuel_1hr_wood_inc_neg(pft,1,m)*
     *            fuel_1hr_wood_inc_neg(pft,nc,m))/
     *            month_length(m)+
     *            fuel_left_minus(pft,2)*
     *            fuel_1hr_wood_left(pft,nc)+
     *            dfuel_update_1hr_wood(pft,1)*
     *            dfuel_update_1hr_wood(pft,nc)/0.45)/
     *            fuel_1hr_wood_left(pft,1)
              END DO
            END IF
	 
C            temp=fuel_1hr_left(pft,1)
C
C            fuel_1hr_left(pft,1)=fuel_1hr_left(pft,1)*
C     *        (1-fire_frac(min(1,d-1)))+
C     *        fuel_1hr_inc_pos(pft,1,m)/month_length(m)+
C     *        fuel_1hr_inc_neg(pft,1,m)/month_length(m)+
C     *        fuel_left_minus(pft,1)+dfuel_leaf(d,pft)+
C     *        dfuel_update_1hr(pft,1)/0.45	!stopc

C            IF (fuel_1hr(pft,1)>0.0) THEN
C              DO nc=2,nco2
C                fuel_1hr_left(pft,nc)=(fuel_1hr_left(pft,nc)
C     *            *temp+(fuel_1hr_inc_pos(pft,1,m)*
C     *            fuel_1hr_inc_pos(pft,nc,m)+
C     *            fuel_1hr_inc_neg(pft,1,m)*
C     *            fuel_1hr_inc_neg(pft,nc,m))/
C     *            month_length(m)+
C     *            fuel_left_minus(pft,1)*
C     *            fuel_1hr_left(pft,nc)+dfuel_leaf(d,pft)*
C     *            fuel_1hr(pft,nc)+dfuel_update_1hr(pft,1)*
C     *            dfuel_update_1hr(pft,nc)/0.45)/
C     *            fuel_1hr_left(pft,1)
C              END DO
C            END IF


            temp=fuel_10hr_left(pft,1)
            fuel_10hr_left(pft,1)=fuel_10hr_left(pft,1)*
     *        (1-fire_frac(max(1,d-1)))+
     *        (fuel_10hr_inc(pft,1,m))/month_length(m)+
     *        fuel_left_minus(pft,3)+
     *        dfuel_update_10hr(pft,1)/0.45

            IF (fuel_10hr(pft,1)>0.0) THEN
              DO nc=2,nco2
                fuel_10hr_left(pft,nc)=(fuel_10hr_left(pft,nc)*
     *            *temp+fuel_10hr_inc(pft,1,m)*
     *            fuel_10hr_inc(pft,nc,m)/month_length(m)+
     *            fuel_left_minus(pft,3)*fuel_10hr_left(pft,nc)
     *            +dfuel_update_10hr(pft,1)*
     *            dfuel_update_10hr(pft,nc)/0.45)/
     *            fuel_10hr_left(pft,1)
               END DO
             END IF


            temp=fuel_100hr_left(pft,1)
            fuel_100hr_left(pft,1)=fuel_100hr_left(pft,1)*
     *        (1-fire_frac(max(1,d-1)))+
     *        (fuel_100hr_inc(pft,1,m))/month_length(m)+
     *        fuel_left_minus(pft,4)+
     *        dfuel_update_100hr(pft,1)/0.45

            IF (fuel_100hr(pft,1)>0.0) THEN
              DO nc=2,nco2
                fuel_100hr_left(pft,nc)=(fuel_100hr_left(pft,nc)*
     *            *temp+fuel_100hr_inc(pft,1,m)*
     *            fuel_100hr_inc(pft,nc,m)/month_length(m)+
     *            fuel_left_minus(pft,4)*fuel_100hr_left(pft,nc)
     *            +dfuel_update_100hr(pft,1)*
     *            dfuel_update_100hr(pft,nc)/0.45)/
     *            fuel_100hr_left(pft,1)
               END DO
             END IF


            temp=fuel_1000hr_left(pft,1)
            fuel_1000hr_left(pft,1)=fuel_1000hr_left(pft,1)*
     *        (1-fire_frac(max(1,d-1)))+
     *        (fuel_1000hr_inc(pft,1,m))/month_length(m)+
     *        fuel_left_minus(pft,5)+
     *        dfuel_update_1000hr(pft,1)/0.45

            IF (fuel_10hr(pft,1)>0.0) THEN
              DO nc=2,nco2
                fuel_1000hr_left(pft,nc)=(fuel_1000hr_left(pft,nc)*
     *            *temp+fuel_1000hr_inc(pft,1,m)*
     *            fuel_1000hr_inc(pft,nc,m)/month_length(m)+
     *            fuel_left_minus(pft,5)*fuel_1000hr_left(pft,nc)
     *            +dfuel_update_1000hr(pft,1)*
     *            dfuel_update_1000hr(pft,nc)/0.45)/
     *            fuel_1000hr_left(pft,1)
               END DO
             END IF

          END IF
        END DO






C    Doug 01/09 MFL: Set any minus fuel load (which is just silly) to zero and record the
C    difference for the next loop. this is not required for fuel_xhr_left as its not used
C    to calculated anything asn is reset at the end of the loop


       DO pft=1,npft
         IF (fuel_1hr_leaf_left(pft,1)<0.0) THEN
           fuel_left_minus(pft,1)=fuel_1hr_leaf_left(pft,1)
           fuel_1hr_leaf_left(pft,1)=0.0
         ELSE
           fuel_left_minus(pft,1)=0
         END IF
		 
         IF (fuel_1hr_wood_left(pft,1)<0.0) THEN
           fuel_left_minus(pft,2)=fuel_1hr_wood_left(pft,1)
           fuel_1hr_wood_left(pft,1)=0.0
         ELSE
           fuel_left_minus(pft,2)=0
         END IF
		 
C         IF (fuel_1hr_left(pft,1)<0.0) THEN
C           fuel_left_minus(pft,1)=fuel_1hr_left(pft,1)
C           fuel_1hr_left(pft,1)=0.0
C         ELSE
C           fuel_left_minus(pft,1)=0
C         END IF
  
         IF (fuel_10hr_left(pft,1)<0.0) THEN
           fuel_left_minus(pft,3)=fuel_10hr_left(pft,1)
           fuel_10hr_left(pft,1)=0.0
         ELSE
           fuel_left_minus(pft,3)=0
         END IF

         IF (fuel_100hr_left(pft,1)<0.0) THEN
           fuel_left_minus(pft,4)=fuel_100hr_left(pft,1)
           fuel_100hr_left(pft,1)=0.0
         ELSE
           fuel_left_minus(pft,4)=0
         END IF

         IF (fuel_1000hr_left(pft,1)<0.0) THEN
           fuel_left_minus(pft,5)=fuel_1000hr_left(pft,1)
           fuel_1000hr_left(pft,1)=0.0
         ELSE
           fuel_left_minus(pft,5)=0
         END IF
       END DO

       fuel_1hr_total=SUM(fuel_1hr_leaf_left(:,1))+
     *   SUM(fuel_1hr_wood_left(:,1))
	 
       fuel_1hr_leaf_total=SUM(fuel_1hr_leaf_left(:,1))
       fuel_1hr_wood_total=SUM(fuel_1hr_wood_left(:,1))

       fuel_10hr_total=0.0
       fuel_100hr_total=0.0
       fuel_1000hr_total=0.0


       DO pft=1,npft
         IF (tree(pft)) THEN
           fuel_10hr_total=fuel_10hr_total+fuel_10hr_left(pft,1)           
           fuel_100hr_total=fuel_100hr_total+fuel_100hr_left(pft,1)
           fuel_1000hr_total=fuel_1000hr_total+fuel_1000hr_left(pft,1)
         END IF
       END DO

        fuel_1hr_leaf_left(:,1)=fuel_1hr_leaf_left(:,1)*0.45
        fuel_1hr_wood_left(:,1)=fuel_1hr_wood_left(:,1)*0.45
        fuel_10hr_left(:,1)=fuel_10hr_left(:,1)*0.45
        fuel_100hr_left(:,1)=fuel_100hr_left(:,1)*0.45
        fuel_1000hr_left(:,1)=fuel_1000hr_left(:,1)*0.45

 

c    FUEL CHARACTERISTICS

c    net fuel load and total amount of dead fuel per fuel class
c    and amount of livegrass
c    ACHTUNG: take grass phenology into account (cf. grass curing) per time step
c             to calculate amount of live grass!!


        dens_fuel_ave=0.0
        livegrass=0.0
		net_fuel=0.0

       fuel_1000hr_total=0.0
        
        do pft=1,npft
          if (present(pft)) then


            if (tree(pft)) then
                fuel_1000hr_total=fuel_1000hr_total	!Doug 03/09: fuel_1000hr_0(pft,1) --> fuel_1000hr_left(pft,1)	(annual to daily)
     *              +fuel_1000hr_left(pft,1)*(1.0-afire_frac)/0.45
            else !grass
c   KIRSTEN:  take proportion of grass leafmass, when green grass leaves are on
c         todays amount of green grass leaves: [gC/m2],influence on ROS only through moist_lg_1hr

              livegrass=livegrass+
     *(lm_ind_0(pft,1)/0.45*nind_0(pft))*(1.0-afire_frac)*dphen(d,pft)  !g/m2

c     *                        (lm_ind(pft)/0.45*nind(pft))  !g/m2


c              used in fire effects section only
               pot_fc_lg(pft,1)=(lm_ind_0(pft,1)/0.45*nind_0(pft))
     *             *(1.0-afire_frac)*dphen(d,pft) !g/m2 
              pot_fc_lg(pft,2)=lm_ind(pft,2)
              pot_fc_lg(pft,3)=lm_ind(pft,3)
c               pot_fc_lg(pft)=(lm_ind(pft)/0.45*nind(pft))  !g/m2 

		pot_fc_lg_0(pft,1)=(lm_ind_0(pft,1)/0.45*nind_0(pft))
     *		*dphen(d,pft) !g/m2 

            endif

          endif
        enddo	!pft


c Doug 05/09: Fire Paradox experiments
C Experiment 1&3: remove 1 hr fuel for 
c 
c see if the day is the 1st day of the growning season
c        IF (year>=5140.AND.year<5190) THEN
c        IF (year>=5160) THEN
c          DO pft=1,npft
c            IF (tree(pft)==.FALSE.) THEN        
c              IF (dphen_change(m,pft)==d.AND.
c     *          test_growing_start) THEN
c                fuel_load_subtraction(:,1)=fuel_1hr_left(:,1)
c                fuel_load_subtraction(:,2)=fuel_10hr_left(:,1)
c                fuel_load_subtraction(:,3)=fuel_100hr_left(:,1)
c                fuel_load_subtraction(:,4)=fuel_1000hr_left(:,1)
c                fuel_1hr_left(:,1)=0 !experiment 1&3
c                fuel_1hr_total=0     !experiment 1&3
c 
c                fuel_10hr_left(:,1)=0 !experiment 3
c                fuel_10hr_total=0     !experiment 3C
c 
c                fuel_100hr_left(:,1)=0 !experiment 3
c                fuel_100hr_total=0     !experiment 3
c 
c 
c                test_growing_start=.FALSE.
c              END IF
c            END IF
c          END DO
c        END IF
c 
c Experiment 2&4
c        IF (year>=5180.AND.year<5190) THEN
c          IF (d==ndaymonth(m).AND.(dphen(d,8)==1.OR.dphen(d,9)==1)) THEN
c            fuel_load_subtraction(:,1)=fuel_load_subtraction(:,1)
c     *        +fuel_1hr_left(:,1) 
c            fuel_1hr_left(:,1)=0       !experiment 2
c            fuel_1hr_total=0           !experiment 2
c          END IF
c        END IF
c 
c          fuel_1000hr_left(:,1)=0    !experiment 4
c          fuel_1000hr_total=0        !experiment 4
c        END IF

        dead_fuel = fuel_1hr_total + fuel_10hr_total
     *     + fuel_100hr_total  !total dead fuel g/m2

c       net fuel load
        net_fuel=0.0
        if (dead_fuel.gt.0.0)
     *    net_fuel=(1.0-MINER_TOT)*(dead_fuel/1000.0)

c    fuel bulk density, weighted per fuel class and fuel load
c    ACHTUNG: WEIGHTING per fpc AND fuel load or per fuel load only? Reg-FIRM1: per FPC       
        do pft=1,npft
          if (present(pft)) then
            if (dead_fuel.gt.0.0) then
              ratio_fbd(pft) = ((fuel_1hr_leaf_left(pft,1) +
     *              fuel_1hr_wood_left(pft,1) +
     *              fbd_a * fuel_10hr_left(pft,1) + 	!Doug 03/09: fuel_xhr--> fuel_xhr_left (annual to daily fuel)
     *                     fbd_b * fuel_100hr_left(pft,1)) / 0.45) 
     *              / dead_fuel
              dens_fuel_ave = dens_fuel_ave+dens_fuel(pft)*
     *              ratio_fbd(pft)    
            else
              dens_fuel_ave=0.0
            endif
          endif
        enddo


c    livegrass

       if (livegrass.gt.0.0) then
c    KIRSTEN: calculate leaf moisture content of grasses from dw1
c    ACHTUNG: change with water scalar value for grasses!!!!
       dlm_lg(d) = max(0.0,((10.0/9.0)*dw1(d)-(1.0/9.0)))
       !dlm_lg(d)=dlm_lg(d)/(1.0-dlm_lg(d))
       if(lon.gt.25.0.and.lon.lt.26.0)then
       !write(1005,*)'dlm_lg,d,lat,lon',dlm_lg,d,lat,lon	!Doug 12/08: Commented out to save space
       !write(1006,*)'dw1(d),d,lat,lon',dw1(d),d,lat,lon		!Doug 12/08: Commented out to save space
       end if
       !Doug 03/13: need to un-hardcode grass pfft number
       ratio_C3_livegrass = pot_fc_lg(8,1) / livegrass
       ratio_C4_livegrass = pot_fc_lg(9,1) / livegrass
       else
         dlm_lg(d)=0.0
       ratio_C3_livegrass = 0.0 
       ratio_C4_livegrass = 0.0
       endif

       if(livegrass.gt.0.0)then
c      influence of livegrass on FBD
c    Kirsten: presence of livegrass additionally lowers FBD, when trees mixed 
c    with grasses
       dens_livegrass_ave =  
     *       fbd_C3_livegrass *  ratio_C3_livegrass +
     *      fbd_C4_livegrass *  ratio_C4_livegrass
         else
           dens_livegrass_ave=0.0
         end if


        if (dead_fuel.gt.0.0 .or. livegrass.gt.0.0) then !K.
         ratio_dead_fuel = dead_fuel  / (dead_fuel + livegrass)
         ratio_live_fuel = livegrass / (dead_fuel + livegrass)
         char_dens_fuel_ave = dens_fuel_ave* ratio_dead_fuel
     *    + dens_livegrass_ave * ratio_live_fuel
        endif
         char_net_fuel = net_fuel +(1.0-MINER_TOT)*livegrass/1000.0 ! in kg biomass
		 if(char_net_fuel.le.0.0)goto 201
		 
         char_net_fuel_0=char_net_fuel

            
           IF(char_dens_fuel_ave<=0.0)then	!Doug 03/09: no longer stops, but reloops. Error is now noted in file fort.10
             !WRITE(10,*),'********************************'
             !WRITE(10,*),'error: char_dens_fuel_ave.le.0.0'
             !WRITE(10,*),'--------------------------------'
             !WRITE(10,*),'lat:',lat,'lon:',lon
             !WRITE(10,*), 'year:', year, 'day:', day
             !WRITE(10,*), 'dens_fuel_ave:', dens_fuel_ave
             !  WRITE(10,*), '    dens_fuel:', dens_fuel
             !  WRITE(10,*), '    ratio_fbd:', ratio_fbd
             !WRITE(10,*),'ratio_dead_fuel:', ratio_dead_fuel
             !WRITE(10,*),'+++++++++++++++++++++++++++++++'

c             stop
             GOTO 201
           end if

           fbdepth=char_net_fuel/char_dens_fuel_ave

        pack_ratio_dead=(pack_ratio_1hr*fuel_1hr_total+
     *                   pack_ratio_10hr*fuel_10hr_total+
     *                   pack_ratio_100hr*fuel_100hr_total)/dead_fuel

       pack_ratio=pack_ratio_dead*ratio_dead_fuel+pack_ratio_lg*
     *            ratio_live_fuel    !pack_ratio if the packing ratio of the fuel bed.

      char_dens_fuelbed=pack_ratio*part_dens  !this is the fuel bulk density.

          if(char_dens_fuelbed.ne.0.0)then
            fbdep=char_net_fuel/char_dens_fuelbed
          else
            fbdep=0.0
          end if


c      moistfactor
c     Kirsten: original method, where moistfactor is a PFT parameter 
c           moistfactor = (moistfactor_1hr * fuel_1hr_total
c     *          + moistfactor_10hr * fuel_10hr_total
c     *          + moistfactor_100hr * fuel_100hr_total)/dead_fuel

           moistfactor = (moistfactor_1hr * fuel_1hr_total
     *          + moistfactor_10hr * fuel_10hr_total
     *          + moistfactor_100hr * fuel_100hr_total)/dead_fuel

c Kirsten: moistfactor for livegrass less than for dead fuel!
       call fire_danger_index1(dlm,dlm_lg,dtemp_min,dtemp_max,dprec,	
     *         d,m,fuel_1hr_total,fuel_10hr_total,fuel_100hr_total,
     *         dead_fuel,ratio_dead_fuel,ratio_live_fuel,dlm_1hr,
     *         dlm_10hr,dlm_100hr,dlm_1000hr,year,ni_acc,
     *         char_alpha_fuel,d_ni,lon,lat,
     *         dlm_1hr_old,dlm_10hr_old,dlm_100hr_old,dlm_1000hr_old,
     *         mlm,dpet,aprec)
       if(livegrass.gt.0.0)then 
        moistfactor_livegrass = 2.9*(1.0/ratio_live_fuel-1.0)*
     *                     (1.0-dlm_1hr(d)/moistfactor_1hr)-0.226

       if(moistfactor_livegrass.lt.moistfactor_1hr)
     *    moistfactor_livegrass=moistfactor_1hr
        else
           moistfactor_livegrass=0.0
        end if

        char_moistfactor = moistfactor *  ratio_dead_fuel
     *       + moistfactor_livegrass * ratio_live_fuel

c    surface-area-to-volume ratio weighted among fuel classes
c     In future, account for SAV diffs between pine needles and
c     broad-leaves.  Previous studies report order
c     of magnitude difference in SAV b/w the two leaf types.

       if (dead_fuel.gt.0.0) then
         sigma=(fuel_1hr_total * sigma_1hr +
     *         fuel_10hr_total * sigma_10hr +
     *         fuel_100hr_total * sigma_100hr) / dead_fuel

c    higher sigma value for livegrass increases fire 
c    WHY should the presence of livegrass increase fire additionally?
c       char_sigma =sigma*ratio_dead_fuel+sigma_livegrass*ratio_live_fuel


      char_sigma =sigma*ratio_dead_fuel+sigma_livegrass*ratio_live_fuel
       
       else
         sigma=0.00001
c          char_sigma=0.00001

          char_sigma=0.00001
       endif
       
       
c Doug 06/09: mask areas of fuel limitation (i.e, where 1hr, 10hr and livegrass
c     is to small for propogation of fire

        IF (fuel_1hr_total+livegrass+fuel_10hr_total<200) THEN
           pfuel_limit=pfuel_limit+1/365 !Doug 12/12: add a day with fuel limitation
           GOTO 201
        END IF

        IF(dead_fuel.le.0.0) THEN
          pfuel_limit=pfuel_limit+1/365 !Doug 12/12: add a day with fuel limitation
          GOTO 201
        END IF
        
        IF (net_fuel.le.0.0) THEN
          pfuel_limit=pfuel_limit+1/365 !Doug 12/12: add a day with fuel limitation
          GOTO 201
        END IF

c influence of wind_speed on rate of forward spread
c wind speed from NCEP reanalysis data is in m/s, ROS is in m/min

        wind_speed=dwindsp(d)*60  !m/min
c  introduce reduction factor for wind speed with respect to forest or grass fpc 
c       ref: Rothermel 1983, Pyne 1996


        wind_speed=(fpc_tree_total*dwindsp(d)*60.0*0.4)+
     *               (fpc_grass_total*dwindsp(d)*60.0*0.6)
c        wind_speed=(fpc_tree_total*1.28*60.0*0.4)+
c     *               (fpc_grass_total*1.28*60.0*0.6)
c      converts wind_speed (m/min) to ft/min
c       for input into Rothermel's formula for phi_wind in the ROS S/R

        wind_forward=3.281*wind_speed       

c        back_ws=wind_speed*exp(-0.05039*wind_speed)     ! backward fire spread
c        wind_backward=3.281*back_ws

        
    
c       weight lb according to presence of grasses and forests
c       wind speed in Canadian formula is in km/h

c      CHECK UNITS!! Check formulae! Max LB ratio = 8
c       check fpc- proportion?


              lb=fpc_tree_total*
     *   (1.0+(8.729*((1.0-(exp(-0.03*0.06*wind_speed)))**2.155)))
     *    +(fpc_grass_total*(1.1+((0.06*wind_speed)**0.0464)))

        if(lb.lt.1.0)lb=1


c Kirsten: check alternatively:
c         else !CHECK THIS
c              lb=min(8,fpc_tree_total*
c     *   (1.0+(8.729*((1.0-(exp(-0.03*0.06*wind_speed)))**2.155)))
c     *    +(fpc_grass_total*(1.1+((0.06*wind_speed)**0.0464))))
c       endif
  

c-----------------------------------------------------------
c      end preparing daily variables - start of simulation
c-----------------------------------------------------------

c       climatic fire danger
c      Assume moisture content of livegrass close to 1.0 ? NO!! Kirsten --> driven by dw1
      

!      call fire_danger_index(d_fdi,dlm,dlm_lg,dtemp_min,dtemp_max,dprec,
!     *    d,moistfactor,fuel_1hr_total,fuel_10hr_total,
!     *    fuel_100hr_total,dead_fuel,
!     *    char_moistfactor, ratio_dead_fuel,ratio_live_fuel,
!     *    dlm_1hr,dlm_10hr, dlm_100hr, dlm_1000hr,year,ni_acc)

        call fire_danger_index(char_moistfactor,char_alpha_fuel,
     *                         ni_acc,d_NI,d_fdi,d)
         
           if (d_fdi(d).gt.0.0) then
           mfdi(m)= mfdi(m)+d_fdi(d)
           an_fdi=an_fdi+d_fdi(d)
           count=count+1
         endif
       
c  number of fires. only ignitions, when dead fuel load - save computation time

         if (net_fuel.gt.0.001) then 

      call fuel_property(dens_fuel_ave,sigma,dlm,char_dens_fuel_ave,
     *                   char_sigma,char_dens_fuelbed,d,bet,q_ig)


       if((0.2388*q_ig).gt.400.0)then
          p_ig=0.0
       else
          p_ig=0.000048*(((400.0-0.2388*q_ig)
     *         /10)**4.3)/50
       end if

          if(p_ig.gt.1.0)p_ig=1.0
 
        if(bet.lt.1.0)p_ig=p_ig*bet
	   dlightn_cons=sum(dlightn)/365	!Doug 12/08, trying out constant lightnng
           dlightn_sca(d)=dlightn(d)*0.15	! 0.15: continous current parameter; 0.2 Cloud-Ground: Total lighting ratio
c	   dlightn_sca(d)=dlightn_cons*0.15*0.2	!Doug 12/08 ""
            dlightn_sca(d)=dlightn_sca(d)*p_ig!*1	!Doug 12/08: play with scaling of ignitions from lightning

           d_numfire(d)=d_fdi(d)*
c     *       (human_ign(popden,a_nd,year)+dlightn_sca(d))*0.000001*
     *       dlightn_sca(d)*0.000001*
     *       area_ha*(1.0-afire_frac)
C Yan
           d_numfire_human(d)=d_fdi(d)*
     *       human_ign(popden,a_nd,year)*0.000001*area_ha
     *       *(1.0-afire_frac)

           d_numfire_lightn(d)=d_fdi(d)*

     *       dlightn_sca(d)*0.000001*area_ha*(1.0-afire_frac)		

             dhuman_ign(d)=human_ign(popden,a_nd,year)


c            continued burning of yesterdays ignitions, not accounting in the # fires statistics - consider later
*             if (d_i_surface(d-1).gt.1000.0.and.dprec(d).le.3.0)
*     *          d_numf_old(d)=d_numfire(d-1)
         else  !not enough fuel
           d_numfire(d)=0.0
C Yan
           d_numfire_human(d)=0.0
           d_numfire_lightn(d)=0.0

         endif

c    area burnt



c rate of spread
!/       if (net_fuel.gt.0.0) then

c     calculation of backward ROS taken out, because of inconsistencies
c     Allan: calculated after Hirsch (1996)
c         call rate_of_spread(U_front,wind_backward,dens_fuel_ave,sigma,
c     *       dlm,d,net_fuel,moistfactor,H, char_dens_fuel_ave, 
c     *             char_sigma, char_net_fuel, char_moistfactor,gamma)
        
c          ros_b(d)=U_front

       call rate_of_spread(U_front,wind_forward,dens_fuel_ave,sigma,dlm,
     *            d, net_fuel,moistfactor,H, char_dens_fuel_ave, 
     *   char_sigma, char_net_fuel, char_moistfactor,gamma,
     *    char_dens_fuelbed)

       ros_f(d)=U_front


 

c      CHECK FORMULAE: HIRSCH (1996)
c     this doubles ros_b(d) compare to the old method to estimate backward wind speed


       ros_b(d) = ros_f(d) * exp(-0.012 * wind_speed) !Can FBP System
           if (ros_b(d).lt. 0.05) ros_b(d) = 0.0 
     

c     check the parameter value!!
c     fire duration as a function of d_fdi


          fire_durat=361.0/(1.0+(((361.0/1.)-1.)*exp(-13.06*d_fdi(d))))
c          fire_durat=241.0/(1.0+(((241.0/1.)-1.)*exp(-11.06*d_fdi(d))))
          db(d)=ros_b(d)*fire_durat ! in min
          df(d)=ros_f(d)*fire_durat

c  area burnt from todays ignitions

         d_area_burnt(d)=d_numfire(d)*
     *        ((pi/(4.0*lb))*((df(d)+db(d))**2.0))/10000.0 ! Can FBP System
        fire_frac(d)=(d_area_burnt(d)*10000.0)/area

C Yan
        d_area_burnt_human(d)=d_numfire_human(d)*
     *        ((pi/(4.0*lb))*((df(d)+db(d))**2.0))/10000.0 ! Can FBP System
        d_area_burnt_lightn(d)=d_numfire_lightn(d)*
     *        ((pi/(4.0*lb))*((df(d)+db(d))**2.0))/10000.0 ! Can FBP System

        fire_frac_human=(d_area_burnt_human(d)*10000.0)/area
        fire_frac_lightn=(d_area_burnt_lightn(d)*10000.0)/area

c     Do not update annual burnt area stats yet. Need to assess FI threshold first.
c     Put info into temp vars for now. KIRSTEN: May not be necessary
c         afire_frac=afire_frac+fire_frac
          afire_frac_temp=afire_frac+fire_frac(d)
		  an_areafires_temp=an_areafires+d_area_burnt(d)

C Yan
          afire_frac_human_temp  
     *   =afire_frac_human +fire_frac_human
          afire_frac_lightn_temp  
     *   =afire_frac_lightn+fire_frac_lightn
          an_areafires_human_temp  
     *   =an_areafires_human+d_area_burnt_human(d)
          an_areafires_lightn_temp 
     *   =an_areafires_lightn+d_area_burnt_lightn(d)

c--------------------------------------------------------
c    prevent to burn the entire cell AND ensure fire_frac is consistent with its later use!
         if ((an_areafires_temp.gt.area_ha).or.
     *      (d_area_burnt(d).gt.area_ha)) then
               d_area_burnt(d)=area_ha-an_areafires_temp+d_area_burnt(d)
               fire_frac(d)=1.0-afire_frac_temp+fire_frac(d)
c               an_areafires_temp=area_ha
c               afire_frac_temp=1.0
         endif



c dead fuel consumption: in g biomass per m?² - not gC/m?² !!

      call fuel_consumption(npft,present,fuel_consum,fire_frac,
!     *  fuel_1hr,fuel_10hr,fuel_100hr,fuel_1000hr,livegrass,pot_fc_lg,
     *  fuel_1hr_leaf_left,fuel_1hr_wood_left,
     *  fuel_10hr_left,fuel_100hr_left,	!Doug 03/09: fuel_xhr_0-->fuel_xhr_left (annual to daily)
     *  fuel_1000hr_left,livegrass,pot_fc_lg_0,
     *  tree,fuel_1hr_total,moistfactor,dlm,d,MINER_TOT,
     *  fc_1hr_leaf,fc_1hr_wood,fc_lg,
     *  fc_10hr,fc_100hr,fc_1000hr,cf,char_moistfactor,
     *  dlm_lg, dlm_1hr, dlm_10hr, dlm_100hr, dlm_1000hr,
     *  moistfactor_livegrass, moistfactor_1hr,moistfactor_10hr,
     *  moistfactor_100hr, moistfactor_1000hr)


c    surface fire intensity
        if ((fire_frac(d).gt.0.0).and.(fuel_consum.gt.0.0)) then

       d_i_surface(d)=H*(fuel_consum/1000.0/fire_frac(d))*
     *     (ros_f(d)/60.0)

c       scorch height - when average over present PFT, better per PFT (see below)
c          sh=f*(d_i_surface(d)**0.667)

c    update of no. fires, area burnt, fire_frac, sh
c    ACHTUNG: there might be two influences of surface fire intensity !!!
c             a) as a stopping rule and b) as criteria for multiple-day burning
c             check references for numbers!
c    Fires may stop but continue to smoulder; or stop and extinguish.
c    Pyne et al (1986): 50 and 300 respectively. CHECK!

        endif !CHECK THIS


c    Fire effects for trees, when surface fire intensity is greater than 50 W/m
c    At the moment it has the same threshold as for extinguishing fires
C    Leilei 12/08 & Doug 12/08: changing intensity threshold for calibraton
           if (d_i_surface(d).le.0.0) then

           d_numfire(d)=0.0
           d_area_burnt(d)=0.0
           fire_frac(d)=0.0
           d_i_surface(d)=0.0
c Yan
           d_numfire_human(d)=0.0
           d_area_burnt_human(d)=0.0
           d_area_burnt_lightn(d)=0.0
           fire_frac_human=0.0

           d_i_surface_human(d)=0.0
           d_numfire_lightn(d)=0.0
           d_i_surface_lightn(d)=0.0
           fire_frac_lightn=0.0

c          sh=0.0
         else
           count_int=count_int+1
         endif

          num_fire(m)=num_fire(m)+d_numfire(d)
          annum_fire=annum_fire+d_numfire(d)

          an_areafires=an_areafires+d_area_burnt(d)

          if(an_areafires.le.0.0000001)an_areafires=0.0		  
          if(an_areafires.lt.0.0)then
          write(*,*)'an_areafires < 0.0'
          stop
          end if


          if(lat.ge.0.0)then !in The North, fire season starts from April 1.
          if(d.gt.90)then
            afire_frac=afire_frac+fire_frac(d)
          else
            afire_frac=afire_frac+afire_frac_afap_old+fire_frac(d)
          end if 
         if(d.eq.90)afire_frac=0.0
          else   !in the South, fire season starts from Sep 1.
           if(d.le.274)then !October 1
              afire_frac=afire_frac+afire_frac_afap_old+fire_frac(d)
            else
              afire_frac=afire_frac+fire_frac(d)
            end if
              if(d.eq.275)afire_frac=0.0
            end if


          m_i_surface(m)=m_i_surface(m)+d_i_surface(d)
          an_i_surface=an_i_surface+d_i_surface(d)
			
          area_burnt(m)=area_burnt(m)+d_area_burnt(d)
C Yan
          num_fire_human(m)=num_fire_human(m)+d_numfire_human(d)
          annum_fire_human=annum_fire_human+d_numfire_human(d)

          an_areafires_human=an_areafires_human+d_area_burnt_human(d)
          afire_frac_human=afire_frac_human+fire_frac_human

          m_i_surface_human(m)=m_i_surface_human(m)+d_i_surface_human(d)
          an_i_surface_human=an_i_surface_human+d_i_surface_human(d)

          area_burnt_human(m)=area_burnt_human(m)+d_area_burnt_human(d)

          num_fire_lightn(m)=num_fire_lightn(m)+d_numfire_lightn(d)
          annum_fire_lightn=annum_fire_lightn+d_numfire_lightn(d)

          an_areafires_lightn=an_areafires_lightn+d_area_burnt_lightn(d)
          afire_frac_lightn=afire_frac_lightn+fire_frac_lightn

          m_i_surface_lightn(m) 
     *  = m_i_surface_lightn(m)+d_i_surface_lightn(d) 
          an_i_surface_lightn
     *  = an_i_surface_lightn+d_i_surface_lightn(d)

          area_burnt_lightn(m)
     *  = area_burnt_lightn(m)+d_area_burnt_lightn(d)


c   FIRE EFFECTS - THIS IS NEW cf REG-FIRM V1: CROWN FIRES
      do pft=1,npft
        
       if (present(pft)) then

c    Fire effects for trees, when surface fire intensity is greater than 50 W/m
c    At the moment it has the same threshold as for extinguishing fires
C    Leilei 12/08 & Doug 12/08: changing intensity threshold for calibraton

       if (d_i_surface(d).gt.00.0) then


c   surface fire: update fuel load per dead fuel class and litter_ag in general

         fuel_1000hr(pft,1)=fuel_1000hr(pft,1)-fc_1000hr(pft,1)
         fuel_1000hr_left(pft,1)=fuel_1000hr_left(pft,1)-
     *     fc_1000hr(pft,1)	!Doug 03/09
	 
	 
         litter_ag_leaf(pft,1)=litter_ag_leaf(pft,1)-
     *     (fc_1hr_leaf(pft,1)*0.45)
	 
         litter_ag_wood(pft,1)=litter_ag_wood(pft,1)-
     *    (fc_1hr_wood(pft,1)*0.45)-
     *    (fc_10hr(pft,1)*0.45)-(fc_100hr(pft,1)*0.45)-fc_1000hr(pft,1)


c         litter_ag(pft,1)=litter_ag(pft,1)-(fc_1hr(pft,1)*0.45)-
c     *   (fc_10hr(pft,1)*0.45)-(fc_100hr(pft,1)*0.45)-fc_1000hr(pft,1)


	 

c        fuel_1hr_left(pft,1)=fuel_1hr_0(pft,1)*(1.0-afire_frac)
c        fuel_10hr_left(pft,1)=fuel_10hr_0(pft,1)*(1.0-afire_frac)
c        fuel_100hr_left(pft,1)=fuel_100hr_0(pft,1)*(1.0-afire_frac)
c        fuel_1000hr_left(pft,1)=fuel_1000hr_0(pft,1)*(1.0-afire_frac)
        livegrass_left(pft,1)=(lm_ind_0(pft,1)/0.45*nind_0(pft))
     *             *dphen(d,pft)*(1.0-afire_frac) ! g/m2

c   start of crown fires
       if (tree(pft)) then

c   scorch height per PFT
          sh(pft)=f(pft)*(d_i_surface(d)**0.667)

c      post-fire mortality from cambial damage
c      tau_r=cf/Gamma, tau_l=tau_r, if not tau_l=2*tau_r
          tau_l(pft)=2.0*(cf(pft)/gamma)

c   uniform distribution of height (5 classes), not needed with real age classes
       do class=0,4

c   crown kill in [%] assuming the crown shape being a cylinder
c   crown height as a fraction of tree height definded per PFT
c   propn of canopy burnt = (SH - (height - cl))/cl = (SH - height + cl)/cl !A


        if (sh(pft).lt.(height_class(class,pft)-cl_t(class,pft)))
     *            ck_t(class,pft)=0.0

        if (sh(pft).ge.(height_class(class,pft)-cl_t(class,pft))
     *       .and.sh(pft).lt.height_class(class,pft))
     *   ck_t(class,pft)=
     *   ((sh(pft)-height_class(class,pft)+cl_t(class,pft))
     *          /cl_t(class,pft))

          if (sh(pft).ge.height_class(class,pft)) ck_t(class,pft)=1.0

          ck(pft)=ck(pft)+ck_t(class,pft)
          
c       post-fire mortality from crown scorching
          pm_ck(pft)=pm_ck(pft)+(r_ck(pft)*(ck_t(class,pft)**p(pft)))

c       Doug 02/13: removed post-fire mortality from cambial damage
c         Allan's version after Peterson&Ryan
c       Doug 02/13: rplaced lots with this. Do a better comment at some point.
          CALL BT_change(dbh_class(class,pft),
     *      BTparam1(pft,1:3),BTparam2(pft,1:3),
     *      tau_l(pft),pm_tau_class)
       
c          if ((tau_l(pft)/tau_c(class,pft)).ge.2.0) then
c               pm_tau_class(class,pft)=1.0
c          else
c             if ((tau_l(pft)/tau_c(class,pft)).gt.0.22) then
c           pm_tau_class(class,pft)=
c     *              (0.563*(tau_l(pft)/tau_c(class,pft)))-0.125
c             else
c               pm_tau_class(class,pft)=0.0
c             endif
c          endif

          pm_tau(pft)=pm_tau(pft)+pm_tau_class(class,pft)

        enddo !height class
        
          ck(pft)=ck(pft)/5.0
          pm_ck(pft)=pm_ck(pft)/5.0
          pm_tau(pft)=pm_tau(pft)/5.0
           
c       Calculate total post-fire mortality from crown scorching AND cambial kill

       postf_mort(pft)=pm_tau(pft)+pm_ck(pft)-(pm_tau(pft)*pm_ck(pft))

c     number of indivs affected by fire in grid cell !A
        nind_fa(pft)=fire_frac(d)*nind_0(pft)

c     ACHTUNG: moisture content of live fuel can be calculated depending from
c     canopy conductance etc.
*              wc_live(pft)=1.0       ! [fraction]

c          amount of LIVE FUEL combusted in a crown fire
c     ALLAN: Assume 100% of leaves, 100% of small twigs, 100% of large twigs
c     and 5% of small branches combusted (Stocks et al 2004 etc).
        disturb_flux(1)=0.0
        disturb_flux(2)=0.0
        disturb_flux(3)=0.0

        ag_non_leaf_ind(pft,1) = sm_ind(pft,1) + hm_ind(pft,1)
        if (ag_non_leaf_ind(pft,1) .gt. 0.0) then
           do nc=2,nco2
              ag_non_leaf_ind(pft,nc)=(sm_ind(pft,1)*sm_ind(pft,nc)
     *        +hm_ind(pft,1)*hm_ind(pft,nc))/ag_non_leaf_ind(pft,1)
           enddo
        endif

c     Kirsten: the fractions of hm_ind and sm_ind to fuel classes should always 
c       add up to 1.0.
c       5% of 1000hr fuel involved in crown kill

       disturb_flux(1) = nind_fa(pft) * ck(pft) * (lm_ind(pft,1) +
     *    0.045 * ag_non_leaf_ind(pft,1) +
     *    0.075 * ag_non_leaf_ind(pft,1) +
     *    0.21 * 0.05 * ag_non_leaf_ind(pft,1)) 
c     *    0.21 * 0.2 * ag_non_leaf_ind(pft)) 
        do nc=2,nco2
           disturb_flux(nc)=(nind_fa(pft) * ck(pft) * (lm_ind(pft,1)
     *          *lm_ind(pft,nc) +
     *    0.045 * ag_non_leaf_ind(pft,1)*ag_non_leaf_ind(pft,nc) +
     *    0.075 * ag_non_leaf_ind(pft,1)*ag_non_leaf_ind(pft,nc) +
     *    0.21 * 0.05 * ag_non_leaf_ind(pft,1)*ag_non_leaf_ind(pft,1)))
     *          /disturb_flux(1)
        enddo

c      check dcflux_fire_pft(d,pft)=disturb
       temp=dcflux_fire_pft(d,pft,1)
       dcflux_fire_pft(d,pft,1)=dcflux_fire_pft(d,pft,1)+disturb_flux(1)
       if (dcflux_fire_pft(d,pft,1) .gt. 0.) then
          do nc=2,nco2
             dcflux_fire_pft(d,pft,nc)=(temp*
     *         dcflux_fire_pft(d,pft,nc)+disturb_flux(1)
     *         *disturb_flux(nc))/dcflux_fire_pft(d,pft,1)             
          enddo
       endif

       fc_crown(d)=fc_crown(d)+disturb_flux(1)

c     Calculate total number of indivs killed for fire affected area!!!

C       nind_kill(pft) = postf_mort(pft) * nind_fa(pft)
C	  Doug 12/12: resprouting pft, fraction of resprouters suriving fire is 
        nind_kill(pft) = postf_mort(pft) * nind_fa(pft) * 
     *    (1-pftpar(pft,59))
 	 
        nind_resp(pft) = postf_mort(pft) * nind_fa(pft) * 
     *    (pftpar(pft,59))


c     Send a/g non-combusted biomass of dead trees to fuel cats & a/g litter !A
c       But make them available for burning in the next year ! Kirsten


 
C       temp=fuel_update_1hr(pft,1)
C       fuel_update_1hr(pft,1) = fuel_update_1hr(pft,1) + 
C     *                (nind_kill(pft) * (1-ck(pft)) *
C     *                (lm_ind(pft,1) + 0.045 * ag_non_leaf_ind(pft,1)))
C	 
C
C	 
C       if (fuel_update_1hr(pft,1) .gt. 0.) then
C          do nc=2,nco2
C             fuel_update_1hr(pft,nc) = (fuel_update_1hr(pft,nc)*temp + 
C     *            (nind_kill(pft) * (1-ck(pft)) *
C     *            (lm_ind(pft,1)*lm_ind(pft,nc) + 0.045 
C     *            * ag_non_leaf_ind(pft,1)*ag_non_leaf_ind(pft,nc))))
C     *            /fuel_update_1hr(pft,1)
C
C             dfuel_update_1hr(pft,nc)=nind_kill(pft) * (1-ck(pft)) *    !Doug 06/09
C     *            (lm_ind(pft,1)*lm_ind(pft,nc) + 0.045 
C     *            * ag_non_leaf_ind(pft,1)*ag_non_leaf_ind(pft,nc))
C          enddo
C       endif
 
       temp=fuel_update_1hr_leaf(pft,1)
       fuel_update_1hr_leaf(pft,1) = fuel_update_1hr_leaf(pft,1) + 
     *    ((nind_kill(pft)+nind_resp(pft)) *   !Doug: All-but killed resprouters loose remainder of leaves
     *    (1-ck(pft)) * lm_ind(pft,1))
	 
       dfuel_update_1hr_leaf(pft,1)=fuel_update_1hr_leaf(pft,1)-temp	!Doug 03/09, find daily update
	 
       if (fuel_update_1hr_leaf(pft,1) .gt. 0.) then
          do nc=2,nco2
             fuel_update_1hr_leaf(pft,nc) =
     *            (fuel_update_1hr_leaf(pft,nc)*temp + 
     *            ((nind_kill(pft)+nind_resp(pft)) * 
     *            (1-ck(pft)) *
     *            (lm_ind(pft,1)*lm_ind(pft,nc))))
     *            /fuel_update_1hr_leaf(pft,1)
          enddo
       endif
 
       temp=fuel_update_1hr_wood(pft,1)
       fuel_update_1hr_wood(pft,1) = fuel_update_1hr_wood(pft,1) + 
     *     (nind_kill(pft) * (1-ck(pft)) *
     *     (0.045 * ag_non_leaf_ind(pft,1)))
	 
       dfuel_update_1hr_wood(pft,1)=fuel_update_1hr_wood(pft,1)-temp	!Doug 03/09, find daily update	 
	 
       if (fuel_update_1hr_wood(pft,1) .gt. 0.) then
          do nc=2,nco2
             fuel_update_1hr_wood(pft,nc) = 
     *            (fuel_update_1hr_wood(pft,nc)*temp + 
     *            (nind_kill(pft) * (1-ck(pft)) *
     *            (0.045 
     *            * ag_non_leaf_ind(pft,1)*ag_non_leaf_ind(pft,nc))))
     *            /fuel_update_1hr_wood(pft,1)
          enddo
       endif
	   

       temp=fuel_update_10hr(pft,1)
       fuel_update_10hr(pft,1) = fuel_update_10hr(pft,1) + 
     *               (nind_kill(pft) * (1-ck(pft)) *
     *               (0.075 * ag_non_leaf_ind(pft,1)))

       dfuel_update_10hr(pft,1)=fuel_update_10hr(pft,1)-temp	!Doug 05/09


       if (fuel_update_10hr(pft,1) .gt. 0.) then
          do nc=2,nco2
             fuel_update_10hr(pft,nc)=(fuel_update_10hr(pft,nc)*temp + 
     *            (nind_kill(pft) * (1-ck(pft)) *
     *            (0.075 * ag_non_leaf_ind(pft,1)*
     *            ag_non_leaf_ind(pft,nc))))/fuel_update_10hr(pft,1)

             dfuel_update_10hr(pft,nc)=nind_kill(pft) * (1-ck(pft)) *	!Doug 06/09
     *            (0.075 * ag_non_leaf_ind(pft,1)*
     *            ag_non_leaf_ind(pft,nc))
          enddo
       endif

c     cambial damage and left-overs from crown scorching

       temp=fuel_update_100hr(pft,1)
       fuel_update_100hr(pft,1)  = fuel_update_100hr(pft,1) +  
     *      (nind_kill(pft) * (1-ck(pft)) * 
     *      (0.21 * ag_non_leaf_ind(pft,1)) +  
     *      nind_kill(pft) * ck(pft) *
     *      (0.95 * 0.21 * ag_non_leaf_ind(pft,1)))
c     *               (0.8 * 0.21 * ag_non_leaf_ind(pft)))

       dfuel_update_100hr(pft,1)=fuel_update_100hr(pft,1)-temp	!Doug 05/09


       if (fuel_update_100hr(pft,1) .gt. 0.) then
          do nc=2,nco2
             fuel_update_100hr(pft,nc)=(fuel_update_100hr(pft,nc)
     *            *temp + (nind_kill(pft) * (1-ck(pft)) * 
     *            (0.21*ag_non_leaf_ind(pft,1)*ag_non_leaf_ind(pft,nc))
     *            + nind_kill(pft)*ck(pft)*(0.95 * 0.21 * 
     *            ag_non_leaf_ind(pft,1)*ag_non_leaf_ind(pft,nc))))
     *            /fuel_update_100hr(pft,1)

             dfuel_update_100hr(pft,nc)=nind_kill(pft) * (1-ck(pft)) * 	!Doug 06/09
     *            (0.21*ag_non_leaf_ind(pft,1)*ag_non_leaf_ind(pft,nc))
     *            + nind_kill(pft)*ck(pft)*(0.95 * 0.21 * 
     *            ag_non_leaf_ind(pft,1)*ag_non_leaf_ind(pft,nc))
          enddo
       endif


c      cambial damage and crown scorching together

       temp=fuel_update_1000hr(pft,1)
       fuel_update_1000hr(pft,1) = fuel_update_1000hr(pft,1) + 

     *                 nind_kill(pft) * (0.67 * ag_non_leaf_ind(pft,1))

       dfuel_update_1000hr(pft,1)=fuel_update_1000hr(pft,1)-temp	!Doug 05/09
 

       if (fuel_update_1000hr(pft,1) .gt. 0.) then
          do nc=2,nco2
             fuel_update_1000hr(pft,nc) = (fuel_update_1000hr(pft,nc)
     *            *temp + nind_kill(pft) * (0.67 * 
     *            ag_non_leaf_ind(pft,1)*ag_non_leaf_ind(pft,nc)))
     *            /fuel_update_1000hr(pft,1)

             dfuel_update_1000hr(pft,nc)=nind_kill(pft) * (0.67 * 	!Doug 06/09
     *            ag_non_leaf_ind(pft,1)*ag_non_leaf_ind(pft,nc))
          enddo
       endif

c     Send roots of dead trees to b/g litter, Kirsten: this can be done daily, no effect on fire spread

       temp=litter_bg(pft,1)
       litter_bg(pft,1)=litter_bg(pft,1)+rm_ind(pft,1)*nind_kill(pft)
       if (litter_bg(pft,1) .gt. 0.0) then
          do nc=2,nco2
             litter_bg(pft,nc)=(litter_bg(pft,nc)*temp+rm_ind(pft,1)
     *            *rm_ind(pft,nc)*nind_kill(pft))/litter_bg(pft,1)
          enddo
       endif


c       update individual density and send dead biomass as a result of postfire mort. to
c       litter and fuel classes, resp. !A
        
c    add DEAD FUEL combustion to fire carbon flux per PFT

        temp=dcflux_fire_pft(d,pft,1)
        dcflux_fire_pft(d,pft,1)=dcflux_fire_pft(d,pft,1)+
     *   ((fc_1hr_leaf(pft,1)+fc_1hr_wood(pft,1)+
     *       fc_10hr(pft,1)+fc_100hr(pft,1))*0.45)
     *       +fc_1000hr(pft,1)
       if (dcflux_fire_pft(d,pft,1) .gt. 0.0) then
          do nc=2,nco2
        dcflux_fire_pft(d,pft,nc)=(dcflux_fire_pft(d,pft,nc)*temp+
     *   ((fc_1hr_leaf(pft,1)*fc_1hr_leaf(pft,nc)
     *   +fc_1hr_wood(pft,1)*fc_1hr_wood(pft,nc)
     *   +fc_10hr(pft,1)*fc_10hr(pft,nc)
     *   +fc_100hr(pft,1)*fc_100hr(pft,nc))*0.45)
     *   +fc_1000hr(pft,1)*fc_1000hr(pft,nc))/dcflux_fire_pft(d,pft,1)
          enddo
       endif


c    total carbon flux from life fuel and dead fuel combustion

       temp=mcflux_fire_pft(m,pft,1)
       mcflux_fire_pft(m,pft,1)=mcflux_fire_pft(m,pft,1)+
     *      dcflux_fire_pft(d,pft,1)
       if (mcflux_fire_pft(m,pft,1) .gt. 0.0) then
          do nc=2,nco2
             mcflux_fire_pft(m,pft,nc)=(mcflux_fire_pft(m,pft,nc)*temp+
     *         dcflux_fire_pft(d,pft,1)*dcflux_fire_pft(d,pft,nc))/
     *         mcflux_fire_pft(m,pft,1)
          enddo
       endif
       
       temp=acflux_fire_pft(pft,1)
       acflux_fire_pft(pft,1)=acflux_fire_pft(pft,1)
     *      +dcflux_fire_pft(d,pft,1)
       if (acflux_fire_pft(pft,1) .gt. 0.0) then
          do nc=2,nco2
             acflux_fire_pft(pft,nc)=(acflux_fire_pft(pft,nc)*temp+
     *         dcflux_fire_pft(d,pft,1)*dcflux_fire_pft(d,pft,nc))/
     *         acflux_fire_pft(pft,1)
          enddo
       endif

c     Number indivs NOT killed by fire
       nind_nkill(pft) = nind(pft) - nind_kill(pft)

c     ALLAN: If at least partial crown kill occurs, and some fire-affected
c       trees survive fire, and 
c     number of trees not killed by fire is non-zero, then redistribute
c     leaf, sapwood, and heartwood and rootmass among survivors

c   Kirsten: error in setting the brackets: the last one wasn't seen by the
c   compiler. Check: performance ok now?

c   Doug 12/12 RS: if resprouter, send suriving trees into normal pft.

      if ((ck(pft).gt.0.0).and.
     *  ((nind_fa(pft)-nind_kill(pft)).gt.0.0) 
     *  .and. (nind_nkill(pft) .gt. 0.0)) then

c       prop_fa_nk(pft) = (nind_fa(pft)-nind_kill(pft))/nind_nkill(pft)
       prop_fa_rs(pft) = nind_resp(pft)/nind_nkill(pft)
       
       prop_fa_nk(pft) = (nind_fa(pft)-nind_kill(pft)-nind_resp(pft))
     *   /nind_nkill(pft)

       lm_ind(pft,1)=lm_ind(pft,1)-prop_fa_nk(pft)*ck(pft)*
     *   lm_ind(pft,1)-prop_fa_rs(pft)*lm_ind(pft,1)

       sm_ind(pft,1) = sm_ind(pft,1) - prop_fa_nk(pft) * ck(pft) *
     *   (0.045 + 0.075 + (0.21 * 0.05)) * sm_ind(pft,1)
     *   -prop_fa_rs(pft)*(0.045 + 0.075 + (0.21 * 0.05)) * 
     *   sm_ind(pft,1)
c     *   (0.045 + 0.075 + (0.21 * 0.2)) * sm_ind(pft)

       hm_ind(pft,1) = hm_ind(pft,1) - prop_fa_nk(pft) * ck(pft) *
     *   (0.045 + 0.075 + (0.21 * 0.05)) * hm_ind(pft,1)
     *   -prop_fa_rs(pft)*(0.045 + 0.075 + (0.21 * 0.05)) * 
     *   hm_ind(pft,1)
c     *   (0.045 + 0.075 + (0.21 * 0.2)) * hm_ind(pft)

cc   KIRSTEN: where is the rootmass??? No rootmass treatment, otherwise carbon balance not closed!!!!!
c       rm_ind(pft) = rm_ind(pft) - prop_fa_nk(pft)*ck(pft)*rm_ind(pft)


       endif 

c     ALLAN:  Assume the proportional decrease in individual root mass due 
c     to stress of surviving a crown-fire among fire-affected indivs = 0%. 
c     So, rm_ind(pft) stays unchanged.  
c     Otherwise, implement proportional reduction (z) = 0.05, say:
c     rm_ind(pft) = rm_ind(pft) -  prop_fa_nk(pft) * ck(pft) * rm_ind(pft) * z
c     litter_bg(pft) = litter_bg(pft) + 
c     *  prop_fa_nk(pft) * ck(pft) * rm_ind(pft) * z


c     Update number of indivs surviving fire, ready for the next day. !A
       nind(pft)=nind_nkill(pft)

      else !grass

c       dead grass

        temp=dcflux_fire_pft(d,pft,1) 
        dcflux_fire_pft(d,pft,1)=dcflux_fire_pft(d,pft,1)
     *       +((fc_1hr_leaf(pft,1)+fc_1hr_wood(pft,1))*0.45)
        if (dcflux_fire_pft(d,pft,1) .gt. 0.0) then
          do nc=2,nco2
             dcflux_fire_pft(d,pft,nc)=(dcflux_fire_pft(d,pft,nc)*temp+
     *         (fc_1hr_leaf(pft,1)*fc_1hr_leaf(pft,nc)
     *         +fc_1hr_wood(pft,1)*fc_1hr_wood(pft,nc))*0.45)/
     *         dcflux_fire_pft(d,pft,1)
          enddo
       endif

c      live grass        

        temp=dcflux_fire_pft(d,pft,1) 
        dcflux_fire_pft(d,pft,1)=dcflux_fire_pft(d,pft,1)+fc_lg(pft,1)
        if (dcflux_fire_pft(d,pft,1) .gt. 0.0) then
          do nc=2,nco2
             dcflux_fire_pft(d,pft,nc)=(dcflux_fire_pft(d,pft,nc)*temp+
     *         fc_lg(pft,1)*fc_lg(pft,nc))/
     *         dcflux_fire_pft(d,pft,1)
          enddo
       endif

c  KIRSTEN: no update of rootmass and litter_bg here, otherwise carbon balance not closed!!
       lm_ind(pft,1)=lm_ind(pft,1)-fc_lg(pft,1)

c    Kirsten: update of livegrass only through lm_ind, due to combination of dphen and livegrass 
c       consumption
c       livegrass=livegrass-fc_lg(pft) !A
c       pot_fc_lg_temp(pft)=lm_ind(pft)-fc_lg(pft)
 
c     Kirsten: never forget to protect variables aganst division by zero!!
c        ratio_lg_reduc=0.0
c        if (pot_fc_lg_temp(pft).gt.0.0)
c     *      ratio_lg_reduc=(lm_ind(pft)/pot_fc_lg_temp(pft))-1.0

c  KIRSTEN: no update of rootmass and litter_bg here, otherwise carbon balance not closed!!
c        rm_ind(pft)=rm_ind(pft)-(ratio_lg_reduc*rm_ind(pft))
c        litter_bg(pft)=litter_bg(pft)+(ratio_lg_reduc*rm_ind(pft))

c       Kirsten: this way we loose all the roots
c       lm_ind(pft)=pot_fc_lg_temp(pft)*0.45 !A
c       litter_bg(pft)=litter_bg(pft) + ratio_lg_reduc * rm_ind(pft)!A
c       rm_ind(pft)=(1.0-ratio_lg_reduc) * rm_ind(pft)!A



c    total carbon flux from life fuel and dead fuel combustion

       temp=mcflux_fire_pft(m,pft,1)
       mcflux_fire_pft(m,pft,1)=mcflux_fire_pft(m,pft,1)+
     *      dcflux_fire_pft(d,pft,1)
       if (mcflux_fire_pft(m,pft,1) .gt. 0.0) then
          do nc=2,nco2
             mcflux_fire_pft(m,pft,nc)=(mcflux_fire_pft(m,pft,nc)*temp+
     *         dcflux_fire_pft(d,pft,1)*dcflux_fire_pft(d,pft,nc))/
     *         mcflux_fire_pft(m,pft,1)
          enddo
       endif
       
       temp=acflux_fire_pft(pft,1)
       acflux_fire_pft(pft,1)=acflux_fire_pft(pft,1)
     *      +dcflux_fire_pft(d,pft,1)
       if (acflux_fire_pft(pft,1) .gt. 0.0) then
          do nc=2,nco2
             acflux_fire_pft(pft,nc)=(acflux_fire_pft(pft,nc)*temp+
     *         dcflux_fire_pft(d,pft,1)*dcflux_fire_pft(d,pft,nc))/
     *         acflux_fire_pft(pft,1)
          enddo
       endif

      endif   !if tree or grass

      endif !i_surface.gt.50 kW/m: combust only when fire spread allowed


      endif     !present

      enddo     !fire effects per PFT

c   KIRSTEN:  take proportion of grass leafmass, when green grass leaves are on
c         todays amount of green grass leaves: [gC/m2],influence on ROS only through moist_lg_1hr

c     *                        (lm_ind(pft)/0.45*nind(pft))  !g/m2

c              used in fire effects section only

c               pot_fc_lg(pft)=(lm_ind(pft)/0.45*nind(pft))  !g/m2 

c       net fuel load
      m_fc_crown(m)=m_fc_crown(m)+fc_crown(d)
      
      an_fc_crown=an_fc_crown+fc_crown(d)

201       continue
c    Kirsten:  adding up of daily to monthly values
         if (d.eq.month(m)) then
            if (count.gt.1) then
               mfdi(m)=mfdi(m)/count
               count_fdi=count_fdi+count
            endif
            if (count_int.gt.1) then
               m_i_surface(m)=m_i_surface(m)/count_int
               count_yr=count_yr+count_int
C Yan
               m_i_surface_human(m)=m_i_surface_human(m)/count_int
               m_i_surface_lightn(m)=m_i_surface_lightn(m)/count_int
            endif
            m=m+1


            count_int=0
            count=0
          elseif (d.le.month(m).and.afire_frac.eq.1.0) then
            mfdi(m)=mfdi(m)/count
            count_fdi=count_fdi+count
            m_i_surface(m)=m_i_surface(m)/count_int
C Yan
            m_i_surface_human(m)=m_i_surface_human(m)/count_int
            m_i_surface_lightn(m)=m_i_surface_lightn(m)/count_int
            count_yr=count_yr+count_int
            count=0
            count_int=0
          endif

        fuel_1hr_leaf_left(:,1)=fuel_1hr_leaf_left(:,1)/0.45			!Doug 03/01
        fuel_1hr_wood_left(:,1)=fuel_1hr_wood_left(:,1)/0.45			!Doug 03/01
        fuel_10hr_left(:,1)=fuel_10hr_left(:,1)/0.45
        fuel_100hr_left(:,1)=fuel_100hr_left(:,1)/0.45
        fuel_1000hr_left(:,1)=fuel_1000hr_left(:,1)/0.45

c     stop daily loop if entire grid cell burnt
         if (afire_frac.eq.1.0) then
           goto 200
          end if

c       pause

          mfuel_1hr_total(min(m,12))=mfuel_1hr_total(min(m,12))+              !Doug 03/09: monthly fuel
     *      (fuel_1hr_total*(1-fire_frac(d)))/month_length(min(m,12)) !for outputtng only
          mfuel_1hr_leaf_total(min(m,12))=
     *      mfuel_1hr_leaf_total(min(m,12))+              !Doug 03/09: monthly fuel
     *      (fuel_1hr_leaf_total*(1-fire_frac(d)))/
     *       month_length(min(m,12)) !for outputtng only
          mfuel_1hr_wood_total(min(m,12))=
     *      mfuel_1hr_wood_total(min(m,12))+              !Doug 03/09: monthly fuel
     *      (fuel_1hr_wood_total*(1-fire_frac(d)))/
     *      month_length(min(m,12)) !for outputtng only

          mfuel_10hr_total(min(m,12))=mfuel_10hr_total(min(m,12))+
     *      (fuel_10hr_total*(1-fire_frac(d)))/month_length(min(m,12))

          mfuel_100hr_total(min(m,12))=mfuel_100hr_total(min(m,12))+
     *      (fuel_100hr_total*(1-fire_frac(d)))/month_length(min(m,12))

          mfuel_1000hr_total(min(m,12))=mfuel_1000hr_total(min(m,12))+
     *      (fuel_1000hr_total*(1-fire_frac(d)))/month_length(min(m,12))

          mlivegrass(min(m,12))=mlivegrass(min(m,12))+
     *      (livegrass*(1-fire_frac(d)))/month_length(min(m,12))

       enddo !daily time step

        fuel_1hr_leaf_left(:,1)=fuel_1hr_leaf_left(:,1)*0.45			!Doug 03/01
        fuel_1hr_wood_left(:,1)=fuel_1hr_wood_left(:,1)*0.45			!Doug 03/01
        fuel_10hr_left(:,1)=fuel_10hr_left(:,1)*0.45
        fuel_100hr_left(:,1)=fuel_100hr_left(:,1)*0.45
        fuel_1000hr_left(:,1)=fuel_1000hr_left(:,1)*0.45

        DO pft=1,npft
          IF (present(pft)) THEN
            temp=fuel_1hr_leaf_left(pft,1)

            fuel_1hr_leaf_left(pft,1)=fuel_1hr_leaf_left(pft,1)*
     *        (1-fire_frac(365))+
     *        fuel_left_minus(pft,1)+
     *        dfuel_update_1hr_leaf(pft,1)/0.45

            IF (fuel_1hr_leaf(pft,1)>0.0) THEN
              DO nc=2,nco2
                fuel_1hr_leaf_left(pft,nc)=(fuel_1hr_leaf_left(pft,nc)
     *            *temp+
     *            fuel_left_minus(pft,1)*
     *            fuel_1hr_leaf_left(pft,nc)+
     *            dfuel_update_1hr_leaf(pft,1)*
     *            dfuel_update_1hr_leaf(pft,nc)/0.45)/
     *            fuel_1hr_leaf_left(pft,1)
              END DO
            END IF
			
            temp=fuel_1hr_wood_left(pft,1)

            fuel_1hr_wood_left(pft,1)=fuel_1hr_wood_left(pft,1)*
     *        (1-fire_frac(365))+
     *        fuel_left_minus(pft,2)+
     *        dfuel_update_1hr_wood(pft,1)/0.45

            IF (fuel_1hr_wood(pft,1)>0.0) THEN
              DO nc=2,nco2
                fuel_1hr_wood_left(pft,nc)=(fuel_1hr_wood_left(pft,nc)
     *            *temp+
     *            fuel_left_minus(pft,2)*
     *            fuel_1hr_wood_left(pft,nc)+
     *            dfuel_update_1hr_wood(pft,1)*
     *            dfuel_update_1hr_wood(pft,nc)/0.45)/
     *            fuel_1hr_wood_left(pft,1)
              END DO
            END IF
	 
C            temp=fuel_1hr_left(pft,1)
C 
C            fuel_1hr_left(pft,1)=fuel_1hr_left(pft,1)*
C     *        (1-fire_frac(365))+
C     *        fuel_left_minus(pft,1)+
C     *        dfuel_update_1hr(pft,1)/0.45

C            IF (fuel_1hr(pft,1)>0.0) THEN
C              DO nc=2,nco2
C                fuel_1hr_left(pft,nc)=(fuel_1hr_left(pft,nc)
C     *            *temp+
C     *            fuel_left_minus(pft,1)*
C     *            fuel_1hr_left(pft,nc)+
C     *            dfuel_update_1hr(pft,1)*
C     *            dfuel_update_1hr(pft,nc)/0.45)/
C     *            fuel_1hr_left(pft,1)
C              END DO
C            END IF


            temp=fuel_10hr_left(pft,1)
            fuel_10hr_left(pft,1)=fuel_10hr_left(pft,1)*
     *        (1-fire_frac(365))+
     *        fuel_left_minus(pft,3)+
     *        dfuel_update_10hr(pft,1)/0.45

            IF (fuel_10hr(pft,1)>0.0) THEN
              DO nc=2,nco2
                fuel_10hr_left(pft,nc)=(fuel_10hr_left(pft,nc)*
     *            *temp+
     *            fuel_left_minus(pft,3)*fuel_10hr_left(pft,nc)
     *            +dfuel_update_10hr(pft,1)*
     *            dfuel_update_10hr(pft,nc)/0.45)/
     *            fuel_10hr_left(pft,1)
               END DO
             END IF


            temp=fuel_100hr_left(pft,1)
            fuel_100hr_left(pft,1)=fuel_100hr_left(pft,1)*
     *        (1-fire_frac(365))+
     *        fuel_left_minus(pft,4)+
     *        dfuel_update_100hr(pft,1)/0.45

            IF (fuel_100hr(pft,1)>0.0) THEN
              DO nc=2,nco2
                fuel_100hr_left(pft,nc)=(fuel_100hr_left(pft,nc)*
     *            *temp+
     *            fuel_left_minus(pft,4)*fuel_100hr_left(pft,nc)
     *            +dfuel_update_100hr(pft,1)*
     *            dfuel_update_100hr(pft,nc)/0.45)/
     *            fuel_100hr_left(pft,1)
               END DO
             END IF


            temp=fuel_1000hr_left(pft,1)
            fuel_1000hr_left(pft,1)=fuel_1000hr_left(pft,1)*
     *        (1-fire_frac(365))+
     *        fuel_left_minus(pft,5)+
     *        dfuel_update_1000hr(pft,1)/0.45

            IF (fuel_10hr(pft,1)>0.0) THEN
              DO nc=2,nco2
                fuel_1000hr_left(pft,nc)=(fuel_1000hr_left(pft,nc)*
     *            *temp+
     *            fuel_left_minus(pft,5)*fuel_1000hr_left(pft,nc)
     *            +dfuel_update_1000hr(pft,1)*
     *            dfuel_update_1000hr(pft,nc)/0.45)/
     *            fuel_1000hr_left(pft,1)
               END DO
             END IF

          END IF
        END DO



200     continue

        if(lat.ge.0.0)then !in the North
          do d=91,365 !fire season starts from April 1.
           afire_frac_afap_old=afire_frac_afap_old+fire_frac(d)
           end do
        else !in the South
          do d=274,365 !fire season starts from Oct. 1
           afire_frac_afap_old=afire_frac_afap_old+fire_frac(d)
          end do
        end if

       afire_frac=an_areafires/area_ha
		
       do m=1,12
         mfire_frac(m)=area_burnt(m)/area_ha
		
         if(mfire_frac(m).lt.0.0001)then
          mfire_frac(m)=0.0
          area_burnt(m)=0.0
         end if
       end do
         
       if(afire_frac.lt.0.0001)then
           afire_frac=0.0
           an_areafires=0.0
       end if


c       Doug 04/09: pasture removed from scaling down of fire, as
c       data-comparision shows this is not the correct mechanism
c       for treating pasture

c             IF(crop>0.0) THEN
               an_areafires=an_areafires*(1.0-crop)

               acflux_fire(1)=acflux_fire(1)*(1.0-crop)

               afire_frac=afire_frac*(1.0-crop)

               livegrass_0=livegrass_0*(1.0-crop)

               dead_fuel_all_0=dead_fuel_all_0*(1.0-crop)

               dead_fuel_0=dead_fuel_0*(1.0-crop)

               fuel_all_0=fuel_all_0*(1.0-crop)

               fuel_1hr_total_0=fuel_1hr_total_0*(1.0-crop)

               fuel_10hr_total_0=fuel_10hr_total_0*(1.0-crop)

               fuel_100hr_total_0=fuel_100hr_total_0*(1.0-crop)

               fuel_1000hr_total_0=fuel_1000hr_total_0*(1.0-crop)


               do l=1,12

                 area_burnt(l)=area_burnt(l)*(1.0-crop)
                 mcflux_fire(l,1)=mcflux_fire(l,1)*(1.0-crop)
                 mfire_frac(l)=mfire_frac(l)*(1.0-crop)
                 num_fire(l)=num_fire(l)*(1.0-crop)
               end do

c             end if

c LeiLei 12/08: experiment in masking out fires on areas fragmented
c by cropland and/or pasture by setting fuel (and thus fires) to zero
c in areas above the specified threshold.
c
c       Doug 04/09: again, poas and crop is now scalar instead of array
c       so the old search method for correct cell has been removed

            IF(crop+pas>10) THEN				

               an_areafires=0.0
               acflux_fire(1)=0.0
               afire_frac=0.0
               livegrass_0=0.0
               dead_fuel_all_0=0.0
               dead_fuel_0=0.0
               fuel_all_0=0.0
               fuel_1hr_total_0=0.0
               fuel_10hr_total_0=0.0
               fuel_100hr_total_0=0.0
               fuel_1000hr_total_0=0.0

               do l=1,12
                 area_burnt(l)=0.0
                 mcflux_fire(l,1)=0.0
                 mfire_frac(l)=0.0
                 num_fire(l)=0.0
               end do

            end if

           if(acflux_fire(1).le.0.1)acflux_fire(1)=0.0



c  send dead biomass from cambial damage to fuel classes at the end of the year
c  to avoid that the dead fuel influence fire spread of the actual fire season (standing dead biomass)
       DO pft=1,npft
         IF (present(pft)) THEN
		 
           temp=fuel_1hr_leaf(pft,1)
           fuel_1hr_leaf(pft,1)=fuel_1hr_leaf_0(pft,1)*(1.0-afire_frac)
     *                  +fuel_update_1hr_leaf(pft,1)
     *                  -(fuel_load_subtraction(pft,1))*
     *                  f1hr_leaf_frac(pft)
	 
         IF (fuel_1hr_leaf(pft,1) .gt. 0.0) THEN
           DO nc=2,nco2
              fuel_1hr_leaf(pft,nc)=(fuel_1hr_leaf(pft,nc)*temp+
     *          fuel_update_1hr_leaf(pft,1)*
     *          fuel_update_1hr_leaf(pft,nc)) !Doug 11/12: needs checking
     *          /fuel_1hr_leaf(pft,1)
           END DO
         END IF
		 
           temp=fuel_1hr_wood(pft,1)
           fuel_1hr_wood(pft,1)=fuel_1hr_wood_0(pft,1)*(1.0-afire_frac)
     *                  +fuel_update_1hr_wood(pft,1)
     *                  -(fuel_load_subtraction(pft,1))*
     *                  (1-f1hr_leaf_frac(pft))
         IF (fuel_1hr_wood(pft,1) .gt. 0.0) THEN
           DO nc=2,nco2
              fuel_1hr_wood(pft,nc)=(fuel_1hr_wood(pft,nc)*temp+
     *          fuel_update_1hr_wood(pft,1)*
     *          fuel_update_1hr_wood(pft,nc)) !Doug 11/12: needs checking
     *          /fuel_1hr_wood(pft,1)
           END DO
         END IF
		
C        temp=fuel_1hr(pft,1)
C        fuel_1hr(pft,1)=fuel_1hr_0(pft,1)*(1.0-afire_frac)
C     *                  +fuel_update_1hr(pft,1)
C     *                 -fuel_load_subtraction(pft,1)
C        if (fuel_1hr(pft,1) .gt. 0.0) then
C           do nc=2,nco2
C              fuel_1hr(pft,nc)=(fuel_1hr(pft,nc)*temp+
C     *          fuel_update_1hr(pft,1)*fuel_update_1hr(pft,nc))
C     *          /fuel_1hr(pft,1)
C           enddo
C        endif
       temp=fuel_10hr(pft,1)
       fuel_10hr(pft,1)=fuel_10hr_0(pft,1)*(1.0-afire_frac)
     *                  +fuel_update_10hr(pft,1)
     *                  -fuel_load_subtraction(pft,2)
	 
       if (fuel_10hr(pft,1) .gt. 0.0) then
          do nc=2,nco2
             fuel_10hr(pft,nc)=(fuel_10hr(pft,nc)*temp+
     *          fuel_update_10hr(pft,1)*fuel_update_10hr(pft,nc))
     *          /fuel_10hr(pft,1)
          enddo
       endif
 
       temp=fuel_100hr(pft,1)
       fuel_100hr(pft,1)=fuel_100hr_0(pft,1)*(1.0-afire_frac)
     *                    +fuel_update_100hr(pft,1)
     *                    -fuel_load_subtraction(pft,3)
       if (fuel_100hr(pft,1) .gt. 0.0) then
          do nc=2,nco2
             fuel_100hr(pft,nc)=(fuel_100hr(pft,nc)*temp+
     *         fuel_update_100hr(pft,1)*fuel_update_100hr(pft,nc))
     *           /fuel_100hr(pft,1)
          enddo
       endif
 
       temp=fuel_1000hr(pft,1)
       fuel_1000hr(pft,1)=fuel_1000hr_0(pft,1)*(1.0-afire_frac)
     *                    +fuel_update_1000hr(pft,1)
     *                    -fuel_load_subtraction(pft,4)
       if (fuel_1000hr(pft,1) .gt. 0.0) then
          do nc=2,nco2
             fuel_1000hr(pft,nc)=(fuel_1000hr(pft,nc)*temp+
     *          fuel_update_1000hr(pft,1)*fuel_update_1000hr(pft,nc))
     *          /fuel_1000hr(pft,1)
          enddo
       endif
	   
       temp=litter_ag_leaf(pft,1)
       litter_ag_leaf(pft,1) = litter_ag_leaf(pft,1) + 
     *             fuel_update_1hr_leaf(pft,1)
     *             -fuel_load_subtraction(pft,1)*
     *             f1hr_leaf_frac(pft)
       IF (litter_ag_leaf(pft,1) .gt. 0.0) THEN
          DO nc=2,nco2
             litter_ag_leaf(pft,nc)= (litter_ag_leaf(pft,nc)*temp+
     *        fuel_update_1hr_leaf(pft,1)*
     *        fuel_update_1hr_leaf(pft,nc))
     *        /litter_ag_leaf(pft,1)
          END DO
       END IF
	   
       temp=litter_ag_wood(pft,1)
       litter_ag_wood(pft,1) = litter_ag_wood(pft,1) +
     *             fuel_update_1hr_wood(pft,1) +
     *             fuel_update_10hr(pft,1) + fuel_update_100hr(pft,1)+
     *             fuel_update_1000hr(pft,1)
     *             -fuel_load_subtraction(pft,1)*
     *             (1-f1hr_leaf_frac(pft))
     *             -SUM(fuel_load_subtraction(pft,2:4))    
       IF (litter_ag_wood(pft,1) .gt. 0.0) THEN
          DO nc=2,nco2
             litter_ag_wood(pft,nc)= (litter_ag_wood(pft,nc)*temp
     *        + fuel_update_1hr_wood(pft,1)*fuel_update_1hr_wood(pft,nc)
     *        + fuel_update_10hr(pft,1)*fuel_update_10hr(pft,nc)
     *        + fuel_update_100hr(pft,1)*fuel_update_100hr(pft,nc)
     *        + fuel_update_1000hr(pft,1)*fuel_update_1000hr(pft,nc))
     *        /litter_ag_wood(pft,1)
          END DO
       END IF
 
c       temp=litter_ag(pft,1)
c       litter_ag(pft,1) = litter_ag(pft,1) + fuel_update_1hr(pft,1) +
c     *             fuel_update_10hr(pft,1) + fuel_update_100hr(pft,1)+
c     *             fuel_update_1000hr(pft,1)
c     *             -SUM(fuel_load_subtraction(pft,:))    
c       if (litter_ag(pft,1) .gt. 0.0) then
c          do nc=2,nco2
c             litter_ag(pft,nc)= (litter_ag(pft,nc)*temp
c     *        + fuel_update_1hr(pft,1)*fuel_update_1hr(pft,nc)
c     *        + fuel_update_10hr(pft,1)*fuel_update_10hr(pft,nc)
c     *        + fuel_update_100hr(pft,1)*fuel_update_100hr(pft,nc)
c     *        + fuel_update_1000hr(pft,1)*fuel_update_1000hr(pft,nc))
c     *        /litter_ag(pft,1)
c          enddo
c       endif
      endif
      enddo

c      fuel_1hr=fuel_1hr_left*0.45
c      fuel_10hr=fuel_10hr_left*0.45
c      fuel_100hr=fuel_100hr_left*0.45
c      fuel_1000hr=fuel_1000hr_left*0.45

c    summation of monthly and annual variables

       do m=1,12 !A
        do pft=1,npft  !A
           temp=mcflux_fire(m,1)
           mcflux_fire(m,1)=mcflux_fire(m,1)+mcflux_fire_pft(m,pft,1)  !A
           if (mcflux_fire(m,1) .gt. 0.0) then
              do nc=2,nco2
                 mcflux_fire(m,nc)=(mcflux_fire(m,nc)*temp+
     *              mcflux_fire_pft(m,pft,1)*mcflux_fire_pft(m,pft,nc)) !A
     *              /mcflux_fire(m,1)
              enddo
           endif
        enddo                   !A
      enddo                     !A

      do pft=1,npft
         temp=acflux_fire(1)
         acflux_fire(1)=acflux_fire(1)+acflux_fire_pft(pft,1) !A
         if (acflux_fire(1) .gt. 0.0) then
            do nc=2,nco2
               acflux_fire(nc)=(acflux_fire(nc)*temp+
     *              acflux_fire_pft(pft,1)*acflux_fire_pft(pft,nc)) !A
     *              /acflux_fire(1)
              enddo
           endif
      enddo

          if(acflux_fire(1).lt.0.1)acflux_fire(1)=0.0

        an_fdi=an_fdi/count_fdi
        an_i_surface=an_i_surface/count_yr

c     Calculate monthly trace gas emission resulting from biomass burning (gSpecies/m?²)


       do m=1,12 !months
        do x=1,6 !trace_species

         do pft=1,npft
         if (present(pft)) then

            if (mcflux_fire_pft(m,pft,1).gt.0.0) then !A 

            mcflux_trace_pft(m,x,pft)=(mcflux_fire_pft(m,pft,1)*
     *                                    ef_trace(pft,x)/0.45)
          mcflux_trace(m,x)=mcflux_trace(m,x)+mcflux_trace_pft(m,x,pft)
           else
               mcflux_trace(m,x)=0.0
           endif

         endif !if present

        enddo
           acflux_trace(x)=acflux_trace(x)+mcflux_trace(m,x)
       enddo
      enddo

      endif !year loop to switch between Glob-FIRM and SPITFIRE


c     ALLAN: For the sake of neatness, account for rounding errors 
c     leading to slightly negative balances.

      DO pft=1,npft
      IF (present(pft)) THEN
        IF (litter_ag_leaf(pft,1).lt.0.0) THEN
          litter_ag_leaf(pft,1)=0.0
          litter_ag_leaf(pft,2)=0.0
          litter_ag_leaf(pft,3)=0.0
        END IF
        IF (litter_ag_wood(pft,1).lt.0.0) THEN
          litter_ag_wood(pft,1)=0.0
          litter_ag_wood(pft,2)=0.0
          litter_ag_wood(pft,3)=0.0
        END IF
		
c        if (litter_ag(pft,1).lt.0.0) then
c          litter_ag(pft,1)=0.0
c          litter_ag(pft,2)=0.0
c          litter_ag(pft,3)=0.0
c        endif
        if (litter_bg(pft,1).lt.0.0) then
          litter_bg(pft,1)=0.0
          litter_bg(pft,2)=0.0
          litter_bg(pft,3)=0.0
        endif
        IF (fuel_1hr_leaf(pft,1).lt.0.0) THEN
          fuel_1hr_leaf(pft,1)=0.0
          fuel_1hr_leaf(pft,2)=0.0
          fuel_1hr_leaf(pft,3)=0.0
        END IF
        IF (fuel_1hr_wood(pft,1).lt.0.0) THEN
          fuel_1hr_wood(pft,1)=0.0
          fuel_1hr_wood(pft,2)=0.0
          fuel_1hr_wood(pft,3)=0.0
        END IF
        if (fuel_10hr(pft,1).lt.0.0) then
          fuel_10hr(pft,1)=0.0
          fuel_10hr(pft,2)=0.0
          fuel_10hr(pft,3)=0.0
        endif
        if (fuel_100hr(pft,1).lt.0.0) then
          fuel_100hr(pft,1)=0.0
          fuel_100hr(pft,2)=0.0
          fuel_100hr(pft,3)=0.0
        endif
        if (fuel_1000hr(pft,1).lt.0.0) then
          fuel_1000hr(pft,1)=0.0
          fuel_1000hr(pft,2)=0.0
          fuel_1000hr(pft,3)=0.0
        endif
        if (lm_ind(pft,1).lt.0.0) then
          lm_ind(pft,1)=0.0
          lm_ind(pft,2)=0.0
          lm_ind(pft,3)=0.0
        endif
        if (sm_ind(pft,1).lt.0.0) then
          sm_ind(pft,1)=0.0
          sm_ind(pft,2)=0.0
          sm_ind(pft,3)=0.0
        endif
        if (hm_ind(pft,1).lt.0.0) then
          hm_ind(pft,1)=0.0
          hm_ind(pft,2)=0.0
          hm_ind(pft,3)=0.0
        endif
        if (rm_ind(pft,1).lt.0.0) then
          rm_ind(pft,1)=0.0
          rm_ind(pft,2)=0.0
          rm_ind(pft,3)=0.0
        endif

       endif
      enddo




         fuel_1hr_total_0=fuel_1hr_total*0.45	! Doug: remove
         fuel_10hr_total_0=fuel_10hr_total*0.45
         fuel_100hr_total_0=fuel_100hr_total*0.45
         fuel_1000hr_total_0=fuel_1000hr_total*0.45


         mfuel_1hr_total=mfuel_1hr_total*.45
         mfuel_1hr_leaf_total=mfuel_1hr_leaf_total*.45
         mfuel_1hr_wood_total=mfuel_1hr_wood_total*.45
         mfuel_10hr_total=mfuel_10hr_total*.45
         mfuel_100hr_total=mfuel_100hr_total*.45
         mfuel_1000hr_total=mfuel_1000hr_total*.45




      return
      end	!FIRE


c----------------------------------------------------------------------
c     Doug 03/09: FUEL 1HR LEAF FALL
C     Redstributes change in 1 hour fuel load according to leaf phenology
C     To do:
C           -remove decay values
c           - include background leaf/twig loss (possibly based on leaf longliverty)


      SUBROUTINE fuel_1hr_redist(dfuel_leaf,
     *    fuel_1hr_leaf_inc_pos,fuel_1hr_leaf_inc_neg,
     *    fuel_1hr_wood_inc_pos,fuel_1hr_wood_inc_neg,
     *    leaf_decayc, dphen, present,
     *    tree, npft, nco2)


      IMPLICIT NONE

c     PARAMETERS
      INTEGER    npft,nco2
      REAL       leaf_decayc(1:npft)

c     ARGUMENTS
      REAL       dfuel_leaf(1:365,1:npft)
      REAL       fuel_1hr_leaf_inc_pos(1:npft,1:nco2,1:12)
      REAL       fuel_1hr_leaf_inc_neg(1:npft,1:nco2,1:12)
      REAL       fuel_1hr_wood_inc_pos(1:npft,1:nco2,1:12)
      REAL       fuel_1hr_wood_inc_neg(1:npft,1:nco2,1:12)
      REAL       dphen(1:365,1:npft)
      LOGICAL    present(1:npft),tree(1:npft)

c     LOCAL VARIABLES
      INTEGER    d,pft,m
      REAL       delta_phen
      REAL       delta_phen_total
      REAL       fuel_1hr_inc_total
      REAL       dfuel_leaf_background(1:365)


      dfuel_leaf(:,:)=0.0

      DO pft=1,npft
        IF (present(pft)) THEN                  !.AND.tree(pft)==.TRUE.
          delta_phen=dphen(1,pft)-dphen(365,pft)	!Doug 03/09, need to refer to previous years dphen?

          fuel_1hr_inc_total=SUM(fuel_1hr_leaf_inc_pos(pft,1,:))/
     *      (1+sum(dphen(:,pft))*(1-exp(leaf_decayc(pft))))

          dfuel_leaf_background=0.0


c          DO m=1,12
c            IF (fuel_1hr_inc(pft,1,m)>0) THEN
c
c              fuel_1hr_inc_total=fuel_1hr_inc_total+
c     *          fuel_1hr_inc(pft,1,m)
c
c            END IF
c          END DO


	  delta_phen_total=0.0

          DO d=1,365
            IF (d>1)  delta_phen=dphen(d,pft)-dphen(d-1,pft)

            IF (dphen(d,pft)>0) THEN

              dfuel_leaf_background(d)=dphen(d,pft)*
     *          fuel_1hr_inc_total*(1-exp(leaf_decayc(pft)))

            END IF

            IF (delta_phen>0) THEN

              delta_phen_total=delta_phen_total+delta_phen

              dfuel_leaf(d,pft)=delta_phen*fuel_1hr_inc_total

            END IF
          END DO

          IF (delta_phen_total>0) THEN

            dfuel_leaf(:,pft)=dfuel_leaf(:,pft)/
     *        delta_phen_total+dfuel_leaf_background(:)

            fuel_1hr_leaf_inc_pos(pft,1,:)=0.0
          END IF

C            DO m=1,12
C              IF (fuel_1hr_inc(pft,1,m)>0) THEN
C                 fuel_1hr_inc(pft,1,m)=0
C              END IF
C            END DO

        END IF
      END DO

      RETURN
      END    !fuel_1hr_redist

c----------------------------------------------------------------------
c     FUNCTION FIRE DANGER INDEX
c     Calculation of the fire danger index using the Nesterov Index
c


!      subroutine fire_danger_index(d_fdi,dlm,dlm_lg,dtemp_min,dtemp_max,
!     *    dprec,d,moistfactor,fuel_1hr_total,fuel_10hr_total,
!     *    fuel_100hr_total,dead_fuel,char_moistfactor,
!     *    ratio_dead_fuel,
!     *    ratio_live_fuel,dlm_1hr,dlm_10hr, dlm_100hr, dlm_1000hr,year,
!     *    ni_acc)

!      implicit none

!      real dlm(1:365),d_fdi(1:365),dlm_lg(1:365)
!      real dtemp_min(1:365),dtemp_max(1:365),dprec(1:365)
!      real moistfactor
!      real fuel_1hr_total,fuel_10hr_total,fuel_100hr_total
!      real dead_fuel
!      integer d,year
!      real char_moistfactor, ratio_dead_fuel,ratio_live_fuel
!      real   dlm_1hr(1:365), dlm_10hr(1:365)  
!      real   dlm_100hr(1:365), dlm_1000hr (1:365)
!      real ni_acc

!      real alpha      !coefficient of fire risk function per fuel class!!
!        parameter (alpha=0.0015)
!      real alpha_1hr, alpha_10hr,alpha_100hr
!      real alpha_livegrass,alpha_1000hr
!        parameter (alpha_1hr=0.001,alpha_10hr=0.00005424)	
!        parameter (alpha_100hr=0.00001485, alpha_1000hr = 0.000001) 
c        parameter (alpha_livegrass=0.0005) 
c     Allan: Alpha values for livegrass and 1000hr dead fuels are particularly
c            subjective 
c     KIRSTEN:i.e. values made up, revisit alpha_livegrass with livegrass
c             moisture driven upper soil moisture

c     LOCAL VARIABLES
c      real dw1(1:365)
!      real d_NI,alpha_fuel,char_alpha_fuel
!      real fdi_spread
!      real wk(1:365),sum_ni_acc
!      integer ii

c     initialise

c     dw1(:)=0.0
!      d_NI=0.0
!      alpha_fuel=0.0
!      char_alpha_fuel=0.0
!      fdi_spread=0.0
!      wk(:)=0.0
!      sum_ni_acc=0.0
!      ii=0

c     calculate Nesterov Index   equ. (2)
!        if (dprec(d).le.3.0.and.(dtemp_min(d)-4.0).ge.0.0) then
c         d_NI=dtemp_max(d)*(dtemp_max(d)-dtemp_min(d)-4.0)     !!! 25/10/07 Yan
!          d_NI=dtemp_max(d)*(dtemp_max(d)-(dtemp_min(d)-4.0))   !!! dtemp_dew=(dtemp_min(d)-4.0)
!          ni_acc=ni_acc+d_NI
!        else
!          d_NI=0.0
!          ni_acc=0.0
!        endif
        
c        if (d_NI.gt.0.0) then
c          ni_acc=ni_acc+d_NI
c        else
c          ni_acc=0.0
c        endif

c   litter moisture index, weighted per dead fuel class and
c   fuel load per dead fuel class; version with 3 alpha now confirmed

c   backcalculation of alpha_livegrass from livegrass moisture
!        if (ni_acc.gt.0.0) then
!          if (dlm_lg(d).gt.0.0) then
!           alpha_livegrass =(log(dlm_lg(d))/ni_acc)*(-1.0)
!         else
!           alpha_livegrass = 0.0
!          endif
!        endif

c    Allan
!        char_alpha_fuel=alpha_fuel * ratio_dead_fuel
!     *          + alpha_livegrass * ratio_live_fuel
!       else
c         alpha_fuel=0.00001
!         char_alpha_fuel=0.00001
!       endif

c        dlm(d)=exp(-alpha_fuel*ni_acc)
!        dlm(d)=exp(-char_alpha_fuel*ni_acc)  ! used for ROS calcs. 1000hr therefore ignored.
       
!       dlm_1hr(d) = exp(-alpha_1hr * ni_acc)
!       dlm_10hr(d) = exp(-alpha_10hr * ni_acc)
!       dlm_100hr(d) = exp(-alpha_100hr * ni_acc)
!       dlm_1000hr(d) = exp(-alpha_1000hr * ni_acc)
      
       
c probability of fire spread
*        if (dlm(d).le.moistfactor) then
*           fdi_spread=1.0-(dlm(d)/moistfactor)
*            fdi_spread=max(0.0,(1.0-2.59*(dlm(d)/moistfactor)+
*     *             5.11*(dlm(d)/moistfactor)**2.0-
*     *             3.52*(dlm(d)/moistfactor)**3.0))
*        else
*           fdi_spread=0.0
*        endif

c     calculate Fire Danger Index equ. (9)
!        if (d_NI.le.0.0) then
!           d_fdi(d)=0.0
!        else
c        d_fdi(d)=max(0.0,
c     *           (1.0-((1.0/moistfactor)*(exp(-alpha_fuel*ni_acc)))))
c
*     *           (1.0-((1.0/moistfactor)*dlm(d))))
!        d_fdi(d)=max(0.0,(1.0-((1.0/char_moistfactor)*
!     *                   (exp(-char_alpha_fuel*ni_acc)))))

!        endif

!      return
!      end	!FIRE DANGER INDEX


c----------------------------------------------------------------------
c     FUNCTION FIRE DANGER INDEX1
c     Calculation of the Nesterov Index and the dead fuel moisture.
c

      subroutine fire_danger_index1(dlm,dlm_lg,dtemp_min,dtemp_max,
     *    dprec,d,m,fuel_1hr_total,fuel_10hr_total,
     *    fuel_100hr_total,dead_fuel,
     *    ratio_dead_fuel,
     *    ratio_live_fuel,dlm_1hr,dlm_10hr, dlm_100hr, dlm_1000hr,year,
     *    ni_acc,char_alpha_fuel,d_NI,lon,lat,
     *    dlm_1hr_old,dlm_10hr_old,dlm_100hr_old,dlm_1000hr_old,
     *    mlm,dpet,aprec)

      implicit none

      real dlm(1:365),d_fdi(1:365),dlm_lg(1:365),mlm(1:12)
      real dtemp_min(1:365),dtemp_max(1:365),dprec(1:365)
      real moistfactor
      real fuel_1hr_total,fuel_10hr_total,fuel_100hr_total
      real dead_fuel
      integer d,year,m
      real char_moistfactor, ratio_dead_fuel,ratio_live_fuel
      real   dlm_1hr(1:365), dlm_10hr(1:365)  
      real   dlm_100hr(1:365), dlm_1000hr (1:365)
      real ni_acc
      
      real alpha      !coefficient of fire risk function per fuel class!!
        parameter (alpha=0.0015)
      real alpha_1hr, alpha_10hr,alpha_100hr
      real alpha_livegrass,alpha_1000hr
        parameter (alpha_1hr=24,alpha_10hr=2.4)	!Doug 08/12: fiddleing with the alphas
        parameter (alpha_100hr=0.24,
     *		alpha_1000hr = 0.024) 

c     Allan: Alpha values for livegrass and 1000hr dead fuels are particularly
c            subjective 
c     KIRSTEN:i.e. values made up, revisit alpha_livegrass with livegrass
c             moisture driven upper soil moisture

c     LOCAL VARIABLES
c      real dw1(1:365)
      real d_NI,alpha_fuel,char_alpha_fuel
      real fdi_spread
      real wk(1:365),sum_ni_acc
      integer ii
      real lon,lat
      real month_length(1:12)
      data (month_length(m),m=1,12)
     *               /31.0,28.0,31.0,30.0,31.0,30.0,31.0,31.0,
     *                30.0,31.0,30.0,31.0/
      
c    Variables for calculating new fuel moisture uisng RH
      REAL rhumid, emc
      REAL rhumid_c1,rhumid_c2 
        PARAMETER (rhumid_c1=17.271,rhumid_c2= 237.7)
      REAL dlm_1hr_old,dlm_10hr_old,dlm_100hr_old,dlm_1000hr_old
      REAL dpet(1:365)
      REAL aprec
      REAL temp_dew
      REAL EF
      REAL dew_c1,dew_c2,dew_c3,dew_c4,dew_c5,dew_c6,dew_c7 
        PARAMETER (dew_c1=-0.127,dew_c2= 1.121,dew_c3= 1.003,
     *      dew_c4=-1.444,dew_c5= 12.312,
     *      dew_c6=-32.766,dew_c7= 0.0006)
c     initialise

c     dw1(:)=0.0
      d_NI=0.0
      alpha_fuel=0.0
      char_alpha_fuel=0.0
      fdi_spread=0.0
      wk(:)=0.0
      sum_ni_acc=0.0
      ii=0

c     calculate Nesterov Index   equ. (2)
        if (dprec(d).le.3.0.and.(dtemp_min(d)-4.0).ge.0.0) then
c         d_NI=dtemp_max(d)*(dtemp_max(d)-dtemp_min(d)-4.0)     !!! 25/10/07 Yan
          d_NI=dtemp_max(d)*(dtemp_max(d)-(dtemp_min(d)-4.0))   !!! dtemp_dew=(dtemp_min(d)-4.0)
          ni_acc=ni_acc+d_NI
        else
          d_NI=0.0
          ni_acc=0.0
        endif
        
c        if (d_NI.gt.0.0) then
c          ni_acc=ni_acc+d_NI
c        else
c          ni_acc=0.0
c        endif

c   litter moisture index, weighted per dead fuel class and
c   fuel load per dead fuel class; version with 3 alpha now confirmed

*        dlm(d)=exp(-alpha*ni_acc)
       if (dead_fuel.gt.0.0) then
        alpha_fuel=(alpha_1hr*fuel_1hr_total+
     *              alpha_10hr*fuel_10hr_total+
     *              alpha_100hr*fuel_100hr_total)/dead_fuel

c   backcalculation of alpha_livegrass from livegrass moisture
        if (ni_acc.gt.0.0) then
          if (dlm_lg(d).gt.0.0) then
           alpha_livegrass =(log(dlm_lg(d))/ni_acc)*(-1.0)
             if(alpha_livegrass.lt.0.0)alpha_livegrass=0.0
         else
           alpha_livegrass = 0.0
          endif
        endif

c    Allan
        char_alpha_fuel=alpha_fuel * ratio_dead_fuel
     *          + alpha_livegrass * ratio_live_fuel
       else
c         alpha_fuel=0.00001
         char_alpha_fuel=0.00001
       endif
        
C Doug 06/13: dlm_xhr calculated using NI replaced with RH 
       IF (dprec(d).le.3.0.and.(dtemp_min(d)-4.0).ge.0.0) THEN
            IF (aprec.EQ.0) THEN
                emc=0.0
            ELSE
                EF=dpet(d)/aprec
                IF (dtemp_max(d)<dtemp_min(d)) THEN
                    dtemp_max(d)=dtemp_min(d)
                ENDIF
            temp_dew=(dtemp_min(d)+273.15)*(dew_c1+dew_c2*
     *          (dew_c3+dew_c4*EF+dew_c5*(EF**2)+
     *           dew_c6*(EF**3))+
     *           dew_c7*(dtemp_max(d)-dtemp_min(d)))
            temp_dew=temp_dew-273.15
            IF (temp_dew<-50) THEN
                rhumid=1
            ELSE
            
                rhumid=100*(exp((rhumid_c1*(temp_dew))/
     *              (rhumid_c2+(temp_dew))))
            
                rhumid=rhumid/exp((rhumid_c1*
     *              dtemp_max(d))/
     *              (rhumid_c2+dtemp_max(d)))
            END IF
            IF (rhumid<0) rhumid=0
            IF (rhumid>100) rhumid=100
            emc=0.942*(rhumid**0.679)+0.000499*exp(0.1*rhumid)+
     *          0.18*(21.1-dtemp_max(d))*
     *          (1-exp(-0.115*rhumid))
            END IF
        ELSE
            emc=100
        ENDIF
        IF (ISNAN(dlm_1hr_old)) THEN
            dlm_1hr(d) = emc/100
        ELSE
            dlm_1hr(d) = emc/100+
     *        (dlm_1hr_old-emc/100)*exp(-alpha_1hr)
        ENDIF
         IF (ISNAN(dlm_10hr_old)) THEN
            dlm_10hr(d) = emc/100
        ELSE 
            dlm_10hr(d) = emc/100+
     *        (dlm_10hr_old-emc/100)*exp(-alpha_10hr)
        ENDIF
        IF (ISNAN(dlm_100hr_old)) THEN
            dlm_100hr(d) = emc/100
        ELSE
            dlm_100hr(d) = emc/100+
     *        (dlm_100hr_old-emc/100)*exp(-alpha_100hr)
        ENDIF
        IF (ISNAN(dlm_1000hr_old)) THEN
            dlm_1000hr(d) = emc/100
        ELSE
            dlm_1000hr(d) = emc/100+
     *        (dlm_1000hr_old-emc/100)*exp(-alpha_1000hr)
        ENDIF
      
        dlm_1hr_old=dlm_1hr(d)
        dlm_10hr_old=dlm_10hr(d)
        dlm_100hr_old=dlm_100hr(d)
        dlm_1000hr_old=dlm_1000hr(d)
      

c Leilei 11/08: dlm= the sum of dlm for each different fuel load/
C Doug 09/13: dlm is just live grass moisture if there in no
C   dead fuel. Not correcting for this maked dlm an NaN if no fuel is left.
       IF (dead_fuel==0) THEN
            dlm(d)=dlm_lg(d)
       ELSE
            dlm(d)=((dlm_1hr(d)*fuel_1hr_total+dlm_10hr(d)*
     *         fuel_10hr_total+
     *         dlm_100hr(d)*fuel_100hr_total)/dead_fuel)
     *         *ratio_dead_fuel+			!Doug 12/08 *1
     *         dlm_lg(d)*ratio_live_fuel		!Doug 12/08 remove
       END IF
        mlm(m)=mlm(m)+dlm(d)/month_length(m)

        !if (dlm(d)>1.05) then
        !    WRITE(11,*), "litter moistures:"
        !    WRITE(11,*), dlm_1hr(d) !StopC
        !    WRITE(11,*), dlm_10hr(d)
        !    WRITE(11,*), dlm_100hr(d)
        !    WRITE(11,*), "grass moisture:"
        !    WRITE(11,*), dlm_lg(d)
        !    WRITE(11,*), "litter amounts:"
        !    WRITE(11,*), fuel_1hr_total
        !    WRITE(11,*), fuel_10hr_total
        !    WRITE(11,*), fuel_100hr_total
        !    WRITE(11,*), "litter-to-grass ratio:"
        !    WRITE(11,*), ratio_dead_fuel
        !    WRITE(11,*), ratio_live_fuel
        !end if
        
        if (isnan(mlm(m))) then
            WRITE(12,*), "dlm:"
            WRITE(12,*), dlm(d)
            WRITE(12,*), "dlms:"
            WRITE(12,*), dlm_1hr(d)
            WRITE(12,*), dlm_10hr(d)
            WRITE(12,*), dlm_100hr(d)
            WRITE(12,*), "dlms_old:"
            WRITE(12,*), dlm_1hr_old
            WRITE(12,*), dlm_10hr_old
            WRITE(12,*), dlm_100hr_old
            WRITE(12,*), "fuel loads"
            WRITE(12,*), fuel_1hr_total
            WRITE(12,*), fuel_10hr_total
            WRITE(12,*), fuel_100hr_total
            WRITE(12,*), dead_fuel
            WRITE(12,*), ratio_dead_fuel
            WRITE(12,*), ratio_live_fuel
            WRITE(12,*), "emc"
            WRITE(12,*), emc
            WRITE(12,*), "rhumid"
            WRITE(12,*), rhumid
            WRITE(12,*), "dtemp:"
            WRITE(12,*), dtemp_max(d)
            WRITE(12,*), dtemp_min(d)
            WRITE(12,*), temp_dew
            WRITE(12,*), "EF"
            WRITE(12,*), EF
            WRITE(12,*), dpet(d)
            WRITE(12,*), aprec
        end if
            
            
c probability of fire spread
*        if (dlm(d).le.moistfactor) then
*           fdi_spread=1.0-(dlm(d)/moistfactor)
*            fdi_spread=max(0.0,(1.0-2.59*(dlm(d)/moistfactor)+
*     *             5.11*(dlm(d)/moistfactor)**2.0-
*     *             3.52*(dlm(d)/moistfactor)**3.0))
*        else
*           fdi_spread=0.0
*        endif

        return
        end


*****************************************************************************************

*******************************************************************************************

       subroutine fire_danger_index(char_moistfactor,char_alpha_fuel,
     *                              ni_acc,d_NI,d_fdi,d)

       implicit none

       real d_fdi(1:365),char_moistfactor,char_alpha_fuel
       real ni_acc,d_NI
       integer d
c     calculate Fire Danger Index equ. (9)
        if (d_NI.le.0.0) then
           d_fdi(d)=0.0
        else
c        d_fdi(d)=max(0.0,
c     *           (1.0-((1.0/moistfactor)*(exp(-alpha_fuel*ni_acc)))))
c
*     *           (1.0-((1.0/moistfactor)*dlm(d))))

        d_fdi(d)=max(0.0,(1.0-((1.0/char_moistfactor)*
     *                   (exp(-char_alpha_fuel*ni_acc)))))

        endif
      return
      end


c-------------------------------------------------------------------------------
c     FUNCTION HUMAN IGNITION
c     Calculation of the numbers of ignition caused by humans, which is
c     determined by population density, see equ. (13), where the value of
c     a(N_d) and k(pop_den) are equal to their mathematical expectation

      real function human_ign(popden,a_nd,year)

      implicit none

      integer year
      real popden,a_nd

      real a_nd_const
      parameter (a_nd_const=0.251)


        integer nran
        real ran,a_ndd



c       if (popden.lt.1.0) popden=1.0
c       if (a_nd.lt.0.11) a_nd=0.165
c      popden=7.5
c       a_nd=1.59
c      human_ign=gamma*popden1**0.43
c       human_ign=6.8*(exp(-0.5*(popden1**0.5)))*a_nd_const*popden1

       nran=5249347
       nran=mod(23*nran,10000001)
       ran=float(nran)*1.E-7
       a_ndd=(-0.005)*log(ran) ! 0.22 = 1/lamda


       if(popden.gt.0.0)then
          human_ign=6.8*(popden**(-0.57))*popden*a_ndd
       else
          human_ign=0.0
       end if
!       human_ign=30.0*(exp(-0.5*(popden**0.5)))*a_nd*popden
      human_ign=0.0
      end
	

c-------------------------------------------------------------------------------
c     FUNCTION LIGHTNING
c     Calculation of lightning strikes. Constant at the moment, in the future
c     depending on latitude and other factors
c
c     Kirsten: function not used with potential lightning ignitions as an input variable
      real function lightn_ign(lat)

      implicit none

      real lat

*      lightn=0.02
c     for case of Brandenburg
       lightn_ign=0.2

      end


c-------------------------------------------------------------------------------
c     subroutine Thermophysical Propertities of the Fuel Array, i.e., bet and q_ig 

      subroutine fuel_property(dens_fuel_ave,sigma,
     *            dlm,char_dens_fuel_ave, 
     * char_sigma,char_dens_fuelbed,d,bet,q_ig)


      implicit none

      integer npft,npftpar
        parameter (npft=13,npftpar=59)
        
      real dens_fuel_ave,sigma,H
      real dlm(1:365)
       real char_dens_fuel_ave,char_sigma
       integer d

c     emissitivity of flames unitless
*      real emiss
*        parameter(emiss=0.3)

      real beta,q_ig
      real bet,beta_op

      real part_dens
        parameter (part_dens=513.0)! this is oven-dry particle density (kg/m3)
		
       real char_dens_fuelbed		

c    initialise
      beta=0.0
      q_ig=0.0
      bet=0.0
      beta_op=0.0

c     start of function

c        bet=dens_fuel_ave*0.200395*(sigma**(-0.8189))/513.0     !=beta/beta_op

       beta = char_dens_fuelbed/part_dens 


c       beta_op = 0.200395*(char_sigma**(-0.8189))
       beta_op = 0.200395*(char_sigma**(-0.8189))

       bet = beta / beta_op

c     heat of pre-ignition
c     ACHTUNG: check influence of litter moisture!!!
       q_ig=581.0+2594.0*dlm(d)

         end
c------------------------------------------------------------------------------------------------------------

c     subroutine RATE OF forward SPREAD

      subroutine rate_of_spread(U_front,base_wind,dens_fuel_ave,sigma,
     *            dlm,d,net_fuel,moistfactor,H, char_dens_fuel_ave, 
     * char_sigma, char_net_fuel, char_moistfactor,gamma,
     *   char_dens_fuelbed)
      implicit none

      integer npft,npftpar
        parameter (npft=13,npftpar=59)
!        parameter (npft=13,npftpar=50)
        
      real U_front,base_wind
      real dens_fuel_ave,sigma,H
      real dlm(1:365)
      integer d
      real net_fuel,moistfactor
       real char_dens_fuel_ave,char_sigma
       real char_net_fuel 
       real char_moistfactor
       real gamma

c     emissitivity of flames unitless
*      real emiss
*        parameter(emiss=0.3)

c     Rothermal fire spread

c      real mw1(1:12),dw1(1:365)
      real dummy,dummy2
      real beta,ir,xi,eps,q_ig,phi_wind
      real gamma_aptr,moist_damp,mw_weight
      real gamma_max,bet,beta_op,a,c,b,e

c     mineral dampening coefficient
      real MINER_DAMP
      real pi
        parameter (pi=3.14159265)
      real part_dens
        parameter (part_dens=513.0)! this is oven-dry particle density (kg/m3)
		
       real char_dens_fuelbed		

c    initialise
c      mw1(:)=0.0
c      dw1(:)=0.0
      dummy=0.0
      dummy2=0.0
      beta=0.0
      ir=0.0
      xi=0.0
      eps=0.0
      q_ig=0.0
      phi_wind=0.0
      gamma_aptr=0.0
      moist_damp=0.0
      mw_weight=0.0
      gamma_max=0.0
      bet=0.0
      beta_op=0.0
      a=0.0
      b=0.0
      c=0.0
      e=0.0

c     start of function

c        bet=dens_fuel_ave*0.200395*(sigma**(-0.8189))/513.0     !=beta/beta_op

       beta = char_dens_fuelbed/part_dens 
c       beta_op = 0.200395*(char_sigma**(-0.8189))
       beta_op = 0.200395*(char_sigma**(-0.8189))

       bet = beta / beta_op

c     heat of pre-ignition
c     ACHTUNG: check influence of litter moisture!!!
       q_ig=581.0+2594.0*dlm(d)
c     effective heating number

!       eps=exp(-4.528/sigma)
c        eps=exp(-4.528/char_sigma)
        eps=exp(-4.528/char_sigma)

c     influence of wind speed

c           b = 0.15988 * (char_sigma**0.54)
c           c = 7.47 * (exp(-0.8711 * (char_sigma**0.55)))
c           e = 0.715 * (exp(-0.01094 * char_sigma))

           b = 0.15988 * (char_sigma**0.54)
           c = 7.47 * (exp(-0.8711 * (char_sigma**0.55)))
           e = 0.715 * (exp(-0.01094 * char_sigma))


       phi_wind=c*(base_wind**b)*(bet**(-e))
c     propagating flux


c       if (char_sigma.le.0.00001) then
       if (char_sigma.le.0.00001) then

        xi=0.0
       else
c        xi=(exp((0.792+3.7597*(sigma**0.5))*
c     *     ((dens_fuel_ave/513.0)+0.1)))/
c     *     (192.0+7.9095*sigma)
	 
	        xi = (exp((0.792 + 3.7597 * (char_sigma**0.5)) *
     *     (beta + 0.1))) / (192 + 7.9095 * char_sigma) 

       endif

c    reaction intensity


c           a = 8.9033 * (char_sigma**(-0.7913))
           a =1.0/(6.71305*(char_sigma**0.1)-7.27)
c           if (char_sigma.le.0.00001) then
           if (char_sigma.le.0.00001) then

            dummy=0.0
           else
            dummy=exp(a*(1.0-bet))
           endif

c           gamma_max = 1.0 / (0.0591 + 2.926 * (char_sigma**(-1.5)))
           gamma_max = 1.0 / (0.0591 + 2.9401 * (char_sigma**(-1.5)))


*        gamma=gamma_max*(bet**a)*exp(a*(1.0-bet))

         gamma_aptr=gamma_max*(bet**a)*dummy

c        if (moistfactor.gt.0.0) then
        if (char_moistfactor.gt.0.0) then

c          mw_weight=dlm(d)/moistfactor
         mw_weight=dlm(d)/char_moistfactor
        else
          mw_weight=0.0
        endif

        moist_damp=max(0.0,(1.0-(2.59*mw_weight)+
     *     (5.11*(mw_weight**2.0))-(3.52*(mw_weight**3.0))))


        MINER_DAMP=0.174*(0.055**(-0.19))

c        ir=gamma*net_fuel*H*moist_damp*MINER_DAMP
         ir=gamma_aptr * char_net_fuel * H * moist_damp * MINER_DAMP

c        for use in postfire mortality (tau_r)
         gamma=gamma_aptr*moist_damp*MINER_DAMP
      
c    reaction intensity end

c        if(dens_fuel_ave.le.0.0.or.eps.le.0.0.or.q_ig.le.0.0) then
         if ((char_dens_fuel_ave.le.0.0).or.
     *      (eps.le.0.0).or.(q_ig.le.0.0)) then
          U_front=0.0
        else
c          U_front=(ir*xi*(1.0+phi_wind))/(dens_fuel_ave*eps*q_ig)
          U_front=(ir * xi * (1.0 + phi_wind)) /
     *            (char_dens_fuelbed * eps * q_ig)

        endif

      return
      end
c-------------------------------------------------------------------------------
c Doug 02/13: Calculates tree mortality and new Bark thickness distribution from cambiol damage
      SUBROUTINE BT_change(DBH,BTparam1,BTparam2,
     *   tau_l,pm_tau_class)
	 
       IMPLICIT NONE
	   
C      Inputs
       REAL DBH
       REAL tau_l
	   
C      INPUTS/OUTPUTS
       REAL BTparam1(1:3),BTparam2(1:3)
       REAL pm_tau_class
	   
C      local vaiables
       REAL W,X,Y,Z
       REAL phi,chi,psi,omg
       REAL s1,s2
       REAL l1,l2,l3,l4,l5,l6,l7,l8
       REAL a,b,AA,BB
       REAL d1,d2
       REAL BTmean,BTmode,BTmode0
       REAL btmode_frac
	   
       pm_tau_class=0.0
	   
       AA=1.125
       BB=0.563*tau_l/2.9
	   
       a=sqrt(abs(BB)/AA);
       b=sqrt(abs(BB)/(AA-1));
	   
       IF (isnan(DBH)) DBH=0.0
       
	   
       s1=BTparam1(1)+DBH*BTparam2(1)
       s2=BTparam1(3)+DBH*BTparam2(3)
	   
       BTmode=BTparam1(2)+DBH*BTparam2(2)
       BTmode0=BTmode
       l1=0.0
       l2=0.0
       l3=0.0
       l4=0.0
       l5=0.0
       l6=0.0
       l7=0.0
       l8=0.0
	  
	   
       IF (a>=s2) THEN
	     pm_tau_class=1.0
         RETURN
       ELSE IF (b<=s1) THEN
	     RETURN
       ELSE IF (BTmode<=a.AND.b>s2) THEN
	     l3=a
	     l4=s2
       ELSE IF (BTmode<=a.AND.b<=s2) THEN
	     l3=a
	     l4=b
	     l7=b
	     l8=s2
       ELSE IF (s1<=a.AND.b>s2) THEN
	     l1=a
	     l2=BTmode
	     l3=BTmode
	     l4=s2
       ELSE IF (s1<=a.AND.b>BTmode) THEN
         l1=a
         l2=BTmode
         l3=BTmode
         l4=b
         l7=b
         l8=s2
       ELSE IF (s1<=a.AND.b<=BTmode) THEN
	     l1=a
	     l2=b
	     l5=b
	     l6=BTmode
	     l7=BTmode
	     l8=s2
       ELSE IF (a<s1.AND.b>s2) THEN
	     l1=s1
	     l2=BTmode
	     l3=BTmode
	     l4=s2
       ELSE IF (a<s2.AND.b>BTmode) THEN
	     l1=s1
	     l2=BTmode
	     l3=BTmode
	     l4=b
	     l7=b
	     l8=s2
       ELSE IF (a<s2.AND.b>s1) THEN
	     l1=s1
	     l2=b
	     l5=b
	     l6=BTmode
	     l7=BTmode
	     l8=s2
       END IF
	   
      W=0.0
      X=0.0
      Y=0.0
      Z=0.0
	  
      phi=0.0
      chi=0.0
      psi=0.0
      omg=0.0
	  
      d1=(s2-s1)*(BTmode-s1)
      d2=(s2-s1)*(s2-BTmode)	   
	   
       IF (l1>=0.0 .AND.l2>0.0) THEN
         W=-(2*AA*((l1**3)-(l2**3))+3*AA*s1*((l2**2)-(l1**2))+
     *    6*BB*s1*(log(l1)-log(l2))+6*BB*(l2-l1))/(3*d1)
     
       IF (W<0.0001 .AND. W >-0.0001) W=0
	 
         phi=-(AA*((l1**2)-(l2**2))+2*AA*s1*(l2-l1)+
     *     2*BB*(log(l2)-log(l1))+2*BB*s1*((1/l2)-(1/l1)))/d1
     
       IF (phi<0.0001 .AND. phi >-0.0001) phi=0
       
       END IF
	   
       IF (l3>=0.0 .AND.l4>0.0) THEN
         X=(2*AA*((l3**3)-(l4**3))+(3*AA*s2*((l4**2)-(l3**2))+
     *     6*BB*s2*(log(l3)-log(l4))+6*BB*(l4-l3)))/(3*d2)
     
         IF (X<0.001 .AND. X >-0.001) X=0.0
	 
	     chi=(AA*((l3**2)-(l4**2))+2*AA*s2*(l4-l3)+
     *     2*BB*(log(l4)-log(l3))+2*BB*s2*((1/l4)-(1/l3)))/d2
     
         IF (chi<0.001 .AND. chi >-0.001) chi=0.0
     
       END IF
	   
       IF (l5>=0.0 .AND. l6>0.0) THEN
	     Y=-(2*((l5**3)-(l6**3))+3*s1*((l5**2)-(l6**2)))/(3*d1)        
         
		 
	     psi=-((l5**2)-(l6**2)+2*s1*(l6-l5))/d1
         
         IF (psi<0.0001 .AND. psi >-0.0001) psi=0.0;
         IF (Y<0.0001 .AND. Y >-0.0001) Y=0.0; psi=0.0
         
       END IF
	   
       IF (l7>=0.0 .AND. l8>0.0) THEN
         Z=(2*((l7**3)-(l8**3))+3*s2*((l8**2)-(l7**2)))/(3*d2)
         
         IF (Z<0.0001 .AND. Z >-0.0001) Z=0.0
	 
	     omg=((l7**2)-(l8**2)+2*s2*(l8-l7))/d2
         
         IF (omg<0.0001 .AND. omg >-0.0001) omg=0.0
         
       END IF
	   
       pm_tau_class=1-(phi+chi+psi+omg)
	   
       IF (pm_tau_class<0) pm_tau_class=0
       IF (pm_tau_class>1) pm_tau_class=1
	
        IF ((phi+chi+psi+omg)<=0) THEN
            BTmean=0
            BTmode_frac=0
        ELSE        
             BTmean=(W+X+Y+Z)/(phi+chi+psi+omg)	   
             BTmode=3*BTmean-s1-s2
             
             IF (BTmode>s2) BTmode=s2*0.99
             BTmode_frac=(BTmode-BTmode0)/(5*(s2-BTmode0))
        END IF
		 
       
       
       IF (BTmode<0 .OR. BTmode_frac>1) THEN
	   !  WRITE(10,*),"+++++++"
       !  WRITE(10,*), "BT change fire error"
       !  
       !  WRITE(10,*), "error: BT medium is lower or higher than"
       !  WRITE(10,*),  "lowest/highest possible limit in fire"
       !  WRITE(10,*), "W, X, Y, Z: ", W, X, Y, Z
       !  WRITE(10,*), "phi, chi, psi, omg", phi, chi, psi, omg
       !  WRITE(10,*), "s1, s2", s1, s2
       !  WRITE(10,*), "BTmode0: ", BTmode0
       !  WRITE(10,*), "BTmean: ", BTmean
       !  WRITE(10,*), "BTmode: ", BTmode
         BTmode_frac=0.0        
	      
       END IF

       BTparam1(2)=BTparam1(2)+BTmode_frac*(BTparam1(3)-BTparam1(2))
       BTparam2(2)=BTparam2(2)+BTmode_frac*(BTparam2(3)-BTparam2(2))
       
       RETURN
	   
       end
            
c-------------------------------------------------------------------------------
c    subroutine for calculation of fuel consumption in the area affected by fire

      subroutine fuel_consumption(npft,present,fuel_consum,fire_frac,
     *  fuel_1hr_leaf_left,fuel_1hr_wood_left,
     *  fuel_10hr_left,fuel_100hr_left,		!Doug 03/09: fuel_xhr_0 --> fuel_xhr_left (annual --> daily)
     *   fuel_1000hr_left,livegrass_left,pot_fc_lg_0,
     *  tree,fuel_1hr_total,moistfactor,dlm,d,MINER_TOT,
     *  fc_1hr_leaf,fc_1hr_wood,fc_lg,
     *  fc_10hr,fc_100hr,fc_1000hr,cf,char_moistfactor,
     *   dlm_lg, dlm_1hr, dlm_10hr, dlm_100hr, dlm_1000hr,
     *   moistfactor_livegrass, moistfactor_1hr,moistfactor_10hr,
     *   moistfactor_100hr, moistfactor_1000hr)

      implicit none

c     PARAMETERS
      integer nco2,nc              ! nc is added for nco2
         parameter (nco2=3)

      integer npft,d
      logical present(1:npft)
      real fuel_consum
      real fuel_consum_all
      real d_fuel_all
      real fire_frac(1:365)
      real pot_fc_lg(1:npft,1:nco2)
      REAL fuel_1hr_leaf(1:npft,1:nco2)
      REAL fuel_1hr_wood(1:npft,1:nco2)
      real fuel_10hr(1:npft,1:nco2)
      real fuel_100hr(1:npft,1:nco2),fuel_1000hr(1:npft,1:nco2)
      real fuel_1hr_total,livegrass
      logical tree(1:npft)
      real moistfactor
      real dlm(1:365)
      real MINER_TOT
      REAL fc_1hr_leaf(1:npft,1:nco2)
      REAL fc_1hr_wood(1:npft,1:nco2)
      REAL fc_lg(1:npft,1:nco2)
      real fc_10hr(1:npft,1:nco2)
      real fc_100hr(1:npft,1:nco2),fc_1000hr(1:npft,1:nco2)
      real cf(1:npft)
      real char_moistfactor
      real  dlm_lg (1:365), dlm_1hr(1:365), dlm_10hr(1:365)  
      real  dlm_100hr(1:365), dlm_1000hr (1:365)
      real moistfactor_livegrass, moistfactor_1hr,moistfactor_10hr
      real moistfactor_100hr, moistfactor_1000hr

c     local variables
      integer pft
      real moist_lg_d1hr,moist_1hr,moist_10_100hr
      real net_moist_lg, net_moist_1hr, net_moist_10hr
      real net_moist_100hr, net_moist_1000hr
      real pot_fc_1hr_total,mw_weight
      real cf_lg, cf_1hr,cf_10hr,cf_100hr,cf_1000hr
c      real dw1(1:365)
      real tau_l(1:npft)

      REAL fuel_1hr_leaf_left(1:npft,1:nco2)
      REAL fuel_1hr_wood_left(1:npft,1:nco2)
      REAL fuel_10hr_left(1:npft,1:nco2)		!Doug 03/09: fuel_xhr_0 --> fuel_xhr_left (annual --> daily)
      real fuel_100hr_left(1:npft,1:nco2)
      real fuel_1000hr_left(1:npft,1:nco2)
      real livegrass_left,pot_fc_lg_0(1:npft,1:nco2)

c     initialise
      pft=0
      moist_lg_d1hr=0.0
      moist_1hr=0.0
      moist_10_100hr=0.0
      net_moist_lg=0.0
      net_moist_1hr=0.0
      net_moist_10hr=0.0
      net_moist_100hr=0.0
      net_moist_1000hr=0.0
      pot_fc_1hr_total=0.0
      mw_weight=0.0
      cf_lg=0.0
      cf_1hr=0.0
      cf_10hr=0.0
      cf_100hr=0.0
      cf_1000hr=0.0
c      dw1(:)=0.0
      tau_l(:)=0.0

      fuel_consum=0.0
      fuel_consum_all=0.0
      d_fuel_all=0.0

c     fuel consumption depending on fuel moisture
c     influence of livefuel on 1hr fuel moisture content
c     ACHTUNG: change with dlm_lg as a function of upper soil moisture
c     CHECK method

       moist_lg_d1hr=dlm_1hr(d)+(dlm_lg(d)*
     *                (livegrass_left/fuel_1hr_total))

      if(moistfactor.gt.0.0) then
        moist_1hr=moist_lg_d1hr/char_moistfactor
        moist_10_100hr=dlm(d)/char_moistfactor
      else
        moist_1hr=1.0
        mw_weight=0.0
        moist_10_100hr=1.0
      endif


      do pft=1,npft
       do nc=1,nco2
        fc_1hr_leaf(pft,nc)=0.0
        fc_1hr_wood(pft,nc)=0.0
        fc_lg(pft,nc)=0.0
        fc_10hr(pft,nc)=0.0
        fc_100hr(pft,nc)=0.0
        fc_1000hr(pft,nc)=0.0
       enddo
c        tau_l(pft)=0.0 !cycle through the dead fuel classes and livegrass


       if (present(pft)) then

c     1hr fuel consumption

      if(moist_1hr.le.0.18) then
         fc_1hr_leaf(pft,1)=1.0*(1.0-MINER_TOT)*
     *     (fuel_1hr_leaf_left(pft,1)/0.45)*fire_frac(d)
         fc_1hr_leaf(pft,2)=fuel_1hr_leaf_left(pft,2)
         fc_1hr_leaf(pft,3)=fuel_1hr_leaf_left(pft,3)
		 
         fc_1hr_wood(pft,1)=1.0*(1.0-MINER_TOT)*
     *     (fuel_1hr_wood_left(pft,1)/0.45)*fire_frac(d)
         fc_1hr_wood(pft,2)=fuel_1hr_wood_left(pft,2)
         fc_1hr_wood(pft,3)=fuel_1hr_wood_left(pft,3)
         cf_1hr=1.0  
      else
         if(moist_1hr.gt.0.18.and.moist_1hr.le.0.73) then
            fc_1hr_leaf(pft,1)=(1.10-0.62*moist_1hr)*(1.0-MINER_TOT)*
     *                   (fuel_1hr_leaf_left(pft,1)/0.45)*fire_frac(d)    !(fuel_1hr_0(pft,1)/0.45)*fire_frac(d)
            fc_1hr_leaf(pft,2)=fuel_1hr_leaf_left(pft,2)
            fc_1hr_leaf(pft,3)=fuel_1hr_leaf_left(pft,3)
			
            fc_1hr_wood(pft,1)=(1.10-0.62*moist_1hr)*(1.0-MINER_TOT)*
     *                   (fuel_1hr_wood_left(pft,1)/0.45)*fire_frac(d)    !(fuel_1hr_0(pft,1)/0.45)*fire_frac(d)
            fc_1hr_wood(pft,2)=fuel_1hr_wood_left(pft,2)
            fc_1hr_wood(pft,3)=fuel_1hr_wood_left(pft,3)
            cf_1hr=1.10-0.62*moist_1hr
         else
           if(moist_1hr.gt.0.73.and.moist_1hr.le.1.0) then
            fc_1hr_leaf(pft,1)=(2.45-2.45*moist_1hr)*(1.0-MINER_TOT)*
     *                   (fuel_1hr_leaf_left(pft,1)/0.45)*fire_frac(d)   !(fuel_1hr_0(pft,1)/0.45)*fire_frac(d)
            fc_1hr_leaf(pft,2)=fuel_1hr_leaf_left(pft,2)
            fc_1hr_leaf(pft,3)=fuel_1hr_leaf_left(pft,3)
			
            fc_1hr_wood(pft,1)=(2.45-2.45*moist_1hr)*(1.0-MINER_TOT)*
     *                   (fuel_1hr_wood_left(pft,1)/0.45)*fire_frac(d)   !(fuel_1hr_0(pft,1)/0.45)*fire_frac(d)
            fc_1hr_wood(pft,2)=fuel_1hr_wood_left(pft,2)
            fc_1hr_wood(pft,3)=fuel_1hr_wood_left(pft,3)
            cf_1hr=2.45-2.45*moist_1hr
           else
            fc_1hr_leaf(pft,1)=0.0
            fc_1hr_leaf(pft,2)=0.0
            fc_1hr_leaf(pft,3)=0.0
			
            fc_1hr_wood(pft,1)=0.0
            fc_1hr_wood(pft,2)=0.0
            fc_1hr_wood(pft,3)=0.0
            cf_1hr=0.0
           endif
         endif
      endif


c    time required for cambial kill, 1hr fuel [gBiomass/cm?²]
c    /1e4 to get from 1m2 to cm2, the units used by P&R(1986) 
c         tau_l(pft)=tau_l(pft)+((fuel_1hr(pft)/0.45/1e4)*
c     *                     (1.0-((1.0-cf_1hr)**0.5)))

c     Cf only for tau_r
          cf(pft)=cf(pft)+cf_1hr


c     livegrass consumption
      if (.not.tree(pft)) then
c      if(moist_1hr.le.0.18) then

c  ACHTUNG can we use dlm_lg directly? Thresholds don't match, we set 
c  dlm_lg to zero at dw1=0.3!!!!!!!!!

      if(dlm_lg(d).le.0.18) then
         fc_lg(pft,1)=1.0*(1.0-MINER_TOT)*pot_fc_lg_0(pft,1)*
     *        fire_frac(d)
         fc_lg(pft,2)=pot_fc_lg_0(pft,2)
         fc_lg(pft,3)=pot_fc_lg_0(pft,3)
      else
cc         if(moist_1hr.gt.0.18.and.moist_1hr.le.0.73) then
         if(dlm_lg(d).gt.0.18.and.dlm_lg(d).le.0.73) then
       fc_lg(pft,1)=(1.10-0.62*dlm_lg(d))*(1.0-MINER_TOT)
     *     *pot_fc_lg_0(pft,1)*fire_frac(d)
           fc_lg(pft,2)=pot_fc_lg_0(pft,2)
           fc_lg(pft,3)=pot_fc_lg_0(pft,3)
        else
cc           if(moist_1hr.gt.0.73.and.moist_1hr.le.1.0) then
           if(dlm_lg(d).gt.0.73.and.dlm_lg(d).le.1.0) then
       fc_lg(pft,1)=(2.45-2.45*dlm_lg(d))*(1.0-MINER_TOT)*
!            fc_lg(pft,3)=pot_fc_lg(pft,3)
     *                 pot_fc_lg_0(pft,1)*fire_frac(d)
            fc_lg(pft,2)=pot_fc_lg_0(pft,2)
            fc_lg(pft,3)=pot_fc_lg_0(pft,3)
           else
            fc_lg(pft,1)=0.0
            fc_lg(pft,2)=0.0
            fc_lg(pft,3)=0.0
           endif
        endif
      endif
      endif !not a tree

c     10hr fuel consumption
      if(moist_10_100hr.le.0.12) then
        fc_10hr(pft,1)=1.0*(1.0-MINER_TOT)*(fuel_10hr_left(pft,1)/0.45)
     *        *fire_frac(d)
         fc_10hr(pft,2)=fuel_10hr_left(pft,2)
         fc_10hr(pft,3)=fuel_10hr_left(pft,3)
        cf_10hr=1.0
      else
       if(moist_10_100hr.gt.0.12.and.moist_10_100hr.le.0.51) then
          fc_10hr(pft,1)=(1.09-0.72*moist_10_100hr)*
     *           (1.0-MINER_TOT)*(fuel_10hr_left(pft,1)/0.45)*
     *           fire_frac(d)
          fc_10hr(pft,2)=fuel_10hr_left(pft,2)
          fc_10hr(pft,3)=fuel_10hr_left(pft,3)
          cf_10hr=1.09-0.72*moist_10_100hr
         else
           if(moist_10_100hr.gt.0.51.and.moist_10_100hr.le.1.0) then
             fc_10hr(pft,1)=(1.47-1.47*moist_10_100hr)*
     *         (1.0-MINER_TOT)*(fuel_10hr_left(pft,1)/0.45)*fire_frac(d)
             fc_10hr(pft,2)=fuel_10hr_left(pft,2)
             fc_10hr(pft,3)=fuel_10hr_left(pft,3)
             cf_10hr=1.47-1.47*moist_10_100hr
           else
             fc_10hr(pft,1)=0.0
             fc_10hr(pft,2)=0.0
             fc_10hr(pft,3)=0.0
             cf_10hr=0.0
           endif
         endif
      endif


c    time required for cambial kill, 10hr fuel
c         tau_l(pft)=tau_l(pft)+((fuel_10hr_left(pft)/0.45/1e4)*
c     *                     (1.0-((1.0-cf_10hr)**0.5)))

      

c     100hr fuel consumption
      if(moist_10_100hr.le.0.38) then
        fc_100hr(pft,1)=(0.98-0.85*moist_10_100hr)*(1.0-MINER_TOT)
     *                *(fuel_100hr_left(pft,1)/0.45)*fire_frac(d)
        fc_100hr(pft,2)=fuel_100hr_left(pft,2)
        fc_100hr(pft,3)=fuel_100hr_left(pft,3)
        cf_100hr=0.98-0.85*moist_10_100hr
      else
        if(moist_10_100hr.gt.0.38.and.moist_10_100hr.le.1.0) then
          fc_100hr(pft,1)=(1.06-1.06*moist_10_100hr)*(1.0-MINER_TOT)
     *               *(fuel_100hr_left(pft,1)/0.45)*fire_frac(d)
          fc_100hr(pft,2)=fuel_100hr_left(pft,2)
          fc_100hr(pft,3)=fuel_100hr_left(pft,3)
          cf_100hr=1.06-1.06*moist_10_100hr
        else
           fc_100hr(pft,1)=0.0
           fc_100hr(pft,2)=0.0
           fc_100hr(pft,3)=0.0
          cf_100hr=0.0
        endif
      endif


c    time required for cambial kill, 100hr fuel
c         tau_l(pft)=tau_l(pft)+((fuel_100hr(pft)/0.45/1e4)*
c     *                     (1.0-((1.0-cf_100hr)**0.5)))

c    calculating Cf from fuel classes, tau_r=Cf/Gamma tau_l=2*tau_r
        cf(pft)=(cf_1hr+cf_10hr+cf_100hr)/3.0


c    Peterson& Ryan: 39.4*tau_l, assuming particle density of 510, 
c    we have 513, thus 39.09*tau_l
c         tau_l(pft)=39.09*tau_l(pft)

c    1000hr fuel consumption, not influencing rate of spread or I_surface (Rothermel 1972)

        fc_1000hr(pft,1)=(-0.8*mw_weight+0.8)*(1.0-MINER_TOT)*
     *                    fuel_1000hr_left(pft,1)*fire_frac(d)
        fc_1000hr(pft,2)=fuel_1000hr_left(pft,2)
        fc_1000hr(pft,3)=fuel_1000hr_left(pft,3)


c    Allan: Approximate form. No data.

c    total fuel consumption (without 1000hr fuel) in g biomass per m?²!!!
c    Used to calculate fire intensity in the FLAMING FRONT.
      fuel_consum=fuel_consum+fc_1hr_leaf(pft,1)+fc_1hr_wood(pft,1)
     *       +fc_10hr(pft,1)+fc_100hr(pft,1)
                 

        fuel_consum_all=fuel_consum_all
     *            +fc_1hr_leaf(pft,1)+fc_1hr_wood(pft,1)+
     *            fc_10hr(pft,1)+fc_100hr(pft,1)+fc_1000hr(pft,1)  !to calculate the total fuel consumed for the further calculation of combustion efficiency.
        d_fuel_all=d_fuel_all+
     *             ((fuel_1hr_leaf_left(pft,1)+
     *             fuel_1hr_wood_left(pft,1)+
     *             fuel_10hr_left(pft,1)+
     *             fuel_100hr_left(pft,1))
     *             /0.45+fuel_1000hr_left(pft,1))*fire_frac(d) 



       endif  !present

      enddo !pft


      return
      end	!fuel_consumption

#ifdef NOT_DEF
c subroutine not used
c TODO : check initialization of locals if enable
c-------------------------------------------------------------------------------
c    subroutine for calculation of fuel consumption in the area affected by fire

      subroutine fuel_consumption_netmoist(npft,present,fuel_consum,
     *  fire_frac,fuel_1hr,fuel_10hr,fuel_100hr,fuel_1000hr,livegrass,
     *  pot_fc_lg,tree,fuel_1hr_total,moistfactor,dlm,d,MINER_TOT,fc_1hr
     *  ,fc_lg,fc_10hr,fc_100hr,fc_1000hr,tau_l,char_moistfactor,
     *   dlm_lg, dlm_1hr, dlm_10hr, dlm_100hr, dlm_1000hr,
     *   moistfactor_livegrass, moistfactor_1hr,moistfactor_10hr,
     *   moistfactor_100hr, moistfactor_1000hr)


      implicit none

c     PARAMETERS
      integer nc         ! nc is added for nco2
      integer nco2
         parameter (nco2=3)

      integer d,npft
      real fuel_1hr_total,livegrass
      real moistfactor,dw1(1:365),fire_frac,dlm(1:365)
      real pot_fc_lg(1:npft,1:nco2),fuel_1hr(1:npft,1:nco2)
      real fuel_10hr(1:npft,1:nco2)
      real fuel_100hr(1:npft,1:nco2),fuel_1000hr(1:npft,1:nco2)
      real MINER_TOT
      real fuel_consum,fc_1hr(1:npft,1:nco2),fc_lg(1:npft,1:nco2)
      real fc_10hr(1:npft,1:nco2)
      real fc_100hr(1:npft,1:nco2),fc_1000hr(1:npft,1:nco2)
      real tau_l(1:npft)
      real char_moistfactor
      real moistfactor_livegrass, moistfactor_1hr,moistfactor_10hr
      real moistfactor_100hr, moistfactor_1000hr
      real  dlm_lg (1:365), dlm_1hr(1:365), dlm_10hr(1:365)  
      real  dlm_100hr(1:365), dlm_1000hr (1:365)
      logical present(1:npft),tree(1:npft)

c     local variables
      integer pft
       real moist_lg_d1hr,moist_1hr,moist_10_100hr
      real net_moist_lg, net_moist_1hr, net_moist_10hr
      real net_moist_100hr, net_moist_1000hr
      real pot_fc_1hr_total,mw_weight
      real cf_lg, cf_1hr,cf_10hr,cf_100hr,cf_1000hr

c     initialise
      fuel_consum=0.0
      moist_lg_d1hr=0.0
      moist_1hr=0.0
      moist_10_100hr=0.0
      mw_weight=0.0
      cf_1hr=0.0
      cf_10hr=0.0
      cf_100hr=0.0

       cf_1000hr=0
       cf_lg=0
c     dlm_lg as a function of upper soil moisture
       net_moist_lg = dlm_lg(d) / moistfactor_livegrass 
       net_moist_1hr = dlm_1hr(d) / moistfactor_1hr
       net_moist_10hr = dlm_10hr(d) / moistfactor_10hr
       net_moist_100hr = dlm_100hr(d) / moistfactor_100hr
       net_moist_1000hr = dlm_1000hr(d) / moistfactor_1000hr

c     fuel consumption depending on fuel moisture
c     ACHTUNG: This is using upper soil moisture until fuel moisture can be calculated
c     per dead fuel class!!!
       

c        if (fuel_1hr(pft).lt.0) print*,'error line 7245 1hr ',
c     *   'fuel_1hr(pft),fc_1hr(pft)',fuel_1hr(pft),fc_1hr(pft)
c        if (fuel_10hr(pft).lt.0) print*,'error line 7245 10hr ',
c     *   'fuel_10hr(pft),fc_10hr(pft)',fuel_10hr(pft),fc_10hr(pft)
c        if (fuel_100hr(pft).lt.0) print*,'error line 7245 100hr ',
c     *   'fuel_100hr(pft),fc_100hr(pft)',fuel_100hr(pft),fc_100hr(pft)
c        if (fuel_1000hr(pft).lt.0) print*,'error line 7245 1000hr ',
c     *   'fuel_1000hr(pft),fc_1000hr(pft)',
c     *   fuel_1000hr(pft),fc_1000hr(pft)

      do pft=1,npft
       do nc=1,nco2
        fc_1hr(pft,nc)=0.0
        fc_lg(pft,nc)=0.0
        fc_10hr(pft,nc)=0.0
        fc_100hr(pft,nc)=0.0
        fc_1000hr(pft,nc)=0.0
       enddo
        tau_l(pft)=0.0 !cycle through the dead fuel classes and livegrass

       if (present(pft)) then

c     1hr fuel consumption
      if(net_moist_1hr.le.0.18) then
         cf_1hr=1.0 
         fc_1hr(pft,1)= cf_1hr *(1.0-MINER_TOT)*(fuel_1hr(pft,1)/0.45)
     *     *fire_frac
         fc_1hr(pft,2)=fuel_1hr(pft,2)
         fc_1hr(pft,3)=fuel_1hr(pft,3)
      else
         if(net_moist_1hr.gt.0.18.and.net_moist_1hr.le.0.73) then
             cf_1hr=1.10-0.62*net_moist_1hr
             fc_1hr(pft,1)=cf_1hr *(1.0-MINER_TOT)*
     *                   (fuel_1hr(pft,1)/0.45)*fire_frac
             fc_1hr(pft,2)=fuel_1hr(pft,2)
             fc_1hr(pft,3)=fuel_1hr(pft,3)
         else
           if(net_moist_1hr.gt.0.73.and.net_moist_1hr.le.1.0) then
            cf_1hr=2.45-2.45*net_moist_1hr
            fc_1hr(pft,1)=cf_1hr*(1.0-MINER_TOT)*
     *                   (fuel_1hr(pft,1)/0.45)*fire_frac
            fc_1hr(pft,2)=fuel_1hr(pft,2)
            fc_1hr(pft,3)=fuel_1hr(pft,3)
           else !includes unrealistic case where Mf > Me !A
            fc_1hr(pft,1)=0.0
            fc_1hr(pft,2)=0.0
            fc_1hr(pft,3)=0.0
            cf_1hr=0.0
           endif
         endif
      endif

c    time required for cambial kill, 1hr fuel [gBiomass/cm?²]
c    /1e4 to get from 1m2 to cm2, the units used by P&R(1986) 
         tau_l(pft)=tau_l(pft)+((fuel_1hr(pft,1)/0.45/1e4)*
     *                     (1.0-((1.0-cf_1hr)**0.5)))

c     livegrass consumption. Same slopes as 1hr dead. A.
      if (.not.tree(pft)) then
      if(net_moist_lg.le.0.18) then
         fc_lg(pft,1)=1.0*(1.0-MINER_TOT)*
     *              pot_fc_lg(pft,1)*fire_frac
         fc_lg(pft,2)=pot_fc_lg(pft,2)
         fc_lg(pft,3)=pot_fc_lg(pft,3)
      else
         if(net_moist_lg.gt.0.18.and.net_moist_lg.le.0.73) then
           fc_lg(pft,1)=(1.10-0.62*net_moist_lg)*(1.0-MINER_TOT)*
     *                 pot_fc_lg(pft,1)*fire_frac
           fc_lg(pft,2)=pot_fc_lg(pft,2)
           fc_lg(pft,3)=pot_fc_lg(pft,3)
         else
           if(net_moist_lg.gt.0.73.and.net_moist_lg.le.1.0) then
            fc_lg(pft,1)=(2.45-2.45*net_moist_lg)*(1.0-MINER_TOT)*
     *                 pot_fc_lg(pft,1)*fire_frac
            fc_lg(pft,2)=pot_fc_lg(pft,2)
            fc_lg(pft,3)=pot_fc_lg(pft,3)
           else  !includes unrealistic case where Mf > Me
            fc_lg(pft,1)=0.0
            fc_lg(pft,2)=0.0
            fc_lg(pft,3)=0.0
           endif
         endif
      endif
      fc_lg(pft,1)= 0.05 * pot_fc_lg(pft,1) ! not quite sure about this and the following statement, MAS, 29/04/07
      fc_lg(pft,2)=pot_fc_lg(pft,2)      
      fc_lg(pft,3)=pot_fc_lg(pft,3)      
      fc_lg(pft,1)=0.0
      fc_lg(pft,2)=0.0
      fc_lg(pft,3)=0.0
      endif

c     10hr fuel consumption
      if(net_moist_10hr.le.0.12) then
        cf_10hr=1.0 
       fc_10hr(pft,1)=cf_10hr*(1.0-MINER_TOT)*
     *              (fuel_10hr(pft,1)/0.45)*fire_frac
      fc_10hr(pft,2)=fuel_10hr(pft,2)
      fc_10hr(pft,3)=fuel_10hr(pft,3)
      else
       if(net_moist_10hr.gt.0.12.and.net_moist_10hr.le.0.51) then
          cf_10hr=1.09-0.72*net_moist_10hr
          fc_10hr(pft,1)=cf_10hr*(1.0-MINER_TOT)*
     *           (fuel_10hr(pft,1)/0.45)*fire_frac
          fc_10hr(pft,2)=fuel_10hr(pft,2)
          fc_10hr(pft,3)=fuel_10hr(pft,3)
         else
           if(net_moist_10hr.gt.0.51.and.net_moist_10hr.le.1.0) then
             cf_10hr=1.47-1.47*net_moist_10hr
             fc_10hr(pft,1)=(1.47-1.47*net_moist_10hr)*
     *             (1.0-MINER_TOT)*(fuel_10hr(pft,1)/0.45)*fire_frac
             fc_10hr(pft,2)=fuel_10hr(pft,2)
             fc_10hr(pft,3)=fuel_10hr(pft,3)
           else !includes unrealistic case where Mf > Me
             fc_10hr(pft,1)=0.0
             fc_10hr(pft,2)=0.0
             fc_10hr(pft,3)=0.0
             cf_10hr=0.0
           endif
         endif
      endif

c    time required for cambial kill, 10hr fuel
         tau_l(pft)=tau_l(pft)+((fuel_10hr(pft,1)/0.45/1e4)*
     *                     (1.0-((1.0-cf_10hr)**0.5)))

c     100hr fuel consumption
      if(net_moist_100hr.le.0.38) then
        cf_100hr=0.98-0.85*net_moist_100hr
        fc_100hr(pft,1)=cf_100hr*(1.0-MINER_TOT)*
     *               (fuel_100hr(pft,1)/0.45)*fire_frac
        fc_100hr(pft,2)=fuel_100hr(pft,2)
        fc_100hr(pft,3)=fuel_100hr(pft,3)
      else
        if(net_moist_100hr.gt.0.38.and.net_moist_100hr.le.1.0) then
           cf_100hr=1.06-1.06*net_moist_100hr
           fc_100hr(pft,1)=cf_100hr*(1.0-MINER_TOT)
     *               *(fuel_100hr(pft,1)/0.45)*fire_frac
           fc_100hr(pft,2)=fuel_100hr(pft,2)
           fc_100hr(pft,3)=fuel_100hr(pft,3)
        else !includes unrealistic case where Mf > Me
           fc_100hr(pft,1)=0.0
           fc_100hr(pft,2)=0.0
           fc_100hr(pft,3)=0.0
           cf_100hr=0.0
        endif
      endif


c    time required for cambial kill, 100hr fuel
         tau_l(pft)=tau_l(pft)+((fuel_100hr(pft,1)/0.45/1e4)*
     *                     (1.0-((1.0-cf_100hr)**0.5)))

c    Peterson& Ryan: 39.4*tau_l, assuming particle density of 510, we have 513, thus 39.09*tau_l
         tau_l(pft)=39.09*tau_l(pft)

c    1000hr fuel consumption, not influencing rate of spread or I_surface (Rothermel 1972)

        if (net_moist_1000hr.le.1) then
         fc_1000hr(pft,1)=(-0.8*net_moist_1000hr+0.8)*(1.0-MINER_TOT)*
     *                    fuel_1000hr(pft,1)*fire_frac !gC/m2
         fc_1000hr(pft,2)=fuel_1000hr(pft,2)
         fc_1000hr(pft,3)=fuel_1000hr(pft,3)
        else
          fc_1000hr(pft,1) = 0
          fc_1000hr(pft,2) = 0
          fc_1000hr(pft,3) = 0
        endif
c    Allan: Approximate form. No data.

c    total fuel consumption (without 1000hr fuel) in g biomass per m?²!!!
c    Used to calculate fire intensity in the FLAMING FRONT.
      fuel_consum=fuel_consum+fc_1hr(pft,1)+fc_10hr(pft,1)
     *       +fc_100hr(pft,1)
                 

       endif  !present

      enddo !pft

c        if (fuel_1hr(pft).lt.0) print*,'error line 7472 1hr ',
c     *   'fuel_1hr(pft),fc_1hr(pft)',fuel_1hr(pft),fc_1hr(pft)
c        if (fuel_10hr(pft).lt.0) print*,'error line 7472 10hr ',
c     *   'fuel_10hr(pft),fc_10hr(pft)',fuel_10hr(pft),fc_10hr(pft)
c        if (fuel_100hr(pft).lt.0) print*,'error line 7472 100hr ',
c     *   'fuel_100hr(pft),fc_100hr(pft)',fuel_100hr(pft),fc_100hr(pft)
c        if (fuel_1000hr(pft).lt.0) print*,'error line 7472 1000hr ',
c     *   'fuel_1000hr(pft),fc_1000hr(pft)',
c     *   fuel_1000hr(pft),fc_1000hr(pft)

      return
      end


#endif


c -----------------------------------------------------------------------------
c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE ESTABLISHMENT
c     Establishment of new individuals (saplings) of woody PFTs,	
c     grass establishment, removal of PFTs not adapted to current climate,
c     update of individual structure and FPC.

      subroutine establishment(pftpar,present,survive,estab,nind,
     *  lm_ind,sm_ind,rm_ind,hm_ind,lm_sapl,sm_sapl,rm_sapl,hm_sapl,
     *  crownarea,fpc_grid,lai_ind,height,dbh,dbh_class,tau_c,cl_t,
     *  BTparam1,BTparam2,BTmode0, ! Doug 03/09
     *  sla,wooddens,
     *  latosa,mprec,reinickerp,
     *  litter_ag_leaf,litter_ag_wood,litter_bg,
     *  fuel_1hr_leaf,fuel_1hr_wood,
     *  fuel_10hr,fuel_100hr,fuel_1000hr,
     *  fuel_1hr_leaf_inc_pos,fuel_1hr_leaf_inc_neg,               !Doug 03/09
     *  fuel_1hr_wood_inc_pos,fuel_1hr_wood_inc_neg, 
     *  fuel_10hr_inc,fuel_100hr_inc,fuel_1000hr_inc,    !Doug 03/09
     *  pas,crop,                                             !Doug 06/09
     *  tree,allom1,allom2,allom3,acflux_estab,
     *  leafondays,leafoffdays,leafon,mnpp,anpp,mnpp_add,anpp_add,
     *  year)
      implicit none

c     PARAMETERS
      integer npft,npftpar,nsoilpar
!        parameter (npft=13,npftpar=50,nsoilpar=7)
        parameter (npft=13,npftpar=59,nsoilpar=7)
      integer nco2
         parameter (nco2=3)
      real pi
        parameter (pi=3.14159265)
      real aprec_min_estab                !minimum annual precipitation for
        parameter (aprec_min_estab=100.0) !establishment (mm)
      real estab_max                      !maximum sapling establishment rate
*        parameter (estab_max=0.12)        !(indiv/m2)
        parameter (estab_max=0.24)
      real nind_min                       !minimum individual density for
        parameter (nind_min=1.0E-10)      !persistence of PFT (indiv/m2)
      real eps
        parameter (eps=1.e-6)

c     ARGUMENTS
      real pftpar(1:npft,1:npftpar)
      logical present(1:npft),survive(1:npft),estab(1:npft)
      real nind(1:npft)
      real lm_ind(1:npft,1:nco2),sm_ind(1:npft,1:nco2)
      real rm_ind(1:npft,1:nco2),hm_ind(1:npft,1:nco2)
      real lm_sapl(1:npft,1:nco2),sm_sapl(1:npft,1:nco2)
      real rm_sapl(1:npft,1:nco2),hm_sapl(1:npft,1:nco2)
      real crownarea(1:npft)
      real fpc_grid(1:npft)
      real lai_ind(1:npft)
      real height(1:npft)
      real dbh(1:npft),tau_c(0:4,1:npft),cl_t(0:4,1:npft)
      REAL BTparam1(1:npft,1:3)        ! Doug 02/13: lower, medium and upper bound describing 1-50-99% quantile p1 in Bark thickness equation BT=p1+p2*DBH
      REAL BTparam2(1:npft,1:3)        ! Doug 02/13: "" for p2
      REAL BTmode0(1:npft,1:2)         ! Doug 02/13: 50% quatile starting point incase of establishment or tree death.  
      real sla(1:npft)
      real wooddens,latosa,reinickerp
      real mprec(1:12)
      REAL litter_ag_leaf(1:npft,1:nco2)
      REAL litter_ag_wood(1:npft,1:nco2)
      REAL litter_bg(1:npft,1:nco2)
      REAL fuel_1hr_leaf(1:npft,1:nco2)
      REAL fuel_1hr_wood(1:npft,1:nco2)
      REAL fuel_10hr(1:npft,1:nco2)
      real fuel_100hr(1:npft,1:nco2)
      real fuel_1000hr(1:npft,1:nco2)

      REAL fuel_1hr_leaf_inc_pos(1:npft,1:nco2,1:12)    !Doug 03/09
      REAL fuel_1hr_leaf_inc_neg(1:npft,1:nco2,1:12)
      REAL fuel_1hr_wood_inc_pos(1:npft,1:nco2,1:12)    !Doug 03/09
      REAL fuel_1hr_wood_inc_neg(1:npft,1:nco2,1:12)
      REAL fuel_10hr_inc(1:npft,1:nco2,1:12)
      REAL fuel_100hr_inc(1:npft,1:nco2,1:12)
      REAL fuel_1000hr_inc(1:npft,1:nco2,1:12)
      REAL temp_fuel(1:5,1:npft,1:nco2)

      REAL pas,crop                                     !Doug 06/09

      logical tree(1:npft)
      real allom1,allom2,allom3
      real acflux_estab(1:nco2)
      integer leafondays(1:npft),leafoffdays(1:npft)
      logical leafon(1:npft)
      real anpp(1:npft,1:nco2),mnpp(1:12,1:npft,1:nco2)
      real mnpp_add(1:12,1:nco2),anpp_add(1:nco2)        !NPPs of PFTs that are killed
      integer year


c     LOCAL VARIABLES
      integer pft,m
      integer npft_estab  !number of regenerating tree PFTs
      real aprec          !annual precipitation (mm)
      real fpc_tree_total,fpc_tree_new !total grid FPC for tree PFTs
      real estab_rate     !sapling establishment rate over area
                          !available for establishment (indiv/m2)
      real estab_grid     !grid-level establishment rate (indiv/m2)
      real nind_old       !number of individuals /m2 before establishment
      real stemdiam       !stem diameter (m)
      real sm_ind_temp    !preliminary sapwood mass (g/m2)
      real fpc_ind        !individual FPC
      real crownarea_max  !maximum crown area (m2)
      real fpc_total,fpc_new      !total grid FPC
      real bare           !gridcell bare ground fraction
      integer ngrass
      real fpc_grass_total,fpc_grass_new
      real bare_max
      integer class
      !real param1(1:npft),param2(1:npft)  !Kirsten: parameter for bark thickness
      real bt(0:4,1:npft),crown(1:npft)      !Kirsten: bark thickness (cm?²)
      real dbh_class(0:4,1:npft),height_class(0:4,1:npft)
      real temp,litter_inc,all
      real fpc_total_new
      integer nc         ! nc is added for nco2

      pft=0
      m=0
      npft_estab=0
      aprec=0.0
      fpc_tree_total=0.0
      fpc_tree_new=0.0
      estab_rate=0.0
      estab_grid=0.0
      nind_old=0.0
      stemdiam=0.0
      sm_ind_temp=0.0
      fpc_ind=0.0
      crownarea_max=0.0
      fpc_total=0.0
      fpc_new=0.0
      bare=0.0
      ngrass=0
      fpc_grass_total=0.0
      fpc_grass_new=0.0
      bare_max=0.0
      class=0
      !param1(:)=0.0
      !param2(:)=0.0
      bt(:,:)=0.0
      crown(:)=0.0
      dbh_class(:,:)=0.0
      height_class(:,:)=0.0
      temp=0.0
      litter_inc=0.0
      all=0.0
      fpc_total_new=0.0
      nc=0

c     Kill PFTs not adapted to current climate, introduce newly "adapted" PFTs

c     Initialize buffers for npps of pfts that go not present
      anpp_add(1)=0.
      anpp_add(2)=0.
      anpp_add(3)=0.
      do m=1,12
         mnpp_add(m,1)=0.
         mnpp_add(m,2)=0.
         mnpp_add(m,3)=0.
      enddo


c     Calculate annual precipitation
      do m=1,12
        aprec=aprec+mprec(m)
      enddo

c     Doug 03/09: Store fuel load values before any changes.
      temp_fuel(1,:,:)=fuel_1hr_leaf
      temp_fuel(2,:,:)=fuel_1hr_wood
      temp_fuel(3,:,:)=fuel_10hr
      temp_fuel(4,:,:)=fuel_100hr
      temp_fuel(5,:,:)=fuel_1000hr

      do pft=1,npft


c         Kirsten: parameter for bark thickness
          crown(pft)=pftpar(pft,44)
c         Doug 02/13: remove BT parameters - now done differently
c          param1(pft)=pftpar(pft,47)
c          param2(pft)=pftpar(pft,48)

        if (present(pft).and.
     *    (.not.survive(pft).or.nind(pft).lt.nind_min)) then     !kill PFT

          present(pft)=.false.

c         Add up NPP of PFTs that are killed in extra balance
          do m=1,12
             temp=mnpp(m,pft,1)
             mnpp_add(m,1)=mnpp_add(m,1)+mnpp(m,pft,1)
             if (mnpp_add(m,1).gt.0.) then
                do nc=2,nco2
                   mnpp_add(m,nc)=(mnpp_add(m,nc)*temp+mnpp(m,pft,1)
     *                  *mnpp(m,pft,nc))/mnpp_add(m,1)
                enddo
             endif
             temp=mnpp(m,pft,1)
             anpp_add(1)=anpp_add(1)+mnpp(m,pft,1)
             if (anpp_add(1).gt.0.) then
                do nc=2,nco2
                   anpp_add(nc)=(anpp_add(nc)*temp+mnpp(m,pft,1)
     *                  *mnpp(m,pft,nc))/anpp_add(1)
                enddo
             endif
          enddo

c         Add killed biomass to litter

          if (tree(pft)) then
		    temp=litter_ag_leaf(pft,1)            
            litter_ag_leaf(pft,1)=litter_ag_leaf(pft,1)+nind(pft)
     *           *lm_ind(pft,1)
            IF (litter_ag_leaf(pft,1) .gt. 0.0) THEN
              DO nc=2,nco2
               litter_inc=(lm_ind(pft,1)*lm_ind(pft,nc))
               litter_ag_leaf(pft,nc)=(litter_ag_leaf(pft,nc)*temp
     *              +litter_inc*nind(pft))/litter_ag_leaf(pft,1)
              END DO
            END IF
			
            temp=litter_ag_wood(pft,1)            
            litter_ag_wood(pft,1)=litter_ag_wood(pft,1)+nind(pft)
     *           *(sm_ind(pft,1)+hm_ind(pft,1))
            IF (litter_ag_wood(pft,1) .gt. 0.0) THEN
              DO nc=2,nco2
               litter_inc=(sm_ind(pft,1)*sm_ind(pft,nc)+hm_ind(pft,1)
     *              *hm_ind(pft,nc))
               litter_ag_wood(pft,nc)=(litter_ag_wood(pft,nc)*temp
     *              +litter_inc*nind(pft))/litter_ag_wood(pft,1)
              END DO
            END IF
		  
c            temp=litter_ag(pft,1)            
c            litter_ag(pft,1)=litter_ag(pft,1)+nind(pft)
c     *           *(sm_ind(pft,1)+hm_ind(pft,1)+lm_ind(pft,1))
c            if (litter_ag(pft,1) .gt. 0.0) then
c            do nc=2,nco2
c               litter_inc=(sm_ind(pft,1)*sm_ind(pft,nc)+hm_ind(pft,1)
c     *              *hm_ind(pft,nc)+lm_ind(pft,1)*lm_ind(pft,nc))
c               litter_ag(pft,nc)=(litter_ag(pft,nc)*temp
c     *              +litter_inc*nind(pft))/litter_ag(pft,1)
c            enddo
c            endif
c         KIRSTEN: per fuel class
            temp=fuel_1hr_leaf(pft,1)
            fuel_1hr_leaf(pft,1)=fuel_1hr_leaf(pft,1)+
     *           lm_ind(pft,1)*nind(pft)
            IF (fuel_1hr_leaf(pft,1) .gt. 0.0) THEN
              DO nc=2,nco2
                litter_inc=lm_ind(pft,1)*lm_ind(pft,nc)
     *              +sm_ind(pft,1)*sm_ind(pft,nc)*0.045
     *              +hm_ind(pft,1)*hm_ind(pft,nc)+0.045
                fuel_1hr_leaf(pft,nc)=(fuel_1hr_leaf(pft,nc)*temp
     *              +litter_inc*nind(pft))/fuel_1hr_leaf(pft,1)
              END DO
            END IF
			
            temp=fuel_1hr_wood(pft,1)
            fuel_1hr_wood(pft,1)=fuel_1hr_wood(pft,1)+
     *           (0.045*sm_ind(pft,1)+ 0.045*hm_ind(pft,1))*nind(pft)
            if (fuel_1hr_wood(pft,1) .gt. 0.0) then
            do nc=2,nco2
               litter_inc=sm_ind(pft,1)*sm_ind(pft,nc)*0.045
     *              +hm_ind(pft,1)*hm_ind(pft,nc)+0.045
               fuel_1hr_wood(pft,nc)=(fuel_1hr_wood(pft,nc)*temp
     *              +litter_inc*nind(pft))/fuel_1hr_wood(pft,1)
            enddo
            endif
			
c            temp=fuel_1hr(pft,1)
c            fuel_1hr(pft,1)=fuel_1hr(pft,1)+(lm_ind(pft,1)+
c     *           0.045*sm_ind(pft,1)+ 0.045*hm_ind(pft,1))*nind(pft)
c            if (fuel_1hr(pft,1) .gt. 0.0) then
c            do nc=2,nco2
c               litter_inc=lm_ind(pft,1)*lm_ind(pft,nc)
c     *              +sm_ind(pft,1)*sm_ind(pft,nc)*0.045
c     *              +hm_ind(pft,1)*hm_ind(pft,nc)+0.045
c               fuel_1hr(pft,nc)=(fuel_1hr(pft,nc)*temp
c     *              +litter_inc*nind(pft))/fuel_1hr(pft,1)
c            enddo
c            endif

            temp=fuel_10hr(pft,1)
            fuel_10hr(pft,1)=fuel_10hr(pft,1)+((0.075*
     *           sm_ind(pft,1)+0.075*hm_ind(pft,1))*nind(pft))
            if (fuel_10hr(pft,1) .gt. 0.0) then
            do nc=2,nco2
               litter_inc=sm_ind(pft,1)*sm_ind(pft,nc)*0.075
     *              +hm_ind(pft,1)*hm_ind(pft,nc)+0.075
               fuel_10hr(pft,nc)=(fuel_10hr(pft,nc)*temp
     *              +litter_inc*nind(pft))/fuel_10hr(pft,1)
            enddo
            endif

            temp=fuel_100hr(pft,1)
            fuel_100hr(pft,1)=fuel_100hr(pft,1)+((0.21*
     *           sm_ind(pft,1) +0.21*hm_ind(pft,1))*nind(pft))
            if (fuel_100hr(pft,1) .gt. 0.0) then
            do nc=2,nco2
               litter_inc=sm_ind(pft,1)*sm_ind(pft,nc)*0.21
     *              +hm_ind(pft,1)*hm_ind(pft,nc)+0.21
               fuel_100hr(pft,nc)=(fuel_100hr(pft,nc)*temp
     *              +litter_inc*nind(pft))/fuel_100hr(pft,1)
            enddo
            endif

            temp=fuel_1000hr(pft,1)
            fuel_1000hr(pft,1)=fuel_1000hr(pft,1)+((0.67
     *           *sm_ind(pft,1)+0.67*hm_ind(pft,1))*nind(pft))
            if (fuel_1000hr(pft,1) .gt. 0.0) then
            do nc=2,nco2
               litter_inc=sm_ind(pft,1)*sm_ind(pft,nc)*0.67
     *              +hm_ind(pft,1)*hm_ind(pft,nc)+0.67
               fuel_1000hr(pft,nc)=(fuel_1000hr(pft,nc)*temp
     *              +litter_inc*nind(pft))/fuel_1000hr(pft,1)
            enddo
            endif

          else  !grasses
             temp=litter_ag_leaf(pft,1)
             litter_ag_leaf(pft,1)=litter_ag_leaf(pft,1)+lm_ind(pft,1)
     *            *nind(pft)
             if (litter_ag_leaf(pft,1) .gt. 0.0) then
             do nc=2,nco2
                litter_inc=lm_ind(pft,1)*lm_ind(pft,nc)
                litter_ag_leaf(pft,nc)=(litter_ag_leaf(pft,nc)*temp
     *               +litter_inc*nind(pft))/litter_ag_leaf(pft,1)
             enddo
             endif
c         KIRSTEN: 1hr fuel class
             temp=fuel_1hr_leaf(pft,1)
             fuel_1hr_leaf(pft,1)=fuel_1hr_leaf(pft,1)+lm_ind(pft,1)
     *            *nind(pft)
             if (fuel_1hr_leaf(pft,1) .gt. 0.0) then
             do nc=2,nco2
                litter_inc=lm_ind(pft,1)*lm_ind(pft,nc)
                fuel_1hr_leaf(pft,nc)=(fuel_1hr_leaf(pft,nc)*temp
     *               +litter_inc*nind(pft))/fuel_1hr_leaf(pft,1)
             enddo
             endif
          endif

          temp=litter_bg(pft,1)
          litter_bg(pft,1)=litter_bg(pft,1)+rm_ind(pft,1)*nind(pft)
          if (litter_bg(pft,1) .gt. 0.0) then
          do nc=2,nco2
             litter_bg(pft,nc)=(litter_bg(pft,nc)*temp+rm_ind(pft,nc)
     *            *rm_ind(pft,1)*nind(pft))/litter_bg(pft,1)
          enddo
          endif

        elseif (.not.present(pft).and.survive(pft).and.estab(pft)
     *    .and.aprec.ge.aprec_min_estab) then

c         Introduce PFT if conditions suitable for establishment

          present(pft)=.true.
          if (tree(pft)) then
            nind(pft)=0.0
          else
            nind(pft)=1.0   !each grass PFT = 1 "individual"
          endif

          do nc=1,nco2
           lm_ind(pft,nc)=0.0
           sm_ind(pft,nc)=0.0
           rm_ind(pft,nc)=0.0
           hm_ind(pft,nc)=0.0
           anpp(pft,nc)=0.0
          enddo
          fpc_grid(pft)=0.0
          do m=1,12
             mnpp(m,pft,1)=0.0
             mnpp(m,pft,2)=0.0
             mnpp(m,pft,3)=0.0
          enddo


          if(.not.tree(pft)) crownarea(pft)=1.0
          leafon(pft)=.true.
          leafondays(pft)=0.0
          leafoffdays(pft)=0.0

        endif

      enddo


c     SAPLING AND GRASS ESTABLISHMENT

c     Calculate total woody FPC and number of woody PFTs present and
c     able to establish

      fpc_tree_total=0.0
      fpc_tree_new=0.0
      fpc_grass_total=0.0
      fpc_grass_new=0.0
      npft_estab=0
      fpc_total=0.0
      fpc_new=0.0
      ngrass=0

      do pft=1,npft
         if (present(pft)) then
            if (tree(pft)) then
              fpc_tree_total=fpc_tree_total+fpc_grid(pft)
              IF (estab(pft)) npft_estab=npft_estab+
     *          (1-0.9*pftpar(pft,59))
                ! Doug 03/13: scale resprouting establishment success. 
                ! (1-0.5*rs) means that 100% rs success pft only has 50% 
                ! of the establistment
            else
              ngrass=ngrass+1
              fpc_grass_total=fpc_grass_total+fpc_grid(pft)
            endif
            fpc_total=fpc_total+fpc_grid(pft)
         endif
      enddo ! pft


      acflux_estab(1)=0.0
      acflux_estab(2)=0.0
      acflux_estab(3)=0.0

c     Prohibit establishment under extreme temperature or water stress.

      do pft=1,npft

      if (aprec.ge.aprec_min_estab.and.
     *  npft_estab.gt.0) then


c       Calculate establishment rate over available space, per tree PFT
c       Maximum establishment rate reduced by shading as tree FPC approaches 1
c       Total establishment rate partitioned equally among regenerating woody
c       PFTs

        estab_rate=estab_max*(1.0-exp(5.0*(fpc_tree_total-1.0)))/
     *    real(npft_estab)

c       Calculate grid-level establishment rate per woody PFT
c       Space available for woody PFT establishment is proportion of grid cell
c       not currently occupied by woody PFTs.
c       Doug 06/09: Grass is more likely to estabilish on pastured land, so
c         reduce establishment rate for trees according to proportion of pasture

        estab_grid=estab_rate*(1.0-fpc_tree_total)!*(1-pas)  !Doug 06/09: *(1-pas) added

      else  !unsuitable climate for establishment

        estab_grid=0.0

      endif

        if (present(pft).and.tree(pft).and.estab(pft)) then
          crownarea_max=pftpar(pft,20)

c         Add new saplings to current population

          nind_old=nind(pft)

          nind(pft)=nind_old+estab_grid*(1-0.9*pftpar(pft,59))
          ! Doug 03/13: scale resprouting establishment success. 
          ! (1-0.5*rs) means that 100% rs success pft only has 50% 
          ! of the establistment

*          if (nind(pft).gt.0.0) then
            temp=lm_ind(pft,1)
            lm_ind(pft,1)=(lm_ind(pft,1)*nind_old+lm_sapl(pft,1)
     *           *estab_grid)/nind(pft)
            if (lm_ind(pft,1).gt.0.0) then
               do nc=2,nco2
                  lm_ind(pft,nc)=(lm_ind(pft,nc)*temp*nind_old+
     *                 lm_sapl(pft,nc)*lm_sapl(pft,1)*estab_grid)/
     *                 (lm_ind(pft,1)*nind(pft))
               enddo
            endif
            temp=sm_ind(pft,1)
            sm_ind_temp=(sm_ind(pft,1)*nind_old+sm_sapl(pft,1)
     *           *estab_grid)/nind(pft)               
            if (sm_ind_temp.gt.0.0) then
               do nc=2,nco2
                  sm_ind(pft,nc)=(sm_ind(pft,nc)*temp*nind_old+
     *                 sm_sapl(pft,nc)*sm_sapl(pft,1)*estab_grid)/
     *                 (sm_ind_temp*nind(pft))
               enddo
            endif
            temp=hm_ind(pft,1)
            hm_ind(pft,1)=(hm_ind(pft,1)*nind_old+hm_sapl(pft,1)
     *           *estab_grid)/nind(pft)
            if (hm_ind(pft,1).gt.0.0) then
               do nc=2,nco2
                  hm_ind(pft,nc)=(hm_ind(pft,nc)*temp*nind_old+
     *                 hm_sapl(pft,nc)*hm_sapl(pft,1)*estab_grid)/                  
     *                 (hm_ind(pft,1)*nind(pft))
               enddo
            endif
            temp=rm_ind(pft,1)
            rm_ind(pft,1)=(rm_ind(pft,1)*nind_old+rm_sapl(pft,1)
     *           *estab_grid)/nind(pft)
            if (rm_ind(pft,1).gt.0.0) then
               do nc=2,nco2
                  rm_ind(pft,nc)=(rm_ind(pft,nc)*temp*nind_old+
     *                 rm_sapl(pft,nc)*rm_sapl(pft,1)*estab_grid)/
     *                 (rm_ind(pft,1)*nind(pft))
               enddo
            endif
*          else
*           lm_ind(pft)=0.0
*           sm_ind(pft)=0.0
*           hm_ind(pft)=0.0
*           rm_ind(pft)=0.0
*          endif


c         Accumulate biomass increment due to sapling establishment
 
            temp=acflux_estab(1)
            all=(lm_sapl(pft,1)+sm_sapl(pft,1)+hm_sapl(pft,1)
     *           +rm_sapl(pft,1))
            if (all*estab_grid.gt.eps) then
               acflux_estab(1)=acflux_estab(1)+all*estab_grid
               if (acflux_estab(1).gt.0.0.and.temp.ge.0.) then
                  do nc=2,nco2
                     acflux_estab(nc)=(acflux_estab(nc)*temp+
     *                  (lm_sapl(pft,1)*lm_sapl(pft,nc)+sm_sapl(pft,1)*
     *                  sm_sapl(pft,nc)+hm_sapl(pft,1)*hm_sapl(pft,nc)+
     *                  rm_sapl(pft,1)*rm_sapl(pft,nc))*estab_grid)/
     *                  acflux_estab(1)
                  enddo
               elseif (acflux_estab(1).gt.0.0.and.temp.lt.0.) then
                  do nc=2,nco2
                     acflux_estab(nc)=(lm_sapl(pft,1)*lm_sapl(pft,nc)
     *                 +sm_sapl(pft,1)*sm_sapl(pft,nc)+hm_sapl(pft,1)
     *                  *hm_sapl(pft,nc)+rm_sapl(pft,1)*rm_sapl(pft,nc))
     *                  /all
                  enddo
               endif
            endif
 
c         Calculate height, diameter and crown area for new average
c         individual such that the basic allometric relationships (A-C below)
c         are satisfied.

c         (A) (leaf area) = latosa * (sapwood xs area)
c                (Pipe Model, Shinozaki et al. 1964a,b; Waring et al 1982)
c         (B) (leaf mass) = lmtorm * (root mass)
c         (C) height = allom2 * (stem diameter)**allom3
c                (source?)
c         (D) (crown area) = min (allom1 * (stem diameter)**reinickerp,
c                                 crownarea_max)

c         From (A),
c          (1) sap_xsa = lm_ind * sla / latosa
c          (2) wooddens = (sm_ind + hm_ind) / stemvolume
c          (3) stemvolume = stem_xsa * height
c         From (1), (2) & (3),
c          (4) stem_xsa = (sm_ind + hm_ind) / wooddens / height
c          (5) stem_xsa = pi * (stemdiam**2) / 4
c         From (5),
c          (6) stemdiam = ( 4 * stem_xsa / pi )**0.5
c         From (4) & (6),
c          (7) stemdiam = ( 4 * (sm_ind + hm_ind) / wooddens / height /
c                pi )**0.5
c         From (C) & (7),
c          (8) stemdiam = ( 4 * (sm_ind + hm_ind) / wooddens /
c                ( allom2 * stemdiam**allom3 ) / pi )**0.5
c         From (8),
c          (9) stemdiam = ( 4 * (sm_ind + hm_ind ) / wooddens / pi /
c                allom2 )**( 1 / (2 + allom3) )

          stemdiam=(4.0*(sm_ind_temp+hm_ind(pft,1))
     *        /wooddens/pi/allom2)**
     *        (1.0/(2.0+allom3))                 !Eqn 9
          height(pft)=allom2*stemdiam**allom3  !Eqn C
          crownarea(pft)=min(crownarea_max,
     *      allom1*stemdiam**reinickerp)       !Eqn D

c         Kirsten recalculate dbh, bt and tau_c for average individual
c          ACHTUNG: Is this needed here? prob. no influence on fire for the next year...
c    KIRSTEN: put in uniform distribution of dbh = 5 classes
            dbh(pft)=stemdiam*100.0 !in cm

           do class=0,4
             dbh_class(class,pft)=(2.0*dbh(pft)+
     *                  (class*0.25*2.0*dbh(pft))- 0.125*2.0*dbh(pft))

c			Doug 02/13: add more detail here later	 
c          bt(class,pft)=(param1(pft)*dbh_class(class,pft)+param2(pft))
c              tau_c(class,pft)=2.9*(bt(class,pft)**2.0)

              height_class(class,pft)=2*height(pft)-
     *            (class*0.25*2.0*height(pft))+ 0.125*2.0*height(pft)   
              cl_t(class,pft)=height_class(class,pft)*crown(pft)  
           enddo !pseudo age class
		   
c		Doug 02/13: add more detail here later 
           IF (tree(pft)) THEN
             !print*, "wow",tree(pft)
		    ! print*, pft, tree
		     !print*, BTparam1(pft,1:3)
		     !print*, BTparam2(pft,1:3)
             CALL bt_establish(BTparam1(pft,1:3),BTmode0(pft,1),
     *         nind_old,estab_grid)
	 
	         CALL bt_establish(BTparam2(pft,1:3),BTmode0(pft,2),
     *         nind_old,estab_grid)

c            bt(pft)=(param1(pft)*dbh(pft)+param2(pft))
c            tau_c(pft)=2.9*(bt(pft)**2.0)
          END IF

c         Recalculate sapwood mass, transferring excess sapwood to heartwood
c         compartment, if necessary to satisfy Eqn A

          sm_ind(pft,1)=lm_ind(pft,1)*height(pft)*wooddens*sla(pft)
     *         /latosa 
          temp=hm_ind(pft,1)
          hm_ind(pft,1)=hm_ind(pft,1)+(sm_ind_temp-sm_ind(pft,1))
          if (hm_ind(pft,1).gt.0.0) then
             do nc=2,nco2
                hm_ind(pft,nc)=(hm_ind(pft,nc)*temp+(sm_ind_temp-
     *               sm_ind(pft,1))*sm_ind(pft,nc))/hm_ind(pft,1)
             enddo
          endif

c          Update LAI and FPC
           if (crownarea(pft).gt.0.0) then
             lai_ind(pft)=(lm_ind(pft,1)*sla(pft))/crownarea(pft)
           else
             lai_ind(pft)=0.0
           endif
		   
           fpc_ind=(1.0-exp(-0.5*lai_ind(pft)))		 
           fpc_grid(pft)=crownarea(pft)*nind(pft)*fpc_ind
       endif
      enddo !pft

      fpc_tree_total=0.0
      fpc_grass_total=0.0
      fpc_total=0.0

      do pft=1,npft
         if (present(pft)) then
            if (tree(pft)) then
              fpc_tree_total=fpc_tree_total+fpc_grid(pft)
            else
              fpc_grass_total=fpc_grass_total+fpc_grid(pft)
            endif
            fpc_total=fpc_total+fpc_grid(pft)
         endif
      enddo ! pft
      if (fpc_tree_total.gt.0.95) then
        do pft=1,npft
          if(tree(pft).and.present(pft)) then
              nind_old=nind(pft)
              nind(pft)=nind(pft)/(fpc_tree_total/0.95)
              fpc_grid(pft)=fpc_grid(pft)/(fpc_tree_total/0.95)
			  
              temp=litter_ag_leaf(pft,1)
              litter_ag_leaf(pft,1)=litter_ag_leaf(pft,1)+lm_ind(pft,1)
     *             *(nind_old-nind(pft))
              IF (litter_ag_leaf(pft,1).gt.0.0) THEN
                 DO nc=2,nco2
                    litter_ag_leaf(pft,nc)=(litter_ag_leaf(pft,nc)*temp+
     *                   ((lm_ind(pft,1)*lm_ind(pft,nc))*
     *                   (nind_old-nind(pft))))/ litter_ag_leaf(pft,1)
                 END DO
              END IF
			  
              temp=litter_ag_wood(pft,1)
              litter_ag_wood(pft,1)=litter_ag_wood(pft,1)+
     *             (sm_ind(pft,1)+hm_ind(pft,1))*(nind_old-nind(pft))
              IF (litter_ag_wood(pft,1).gt.0.0) THEN
                 DO nc=2,nco2
                    litter_ag_wood(pft,nc)=(litter_ag_wood(pft,nc)*temp+
     *                   ((hm_ind(pft,1)*hm_ind(pft,nc)+sm_ind(pft,1)*
     *                   sm_ind(pft,nc)*(nind_old-nind(pft)))))/
     *                   litter_ag_wood(pft,1)
                 END DO
              END IF
			  
c              temp=litter_ag(pft,1)
c              litter_ag(pft,1)=litter_ag(pft,1)+(lm_ind(pft,1)+
c     *             sm_ind(pft,1)+hm_ind(pft,1))*(nind_old-nind(pft))
c              if (litter_ag(pft,1).gt.0.0) then
c                 do nc=2,nco2
c                    litter_ag(pft,nc)=(litter_ag(pft,nc)*temp+
c     *                   ((hm_ind(pft,1)*hm_ind(pft,nc)+sm_ind(pft,1)*
c     *                   sm_ind(pft,nc)+lm_ind(pft,1)*lm_ind(pft,nc))*
c     *                   (nind_old-nind(pft))))/ litter_ag(pft,1)
c                 enddo
c              endif
c      Kirsten: fuel class
            temp=fuel_1hr_leaf(pft,1)
            fuel_1hr_leaf(pft,1)=fuel_1hr_leaf(pft,1)+
     *           lm_ind(pft,1)*nind(pft)
            IF (fuel_1hr_leaf(pft,1) .gt. 0.0) THEN
              DO nc=2,nco2
                litter_inc=sm_ind(pft,1)*sm_ind(pft,nc)*0.045
     *              +hm_ind(pft,1)*hm_ind(pft,nc)*0.045
                fuel_1hr_leaf(pft,nc)=(fuel_1hr_leaf(pft,nc)*temp
     *              +litter_inc*nind(pft))/fuel_1hr_leaf(pft,1)
              END DO
            END IF
			
            temp=fuel_1hr_wood(pft,1)
            fuel_1hr_wood(pft,1)=fuel_1hr_wood(pft,1)+
     *           (0.045*sm_ind(pft,1)+ 0.045*hm_ind(pft,1))*nind(pft)
            IF (fuel_1hr_wood(pft,1) .gt. 0.0) THEN
              DO nc=2,nco2
                litter_inc=sm_ind(pft,1)*sm_ind(pft,nc)*0.045
     *              +hm_ind(pft,1)*hm_ind(pft,nc)*0.045
                fuel_1hr_wood(pft,nc)=(fuel_1hr_wood(pft,nc)*temp
     *              +litter_inc*nind(pft))/fuel_1hr_wood(pft,1)
              END DO
            END IF
			
C            temp=fuel_1hr(pft,1)
C            fuel_1hr(pft,1)=fuel_1hr(pft,1)+(lm_ind(pft,1)+
C     *           0.045*sm_ind(pft,1)+ 0.045*hm_ind(pft,1))*nind(pft)
C            if (fuel_1hr(pft,1) .gt. 0.0) then
C            do nc=2,nco2
C               litter_inc=lm_ind(pft,1)*lm_ind(pft,nc)
C     *              +sm_ind(pft,1)*sm_ind(pft,nc)*0.045
C     *              +hm_ind(pft,1)*hm_ind(pft,nc)+0.045
C               fuel_1hr(pft,nc)=(fuel_1hr(pft,nc)*temp
C     *              +litter_inc*nind(pft))/fuel_1hr(pft,1)
C            enddo
C            endif


            temp=fuel_10hr(pft,1)
            fuel_10hr(pft,1)=fuel_10hr(pft,1)+((0.075*
     *           sm_ind(pft,1)+0.075*hm_ind(pft,1))*nind(pft))
            if (fuel_10hr(pft,1) .gt. 0.0) then
            do nc=2,nco2
               litter_inc=sm_ind(pft,1)*sm_ind(pft,nc)*0.075
     *              +hm_ind(pft,1)*hm_ind(pft,nc)+0.075
               fuel_10hr(pft,nc)=(fuel_10hr(pft,nc)*temp
     *              +litter_inc*nind(pft))/fuel_10hr(pft,1)
            enddo
            endif

            temp=fuel_100hr(pft,1)
            fuel_100hr(pft,1)=fuel_100hr(pft,1)+((0.21*
     *           sm_ind(pft,1) +0.21*hm_ind(pft,1))*nind(pft))
            if (fuel_100hr(pft,1) .gt. 0.0) then
            do nc=2,nco2
               litter_inc=sm_ind(pft,1)*sm_ind(pft,nc)*0.21
     *              +hm_ind(pft,1)*hm_ind(pft,nc)+0.21
               fuel_100hr(pft,nc)=(fuel_100hr(pft,nc)*temp
     *              +litter_inc*nind(pft))/fuel_100hr(pft,1)
            enddo
            endif

            temp=fuel_1000hr(pft,1)
            fuel_1000hr(pft,1)=fuel_1000hr(pft,1)+((0.67
     *           *sm_ind(pft,1)+0.67*hm_ind(pft,1))*nind(pft))
            if (fuel_1000hr(pft,1) .gt. 0.0) then
            do nc=2,nco2
               litter_inc=sm_ind(pft,1)*sm_ind(pft,nc)*0.67
     *              +hm_ind(pft,1)*hm_ind(pft,nc)+0.67
               fuel_1000hr(pft,nc)=(fuel_1000hr(pft,nc)*temp
     *              +litter_inc*nind(pft))/fuel_1000hr(pft,1)
            enddo
            endif

              temp=litter_bg(pft,1)
              litter_bg(pft,1)=litter_bg(pft,1)+rm_ind(pft,1)
     *             *(nind_old- nind(pft))
              if (litter_bg(pft,1).gt.0.0) then
                 do nc=2,nco2
                    litter_bg(pft,nc)=(litter_bg(pft,nc)*temp
     *                   +rm_ind(pft,1)*rm_ind(pft,nc)*
     *                   (nind_old-nind(pft)))/litter_bg(pft,1)
                 enddo
              endif
          endif

        enddo

        fpc_total=fpc_total-(fpc_tree_total-0.95)
        if (fpc_total.gt.1.0) fpc_total=1.0
        fpc_tree_total=0.95

      endif     !fpc_tree_total gt 0.95


c     SECTION FOR GRASSES

      do pft=1,npft
    
        if (present(pft).and..not.tree(pft)) then
           if (estab(pft)) then
c          Grasses can establish in non-vegetated areas
              if (ngrass.gt.0) then
                 bare=(1.0-fpc_total)/real(ngrass)
              else
                 bare=0.0
              endif

              bare_max=( (-2.0*crownarea(pft)*alog(1.0-bare-
     *             fpc_grid(pft))/sla(pft))-
     *             lm_ind(pft,1) )/lm_sapl(pft,1)

c         Accumulate biomass increment due to grass establishment

             if (bare.le.0.0) then
                litter_bg(pft,1)=litter_bg(pft,1)-(bare*rm_sapl(pft,1))
                litter_ag_leaf(pft,1)=litter_ag_leaf(pft,1)-
     *             (bare*lm_sapl(pft,1))
                
c            KIRSTEN: per fuel class
              fuel_1hr_leaf(pft,1)=fuel_1hr_leaf(pft,1)-
     *          (bare*lm_sapl(pft,1))
             else
               if(bare.gt.bare_max) bare=bare_max
               temp=lm_ind(pft,1)
               lm_ind(pft,1)=lm_ind(pft,1)+bare*lm_sapl(pft,1)
               if (lm_ind(pft,1).ne.0.0) then
                  do nc=2,nco2
                     lm_ind(pft,nc)=(lm_ind(pft,nc)*temp+bare*
     *                    lm_sapl(pft,1)*lm_sapl(pft,nc))/lm_ind(pft,1)
                  enddo
               endif
               temp=rm_ind(pft,1)
               rm_ind(pft,1)=rm_ind(pft,1)+bare*rm_sapl(pft,1)
               if (rm_ind(pft,1).ne.0.0) then
                  do nc=2,nco2
                     rm_ind(pft,nc)=(rm_ind(pft,nc)*temp+bare*
     *                    rm_sapl(pft,1)*rm_sapl(pft,nc))/rm_ind(pft,1)
                  enddo
               endif     
             endif

             temp=acflux_estab(1)
             all=bare*(lm_sapl(pft,1)+rm_sapl(pft,1))
             if (all*crownarea(pft).gt.eps) then
                acflux_estab(1)=acflux_estab(1)+all*crownarea(pft)
                if (acflux_estab(1).gt.0.0.and.temp.ge.0.) then
                   do nc=2,nco2
                      acflux_estab(nc)=(acflux_estab(nc)*temp+bare*
     *                     (lm_sapl(pft,1)*lm_sapl(pft,nc)+
     *                     rm_sapl(pft,1)*rm_sapl(pft,nc))*
     *                     crownarea(pft))/acflux_estab(1)
                   enddo
                elseif (acflux_estab(1).gt.0.0.and.temp.lt.0.) then
                   do nc=2,nco2
                      acflux_estab(nc)=bare*(lm_sapl(pft,1)
     *                     *lm_sapl(pft,nc)+rm_sapl(pft,1)
     *                     *rm_sapl(pft,nc))/all
                   enddo
                endif
             endif             
           endif

           if (lm_ind(pft,1).le.0.0) then
            present(pft)=.false.
            do m=1,12
               temp=mnpp(m,pft,1)
               mnpp_add(m,1)=mnpp_add(m,1)+mnpp(m,pft,1)
               if (mnpp_add(m,1).gt.0.) then
                  do nc=2,nco2
                     mnpp_add(m,nc)=(mnpp_add(m,nc)*temp+mnpp(m,pft,1)
     *                    *mnpp(m,pft,nc))/mnpp_add(m,1)
                  enddo
               endif
               temp=anpp_add(1)
               anpp_add(1)=anpp_add(1)+mnpp(m,pft,1)
               if (anpp_add(1).gt.0.) then
                  do nc=2,nco2
                     anpp_add(nc)=(anpp_add(nc)*temp+mnpp(m,pft,1)
     *                    *mnpp(m,pft,nc))/anpp_add(1)
                  enddo
               endif
            enddo
            temp=litter_bg(pft,1)
            litter_bg(pft,1)=litter_bg(pft,1)+rm_ind(pft,1)*nind(pft)
            if (litter_bg(pft,1).gt.0.0) then
               do nc=2,nco2
                  litter_bg(pft,nc)=(litter_bg(pft,nc)*temp+
     *                 (rm_ind(pft,1)*rm_ind(pft,nc)*nind(pft)))/
     *                 litter_bg(pft,1)
               enddo
            endif
         endif
       endif
      enddo

c     recalculate fpc's

      do pft=1,npft
         if (present(pft)) then
           if (crownarea(pft).gt.0.0) then
             lai_ind(pft)=(lm_ind(pft,1)*sla(pft))/crownarea(pft)
           else
             lai_ind(pft)=0.0
           endif

           fpc_ind=(1.0-exp(-0.5*lai_ind(pft)))
           fpc_grid(pft)=crownarea(pft)*nind(pft)*fpc_ind
         endif
      enddo


c     Doug 03/09: update changein fuel loads
      DO pft=1,npft
        IF (fuel_1hr_leaf(pft,1)-temp_fuel(1,pft,1)>0) THEN
          DO nc=1,nco2
            fuel_1hr_leaf_inc_pos(pft,nc,:)=
     *         fuel_1hr_leaf_inc_pos(pft,nc,:)+
     *        (fuel_1hr_leaf(pft,1)-temp_fuel(1,pft,1))/12

          END DO
        ELSE
          DO nc=1,nco2
            fuel_1hr_leaf_inc_neg(pft,nc,:)=
     *         fuel_1hr_leaf_inc_neg(pft,nc,:)+
     *        (fuel_1hr_leaf(pft,1)-temp_fuel(1,pft,1))/12
          END DO
        END IF
		
        IF (fuel_1hr_wood(pft,1)-temp_fuel(2,pft,1)>0) THEN
          DO nc=1,nco2
            fuel_1hr_wood_inc_pos(pft,nc,:)=
     *         fuel_1hr_wood_inc_pos(pft,nc,:)+
     *        (fuel_1hr_wood(pft,1)-temp_fuel(2,pft,1))/12

          END DO
        ELSE
          DO nc=1,nco2
            fuel_1hr_wood_inc_neg(pft,nc,:)=
     *         fuel_1hr_wood_inc_neg(pft,nc,:)+
     *        (fuel_1hr_wood(pft,1)-temp_fuel(2,pft,1))/12
          END DO
        END IF

        DO nc=1,nco2
          fuel_10hr_inc(pft,nc,:)=fuel_10hr_inc(pft,nc,:)+
     *      (fuel_10hr(pft,1)-temp_fuel(3,pft,1))/12

          fuel_100hr_inc(pft,nc,:)=fuel_100hr_inc(pft,nc,:)+
     *      (fuel_100hr(pft,1)-temp_fuel(4,pft,1))/12

          fuel_1000hr_inc(pft,nc,:)=fuel_1000hr_inc(pft,nc,:)+
     *      (fuel_1000hr(pft,1)-temp_fuel(5,pft,1))/12
        END DO
      END DO

      return
      end
c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE bt_establish
c Doug 02/13: Calculates new bark thickness spread after establishment
c     SUBROUTINE bt_establish
c Doug 02/13: Calculates new bark thickness spread after establishment
       SUBROUTINE bt_establish(BTparam,BTmode0,nind0,estab_grid)
	 
       IMPLICIT NONE
	   
C      Inputs
       REAL BTmode0
       REAL nind0, estab_grid
	   
C      INPUTS/OUTPUTS
       REAL BTparam(1:3)
	   
C      local vaiables
       REAL A,B,C,D
       !REAL alpha,beta,gamma,delta
       REAL s1,s2
       REAL c1,c2
       REAL d1,d2
       REAL m,n
       REAL BTmean,BTmode,BTmode_old,BTmode_frac
	   
       IF (nind0==0) THEN
         BTparam(2)=BTmode0
         RETURN
       END IF
	   !PRINT*, "EST"
	   !PRINT*, "@@@@@"
	   !PRINT*, BTparam
       s1=BTparam(1)
       s2=BTparam(3)
       c1=BTmode0;
       c2=BTparam(2)
       m=estab_grid
       n=nind0
       BTmode_old=BTparam(2)
	   
	   !PRINT*, "EST"
	   !print*, "!!!!!!!"
	   !print*, s1,s2,c1,c2,m,n
	   
       A=-(2*c2**3-3*s1*c2**2+s1**3)/(3*(s1-s2)*(c2-s1))
       B=(2*c2**3-3*s2*c2**2+s2**3)/(3*(s1-s2)*(c2-s2))
       C=-m*(2*c1**3-3*s1*c1**2+s1**3)/(3*n*(s1-s2)*(c1-s1))
       D=m*(2*c1**3-3*s2*c1**2+s2**3)/(3*n*(s1-s2)*(c1-s2))
	   
       !alpha=-(c2**2-2*c2*s1+s1**2)/((s1-s2)*(c2-s1))
       !beta=(c2**2-2*c2*s2+s2**2)/((s1-s2)*(c2-s2))
       !gamma=-m*(c1**2-2*c1*s1+s1**2)/(n*(s1-s2)*(c1-s1))
       !delta=m*(c1**2-2*c1*s1+s1**2)/(n*(s1-s2)*(c1-s2))
	   
       BTmean=(A+B+C+D)/(1+m/n);	   
       BTmode=3*BTmean-s1-s2; 
	   
       BTmode_frac=(BTmode-s1)/(s2-s1)
       BTparam(2)=BTparam(1)+BTmode_frac*(BTparam(3)-BTparam(1))
	   

	   
	   
       IF (BTparam(2)<BTmode0 .OR. BTparam(2)>BTmode_old) THEN
         IF (abs(BTparam(2)-BTmode0)<0.0001) THEN
           BTparam(2)=BTmode0
         ELSEIF (abs(BTparam(2)-BTmode_old)<0.0001) THEN
           BTparam(2)=BTmode_old
         ELSE
           !WRITE(10,*), "EST BT PROBLEM"
           !WRITE(10,*), "error: BT medium is lower or higher than"
           !WRITE(10,*), "lowest/highest possible limit"
           !WRITE(10,*), "EST BT PROBLEM"
	       !WRITE(10,*), "+_+_+_+_+_+_+"
	       !WRITE(10,*), "A,B,C,D: ", A,B,C,D
	       !WRITE(10,*), "BTmean: ", BTmean
	       !WRITE(10,*), "BTmode: ", BTmode
	       !WRITE(10,*), "BTparam(2): ", BTparam(2)
	       !WRITE(10,*), "BTmode0: ", BTmode0
	       !WRITE(10,*), "BTmode_old: ", BTmode_old
	       !WRITE(10,*), "*/*/*/*/"
           IF (BTparam(2)<BTmode0) BTparam(2)=BTmode0
           IF (BTparam(2)>BTmode_old) BTparam(2)=BTmode_old
         END IF
       END IF
	   
       RETURN
	   
       end
           
	   
c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE decay
c     Calculation of radioactive 14C decay

      subroutine decay(present,lm_ind,sm_ind,hm_ind,rm_ind,
     *     litter_ag_leaf,litter_ag_wood,
     *     litter_bg,cpool_fast,cpool_slow,
     *     fuel_1hr_leaf,fuel_1hr_wood,fuel_10hr,fuel_100hr,fuel_1000hr)

      implicit none

c     PARAMETERS
      integer npft
      parameter (npft=13)
      integer nco2
         parameter (nco2=3)
      real lambda
      parameter (lambda=8267) 

c     ARGUMENTS
      logical present(1:npft)
      real lm_ind(1:npft,1:nco2),rm_ind(1:npft,1:nco2)
      real hm_ind(1:npft,1:nco2),sm_ind(1:npft,1:nco2)
      REAL litter_ag_leaf(1:npft,1:nco2)
      REAL litter_ag_wood(1:npft,1:nco2)
      real litter_bg(1:npft,1:nco2)
      real cpool_fast(1:nco2),cpool_slow(1:nco2)
      REAL fuel_1hr_leaf(1:npft,1:nco2)      !1hr dead fuel: dead grass,"cured" grass,shed tree leaves, small twigs
      REAL fuel_1hr_wood(1:npft,1:nco2)      !1hr dead fuel: dead grass,"cured" grass,shed tree leaves, small twigs
      real fuel_10hr(1:npft,1:nco2)          !10hr dead fuel: large twigs
      real fuel_100hr(1:npft,1:nco2)         !100hr dead fuel: small branches
      real fuel_1000hr(1:npft,1:nco2)        !1000hr dead fuel: logs, bole, large branches
      

c     LOCAL VARIABLES 
      integer i,j,pft,c
      real dfac
      
      i=0
      j=0
      pft=0
      c=3
      dfac=exp(-1./lambda)

      cpool_slow(c)=(((cpool_slow(c)/1000.+1.)*dfac)-1.)*1000.
      cpool_fast(c)=(((cpool_fast(c)/1000.+1.)*dfac)-1.)*1000.

      
      
      do pft=1,npft
         lm_ind(pft,c)=(((lm_ind(pft,c)/1000.+1.)*dfac)-1.)*1000.
         sm_ind(pft,c)=(((sm_ind(pft,c)/1000.+1.)*dfac)-1.)*1000.
         hm_ind(pft,c)=(((hm_ind(pft,c)/1000.+1.)*dfac)-1.)*1000.
         rm_ind(pft,c)=(((rm_ind(pft,c)/1000.+1.)*dfac)-1.)*1000.
         litter_ag_leaf(pft,c)=
     *       (((litter_ag_leaf(pft,c)/1000.+1.)*dfac)-1.)*1000.
         litter_ag_wood(pft,c)=
     *       (((litter_ag_wood(pft,c)/1000.+1.)*dfac)-1.)*1000.
         litter_bg(pft,c)=(((litter_bg(pft,c)/1000.+1.)*dfac)-1.)*1000.
         fuel_1hr_leaf(pft,c)=
     *     (((fuel_1hr_leaf(pft,c)/1000.+1.)*dfac)-1.)*1000.
         fuel_1hr_wood(pft,c)=
     *     (((fuel_1hr_wood(pft,c)/1000.+1.)*dfac)-1.)*1000.
         fuel_10hr(pft,c)=(((fuel_10hr(pft,c)/1000.+1.)*dfac)-1.)*1000.
         fuel_100hr(pft,c)=(((fuel_100hr(pft,c)/1000.+1.)
     *                     *dfac)-1.)*1000.
         fuel_1000hr(pft,c)=(((fuel_1000hr(pft,c)/1000.+1.)
     *                     *dfac)-1.)*1000.
      enddo       
      

      return

      end



***************************************************************************      c
c         subroutine calculates the area of the pixel(0.5x0.5 or other)degress
c         alberte's method

      subroutine pixelarea(lat,area)
      implicit none
      real lat,ilat
      double precision pie,re,dip,equ
c     >     sum,earth
       real area

*       re = 6.378E6
       pie = 4. * ATAN(1.)
       dip = pie/180.

*       area=((111.0*10**3*0.167)**2)*cos(lat*dip)
        area=((111.0*10**3*0.5)**2)*cos(lat*dip)



       return
       end

c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     Stephen & Doug 05/09: land use change, agriculture subroutine
c     clearance = conversion from natural gridcell - agricultural 
c     abandonment= reset the state variables, start natural regeneration
c                  from bare ground.
c
      subroutine agriculture(clear,mtemp_min20,tree,agri_burn, 
     *  lm_ind,sm_ind,hm_ind,rm_ind,nind,convflux,agri_litter_ag,
     *  agri_litter_bg,litter_ag,litter_bg,abandon,prod10,prod100,
     *  fpc_grid,present)

      implicit none  

c     PARAMETERS
      INTEGER npft
          PARAMETER(npft=13)
      
c     ARGUMENTS
      logical clear,tree(1:npft),abandon,present(1:npft)
      real mtemp_min20,agri_burn 
      real lm_ind(1:npft),sm_ind(1:npft)
      real hm_ind(1:npft),rm_ind(1:npft),nind(1:npft) 
      real convflux,litter_ag(1:npft),litter_bg(1:npft)
      real prod10(0:10),prod100(0:100)
      real fpc_grid(1:npft)
      real agri_litter_ag,agri_litter_bg        

c     LOCAL VARIABLES 
      real above
      integer pft
      

c     intialisations
      above=0.0   
      agri_burn=0.0
      prod10(0)=0.0
      prod100(0)=0.0
      convflux=0.0      

c     Grid cell cleared in this year for agriculture
      if (clear) then
         if (mtemp_min20.gt.15.5) then    ! tropical ecosystem
            do pft=1,npft   
               if (.not.tree(pft)) then  ! tropical grass
                 agri_burn=agri_burn+(lm_ind(pft)*nind(pft))+
     *                      litter_ag(pft)
                 convflux=convflux+(lm_ind(pft)*nind(pft))+
     *                      litter_ag(pft)
                 agri_litter_bg=agri_litter_bg+litter_bg(pft)
     *                           + (rm_ind(pft)*nind(pft))
               else  ! woody pfts
                 above=(lm_ind(pft)+sm_ind(pft)+hm_ind(pft))*nind(pft)
                 agri_burn=agri_burn+ ((40.0/67.0)*above)+litter_ag(pft)
                 convflux=convflux+ ((40.0/67.0)*above)+litter_ag(pft)
                 prod10(0)=prod10(0)+(27.0/67.0)*above
                 agri_litter_bg=agri_litter_bg+litter_bg(pft)
     *               +(rm_ind(pft)*nind(pft))
               endif
             enddo 
             agri_litter_ag=0.0

           else  ! temperate or boreal
             do pft=1,npft
               if (.not.tree(pft)) then  ! grass
                 convflux=convflux+(lm_ind(pft)*nind(pft))+
     *                    litter_ag(pft)   ! fodder or burnt  
                 agri_litter_bg=agri_litter_bg+litter_bg(pft)
     *                          + (rm_ind(pft)*nind(pft))
               else ! woody pfts
                 above=(lm_ind(pft)+sm_ind(pft)+hm_ind(pft))*nind(pft)
                 convflux=convflux+(40.0/67.0)*above
                 prod10(0)=prod10(0)+(20.0/67.0)*above
                 prod100(0)=prod100(0)+(7.0/67.0)*above
                 agri_litter_ag=agri_litter_ag+litter_ag(pft)
                 agri_litter_bg=agri_litter_bg+litter_bg(pft)
     *               +(rm_ind(pft)*nind(pft))
               endif
             enddo
           endif

       endif  ! gridcell clearance

c      abandonment of previous agricultural land to natural 
c      vegetation
         
       if (abandon) then
c
c         set all state variables to zero: build up from bare ground
c                                                                   
         do pft=1,npft
            nind(pft)=0.0
            lm_ind(pft)=0.0
            sm_ind(pft)=0.0
            hm_ind(pft)=0.0 
            rm_ind(pft)=0.0
            litter_ag(pft)=0.0
            litter_bg(pft)=0.0
            fpc_grid(pft)=0.0
            present(pft)=.false.
         enddo
       endif   

       return
       end

c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     Stephen & Doug 05/09: land use change,
c     agricultural production and product decomposition
c                                                    
      subroutine agriprod(present,agri,anpp,cflux_prod_total,
     *   prod10_total,prod100_total,prod1,anpp_agri,rap,x10,
     *   x100,agri_litter_ag,agri_litter_bg,prod10,prod100)

      implicit none
 
c     PARAMETERS
      real above,below
        parameter (above=0.4,below=0.6)
      INTEGER npft
        PARAMETER(npft=13)

c     ARGUMENTS
      logical present(1:npft),agri
      real anpp(1:npft)
      real cflux_prod_total
      real prod10_total,prod100_total,prod1
      real anpp_agri,rap
      real x10(0:10),x100(0:100)
      real agri_litter_ag,agri_litter_bg
      real prod10(0:10),prod100(0:100)

c     LOCAL VARIABLES 
      integer i,j,pft       
      real cflux_prod1,cflux_prod10,cflux_prod100 
      real anpp_total
     

      cflux_prod_total=0.0
      cflux_prod1=0.0
      cflux_prod10=0.0 
      cflux_prod100=0.0
      prod10_total=0.0
      prod100_total=0.0
      anpp_total=0.0
      anpp_agri=0.0 
       
      do pft=1,npft
        if (present(pft)) anpp_total=anpp_total+anpp(pft)
      enddo
                               
c     calculate agricultural npp
      if (agri) anpp_agri=rap*anpp_total  

c     update the product pools: as a linear decay of the 
c     initial carbon inputs, e.g. the annual release from
c     the prod10 pool represents 10% of the initial carbon
c     entering the prod10 pool during the previous 10 years
c                                          
c     x10 represents 10% of the initial conversion pool with
c         10 year turnover. 
c     x100 represents 1% of the initial conversion pool with
c         100 year turnover 
c 
c     Each year the pool decomposes by x10, x100 respectively.
c
      do i=0,8            ! check
       j=10-i
       cflux_prod10=cflux_prod10+x10(j-1)
       prod10(j)=prod10(j-1)-x10(j-1)
       prod10_total=prod10_total+prod10(j) 
       x10(j)=x10(j-1)
       if (prod10(j).lt.1.0) prod10(j)=0.0
      enddo  

      x10(1)=0.1*prod10(0)
      prod10(1)=prod10(0)
      prod10_total=prod10_total+prod10(1)         
      
      do i=0,98
       j=100-i
       cflux_prod100=cflux_prod100+x100(j-1)
       prod100(j)=prod100(j-1)-x100(j-1)
       prod100_total=prod100_total+prod100(j)
       x100(j)=x100(j-1)
       if (prod100(j).lt.1.0) prod100(j)=0.0
      enddo
      x100(1)=0.01*prod100(0)
      prod100(1)=prod100(0)
      prod100_total=prod100_total+prod100(1)

c
c     assume harvest is consumed in one year
c     the crop is partitioned into 40% above ground biomass
c                                  60% below ground biomass
c     assume above ground biomass goes into the prod1 pool.
c
      prod1=above*anpp_agri
      cflux_prod1=prod1
      agri_litter_bg=agri_litter_bg+(below*anpp_agri)

c
c     product flux
c                 
      cflux_prod_total=cflux_prod1+cflux_prod10+cflux_prod100

      return
      end
c -----------------------------------------------------------------------------
c                                 REFERENCES
c -----------------------------------------------------------------------------

c Carslaw, HS & Jaeger JC 1959 Conduction of Heat in Solids, Oxford University
c   Press, London
c Collatz, GJ, Ball, JT, Grivet C & Berry, JA 1991 Physiological and
c   environmental regulation of stomatal conductance, photosynthesis and
c   transpiration: a model that includes a laminar boundary layer. Agricultural
c   and Forest Meteorology 54: 107-136
c Collatz, GJ, Ribas-Carbo, M & Berry, JA 1992 Coupled photosynthesis-stomatal
c   conductance models for leaves of C4 plants. Australian Journal of Plant
c   Physiology 19: 519-538
c Farquhar GD & von Caemmerer 1982 Modelling of photosynthetic response to
c   environmental conditions. In: Lange, OL, Nobel PS, Osmond CB, Ziegler H
c   (eds) Physiological Plant Ecology II: Water Relations and Carbon
c   Assimilation, Vol 12B. Springer, Berlin, pp 549-587.
c Foley J A 1995 An equilibrium model of the terrestrial carbon budget
c   Tellus (1995), 47B, 310-319
c Harper JL 1977 Population Biology of Plants, Academic Press, London
c Haxeltine A & Prentice IC 1996 BIOME3: an equilibrium terrestrial biosphere
c   model based on ecophysiological constraints, resource availability, and
c   competition among plant functional types. Global Biogeochemical Cycles 10:
c    693-709
c Haxeltine A & Prentice IC 1996a A general model for the light-use efficiency
c   of primary production. Functional Ecology 10: 551-561
c Henderson-Sellers, A & Robinson, PJ 1986 Contemporary Climatology. Longman,
c   Essex.
c Jarvis, PG & McNaughton KG 1986 Stomatal control of transpiration: scaling up
c   from leaf to region. Advances in Ecological Research 15: 1-49
c Jury WA, Gardner WR & Gardner WH 1991 Soil Physics 5th ed, John Wiley, NY
c Larcher W 1983 Physiological Plant Ecology, 2nd ed, Springer-Verlag, Berlin
c Lloyd, J & Taylor JA 1994 On the temperature dependence of soil respiration
c   Functional Ecology 8: 315-323
c Monsi, M & Saeki, T 1953 Ueber den Lichtfaktor in den Pflanzengesellschaften
c   und seine Bedeutung fuer die Stoffproduktion. Japanese Journal of Botany
c   14: 22-52
c Monteith, JL & Unsworth, MH 1990 Principles of Environmental Physics, 2nd ed,
c   Arnold, London
c Prentice, IC, Sykes, MT & Cramer W 1993 A simulation model for the transient
c   effects of climate change on forest landscapes. Ecological Modelling 65:
c   51-70.
c Press, WH, Teukolsky, SA, Vetterling, WT & Flannery, BT. 1986. Numerical
c   Recipes in FORTRAN, 2nd ed. Cambridge University Press, Cambridge
c Reich, PB, Walters, MB & Ellsworth, DS, 1997. From tropics to tundra: global
c   convergence in plant functioning. Proceedings of the National Academy of
c   Sciences USA 94: 13730-13734.
c Ryan, GR 1991 Effects of climate change on plant respiration. Ecological
c   applications 1: 157-167
c Shinozaki, K, Yoda, K, Hozumi, K & Kira, T 1964 A quantitative analysis of
c   plant form - the pipe model theory. I. basic analyses. Japanese Journal of
c   Ecology 14: 97-105
c Shinozaki, K, Yoda, K, Hozumi, K & Kira, T 1964 A quantitative analysis of
c   plant form - the pipe model theory. II. further evidence of the theory and
c   its application in forest ecology. Japanese Journal of Ecology 14: 133-139
c Sprugel, DG, Ryan, MG, Brooks, JR, Vogt, KA, Martin, TA, Respiration from the
c   organ level to the stand (in press '93 ---> look up)
c van Duin, RHA 1963 The influence of soil management on the temperature
c   wave near the surface. Tech Bull 29 Inst for Land and Water Management
c   Research, Wageningen, Netherlands
c Waring, RH Schroeder, PE & Oren, R 1982 Application of the pipe model theory
c   to predict canopy leaf area. Canadian Journal of Forest Research 12:
c   556-560

c -----------------------------------------------------------------------------
