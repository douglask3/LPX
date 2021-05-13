# LPX

LPX is a fire enabled  DGVM in the LPJ (Sitch et al. 2003) suite of models. It is a development of LPJ-SPITFIRE (Thonicke et al., 2010), with some improvement in the fire component, documented in Prentice et al (2011). 

The configuration here (in lpj.cfg) is from Sato et al. 2021. LPX  was also used to simulate current (Murray et al., 2012) & future  (Murray et al., 2013) global hydrological cycles; carbon stores at the last ice age (Ciais et al., 2012); and in the development of a DGVM benchmarking system (Kelley et al., 2013).

## How to run
 To run the model, you will have to have the following installed
 
1.  fortran compiler - to compile the main, science bit of the model. e.g gfortran
2. c++ compiler - to compile the driver (I/O) bit. We use GCC c++
3. Netcdf for c++. Specifically, -lnetcdf_c++ -lnetcdf libraries
4. A little file that helps fortran and c++ talk to each other, called something along the lines of
 lstdc++

### Makefile
Before compiling, you need to alter the "makefile" compilers and paths. Look for and alter the following lines:
```
NETCDF = <netcdf binary folder>
FC = <your fortran compiler>
FCOPTIONS = <relevant fortran options, including option for pre-compiling. For gfortran, these 
        are typically -x -f77 -cpp-input>
CPPOPTIONS= <you c++ libraries>
cxx = <you c++ compiler>
cxxoptions = <relevant  c++ options> Typically -o3
you will finally need to change the 2nd path listed on the line "LDLIBS" to point to your talky 
c++ Fortran file (lstdc++)
```
### Compile
Once you've done this, it’s time to compile. LPX has 3 compiling options. Type either:

1. ```> make one_step``` makes a file for running the model in one go. Useful for equilibrium climate runs.
2. ```> make two_step``` makes two running options: a spin-up; and a run. Spin-up is important for  
         dynamic runs, as you need to get the plants growing and various carbon pools filled before  
         you can do a proper run. This is the option most people should use
3. ```> make three_step ``` makes three options: a spin-up; a ramp-up; and a run. This is useful for 
        anyone doing futures runs, or for anyone unsure about how long their spin-up should take. 
        Spin up gets the plants and carbon pools going, ramp-up can either continue the spin-up (i.e, if 
        the model is not in equilibrium – see spinup protocol on wiki); or run the model dynamically to a 
        common point (e.g. to the present day), where all data is saved, and you can start new runs 
        (e.g. each future run) from that one point. 
  
There is also an option call "make all", but I don’t think you'll need this unless you’re doing some funny development stuff. This makes all the options from all the other make commands.
 
### Inputs
Before running the model, you need your climate inputs. The historical CRU inputs most of us 
use are on terrafirma server at Macquarie Uni. Server address is sceince.terrafirma.mq.edu.au. 
Hopefully at some point, we will make some inputs public on www.bcd.mq.edu.au , but for now, 
you'll need a terrafirma login to get data. 
LPX takes inputs in netcdf format. If you want to provide you own, then visit our wiki for 
information about what you need.


 For most runs, you will need spin-up (typically detrended) data for the spin up, and then real (or 
 transient) data from the actual run. See configuration help file for full list of input data.
 
 ### Configure
 The configuration file (lpj.cfg) lists all the information to set up the model. To set the model 
 properly, change each file path in the config file (paths ending in .nc) to the relevant inputs 
 (examples given in the file already):
 
* ```TEMP:``` monthly gridded mean temperature over the period of the run
* ```PREC:``` monthly gridded precipitation over the period of the run
* ```WETDAYS:``` monthly gridded number of wetdays over the period of the run
* ```SUN:``` monthly gridded cloud over the period of the run
* ```TMIN:``` monthly gridded mean minimum temperature over the period of the run
* ```TMAX:``` monthly gridded mean maximum temperature over the period of the run
* ```WINDSPEED:``` monthly gridded windspeed over the period of the run
* ```MASK:``` A landmask (no time dimension) defining which cells should be simulated. All other 
		inputs must be on the same sized grid as this (i.e, if this is 0.5 degree, so should
		everything else). If there is data in the other files that is not in the defined mask from here
		(i.e, if any of the other inputs have ocean values), they will simply be ignored.
* ```MASK_VAR:``` name of the variable in the landmask file that corresponds to the mask.
* ```SOIL_TYPE:``` gridded soil type (no time dimension). Soils should be listed as of Sitch et al 
		(2003).
* ```SOIL_TYPE_VAR:``` name of variable in soiltype file corresponding to the soil.
* ```POPDENS:``` annual gridded population density. Not currently used by the model, but we 
		haven’t taken it off the input list cos we might use it again one day
* ```Crop:``` annual gridded cropland fraction in each cell
* ```Pasture:``` annual gridded paster for each cell. Only used on some runs, but still there’s a line 
		for it in the file.
* ```Lightn:``` gridded lighting climatology (no time dimension…. yet).
* ```A_ND:``` no longer used, but still listed as an input. Maybe I’ll find time to remove this at some 
		point. Annual, gridded input of whatever you want, cos it makes no difference.
* ```CO2_FILE:``` Annual (not gridded) atmospheric CO2 concentration
* ```C13_FILE:``` Annual (not gridded) atmospheric delta C13 values.
* ```C14_NORTH_FILE:``` Annual (not gridded) atmospheric delta C14 values for the northern 
		hemisphere.
* ```C14_EQUATOR_FILE:``` Annual (not gridded) atmospheric delta C14 values for the equator
* ```C14_SOUTH_FILE:``` Annual (not gridded) atmospheric delta C14 values for the southern 
		hemisphere.
* ```CO2_VAR:``` CO2 variable in CO2 file
* ```C13_VAR:``` c13 variable in delta C13 file
* ```C14_VAR:``` c14 variable in delat C14 file

 For the spin-up, these line should point to your detrended data. Once the spin-up is complete, this should be 
 your real data. 

 Then, set the following:
* ``` SUN_IS_CLOUD:``` ```TRUE``` If you are using cloud data as inputs, else ```FALSE```
* ``` SPINUP_YEARS:``` Set to the number of years you wish to do a spin-up. Different applications 
		have different protocols for the length of the spin-up. There will be some listed on 
		the wiki  at some point.
* ``` RUN_YEARS:``` The number of years the model needs to run past spin-up. E.g, for spin ups 
		ending in 1850, set to 150 to run the model onto 2000.
* ``` RUN_OUT_YEARS:``` The year from which you want LPX to start giving you outputs. Ie., staying 
		with the same example, if you want LPX to give you outputs from 1991, set to 141.
* ``` RUN_OUT_FREQ:``` How often you want the model to output. The model averages the 
		variables over this time. i.e, if you want the average output over the period 1991-
		2000, set to 10.
* ``` SPINUP_FILE:``` This is the file where the state of the model after the spin-up will be stored. 
* ``` RAMPUP_FILE:``` This is the file where the state of the model after the ramp-up will be stored.
* ``` OUTPUT_FILE:``` This is the netcdf file where defined outputs will be stored.
* ``` OUTPUT_VAR:``` A variable you’d like to output. You will need an OUPUT_VAR line for each variable 
       your outputting. You’ll have to look in the lpjmain code and search through the available 
       variables for now. At some point, I will put an OUPUT key together.

### Run
Once this is all set up, you’re ready to run. Run the model with the following commands:

1. ```> motif-lpj :``` obtained by the make one_step command. Runs the model in 1 step, useful for 
        equilibrium runs
2. ```> motif-lpj-step1a :```  obtained by the make two_step & make three_step commands. Runs the 
        spinup component and outputs the model state at the end of the spin-up into the file listed 
        after ```SPINUP_FILE:``` in lpj.cfg
3, ```> motif-lpj-step1b : <file?``` obtained by the make three_step command Runs the rampup component. 
        Reads state of model from ```SPINUP_FILE: <file>``` , runs through rampup, and outputs state of 
        model into ```RAMPUP_FILE: <file>```
4.  ```> motif-lpj-step2  :``` obtained by the make two_step & make three_step commands. Runs the 
        final step of the model. Read state of model form SPINUP_FILE: <file> and outputs defined 
        inputs into netcdf file defined by ```OUTPUT_FILE: <file>```.
   
 The model outputs the defined variables (listed after ```OUTPUT_VAR``` in the config file) into a netcdf 
 file. The output file will be named after the ```OUTPUT_FILE: <file>``` line, and will have the year the 
 output related to after it. This year is in model years, so if you run a 7000 year spin up to pre-
 industrial 1850, and you want to find your output for the year 1990, it will be referred to as 7000 
 (number of spin-up years) +140 (number of years after pre-industrial) = 7140. So the output will be 
 called <name>-7140.nc . If you have asked for an output freqancy of e.g. 10 years, the output file in 
 this example, 7140, will cover 1981-1990 (i.e, the previus 9 years, plus the year its outputting on).
 And thats how to setup and run the model.
 
 Goodluck!

 Doug (douglas.i.kelley@gmail.com)





## references
Ciais P, Tagliabue A, Cuntz M, Bopp L, Scholze M, Hoffmann G, Lourantou A, Harrison SP, Prentice IC, Kelley DI, Koven C. Large inert carbon pool in the terrestrial biosphere during the Last Glacial Maximum. Nature Geoscience. 2012 Jan;5(1):74-9.

Kelley DI, Prentice IC, Harrison SP, Wang H, Simard M, Fisher JB, Willis KO. A comprehensive benchmarking system for evaluating global vegetation models. Biogeosciences. 2013 May 17;10(5):3313-40.

Murray SJ, Watson IM, Prentice IC. The use of dynamic global vegetation models for simulating hydrology and the potential integration of satellite observations. Progress in Physical Geography. 2013 Feb;37(1):63-97.

Murray SJ, Watson IM, Prentice IC. The use of dynamic global vegetation models for simulating hydrology and the potential integration of satellite observations. Progress in Physical Geography. 2013 Feb;37(1):63-97.

Prentice IC, Kelley DI, Foster PN, Friedlingstein P, Harrison SP, Bartlein PJ. Modeling fire and the terrestrial carbon balance. Global Biogeochemical Cycles. 2011 Sep;25(3).
 
Sato H, Kelley DI, Mayor S, Calvo MM, Cowling S, Prentice IC, Amazonian Dry Corridors Opened by Fire and CO2 Deprivation during the Last Glacial Maximum. 2021, accepted in Nature Geoscience 

Sitch S, Smith B, Prentice IC, Arneth A, Bondeau A, Cramer W, Kaplan JO, Levis S, Lucht W, Sykes MT, Thonicke K. Evaluation of ecosystem dynamics, plant geography and terrestrial carbon cycling in the LPJ dynamic global vegetation model. Global change biology. 2003 Feb;9(2):161-85.

Thonicke K, Spessa A, Prentice IC, Harrison SP, Dong L, Carmona-Moreno C. The influence of vegetation, fire spread and fire behaviour on biomass burning and trace gas emissions: results from a process-based model. Biogeosciences. 2010 Jun 23;7(6):1991-2011.
