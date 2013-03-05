# To run the model, you will have to have the following installed:
#   1) fortran compiler - to compile the main, science bit of the model. At Macquarie, we use gfortran
#   2) c++ compiler - to compile the driver (I/O) bit. We use GCC c++
#   3) Netcdf for c++. Specifically, -lnetcdf_c++ -lnetcdf libraries
#   4) A little file that helps fortran and c++ talk to each other. e.g lstdc++
# 
# Before compiling, you need to alter the "makefile" compilers and paths. Look for and alter the following lines. If you are running on the holocene, lgm or eocene servers at Bristol, you might be able to get away with skipping this bit:
#   1) NETCDF = <netcdf binary folder>>
#   2) FC = <your fortran compiler>
#   3) FCOPTIONS = <relevant fortran options, including option for pre-compiling>. For gfortran, these are typically -x -f77 -cpp-input
#   4) CPPOPTIONS: <you c++ libraries>
#   5) cxx = <you c++ compiler>
#   6) cxxoptions = <relevet c++ options. Typically -o3>
#   7) You will finally need to change the 2nd path listed on the line "LDLIBS" to point to your talky c++ Fortran file
#  
# Once you've done this, its time to compile. LPX has 3 compiling options. Type either:
#   1) > make one_step : makes a file you can run in one go. Useful for equilibrium climate runs.
#   2) > make two_step: makes two running options: a spin-up; and a run. Spin-up is important for dynamic runs, as you need to get the plants agrowing and various carbon pools filled before you can do a proper run. This is the option most people should use
#   3) > make three_step: makes three options: a spin-up; a ramp-up; and a run. This is useful for anyone doing futures runs. spin up gets the plants and carbon pools going, ramp-up can e.g. run the model to the present day, where all data is saved, and you can start each future run from one point, rather then running the historic for each future run.
#  
#   There is also an option call "make all"., but I dont think you'll need this unless your doing some funny development stuff
#  
#  
# Now, before running the model, you need your climate inputs. The typical historical CRY inouts most of us use are on terrafirma server at Macquarie Uni. Server address is sceince.terrafirma.mq.edu.au. Hopefully at some point, we will make some inputs public on www.bcd.mq.edu.au , but for now, you'll ned a terrafirma login to get data. EMail Rhys  on rhys.whitely@mq.edu.au for more info.
# For most runs, you will need spin-up (typically detrended) data for the spin up, and then real data fro the actual run.
# Data required is listed in the config file, and is likley to change with new developments. As of Prentice et al (2011), see input helpfile for me info.
# 
# Set each path in the config file (lpj.cfg) to the relevent inputs (examples given in the file already). For the spin-up, this shiuld be your detrended data. Once the spin-up is complete, this should be your real data. There are some other opitions in the cfg file that I'll write about later.
# Run the model with the follwing commands:
# 
#   1) > motif-lpj : runs the model in 1 step, useful for equlibrium runs
#   2) > motif-lpj-step1a : runs the spinup component, and ouputs the model state at the end of the spinup into the file listed after "SPINUP_FILE:" in lpj.cfg
#   3) > motif-lpj-step1b : runs the rampup component. Reads state of model from SPINUP_FILE: <file> , runs through rampup, and ouputs state of model into RAMPUp_FILE: <file>
#   4) > motif-lpj-step2  : runs the final step of the model, and ouputs defined inputs into netcdf file defined by OUTPUT_FILE: <file>
#   
# And thats how to setup and run the model.
# 
# Goodluck!

Doug (douglas.kelley@mq.edu.au)
