If you are a member of the Yao group, please the Odyssey specific install instructions (much simpler) in the examples folder in /n/home06/jkemp/SharedDMT .

# Install Instructions #

First, you must install iTensor. Download it here:

	git clone https://github.com/ITensor/ITensor

Now follow the install instructions in the INSTALL file (alternatively you can read them online [here](http://www.itensor.org/docs.cgi?page=install&vers=cppv3)). Notice you need a C++ compiler which supports C++17, so if, for example, you have an ancient LTS version of Ubuntu, you may need to update gcc (using apt-get). You will also need to have development LAPACK and BLAS libraries. (Ulimately I recommend [Intel MKL] (https://software.intel.com/content/www/us/en/develop/documentation/get-started-with-mkl-for-dpcpp/top.html).

Then clone this repository:

	git clone https://github.com/Jack-Kemp/dmt.git

Now copy the examples folder into a new directory for your first project (rename it to something helpful!):

	cp -r <THIS_DIR>/examples <YOUR-DIR>
	
	cd <YOUR-DIR>
	
Edit the Makefile in your new directory so that the following lines point to the correct location:

	LIBRARY_DIR=<YOUR-ITENSOR-INSTALL>
	DMT_DIR=<YOUR-DMT-INSTALL>

Finally attempt to compile:

	make
	
Now if it makes successfully you can run using:

	export OMP_NUM_THREADS=1 # or number of threads desired.

	./xyz xyz_example_input_params.txt
	
If you haven't changed the input parameters, the output should go into ./out. The python notebook Analysis.ipynb has an example of how you can read in and plot this output data.


# Compile Options #


If you wish to compile your own code, use 
	
	make APP=your_name

or change the hardcoded APP variable in the Makefile.
	
	make debug
	
will compile the debug version of the code, which has a "-g" appended to it.

You can then use gdb to debug (or run under valgrind, heaptrack etc.).

# Examples #

xyz.cc represents the most basic usage of DMT. It runs DMT for the XYZ spin chain (nearest neighbour XX + YY + ZZ interactions) and a small number of observables.

xyzAdvanced.cc adds use of checkpointing, long range interactions, and more complex observables such as local energy density. All of the features of the DMT code currently are featured in xyzAdvanced.cc.

xyzladderRung.cc and xyzLadderSwap.cc both implement an XYZ ladder; i.e. a pseudo 2-dimensional systems. xyzladderRung.cc does this by combining two spins on a rung of a ladder into a single site, as shown in RungBasis.h. xyzladderSwap.cc uses the inbuilt long-range interactions of the DMT code, which implements swap gates. What this means is that gates are still applied between two neighbouring sites, but swap gates are used to move non-neighbouring sites together.

# Checkpointing #

Writing Checkpoint = True and CheckpointTime = N will cause the code to write the output and the density matrix to a checkpoint file every N hours (use the code in xyzAdvanced.cc as an example for how to setup checkpointing for your simulation). By default these files will be named "OutputName\_t\_X.(dat/dmt)" where X is the simulation time the snapshot was taken. Using the .dmt file, one can then restart the simulation at the same point you left off, by simply using the same input parameter file (the start time will change automatically). This is useful if a simulation runs out of time on the cluster (it also works with the automatic restart that serial_requeue sometimes performs), or if you want to extend the simulation time by using the same input parameter file but increasing tTotal.

When a simulation is restarted from a checkpoint like this, the output files will begin from at t=tstart, not t=0. To concatenate all of the output files, possibly from multiple restarts, use the script checkpoint\_glob.sh in examples (you may need to change the sys.path at the top of the script to point to your DMT install location) like this:

    python checkpoint_cat.py example_filename_t_*.dat

Which will output a file example\_filename\_t\_*.dat\_globbed with the entire continous time series. Use checkpoint\_cat\_all.py to do the same for an entire directory of different checkpoint files.










	
	
	
	
	
	
