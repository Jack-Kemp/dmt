# Install Instructions #

Copy this folder to your desired directory (rename it to something helpful!), then source the required modules and compile.

	cp -r /n/home06/jkemp/SharedDMT/examples <YOUR-DIR>
	
	cd <YOUR-DIR>

	source load_modules.sh
	
	make
	
Now if it makes successfully you can run using (on a test node if you can get it to work -- they fail to run anything at all for me):

	export OMP_NUM_THREADS=1 # or number of threads desired.

	./xyz xyz_example_input_params.txt
	
If you haven't changed the input parameters, the output should go into ./out. The python notebook Analysis.ipynb has an example of how you can read in and plot this output data.


# Submitting to the cluster #

The python script submit.py can submit to Odyssey. Usage:

	python submit.py ./xyz xyz_example_input_params.txt
	
submit_batch.py can submit a batch of jobs to the cluster, sweeping over the value of an input variable in a closed interval, for example:

	python submit_batch.py ./xyz xyz_example_input_params.txt Jx 6 0 1.0 2

will submit for 2 jobs handling cases Jx = (0, 0.2, 0.4) and (0.6, 0.8, 1.0) respectively.

To change memory, time allowance, queue etc. just edit the scripts directly, making sure also to change

	export OMP_NUM_THREADS=
	
also.


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

When a simulation is restarted from a checkpoint like this, the output files will begin from at t=tstart, not t=0. To concatenate all of the output files, possibly from multiple restarts, use the script checkpoint\_glob.sh in examples/util like this:

    python checkpoint_cat.py example_filename_t_*.dat

Which will output a file example\_filename\_t\_*.dat\_globbed with the entire continous time series. Use checkpoint\_cat\_all.py to do the same for an entire directory of different checkpoint files.










	
	
	
	
	
	
