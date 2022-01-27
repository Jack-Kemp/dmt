# Summary #

A library implementing Density Matrix Truncation (DMT) for quantum dynamics simulations in one dimension, based on the work of [Chris White *et al.*](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.97.035127), using the [ITensor](http://www.itensor.org) library. The algorithm can be thought of as a variant of the standard time-evolving block decimation (TEBD) for matrix product state wavefunctions, but applied to the density matrix instead. The truncation method is chosen to preserve local observables, making it particularly well-suited for calculating hyrdrodynamic properties such as transport coefficients.

Features:

  * Arbitrary one-dimensional Hamiltonians over arbitrary Hilbert spaces with arbitrary initial condition density matrices can be easily defined and time-evolved.
  * Long-range interactions are supported via SWAP gates.
  * Automatic checkpointing and output parsing for easy analysis in python Jupyter notebooks.

# Install Instructions #
If you are a member of the Yao group, please the Odyssey specific install instructions (much simpler) in the examples folder in /n/home06/jkemp/SharedDMT .

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

N.B. Notice the Makefile has a section for headers files:

HEADERS=RungBasis.h

If you delete RungBasis.h from your own folder, or add your own header files, you **must** change this section fo the Makefile accordingly.

# Examples #

xyz.cc represents the most basic usage of DMT. It runs DMT for the XYZ spin chain (nearest neighbour XX + YY + ZZ interactions) and a small number of observables.

xyzAdvanced.cc adds use of checkpointing, long range interactions, and more complex observables such as local energy density. All of the features of the DMT code currently are featured in xyzAdvanced.cc.

xyzladderRung.cc and xyzLadderSwap.cc both implement an XYZ ladder; i.e. a pseudo 2-dimensional systems. xyzladderRung.cc does this by combining two spins on a rung of a ladder into a single site, as shown in RungBasis.h. xyzladderSwap.cc uses the inbuilt long-range interactions of the DMT code, which implements swap gates. What this means is that gates are still applied between two neighbouring sites, but swap gates are used to move non-neighbouring sites together.

# Checkpointing #

Writing Checkpoint = True and CheckpointTime = N will cause the code to write the output and the density matrix to a checkpoint file every N hours (use the code in xyzAdvanced.cc as an example for how to setup checkpointing for your simulation). By default these files will be named "OutputName\_t\_X.(dat/dmt)" where X is the simulation time the snapshot was taken. Using the .dmt file, one can then restart the simulation at the same point you left off, by simply using the same input parameter file (the start time will change automatically). This is useful if a simulation runs out of time on the cluster (it also works with the automatic restart that serial_requeue sometimes performs), or if you want to extend the simulation time by using the same input parameter file but increasing tTotal.

When a simulation is restarted from a checkpoint like this, the output files will begin from at t=tstart, not t=0. To concatenate all of the output files, possibly from multiple restarts, use the script checkpoint\_glob.sh in examples (you may need to change the sys.path at the top of the script to point to your DMT install location) like this:

    python checkpoint_cat.py example_filename_t_*.dat

Which will output a file example\_filename\_t\_*.dat\_globbed with the entire continous time series. Use checkpoint\_cat\_all.py to do the same for an entire directory of different checkpoint files.

# Input Parameters #

## Required Parameters ##

	tStep
	tTotal
	nSweeps
The time step of each evolution, the total evolution time, and the number of update, DMRG-style sweeps per time step that the Hamiltonian gates are applied to the density matrix.

	MaxDim
The maximum bond dimension for the MPDO representing the density matrix.

	OutputName
The root name for the output files. `OutputName.dat` contains measurements of observables over the evolution, while `OutputName.dmt` gives the final output density matrix itself.

## Basic Optional Parameters ##

	OutputDir = ./
	Verbose = true
	Checkpoint = false
	CheckpointTime = 10.0
	CheckpointName = OutputName
	tStart = 0

The directory for output files, verbosity of output, and checkpointing options, see explanation of checkpointing above. tStart labels the start time of the evolution and is mainly used for checkpointing.

	PresRadius = 1
	
The 'preservation radius' of the DMT's truncation scheme: the number of sites away from the bond being updated in one direction for which expectation values of local observables remain unchanged by the truncation. Because the MPDO must represent the local density matrix within the preservation radius exactly at each truncation, this limits maximum bond dimension which can be truncated, and thus larger preservation radius requires larger bond dimension.

	Cutoff = 1e-16
	AbsoluteCutoff = true

The cutoff for Schmidt values in the truncation step. If AbsoluteCutoff = true, any values smaller than cutoff are truncated. Because this is set to machine precision by default, the DMT algorithm will immediately saturate to MaxDim. Alternatively, one could set MaxDim very large and control the accuracy instead via Cutoff. Setting AbsoluteCutoff = false, Cutoff instead controls the truncation error, the sum of total Schmidt values truncated.

	WriteX
	
Whether to measure and write to the output .dat file observable 'X' as the evolution proceeds.

## Advanced Optional Parameters ##

	SVDMethod = "gesdd"
	
Which underlying SVDMethod to tell ITensor to use when truncating. "gesdd" is the fastest and most stable but occasionally fails due to a bug in all LAPACK implementations, "gesvd" is more stable but slower, and "ITensor" is significantly slower and less accurate but will not fail. "automatic" will attempt to perform "gesdd" and then "gesvd".


	FirstSVDCutoff = 1e-16
	ThirdSVDCutoff = 1e-10
	AbsolutePresCutoff = true

There are actually three SVDs in each DMT trucation step, the main truncating SVD and two additional SVDs sandwiching the actual truncation which are required before and after changing basis to the preferred truncation basis for DMT. These cutoff parameters control these SVDs like the Cutoff parameter controls the basic SVD -- it is recommended you do not change these settings.

	Vectorize = true
	HermitianBasis = true
	
If Vectorize = true then the density matrix MPDO, which has two physical legs per site, is flattened into an MPS, with one physical leg per site. The advantage of this is that the HermitianBasis can then be chosen, which is the basis in which the density matrix is real, which is computationally advantageous. However, conceptually, it can be advantageous to leave the density matrix in MPDO form, in which case Vectorize = false may be chosen.

	DoNormalize = true
	CacheTrace = true
	
Require density matrix normalization and cache partial traces of the density matrix intelligently to speed up calculations. Recommended that both settings are left unchanged except for debugging purposes.

## Experimental Unsupported Parameters ##

	ConserveQNs = false
	OnlyPreserveEnergyDensity = false





	







	









	
	
	
	
	
	
