# Install Instructions #

Copy this folder to your desired directory (rename it to something helpful!), then source the required modules and compile.

	cp /n/home06/jkemp/SharedDMT/examples <YOUR-DIR>
	
	cd <YOUR-DIR>

	source load_modules.sh
	
	make
	
Now if it makes successfully you can run using (on a test node!):

	./xyz xyz_example_input_params.txt
	
The output should go into ./out. The python notebook Analysis.ipynb has an example of how you can read in and plot this output data.


# Submitting to the cluster #

The python script submit.py can sumbit to Odyssey. Usage:

	python submit.py ./xyz xyz_example_input_params.txt
	
submit_batch.py can submit a batch of jobs to the cluster, sweeping over the value of an input variable in a closed interval, for example:

	python submit_batch.py ./xyz xyz_example_input_params.txt Jx 6 0 1.0 2

will submit for 2 jobs handling cases Jx = (0, 0.2, 0.4) and (0.6, 0.8, 1.0) respectively.

To change memory, time allowance, queue etc. just edit the scripts directly.


# Compile Options #


If you wish to compile your own code, use 
	
	make APP=your_name

or change the hardcoded APP variables in the Makefile.
	
	make debug
	
will compile the debug version of the code, which has a "-g" appended to it.

You can then use gdb to debug (or run under valgrind, etc.).






	
	
	
	
	
	
