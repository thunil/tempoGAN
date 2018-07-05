Here's a quick guide how to run this example:

First generate simulation data with the following command. Note that the
parameters in this scene will lead to large velocities in the UI; this is
correct, and no sign of a solver explosion, or so.
	>>> manta manta_flip.py

The first step only generates simulation data which is not directly usable for the NN. 
the training data for tensorflow is prepared with:
	>>> manta manta_gendata.py

Then use it to train a first model with:
	>>>python tf_train.py ../data/manta-flip/training_data/ 
the last parameter contains the extracted data from the previous step, you can also enter multiple ones here.

Once the model is trained, you can use it in a new sim with
	>>> manta manta_mlflip.py

Note that all data is stored in ../data/ by default.  You can pass a directory
to load from in manta_mlflip with the "--load" option.

