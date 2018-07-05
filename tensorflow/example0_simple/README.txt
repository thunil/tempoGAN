This is an "as simple as possible" example using mantaflow data in tensorflow.
It's a simple 3 layer network resembling a primitive auto-encoder (i.e., it
tries to restore the input from a small hidden inner layer, the latent space).

This example has no parameters. Simply run the generation script a few times
(ca. 10 times) to generate ../data/simSimple_1xxx directories with data for
training:
	>>> manta manta_genSimSimple.py

Then call the training script to run tensorflow:
	>>> python tf_simple.py 

This will load all possible uni files, use 10% as validation data, and output a
sequence of re-generated fields as out_XXX.png. Note that the result most
likely will be quite flickery due to the small latent space and the sub-optimal
network.

(Try replacing the fully connected layers with convolutional ones to improve
the results - or check out the next example.)

