Maya plugins for loading mesh surfaces and density grids
--------------------------------------------------------------

These plugins provide Maya DAG objects which automatically load a 
frame from .bobj.gz or .uni file sequence, corresponding to the current animation frame.
This means that only the current frame has to be held in memory.
The building process has been tested with Maya 2010-2012 under Windows and Linux,
but should also work with different versions.

Installation
-------------

1. Linux:

> make && make install
This should compile and install the modules in the user's local maya directory.
If maya is installed in a nonstandard directory, you might need to set the MAYA_LOCATION
environment variable or edit buildconfig.linux
This will build the plugins and copy them in the user plugin directory

2. Windows

For Maya2012 x64, you can use the prebuilt .mll files in this directory. For other setups,
you will need to build them using the provided MSVC projects (you will need zlib to compile).
Copy The .mll files to a plugin path (e.g. "My Documents\maya\plug-ins"), and the .mel files
to a script path (e.g. "My Documents\maya\scripts"). 


Usage
------

1. In Maya, you can activate the plugins in Windows->Settings->Plugin Manager.
   If the modules are installed correctly, 'densityloader' and 'bobjloader' should appear.
   Check the load, autoload boxes.
2. To use the objects in a Maya scene, input the following MEL commands:
   > source createBobjLoader; 
   or
   > source createDensityLoader; 
  
3. If you open the attribute editor, you should now be editing the attributes of
   bobjFluidObjectNode1 or gridFluidObjectNode1. 
   Next, you'll need to specify the input file mask. This 
   is the filename of the .bobj.gz-files, containing one printf-like decimal
   placeholder for the frame number (e.g. "C:\data\surface%04d.bobj.gz"). Note
   that the attribute editor script will automatically replace the string "0000"
   with "%04d", so that if you pick the first file of a sequence in the file
   browser, you shouldn't have to make any changes to the file mask.

4. The index offset attribute will be added to the current Maya animation frame
   number to determine the actual file index. The default is -1, because
   .bobj.gz-sequences may be 0-based and Maya frames are 1-based.
   