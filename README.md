# SingleMicrotubuleBuckling
This repository contains the source code for simulation of buckling of a single microtubule under the force of dynein molecular motors.

Compiling the code:
$ make all
You can change the compiler in the Makefile. Please check the Makefile before compiling the code.

Output of the code:
Once you run
$make all
the Makefile will create an executable named mt.
The output of mt will be the x,y,z coordinates of each element of the microtubule at each timepoint of the simulation. These coordinates are best viewed in VMD, which can be downloaded from here:
http://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD
VMD needs to know the connectivity between differnt segments to show them correctly. This information is stored in a.psf file, which can be viewed and edited in a text editor. It shows the number of segments and the connectivity between them.

Big picture information:
The unit.xlsx (and unit.ods) file contains the units of all physical quantities employed in the model.
The main.c file employs slightly different units that the rod.c file. However, all unit conversions are taken care of within the codes.
the unit of distance is 1 micrometer in main.c, while it is 1 nanometer in rod.c.
the unit of time is 1 millisecond in main.c, while it is 1 microsecond in rod.c.
the unit of force is 1 femtonewton in main.c, it is 1 femtonewton in rod.c as well.
