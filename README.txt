             The "Z2Z4-Additive Codes" Package
             ---------------------------------

"Z2Z4-Additive Codes"  is  a  package  that  supports the basic
facilities  to  work  with Z2Z4-additive codes (see [1],[2]) in
Magma (see [3]). With  this  package, this kind of codes can be
created and manipulated and  information about  these codes can
be calculated.

"Z2Z4-Additive Codes"  consists  of five  files  written in the 
Magma language. Please send your  bug reports to Combinatorics, 
Coding and  Security Group (CCSG)  at support-ccsg@deic.uab.cat 
or if it is a Magma problem to magma-trouble (magma@maths.usyd.
edu.au). See the section below.

"Z2Z4-Additive Codes" was originally written in Magma by Bernat
Gastón  (supervised by Mercè Villanueva and Jaume Pujol)  as  a
support for a research project on Z2Z4-additive codes developed
by  the  Combinatorics, Coding and Security Group (CCSG) within
the  Department of  Information  and Communications Engineering
(dEIC) at the Autonomous University of Barcelona (UAB).

This version of the package has been developped in Magma version
2.26-8.


                    Composition of the package
                        ------------------

The   "Z2Z4-Additive  Codes"   package  is  composed  of  five
directories:

/src: The files to attach to Magma “Z2Z4AdditiveCodes.m”,
      "ZZZ4StandardForm.m", "Z2Z4MinimumWeight.m", "Z2Z4Cyclic.m", 
      "Z2Z4CoveringRadius.m", "Z2Z4CodeConstructions" and
      "Z2Z4Decode.m".
/license: The license of the package.
/doc: The manual to use the package in pdf format.
/examples: Examples  from  the  manual.  They can  be loaded in
           Magma as soon as the package is attached.
/test: Some test files that can be used to check the package.

            Using/Installing "Z2Z4-Additive Codes"
             ---------------------------------

To use  "Z2Z4-Additive Codes"  temporally  (as a Magma Package)
unpack  the  archive  file in a directory.   Enter to the ./src
directory. Call Magma and then write:

   AttachSpec(“Z2Z4AdditiveCodes.spec”);

in the Magma command line.

To install "Z2Z4-Additive Codes" permanent (as a Magma Package):

1. Unpack the archive file in a directory.

2. Enter  to  the  directory  where  Magma  is installed, go to
   package directory    $PATHMAGMA/package/    and create a new
   directory.

     mkdir Z2Z4AdditiveCodes

3. Copy $PATH/src/Z2Z4AdditiveCodes.spec $PATH/src/Z2Z4AdditiveCodes.m 
   $PATH/src/Z2Z4StandardForm.m
   $PATH/src/Z2Z4MinimumWeight.m $PATH/src/Z2Z4Cyclic.m,
   $PATH/src/Z2Z4CoveringRadius.m, $PATH/src/Z2Z4CodeConstructions.m
   and $PATH/src/Z2Z4Decode.m
   to the Magma directory $PATHMAGMA/package/Z2Z4AdditiveCodes/

   cp $PATH/src/Z2Z4AdditiveCodes.spec
			  $PATHMAGMA/package/Z2Z4AdditiveCodes/

   cp $PATH/src/Z2Z4AdditiveCodes.m
			  $PATHMAGMA/package/Z2Z4AdditiveCodes/
			
   cp $PATH/src/Z2Z4StandardForm.m
			  $PATHMAGMA/package/Z2Z4AdditiveCodes/

   cp $PATH/src/Z2Z4MinimumWeight.m
			  $PATHMAGMA/package/Z2Z4AdditiveCodes/

   cp $PATH/src/Z2Z4Cyclic.m
			  $PATHMAGMA/package/Z2Z4AdditiveCodes/

   cp $PATH/src/Z2Z4CoveringRadius.m
			  $PATHMAGMA/package/Z2Z4AdditiveCodes/

   cp $PATH/src/Z2Z4CodeConstructions.m
			  $PATHMAGMA/package/Z2Z4AdditiveCodes/

   cp $PATH/src/Z2Z4Decode.m
			  $PATHMAGMA/package/Z2Z4AdditiveCodes/

4. Edit the file      $PATHMAGMA/package/spec     and write the
   following lines at the end:

      Z2Z4AdditiveCodes
      {
         +Z2Z4AdditiveCodes.spec
      }


                             Bug reports
                             -----------

When  sending a  bug  report to support-ccsg@deic.uab.cat or to
magma@maths.usyd.edu.au,    remember we will need to be able to
reproduce the problem; so please include:

 * The  version  of  Magma  you  are  using; look at the header 
   when you start up Magma.
 * The  operating  system you are using e.g. Linux, SunOS 5.8 =
   Solaris 2.8, IRIX 6.5, Windows, ...
 * A script that demonstrates the bug, along with a description
   of why it's a bug (e.g.  by  adding comments to  the  script
   - recall  comments  in Magma  begin  with  a  //  or between
   /*  */).


                             Bibliography
                               --------

[1]  J. Borges, C. Fernández, J. Pujol, J. Rifà and M. Villanueva,
   "Z2Z4-linear codes: generator matrices  and  duality", Designs,
   Codes and Cryptography, vol 54, no. 2, pp. 167-179, 2010.

[2] C. Fernández, J. Pujol, and M. Villanueva, "Z2Z4-linear codes:
   rank and kernel", Designs, Codes and Cryptography,  vol 56, no.
   1, pp. 43-59, 2010.

[3] J.J. Cannon and W. Bosma (Eds.) "Handbook of Magma Functions",
   Edition 2.13, 4350 pages, 2006.


July 28, 2022