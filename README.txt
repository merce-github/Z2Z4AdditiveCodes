             The "Z2Z4-Additive Codes" Package
             ---------------------------------

"Z2Z4-Additive Codes"  is  a  package  that  supports the basic
facilities  to  work  with Z2Z4-additive codes (see [1],[2]) in
Magma (see [3]). With  this  package, this kind of codes can be
created and manipulated and  information about  these codes can
be calculated.

"Z2Z4-Additive Codes"  consists  of five  files  written in the 
Magma language. Please send your  bug reports to Combinatorics, 
Coding and  Security Group (CCSG)  at support.ccsg@uab.cat 
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

The   "Z2Z4-Additive  Codes"   package  is  composed  of  four
directories:

/src: The files to attach to Magma “Z2Z4AdditiveCodes.m”,
      "ZZZ4StandardForm.m", "Z2Z4MinimumWeight.m", "Z2Z4Cyclic.m",
      "Z2Z4CoveringRadius.m", "Z2Z4CodeConstructions" and
      "Z2Z4Decode.m".
/docs: The manual to use the package, in pdf format.
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

1. Unpack the archive file in a directory $DIR.

2. If you do not have a directory to store user-defined packages,
   create one in your preferred location, for instance $HOME:
	
     mkdir UserPackages

   If you already have such a directory, proceed with the instructions
   changing $HOME/UserPackages by the path to your directory.

3. Create a new directory in the $HOME/UserPackages directory:

     cd UserPackages
     mkdir Z2Z4AdditiveCodes

4. Copy all files in $DIR/src/ into this new directory Z2Z4AdditiveCodes:

     cp $DIR/src/* $HOME/UserPackages/Z2Z4AdditiveCodes/

5. Create a file named spec in the directory $HOME/UserPackages:

     touch spec

   Edit the spec file and add the following content: 

     Z2Z4AdditiveCodes
     {
        +Z2Z4AdditiveCodes.spec
     }

   In case that the spec file already exists, add the lines above at
   the end of the old spec file.

6. Ensure that all files have the correct permissions:

     chmod -R a+rX .

7. Set the environment variable MAGMA_USER_SPEC to the spec file.
   Change to the directory where Magma is installed and edit the magma script.
   Locate the line

     export MAGMA_SYSTEM_SPEC

   and add the following lines just after that:

     MAGMA_USER_SPEC="$HOME/UserPackages/spec"
     export MAGMA_USER_SPEC

8. In order to check that the package has been installed correctly,
   run magma in a terminal window and try to run the following lines:
   
     Z2Z4AdditiveCodes_version();
     Z2Z4CodeConstructions_version();
     Z2Z4CoveringRadius_version();
     Z2Z4Cyclic_version();
     Z2Z4Decode_version();
     Z2Z4MinimumWeight_version();
     Z2Z4StandardForm_version();

   If the installation has been successful, magma should return the
   following lines, one for each function, respectively:

     [4,6]
     [1,3]
     [1,5]
     [1,5]
     [2,2]
     [2,4]
     [2,1]
  
   If the numbers appear but are different from the ones shown above,
   then the respective files have not been installed correctly and
   they may correspond to a previous version.


                             Bug reports
                             -----------

When  sending a  bug  report to support.ccsg@uab.cat or to
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


September 29, 2022
