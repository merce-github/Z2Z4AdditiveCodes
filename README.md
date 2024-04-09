# The "Z2Z4-Additive Codes" Package

"Z2Z4-Additive Codes"  is  a  package  that  supports the basic
facilities  to  work  with Z2Z4-additive codes (see [[1]](#1), [[2]](#2)) in
Magma (see [[3]](#3)). With  this  package, this kind of codes can be
created and manipulated and  information about  these codes can
be calculated.

"Z2Z4-Additive Codes"  consists  of five  files  written in the 
Magma language. Please send your  bug reports to Combinatorics, 
Coding and  Security Group (CCSG)  at [support.ccsg@uab.cat](mailto:support.ccsg@uab.cat)
or if it is a Magma problem to magma-trouble (magma@maths.usyd.edu.au). See section [Bug reports](#bug-reports) below.

"Z2Z4-Additive Codes" was originally written in Magma by Bernat
Gastón  (supervised by Mercè Villanueva and Jaume Pujol)  as  a
support for a research project on Z2Z4-additive codes developed
by  the  Combinatorics, Coding and Security Group (CCSG) within
the  Department of  Information  and Communications Engineering
(dEIC) at the Autonomous University of Barcelona (UAB).

This version of the package has been developed in Magma version
2.26-8.


## Composition of the package

The   "Z2Z4-Additive  Codes"   package  is  composed  of  four
directories:

* [/src](src): The files to attach to Magma “Z2Z4AdditiveCodes.m”,
      "ZZZ4StandardForm.m", "Z2Z4MinimumWeight.m", "Z2Z4Cyclic.m", 
      "Z2Z4CoveringRadius.m", "Z2Z4CodeConstructions" and
      "Z2Z4Decode.m".
* [/docs](docs): The manual to use the package, in pdf format.
* [/examples](examples): Examples  from  the  manual.  They can  be loaded in
           Magma as soon as the package is attached.
* [/test](test): Some test files that can be used to check the package.


## Using/Installing "Z2Z4-Additive Codes"

To use  "Z2Z4-Additive Codes"  temporally  (as a Magma Package)
unpack  the  archive  file in a directory.   Enter to the ./src
directory. Call Magma and then write:
```
   AttachSpec(“Z2Z4AdditiveCodes.spec”);
```
in the Magma command line.

To install "Z2Z4-Additive Codes" permanent (as a Magma Package) on Linux:

1. Unpack the archive file in a directory <code>$DIR</code>.

2. If you do not have a directory to store user-defined packages, create one in your preferred location, for instance <code>$HOME</code>:

   ```
      mkdir UserPackages
   ```

   If you already have such a directory, proceed with the instructions changing <code>$HOME/UserPackages</code> by the path to your directory.

3. Create a new directory in the <code>$HOME/UserPackages</code> directory:

   ```
      cd UserPackages
      mkdir Z2Z4AdditiveCodes
   ```

4. Copy all files in <code>$DIR/src/</code> into this new directory <code>Z2Z4AdditiveCodes</code>:

   ```
      cp $DIR/src/* $HOME/UserPackages/Z2Z4AdditiveCodes/
   ```

5. Create a file named <code>spec</code> in the directory <code>$HOME/UserPackages</code>:

   ```
      touch spec
   ```

   Edit the <code>spec</code> file and add the following content:

   ```
      Z2Z4AdditiveCodes
      {
         +Z2Z4AdditiveCodes.spec
      }
   ```

   In case that the <code>spec</code> file already exists, add the lines above at the end of the old <code>spec</code> file.

6. Ensure that all files have the correct permissions:

   ```
      chmod -R a+rX .
   ```

7. Set the environment variable <code>MAGMA_USER_SPEC</code> to the <code>spec</code> file. Change to the directory where Magma is installed and edit the <code>magma</code> script. Locate the line

   ```
      export MAGMA_SYSTEM_SPEC
   ```

   and add the following lines just after that:

   ```
      MAGMA_USER_SPEC="$HOME/UserPackages/spec"
      export MAGMA_USER_SPEC
   ```

To do the installation on Windows OS, follow the above items from 1 to 5 (in item 4, use "copy" instead of "cp"; and in item 5, use "type nul >>" instead of "touch" if a spec file does not exist). Finally, edit the system environment variables by adding MAGMA_USER_SPEC and set it to the spec file.

In order to check that the package has been installed correctly, run Magma in a terminal window and try to run the following lines:

```
   Z2Z4AdditiveCodes_version();
   Z2Z4CodeConstructions_version();
   Z2Z4CoveringRadius_version();
   Z2Z4Cyclic_version();
   Z2Z4Decode_version();
   Z2Z4MinimumWeight_version();
   Z2Z4StandardForm_version();
```

If the installation has been successful, Magma should return the following lines, one for each function, respectively:

```
   [4,6]
   [1,3]
   [1,5]
   [1,5]
   [2,2]
   [2,4]
   [2,1]
```

If the numbers appear but are different from the ones shown above, then the respective files have not been installed correctly and they may correspond to a previous version.

## Bug reports

When  sending a  bug  report to [support.ccsg@uab.cat](mailto:support.ccsg@uab.cat) or to
magma@maths.usyd.edu.au,    remember we will need to be able to
reproduce the problem; so please include:

 * The  version  of  Magma  you  are  using; look at the
   header when you start up Magma.
 * The  operating  system you are using e.g. Linux, SunOS 5.8 =
   Solaris 2.8, IRIX 6.5, Windows, ...
 * A script that demonstrates the bug, along with a description
   of why it's a bug (e.g.  by  adding comments to  the  script
   _-_ recall  comments  in Magma  begin  with  a  //  or between
   /*  */).


## Bibliography

<a id="1">[1]</a>  J. Borges, C. Fernández, J. Pujol, J. Rifà and M. Villanueva,
   "Z2Z4-linear codes: generator matrices  and  duality", Designs,
   Codes and Cryptography, vol 54, no. 2, pp. 167-179, 2010.

<a id="2">[2]</a> C. Fernández, J. Pujol, and M. Villanueva, "Z2Z4-linear codes:
   rank and kernel", Designs, Codes and Cryptography,  vol 56, no.
   1, pp. 43-59, 2010.

<a id="3">[3]</a> J.J. Cannon and W. Bosma (Eds.) "Handbook of Magma Functions",
   Edition 2.13, 4350 pages, 2006.


October 4, 2022
