/*************************************************************/
/*                                                           */
/* Package or project name: Z2Z4AdditiveCodes package        */
/* Test file name: Z2Z4NumInvariants_BB_test.m               */
/*                                                           */
/* Comments: Black-box tests for the intrinsic functions     */
/*           Length, BinaryLength, Z2Z4Type                  */
/*           Ngens, '#', InformationRate                     */
/*           included in the Z2Z4AdditiveCodes.m file        */
/*                                                           */
/* Authors: M. Villanueva                                    */
/*                                                           */
/* Revision version and last date: v1.0, 2016/02/19          */
/*                                 v1.1, 2018/01/15          */
/*         user defined type       v1.2  2019/01/30          */
/*                                                           */
/*************************************************************/

//needs Z2Z4AdditiveCode file

SetAssertions(true);
Alarm(30*60);

/*************************************************************
	GLOBAL VARIABLES
*************************************************************/	
Z4 := Integers(4);

/****************************************************************/
/*                                                              */
/* Function name: Length                                        */
/* Parameters: C                                                */
/* Function description: Given a Z2Z4-additive code C of type   */
/*   (alpha, beta; gamma, delta; kappa), return the length      */
/*   n = alpha + beta, and the sequence [alpha, beta] with the  */
/*   number of coordinates over Z2 and the number of            */
/*   coordinates over Z4, respectively.                         */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - The length of the Z2Z4-additive code C                   */
/*   - The sequence [alpha, beta]                               */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> RngIntElt, SeqEnum              */
/*                                                              */
/****************************************************************/
/****************************************************************/
/*                                                              */
/* Function name: BinaryLength                                  */
/* Parameters: C                                                */
/* Function description: Given a Z2Z4-additive code C of type   */
/*   (alpha, beta; gamma, delta; kappa), return the binary      */
/*   length nbin = alpha + 2 * beta of C, which corresponds to  */
/*   the length of the binary code Cbin = Phi(C), where Phi is  */
/*   the Gray map.                                              */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - The length of the the binary code Cbin = Phi(C)          */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> RngIntElt                       */
/*                                                              */
/****************************************************************/
/****************************************************************/
/*                                                              */
/* Function name: Ngens                                         */
/* Parameters: C                                                */
/* Function description: The number of generators (which equals */
/*   the pseudo-dimension k) of the Z2Z4-additive code C as a   */
/*   quaternary linear code.                                    */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - The number of the generators of the code C               */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> RngIntElt                       */
/*                                                              */
/****************************************************************/
/****************************************************************/
/*                                                              */
/* Function name: Z2Z4Type                                      */
/* Parameters: C                                                */
/* Function description: Given a Z2Z4-additive code C, return a */
/*   sequence with the parameters [alpha,beta,gamma,delta,kappa],*/
/*   that is, the type of the code.                             */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - The parameters [alpha, beta, gamma, delta, kappa]        */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> RngIntElt                       */
/*                                                              */
/****************************************************************/
/****************************************************************/
/*                                                              */
/* Function name: #                                             */
/* Parameters: C                                                */
/* Function description: Given a Z2Z4-additive code C of type   */
/*   (alpha, beta; gamma, delta; kappa), return the number of   */
/*   codewords belonging to C, that is 2^gamma 4^delta.         */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - The number of codewords of C                             */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> RngIntElt                       */
/*                                                              */
/****************************************************************/
print "BASIC NUMERICAL INVARIANTS";

print "test 1: Trivial Z2Z4-additive zero code
               alpha = 2, beta = 4, gamma = 0, delta = 0, kappa = 0, length = 6, #C = 1";
C := Z2Z4AdditiveZeroCode(2, 4);

expectedOutputLength := 6;
expectedOutputBinaryLength := 10;
expectedOutputNgens := 0;
expectedOutputType := [2,4,0,0,0];
expectedOutputCardinal := 1;
expectedOutputInformationRate := 0;

assert Length(C) eq expectedOutputLength;
assert BinaryLength(C) eq expectedOutputBinaryLength;
assert Ngens(C) eq expectedOutputNgens;
assert Z2Z4Type(C) eq expectedOutputType;
assert #(C) eq expectedOutputCardinal;
assert InformationRate(C) eq expectedOutputInformationRate;

/************************************************************/
print "test 2: Trivial Z2Z4-additive universe code
               alpha = 2, beta = 4, gamma = 2, delta = 4, kappa = 2, length = 6, #C = 1024";
C := Z2Z4AdditiveUniverseCode(2, 4);

expectedOutputLength := 6;
expectedOutputBinaryLength := 10;
expectedOutputNgens := 6;
expectedOutputType := [2,4,2,4,2];
expectedOutputCardinal := 1024;
expectedOutputInformationRate := 1;

assert  Length(C) eq expectedOutputLength;
assert  BinaryLength(C) eq expectedOutputBinaryLength;
assert  Ngens(C) eq expectedOutputNgens;
assert Z2Z4Type(C) eq expectedOutputType;
assert #(C) eq expectedOutputCardinal;
assert  InformationRate(C) eq expectedOutputInformationRate;

/************************************************************/
print "test 3: Repetition Z2Z4-additive code
               alpha = 4, beta = 8, gamma = 1, delta = 0, kappa = 1, length = 12, #C = 2";
C := Z2Z4AdditiveRepetitionCode(4, 8);

expectedOutputLength := 12;
expectedOutputBinaryLength := 20;
expectedOutputNgens := 1;
expectedOutputType := [4,8,1,0,1];
expectedOutputCardinal := 2;
expectedOutputInformationRate := 0.05;

assert  Length(C) eq expectedOutputLength;
assert  BinaryLength(C) eq expectedOutputBinaryLength;
assert  Ngens(C) eq expectedOutputNgens;
assert Z2Z4Type(C) eq expectedOutputType;
assert #(C) eq expectedOutputCardinal;
assert  InformationRate(C) eq expectedOutputInformationRate;

/************************************************************/
print "test 4: A Z2Z4-additive code with alpha=0
               alpha = 0, beta = 19, gamma = 7, delta = 2, kappa = 0, length = 19, #C = 2048";
R := RSpace(Z4, 19);
L := [R![0,0,2,0,0,0,2,0,2,0,0,0,0,0,2,2,2,2,0],
      R![0,0,0,0,0,0,2,0,2,0,2,0,2,0,2,2,2,2,0],
      R![0,0,0,0,2,2,0,0,2,0,0,0,0,2,0,2,2,0,0],
      R![0,2,0,0,0,0,2,0,2,0,2,0,0,2,2,2,0,2,0],
      R![0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,2,0,0],
      R![0,0,0,0,2,0,2,2,2,2,2,0,0,0,0,0,2,2,0],
      R![0,0,0,0,2,0,0,0,0,0,2,0,0,0,0,2,0,2,2],
      R![0,0,0,1,2,0,2,1,0,3,0,1,1,0,2,0,2,0,1],
      R![1,1,1,0,0,1,2,0,2,3,0,1,1,2,0,2,2,3,1]];
C := Z2Z4AdditiveCode(L, 0);

expectedOutputLength := 19;
expectedOutputBinaryLength := 38;
expectedOutputNgens := 9;
expectedOutputType := [0,19,7,2,0];
expectedOutputCardinal := 2048;
expectedOutputInformationRate := 11/38;

assert  Length(C) eq expectedOutputLength;
assert  BinaryLength(C) eq expectedOutputBinaryLength;
assert  Ngens(C) eq expectedOutputNgens;
assert Z2Z4Type(C) eq expectedOutputType;
assert #(C) eq expectedOutputCardinal;
assert  InformationRate(C) eq expectedOutputInformationRate;

/************************************************************/
print "test 5: A Z2Z4-additive code with beta=0
               alpha = 10, beta = 0, gamma = 4, delta = 0, kappa = 4, length = 10, #C = 16";
R := RSpace(Z4, 10);
L := [R![0,0,1,0,1,0,1,0,1,0],
      R![1,0,0,0,0,0,1,1,1,1],
      R![0,0,1,0,0,0,0,1,0,0],
      R![0,1,0,1,0,1,0,1,1,1]];
C := Z2Z4AdditiveCode(L, 10);

expectedOutputLength := 10;
expectedOutputBinaryLength := 10;
expectedOutputNgens := 4;
expectedOutputType := [10,0,4,0,4];
expectedOutputCardinal := 16;
expectedOutputInformationRate := 0.4;

assert  Length(C) eq expectedOutputLength;
assert  BinaryLength(C) eq expectedOutputBinaryLength;
assert  Ngens(C) eq expectedOutputNgens;
assert Z2Z4Type(C) eq expectedOutputType;
assert #(C) eq expectedOutputCardinal;
assert  InformationRate(C) eq expectedOutputInformationRate;

/************************************************************/
print "test 6: A Z2Z4-additive code with gamma=0
               alpha = 10, beta = 19, gamma = 0, delta = 7, kappa = 0, length = 19, #C = 16384";
R := RSpace(Z4, 29);
M := Matrix(Z4,[[2,2,0,0,2,2,0,2,2,2,0,0,0,1,0,0,0,0,2,3,2,1,1,3,3,1,2,3,0],
                [0,2,0,0,0,0,2,0,0,0,0,0,0,0,0,2,0,1,0,3,0,2,2,1,3,1,1,0,1],
                [0,0,2,2,0,0,2,2,2,2,0,1,0,0,0,0,0,0,2,3,1,1,0,2,1,0,0,0,1],
                [2,0,0,0,2,0,2,0,0,0,0,0,0,0,0,0,1,0,3,3,1,3,1,0,0,3,1,0,2],
                [2,2,2,2,2,0,0,0,0,0,1,0,0,0,0,1,0,0,0,2,3,0,0,1,1,0,2,2,1],
                [2,0,2,2,2,0,2,2,0,2,0,0,1,0,0,3,0,0,0,1,2,3,1,3,0,2,1,3,2],
                [0,2,2,0,2,0,2,2,0,0,0,0,0,0,1,1,0,0,0,3,3,2,0,2,0,3,1,1,3]]);
L := [R!v : v in RowSequence(M)];
C := Z2Z4AdditiveCode(L, 10);

expectedOutputLength := 29;
expectedOutputBinaryLength := 48;
expectedOutputNgens := 14;
expectedOutputType := [10,19,0,7,0];
expectedOutputCardinal := 16384;
expectedOutputInformationRate := 14.0/48;

assert  Length(C) eq expectedOutputLength;
assert  BinaryLength(C) eq expectedOutputBinaryLength;
assert  Ngens(C) eq expectedOutputNgens;
assert Z2Z4Type(C) eq expectedOutputType;
assert #(C) eq expectedOutputCardinal;
assert  InformationRate(C) eq expectedOutputInformationRate;

/************************************************************/
print "test 7: A Z2Z4-additive code with delta=0
               alpha = 10, beta = 19, gamma = 7, delta = 0, kappa = 4, length = 19, #C = 128";
R := RSpace(Z4, 29);
M := Matrix(Z4,[[2,0,0,0,0,0,2,2,2,0,0,0,0,0,0,2,2,0,0,0,0,0,0,2,0,2,0,2,2],
                [0,2,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,2,2,2,2],
                [0,0,2,0,0,2,2,0,2,0,0,0,0,0,2,0,2,2,0,2,0,0,0,2,2,2,2,0,0],
                [0,0,0,2,0,2,0,0,2,0,0,0,0,0,2,0,0,0,2,0,2,2,0,0,2,0,2,0,2],
                [0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,0,2,0,2,0,0,2,0,2,0,0,0,0,2],
                [0,0,0,0,0,0,0,0,0,0,0,2,0,2,0,2,2,2,0,2,0,2,0,2,2,2,2,0,2],
                [0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,2,2,0,0,0,0,0,2,2,0,2,2]]);
L := [R!v : v in RowSequence(M)];
C := Z2Z4AdditiveCode(L , 10);

expectedOutputLength := 29;
expectedOutputBinaryLength := 48;
expectedOutputNgens := 7;
expectedOutputType := [10,19,7,0,4];
expectedOutputCardinal := 128;
expectedOutputInformationRate := 7.0/48;

assert  Length(C) eq expectedOutputLength;
assert  BinaryLength(C) eq expectedOutputBinaryLength;
assert  Ngens(C) eq expectedOutputNgens;
assert Z2Z4Type(C) eq expectedOutputType;
assert #(C) eq expectedOutputCardinal;
assert  InformationRate(C) eq expectedOutputInformationRate;					

/************************************************************/
print "test 8: A Z2Z4-additive code with kappa=0
               alpha = 10, beta = 19, gamma = 7, delta = 2, kappa = 0, length = 29, #C = 2048";
R := RSpace(Z4, 29);
M := Matrix(Z4,[[0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,0,2,0,0,0,0,0,2,2,2,2,0],
                [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,2,0,2,0,2,0,2,2,2,2,0],
                [0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,0,0,2,0,0,0,0,2,0,2,2,0,0],
                [0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,2,0,2,0,2,0,0,2,2,2,0,2,0],
                [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,2,0,0],
                [0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,2,2,2,2,2,0,0,0,0,0,2,2,0],
                [0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2,0,0,0,0,2,0,2,2],
                [0,0,2,0,0,0,2,0,0,0,0,0,0,1,2,0,2,1,0,3,0,1,1,0,2,0,2,0,1],
                [0,2,0,2,0,2,0,2,2,2,1,1,1,0,0,1,2,0,2,3,0,1,1,2,0,2,2,3,1]]);
L := [R!v : v in RowSequence(M)];
C := Z2Z4AdditiveCode(L, 10);

expectedOutputLength := 29;
expectedOutputBinaryLength := 48;
expectedOutputNgens := 11;
expectedOutputType := [10,19,7,2,0];
expectedOutputCardinal := 2048;
expectedOutputInformationRate := 11.0/48;

assert  Length(C) eq expectedOutputLength;
assert  BinaryLength(C) eq expectedOutputBinaryLength;
assert  Ngens(C) eq expectedOutputNgens;
assert Z2Z4Type(C) eq expectedOutputType;
assert #(C) eq expectedOutputCardinal;
assert  InformationRate(C) eq expectedOutputInformationRate;

/************************************************************/
print "test 9: A Z2Z4-additive code with kappa=alpha
               alpha = 7, beta = 12, gamma = 8, delta = 2, kappa = 7, length = 19, #C = 4096";
R := RSpace(Z4, 19);
M := Matrix(Z4,[[0,0,0,2,0,2,2,0,0,0,2,2,0,0,2,2,0,0,0],
                [0,0,0,2,0,0,0,0,0,0,2,0,0,0,0,2,2,0,0],
                [0,0,2,0,2,2,2,0,0,0,2,2,0,0,0,0,0,0,0],
                [0,0,0,2,0,0,0,0,0,2,2,2,0,0,0,2,0,0,0],
                [2,0,0,0,2,2,0,0,0,0,2,0,0,0,0,0,0,0,0],
                [0,2,0,0,2,2,0,0,0,0,0,2,0,2,0,0,0,0,0],
                [0,2,0,2,2,0,0,0,2,0,2,2,0,0,0,0,0,0,2],
                [0,2,0,2,0,0,0,0,0,0,2,2,0,0,0,0,0,2,2],
                [0,2,0,0,0,0,2,0,0,0,0,0,1,1,1,2,1,1,2],
                [0,2,0,2,0,0,2,1,1,1,2,0,0,1,0,0,1,0,3]]);
L := [R!v : v in RowSequence(M)];
C := Z2Z4AdditiveCode(L, 7);

expectedOutputLength := 19;
expectedOutputBinaryLength := 31;
expectedOutputNgens := 11;
expectedOutputType := [7,12,8,2,7];
expectedOutputCardinal := 4096;
expectedOutputInformationRate := 12.0/31;

assert  Length(C) eq expectedOutputLength;
assert  BinaryLength(C) eq expectedOutputBinaryLength;
assert  Ngens(C) eq expectedOutputNgens;
assert Z2Z4Type(C) eq expectedOutputType;
assert #(C) eq expectedOutputCardinal;
assert  InformationRate(C) eq expectedOutputInformationRate;

/************************************************************/
print "test 10: A Z2Z4-additive code with kappa<>gamma
                alpha = 9, beta = 12, gamma = 6, delta = 3, kappa = 4, length = 21, #C = 4096";
R := RSpace(Z4, 21);
M := Matrix(Z4,[[2,0,0,0,2,2,2,0,2,0,0,0,0,2,0,0,0,0,0,0,2],
                [2,0,0,2,0,0,2,0,2,0,0,2,0,2,0,0,2,0,0,0,0],
                [0,0,2,2,0,2,2,0,0,0,0,0,0,2,2,0,0,0,2,0,0],
                [0,0,2,2,0,2,2,0,0,0,0,0,0,0,0,2,2,0,2,0,0],
                [0,0,0,0,2,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0],
                [2,0,0,2,0,0,2,0,2,0,0,0,0,0,0,0,0,0,0,0,0],
                [0,0,0,0,0,2,2,2,2,0,0,0,1,0,1,1,0,0,1,0,1],
                [0,0,2,0,2,2,0,0,2,1,0,1,0,2,0,1,0,0,0,0,0],
                [0,0,0,2,0,2,2,0,2,0,0,0,1,0,1,1,0,0,1,0,0]]);
L := [R!v : v in RowSequence(M)];
C := Z2Z4AdditiveCode(L, 9);

expectedOutputLength := 21;
expectedOutputBinaryLength := 33;
expectedOutputNgens := 12;
expectedOutputType := [9,12,6,3,4];
expectedOutputCardinal := 4096;
expectedOutputInformationRate := 12.0/33;

assert  Length(C) eq expectedOutputLength;
assert  BinaryLength(C) eq expectedOutputBinaryLength;
assert  Ngens(C) eq expectedOutputNgens;
assert Z2Z4Type(C) eq expectedOutputType;
assert #(C) eq expectedOutputCardinal;
assert  InformationRate(C) eq expectedOutputInformationRate;					


