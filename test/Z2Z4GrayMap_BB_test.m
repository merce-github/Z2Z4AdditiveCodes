/*************************************************************/
/*                                                           */
/* Package or project name: Z2Z4AdditiveCodes package        */
/* Test file name: Z2Z4GrayMap_BB_test.m                     */
/*                                                           */
/* Comments: Black-box tests for the intrinsic functions     */
/*           GrayMap, GrayMapImage,                          */
/*           HasLinearGrayMapImage                           */
/*           included in the Z2Z4AdditiveCodes.m file        */
/*                                                           */
/* Authors: M. Villanueva                                    */
/*                                                           */
/* Revision version and last date: v1.0, 2016/02/19          */
/*                                 v1.1, 2018/01/15          */
/*                                 v1.2, 2018/10/08          */
/*         user defined type       v1.3  2019/01/29          */
/*                                                           */
/*************************************************************/

//needs Z2Z4AdditiveCode file

SetAssertions(true);
Alarm(30*60);

/*************************************************************
    GLOBAL VARIABLES
*************************************************************/    
Z4 := Integers(4);

/************************************************************/
/*                                                          */
/* Function name: GrayMap                                   */
/* Parameters: C                                            */
/* Function description: Given a Z2Z4-additive code C,      */
/*   return the Gray map for C. This is the map Phi from C  */ 
/*   to Phi(C), as defined above.                           */
/* Input parameters description:                            */
/*   - C: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - The Gray map for C                                   */
/*                                                          */
/* Signature: (<Z2Z4Code> C) -> Map                         */
/*                                                          */
/************************************************************/
/************************************************************/
/*                                                          */
/* Function name: GrayMapImage                              */
/* Parameters: C                                            */
/* Function description: Given a Z2Z4-additive code C,      */
/*   return the image of C under the Gray map as a sequence */
/*   of vectors in Z2^(alpha+2*beta). As the resulting      */
/*   image may not be a binary linear code, a sequence of   */
/*   vectors is returned rather than a code.                */
/* Input parameters description:                            */
/*   - C: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - The image of C under the Gray map as a sequence      */
/*     of vectors                                           */
/*                                                          */
/* Signature: (<Z2Z4Code> C) -> SeqEnum                     */
/*                                                          */
/************************************************************/
/************************************************************/
/*                                                          */
/* Function name: HasLinearGrayMapImage                     */
/* Parameters: C                                            */
/* Function description: Given a Z2Z4-additive code C,      */
/*   return true if and only if the image of C under the    */
/*   Gray map is a binary linear code. If so, the function  */
/*   also returns the image B as a binary linear code,      */
/*   together with the bijection Phi: C -> B.               */
/* Input parameters description:                            */
/*   - C: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - true if and only if the image of C under the         */
/*     Gray map is a binary linear code                     */
/*   - The binary linear code if return true                */
/*   - The bijection Phi: C -> B if return true             */ 
/*                                                          */
/* Signature: (<Z2Z4Code> C) -> BoolElt, CodeLinFld, Map    */
/*                                                          */
/************************************************************/
print "THE GRAY MAP";

print "test 1: Trivial Z2Z4-additive zero code
               alpha = 2, beta = 4, gamma = 0, delta = 0, kappa = 0, length = 6, #C = 1";
C := Z2Z4AdditiveZeroCode(2, 4);

grayMap := GrayMap(C);
expectedOutputGrayMapImage := [grayMap(c) : c in C`Code];
expectedOutputHasLinearGrayMap := true;
expectedOutputLinearCode := ZeroCode(GF(2), 10);
OutputHasLinearGrayMap, OutputLinearCode := HasLinearGrayMapImage(C);

assert GrayMapImage(C) eq expectedOutputGrayMapImage;
assert OutputHasLinearGrayMap eq expectedOutputHasLinearGrayMap;
assert OutputLinearCode eq expectedOutputLinearCode;

/************************************************************/
print "test 2: Trivial Z2Z4-additive universe code
               alpha = 2, beta = 4, gamma = 2, delta = 4, kappa = 2, length = 6, #C = 1024";
C := Z2Z4AdditiveUniverseCode(2, 4);

grayMap := GrayMap(C);
expectedOutputGrayMapImage := [grayMap(c) : c in C`Code];
expectedOutputHasLinearGrayMap := true;
expectedOutputLinearCode := UniverseCode(GF(2), 10);
OutputHasLinearGrayMap, OutputLinearCode := HasLinearGrayMapImage(C);

assert GrayMapImage(C) eq expectedOutputGrayMapImage;
assert OutputHasLinearGrayMap eq expectedOutputHasLinearGrayMap;
assert OutputLinearCode eq expectedOutputLinearCode;

/************************************************************/
print "test 3a: Repetition Z2Z4-additive code
               alpha = 4, beta = 8, gamma = 1, delta = 0, kappa = 1, length = 12, #C = 2";
C := Z2Z4AdditiveRepetitionCode(4, 8);

grayMap := GrayMap(C);
expectedOutputGrayMapImage := [grayMap(c) : c in C`Code];
expectedOutputHasLinearGrayMap := true;
expectedOutputLinearCode := RepetitionCode(GF(2), 20);
OutputHasLinearGrayMap, OutputLinearCode := HasLinearGrayMapImage(C);

assert GrayMapImage(C) eq expectedOutputGrayMapImage;
assert OutputHasLinearGrayMap eq expectedOutputHasLinearGrayMap;
assert OutputLinearCode eq expectedOutputLinearCode;

/************************************************************/
print "test 3b: Even weight Z2Z4-additive code
               alpha = 2, beta = 5, gamma = 1, delta = 5, kappa = 1, length = 7, #C = 2048";
C := Z2Z4AdditiveEvenWeightCode(2, 5);

grayMap := GrayMap(C);
expectedOutputGrayMapImage := [grayMap(c) : c in C`Code];
expectedOutputHasLinearGrayMap := true;
expectedOutputLinearCode := EvenWeightCode(12);
OutputHasLinearGrayMap, OutputLinearCode := HasLinearGrayMapImage(C);

assert GrayMapImage(C) eq expectedOutputGrayMapImage;
assert OutputHasLinearGrayMap eq expectedOutputHasLinearGrayMap;
assert OutputLinearCode eq expectedOutputLinearCode;

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

grayMap := GrayMap(C);
expectedOutputGrayMapImage := [grayMap(c) : c in C`Code];
expectedOutputHasLinearGrayMap := false;

assert GrayMapImage(C) eq expectedOutputGrayMapImage;
assert HasLinearGrayMapImage(C) eq expectedOutputHasLinearGrayMap;

/************************************************************/
print "test 5: A Z2Z4-additive code with beta=0
               alpha = 10, beta = 0, gamma = 4, delta = 0, kappa = 4, length = 10, #C = 16";
R := RSpace(Z4, 10);
L := [R![0,0,1,0,1,0,1,0,1,0],
      R![1,0,0,0,0,0,1,1,1,1],
      R![0,0,1,0,0,0,0,1,0,0],
      R![0,1,0,1,0,1,0,1,1,1]];                
C := Z2Z4AdditiveCode(L, 10);

grayMap := GrayMap(C);
expectedOutputGrayMapImage := [grayMap(c) : c in C`Code];
expectedOutputHasLinearGrayMap := true;

assert GrayMapImage(C) eq expectedOutputGrayMapImage;
assert HasLinearGrayMapImage(C) eq expectedOutputHasLinearGrayMap;

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

grayMap := GrayMap(C);
expectedOutputGrayMapImage := [grayMap(c) : c in C`Code];
expectedOutputHasLinearGrayMap := false;

assert GrayMapImage(C) eq expectedOutputGrayMapImage;
assert HasLinearGrayMapImage(C) eq expectedOutputHasLinearGrayMap;

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
C := Z2Z4AdditiveCode(L, 10);

grayMap := GrayMap(C);
expectedOutputGrayMapImage := [grayMap(c) : c in C`Code];
expectedOutputHasLinearGrayMap := true;

assert GrayMapImage(C) eq expectedOutputGrayMapImage;
assert HasLinearGrayMapImage(C) eq expectedOutputHasLinearGrayMap;            

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

grayMap := GrayMap(C);
expectedOutputGrayMapImage := [grayMap(c) : c in C`Code];
expectedOutputHasLinearGrayMap := false;

assert GrayMapImage(C) eq expectedOutputGrayMapImage;
assert HasLinearGrayMapImage(C) eq expectedOutputHasLinearGrayMap;

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

grayMap := GrayMap(C);
expectedOutputGrayMapImage := [grayMap(c) : c in C`Code];
expectedOutputHasLinearGrayMap := false;

assert GrayMapImage(C) eq expectedOutputGrayMapImage;
assert HasLinearGrayMapImage(C) eq expectedOutputHasLinearGrayMap;

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

grayMap := GrayMap(C);
expectedOutputGrayMapImage := [grayMap(c) : c in C`Code];
expectedOutputHasLinearGrayMap := false;

assert GrayMapImage(C) eq expectedOutputGrayMapImage;
assert HasLinearGrayMapImage(C) eq expectedOutputHasLinearGrayMap;

