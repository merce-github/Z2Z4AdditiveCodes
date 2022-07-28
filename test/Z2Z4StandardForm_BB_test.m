/*************************************************************/
/*                                                           */
/* Package or project name: Z2Z4AdditiveCodes package        */
/* Test file name: Z2Z4StandardForm_BB_test.m                */
/*                                                           */
/* Comments: Black-box tests for the intrinsic function      */
/*           StandardForm and IsStandardFormMatrix included  */
/*           in the Z2Z4StandardForm.m file                  */
/*                                                           */
/* Authors: J. Pujol                                         */
/*                                                           */
/* Revision version and last date: v1.0, 2012/06/21          */
/*                                 v1.1, 2015/02/11          */
/*                                 v1.2, 2016/02/19          */
/*                                 v1.3, 2018/01/15          */
/*                                 v1.4, 2018/03/06          */
/*         user defined type       v1.5  2019/01/29          */
/*                                                           */
/*************************************************************/

//needs Z2Z4AdditiveCode file
//needs Z2Z4StandardForm file

SetAssertions(true);
Alarm(30*60);

/*************************************************************
    GLOBAL VARIABLES
*************************************************************/    
Z4:=IntegerRing(4);

/************************************************************/
/*                                                          */
/* Function name: StandardForm                              */
/* Parameters: C                                            */
/* Function description: Given any Z2Z4-additive code C,    */
/*   return a permutation-equivalent Z2Z4-additive code Csf */
/*   in standard form, together with the corresponding      */
/*   isomorphism from C to Csf, the generator matrix in     */
/*   standard form, and the coordinate permutation used to  */
/*   define the isomorphism.                                */
/* Input parameters description:                            */
/*   - C: a Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - a permutation-equivalent Z2Z4-additive code Csf      */
/*   - the corresponding isomorphism f from C to Csf        */
/*   - the generator matrix Gsf in standard form            */
/*   - the coordinate permutation p used to define the      */
/*     isomorphism                                          */
/*                                                          */
/* Signature:(<Z2Z4Code> C) -> Z2Z4Code, Map,               */
/*                             ModMatRngElt, GrpPermElt     */
/*                                                          */
/************************************************************/
/************************************************************/
/*                                                          */
/* Function name: IsStandardFormMatrix                      */
/* Parameters: M, L                                         */
/* Function description: Return true if and only if the     */
/*   matrix M over Z4, where the ones in the first alpha    */
/*   coordinates are represented by twos, is a generator    */
/*   matrix of a Z2Z4-additive code of type (alpha, beta;   */
/*   gamma, delta; kappa) in standard form.                 */
/* Input parameters description:                            */
/*   - M: a matrix over Z4                                  */
/*   - L: a sequence [alpha, beta, gamma, delta, kappa]     */
/* Output parameters description:                           */
/*   - A boolean, true if and only if M is in standard form */
/*                                                          */
/* Signature:(<ModMatRngElt> M, <SeqEnum> L) -> BoolElt     */
/* Signature:(<AlgMatRngElt> M, <SeqEnum> L) -> BoolElt     */
/*                                                          */
/************************************************************/
print "THE STANDARD FORM";

print "test 1: Trivial Z2Z4-additive zero code 
             alpha = 4, beta = 8, gamma = 0, delta = 0, kappa = 0, length = 12, #C = 1"; 
//C := Z2Z4AdditiveZeroCode(4, 8);
//expectedOutput := "C can not be the zero code";
//print expectedOutput;
//Output := Z2Z4StandardForm(C);

/************************************************************/
print "test 2: Trivial Z2Z4-additive universe code 
               alpha = 4, beta = 8, gamma = 4, delta = 8, kappa = 4, length = 12, #C = 1048576";
C := Z2Z4AdditiveUniverseCode(4, 8);
expectedOutputGsf := Matrix(Z4,[[2,0,0,0,0,0,0,0,0,0,0,0],
                                [0,2,0,0,0,0,0,0,0,0,0,0],
                                [0,0,2,0,0,0,0,0,0,0,0,0],
                                [0,0,0,2,0,0,0,0,0,0,0,0],
                                [0,0,0,0,1,0,0,0,0,0,0,0],
                                [0,0,0,0,0,1,0,0,0,0,0,0],
                                [0,0,0,0,0,0,1,0,0,0,0,0],
                                [0,0,0,0,0,0,0,1,0,0,0,0],
                                [0,0,0,0,0,0,0,0,1,0,0,0],
                                [0,0,0,0,0,0,0,0,0,1,0,0],
                                [0,0,0,0,0,0,0,0,0,0,1,0],
                                [0,0,0,0,0,0,0,0,0,0,0,1]]);
expectedOutputC := C;

OutputCsf, Outputf, OutputGsf, Outputp := StandardForm(C); 
OutputCsf2 := Z2Z4AdditiveCode(GeneratorMatrix(C`Code)*PermutationMatrix(Z4, Outputp), 4);
OutputC := Z2Z4AdditiveCode(OutputGsf*PermutationMatrix(Z4, Outputp^(-1)), 4);

assert expectedOutputGsf eq OutputGsf;
assert expectedOutputC eq OutputC;
assert OutputCsf eq OutputCsf2;
assert IsStandardFormMatrix(OutputGsf, [4,8,4,8,4]);
//assert {Outputf(c) : c in C`Code} eq Set(OutputCsf`Code);   //Too slow

/************************************************************/
print "test 3: Repetition Z2Z4-additive code
               alpha = 4, beta = 8, gamma = 1, delta = 0, kappa = 1, length = 12, #C = 2";
C := Z2Z4AdditiveRepetitionCode(4, 8);
expectedOutputGsf := Matrix(Z4, [[2,2,2,2,2,2,2,2,2,2,2,2]]);
expectedOutputC := C;

OutputCsf, Outputf, OutputGsf, Outputp := StandardForm(C); 
OutputCsf2 := Z2Z4AdditiveCode(GeneratorMatrix(C`Code)*PermutationMatrix(Z4, Outputp), 4);
OutputC := Z2Z4AdditiveCode(OutputGsf*PermutationMatrix(Z4, Outputp^(-1)), 4);

assert expectedOutputGsf eq OutputGsf;
assert expectedOutputC eq OutputC;
assert OutputCsf eq OutputCsf2;
assert IsStandardFormMatrix(OutputGsf, [4,8,1,0,1]);
assert {Outputf(c) : c in C`Code} eq Set(OutputCsf`Code);

/************************************************************/
print "test 4: A Z2Z4-additive code with alpha, beta even, p=Id.
               alpha = 10, beta = 20, gamma = 6, delta = 2, kappa = 1, length = 30, #C = 1024";
M := Matrix(Z4,[[2,0,2,0,0,0,2,0,2,2,0,0,0,1,0,0,0,3,2,1,2,0,1,2,0,1,2,1,0,3],
                [0,2,2,0,0,2,0,0,2,2,1,1,1,1,0,1,0,3,2,0,0,2,2,0,1,0,3,2,1,1],
                [0,0,0,0,2,0,2,2,2,0,0,0,0,1,0,0,0,1,2,3,0,0,1,0,2,3,0,3,2,3],
                [0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,2,2,0,0,2,0,2,2,0,0,2,2],
                [0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2,0,2,2,2,2,2,0,0,2,2,2,0],
                [0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,2,0,2,0,0,0,2,2,0,2,2,2,2],
                [0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,0,2,0,0,2,0,0,2,0,2,0,2],
                [0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,2,0,2,0,0,2,2,0,2,0,0,2,2],
                [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,2,0,2,2,2,0,2,0,2,2,0,0],
                [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,0,0,0,0,0,0,2,2,0,2,2,2]]);
C := Z2Z4AdditiveCode(M, 10);
expectedOutputC := C;

OutputCsf, Outputf, OutputGsf, Outputp := StandardForm(C); 
OutputCsf2 := Z2Z4AdditiveCode(GeneratorMatrix(C`Code)*PermutationMatrix(Z4, Outputp), 10);
OutputC := Z2Z4AdditiveCode(OutputGsf*PermutationMatrix(Z4, Outputp^(-1)), 10);

assert expectedOutputC eq OutputC;
assert OutputCsf eq OutputCsf2;
assert IsStandardFormMatrix(OutputGsf, [10,20,6,2,1]);
assert {Outputf(c) : c in C`Code} eq Set(OutputCsf`Code);

/************************************************************/
print "test 5: A Z2Z4-additive code with alpha even, beta odd, p <> Id in delta quaternary part
               alpha = 10, beta = 19, gamma = 7, delta = 2, kappa = 3, length = 29, #C = 2048";
M := Matrix(Z4,[[2,0,2,0,0,0,2,0,2,2,0,0,0,1,0,0,0,3,2,1,2,0,1,2,0,1,2,1,0],
                [0,2,2,0,0,2,0,0,2,2,1,1,1,1,0,1,0,3,2,0,0,2,2,0,2,0,3,2,2],
                [0,0,0,0,2,0,2,2,2,0,0,0,0,1,0,0,0,1,2,3,0,0,1,0,2,3,0,3,2],
                [0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,2,2,0,0,2,0,2,2,0,0,2],
                [0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2,0,2,2,2,2,2,0,0,2,2,2],
                [0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,2,0,2,0,0,0,2,2,0,2,2,2],
                [0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,0,2,0,0,2,0,0,2,0,2,0],
                [0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,2,0,2,0,0,2,2,0,2,0,0,2],
                [0,2,0,2,0,0,0,0,0,0,0,0,0,0,0,2,0,0,2,0,2,2,2,0,2,0,2,2,0],
                [0,0,2,0,2,0,2,0,2,0,2,0,0,0,0,0,2,2,0,0,0,0,0,0,2,2,0,2,2]]);
C := Z2Z4AdditiveCode(M, 10);
expectedOutputC := C;

OutputCsf, Outputf, OutputGsf, Outputp := StandardForm(C); 
OutputCsf2 := Z2Z4AdditiveCode(GeneratorMatrix(C`Code)*PermutationMatrix(Z4, Outputp), 10);
OutputC := Z2Z4AdditiveCode(OutputGsf*PermutationMatrix(Z4, Outputp^(-1)), 10);

assert expectedOutputC eq OutputC;
assert OutputCsf eq OutputCsf2;
assert IsStandardFormMatrix(OutputGsf, [10,19,7,2,3]);
assert {Outputf(c) : c in C`Code} eq Set(OutputCsf`Code);    
    
/************************************************************/                      
print "test 6: A Z2Z4-additive code with alpha even, beta odd, p <> Id in non delta quaternary part.    
               alpha = 10, beta = 19, gamma = 7, delta = 2, kappa = 3, length = 29, #C = 2048";
M := Matrix(Z4,[[2,0,0,0,0,0,2,2,2,2,0,0,2,2,0,0,2,0,2,0,2,0,0,0,0,0,0,0,0],
                [0,2,0,2,0,0,0,0,0,0,0,0,2,0,2,2,0,0,2,0,2,2,2,0,0,2,0,0,0],
                [0,0,2,0,2,0,2,0,2,0,0,0,0,2,0,0,2,0,2,0,0,0,2,0,0,0,0,0,0],
                [0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,0,2,2,2,0,0,0,2,0,0,0,0,0],
                [0,0,0,0,0,0,0,0,0,0,2,2,0,0,2,2,0,0,0,2,0,0,0,0,2,0,0,0,0],
                [0,0,0,0,0,0,0,0,0,0,2,0,2,0,2,2,0,0,0,2,2,2,0,0,0,0,0,0,0],
                [0,0,0,0,0,0,0,0,0,0,2,2,2,0,0,0,0,0,2,2,2,2,2,0,0,0,2,0,0],
                [0,0,0,2,2,2,2,0,0,2,3,3,3,3,0,3,2,3,0,2,0,2,0,0,0,0,0,1,0],
                [0,0,0,0,2,0,2,2,2,0,2,0,0,3,0,0,0,3,0,3,0,0,3,0,0,1,0,0,1]]);
C := Z2Z4AdditiveCode(M, 10);                  
expectedOutputC := C;

OutputCsf, Outputf, OutputGsf, Outputp := StandardForm(C);
OutputCsf2 := Z2Z4AdditiveCode(GeneratorMatrix(C`Code)*PermutationMatrix(Z4, Outputp), 10);
OutputC := Z2Z4AdditiveCode(OutputGsf*PermutationMatrix(Z4, Outputp^(-1)), 10);

assert expectedOutputC eq OutputC;    
assert OutputCsf eq OutputCsf2;
assert IsStandardFormMatrix(OutputGsf, [10,19,7,2,3]);    
assert {Outputf(c) : c in C`Code} eq Set(OutputCsf`Code);                          

/************************************************************/ 
print "test 7: A Z2Z4-additive code with alpha even, beta odd, p <> Id in alpha part.
               alpha = 10, beta = 19, gamma = 7, delta = 2, kappa = 3, length = 29, #C = 2048";
M := Matrix(Z4,[[0,0,2,0,2,0,2,0,2,0,0,0,2,0,0,0,2,0,2,0,0,0,0,0,2,2,2,2,0],
                [2,0,0,0,0,0,2,2,2,2,0,0,0,0,0,0,2,0,2,0,2,0,2,0,2,2,2,2,0],
                [2,0,2,0,2,0,0,2,0,2,0,0,0,0,2,2,0,0,2,0,0,0,0,2,0,2,2,0,0],
                [0,0,2,0,2,0,2,0,2,0,0,2,0,0,0,0,2,0,2,0,2,0,0,2,2,2,0,2,0],
                [0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,2,0,0],
                [2,0,0,0,0,0,2,2,2,2,0,0,0,0,2,0,2,2,2,2,2,0,0,0,0,0,2,2,0],
                [2,0,2,2,2,0,0,2,0,2,0,0,0,0,2,0,0,0,0,0,2,0,0,0,0,2,0,2,2],
                [0,0,2,0,0,0,0,2,0,0,0,0,0,1,2,0,2,1,0,3,0,1,1,0,2,0,2,0,1],
                [0,2,0,2,0,2,0,2,2,2,1,1,1,0,0,1,2,0,2,3,0,1,1,2,0,2,2,3,1]]);
C := Z2Z4AdditiveCode(M, 10);
expectedOutputC := C;

OutputCsf, Outputf, OutputGsf, Outputp := StandardForm(C); 
OutputCsf2 := Z2Z4AdditiveCode(GeneratorMatrix(C`Code)*PermutationMatrix(Z4, Outputp), 10);
OutputC := Z2Z4AdditiveCode(OutputGsf*PermutationMatrix(Z4, Outputp^(-1)), 10);

assert expectedOutputC eq OutputC;    
assert OutputCsf eq OutputCsf2;
assert IsStandardFormMatrix(OutputGsf, [10,19,7,2,3]);    
assert {Outputf(c) : c in C`Code} eq Set(OutputCsf`Code);                       

/************************************************************/
print "test 8: A Z2Z4-additive code with alpha=0
               alpha = 0, beta = 19, gamma = 7, delta = 2, kappa = 0, length = 19, #C = 2048";
M := Matrix(Z4,[[0,0,2,0,0,0,2,0,2,0,0,0,0,0,2,2,2,2,0],
                [0,0,0,0,0,0,2,0,2,0,2,0,2,0,2,2,2,2,0],
                [0,0,0,0,2,2,0,0,2,0,0,0,0,2,0,2,2,0,0],
                [0,2,0,0,0,0,2,0,2,0,2,0,0,2,2,2,0,2,0],
                [0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,2,0,0],
                [0,0,0,0,2,0,2,2,2,2,2,0,0,0,0,0,2,2,0],
                [0,0,0,0,2,0,0,0,0,0,2,0,0,0,0,2,0,2,2],
                [0,0,0,1,2,0,2,1,0,3,0,1,1,0,2,0,2,0,1],
                [1,1,1,0,0,1,2,0,2,3,0,1,1,2,0,2,2,3,1]]);
C := Z2Z4AdditiveCode(M, 0);
expectedOutputC := C;

OutputCsf, Outputf, OutputGsf, Outputp := StandardForm(C);  
OutputCsf2 := Z2Z4AdditiveCode(GeneratorMatrix(C`Code)*PermutationMatrix(Z4, Outputp), 0);
OutputC := Z2Z4AdditiveCode(OutputGsf*PermutationMatrix(Z4, Outputp^(-1)), 0);

assert expectedOutputC eq OutputC;    
assert OutputCsf eq OutputCsf2;
assert IsStandardFormMatrix(OutputGsf, [0,19,7,2,0]);    
assert {Outputf(c) : c in C`Code} eq Set(OutputCsf`Code);                       

/************************************************************/
print "test 9: A Z2Z4-additive code with beta=0
                alpha = 10, beta = 0, gamma = 5, delta = 0, kappa = 5, length = 10, #C = 32";
M := Matrix(Z4,[[0,0,2,0,2,0,2,0,2,0],
                [2,0,0,0,0,0,2,2,2,2],
                [2,0,2,0,2,0,0,2,0,2],
                [0,0,2,0,2,0,2,0,2,0],
                [0,0,0,2,0,0,0,0,0,0],
                [2,0,0,0,0,0,2,2,2,2],
                [2,0,2,2,2,0,0,2,0,2],
                [0,0,2,0,0,0,0,2,0,0],
                [0,2,0,2,0,2,0,2,2,2]]);
C := Z2Z4AdditiveCode(M, 10);                  
expectedOutputC := C;

OutputCsf, Outputf, OutputGsf, Outputp := StandardForm(C);
OutputCsf2 := Z2Z4AdditiveCode(GeneratorMatrix(C`Code)*PermutationMatrix(Z4, Outputp), 10);
OutputC := Z2Z4AdditiveCode(OutputGsf*PermutationMatrix(Z4, Outputp^(-1)), 10);

assert expectedOutputC eq OutputC;    
assert OutputCsf eq OutputCsf2;
assert IsStandardFormMatrix(OutputGsf, [10,0,5,0,5]);
assert {Outputf(c) : c in C`Code} eq Set(OutputCsf`Code);

/************************************************************/
print "test 10: A Z2Z4-additive code with gamma=0
                alpha = 10, beta = 19, gamma = 0, delta = 7, kappa = 0, length = 19, #C = 16384";
M := Matrix(Z4,[[2,2,0,0,2,2,0,2,2,2,0,0,0,1,0,0,0,0,2,3,2,1,1,3,3,1,2,3,0],
                [0,2,0,0,0,0,2,0,0,0,0,0,0,0,0,2,0,1,0,3,0,2,2,1,3,1,1,0,1],
                [0,0,2,2,0,0,2,2,2,2,0,1,0,0,0,0,0,0,2,3,1,1,0,2,1,0,0,0,1],
                [2,0,0,0,2,0,2,0,0,0,0,0,0,0,0,0,1,0,3,3,1,3,1,0,0,3,1,0,2],
                [2,2,2,2,2,0,0,0,0,0,1,0,0,0,0,1,0,0,0,2,3,0,0,1,1,0,2,2,1],
                [2,0,2,2,2,0,2,2,0,2,0,0,1,0,0,3,0,0,0,1,2,3,1,3,0,2,1,3,2],
                [0,2,2,0,2,0,2,2,0,0,0,0,0,0,1,1,0,0,0,3,3,2,0,2,0,3,1,1,3]]);
C := Z2Z4AdditiveCode(M, 10);
expectedOutputC := C;

OutputCsf, Outputf, OutputGsf, Outputp := StandardForm(C); 
OutputCsf2 := Z2Z4AdditiveCode(GeneratorMatrix(C`Code)*PermutationMatrix(Z4, Outputp), 10);
OutputC := Z2Z4AdditiveCode(OutputGsf*PermutationMatrix(Z4, Outputp^(-1)), 10);

assert expectedOutputC eq OutputC;    
assert OutputCsf eq OutputCsf2;
assert IsStandardFormMatrix(OutputGsf, [10,19,0,7,0]);
assert {Outputf(c) : c in C`Code} eq Set(OutputCsf`Code);

/************************************************************/
print "test 11: A Z2Z4-additive code with delta=0
                alpha = 10, beta = 19, gamma = 7, delta = 0, kappa = 4, length = 19, #C = 128";
M := Matrix(Z4,[[2,0,0,0,0,0,2,2,2,0,0,0,0,0,0,2,2,0,0,0,0,0,0,2,0,2,0,2,2],
                [0,2,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,2,2,2,2],
                [0,0,2,0,0,2,2,0,2,0,0,0,0,0,2,0,2,2,0,2,0,0,0,2,2,2,2,0,0],
                [0,0,0,2,0,2,0,0,2,0,0,0,0,0,2,0,0,0,2,0,2,2,0,0,2,0,2,0,2],
                [0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,0,2,0,2,0,0,2,0,2,0,0,0,0,2],
                [0,0,0,0,0,0,0,0,0,0,0,2,0,2,0,2,2,2,0,2,0,2,0,2,2,2,2,0,2],
                [0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,2,2,0,0,0,0,0,2,2,0,2,2]]);
C := Z2Z4AdditiveCode(M, 10);
expectedOutputC := C;

OutputCsf, Outputf, OutputGsf, Outputp := StandardForm(C);
OutputCsf2 := Z2Z4AdditiveCode(GeneratorMatrix(C`Code)*PermutationMatrix(Z4, Outputp), 10);
OutputC := Z2Z4AdditiveCode(OutputGsf*PermutationMatrix(Z4, Outputp^(-1)), 10);

assert expectedOutputC eq OutputC;    
assert OutputCsf eq OutputCsf2;
assert IsStandardFormMatrix(OutputGsf, [10,19,7,0,4]);    
assert {Outputf(c) : c in C`Code} eq Set(OutputCsf`Code);                           

/************************************************************/
print "test 12: A Z2Z4-additive code with kappa=0
                alpha = 10, beta = 19, gamma = 7, delta = 2, kappa = 0, length = 29, #C = 2048";
M := Matrix(Z4,[[0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,0,2,0,0,0,0,0,2,2,2,2,0],
                [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,2,0,2,0,2,0,2,2,2,2,0],
                [0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,0,0,2,0,0,0,0,2,0,2,2,0,0],
                [0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,2,0,2,0,2,0,0,2,2,2,0,2,0],
                [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,2,0,0],
                [0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,2,2,2,2,2,0,0,0,0,0,2,2,0],
                [0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2,0,0,0,0,2,0,2,2],
                [0,0,2,0,0,0,2,0,0,0,0,0,0,1,2,0,2,1,0,3,0,1,1,0,2,0,2,0,1],
                [0,2,0,2,0,2,0,2,2,2,1,1,1,0,0,1,2,0,2,3,0,1,1,2,0,2,2,3,1]]);
C := Z2Z4AdditiveCode(M, 10);
expectedOutputC := C;

OutputCsf, Outputf, OutputGsf, Outputp := StandardForm(C); 
OutputCsf2 := Z2Z4AdditiveCode(GeneratorMatrix(C`Code)*PermutationMatrix(Z4, Outputp), 10);
OutputC := Z2Z4AdditiveCode(OutputGsf*PermutationMatrix(Z4, Outputp^(-1)), 10);

assert expectedOutputC eq OutputC;
assert OutputCsf eq OutputCsf2;
assert IsStandardFormMatrix(OutputGsf, [10,19,7,2,0]);
assert {Outputf(c) : c in C`Code} eq Set(OutputCsf`Code);

/************************************************************/
print "test 13: A Z2Z4-additive code with kappa=alpha
                alpha = 7, beta = 12, gamma = 8, delta = 2, kappa = 7, length = 19, #C = 4096";
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
C := Z2Z4AdditiveCode(M, 7);
expectedOutputC := C;

OutputCsf, Outputf, OutputGsf, Outputp := StandardForm(C); 
OutputCsf2 := Z2Z4AdditiveCode(GeneratorMatrix(C`Code)*PermutationMatrix(Z4, Outputp), 7);
OutputC := Z2Z4AdditiveCode(OutputGsf*PermutationMatrix(Z4, Outputp^(-1)), 7);

assert expectedOutputC eq OutputC;
assert OutputCsf eq OutputCsf2;
assert IsStandardFormMatrix(OutputGsf, [7,12,8,2,7]);
assert {Outputf(c) : c in C`Code} eq Set(OutputCsf`Code);
                           
/************************************************************/
print "test 14: A Z2Z4-additive code with kappa=gamma
                alpha = 9, beta = 12, gamma = 8, delta = 2, kappa = 8, length = 21, #C = 4096";
M := Matrix(Z4,[[0,0,0,2,0,2,2,2,2,0,0,2,0,0,0,0,2,0,2,0,0],
                [0,0,0,0,2,2,2,0,2,0,0,0,0,0,0,2,2,0,0,0,0],
                [2,2,2,0,0,0,2,2,2,0,0,0,0,0,2,0,0,0,0,0,0],
                [0,0,0,0,2,0,2,2,0,0,0,0,2,0,0,0,2,0,0,0,0],
                [0,0,0,2,0,0,0,2,0,0,0,0,0,0,0,0,0,0,2,2,0],
                [2,0,0,0,0,0,0,2,2,0,0,0,0,2,0,0,0,0,0,0,0],
                [2,0,2,2,2,0,2,0,2,0,0,0,0,0,0,0,0,2,0,0,0],
                [2,0,0,2,0,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,2],
                [2,2,0,2,0,2,0,0,2,0,1,0,0,1,1,0,2,1,2,1,1],
                [2,2,0,0,0,2,0,2,0,1,0,1,1,1,0,1,0,1,2,1,0]]);
C := Z2Z4AdditiveCode(M, 9);
expectedOutputC := C;

OutputCsf, Outputf, OutputGsf, Outputp := StandardForm(C); 
OutputCsf2 := Z2Z4AdditiveCode(GeneratorMatrix(C`Code)*PermutationMatrix(Z4, Outputp), 9);
OutputC := Z2Z4AdditiveCode(OutputGsf*PermutationMatrix(Z4, Outputp^(-1)), 9);

assert expectedOutputC eq OutputC;
assert OutputCsf eq OutputCsf2;
assert IsStandardFormMatrix(OutputGsf, [9,12,8,2,8]);
assert {Outputf(c) : c in C`Code} eq Set(OutputCsf`Code);

/************************************************************/
print "test 15: A Z2Z4-additive code with a zero column in alpha part.
                alpha = 9, beta = 12, gamma = 6, delta = 3, kappa = 4, length = 21, #C = 4096";
M := Matrix(Z4,[[2,0,0,0,2,2,2,0,2,0,0,0,0,2,0,0,0,0,0,0,2],
                [2,0,0,2,0,0,2,0,2,0,0,2,0,2,0,0,2,0,0,2,0],
                [0,0,2,2,0,2,2,0,0,0,0,0,0,2,2,0,0,0,2,2,0],
                [0,0,2,2,0,2,2,0,0,0,0,0,0,0,0,2,2,0,2,2,0],
                [0,0,0,0,2,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0],
                [2,0,0,2,0,0,2,0,2,0,0,0,0,0,0,0,0,2,0,2,0],
                [0,0,0,0,0,2,2,2,2,0,0,0,1,0,1,1,0,0,1,2,1],
                [0,0,2,0,2,2,0,0,2,1,0,1,0,2,0,1,0,1,0,3,0],
                [0,0,0,2,0,2,2,0,2,0,1,1,0,1,0,0,2,0,1,0,0]]);
C := Z2Z4AdditiveCode(M, 9);
expectedOutputC := C;

OutputCsf, Outputf, OutputGsf, Outputp := StandardForm(C); 
OutputCsf2 := Z2Z4AdditiveCode(GeneratorMatrix(C`Code)*PermutationMatrix(Z4, Outputp), 9);
OutputC := Z2Z4AdditiveCode(OutputGsf*PermutationMatrix(Z4, Outputp^(-1)), 9);

assert expectedOutputC eq OutputC;
assert OutputCsf eq OutputCsf2;
assert IsStandardFormMatrix(OutputGsf, [9,12,6,3,4]);
assert {Outputf(c) : c in C`Code} eq Set(OutputCsf`Code);

/************************************************************/
print "test 16: A Z2Z4-additive code with a zero column in beta part.
                alpha = 9, beta = 12, gamma = 6, delta = 3, kappa = 4, length = 21, #C = 4096";
M := Matrix(Z4,[[2,0,0,0,2,2,2,0,2,0,0,0,0,2,0,0,0,0,0,0,2],
                [2,0,0,2,0,0,2,0,2,0,0,2,0,2,0,0,2,0,0,2,0],
                [0,0,2,2,0,2,2,0,0,0,0,0,0,2,2,0,0,0,2,2,0],
                [0,0,2,2,0,2,2,0,0,0,0,0,0,0,0,2,2,0,2,2,0],
                [0,0,0,0,2,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0],
                [2,0,0,2,0,0,2,0,2,0,0,0,0,0,0,0,0,0,0,2,0],
                [0,0,0,0,0,2,2,2,2,0,0,0,1,0,1,1,0,0,1,2,1],
                [0,0,2,0,2,2,0,0,2,1,0,1,0,2,0,1,0,0,0,3,0],
                [0,0,0,2,0,2,2,0,2,0,1,1,0,1,0,0,2,0,1,0,0]]);
C := Z2Z4AdditiveCode(M, 9);
expectedOutputC := C;

OutputCsf, Outputf, OutputGsf, Outputp := StandardForm(C); 
OutputCsf2 := Z2Z4AdditiveCode(GeneratorMatrix(C`Code)*PermutationMatrix(Z4, Outputp), 9);
OutputC := Z2Z4AdditiveCode(OutputGsf*PermutationMatrix(Z4, Outputp^(-1)), 9);

assert expectedOutputC eq OutputC;
assert OutputCsf eq OutputCsf2;
assert IsStandardFormMatrix(OutputGsf, [9,12,6,3,4]);
assert {Outputf(c) : c in C`Code} eq Set(OutputCsf`Code);
    
/************************************************************/                           
print "test 17: Another Z2Z4-additive code with two zero columns in beta part.
                alpha = 9, beta = 12, gamma = 6, delta = 3, kappa = 4, length = 21, #C = 4096";
M := Matrix(Z4,[[2,0,0,0,2,2,2,0,2,0,0,0,0,2,0,0,0,0,0,0,2],
                [2,0,0,2,0,0,2,0,2,0,0,2,0,2,0,0,2,0,0,0,0],
                [0,0,2,2,0,2,2,0,0,0,0,0,0,2,2,0,0,0,2,0,0],
                [0,0,2,2,0,2,2,0,0,0,0,0,0,0,0,2,2,0,2,0,0],
                [0,0,0,0,2,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0],
                [2,0,0,2,0,0,2,0,2,0,0,0,0,0,0,0,0,0,0,0,0],
                [0,0,0,0,0,2,2,2,2,0,0,0,1,0,1,1,0,0,1,0,1],
                [0,0,2,0,2,2,0,0,2,1,0,1,0,2,0,1,0,0,0,0,0],
                [0,0,0,2,0,2,2,0,2,0,0,0,1,0,1,1,0,0,1,0,0]]);
C := Z2Z4AdditiveCode(M, 9);
expectedOutputC := C;

OutputCsf, Outputf, OutputGsf, Outputp := StandardForm(C); 
OutputCsf2 := Z2Z4AdditiveCode(GeneratorMatrix(C`Code)*PermutationMatrix(Z4, Outputp), 9);
OutputC := Z2Z4AdditiveCode(OutputGsf*PermutationMatrix(Z4, Outputp^(-1)), 9);

assert expectedOutputC eq OutputC;
assert OutputCsf eq OutputCsf2;
assert IsStandardFormMatrix(OutputGsf, [9,12,6,3,4]);
assert {Outputf(c) : c in C`Code} eq Set(OutputCsf`Code);

