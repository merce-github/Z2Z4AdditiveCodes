/*************************************************************/
/*                                                           */
/* Package or project name: Z2Z4AdditiveCodes package        */
/* Test file name: Z2Z4AdditiveCodes_BB_test.m               */
/*                                                           */
/* Comments: Black-box tests for the intrinsic functions     */
/*           Z2Z4AdditiveCode, Z2Z4AdditiveZeroCode,         */
/*           Z2Z4AdditiveUniverseCode, Z2Z4RepetitionCode,   */
/*           Z2Z4HadamardCode, Z2Z4ExtendedPerfectCode and   */
/*           RandomZ2Z4AdditiveCode                          */
/*           included in the Z2Z4AdditiveCodes.m file        */
/*                                                           */
/* Authors: Jaume Pujol and Mercè Villanueva                 */
/*                                                           */
/* Revision version and last date: v1.0, 2015/02/11          */
/*                                 v2.0, 2016/02/19          */
/*                                 v2.1, 2018/01/13          */
/*                                 v2.2, 2018/02/23          */
/*                                 v2.3, 2018/10/08          */
/*                                 v2.4, 2019/01/26          */
/*         user defined type       v2.5  2019/01/29          */
/*                                                           */
/*************************************************************/

//needs Z2Z4AdditiveCode file

SetAssertions(true);
Alarm(30*60);

/****************************************************************
    GLOBAL VARIABLES
*****************************************************************/    
Z4 := Integers(4);

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4AdditiveCode                              */
/* Parameters: L                                                */
/* Function description: Create a Z2Z4-additive code C.         */
/*   Alpha is the length of the binary part of the Z2Z4-additive*/
/*   code, and Code is the quaternary linear code equal to the  */
/*   Z2Z4-additive code, where the ones in the first Alpha      */
/*   coordinates are represented by twos.                       */  
/* Input parameters description:                                */
/*   -    a sequence of elements of Z4^n, or a subspace         */
/*        of Z4^n, or a m x n matrix M over the ring Z4, or     */
/*        a quaternary linear code, or a binary linear code,    */
/*        or a sequence of elements of Z2^alpha x Z4^beta,      */
/*        or two matrices M over Z2^alpha and N over Z4^beta,   */
/*        or two codes D, a binary linear code, and E, a        */
/*        quaternary linear code.                               */
/* Output parameters description:                               */
/*   - a Z2Z4-additive code                                     */
/*                                                              */
/* Signature: (<ModMatRngElt> M, <RngIntElt> alpha) -> Z2Z4Code */
/*            (<AlgMatElt> M, <RngIntElt> alpha) -> Z2Z4Code    */
/*            (<SeqEnum> L, <RngIntElt> alpha) -> Z2Z4Code      */
/*            (<ModTupRng> V, <RngIntElt> alpha) -> Z2Z4Code    */
/*            (<CodeLinRng> C) -> Z2Z4Code                      */
/*            (<CodeLinRng> C, <RngIntElt> alpha) -> Z2Z4Code   */
/*            (<CodeLinFld> C) -> Z2Z4Code                      */
/*            (<SeqEnum> T) -> Z2Z4Code                         */
/*            (<ModMatRngElt> M, <ModMatRngElt> N) -> Z2Z4Code  */
/*            (<AlgMatElt> M, <AlgMatElt> N) -> Z2Z4Code        */
/*            (<CodeLinFld> D, <CodeLinRng> E) -> Z2Z4Code      */
/*                                                              */
/****************************************************************/
/****************************************************************/
/*                                                              */
/* Function name: Z2Z4HadamardCode                              */
/* Parameters:  delta, m                                        */
/* Function description: The function returns a Z2Z4-additive   */
/*   Hadamard code of type (alpha, beta; gamma, delta; kappa).  */
/*   The parameter OverZ4 specifies whether the code is over Z4,*/
/*   that is alpha=0, or, otherwise, alpha<>0. The default value*/
/*   is false. When OverZ4 is true, given an integer m >= 1 and */
/*   an integer delta such that 1 <= delta <= Floor((m+1)/2),   */
/*   return a Z2Z4-additive Hadamard code of type (0,beta;gamma,*/
/*   delta;0), where beta=2^{m-1} and gamma=m+1-2*delta. When   */
/*   OverZ4 is false, given an integer m >= 1 and an integer    */
/*   delta such that 0 <= delta <= Floor(m/2), return a Z2Z4-   */
/*   additive Hadamard code of type (alpha, beta; gamma, delta; */
/*   gamma), where alpha=2^{m-delta}, beta=2^{m-1}-2^{m-delta-1}*/
/*   and gamma=m+1-2*delta. Moreover, return a generator matrix */
/*   with gamma+delta rows constructed in a recursive way, where*/
/*   the ones in the first alpha coordinates are represented by */
/*   twos.                                                      */
/* Input parameters description:                                */
/*   - delta : A positive integer                               */
/*   - m :  A positive integer                                  */
/* Output parameters description:                               */
/*   -  A record representing a Z2Z4-additive code              */
/*                                                              */
/* Signature: (<RngIntElt> delta, <RngIntElt> m ) -> Z2Z4Code   */
/*                                                              */
/****************************************************************/ 
/****************************************************************/
/*                                                              */
/* Function name: Z2Z4ExtendedPerfectCode                       */
/* Parameters:  delta, m                                        */
/* Function description: The function returns a Z2Z4-additive   */
/*   extended perfect code of type alpha, beta; gamma, delta;   */
/*   kappa). The parameter OverZ4 specifies whether the code is */
/*   over Z4, that is alpha=0, or, otherwise, alpha<>0. The     */
/*   default value is false. When OverZ4 is true, given an      */
/*   integer m >= 1 and an integer delta such that 1 <=delta<=  */
/*   Floor((m+1)/2), return a Z2Z4-additive extended perfect    */
/*   code, such that its additive dual code is of type of type  */
/*   (0, beta; gamma, delta; 0), where beta=2^{m-1} and gamma=  */
/*   m+1-2*delta. When OverZ4 is false, given an integer m >= 1 */
/*   and an integer delta such that 0 <= delta <= Floor(m/2),   */
/*   return a Z2Z4-additive extended perfect code, such that its*/
/*   additive dual code is of type (alpha, beta; gamma, delta;  */
/*   gamma), where alpha=2^{m-delta}, beta=2^{m-1}-2^{m-delta-1}*/
/*   and gamma=m+1-2*delta. Moreover, return a parity check     */
/*   matrix with gamma+delta rows constructed in a recursive    */
/*   way, where the ones in the first alpha coordinates are     */
/*   represented by twos.                                       */
/* Input parameters description:                                */
/*   - delta : A positive integer                               */
/*   - m :  A positive integer                                  */
/* Output parameters description:                               */
/*   -  A record representing a Z2Z4-additive code              */
/*                                                              */
/* Signature: (<RngIntElt> delta, <RngIntElt> m ) -> Z2Z4Code   */
/*                                                              */
/****************************************************************/ 
print "CONSTRUCTION OF GENERAL Z2Z4-ADDITIVE CODES AND SOME TRIVIAL ONES";

print "test 1: Trivial Z2Z4-additive zero code 
               alpha = 2, beta = 4, gamma = 0, delta = 0, kappa = 0, length = 6, #C = 1"; 
R := RSpace(Z4, 6);
L := [R!0];
C1 := Z2Z4AdditiveCode(L, 2);
C2 := Z2Z4AdditiveCode(L, 2);
C3 := Z2Z4AdditiveCode(sub<R|L>, 2);   
C4 := Z2Z4AdditiveCode(Matrix(L), 2);
C5 := Z2Z4AdditiveZeroCode(2, 4);

T := FromZ4toZ2Z4(L, 2);
C6 := Z2Z4AdditiveCode(T);

M := Matrix(T[1][1]);
N := Matrix(T[1][2]);
C7 := Z2Z4AdditiveCode(M, N);

D := ZeroCode(GF(2), 2);
E := ZeroCode(Z4, 4);
C8 := Z2Z4AdditiveCode(D, E);

expectedOutputCode := LinearCode<Z4, 6 | L>;
expectedAlpha := 2;

assert C1`Code eq expectedOutputCode;
assert C1`Alpha eq expectedAlpha; 
assert C2`Code eq expectedOutputCode;
assert C2`Alpha eq expectedAlpha; 
assert C3`Code eq expectedOutputCode;
assert C3`Alpha eq expectedAlpha;
assert C4`Code eq expectedOutputCode;
assert C4`Alpha eq expectedAlpha;
assert C5`Code eq expectedOutputCode;
assert C5`Alpha eq expectedAlpha;
assert C1 eq C2;
assert C2 eq C3;
assert C3 eq C4;
assert C4 eq C5;
assert C5 eq C6;
assert C6 eq C7;
assert C7 eq C8;
assert IsSeparable(C1);

/************************************************************/
print "test 2: Trivial Z2Z4-additive universe code 
               alpha = 2, beta = 4, gamma = 2, delta = 4, kappa = 2, length = 6, #C = 1024";
R := RSpace(Z4, 6);
L1 := Basis(R); 
C1 := Z2Z4AdditiveCode(L1, 2);
L2 := L1;  
L2[1][1] := 2;
L2[2][2] := 2;
C2 := Z2Z4AdditiveCode(L2, 2);
C3 := Z2Z4AdditiveCode(sub<R|L2>, 2);
C4 := Z2Z4AdditiveCode(Matrix(L2), 2);
C5 := Z2Z4AdditiveUniverseCode(2, 4);

T := FromZ4toZ2Z4(L2, 2);
C6 := Z2Z4AdditiveCode(T);

M := Matrix([T[i][1] : i in [1..6]]);
N := Matrix([T[i][2] : i in [1..6]]);
C7 := Z2Z4AdditiveCode(M, N);

D := UniverseCode(GF(2), 2);
E := UniverseCode(Z4, 4);
C8 := Z2Z4AdditiveCode(D, E);

expectedOutputCode := LinearCode<Z4, 6 | L2>;
expectedAlpha := 2;

assert C1`Code eq expectedOutputCode;
assert C1`Alpha eq expectedAlpha; 
assert C2`Code eq expectedOutputCode;
assert C2`Alpha eq expectedAlpha;
assert C3`Code eq expectedOutputCode;
assert C3`Alpha eq expectedAlpha;
assert C4`Code eq expectedOutputCode;
assert C4`Alpha eq expectedAlpha;
assert C5`Code eq expectedOutputCode;
assert C5`Alpha eq expectedAlpha;
assert C1 eq C2;
assert C2 eq C3;
assert C3 eq C4;
assert C4 eq C5;
assert C5 eq C6;
assert C6 eq C7;
assert C7 eq C8;
assert IsSeparable(C1);

/************************************************************/
print "test 3a: Repetition Z2Z4-additive code
               alpha = 4, beta = 8, gamma = 1, delta = 0, kappa = 1, length = 12, #C = 2";
R := RSpace(Z4, 12);
L1 := [R![1,1,1,1,2,2,2,2,2,2,2,2]];
C1 := Z2Z4AdditiveCode(L1, 4);
L2 := [R![2,2,2,2,2,2,2,2,2,2,2,2]];
C2 := Z2Z4AdditiveCode(L2, 4);
C3 := Z2Z4AdditiveCode(sub<R|L2>, 4);
C4 := Z2Z4AdditiveCode(Matrix(L2), 4);
C5 := Z2Z4AdditiveRepetitionCode(4, 8);

T := FromZ4toZ2Z4(L2, 4);
C6 := Z2Z4AdditiveCode(T);

M := Matrix(T[1][1]);
N := Matrix(T[1][2]);
C7 := Z2Z4AdditiveCode(M, N);

expectedOutputCode := LinearCode<Z4, 12 | L2>;
expectedAlpha := 4;

assert C1`Code eq expectedOutputCode;
assert C1`Alpha eq expectedAlpha; 
assert C2`Code eq expectedOutputCode;
assert C2`Alpha eq expectedAlpha;
assert C3`Code eq expectedOutputCode;
assert C3`Alpha eq expectedAlpha;
assert C4`Code eq expectedOutputCode;
assert C4`Alpha eq expectedAlpha;
assert C5`Code eq expectedOutputCode;
assert C5`Alpha eq expectedAlpha;
assert C1 eq C2;
assert C2 eq C3; 
assert C3 eq C4;
assert C4 eq C5; 
assert C5 eq C6;
assert C6 eq C7;
assert not IsSeparable(C1);    

/************************************************************/
print "test 3b: Even weight Z2Z4-additive code
               alpha = 2, beta = 5, gamma = 1, delta = 5, kappa = 1, length = 7, #C = 2048";
R := RSpace(Z4, 7);
L1 := Basis(R);
L1[7][7] := 2;
for i in [1..6] do L1[i][7] := 1; end for;
C1 := Z2Z4AdditiveCode(L1, 2);
L2 := L1; L2[1][1] := 2; L2[2][2] := 2;
C2 := Z2Z4AdditiveCode(L2, 2);
C3 := Z2Z4AdditiveCode(sub<R|L2>, 2);
C4 := Z2Z4AdditiveCode(Matrix(L2), 2);
C5 := Z2Z4AdditiveEvenWeightCode(2, 5);

T := FromZ4toZ2Z4(L2, 2);
C6 := Z2Z4AdditiveCode(T);

M := Matrix([T[i][1] : i in [1..7]]);
N := Matrix([T[i][2] : i in [1..7]]);
C7 := Z2Z4AdditiveCode(M, N);

expectedOutputCode := LinearCode<Z4, 7 | L2>;
expectedAlpha := 2;

assert C1`Code eq expectedOutputCode;
assert C1`Alpha eq expectedAlpha;
assert C2`Code eq expectedOutputCode;
assert C2`Alpha eq expectedAlpha;
assert C3`Code eq expectedOutputCode;
assert C3`Alpha eq expectedAlpha;
assert C4`Code eq expectedOutputCode;
assert C4`Alpha eq expectedAlpha;
assert C5`Code eq expectedOutputCode;
assert C5`Alpha eq expectedAlpha;
assert C1 eq C2;
assert C2 eq C3;
assert C3 eq C4;
assert C4 eq C5; 
assert C5 eq C6;
assert C6 eq C7;
assert not IsSeparable(C1);           
 
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
C1 := Z2Z4AdditiveCode(L, 0);
C2 := Z2Z4AdditiveCode(L, 0);
C3 := Z2Z4AdditiveCode(sub<R|L>, 0);
C4 := Z2Z4AdditiveCode(Matrix(L), 0);
C5 := Z2Z4AdditiveCode(LinearCode(Matrix(L)));

 T := FromZ4toZ2Z4(L, 0);
 C6 := Z2Z4AdditiveCode(T);
 
M := Matrix([T[i][1] : i in [1..9]]);
N := Matrix([T[i][2] : i in [1..9]]);
C7 := Z2Z4AdditiveCode(M, N);

expectedOutputCode := LinearCode<Z4, 19 | L>;
expectedAlpha := 0;

assert C1`Code eq expectedOutputCode;
assert C1`Alpha eq expectedAlpha; 
assert C2`Code eq expectedOutputCode;
assert C2`Alpha eq expectedAlpha;
assert C3`Code eq expectedOutputCode;
assert C3`Alpha eq expectedAlpha;
assert C4`Code eq expectedOutputCode;
assert C4`Alpha eq expectedAlpha;
assert C5`Code eq expectedOutputCode;
assert C5`Alpha eq expectedAlpha;
assert C1 eq C2;
assert C2 eq C3;
assert C3 eq C4;   
assert C4 eq C5;
assert C5 eq C6;
assert C6 eq C7;
assert IsSeparable(C1);               

/************************************************************/
print "test 5: A Z2Z4-additive code with beta=0
               alpha = 10, beta = 0, gamma = 4, delta = 0, kappa = 4, length = 10, #C = 16";
R := RSpace(Z4, 10);
L1 := [R![0,0,1,0,1,0,1,0,1,0],
       R![1,0,0,0,0,0,1,1,1,1],
       R![0,0,1,0,0,0,0,1,0,0],
       R![0,1,0,1,0,1,0,1,1,1]];
C1 := Z2Z4AdditiveCode(L1, 10);
L2 := [R![0,0,2,0,2,0,2,0,2,0],
       R![2,0,0,0,0,0,2,2,2,2],
       R![0,0,2,0,0,0,0,2,0,0],
       R![0,2,0,2,0,2,0,2,2,2]];
C2 := Z2Z4AdditiveCode(L2, 10);
C3 := Z2Z4AdditiveCode(sub<R|L2>, 10);
C4 := Z2Z4AdditiveCode(Matrix(L2), 10);
C5 := Z2Z4AdditiveCode(LinearCode(ChangeRing( Matrix(L1), GF(2) )));

T := FromZ4toZ2Z4(L2, 10);
C6 := Z2Z4AdditiveCode(T);

M := Matrix([T[i][1] : i in [1..4]]);
N := Matrix([T[i][2] : i in [1..4]]);
C7 := Z2Z4AdditiveCode(M, N);

expectedOutputCode := LinearCode<Z4, 10 | L2>;
expectedAlpha := 10;

assert C1`Code eq expectedOutputCode;
assert C1`Alpha eq expectedAlpha; 
assert C2`Code eq expectedOutputCode;
assert C2`Alpha eq expectedAlpha;
assert C3`Code eq expectedOutputCode;
assert C3`Alpha eq expectedAlpha;
assert C4`Code eq expectedOutputCode;
assert C4`Alpha eq expectedAlpha;
assert C5`Code eq expectedOutputCode;
assert C5`Alpha eq expectedAlpha;
assert C1 eq C2;
assert C2 eq C3;
assert C3 eq C4;   
assert C4 eq C5;
assert C5 eq C6;
assert C6 eq C7;
assert IsSeparable(C1);

/************************************************************/
print "test 6: A Z2Z4-additive code with gamma=0
               alpha = 10, beta = 19, gamma = 0, delta = 7, kappa = 0, length = 19, #C = 16384";
R := RSpace(Z4, 29);
L1 := [R![1,1,0,0,1,1,0,1,1,1,0,0,0,1,0,0,0,0,2,3,2,1,1,3,3,1,2,3,0],
       R![0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,2,0,1,0,3,0,2,2,1,3,1,1,0,1],
       R![0,0,1,1,0,0,1,1,1,1,0,1,0,0,0,0,0,0,2,3,1,1,0,2,1,0,0,0,1],
       R![1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1,0,3,3,1,3,1,0,0,3,1,0,2],
       R![1,1,1,1,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,2,3,0,0,1,1,0,2,2,1],
       R![1,0,1,1,1,0,1,1,0,1,0,0,1,0,0,3,0,0,0,1,2,3,1,3,0,2,1,3,2],
       R![0,1,1,0,1,0,1,1,0,0,0,0,0,0,1,1,0,0,0,3,3,2,0,2,0,3,1,1,3]];
C1 := Z2Z4AdditiveCode(L1, 10);
L2 := [R![2,2,0,0,2,2,0,2,2,2,0,0,0,1,0,0,0,0,2,3,2,1,1,3,3,1,2,3,0],
       R![0,2,0,0,0,0,2,0,0,0,0,0,0,0,0,2,0,1,0,3,0,2,2,1,3,1,1,0,1],
       R![0,0,2,2,0,0,2,2,2,2,0,1,0,0,0,0,0,0,2,3,1,1,0,2,1,0,0,0,1],
       R![2,0,0,0,2,0,2,0,0,0,0,0,0,0,0,0,1,0,3,3,1,3,1,0,0,3,1,0,2],
       R![2,2,2,2,2,0,0,0,0,0,1,0,0,0,0,1,0,0,0,2,3,0,0,1,1,0,2,2,1],
       R![2,0,2,2,2,0,2,2,0,2,0,0,1,0,0,3,0,0,0,1,2,3,1,3,0,2,1,3,2],
       R![0,2,2,0,2,0,2,2,0,0,0,0,0,0,1,1,0,0,0,3,3,2,0,2,0,3,1,1,3]];
C2 := Z2Z4AdditiveCode(L2, 10);
C3 := Z2Z4AdditiveCode(sub<R|L2>, 10);
C4 := Z2Z4AdditiveCode(Matrix(L2), 10);
C5 := Z2Z4AdditiveCode(LinearCode(Matrix(L2)), 10);

T := FromZ4toZ2Z4(L2, 10);
C6 := Z2Z4AdditiveCode(T);

M := Matrix([T[i][1] : i in [1..7]]);
N := Matrix([T[i][2] : i in [1..7]]);
C7 := Z2Z4AdditiveCode(M, N);

expectedOutputCode := LinearCode<Z4, 29 | L2>;
expectedAlpha := 10;

assert C1`Code eq expectedOutputCode;
assert C1`Alpha eq expectedAlpha; 
assert C2`Code eq expectedOutputCode;
assert C2`Alpha eq expectedAlpha;
assert C3`Code eq expectedOutputCode;
assert C3`Alpha eq expectedAlpha;
assert C4`Code eq expectedOutputCode;
assert C4`Alpha eq expectedAlpha;
assert C5`Code eq expectedOutputCode;
assert C5`Alpha eq expectedAlpha;
assert C1 eq C2;
assert C2 eq C3;
assert C3 eq C4;   
assert C4 eq C5;
assert C5 eq C6;
assert C6 eq C7;
assert not IsSeparable(C1);

/************************************************************/
print "test 7: A Z2Z4-additive code with delta=0
               alpha = 10, beta = 19, gamma = 7, delta = 0, kappa = 4, length = 19, #C = 128";
R := RSpace(Z4, 29);
L1 := [R![1,0,0,0,0,0,1,1,1,0,0,0,0,0,0,2,2,0,0,0,0,0,0,2,0,2,0,2,2],
       R![0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,2,2,2,2],
       R![0,0,1,0,0,1,1,0,1,0,0,0,0,0,2,0,2,2,0,2,0,0,0,2,2,2,2,0,0],
       R![0,0,0,1,0,1,0,0,1,0,0,0,0,0,2,0,0,0,2,0,2,2,0,0,2,0,2,0,2],
       R![0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,0,2,0,2,0,0,2,0,2,0,0,0,0,2],
       R![0,0,0,0,0,0,0,0,0,0,0,2,0,2,0,2,2,2,0,2,0,2,0,2,2,2,2,0,2],
       R![0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,2,2,0,0,0,0,0,2,2,0,2,2]];      
C1 := Z2Z4AdditiveCode(L1, 10);
L2 := [R![2,0,0,0,0,0,2,2,2,0,0,0,0,0,0,2,2,0,0,0,0,0,0,2,0,2,0,2,2],
       R![0,2,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,2,2,2,2],
       R![0,0,2,0,0,2,2,0,2,0,0,0,0,0,2,0,2,2,0,2,0,0,0,2,2,2,2,0,0],
       R![0,0,0,2,0,2,0,0,2,0,0,0,0,0,2,0,0,0,2,0,2,2,0,0,2,0,2,0,2],
       R![0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,0,2,0,2,0,0,2,0,2,0,0,0,0,2],
       R![0,0,0,0,0,0,0,0,0,0,0,2,0,2,0,2,2,2,0,2,0,2,0,2,2,2,2,0,2],
       R![0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,2,2,0,0,0,0,0,2,2,0,2,2]];
C2 := Z2Z4AdditiveCode(L2, 10);
C3 := Z2Z4AdditiveCode(sub<R|L2>, 10);
C4 := Z2Z4AdditiveCode(Matrix(L2), 10);
C5 := Z2Z4AdditiveCode(LinearCode(Matrix(L2)), 10);

T := FromZ4toZ2Z4(L2, 10);
C6 := Z2Z4AdditiveCode(T);

M := Matrix([T[i][1] : i in [1..7]]);
N := Matrix([T[i][2] : i in [1..7]]);
C7 := Z2Z4AdditiveCode(M, N);

expectedOutputCode := LinearCode<Z4, 29 | L2>;
expectedAlpha := 10;

assert C1`Code eq expectedOutputCode;
assert C1`Alpha eq expectedAlpha; 
assert C2`Code eq expectedOutputCode;
assert C2`Alpha eq expectedAlpha;
assert C3`Code eq expectedOutputCode;
assert C3`Alpha eq expectedAlpha;
assert C4`Code eq expectedOutputCode;
assert C4`Alpha eq expectedAlpha;
assert C5`Code eq expectedOutputCode;
assert C5`Alpha eq expectedAlpha;
assert C1 eq C2;
assert C2 eq C3;
assert C3 eq C4;   
assert C4 eq C5;
assert C5 eq C6;
assert C6 eq C7;
assert not IsSeparable(C1);                       

/************************************************************/
print "test 8: A Z2Z4-additive code with kappa=0
               alpha = 10, beta = 19, gamma = 7, delta = 2, kappa = 0, length = 29, #C = 2048";
R := RSpace(Z4, 29);
L1 := [R![0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,0,2,0,0,0,0,0,2,2,2,2,0],
       R![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,2,0,2,0,2,0,2,2,2,2,0],
       R![0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,0,0,2,0,0,0,0,2,0,2,2,0,0],
       R![0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,2,0,2,0,2,0,0,2,2,2,0,2,0],
       R![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,2,0,0],
       R![0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,2,2,2,2,2,0,0,0,0,0,2,2,0],
       R![0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2,0,0,0,0,2,0,2,2],
       R![0,0,1,0,0,0,1,0,0,0,0,0,0,1,2,0,2,1,0,3,0,1,1,0,2,0,2,0,1],
       R![0,1,0,1,0,1,0,1,1,1,1,1,1,0,0,1,2,0,2,3,0,1,1,2,0,2,2,3,1]];
C1 := Z2Z4AdditiveCode(L1, 10);
L2 := [R![0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,0,2,0,0,0,0,0,2,2,2,2,0],
       R![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,2,0,2,0,2,0,2,2,2,2,0],
       R![0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,0,0,2,0,0,0,0,2,0,2,2,0,0],
       R![0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,2,0,2,0,2,0,0,2,2,2,0,2,0],
       R![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,2,0,0],
       R![0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,2,2,2,2,2,0,0,0,0,0,2,2,0],
       R![0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2,0,0,0,0,2,0,2,2],
       R![0,0,2,0,0,0,2,0,0,0,0,0,0,1,2,0,2,1,0,3,0,1,1,0,2,0,2,0,1],
       R![0,2,0,2,0,2,0,2,2,2,1,1,1,0,0,1,2,0,2,3,0,1,1,2,0,2,2,3,1]];
C2 := Z2Z4AdditiveCode(L2, 10);
C3 := Z2Z4AdditiveCode(sub<R|L2>, 10);
C4 := Z2Z4AdditiveCode(Matrix(L2), 10);
C5 := Z2Z4AdditiveCode(LinearCode(Matrix(L2)), 10);

T := FromZ4toZ2Z4(L2, 10);
C6 := Z2Z4AdditiveCode(T);

M := Matrix([T[i][1] : i in [1..9]]);
N := Matrix([T[i][2] : i in [1..9]]);
C7 := Z2Z4AdditiveCode(M, N);

expectedOutputCode := LinearCode<Z4, 29 | L2>;
expectedAlpha := 10;

assert C1`Code eq expectedOutputCode;
assert C1`Alpha eq expectedAlpha; 
assert C2`Code eq expectedOutputCode;
assert C2`Alpha eq expectedAlpha;
assert C3`Code eq expectedOutputCode;
assert C3`Alpha eq expectedAlpha;
assert C4`Code eq expectedOutputCode;
assert C4`Alpha eq expectedAlpha;
assert C5`Code eq expectedOutputCode;
assert C5`Alpha eq expectedAlpha;
assert C1 eq C2;
assert C2 eq C3;
assert C3 eq C4;   
assert C4 eq C5;
assert C5 eq C6;
assert C6 eq C7;
assert not IsSeparable(C1);            

/************************************************************/
print "test 9: A Z2Z4-additive code with kappa=alpha
               alpha = 7, beta = 12, gamma = 8, delta = 2, kappa = 7, length = 19, #C = 4096";
R := RSpace(Z4, 19);
L1 := [R![0,0,0,1,0,1,1,0,0,0,2,2,0,0,2,2,0,0,0],
       R![0,0,0,1,0,0,0,0,0,0,2,0,0,0,0,2,2,0,0],
       R![0,0,1,0,1,1,1,0,0,0,2,2,0,0,0,0,0,0,0],
       R![0,0,0,1,0,0,0,0,0,2,2,2,0,0,0,2,0,0,0],
       R![1,0,0,0,1,1,0,0,0,0,2,0,0,0,0,0,0,0,0],
       R![0,1,0,0,1,1,0,0,0,0,0,2,0,2,0,0,0,0,0],
       R![0,1,0,1,1,0,0,0,2,0,2,2,0,0,0,0,0,0,2],
       R![0,1,0,1,0,0,0,0,0,0,2,2,0,0,0,0,0,2,2],
       R![0,1,0,0,0,0,1,0,0,0,0,0,1,1,1,2,1,1,2],
       R![0,1,0,1,0,0,1,1,1,1,2,0,0,1,0,0,1,0,3]]; 
C1 := Z2Z4AdditiveCode(L1, 7);
L2 := [R![0,0,0,2,0,2,2,0,0,0,2,2,0,0,2,2,0,0,0],
       R![0,0,0,2,0,0,0,0,0,0,2,0,0,0,0,2,2,0,0],
       R![0,0,2,0,2,2,2,0,0,0,2,2,0,0,0,0,0,0,0],
       R![0,0,0,2,0,0,0,0,0,2,2,2,0,0,0,2,0,0,0],
       R![2,0,0,0,2,2,0,0,0,0,2,0,0,0,0,0,0,0,0],
       R![0,2,0,0,2,2,0,0,0,0,0,2,0,2,0,0,0,0,0],
       R![0,2,0,2,2,0,0,0,2,0,2,2,0,0,0,0,0,0,2],
       R![0,2,0,2,0,0,0,0,0,0,2,2,0,0,0,0,0,2,2],
       R![0,2,0,0,0,0,2,0,0,0,0,0,1,1,1,2,1,1,2],
       R![0,2,0,2,0,0,2,1,1,1,2,0,0,1,0,0,1,0,3]]; 
C2 := Z2Z4AdditiveCode(L2, 7);
C3 := Z2Z4AdditiveCode(sub<R|L2>, 7);
C4 := Z2Z4AdditiveCode(Matrix(L2), 7);
C5 := Z2Z4AdditiveCode(LinearCode(Matrix(L2)), 7);

T := FromZ4toZ2Z4(L2, 7);
C6 := Z2Z4AdditiveCode(T);

M := Matrix([T[i][1] : i in [1..10]]);
N := Matrix([T[i][2] : i in [1..10]]);
C7 := Z2Z4AdditiveCode(M, N);

expectedOutputCode := LinearCode<Z4, 19 | L2>;
expectedAlpha := 7;

assert C1`Code eq expectedOutputCode;
assert C1`Alpha eq expectedAlpha; 
assert C2`Code eq expectedOutputCode;
assert C2`Alpha eq expectedAlpha;
assert C3`Code eq expectedOutputCode;
assert C3`Alpha eq expectedAlpha;
assert C4`Code eq expectedOutputCode;
assert C4`Alpha eq expectedAlpha;
assert C5`Code eq expectedOutputCode;
assert C5`Alpha eq expectedAlpha;
assert C1 eq C2;
assert C2 eq C3;
assert C3 eq C4;   
assert C4 eq C5;
assert C5 eq C6;
assert C6 eq C7;    
assert not IsSeparable(C1);        
                           
/************************************************************/
print "test 10: A Z2Z4-additive code with kappa<>gamma
                alpha = 9, beta = 12, gamma = 6, delta = 3, kappa = 4, length = 21, #C = 4096";
R := RSpace(Z4, 21);
L1 := [R![1,0,0,0,1,1,1,0,1,0,0,0,0,2,0,0,0,0,0,0,2],
       R![1,0,0,1,0,0,1,0,1,0,0,2,0,2,0,0,2,0,0,0,0],
       R![0,0,1,1,0,1,1,0,0,0,0,0,0,2,2,0,0,0,2,0,0],
       R![0,0,1,1,0,1,1,0,0,0,0,0,0,0,0,2,2,0,2,0,0],
       R![0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
       R![1,0,0,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
       R![0,0,0,0,0,1,1,1,1,0,0,0,1,0,1,1,0,0,1,0,1],
       R![0,0,1,0,1,1,0,0,1,1,0,1,0,2,0,1,0,0,0,0,0],
       R![0,0,0,1,0,1,1,0,1,0,0,0,1,0,1,1,0,0,1,0,0]];
C1 := Z2Z4AdditiveCode(L1, 9);
L2 := [R![2,0,0,0,2,2,2,0,2,0,0,0,0,2,0,0,0,0,0,0,2],
       R![2,0,0,2,0,0,2,0,2,0,0,2,0,2,0,0,2,0,0,0,0],
       R![0,0,2,2,0,2,2,0,0,0,0,0,0,2,2,0,0,0,2,0,0],
       R![0,0,2,2,0,2,2,0,0,0,0,0,0,0,0,2,2,0,2,0,0],
       R![0,0,0,0,2,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0],
       R![2,0,0,2,0,0,2,0,2,0,0,0,0,0,0,0,0,0,0,0,0],
       R![0,0,0,0,0,2,2,2,2,0,0,0,1,0,1,1,0,0,1,0,1],
       R![0,0,2,0,2,2,0,0,2,1,0,1,0,2,0,1,0,0,0,0,0],
       R![0,0,0,2,0,2,2,0,2,0,0,0,1,0,1,1,0,0,1,0,0]];
C2 := Z2Z4AdditiveCode(L2, 9);
C3 := Z2Z4AdditiveCode(sub<R|L2>, 9);
C4 := Z2Z4AdditiveCode(Matrix(L2), 9);
C5 := Z2Z4AdditiveCode(LinearCode(Matrix(L2)), 9);

T := FromZ4toZ2Z4(L2, 9);
C6 := Z2Z4AdditiveCode(T);

M := Matrix([T[i][1] : i in [1..9]]);
N := Matrix([T[i][2] : i in [1..9]]);
C7 := Z2Z4AdditiveCode(M, N);

expectedOutputCode := LinearCode<Z4, 21 | L2>;
expectedAlpha := 9;

assert C1`Code eq expectedOutputCode;
assert C1`Alpha eq expectedAlpha; 
assert C2`Code eq expectedOutputCode;
assert C2`Alpha eq expectedAlpha;
assert C3`Code eq expectedOutputCode;
assert C3`Alpha eq expectedAlpha;
assert C4`Code eq expectedOutputCode;
assert C4`Alpha eq expectedAlpha;
assert C5`Code eq expectedOutputCode;
assert C5`Alpha eq expectedAlpha;
assert C1 eq C2;
assert C2 eq C3;
assert C3 eq C4;   
assert C4 eq C5;
assert C5 eq C6;
assert C6 eq C7;
assert not IsSeparable(C1);            

/************************************************************/
print "test 11: A Z2Z4-additive Hadamard and extended perfect code with alpha=0";

m := 5;
delta := 3;
C, G1 := Z2Z4HadamardCode(delta, m : OverZ4 := true);
D, G2 := Z2Z4ExtendedPerfectCode(delta, m : OverZ4 := true);
 
expectedAlpha := 0;
expectedGamma := m+1-2*delta;
expectedDelta := delta;
expectedBinaryLength := 2^m;
expectedCardinal := 2*expectedBinaryLength; 
expectedMinimumWeight := expectedBinaryLength/2;
expectedWeightDistribution := [<0, 1>, <2^(m-1), 2^(m+1)-2>, <2^m, 1>];

assert C`Code eq HadamardCodeZ4(delta, m); 
assert C`Code eq LinearCode(G1);
assert C`Alpha eq expectedAlpha;
assert Z2Z4Type(C)[3] eq expectedGamma;
assert Z2Z4Type(C)[4] eq expectedDelta;
assert LeeWeightDistribution(C) eq expectedWeightDistribution;
assert BinaryLength(C) eq expectedBinaryLength;
assert #C eq expectedCardinal;
assert Z2Z4MinimumLeeWeight(C) eq expectedMinimumWeight; 

assert D eq Dual(C);
assert G1 eq G2;
assert D`Code eq Dual(LinearCode(G2));
assert D`Alpha eq expectedAlpha;
assert BinaryLength(D) eq expectedBinaryLength;
assert #D eq 2^(expectedBinaryLength-m-1);
assert Z2Z4MinimumLeeWeight(D) eq 4;
assert not IsSeparable(C1);

/************************************************************/
print "test 12: A Z2Z4-additive Hadamard and extended perfect code with alpha<>0";

m := 4;
delta := 2;
C, G1 := Z2Z4HadamardCode(delta, m);
D, G2 := Z2Z4ExtendedPerfectCode(delta, m);

expectedAlpha := 4;
expectedGamma := m+1-2*delta;
expectedDelta := delta;
expectedBinaryLength := 2^m;
expectedCardinal := 2*expectedBinaryLength; 
expectedMinimumWeight := expectedBinaryLength/2;
expectedWeightDistribution := [<0, 1>, <2^(m-1), 2^(m+1)-2>, <2^m, 1>];

assert C`Code eq LinearCode(G1);
assert C`Alpha eq expectedAlpha;
assert Z2Z4Type(C)[3] eq expectedGamma;
assert Z2Z4Type(C)[4] eq expectedDelta;
assert LeeWeightDistribution(C) eq expectedWeightDistribution;
assert BinaryLength(C) eq expectedBinaryLength;
assert #C eq expectedCardinal;
assert Z2Z4MinimumLeeWeight(C) eq expectedMinimumWeight; 

assert D eq Dual(C);
assert G1 eq G2;
assert D`Alpha eq expectedAlpha;
assert BinaryLength(D) eq expectedBinaryLength;
assert #D eq 2^(expectedBinaryLength-m-1);
assert Z2Z4MinimumLeeWeight(D) eq 4;
assert not IsSeparable(C1);

/************************************************************/
print "test 13: A random Z2Z4-additive code with alpha=0";

alpha := 0;
beta := 5;
C := RandomZ2Z4AdditiveCode(alpha, beta);

assert IsZ2Z4AdditiveCode(C);
assert Z2Z4Type(C)[1] eq alpha;
assert Z2Z4Type(C)[2] eq beta;

/************************************************************/
print "test 14: A random Z2Z4-additive code with beta=0";

alpha := 4;
beta := 0;
C := RandomZ2Z4AdditiveCode(alpha, beta);

assert IsZ2Z4AdditiveCode(C);
assert Z2Z4Type(C)[1] eq alpha;
assert Z2Z4Type(C)[2] eq beta;

/************************************************************/
print "test 15: A random Z2Z4-additive code with gamma=0";

alpha := 4;
beta := 10;
gamma := 0;
delta := 6;
C := RandomZ2Z4AdditiveCode(alpha, beta, gamma, delta);

assert IsZ2Z4AdditiveCode(C);
assert Z2Z4Type(C)[1] eq alpha;
assert Z2Z4Type(C)[2] eq beta;
assert Z2Z4Type(C)[3] eq gamma;
assert Z2Z4Type(C)[4] eq delta;

/************************************************************/
print "test 16: A random Z2Z4-additive code with delta=0";

alpha := 4;
beta := 10;
gamma := 6;
delta := 0;
C := RandomZ2Z4AdditiveCode(alpha, beta, gamma, delta);

assert IsZ2Z4AdditiveCode(C);
assert Z2Z4Type(C)[1] eq alpha;
assert Z2Z4Type(C)[2] eq beta;
assert Z2Z4Type(C)[3] eq gamma;
assert Z2Z4Type(C)[4] eq delta;

/************************************************************/
print "test 17: A random Z2Z4-additive code with kappa=0";

alpha := 4;
beta := 10;
gamma := 6;
delta := 3;
kappa := 0;
C := RandomZ2Z4AdditiveCode(alpha, beta, gamma, delta, kappa);

assert IsZ2Z4AdditiveCode(C);
assert Z2Z4Type(C)[1] eq alpha;
assert Z2Z4Type(C)[2] eq beta;
assert Z2Z4Type(C)[3] eq gamma;
assert Z2Z4Type(C)[4] eq delta;
assert Z2Z4Type(C)[5] eq kappa;

/************************************************************/
print "test 18: A random Z2Z4-additive code given alpha and beta";

alpha := 3;
beta := 10;
C := RandomZ2Z4AdditiveCode(alpha, beta);

assert IsZ2Z4AdditiveCode(C);
assert Z2Z4Type(C)[1] eq alpha;
assert Z2Z4Type(C)[2] eq beta;

/************************************************************/
print "test 19: A random Z2Z4-additive code given alpha, beta and gamma";

alpha := 5;
beta := 12;
gamma := 2;
C := RandomZ2Z4AdditiveCode(alpha, beta, gamma);

assert IsZ2Z4AdditiveCode(C);
assert Z2Z4Type(C)[1] eq alpha;
assert Z2Z4Type(C)[2] eq beta;
assert Z2Z4Type(C)[3] eq gamma;

/************************************************************/
print "test 20: A random Z2Z4-additive code given alpha, beta, gamma and delta";

alpha := 5;
beta := 12;
gamma := 2;
delta := 4;
C := RandomZ2Z4AdditiveCode(alpha, beta, gamma, delta);

assert IsZ2Z4AdditiveCode(C);
assert Z2Z4Type(C)[1] eq alpha;
assert Z2Z4Type(C)[2] eq beta;
assert Z2Z4Type(C)[3] eq gamma;
assert Z2Z4Type(C)[4] eq delta;

/************************************************************/
print "test 21: A random Z2Z4-additive code given alpha, beta, gamma, delta and kappa";

alpha := 5;
beta := 12;
gamma := 2;
delta := 4;
kappa := 2;
C := RandomZ2Z4AdditiveCode(alpha, beta, gamma, delta, kappa);

assert IsZ2Z4AdditiveCode(C);
assert Z2Z4Type(C)[1] eq alpha;
assert Z2Z4Type(C)[2] eq beta;
assert Z2Z4Type(C)[3] eq gamma;
assert Z2Z4Type(C)[4] eq delta;
assert Z2Z4Type(C)[5] eq kappa;

/************************************************************/
print "test 22: A random Z2Z4-additive code given alpha, beta, gamma, delta and kappa";

alpha := 5;
beta := 12;
gamma := 2;
delta := 4;
kappa := 1;
C := RandomZ2Z4AdditiveCode(alpha, beta, gamma, delta, kappa);

assert IsZ2Z4AdditiveCode(C);
assert Z2Z4Type(C)[1] eq alpha;
assert Z2Z4Type(C)[2] eq beta;
assert Z2Z4Type(C)[3] eq gamma;
assert Z2Z4Type(C)[4] eq delta;
assert Z2Z4Type(C)[5] eq kappa;

/************************************************************/
print "test 23: A random Z2Z4-additive code given alpha, beta, gamma, delta and kappa=alpha";

alpha := 5;
beta := 12;
gamma := 6;
delta := 4;
kappa := 5;
C := RandomZ2Z4AdditiveCode(alpha, beta, gamma, delta, kappa);

assert IsZ2Z4AdditiveCode(C);
assert Z2Z4Type(C)[1] eq alpha;
assert Z2Z4Type(C)[2] eq beta;
assert Z2Z4Type(C)[3] eq gamma;
assert Z2Z4Type(C)[4] eq delta;
assert Z2Z4Type(C)[5] eq kappa;

/************************************************************/
print "test 24: A random Z2Z4-additive code given alpha, beta, gamma, delta and kappa=gamma";

alpha := 5;
beta := 12;
gamma := 2;
delta := 4;
kappa := 2;
C := RandomZ2Z4AdditiveCode(alpha, beta, gamma, delta, kappa);

assert IsZ2Z4AdditiveCode(C);
assert Z2Z4Type(C)[1] eq alpha;
assert Z2Z4Type(C)[2] eq beta;
assert Z2Z4Type(C)[3] eq gamma;
assert Z2Z4Type(C)[4] eq delta;
assert Z2Z4Type(C)[5] eq kappa;

/************************************************************/
print "test 25: A random Z2Z4-additive code given alpha, beta, gamma, delta and kappa, 
                with gamma+delta = beta+kappa";

alpha := 5;
beta := 8;
gamma := 3;
delta := 7;
kappa := 2;
C := RandomZ2Z4AdditiveCode(alpha, beta, gamma, delta, kappa);

assert IsZ2Z4AdditiveCode(C);
assert Z2Z4Type(C)[1] eq alpha;
assert Z2Z4Type(C)[2] eq beta;
assert Z2Z4Type(C)[3] eq gamma;
assert Z2Z4Type(C)[4] eq delta;
assert Z2Z4Type(C)[5] eq kappa;

/************************************************************/
print "test 26: A Z2Z4-additive code 
               alpha = 2, beta = 1, gamma = 2, delta = 0, kappa = 1, length = 3, #C = 4";


D := LinearCode<GF(2), 2 | [1,1]>;
E := LinearCode<Z4, 1| [2]>;
R := RSpace(Z4, 3);
L := [R![2, 2, 0],
      R![0, 0, 2],
      R![2, 2, 2],
      R![0, 0, 0]];

expectedOutputCode := LinearCode(Matrix(L));
expectedOutputAlpha := 2;

C := Z2Z4AdditiveCode(D, E);

assert IsZ2Z4AdditiveCode(C);
assert C`Code eq expectedOutputCode;
assert C`Alpha eq expectedOutputAlpha;
assert IsSeparable(C);

/************************************************************/
print "test 27: A Z2Z4-additive code 
               alpha = 8, beta = 4, gamma = 6, delta = 1, kappa = 4, length = 12, #C = 256";
  
D := ExtendCode(HammingCode(GF(2),3));
E := LinearCode<Z4, 4 | [[2,2,0,0],[2,0,2,0],[1,1,1,1]]>;
R := RSpace(Z4, 12);
L := [R![2, 0, 0, 0, 2, 2, 0, 2, 0, 0, 0, 0],
      R![0, 2, 0, 0, 0, 2, 2, 2, 0, 0, 0, 0],
      R![0, 0, 2, 0, 2, 2, 2, 0, 0, 0, 0, 0],
      R![0, 0, 0, 2, 2, 0, 2, 2, 0, 0, 0, 0],
      R![0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
      R![0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2],
      R![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2]];

expectedOutputCode := LinearCode(Matrix(L));
expectedOutputAlpha := 8;

C := Z2Z4AdditiveCode(D, E);

assert IsZ2Z4AdditiveCode(C);
assert C`Code eq expectedOutputCode;
assert C`Alpha eq expectedOutputAlpha;
assert IsSeparable(C);
assert IsSelfOrthogonal(C);
assert IsSelfDual(C);