/*************************************************************/
/*                                                           */
/* Package or project name: Z2Z4AdditiveCodes package        */
/* Test file name: Z2Z4CoveringRadius_BB_test.m      	     */
/*                                                           */
/* Comments: Black-box tests for the intrinsic functions     */
/*        	 CoveringRadius, CoveringRadiusBounds and        */
/*           IsCompletelyRegular included in the             */
/*           Z2Z4CoveringRadius.m file                       */
/*                                                           */
/* Authors: Jaume Pujol and Mercè Villanueva                 */
/*                                                           */
/* Revision version and last date: v1.0   2019/03/28         */                  
/*                                                           */
/*************************************************************/

SetAssertions(true);
Alarm(30*60); 

/*************************************************************
	GLOBAL VARIABLES
*************************************************************/	
Z4 := Integers(4);


/****************************************************************/
/*                                                              */
/* Function name: CoveringRadiusBounds                          */
/* Parameters: C                                                */
/* Function description: Return a lower and upper bounds on the */
/*   covering radius of the Z2Z4-additive code C. The lower     */
/*   bound is given by the maximum between the error correcting */
/*   capability of C and the covering radius of the linear span */
/*   of Cbin, where Cbin=Phi(C) and Phi is the Gray map. The    */
/*   upper bound is given by the minimum between the external   */
/*   distance of C and the covering radius of Kbin=Phi(K_C),    */
/*   where K_C is the kernel of C. Note that this function is   */
/*   only applicable when C is small.                           */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/* Output parameters description:                               */
/*   - Two integers with a lower and upper bound for the        */
/*     covering radius of C                                     */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> RngIntElt, RngIntElt            */
/*                                                              */
/****************************************************************/
/****************************************************************/
/*                                                              */
/* Function name: CoveringRadius                                */
/* Parameters: C                                                */
/* Function description: The covering radius of the Z2Z4-       */
/*   additive code C, which  is the smallest radius rho         */
/*   (considering the Lee distance) such that the spheres of    */
/*   radius rho centered at all codewords of C cover the ambient*/
/*   space V=Z2^alpha x Z4^beta. Note that this function is only*/
/*   applicable when C is small. The parameter MaximumTime sets */
/*   a time limit (in seconds of ``user time'') after which the */
/*   calculation is aborted. The default value is infinity,     */
/*   when there is no restriction on time.                      */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/* Output parameters description:                               */
/*   - An integer with the covering radius of C                 */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> RngIntElt                       */
/*                                                              */
/****************************************************************/
/****************************************************************/
/*                                                              */
/* Function name: IsCompletelyRegular                           */
/* Parameters: C, MaximumCosetWeight, MaximumTime               */
/* Function description: Given a Z2Z4-additive code C of        */
/*    type (alpha, beta; gamma, delta; kappa), return true      */
/*    if and only if Cbin is completely regular, where          */
/*    Cbin=Phi(C) and Phi is the Gray map. If return true,      */
/*    then it also returns the intersection array, a            */
/*    sequence with the number of coset leaders of Cbin of      */
/*    each weight, and the covering radius. Note that, if       */
/*    no optional parameters are specified, then this           */
/*    function is only applicable when C is small.              */
/*    A code Cbin subset of Z2^(alpha+2beta) is                 */
/*    called completely regular if, for all i>=0, every         */
/*    vector x in C(i) has the same number c(i) of              */
/*    neighbours in C(i-1) and the same number b(i) of          */
/*    neighbours in C(i+1), where                               */
/*    C(i)={x in Z2^(alpha+2beta) | d(x,Cbin)=i}. The           */
/*    intersection array is given as a seq. with the            */
/*    parameters [ [b(0),...,b(rho-1)], [c(1),...,c(rho)] ],    */
/*    where rho is the covering radius of Cbin.                 */
/* Input parameters description:                                */
/*   - C:   Z2Z4-additive code                                  */
/*   - MaximumCosetWeight: Integer that sets the maximum        */
/*     weight of the coset leaders which are used in the        */
/*     calculation to check whether the code is completely      */
/*     regular or not. The default value is rho.                */
/*   - MaximumTime: Integer that sets a time limit (in          */
/*     seconds of "user time") after which the calculation      */
/*     is aborted. The default value is infinitive, when        */
/*     there is no restriction on time.                         */
/* Output parameters description:                               */
/*   - Whether the image of C under the Gray map is a           */
/*     binary completely regular code                           */
/*   - The intersection array                                   */
/*   - List with the number of coset leaders of each weight     */
/*   - The covering radius of the image of C under Gray map     */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> BoolElt, SeqEnum,               */
/*                              SeqEnum, RngIntElt              */
/*                                                              */
/****************************************************************/
print "THE COVERING RADIUS";

print "test 1: Trivial Z2Z4-additive zero code 
               alpha = 2, beta = 4, gamma = 0, delta = 0, kappa = 0, length = 6, #C = 1"; 
alpha := 2;
beta:= 4;
gamma := 0;
delta := 0;
C := Z2Z4AdditiveCode([RSpace(Z4, 6)!0], 2);

expectedExternalDistance := alpha + 2*beta;
expectedCoveringRadius := alpha + 2*beta;
expectedIsCompletelyRegular := true;
expectedIntersectionArray := [[10,9,8,7,6,5,4,3,2,1],[1..10]];
expectedCosetDistribution := [10, 45, 120, 210, 252, 210, 120, 45, 10, 1];  

extDist := ExternalDistance(C);
rho := CoveringRadius(C);
assert extDist eq expectedExternalDistance;
assert rho eq expectedCoveringRadius;
//rhoLowerBound, rhoUpperBound := CoveringRadiusBounds(C);
//assert (rhoLowerBound le rho) and (rho le rhoUpperBound);

isCR, InterArray, CosetDist, CR := IsCompletelyRegular(C);
assert isCR eq expectedIsCompletelyRegular;
assert InterArray eq expectedIntersectionArray;
assert CosetDist eq expectedCosetDistribution;
assert CR eq expectedCoveringRadius;

/************************************************************/
print "test 2: Trivial Z2Z4-additive universe code 
               alpha = 2, beta = 6, gamma = 2, delta = 6, kappa = 2, length = 8, #C = 16384";
alpha := 2;
beta:= 6;
gamma := 2;
delta := 6;
C := Z2Z4AdditiveUniverseCode(alpha, beta);

expectedExternalDistance := 0;
expectedCoveringRadius := 0;
expectedIsCompletelyRegular := true;
expectedIntersectionArray := [[],[]];
expectedCosetDistribution := [];  

extDist := ExternalDistance(C);
rho := CoveringRadius(C);
assert extDist eq expectedExternalDistance;
assert rho eq expectedCoveringRadius;
rhoLowerBound, rhoUpperBound := CoveringRadiusBounds(C);
assert (rhoLowerBound le rho) and (rho le rhoUpperBound);

isCR, InterArray, CosetDist, CR := IsCompletelyRegular(C);
assert isCR eq expectedIsCompletelyRegular;
assert InterArray eq expectedIntersectionArray;
assert CosetDist eq expectedCosetDistribution;
assert CR eq expectedCoveringRadius;

/************************************************************/
print "test 3: Repetition Z2Z4-additive code
               alpha = 1, beta = 3, gamma = 1, delta = 0, kappa = 1, length = 4, #C = 2";
alpha := 1;
beta:= 3;
gamma := 1;
delta := 0;
C := Z2Z4AdditiveRepetitionCode(alpha, beta);

expectedExternalDistance := 3;
expectedCoveringRadius := 3;
expectedIsCompletelyRegular := true;
expectedIntersectionArray := [[7,6,5], [1,2,3]];
expectedCosetDistribution := [7, 21, 35];  

extDist := ExternalDistance(C);
rho := CoveringRadius(C);
assert extDist eq expectedExternalDistance;
assert rho eq expectedCoveringRadius;
rhoLowerBound, rhoUpperBound := CoveringRadiusBounds(C);
assert (rhoLowerBound le rho) and (rho le rhoUpperBound);

isCR, InterArray, CosetDist, CR := IsCompletelyRegular(C);
assert isCR eq expectedIsCompletelyRegular;
assert InterArray eq expectedIntersectionArray;
assert CosetDist eq expectedCosetDistribution;
assert CR eq expectedCoveringRadius;

/************************************************************/
print "test 4: A Z2Z4-additive code from a quaternary Kerdock code 
               alpha = 0, beta = 8, gamma = 0, delta = 4, kappa = 0, length = 8, #C = 256";
C := Z2Z4AdditiveCode(KerdockCode(3));

expectedExternalDistance := 4;
expectedCoveringRadius := 4;
expectedIsCompletelyRegular := true;
expectedIntersectionArray := [[16,15,14,1], [1,2,15,16]];
expectedCosetDistribution := [16, 120, 112, 7];  

extDist := ExternalDistance(C);
rho := CoveringRadius(C);
assert extDist eq expectedExternalDistance;
assert rho eq expectedCoveringRadius;
rhoLowerBound, rhoUpperBound := CoveringRadiusBounds(C);
assert (rhoLowerBound le rho) and (rho le rhoUpperBound);

isCR, InterArray, CosetDist, CR := IsCompletelyRegular(C);
assert isCR eq expectedIsCompletelyRegular;
assert InterArray eq expectedIntersectionArray;
assert CosetDist eq expectedCosetDistribution;
assert CR eq expectedCoveringRadius;

/************************************************************/
print "test 5: A Z2Z4-additive code from a quaternary Preparata code 
               alpha = 0, beta = 8, gamma = 0, delta = 4, kappa = 0, length = 8, #C = 256";               
C := Z2Z4AdditiveCode(PreparataCode(4));

expectedExternalDistance := 4;
expectedCoveringRadius := 4;
expectedIsCompletelyRegular := false;

extDist := ExternalDistance(C);
rho := CoveringRadius(C);
assert extDist eq expectedExternalDistance;
assert rho eq expectedCoveringRadius;
rhoLowerBound, rhoUpperBound := CoveringRadiusBounds(C);
assert (rhoLowerBound le rho) and (rho le rhoUpperBound);

isCR, InterArray, CosetDist, CR := IsCompletelyRegular(C);
assert isCR eq expectedIsCompletelyRegular;

/************************************************************/
print "test 6: A Z2Z4-additive code from a quaternary Hadamard code 
               alpha = 0, beta = 8, gamma = 0, delta = 4, kappa = 0, length = 8, #C = 256";             
C := Z2Z4AdditiveCode(HadamardCodeZ4(2, 4));

expectedExternalDistance := 6;
expectedCoveringRadius := 6;
expectedIsCompletelyRegular := false;

extDist := ExternalDistance(C);
rho := CoveringRadius(C);
assert extDist eq expectedExternalDistance;
assert rho eq expectedCoveringRadius;
rhoLowerBound, rhoUpperBound := CoveringRadiusBounds(C);
assert (rhoLowerBound le rho) and (rho le rhoUpperBound);

isCR, InterArray, CosetDist, CR := IsCompletelyRegular(C);
assert isCR eq expectedIsCompletelyRegular;

/************************************************************/
print "test 7: A Z2Z4-additive code from a quaternary extended 1-perfect code 
               alpha = 0, beta = 8, gamma = 0, delta = 4, kappa = 0, length = 8, #C = 256";
C := Z2Z4AdditiveCode(Dual(HadamardCodeZ4(3, 5)));

expectedExternalDistance := 2;
expectedCoveringRadius := 2;
expectedIsCompletelyRegular := true;
expectedIntersectionArray := [[32, 31], [1,32]];
expectedCosetDistribution := [32, 31];  

extDist := ExternalDistance(C);
rho := CoveringRadius(C);
assert extDist eq expectedExternalDistance;
assert rho eq expectedCoveringRadius;
rhoLowerBound, rhoUpperBound := CoveringRadiusBounds(C);
assert (rhoLowerBound le rho) and (rho le rhoUpperBound);

isCR, InterArray, CosetDist, CR := IsCompletelyRegular(C);
assert isCR eq expectedIsCompletelyRegular;
assert InterArray eq expectedIntersectionArray;
assert CosetDist eq expectedCosetDistribution;
assert CR eq expectedCoveringRadius;

/************************************************************/
print "test 8: A Z2Z4-additive 1-perfect code 
               (9, 2048) Z2Z4-additive code of type (3, 6; 3, 4; 3)";
C := PunctureCode(Z2Z4ExtendedPerfectCode(2, 4), 1);

expectedExternalDistance := 1;
expectedCoveringRadius := 1;
expectedIsCompletelyRegular := true;
expectedIntersectionArray := [[15], [1]];
expectedCosetDistribution := [15];

extDist := ExternalDistance(C);
rho := CoveringRadius(C);
assert extDist eq expectedExternalDistance;
assert rho eq expectedCoveringRadius;
rhoLowerBound, rhoUpperBound := CoveringRadiusBounds(C);
assert (rhoLowerBound le rho) and (rho le rhoUpperBound);

isCR, InterArray, CosetDist, CR := IsCompletelyRegular(C);
assert isCR eq expectedIsCompletelyRegular;
assert InterArray eq expectedIntersectionArray;
assert CosetDist eq expectedCosetDistribution;
assert CR eq expectedCoveringRadius;

/************************************************************/
print "test 9: A Z2Z4-additive 1-perfect code 
               (9, 2048) Z2Z4-additive code of type (2, 6; 3, 4; 2)";
C := PunctureCode(Z2Z4ExtendedPerfectCode(2, 4), {1, 2});

expectedExternalDistance := 1;
expectedCoveringRadius := 1;
expectedIsCompletelyRegular := true;
expectedIntersectionArray := [[14], [2]];
expectedCosetDistribution := [7];

extDist := ExternalDistance(C);
rho := CoveringRadius(C);
assert extDist eq expectedExternalDistance;
assert rho eq expectedCoveringRadius;
rhoLowerBound, rhoUpperBound := CoveringRadiusBounds(C);
assert (rhoLowerBound le rho) and (rho le rhoUpperBound);

isCR, InterArray, CosetDist, CR := IsCompletelyRegular(C);
assert isCR eq expectedIsCompletelyRegular;
assert InterArray eq expectedIntersectionArray;
assert CosetDist eq expectedCosetDistribution;
assert CR eq expectedCoveringRadius;

/************************************************************/
print "test 10: A Z2Z4-additive 1-perfect code 
               (9, 2048) Z2Z4-additive code of type (2, 5; 2, 4; 2)";
C := PunctureCode(Z2Z4ExtendedPerfectCode(2, 4), {1, 2, 9});

expectedExternalDistance := 1;
expectedCoveringRadius := 1;
expectedIsCompletelyRegular := true;
expectedIntersectionArray := [[12], [4]];
expectedCosetDistribution := [3];

extDist := ExternalDistance(C);
rho := CoveringRadius(C);
assert extDist eq expectedExternalDistance;
assert rho eq expectedCoveringRadius;
rhoLowerBound, rhoUpperBound := CoveringRadiusBounds(C);
assert (rhoLowerBound le rho) and (rho le rhoUpperBound);

isCR, InterArray, CosetDist, CR := IsCompletelyRegular(C);
assert isCR eq expectedIsCompletelyRegular;
assert InterArray eq expectedIntersectionArray;
assert CosetDist eq expectedCosetDistribution;
assert CR eq expectedCoveringRadius;

/************************************************************/
print "test 11: A Z2Z4-additive 1-perfect code 
               (9, 2048) Z2Z4-additive code of type (3, 6; 3, 4; 3)";
C := ShortenCode(PunctureCode(Z2Z4ExtendedPerfectCode(2, 4), {1}), 2);

expectedExternalDistance := 2;
expectedCoveringRadius := 2;
expectedIsCompletelyRegular := true;
expectedIntersectionArray := [[14, 1], [1, 14]];
expectedCosetDistribution := [14, 1];

extDist := ExternalDistance(C);
rho := CoveringRadius(C);
assert extDist eq expectedExternalDistance;
assert rho eq expectedCoveringRadius;
rhoLowerBound, rhoUpperBound := CoveringRadiusBounds(C);
assert (rhoLowerBound le rho) and (rho le rhoUpperBound);

isCR, InterArray, CosetDist, CR := IsCompletelyRegular(C);
assert isCR eq expectedIsCompletelyRegular;
assert InterArray eq expectedIntersectionArray;
assert CosetDist eq expectedCosetDistribution;
assert CR eq expectedCoveringRadius;

/************************************************************/
print "test 12: A Z2Z4-additive 1-perfect code 
               (9, 2048) Z2Z4-additive code of type (3, 6; 3, 4; 3)";
C := ShortenCode(PunctureCode(Z2Z4ExtendedPerfectCode(2, 4), {1}), 6);

expectedExternalDistance := 3;
expectedCoveringRadius := 2;
expectedIsCompletelyRegular := false;

extDist := ExternalDistance(C);
rho := CoveringRadius(C);
assert extDist eq expectedExternalDistance;
assert rho eq expectedCoveringRadius;
rhoLowerBound, rhoUpperBound := CoveringRadiusBounds(C);
assert (rhoLowerBound le rho) and (rho le rhoUpperBound);

isCR, InterArray, CosetDist, CR := IsCompletelyRegular(C);
assert isCR eq expectedIsCompletelyRegular;

/************************************************************/
print "test 13: A Z2Z4-additive code from a quaternary Klemm code 
               alpha = 0, beta = 8, gamma = 6, delta = 1, kappa = 0, length = 8, #C = 256";
C := KlemmCodeZ4(2);

expectedExternalDistance := 4;
expectedCoveringRadius := 4;
expectedIsCompletelyRegular := false;

extDist := ExternalDistance(C);
rho := CoveringRadius(C);
assert extDist eq expectedExternalDistance;
assert rho eq expectedCoveringRadius;
rhoLowerBound, rhoUpperBound := CoveringRadiusBounds(C);
assert (rhoLowerBound le rho) and (rho le rhoUpperBound);

isCR, InterArray, CosetDist, CR := IsCompletelyRegular(C);
assert isCR eq expectedIsCompletelyRegular;

/************************************************************/
print "test 14: A Z2Z4-additive code from a quaternary Klemm code 
               alpha = 1, beta = 7, gamma = 6, delta = 1, kappa = 0, length = 8, #C = 256";
_, G := KlemmCodeZ4(2);
C := Z2Z4AdditiveCode(G, 1);

expectedExternalDistance := 7;
expectedCoveringRadius := 4;
expectedIsCompletelyRegular := false;

extDist := ExternalDistance(C);
rho := CoveringRadius(C);
assert extDist eq expectedExternalDistance;
assert rho eq expectedCoveringRadius;
rhoLowerBound, rhoUpperBound := CoveringRadiusBounds(C);
assert (rhoLowerBound le rho) and (rho le rhoUpperBound);

isCR, InterArray, CosetDist, CR := IsCompletelyRegular(C);
assert isCR eq expectedIsCompletelyRegular;
