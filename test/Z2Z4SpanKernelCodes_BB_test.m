/*************************************************************/
/*                                                           */
/* Package or project name: Z2Z4AdditiveCodes package        */
/* Test file name: Z2Z4SpanKernelCodes_BB_test.m             */
/*                                                           */
/* Comments: Black-box tests for the intrinsic functions     */
/*            SpanZ2Code,  KernelZ2Code,                     */
/*            DimensionOfSpanZ2,  DimensionOfKernelZ2,       */
/*            KernelCosetRepresentatives, and                */
/*            CosetRepresentatives                           */  
/*            included in the Z2Z4AdditiveCodes.m file       */
/*                                                           */
/* Authors: M. Villanueva                                    */
/*                                                           */
/* Revision version and last date: v1.0, 2016/08/17          */
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

/************************************************************/
/*                                                          */
/* Function name: SpanZ2Code                                */
/* Parameters: C                                            */
/* Function description: Given a Z2Z4-additive code C of    */
/*   type (alpha, beta; gamma, delta; kappa), return        */
/*   S_C=Phi^(-1)(Sbin) as a Z2Z4-additive code, and        */
/*   Sbin=<Cbin>, that is the linear span of Cbin, as a     */
/*   binary linear code of length alpha + 2*beta, where     */
/*   Cbin=Phi(C) and Phi is the Gray map.                   */
/* Input parameters description:                            */
/*   - C: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - The Phi^(-1) of the linear span of Cbin as a         */
/*     Z2Z4-additive code                                   */
/*   - The linear span of Cbin as a binary linear code      */
/*                                                          */
/* Signature: (<Z2Z4Code> C) -> Z2Z4Code, CodeLinFld        */
/*                                                          */
/************************************************************/
/************************************************************/
/*                                                          */
/* Function name: DimensionOfSpanZ2                         */
/* Parameters: C                                            */
/* Function description: Given a Z2Z4-additive code C,      */
/*   return the dimension of the linear span of Cbin, that  */
/*   is, the dimension of <Cbin>, where Cbin=Phi(C) and Phi */
/*   is the Gray map.                                       */
/* Input parameters description:                            */
/*   - C: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - The dimension of the linear span of Cbin             */
/*                                                          */
/* Signature: (<Z2Z4Code> C) -> RngIntElt                   */
/*                                                          */
/************************************************************/
/************************************************************/
/*                                                          */
/* Function name: KernelZ2Code                              */
/* Parameters: C                                            */
/* Function description: Given a Z2Z4-additive code C of    */
/*   type (alpha, beta; gamma, delta; kappa), return its    */
/*   kernel K_C as a Z2Z4-additive subcode of C, and        */
/*   Kbin=Phi(K_C) as a binary linear subcode of Cbin of    */
/*   length alpha + 2 beta, where Cbin=Phi(C) and Phi is    */
/*   the Gray map.                                          */
/*   The kernel K_C contains the codewords v such that      */
/*   2v*u in C for all u in C, where * denotes the          */
/*   component-wise product. Equivalently, the kernel       */
/*   Kbin=Phi(K_C) contains the codewords c in Cbin such    */
/*   that c+Cbin=Cbin, where Cbin=Phi(C) and Phi is the     */
/*   Gray map.                                              */
/* Input parameters description:                            */
/*   - C: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - The kernel K_C as a Z2Z4-additive subcode of C       */
/*   - The kernel Phi(K_C) as a binary linear subcode of Cbin*/
/*                                                          */
/* Signature: (<Z2Z4Code> C) -> Z2Z4Code, CodeLinFld        */
/*                                                          */
/************************************************************/
/************************************************************/
/*                                                          */
/* Function name: KernelCosetRepresentatives                */
/* Parameters: C                                            */
/* Function description: Given a Z2Z4-additive code C of    */
/*   type (alpha, beta; gamma, delta; kappa), return the    */
/*   coset representatives [c_1,..., c_t] as a sequence of  */
/*   codewords of C, such that C = KC U (U_(i=1)^t (KC + c_i)),*/
/*   where KC is the kernel of C as a Z2Z4-additive subcode.*/
/*   It also returns the coset representatives of the       */
/*   corresponding binary code Cbin = phi(C) as a sequence  */
/*   of binary codewords [phi(c_1),..., phi(c_t)], such that*/
/*   Cbin = Kbin U (U_(i=1)^t (Kbin + phi(c_i))) where      */
/*   Kbin = phi(KC) and phi is the Gray map.                */
/* Input parameters description:                            */
/*   - C: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - The coset representatives [c_1,..., c_t]             */
/*   - The coset representatives of the binary code Cbin    */
/*                                                          */
/* Signature: (<Z2Z4Code> C) -> SeqEnum, SeqEnum            */
/*                                                          */
/************************************************************/
/************************************************************/
/*                                                          */
/* Function name: DimensionOfKernelZ2                       */
/* Parameters: C                                            */
/* Function description: Given a Z2Z4-additive code C,      */
/*   return the dimension of the Gray map image of its      */
/*   Z2Z4-additive kernel K_C, that is the dimension of     */
/*   Kbin = Phi(K_C), where Phi is the Gray map. Note that  */
/*   Kbin is always a binary linear code.                   */
/* Input parameters description:                            */
/*   - C: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - The dimension of the Gray map image of the kernel    */
/*                                                          */
/* Signature: (<Z2Z4Code> C) -> RngIntElt                   */
/*                                                          */
/************************************************************/
print "SPAN AND KERNEL CODES";
print "COSET REPRESENTATIVES";

print "test 1: Trivial Z2Z4-additive zero code
               alpha = 2, beta = 4, gamma = 0, delta = 0, kappa = 0, length = 6, #C = 1";
C := Z2Z4AdditiveZeroCode(2, 4);
V := VectorSpace(GF(2), 10);

expectedOutputSpan := C;
_, expectedOutputSpanZ2 := HasLinearGrayMapImage(C);
expectedOutputKernel := C;
expectedOutputKernelZ2 := expectedOutputSpanZ2;
expectedOutputSpanDimension := Dimension(expectedOutputSpanZ2);
expectedOutputKernelDimension := Dimension(expectedOutputKernelZ2);
expectedOutputKernelRepresentatives := [ ];
expectedOutputKernelRepresentativesZ2 := [ ]; 

outputSpan, outputSpanZ2 :=  SpanZ2Code(C);
outputKernel, outputKernelZ2 :=  KernelZ2Code(C);
outputv1KernelCosets, outputv1KernelCosetsZ2 :=  KernelCosetRepresentatives(C);
outputv2KernelCosets, outputv2KernelCosetsZ2 := 
                                 CosetRepresentatives(C, expectedOutputKernel); 
outputCosets, outputCosetsZ2 :=  CosetRepresentatives(C, C);

assert outputSpan eq expectedOutputSpan;
assert outputSpanZ2 eq expectedOutputSpanZ2;
assert outputKernel eq expectedOutputKernel;
assert outputKernelZ2 eq expectedOutputKernelZ2;
assert  DimensionOfSpanZ2(C) eq expectedOutputSpanDimension;
assert  DimensionOfKernelZ2(C) eq expectedOutputKernelDimension;
assert outputv1KernelCosets eq expectedOutputKernelRepresentatives;
assert outputv1KernelCosetsZ2 eq expectedOutputKernelRepresentativesZ2;
assert Append(outputv1KernelCosets, C`Code!0) eq Setseq(outputv2KernelCosets); 
assert Append(outputv1KernelCosetsZ2, V!0) eq Setseq(outputv2KernelCosetsZ2); 
assert outputCosets eq {@ C`Code!0 @};
assert outputCosetsZ2 eq {@ V!0 @};

/************************************************************/
print "test 2: Trivial Z2Z4-additive universe code
               alpha = 2, beta = 4, gamma = 2, delta = 4, kappa = 2, length = 6, #C = 1024";
C := Z2Z4AdditiveUniverseCode(2, 4);
V := VectorSpace(GF(2), 10);
R := RSpace(Z4, 6);
S1 := Z2Z4AdditiveZeroCode(2, 4);
S2 := Z2Z4AdditiveCode(Matrix(Z4, [[0,0,1,0,0,0], 
                                   [0,0,0,1,0,1],
                                   [0,0,0,0,1,0],
                                   [0,0,0,0,0,1]]), 2); 

expectedOutputSpan := C;
_, expectedOutputSpanZ2 := HasLinearGrayMapImage(C);
expectedOutputKernel := C;
expectedOutputKernelZ2 := expectedOutputSpanZ2;
expectedOutputSpanDimension := Dimension(expectedOutputSpanZ2);
expectedOutputKernelDimension := Dimension(expectedOutputKernelZ2);
expectedOutputKernelRepresentatives := [];
expectedOutputKernelRepresentativesZ2 := [];
expectedOutputS1Cosets := SetToIndexedSet(Set(C`Code));
expectedOutputS1CosetsZ2 := SetToIndexedSet({ V!v : v in  GrayMapImage(C) });
expectedOutputS2Cosets := {@ R!0, R![0,2,0,0,0,0], 
                                  R![2,0,0,0,0,0], 
                                  R![2,2,0,0,0,0] @};
expectedOutputS2CosetsZ2 := {@ V!0, V![0,1,0,0,0,0,0,0,0,0], 
                                  V![1,0,0,0,0,0,0,0,0,0], 
                                  V![1,1,0,0,0,0,0,0,0,0] @};

outputSpan, outputSpanZ2 :=  SpanZ2Code(C);
outputKernel, outputKernelZ2 :=  KernelZ2Code(C);
outputv1KernelCosets, outputv1KernelCosetsZ2 :=  KernelCosetRepresentatives(C);
outputv2KernelCosets, outputv2KernelCosetsZ2 := 
                            CosetRepresentatives(C, expectedOutputKernel); 
outputCosets, outputCosetsZ2 :=  CosetRepresentatives(C, C);
outputS1Cosets, outputS1CosetsZ2 :=  CosetRepresentatives(C, S1);
outputS2Cosets, outputS2CosetsZ2 :=  CosetRepresentatives(C, S2);

assert outputSpan eq expectedOutputSpan;
assert outputSpanZ2 eq expectedOutputSpanZ2;
assert outputKernel eq expectedOutputKernel;
assert outputKernelZ2 eq expectedOutputKernelZ2;
assert  DimensionOfSpanZ2(C) eq expectedOutputSpanDimension;
assert  DimensionOfKernelZ2(C) eq expectedOutputKernelDimension;
assert outputv1KernelCosets eq expectedOutputKernelRepresentatives;
assert outputv1KernelCosetsZ2 eq expectedOutputKernelRepresentativesZ2;
assert Append(outputv1KernelCosets, C`Code!0) eq Setseq(outputv2KernelCosets); 
assert Append(outputv1KernelCosetsZ2, V!0) eq Setseq(outputv2KernelCosetsZ2); 
assert outputCosets eq {@ C`Code!0 @};
assert outputCosetsZ2 eq {@ V!0 @};
assert outputS1Cosets eq expectedOutputS1Cosets;
assert outputS1CosetsZ2 eq expectedOutputS1CosetsZ2;
assert outputS2Cosets eq expectedOutputS2Cosets;
assert outputS2CosetsZ2 eq expectedOutputS2CosetsZ2;

/************************************************************/
print "test 3: Repetition Z2Z4-additive code
               alpha = 4, beta = 8, gamma = 1, delta = 0, kappa = 1, length = 12, #C = 2";
C := Z2Z4AdditiveRepetitionCode(4, 8);
V := VectorSpace(GF(2), 20);
S1 := Z2Z4AdditiveZeroCode(4, 8);

expectedOutputSpan := C;
_, expectedOutputSpanZ2 := HasLinearGrayMapImage(C);
expectedOutputKernel := C;
expectedOutputKernelZ2 := expectedOutputSpanZ2;
expectedOutputSpanDimension := Dimension(expectedOutputSpanZ2);
expectedOutputKernelDimension := Dimension(expectedOutputKernelZ2);
expectedOutputKernelRepresentatives := [ ];
expectedOutputKernelRepresentativesZ2 := [ ];
expectedOutputS1Cosets := SetToIndexedSet(Set(C`Code));
expectedOutputS1CosetsZ2 := SetToIndexedSet({ V!v : v in  GrayMapImage(C) }); 

outputSpan, outputSpanZ2 :=  SpanZ2Code(C);
outputKernel, outputKernelZ2 :=  KernelZ2Code(C);
outputv1KernelCosets, outputv1KernelCosetsZ2 :=  KernelCosetRepresentatives(C);
outputv2KernelCosets, outputv2KernelCosetsZ2 := 
                            CosetRepresentatives(C, expectedOutputKernel); 
outputCosets, outputCosetsZ2 :=  CosetRepresentatives(C, C);
outputS1Cosets, outputS1CosetsZ2 :=  CosetRepresentatives(C, S1);

assert outputSpan eq expectedOutputSpan;
assert outputSpanZ2 eq expectedOutputSpanZ2;
assert outputKernel eq expectedOutputKernel;
assert outputKernelZ2 eq expectedOutputKernelZ2;
assert  DimensionOfSpanZ2(C) eq expectedOutputSpanDimension;
assert  DimensionOfKernelZ2(C) eq expectedOutputKernelDimension;
assert outputv1KernelCosets eq expectedOutputKernelRepresentatives;
assert outputv1KernelCosetsZ2 eq expectedOutputKernelRepresentativesZ2;
assert Append(outputv1KernelCosets, C`Code!0) eq Setseq(outputv2KernelCosets); 
assert Append(outputv1KernelCosetsZ2, V!0) eq Setseq(outputv2KernelCosetsZ2); 
assert outputCosets eq {@ C`Code!0 @};
assert outputCosetsZ2 eq {@ V!0 @};
assert outputS1Cosets eq expectedOutputS1Cosets;
assert outputS1CosetsZ2 eq expectedOutputS1CosetsZ2;

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
V := VectorSpace(GF(2), 19*2);
R := RSpace(Z4, 19);
S1 := Z2Z4AdditiveZeroCode(0, 19);

expectedOutputSpan := Z2Z4AdditiveCode(VerticalJoin(GeneratorMatrix(C`Code), 
                       Matrix(Z4,[[0,0,0,0,0,0,0,0,0,2,0,0,2,0,0,2,2,0,2]])), 0);
_, expectedOutputSpanZ2 := HasLinearGrayMapImage(expectedOutputSpan);
expectedOutputKernel :=  OrderTwoSubcode(C);
_, expectedOutputKernelZ2 := HasLinearGrayMapImage(expectedOutputKernel);
expectedOutputSpanDimension := Dimension(expectedOutputSpanZ2);
expectedOutputKernelDimension := Dimension(expectedOutputKernelZ2);
expectedOutputKernelRepresentatives := [ R![0,0,0,3,0,0,0,3,2,1,0,3,1,0,0,0,0,0,1],
                                         R![1,1,1,0,0,1,0,0,0,3,2,1,3,2,2,0,0,1,1],
                                         R![1,1,1,3,0,1,0,3,2,0,2,0,0,2,2,0,0,1,2]];
expectedOutputKernelRepresentativesZ2 := [ V! GrayMap(C)(v) : v in expectedOutputKernelRepresentatives];
expectedOutputS1Cosets := SetToIndexedSet(Set(C`Code));
expectedOutputS1CosetsZ2 := SetToIndexedSet({ V!v : v in  GrayMapImage(C) });

outputSpan, outputSpanZ2 :=  SpanZ2Code(C);
outputKernel, outputKernelZ2 :=  KernelZ2Code(C);
outputv1KernelCosets, outputv1KernelCosetsZ2 :=  KernelCosetRepresentatives(C);
outputv2KernelCosets, outputv2KernelCosetsZ2 := 
                            CosetRepresentatives(C, expectedOutputKernel); 
outputCosets, outputCosetsZ2 :=  CosetRepresentatives(C, C);
outputS1Cosets, outputS1CosetsZ2 :=  CosetRepresentatives(C, S1);

assert  outputSpan eq expectedOutputSpan;
assert outputSpanZ2 eq expectedOutputSpanZ2;
assert  outputKernel eq expectedOutputKernel;
assert outputKernelZ2 eq expectedOutputKernelZ2;
assert  DimensionOfSpanZ2(C) eq expectedOutputSpanDimension;
assert  DimensionOfKernelZ2(C) eq expectedOutputKernelDimension;
assert outputv1KernelCosets eq expectedOutputKernelRepresentatives;
assert outputv1KernelCosetsZ2 eq expectedOutputKernelRepresentativesZ2;
assert Set(Append(outputv1KernelCosets, C`Code!0)) eq outputv2KernelCosets; 
assert Set(Append(outputv1KernelCosetsZ2, V!0)) eq outputv2KernelCosetsZ2; 
assert outputCosets eq {@ C`Code!0 @};
assert outputCosetsZ2 eq {@ V!0 @};
assert outputS1Cosets eq expectedOutputS1Cosets;
assert outputS1CosetsZ2 eq expectedOutputS1CosetsZ2;

/************************************************************/
print "test 5: A Z2Z4-additive code with beta=0
               alpha = 10, beta = 0, gamma = 4, delta = 0, kappa = 4, length = 10, #C = 16";
R := RSpace(Z4, 10);
L := [R![0,0,1,0,1,0,1,0,1,0],
       R![1,0,0,0,0,0,1,1,1,1],
       R![0,0,1,0,0,0,0,1,0,0],
       R![0,1,0,1,0,1,0,1,1,1]];
C := Z2Z4AdditiveCode(L, 10);
V := VectorSpace(GF(2), 10);
R := RSpace(Z4, 10);
S1 := Z2Z4AdditiveZeroCode(10, 0);

expectedOutputSpan := C;
_, expectedOutputSpanZ2 := HasLinearGrayMapImage(expectedOutputSpan);
expectedOutputKernel := C;
_, expectedOutputKernelZ2 := HasLinearGrayMapImage(expectedOutputKernel);
expectedOutputSpanDimension := Dimension(expectedOutputSpanZ2);
expectedOutputKernelDimension := Dimension(expectedOutputKernelZ2);
expectedOutputKernelRepresentatives := [ ];
expectedOutputKernelRepresentativesZ2 := [ ];
expectedOutputS1Cosets := SetToIndexedSet(Set(C`Code));
expectedOutputS1CosetsZ2 := SetToIndexedSet({ V!v : v in  GrayMapImage(C) });

outputSpan, outputSpanZ2 :=  SpanZ2Code(C);
outputKernel, outputKernelZ2 :=  KernelZ2Code(C);
outputv1KernelCosets, outputv1KernelCosetsZ2 :=  KernelCosetRepresentatives(C);
outputv2KernelCosets, outputv2KernelCosetsZ2 := 
                            CosetRepresentatives(C, expectedOutputKernel); 
outputCosets, outputCosetsZ2 :=  CosetRepresentatives(C, C);
outputS1Cosets, outputS1CosetsZ2 :=  CosetRepresentatives(C, S1);

assert  outputSpan eq expectedOutputSpan;
assert outputSpanZ2 eq expectedOutputSpanZ2;
assert  outputKernel eq expectedOutputKernel;
assert outputKernelZ2 eq expectedOutputKernelZ2;
assert  DimensionOfSpanZ2(C) eq expectedOutputSpanDimension;
assert  DimensionOfKernelZ2(C) eq expectedOutputKernelDimension;
assert outputv1KernelCosets eq expectedOutputKernelRepresentatives;
assert outputv1KernelCosetsZ2 eq expectedOutputKernelRepresentativesZ2;
assert Set(Append(outputv1KernelCosets, C`Code!0)) eq outputv2KernelCosets; 
assert Set(Append(outputv1KernelCosetsZ2, V!0)) eq outputv2KernelCosetsZ2; 
assert outputCosets eq {@ C`Code!0 @};
assert outputCosetsZ2 eq {@ V!0 @};
assert outputS1Cosets eq expectedOutputS1Cosets;
assert outputS1CosetsZ2 eq expectedOutputS1CosetsZ2;

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
V := VectorSpace(GF(2), 48);
R := RSpace(Z4, 29);
S1 := Z2Z4AdditiveZeroCode(10, 19);

expectedOutputSpan := Z2Z4AdditiveCode(SpanZ2CodeZ4(C`Code), 10); 
_, expectedOutputSpanZ2 := HasLinearGrayMapImage(expectedOutputSpan);
expectedOutputKernel :=  OrderTwoSubcode(C);
_, expectedOutputKernelZ2 := HasLinearGrayMapImage(expectedOutputKernel);
expectedOutputSpanDimension := Dimension(expectedOutputSpanZ2);
expectedOutputKernelDimension := Dimension(expectedOutputKernelZ2);
expectedOutputKernelRepresentatives := KernelCosetRepresentatives(C`Code);
expectedOutputKernelRepresentativesZ2 := [ V! GrayMap(C)(v) : v in expectedOutputKernelRepresentatives];
expectedOutputS1Cosets := SetToIndexedSet(Set(C`Code));
expectedOutputS1CosetsZ2 := SetToIndexedSet({ V!v : v in  GrayMapImage(C) });

outputSpan, outputSpanZ2 :=  SpanZ2Code(C);
outputKernel, outputKernelZ2 :=  KernelZ2Code(C);
outputv1KernelCosets, outputv1KernelCosetsZ2 :=  KernelCosetRepresentatives(C);
outputv2KernelCosets, outputv2KernelCosetsZ2 := 
                            CosetRepresentatives(C, expectedOutputKernel); 
outputCosets, outputCosetsZ2 :=  CosetRepresentatives(C, C);
outputS1Cosets, outputS1CosetsZ2 :=  CosetRepresentatives(C, S1);

assert  outputSpan eq expectedOutputSpan;
assert outputSpanZ2 eq expectedOutputSpanZ2;
assert  outputKernel eq expectedOutputKernel;
assert outputKernelZ2 eq expectedOutputKernelZ2;
assert  DimensionOfSpanZ2(C) eq expectedOutputSpanDimension;
assert  DimensionOfKernelZ2(C) eq expectedOutputKernelDimension;
assert outputv1KernelCosets eq expectedOutputKernelRepresentatives;
assert outputv1KernelCosetsZ2 eq expectedOutputKernelRepresentativesZ2;
assert Set(Append(outputv1KernelCosets, C`Code!0)) eq outputv2KernelCosets; 
assert Set(Append(outputv1KernelCosetsZ2, V!0)) eq outputv2KernelCosetsZ2; 
assert outputCosets eq {@ C`Code!0 @};
assert outputCosetsZ2 eq {@ V!0 @};
assert outputS1Cosets eq expectedOutputS1Cosets;
assert outputS1CosetsZ2 eq expectedOutputS1CosetsZ2;

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
V := VectorSpace(GF(2), 48);
R := RSpace(Z4, 29);
S1 := Z2Z4AdditiveZeroCode(10, 19);

expectedOutputSpan := C;
_, expectedOutputSpanZ2 := HasLinearGrayMapImage(expectedOutputSpan);
expectedOutputKernel := C;
_, expectedOutputKernelZ2 := HasLinearGrayMapImage(expectedOutputKernel);
expectedOutputSpanDimension := Dimension(expectedOutputSpanZ2);
expectedOutputKernelDimension := Dimension(expectedOutputKernelZ2);
expectedOutputKernelRepresentatives := [ ];
expectedOutputKernelRepresentativesZ2 := [ ];
expectedOutputS1Cosets := SetToIndexedSet(Set(C`Code));
expectedOutputS1CosetsZ2 := SetToIndexedSet({ V!v : v in  GrayMapImage(C) });

outputSpan, outputSpanZ2 :=  SpanZ2Code(C);
outputKernel, outputKernelZ2 :=  KernelZ2Code(C);
outputv1KernelCosets, outputv1KernelCosetsZ2 :=  KernelCosetRepresentatives(C);
outputv2KernelCosets, outputv2KernelCosetsZ2 := 
                            CosetRepresentatives(C, expectedOutputKernel); 
outputCosets, outputCosetsZ2 :=  CosetRepresentatives(C, C);
outputS1Cosets, outputS1CosetsZ2 :=  CosetRepresentatives(C, S1);

assert  outputSpan eq expectedOutputSpan;
assert outputSpanZ2 eq expectedOutputSpanZ2;
assert  outputKernel eq expectedOutputKernel;
assert outputKernelZ2 eq expectedOutputKernelZ2;
assert  DimensionOfSpanZ2(C) eq expectedOutputSpanDimension;
assert  DimensionOfKernelZ2(C) eq expectedOutputKernelDimension;
assert outputv1KernelCosets eq expectedOutputKernelRepresentatives;
assert outputv1KernelCosetsZ2 eq expectedOutputKernelRepresentativesZ2;
assert Set(Append(outputv1KernelCosets, C`Code!0)) eq outputv2KernelCosets; 
assert Set(Append(outputv1KernelCosetsZ2, V!0)) eq outputv2KernelCosetsZ2; 
assert outputCosets eq {@ C`Code!0 @};
assert outputCosetsZ2 eq {@ V!0 @};
assert outputS1Cosets eq expectedOutputS1Cosets;
assert outputS1CosetsZ2 eq expectedOutputS1CosetsZ2;

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
V := VectorSpace(GF(2), 48);
R := RSpace(Z4, 29);
S1 := Z2Z4AdditiveZeroCode(10, 19);

expectedOutputSpan := Z2Z4AdditiveCode(SpanZ2CodeZ4(C`Code), 10);
_, expectedOutputSpanZ2 := HasLinearGrayMapImage(expectedOutputSpan);
expectedOutputKernel :=  OrderTwoSubcode(C);
_, expectedOutputKernelZ2 := HasLinearGrayMapImage(expectedOutputKernel);
expectedOutputSpanDimension := Dimension(expectedOutputSpanZ2);
expectedOutputKernelDimension := Dimension(expectedOutputKernelZ2);
expectedOutputKernelRepresentatives := KernelCosetRepresentatives(C`Code);
expectedOutputKernelRepresentativesZ2 := [ V! GrayMap(C)(v) : v in expectedOutputKernelRepresentatives];
expectedOutputS1Cosets := SetToIndexedSet(Set(C`Code));
expectedOutputS1CosetsZ2 := SetToIndexedSet({ V!v : v in  GrayMapImage(C) });

outputSpan, outputSpanZ2 :=  SpanZ2Code(C);
outputKernel, outputKernelZ2 :=  KernelZ2Code(C);
outputv1KernelCosets, outputv1KernelCosetsZ2 :=  KernelCosetRepresentatives(C);
outputv2KernelCosets, outputv2KernelCosetsZ2 := 
                            CosetRepresentatives(C, expectedOutputKernel); 
outputCosets, outputCosetsZ2 :=  CosetRepresentatives(C, C);
outputS1Cosets, outputS1CosetsZ2 :=  CosetRepresentatives(C, S1);

assert  outputSpan eq expectedOutputSpan;
assert outputSpanZ2 eq expectedOutputSpanZ2;
assert  outputKernel eq expectedOutputKernel;
assert outputKernelZ2 eq expectedOutputKernelZ2;
assert  DimensionOfSpanZ2(C) eq expectedOutputSpanDimension;
assert  DimensionOfKernelZ2(C) eq expectedOutputKernelDimension;
assert outputv1KernelCosets eq expectedOutputKernelRepresentatives;
assert outputv1KernelCosetsZ2 eq expectedOutputKernelRepresentativesZ2;
assert Set(Append(outputv1KernelCosets, C`Code!0)) eq outputv2KernelCosets; 
assert Set(Append(outputv1KernelCosetsZ2, V!0)) eq outputv2KernelCosetsZ2; 
assert outputCosets eq {@ C`Code!0 @};
assert outputCosetsZ2 eq {@ V!0 @};
assert outputS1Cosets eq expectedOutputS1Cosets;
assert outputS1CosetsZ2 eq expectedOutputS1CosetsZ2;

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
V := VectorSpace(GF(2), 31);
S1 := Z2Z4AdditiveZeroCode(7, 12);

expectedOutputSpan := Z2Z4AdditiveCode(SpanZ2CodeZ4(C`Code), 7);
_, expectedOutputSpanZ2 := HasLinearGrayMapImage(expectedOutputSpan);
expectedOutputKernel :=  OrderTwoSubcode(C);
_, expectedOutputKernelZ2 := HasLinearGrayMapImage(expectedOutputKernel);
expectedOutputSpanDimension := Dimension(expectedOutputSpanZ2);
expectedOutputKernelDimension := Dimension(expectedOutputKernelZ2);
expectedOutputKernelRepresentatives := KernelCosetRepresentatives(C`Code);
expectedOutputKernelRepresentativesZ2 := [ V! GrayMap(C)(v) : v in expectedOutputKernelRepresentatives];
expectedOutputS1Cosets := SetToIndexedSet(Set(C`Code));
expectedOutputS1CosetsZ2 := SetToIndexedSet({ V!v : v in  GrayMapImage(C) });

outputSpan, outputSpanZ2 :=  SpanZ2Code(C);
outputKernel, outputKernelZ2 :=  KernelZ2Code(C);
outputv1KernelCosets, outputv1KernelCosetsZ2 :=  KernelCosetRepresentatives(C);
outputv2KernelCosets, outputv2KernelCosetsZ2 := 
                            CosetRepresentatives(C, expectedOutputKernel); 
outputCosets, outputCosetsZ2 :=  CosetRepresentatives(C, C);
outputS1Cosets, outputS1CosetsZ2 :=  CosetRepresentatives(C, S1);

assert  outputSpan eq expectedOutputSpan;
assert outputSpanZ2 eq expectedOutputSpanZ2;
assert  outputKernel eq expectedOutputKernel;
assert outputKernelZ2 eq expectedOutputKernelZ2;
assert  DimensionOfSpanZ2(C) eq expectedOutputSpanDimension;
assert  DimensionOfKernelZ2(C) eq expectedOutputKernelDimension;
assert outputv1KernelCosets eq expectedOutputKernelRepresentatives;
assert outputv1KernelCosetsZ2 eq expectedOutputKernelRepresentativesZ2;
assert Set(Append(outputv1KernelCosets, C`Code!0)) eq outputv2KernelCosets; 
assert Set(Append(outputv1KernelCosetsZ2, V!0)) eq outputv2KernelCosetsZ2; 
assert outputCosets eq {@ C`Code!0 @};
assert outputCosetsZ2 eq {@ V!0 @};
assert outputS1Cosets eq expectedOutputS1Cosets;
assert outputS1CosetsZ2 eq expectedOutputS1CosetsZ2;

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
V := VectorSpace(GF(2), 33);
S1 := Z2Z4AdditiveZeroCode(9, 12);

expectedOutputSpan := Z2Z4AdditiveCode(SpanZ2CodeZ4(C`Code), 9);
_, expectedOutputSpanZ2 := HasLinearGrayMapImage(expectedOutputSpan);
expectedOutputKernel :=  Z2Z4AdditiveCode(KernelZ2CodeZ4(C`Code), 9);
_, expectedOutputKernelZ2 := HasLinearGrayMapImage(expectedOutputKernel);
expectedOutputSpanDimension := Dimension(expectedOutputSpanZ2);
expectedOutputKernelDimension := Dimension(expectedOutputKernelZ2);
expectedOutputKernelRepresentatives := KernelCosetRepresentatives(C`Code);
expectedOutputKernelRepresentativesZ2 := [ V! GrayMap(C)(v) : v in expectedOutputKernelRepresentatives];
expectedOutputS1Cosets := SetToIndexedSet(Set(C`Code));
expectedOutputS1CosetsZ2 := SetToIndexedSet({ V!v : v in  GrayMapImage(C) });

outputSpan, outputSpanZ2 :=  SpanZ2Code(C);
outputKernel, outputKernelZ2 :=  KernelZ2Code(C);
outputv1KernelCosets, outputv1KernelCosetsZ2 :=  KernelCosetRepresentatives(C);
outputv2KernelCosets, outputv2KernelCosetsZ2 := 
                            CosetRepresentatives(C, expectedOutputKernel); 
outputCosets, outputCosetsZ2 :=  CosetRepresentatives(C, C);
outputS1Cosets, outputS1CosetsZ2 :=  CosetRepresentatives(C, S1);

assert  outputSpan eq expectedOutputSpan;
assert outputSpanZ2 eq expectedOutputSpanZ2;
assert  outputKernel eq expectedOutputKernel;
assert outputKernelZ2 eq expectedOutputKernelZ2;
assert  DimensionOfSpanZ2(C) eq expectedOutputSpanDimension;
assert  DimensionOfKernelZ2(C) eq expectedOutputKernelDimension;
assert outputv1KernelCosets eq expectedOutputKernelRepresentatives;
assert outputv1KernelCosetsZ2 eq expectedOutputKernelRepresentativesZ2;
assert Set(Append(outputv1KernelCosets, C`Code!0)) eq outputv2KernelCosets; 
assert Set(Append(outputv1KernelCosetsZ2, V!0)) eq outputv2KernelCosetsZ2; 
assert outputCosets eq {@ C`Code!0 @};
assert outputCosetsZ2 eq {@ V!0 @};
assert outputS1Cosets eq expectedOutputS1Cosets;
assert outputS1CosetsZ2 eq expectedOutputS1CosetsZ2;

