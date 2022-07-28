/*************************************************************/
/*                                                           */
/* Package or project name: Z2Z4AdditiveCodes package        */
/* Test file name: Z2Z4MinimumWeight_BB_test.m      	     */
/*                                                           */
/* Comments: Black-box tests for the intrinsic functions     */
/*        	 Z2Z4MinimumLeeWeight, Z2Z4MinimumLeeDistance,   */
/*           MinimumWord, MinimumWords,                      */
/*           LeeWeightDistribution,                          */
/*           LeeWeightDistributionDual and                   */
/*           ExternalDistance included in the                */
/*           Z2Z4MinimumWeight.m file                        */
/*                                                           */
/* Authors: Cristina Diéguez, Cristina Fernández-Córdoba,    */
/*          Jaume Pujol and Mercè Villanueva                 */
/*                                                           */
/* Revision version and last date: v1.0   2016/06/24         */                  
/*                                 v1.1   2016/08/15         */
/*                                 v1.2   2018/02/12         */
/*         user defined type       v1.3   2019/01/30         */
/*                                                           */
/*************************************************************/

//needs Z2Z4AdditiveCode package

SetAssertions(true);
Alarm(30*60); 

/*************************************************************
	GLOBAL VARIABLES
*************************************************************/	
Z4 := Integers(4);

SetVerbose("IgnoreWeightAttributes", 1);

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4MinimumLeeWeight                          */
/* Parameters: C                                                */
/* Function description: Given a Z2Z4-additive code C, return   */
/*   the minimum Lee weight of the codewords belonging to the   */
/*   code C, which is also the minimum Lee distance between any */
/*   two codewords.                                             */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/* Output parameters description:                               */
/*   - The minimum Lee weight of the code                       */
/*                                                              */
/* Function first developed by C. Diéguez-Martí and             */
/*   updated by A. Torres                                       */
/*                                                              */
/* Signature: (Z2Z4Code C) -> RngIntElt                         */
/*                                                              */
/****************************************************************/
/****************************************************************/
/*                                                              */
/* Function name: Z2Z4MinimumLeeDistance                        */
/* Parameters:  C                                               */
/* Function description: Given a Z2Z4-additive code C, return   */
/*   the minimum Lee weight of the codewords belonging to the   */
/*   code C, which is also the minimum Lee distance between any */
/*   two codewords.                                             */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/* Output parameters description:                               */
/*   - The minimum Lee distance of the code                     */
/*                                                              */
/* Function first developed by C. Diéguez-Martí and             */
/*   updated by A. Torres                                       */
/*                                                              */
/* Signature: (Z2Z4Code C) -> RngIntElt                         */
/*                                                              */
/****************************************************************/
/****************************************************************/
/*                                                              */
/* Function name: MinimumWord                                   */
/* Parameters: C                                                */
/* Function description: Given a Z2Z4-additive code C, return   */
/*   one codeword of the code C having minimum Lee weight.      */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/* Output parameters description:                               */
/*   - A codeword with minimum Lee weight                       */
/*                                                              */
/* Function first developed by (?) and updated by A. Torres     */
/*                                                              */
/* Signature: (Z2Z4Code C) -> ModTupFldElt                      */
/*                                                              */
/****************************************************************/
/****************************************************************/
/*                                                              */
/* Function name: MinimumWords                                  */
/* Parameters: C                                                */
/* Function description: Given a Z2Z4-additive code C, return   */
/*   the set of all codewords of C having minimum Lee weight.   */
/*   If NumWords is set to a non-negative integer, then the     */
/*   algorithm will terminate after at least that total of      */
/*   codewords have been found.                                 */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/* Output parameters description:                               */
/*   - A set with all codewords of minimum Lee weight           */
/*                                                              */
/* Function first developed by C. Diéguez-Martí and             */
/*   updated by A. Torres                                       */
/*                                                              */
/* Signature: (Z2Z4Code C) -> Set                               */
/*                                                              */
/****************************************************************/
/****************************************************************/
/*                                                              */
/* Function name: LeeWeightDistribution                         */
/* Parameters: C                                                */
/* Function description: Determine the Lee weight distribution  */
/*   for the Z2Z4-additive code C. The distribution is returned */
/*   in the form of a sequence of tuples, where the i-th tuple  */
/*   contains the i-th weight, wi say, and the number of        */
/*   codewords having weight wi.                                */
/* Input parameters description:                                */
/*   - C: a Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - Sequence of tuples <Lee weight, number of codewords>     */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> SeqEnum                         */
/*                                                              */
/****************************************************************/
/****************************************************************/
/*                                                              */
/* Function name: DualLeeWeightDistribution                     */
/* Parameters: C                                                */
/* Function description: Determine the Lee weight distribution  */
/*   for the additive dual code of C. The distribution is       */
/*   returned in the form of a sequence of tuples, where the    */
/*   i-th tuple contains the i-th weight, wi say, and the       */
/*   number of codewords having weight wi.                      */      
/* Input parameters description:                                */
/*   - C: a Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - Enumerated sequence of tuples with the weight            */
/*     distribution                                             */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> SeqEnum                         */
/*                                                              */
/****************************************************************/
/****************************************************************/
/*                                                              */
/* Function name: ExternalDistance                              */
/* Parameters: C                                                */
/* Function description: Determine the external distance for    */
/*   the Z2Z4-additive code C. The external distance of a code  */
/*   is the number of different nonzero weights of the dual     */
/*   code of C.                                                 */      
/* Input parameters description:                                */
/*   - C: a Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - Integer with the external distance                       */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> RngIntElt                       */
/*                                                              */
/****************************************************************/
print "THE WEIGHT DISTRIBUTION";

print "test 1: Trivial Z2Z4-additive zero code 
               alpha = 2, beta = 4, gamma = 0, delta = 0, kappa = 0, length = 6, #C = 1"; 
alpha := 2;
beta:= 4;
gamma := 0;
delta := 0;
C := Z2Z4AdditiveCode([RSpace(Z4, 6)!0], 2);

expectedOutputWeigthDistribution := [ <0, 1> ];
expectedOutputWeigthDistributionDual := 
    MacWilliamsTransform(alpha+2*beta, gamma+2*delta, 2, expectedOutputWeigthDistribution);
expectedOutputExternalDistance := 10; 

assert LeeWeightDistribution(C) eq expectedOutputWeigthDistribution;
assert DualLeeWeightDistribution(C) eq expectedOutputWeigthDistributionDual;
assert ExternalDistance(C) eq expectedOutputExternalDistance;

assert LeeWeightDistribution(C : Method := "Distribution") eq expectedOutputWeigthDistribution;

assert LeeWeightDistribution(C : Method := "KernelCosets") eq expectedOutputWeigthDistribution;

/************************************************************/
print "test 2: Trivial Z2Z4-additive universe code 
               alpha = 2, beta = 6, gamma = 2, delta = 6, kappa = 2, length = 8, #C = 16384";
alpha := 2;
beta:= 6;
gamma := 2;
delta := 6;
C := Z2Z4AdditiveUniverseCode(alpha, beta);
delete C`MinimumLeeWeight;

expectedOutputMinimumLeeDistance :=  1;
expectedOutputMinimumWords := MinimumWords(C : Method := "Distribution");
expectedOutputWeigthDistribution := [ <0, 1>, <1, 14>, <2, 91>, <3, 364>, <4, 1001>, 
   <5, 2002>, <6, 3003>, <7, 3432>, <8, 3003>, <9, 2002>, <10, 1001>, <11, 364>, 
   <12, 91>, <13, 14>, <14, 1> ];
expectedOutputWeigthDistributionDual := 
    MacWilliamsTransform(alpha+2*beta, gamma+2*delta, 2, expectedOutputWeigthDistribution);
expectedOutputExternalDistance := 0; 

assert Z2Z4MinimumLeeWeight(C) eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C) eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C) in expectedOutputMinimumWords;
assert MinimumWords(C) eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C) eq expectedOutputWeigthDistribution;
assert DualLeeWeightDistribution(C) eq expectedOutputWeigthDistributionDual;
assert ExternalDistance(C) eq expectedOutputExternalDistance;

assert Z2Z4MinimumLeeWeight(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Distribution") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Distribution") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "Distribution") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "KernelCosets") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "KernelCosets") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "KernelCosets") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Brouwer") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Brouwer") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Zimmermann") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Zimmermann") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Quaternary") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Quaternary") eq expectedOutputMinimumWords;

/************************************************************/
print "test 3: Repetition Z2Z4-additive code
               alpha = 4, beta = 8, gamma = 1, delta = 0, kappa = 1, length = 12, #C = 2";
alpha := 4;
beta:= 8;
gamma := 1;
delta := 0;
C := Z2Z4AdditiveRepetitionCode(alpha, beta);
delete C`MinimumLeeWeight;

expectedOutputMinimumLeeDistance :=  20;
expectedOutputMinimumWords := MinimumWords(C : Method := "Distribution");
expectedOutputWeigthDistribution := [ <0, 1>, <20, 1> ];
expectedOutputWeigthDistributionDual := 
    MacWilliamsTransform(alpha+2*beta, gamma+2*delta, 2, expectedOutputWeigthDistribution);
//[ <0, 1>, <2, 190>, <4, 4845>, <6, 38760>, <8, 125970>, 
//<10, 184756>, <12, 125970>, <14, 38760>, <16, 4845>, <18, 190>, <20, 1> ];
expectedOutputExternalDistance := 10; 

assert Z2Z4MinimumLeeWeight(C) eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C) eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C) in expectedOutputMinimumWords;
assert MinimumWords(C) eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C) eq expectedOutputWeigthDistribution;
assert DualLeeWeightDistribution(C) eq expectedOutputWeigthDistributionDual;
assert ExternalDistance(C) eq expectedOutputExternalDistance;

assert Z2Z4MinimumLeeWeight(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Distribution") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Distribution") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "Distribution") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "KernelCosets") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "KernelCosets") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "KernelCosets") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Brouwer") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Brouwer") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Zimmermann") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Zimmermann") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Quaternary") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Quaternary") eq expectedOutputMinimumWords;

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
alpha := 10;
beta := 20;
gamma := 6;
delta := 2;
C := Z2Z4AdditiveCode(M, alpha);

expectedOutputMinimumLeeDistance :=  10;
expectedOutputMinimumWords :=  MinimumWords(C : Method := "Distribution");
expectedOutputWeigthDistribution := [ <0, 1>, <10, 1>, <12, 5>, <13, 3>, <14, 10>, 
<15, 4>, <16, 14>, <17, 15>, <18, 26>, <19, 31>, <20, 40>, <21, 55>, <22, 72>, 
<23, 118>, <24, 98>, <25, 126>, <26, 85>, <27, 72>, <28, 75>, <29, 53>, <30, 62>, 
<31, 30>, <32, 23>, <33, 3>, <35, 1>, <37, 1> ];
expectedOutputWeigthDistributionDual := 
    MacWilliamsTransform(alpha+2*beta, gamma+2*delta, 2, expectedOutputWeigthDistribution);
expectedOutputExternalDistance := 48;

assert Z2Z4MinimumLeeWeight(C) eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C) eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C) in expectedOutputMinimumWords;
assert MinimumWords(C) eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C) eq expectedOutputWeigthDistribution;
assert DualLeeWeightDistribution(C) eq expectedOutputWeigthDistributionDual;
assert ExternalDistance(C) eq expectedOutputExternalDistance;

assert Z2Z4MinimumLeeWeight(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Distribution") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Distribution") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "Distribution") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "KernelCosets") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "KernelCosets") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "KernelCosets") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Brouwer") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Brouwer") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Zimmermann") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Zimmermann") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Quaternary") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Quaternary") eq expectedOutputMinimumWords;

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
alpha := 10;
beta := 19;
gamma := 7;
delta := 2;
C := Z2Z4AdditiveCode(M, alpha);

expectedOutputMinimumLeeDistance :=  6;
expectedOutputMinimumWords :=  MinimumWords(C : Method := "Distribution");
expectedOutputWeigthDistribution := [ <0, 1>, <6, 1>, <10, 1>, <11, 1>, <12, 13>, 
<14, 21>, <15, 3>, <16, 31>, <17, 41>, <18, 53>, <19, 92>, <20, 109>, <21, 122>, 
<22, 197>, <23, 204>, <24, 224>, <25, 220>, <26, 189>, <27, 135>, <28, 117>, 
<29, 118>, <30, 45>, <31, 73>, <32, 16>, <33, 11>, <34, 5>, <35, 4>, <36, 1> ];
expectedOutputWeigthDistributionDual := 
    MacWilliamsTransform(alpha+2*beta, gamma+2*delta, 2, expectedOutputWeigthDistribution);
expectedOutputExternalDistance := 46;

assert Z2Z4MinimumLeeWeight(C) eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C) eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C) in expectedOutputMinimumWords;
assert MinimumWords(C) eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C) eq expectedOutputWeigthDistribution;
assert DualLeeWeightDistribution(C) eq expectedOutputWeigthDistributionDual;
assert ExternalDistance(C) eq expectedOutputExternalDistance;

assert Z2Z4MinimumLeeWeight(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Distribution") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Distribution") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "Distribution") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "KernelCosets") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "KernelCosets") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "KernelCosets") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Brouwer") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Brouwer") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Zimmermann") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Zimmermann") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Quaternary") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Quaternary") eq expectedOutputMinimumWords;

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
alpha := 10;
beta := 19;
gamma := 7;
delta := 2;
C := Z2Z4AdditiveCode(M, alpha);	

expectedOutputMinimumLeeDistance :=  8;
expectedOutputMinimumWords :=  MinimumWords(C : Method := "Distribution");
expectedOutputWeigthDistribution := [ <0, 1>, <8, 1>, <10, 1>, <11, 1>, <12, 16>, 
<13, 1>, <14, 17>, <15, 1>, <16, 33>, <17, 48>, <18, 52>, <19, 82>, <20, 111>, 
<21, 124>, <22, 192>, <23, 200>, <24, 225>, <25, 226>, <26, 191>, <27, 147>, 
<28, 112>, <29, 99>, <30, 55>, <31, 79>, <32, 12>, <33, 14>, <34, 4>, <35, 2>, 
<36, 1> ];
expectedOutputWeigthDistributionDual := 
    MacWilliamsTransform(alpha+2*beta, gamma+2*delta, 2, expectedOutputWeigthDistribution);
expectedOutputExternalDistance := 46;

assert Z2Z4MinimumLeeWeight(C) eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C) eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C) in expectedOutputMinimumWords;
assert MinimumWords(C) eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C) eq expectedOutputWeigthDistribution;
assert DualLeeWeightDistribution(C) eq expectedOutputWeigthDistributionDual;
assert ExternalDistance(C) eq expectedOutputExternalDistance;

assert Z2Z4MinimumLeeWeight(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Distribution") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Distribution") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "Distribution") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "KernelCosets") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "KernelCosets") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "KernelCosets") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Brouwer") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Brouwer") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Zimmermann") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Zimmermann") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Quaternary") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Quaternary") eq expectedOutputMinimumWords;

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
alpha := 10;
beta := 19;
gamma := 7;
delta := 2;
C := Z2Z4AdditiveCode(M, alpha);

expectedOutputMinimumLeeDistance :=  7;
expectedOutputMinimumWords :=  MinimumWords(C : Method := "Distribution");
expectedOutputWeigthDistribution := [ <0, 1>, <7, 1>, <10, 1>, <11, 6>, <12, 11>, 
<13, 8>, <14, 10>, <15, 11>, <16, 27>, <17, 41>, <18, 52>, <19, 91>, <20, 111>, 
<21, 138>, <22, 176>, <23, 193>, <24, 222>, <25, 220>, <26, 189>, <27, 158>, 
<28, 109>, <29, 94>, <30, 78>, <31, 51>, <32, 30>, <33, 11>, <34, 6>, <35, 1>, 
<36, 1> ];
expectedOutputWeigthDistributionDual := 
    MacWilliamsTransform(alpha+2*beta, gamma+2*delta, 2, expectedOutputWeigthDistribution);
expectedOutputExternalDistance := 45;

assert Z2Z4MinimumLeeWeight(C) eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C) eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C) in expectedOutputMinimumWords;
assert MinimumWords(C) eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C) eq expectedOutputWeigthDistribution;
assert DualLeeWeightDistribution(C) eq expectedOutputWeigthDistributionDual;
assert ExternalDistance(C) eq expectedOutputExternalDistance;

assert Z2Z4MinimumLeeWeight(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Distribution") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Distribution") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "Distribution") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "KernelCosets") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "KernelCosets") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "KernelCosets") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Brouwer") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Brouwer") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Zimmermann") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Zimmermann") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Quaternary") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Quaternary") eq expectedOutputMinimumWords;

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
alpha := 0;
beta := 19;
gamma := 7;
delta := 2;
C := Z2Z4AdditiveCode(M, alpha);

expectedOutputMinimumLeeDistance :=  6;
expectedOutputMinimumWords :=  MinimumWords(C : Method := "Distribution");
expectedOutputWeigthDistribution := [ <0, 1>, <6, 2>, <8, 8>, <10, 11>, <11, 12>, 
<12, 47>, <13, 60>, <14, 94>, <15, 116>, <16, 134>, <17, 204>, <18, 206>, 
<19, 252>, <20, 206>, <21, 180>, <22, 150>, <23, 124>, <24, 104>, <25, 68>, 
<26, 47>, <27, 8>, <28, 11>, <30, 2>, <32, 1> ];
expectedOutputWeigthDistributionDual := 
    MacWilliamsTransform(alpha+2*beta, gamma+2*delta, 2, expectedOutputWeigthDistribution);
expectedOutputExternalDistance := 35;

assert Z2Z4MinimumLeeWeight(C) eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C) eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C) in expectedOutputMinimumWords;
assert MinimumWords(C) eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C) eq expectedOutputWeigthDistribution;
assert DualLeeWeightDistribution(C) eq expectedOutputWeigthDistributionDual;
assert ExternalDistance(C) eq expectedOutputExternalDistance;

assert Z2Z4MinimumLeeWeight(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Distribution") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Distribution") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "Distribution") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "KernelCosets") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "KernelCosets") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "KernelCosets") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Brouwer") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Brouwer") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Zimmermann") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Zimmermann") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Quaternary") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Quaternary") eq expectedOutputMinimumWords;

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
alpha := 10;
beta := 0;
gamma := 5;
delta := 0;
C := Z2Z4AdditiveCode(M, alpha);

expectedOutputMinimumLeeDistance :=  1;
expectedOutputMinimumWords :=  MinimumWords(C : Method := "Distribution");
expectedOutputWeigthDistribution := [ <0, 1>, <1, 1>, <2, 1>, <3, 2>, <4, 4>, 
<5, 9>, <6, 9>, <7, 4>, <8, 1> ];
expectedOutputWeigthDistributionDual := 
    MacWilliamsTransform(alpha+2*beta, gamma+2*delta, 2, expectedOutputWeigthDistribution);
expectedOutputExternalDistance := 7;

assert Z2Z4MinimumLeeWeight(C) eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C) eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C) in expectedOutputMinimumWords;
assert MinimumWords(C) eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C) eq expectedOutputWeigthDistribution;
assert DualLeeWeightDistribution(C) eq expectedOutputWeigthDistributionDual;
assert ExternalDistance(C) eq expectedOutputExternalDistance;

assert Z2Z4MinimumLeeWeight(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Distribution") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Distribution") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "Distribution") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "KernelCosets") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "KernelCosets") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "KernelCosets") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Brouwer") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Brouwer") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Zimmermann") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Zimmermann") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Quaternary") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Quaternary") eq expectedOutputMinimumWords;

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
alpha := 10;
beta := 19;
gamma := 0;
delta := 7;
C := Z2Z4AdditiveCode(M, alpha);

expectedOutputMinimumLeeDistance :=  10;
expectedOutputMinimumWords :=  MinimumWords(C : Method := "Distribution");
expectedOutputWeigthDistribution := [ <0, 1>, <10, 1>, <12, 9>, <13, 16>, <14, 33>, 
<15, 64>, <16, 132>, <17, 264>, <18, 404>, <19, 650>, <20, 936>, <21, 1298>, 
<22, 1636>, <23, 1780>, <24, 1834>, <25, 1778>, <26, 1691>, <27, 1358>, <28, 951>, 
<29, 694>, <30, 419>, <31, 226>, <32, 97>, <33, 46>, <34, 40>, <35, 16>, <36, 8>, 
<39, 2> ];
expectedOutputWeigthDistributionDual := 
    MacWilliamsTransform(alpha+2*beta, gamma+2*delta, 2, expectedOutputWeigthDistribution);
expectedOutputExternalDistance := 43;

assert Z2Z4MinimumLeeWeight(C) eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C) eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C) in expectedOutputMinimumWords;
assert MinimumWords(C) eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C) eq expectedOutputWeigthDistribution;
assert DualLeeWeightDistribution(C) eq expectedOutputWeigthDistributionDual;
assert ExternalDistance(C) eq expectedOutputExternalDistance;

assert Z2Z4MinimumLeeWeight(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Distribution") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Distribution") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "Distribution") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "KernelCosets") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "KernelCosets") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "KernelCosets") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Brouwer") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Brouwer") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Zimmermann") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Zimmermann") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Quaternary") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Quaternary") eq expectedOutputMinimumWords;

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
alpha := 10;
beta := 19;
gamma := 7;
delta := 0;
C := Z2Z4AdditiveCode(M, alpha);

expectedOutputMinimumLeeDistance :=  12;
expectedOutputMinimumWords :=  MinimumWords(C : Method := "Distribution");
expectedOutputWeigthDistribution := [ <0, 1>, <12, 1>, <14, 4>, <15, 3>, <16, 5>, 
<17, 2>, <18, 6>, <19, 6>, <20, 12>, <21, 9>, <22, 11>, <23, 7>, <24, 7>, <25, 14>, 
<26, 8>, <27, 12>, <28, 5>, <29, 1>, <30, 3>, <31, 4>, <32, 1>, <33, 6> ];
expectedOutputWeigthDistributionDual := 
    MacWilliamsTransform(alpha+2*beta, gamma+2*delta, 2, expectedOutputWeigthDistribution);
expectedOutputExternalDistance := 47;

assert Z2Z4MinimumLeeWeight(C) eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C) eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C) in expectedOutputMinimumWords;
assert MinimumWords(C) eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C) eq expectedOutputWeigthDistribution;
assert DualLeeWeightDistribution(C) eq expectedOutputWeigthDistributionDual;
assert ExternalDistance(C) eq expectedOutputExternalDistance;

assert Z2Z4MinimumLeeWeight(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Distribution") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Distribution") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "Distribution") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "KernelCosets") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "KernelCosets") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "KernelCosets") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Brouwer") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Brouwer") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Zimmermann") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Zimmermann") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Quaternary") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Quaternary") eq expectedOutputMinimumWords;

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
alpha := 10;
beta := 19;
gamma := 7;
delta := 2;
C := Z2Z4AdditiveCode(M, alpha);

expectedOutputMinimumLeeDistance :=  6;
expectedOutputMinimumWords :=  MinimumWords(C : Method := "Distribution");
expectedOutputWeigthDistribution := [ <0, 1>, <6, 2>, <8, 2>, <10, 17>, <12, 29>, 
<14, 64>, <16, 114>, <18, 162>, <19, 44>, <20, 210>, <21, 92>, <22, 178>, <23, 148>, 
<24, 130>, <25, 236>, <26, 85>, <27, 220>, <28, 25>, <29, 148>, <30, 4>, <31, 92>, 
<32, 1>, <33, 36>, <35, 8> ];
expectedOutputWeigthDistributionDual := 
    MacWilliamsTransform(alpha+2*beta, gamma+2*delta, 2, expectedOutputWeigthDistribution);
expectedOutputExternalDistance := 47;

assert Z2Z4MinimumLeeWeight(C) eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C) eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C) in expectedOutputMinimumWords;
assert MinimumWords(C) eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C) eq expectedOutputWeigthDistribution;
assert DualLeeWeightDistribution(C) eq expectedOutputWeigthDistributionDual;
assert ExternalDistance(C) eq expectedOutputExternalDistance;

assert Z2Z4MinimumLeeWeight(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Distribution") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Distribution") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "Distribution") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "KernelCosets") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "KernelCosets") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "KernelCosets") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Brouwer") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Brouwer") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Zimmermann") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Zimmermann") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Quaternary") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Quaternary") eq expectedOutputMinimumWords;

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
alpha := 7;
beta := 12;
gamma := 8;
delta := 2;
C := Z2Z4AdditiveCode(M, alpha);

expectedOutputMinimumLeeDistance :=  5;
expectedOutputMinimumWords :=  MinimumWords(C : Method := "Distribution");
expectedOutputWeigthDistribution := [ <0, 1>, <5, 3>, <6, 6>, <7, 13>, <8, 19>, 
<9, 53>, <10, 103>, <11, 178>, <12, 280>, <13, 366>, <14, 477>, <15, 531>, 
<16, 542>, <17, 513>, <18, 379>, <19, 273>, <20, 171>, <21, 82>, <22, 57>, 
<23, 28>, <24, 10>, <25, 6>, <26, 2>, <27, 1>, <28, 1>, <29, 1> ];
expectedOutputWeigthDistributionDual := 
    MacWilliamsTransform(alpha+2*beta, gamma+2*delta, 2, expectedOutputWeigthDistribution);
expectedOutputExternalDistance := 24;

assert Z2Z4MinimumLeeWeight(C) eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C) eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C) in expectedOutputMinimumWords;
assert MinimumWords(C) eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C) eq expectedOutputWeigthDistribution;
assert DualLeeWeightDistribution(C) eq expectedOutputWeigthDistributionDual;
assert ExternalDistance(C) eq expectedOutputExternalDistance;

assert Z2Z4MinimumLeeWeight(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Distribution") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Distribution") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "Distribution") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "KernelCosets") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "KernelCosets") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "KernelCosets") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Brouwer") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Brouwer") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Zimmermann") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Zimmermann") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Quaternary") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Quaternary") eq expectedOutputMinimumWords;

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
alpha := 9;
beta := 12;
gamma := 8;
delta := 2;
C := Z2Z4AdditiveCode(M, alpha);

expectedOutputMinimumLeeDistance :=  5;
expectedOutputMinimumWords :=  MinimumWords(C : Method := "Distribution");
expectedOutputWeigthDistribution := [ <0, 1>, <5, 1>, <6, 1>, <7, 6>, <8, 8>, 
<9, 20>, <10, 47>, <11, 105>, <12, 206>, <13, 278>, <14, 356>, <15, 483>, 
<16, 519>, <17, 526>, <18, 528>, <19, 391>, <20, 258>, <21, 157>, <22, 83>, 
<23, 71>, <24, 32>, <25, 10>, <26, 9> ];
expectedOutputWeigthDistributionDual := 
    MacWilliamsTransform(alpha+2*beta, gamma+2*delta, 2, expectedOutputWeigthDistribution);
expectedOutputExternalDistance := 27;

assert Z2Z4MinimumLeeWeight(C) eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C) eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C) in expectedOutputMinimumWords;
assert MinimumWords(C) eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C) eq expectedOutputWeigthDistribution;
assert DualLeeWeightDistribution(C) eq expectedOutputWeigthDistributionDual;
assert ExternalDistance(C) eq expectedOutputExternalDistance;

assert Z2Z4MinimumLeeWeight(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Distribution") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Distribution") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "Distribution") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "KernelCosets") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "KernelCosets") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "KernelCosets") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Brouwer") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Brouwer") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Zimmermann") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Zimmermann") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Quaternary") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Quaternary") eq expectedOutputMinimumWords;

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
alpha := 9;
beta := 12;
gamma := 6;
delta := 3;
C := Z2Z4AdditiveCode(M, alpha);

expectedOutputMinimumLeeDistance :=  2;
expectedOutputMinimumWords :=  MinimumWords(C : Method := "Distribution");
expectedOutputWeigthDistribution := [ <0, 1>, <2, 1>, <6, 1>, <8, 12>, <9, 20>, 
<10, 68>, <11, 156>, <12, 240>, <13, 324>, <14, 414>, <15, 508>, <16, 553>, 
<17, 540>, <18, 469>, <19, 340>, <20, 200>, <21, 140>, <22, 65>, <23, 20>, 
<24, 18>, <26, 6> ];
expectedOutputWeigthDistributionDual := 
    MacWilliamsTransform(alpha+2*beta, gamma+2*delta, 2, expectedOutputWeigthDistribution);
expectedOutputExternalDistance := 31;

assert Z2Z4MinimumLeeWeight(C) eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C) eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C) in expectedOutputMinimumWords;
assert MinimumWords(C) eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C) eq expectedOutputWeigthDistribution;
assert DualLeeWeightDistribution(C) eq expectedOutputWeigthDistributionDual;
assert ExternalDistance(C) eq expectedOutputExternalDistance;

assert Z2Z4MinimumLeeWeight(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Distribution") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Distribution") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "Distribution") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "KernelCosets") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "KernelCosets") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "KernelCosets") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Brouwer") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Brouwer") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Zimmermann") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Zimmermann") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Quaternary") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Quaternary") eq expectedOutputMinimumWords;

/************************************************************/
print "test 16: A Z2Z4-additive code with a zero column in beta part.
                alpha = 9, beta = 12, gamma = 6, delta = 3, kappa = 4, length = 21, #C = 4096";
M := Matrix(Z4,[[2,0,0,0,2,2,2,0,2,0,0,0,0,2,0,0,0,0,0,0,2],
                [2,0,0,2,0,0,2,0,2,0,0,2,0,2,0,0,2,0,0,2,0],
                [0,0,2,2,0,2,2,0,0,0,0,0,0,2,2,0,0,0,2,2,0],
                [0,0,2,2,0,2,2,0,0,0,0,0,0,0,0,2,2,0,2,2,0],
                [0,2,0,0,2,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0],
                [2,0,0,2,0,0,2,0,2,0,0,0,0,0,0,0,0,0,0,2,0],
                [0,0,0,0,0,2,2,2,2,0,0,0,1,0,1,1,0,0,1,2,1],
                [0,0,2,0,2,2,0,0,2,1,0,1,0,2,0,1,0,0,0,3,0],
                [0,0,0,2,0,2,2,0,2,0,1,1,0,1,0,0,2,0,1,0,0]]);
alpha := 9;
beta := 12;
gamma := 6;
delta := 3;
C := Z2Z4AdditiveCode(M, alpha);

expectedOutputMinimumLeeDistance :=  3;
expectedOutputMinimumWords :=  MinimumWords(C : Method := "Distribution");
expectedOutputWeigthDistribution := [ <0, 1>, <3, 1>, <6, 5>, <7, 1>, <8, 18>, 
<9, 33>, <10, 90>, <11, 185>, <12, 294>, <13, 422>, <14, 420>, <15, 508>, 
<16, 603>, <17, 536>, <18, 406>, <19, 245>, <20, 170>, <21, 90>, <22, 39>, 
<23, 19>, <24, 2>, <25, 7>, <27, 1> ];
expectedOutputWeigthDistributionDual := 
    MacWilliamsTransform(alpha+2*beta, gamma+2*delta, 2, expectedOutputWeigthDistribution);
expectedOutputExternalDistance := 30;

assert Z2Z4MinimumLeeWeight(C) eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C) eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C) in expectedOutputMinimumWords;
assert MinimumWords(C) eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C) eq expectedOutputWeigthDistribution;
assert DualLeeWeightDistribution(C) eq expectedOutputWeigthDistributionDual;
assert ExternalDistance(C) eq expectedOutputExternalDistance;

assert Z2Z4MinimumLeeWeight(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Distribution") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Distribution") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "Distribution") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "KernelCosets") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "KernelCosets") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "KernelCosets") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Brouwer") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Brouwer") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Zimmermann") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Zimmermann") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Quaternary") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Quaternary") eq expectedOutputMinimumWords;

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
alpha := 9;
beta := 12;
gamma := 6;
delta := 3;
C := Z2Z4AdditiveCode(M, alpha);

expectedOutputMinimumLeeDistance :=  2;
expectedOutputMinimumWords :=  MinimumWords(C : Method := "Distribution");
expectedOutputWeigthDistribution := [ <0, 1>, <2, 2>, <3, 4>, <4, 5>, <5, 6>, 
<6, 18>, <7, 40>, <8, 114>, <9, 188>, <10, 340>, <11, 480>, <12, 554>, <13, 632>, 
<14, 540>, <15, 440>, <16, 317>, <17, 196>, <18, 122>, <19, 60>, <20, 33>, <21, 2>, 
<22, 2> ];
expectedOutputWeigthDistributionDual := 
    MacWilliamsTransform(alpha+2*beta, gamma+2*delta, 2, expectedOutputWeigthDistribution);
expectedOutputExternalDistance := 31;

assert Z2Z4MinimumLeeWeight(C) eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C) eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C) in expectedOutputMinimumWords;
assert MinimumWords(C) eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C) eq expectedOutputWeigthDistribution;
assert DualLeeWeightDistribution(C) eq expectedOutputWeigthDistributionDual;
assert ExternalDistance(C) eq expectedOutputExternalDistance;

assert Z2Z4MinimumLeeWeight(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Distribution") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Distribution") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "Distribution") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "KernelCosets") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "KernelCosets") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "KernelCosets") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Brouwer") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Brouwer") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Zimmermann") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Zimmermann") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Quaternary") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Quaternary") eq expectedOutputMinimumWords;

/************************************************************/						   
print "test 18: A Z2Z4-additive code 
                alpha = 3, beta = 5, gamma = 1, delta = 2, kappa = 0, length = 8, #C = 32";
M := Matrix(Z4,[[2,0,2,1,1,1,2,1],
                [0,0,0,2,0,1,1,1],
                [0,0,0,0,2,1,3,1],
                [0,0,0,0,0,2,2,0],
                [0,0,0,0,0,0,0,2]]);
alpha := 3;
beta := 5;
gamma := 1;
delta := 2;
C := Z2Z4AdditiveCode(M, alpha);

expectedOutputMinimumLeeDistance :=  2;
expectedOutputMinimumWords :=  MinimumWords(C : Method := "Distribution");
expectedOutputWeigthDistribution := [ <0, 1>, <2, 1>, <4, 1>, <5, 10>, <6, 7>, <7, 4>, <8, 6>, <9, 2> ];
expectedOutputWeigthDistributionDual := 
    MacWilliamsTransform(alpha+2*beta, gamma+2*delta, 2, expectedOutputWeigthDistribution);
expectedOutputExternalDistance := 11;

assert Z2Z4MinimumLeeWeight(C) eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C) eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C) in expectedOutputMinimumWords;
assert MinimumWords(C) eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C) eq expectedOutputWeigthDistribution;
assert DualLeeWeightDistribution(C) eq expectedOutputWeigthDistributionDual;
assert ExternalDistance(C) eq expectedOutputExternalDistance;

assert Z2Z4MinimumLeeWeight(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Distribution") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Distribution") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "Distribution") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "KernelCosets") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "KernelCosets") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "KernelCosets") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Brouwer") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Brouwer") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Zimmermann") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Zimmermann") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Quaternary") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Quaternary") eq expectedOutputMinimumWords;

/************************************************************/						   
print "test 19: A Z2Z4-additive code 
                alpha = 2, beta = 5, gamma = 2, delta = 4, kappa = 1, length = 7, #C = 1024";
M := Matrix(Z4,[[2,0,0,0,0,0,0],
                [0,2,0,0,0,0,1],
                [0,0,1,0,0,1,0],
                [0,0,0,1,0,1,0],
                [0,0,0,0,1,0,1],
                [0,0,0,0,0,2,0],
                [0,0,0,0,0,0,2]]);
alpha := 2;
beta := 5;
gamma := 2;
delta := 4;
C := Z2Z4AdditiveCode(M, alpha);

expectedOutputMinimumLeeDistance :=  1;
expectedOutputMinimumWords :=  MinimumWords(C : Method := "Distribution");
expectedOutputWeigthDistribution := [ <0, 1>, <1, 1>, <2, 25>, <3, 25>, <4, 170>, <5, 170>, 
<6, 226>, <7, 226>, <8, 85>, <9, 85>, <10, 5>, <11, 5> ];
expectedOutputWeigthDistributionDual := 
    MacWilliamsTransform(alpha+2*beta, gamma+2*delta, 2, expectedOutputWeigthDistribution);
expectedOutputExternalDistance := 3;

assert Z2Z4MinimumLeeWeight(C) eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C) eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C) in expectedOutputMinimumWords;
assert MinimumWords(C) eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C) eq expectedOutputWeigthDistribution;
assert DualLeeWeightDistribution(C) eq expectedOutputWeigthDistributionDual;
assert ExternalDistance(C) eq expectedOutputExternalDistance;

assert Z2Z4MinimumLeeWeight(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Distribution") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Distribution") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "Distribution") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "KernelCosets") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "KernelCosets") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "KernelCosets") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Brouwer") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Brouwer") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Zimmermann") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Zimmermann") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Quaternary") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Quaternary") eq expectedOutputMinimumWords;

/************************************************************/						   
print "test 20: A Z2Z4-additive code 
                alpha = 1, beta = 9, gamma = 3, delta = 4, kappa = 1, length = 10, #C = 2048";
M := Matrix(Z4,[[2,0,0,0,0,1,0,3,1,3],
                [0,1,0,0,0,0,0,0,0,1],
                [0,0,2,0,0,1,0,1,1,3],
                [0,0,0,1,0,1,0,0,1,0],
                [0,0,0,0,2,1,0,1,1,1],
                [0,0,0,0,0,2,0,2,0,2],
                [0,0,0,0,0,0,1,0,1,2],
                [0,0,0,0,0,0,0,0,2,0]]);
alpha := 1;
beta := 9;
gamma := 3;
delta := 4;
C := Z2Z4AdditiveCode(M, alpha);

expectedOutputMinimumLeeDistance :=  2;
expectedOutputMinimumWords :=  MinimumWords(C : Method := "Distribution");
expectedOutputWeigthDistribution := [ <0, 1>, <2, 4>, <3, 4>, <4, 27>, <5, 59>, <6, 104>, <7, 208>, 
<8, 293>, <9, 355>, <10, 316>, <11, 236>, <12, 217>, <13, 129>, <14, 56>, <15, 32>, <16, 6>, <17, 1> ];
expectedOutputWeigthDistributionDual := 
    MacWilliamsTransform(alpha+2*beta, gamma+2*delta, 2, expectedOutputWeigthDistribution);
expectedOutputExternalDistance := 14;

assert Z2Z4MinimumLeeWeight(C) eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C) eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C) in expectedOutputMinimumWords;
assert MinimumWords(C) eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C) eq expectedOutputWeigthDistribution;
assert DualLeeWeightDistribution(C) eq expectedOutputWeigthDistributionDual;
assert ExternalDistance(C) eq expectedOutputExternalDistance;

assert Z2Z4MinimumLeeWeight(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Distribution") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Distribution") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Distribution") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "Distribution") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "KernelCosets") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "KernelCosets") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "KernelCosets") eq expectedOutputMinimumWords;
assert LeeWeightDistribution(C : Method := "KernelCosets") eq expectedOutputWeigthDistribution;

assert Z2Z4MinimumLeeWeight(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Brouwer") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Brouwer") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Brouwer") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Zimmermann") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Zimmermann") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Zimmermann") eq expectedOutputMinimumWords;

assert Z2Z4MinimumLeeWeight(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert Z2Z4MinimumLeeDistance(C : Method := "Quaternary") eq expectedOutputMinimumLeeDistance;
assert MinimumWord(C : Method := "Quaternary") in expectedOutputMinimumWords;
assert MinimumWords(C : Method := "Quaternary") eq expectedOutputMinimumWords;

SetVerbose("IgnoreWeightAttributes", 0);