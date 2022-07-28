/*************************************************************/
/*                                                           */
/* Package or project name: Z2Z4AdditiveCodes package        */
/* Test file name: Z2Z4Families_BB_test.m      	             */
/*                                                           */
/* Comments: Black-box tests for the intrinsic functions     */
/*        	 Z2Z4ReedMullerCode, Z2Z4ExtendedPerfectCode,    */
/*           Z2Z4HadamardCode and Z2Z4ReedMullerCode included*/
/*           in the Z2Z4CodeConstructions.m file.            */
/*                                                           */
/* Authors: Jaume Pujol and Mercè Villanueva                 */
/*                                                           */
/* Revision version and last date: v1.0   2022/02/16         */
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
/* Function name: Z2Z4ExtendedPerfectCode                       */
/* Parameters: delta, m                                         */
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
/*   - m : A positive integer                                   */
/* Output parameters description:                               */
/*   - The Z2Z4-extended perfect additive code                  */
/*                                                              */
/* Signature: (<RngIntElt> delta, <RngIntElt> m ) -> Z2Z4Code   */
/*                                                              */
/****************************************************************/ 
/****************************************************************/
/*                                                              */
/* Function name: Z2Z4HadamardCode                              */
/* Parameters: delta, m                                         */
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
/*   - m : A positive integer                                   */
/* Output parameters description:                               */
/*   - The Z2Z4-Hadamard additive code                          */
/*                                                              */
/* Signature: (<RngIntElt> delta, <RngIntElt> m ) -> Z2Z4Code   */
/*                                                              */
/****************************************************************/ 
/****************************************************************/
/*                                                              */
/* Function name: Z2Z4ReedMullerCode                            */
/* Parameters: s, r, m                                          */
/* Function description: The function returns a Z2Z4-additive   */
/*   Reed-Muller code of type (alpha, beta; gamma, delta;kappa).*/
/*   The parameter OverZ4 specifies whether the code is over Z4,*/
/*   that is alpha=0, or, otherwise, alpha<>0. The default value*/
/*   is false. When OverZ4 is true, given an integer m>=1, an   */
/*   integer s such that 0 <= s <= Floor((m-1)/2) and an integer*/
/*   0 <= r <= m, return an r-th order Z2Z4-Additive Reed-Muller*/
/*   Code RM_s(r,m), with alpha = 0 and beta = 2^(m-1). When    */
/*   OverZ4 is false, given m>=1, 0 <= s <= Floor(m/2) and      */
/*   0 <= r <= m, return an r-th Z2Z4-Additive Reed-MullerCode  */
/*   ARM_s(r,m), with alpha = 2^(m-s) and beta =                */
/*   2^(m-1)-2^(m-s-1). In both cases, it returns a generator   */
/*   matrix with gamma+delta rows constructed in a recursive    */
/*   way, where the ones in the first alpha coordinates are     */
/*   represented by twos.                                       */
/* Input parameters description:                                */
/*   - s : An integer between 0 and Floor(m/2)                  */
/*   - r : An integer between 0 and m                           */
/*   - m : A positive integer                                   */
/* Output parameters description:                               */
/*   - The Z2Z4-Additive Reed-Muller code                       */
/*   - A generator matrix                                       */
/*                                                              */
/* Signature: (<RngIntElt> s, <RngIntElt> r, <RngIntElt> m)     */
/*                                -> Z2Z4Code, ModMatRngElt     */
/*                                                              */
/****************************************************************/
/****************************************************************/
/*                                                              */
/* Function name: Z2Z4ReedMullerCodes                           */
/* Parameters: s, m                                             */
/* Function description: The function returns a sequence        */
/*   containing a family of Z2Z4-additive Reed-Muller codes. The*/
/*   parameter OverZ4 specifies whether the codes are over Z4,  */
/*   that is alpha=0, or, otherwise, alpha<>0. The default value*/
/*   is false. When OverZ4 is true, return the codes RM_s(r,m), */
/*   given any m>=1 and 0<=s<=Floor((m-1)/2). When OverZ4 is    */
/*   false, return the codes ARM_s(r,m), given any m>=1 and     */
/*   0<=s<=Floor(m/2).                                          */
/* Input parameters description:                                */
/*   - s : An integer between 0 and Floor(m/2)                  */
/*   - r : An integer between 0 and m                           */
/*   - m : A positive integer                                   */
/* Output parameters description:                               */
/*   - The Z2Z4-Additive Reed-Muller code                       */
/*   - A generator matrix                                       */
/*                                                              */
/* Signature: (<RngIntElt> s, <RngIntElt> m ) -> Z2Z4Code,      */
/*                                            ModMatRngElt      */
/*                                                              */
/****************************************************************/
print "ADDITIVE REED-MULLER CODES";

CheckARM := procedure(s,m)
    ARM := [];

    for i in [1..(m+1)] do
        ARM[i] := Z2Z4ReedMullerCode(s,i-1,m);
    end for;

    for i in [2..(m+1)] do
        assert ARM[i-1] subset ARM[i];
    end for;

    for i in [1..(m+1)] do

        expectedBinaryLength := 2^m;
        expectedSize := 2^&+[Binomial(m,r) : r in [0..(i-1)]];
        expectedMinimumLeeWeight := 2^(m-i+1);

        _, alphabeta := Length(ARM[i]);
        binaryLength := alphabeta[1]+2*alphabeta[2];
        size := #ARM[i];
        minimumLeeWeight := Z2Z4MinimumLeeWeight(ARM[i]);

        assert binaryLength eq expectedBinaryLength;
        assert size eq expectedSize;
        assert minimumLeeWeight eq expectedMinimumLeeWeight;

    end for;

    assert ARM[1] eq Z2Z4AdditiveRepetitionCode(alphabeta[1], alphabeta[2]);

    hadamardDistribution := LeeWeightDistribution(Z2Z4HadamardCode(s, m));
    distribution := LeeWeightDistribution(ARM[2]);
    assert distribution eq hadamardDistribution;
     
    extendedPerfectDistribution := LeeWeightDistribution(Z2Z4ExtendedPerfectCode(s, m));
    distribution := LeeWeightDistribution(ARM[m-1]);
    assert distribution eq extendedPerfectDistribution;

    assert ARM[m] eq Z2Z4AdditiveEvenWeightCode(alphabeta[1], alphabeta[2]);
    
    assert ARM[m+1] eq Z2Z4AdditiveUniverseCode(alphabeta[1], alphabeta[2]);


end procedure;

print "test 1: ARM_s(r,m) family with s=1, m=5 and 0<=r<=m";

s:=1;
m:=5;

CheckARM(s,m);

print "test 2: ARM_s(r,m) family with s=3, m=6 and 0<=r<=m";

s:=3;
m:=6;

CheckARM(s, m);

print "test 3: ARM_s(r,m) family with s=0, m=8 and 0<=r<=m";

s:=0;
m:=8;

CheckARM(s, m);

print "QUATERNARY REED-MULLER CODES";

CheckRM := procedure(s,m)
    RM := [];

    for i in [1..(m+1)] do
        RM[i] := Z2Z4ReedMullerCode(s,i-1,m : OverZ4 := true);
    end for;

    for i in [2..(m+1)] do
        assert RM[i-1] subset RM[i];
    end for;

    for i in [1..(m+1)] do

        expectedBinaryLength := 2^m;
        expectedSize := 2^&+[Binomial(m,r) : r in [0..(i-1)]];
        expectedMinimumLeeWeight := 2^(m-i+1);

        _, alphabeta := Length(RM[i]);
        binaryLength := alphabeta[1]+2*alphabeta[2];
        size := #RM[i];
        minimumLeeWeight := Z2Z4MinimumLeeWeight(RM[i]);

        assert alphabeta[1] eq 0;
        assert binaryLength eq expectedBinaryLength;
        assert size eq expectedSize;
        assert minimumLeeWeight eq expectedMinimumLeeWeight;

    end for;

    assert RM[1] eq Z2Z4AdditiveRepetitionCode(alphabeta[1], alphabeta[2]);

    hadamardDistribution := LeeWeightDistribution(Z2Z4HadamardCode(s, m));
    distribution := LeeWeightDistribution(RM[2]);
    assert distribution eq hadamardDistribution;
    
    extendedPerfectDistribution := LeeWeightDistribution(Z2Z4ExtendedPerfectCode(s, m));
    distribution := LeeWeightDistribution(RM[m-1]);
    assert distribution eq extendedPerfectDistribution;
    
    assert RM[m] eq Z2Z4AdditiveEvenWeightCode(alphabeta[1], alphabeta[2]);
    
    assert RM[m+1] eq Z2Z4AdditiveUniverseCode(alphabeta[1], alphabeta[2]);

end procedure;

print "test 4: RM_s(r,m) (alpha=0) family with s=1, m=5 and 0<=r<=m";

s:=0;
m:=5;

CheckRM(s,m);

print "test 5: RM_s(r,m) (alpha=0) family with s=2, m=5 and 0<=r<=m";

s:=2;
m:=5;

CheckRM(s,m);

print "test 6: RM_s(r,m) (alpha=0) family with s=3, m=7 and 0<=r<=m";

s:=3;
m:=7;

CheckRM(s,m);

print "test 7: RM_s(r,m) (alpha=0) family with s=1, m=7 and 0<=r<=m";

s:=1;
m:=7;

CheckRM(s,m);

print "FAMILY OF ADDITIVE REED-MULLER CODES";

CheckARMs := procedure(s,m)

    ARM := Z2Z4ReedMullerCodes(s,m);

    for i in [2..(m+1)] do
        assert ARM[i-1] subset ARM[i];
    end for;

    for i in [1..(m+1)] do

        expectedBinaryLength := 2^m;
        expectedSize := 2^&+[Binomial(m,r) : r in [0..(i-1)]];
        expectedMinimumLeeWeight := 2^(m-i+1);

        _, alphabeta := Length(ARM[i]);
        binaryLength := alphabeta[1]+2*alphabeta[2];
        size := #ARM[i];
        minimumLeeWeight := Z2Z4MinimumLeeWeight(ARM[i]);

        assert binaryLength eq expectedBinaryLength;
        assert size eq expectedSize;
        assert minimumLeeWeight eq expectedMinimumLeeWeight;

    end for;

    assert ARM[1] eq Z2Z4AdditiveRepetitionCode(alphabeta[1], alphabeta[2]);

    hadamardDistribution := LeeWeightDistribution(Z2Z4HadamardCode(s, m));
    distribution := LeeWeightDistribution(ARM[2]);
    assert distribution eq hadamardDistribution;
     
    extendedPerfectDistribution := LeeWeightDistribution(Z2Z4ExtendedPerfectCode(s, m));
    distribution := LeeWeightDistribution(ARM[m-1]);
    assert distribution eq extendedPerfectDistribution;

    assert ARM[m] eq Z2Z4AdditiveEvenWeightCode(alphabeta[1], alphabeta[2]);
    
    assert ARM[m+1] eq Z2Z4AdditiveUniverseCode(alphabeta[1], alphabeta[2]);


end procedure;

print "test 8: ARM_s(r,m) family with s=2, m=5";

s:=2;
m:=5;

CheckARMs(s,m);

print "test 9: ARM_s(r,m) family with s=0, m=6";

s:=3;
m:=6;

CheckARMs(s, m);

print "FAMILY OF QUATERNARY REED-MULLER CODES";

CheckRMs := procedure(s,m)
    RM := Z2Z4ReedMullerCodes(s,m : OverZ4 := true);

    for i in [2..(m+1)] do
        assert RM[i-1] subset RM[i];
    end for;

    for i in [1..(m+1)] do

        expectedBinaryLength := 2^m;
        expectedSize := 2^&+[Binomial(m,r) : r in [0..(i-1)]];
        expectedMinimumLeeWeight := 2^(m-i+1);

        _, alphabeta := Length(RM[i]);
        binaryLength := alphabeta[1]+2*alphabeta[2];
        size := #RM[i];
        minimumLeeWeight := Z2Z4MinimumLeeWeight(RM[i]);

        assert alphabeta[1] eq 0;
        assert binaryLength eq expectedBinaryLength;
        assert size eq expectedSize;
        assert minimumLeeWeight eq expectedMinimumLeeWeight;

    end for;

    assert RM[1] eq Z2Z4AdditiveRepetitionCode(alphabeta[1], alphabeta[2]);

    hadamardDistribution := LeeWeightDistribution(Z2Z4HadamardCode(s, m));
    distribution := LeeWeightDistribution(RM[2]);
    assert distribution eq hadamardDistribution;
    
    extendedPerfectDistribution := LeeWeightDistribution(Z2Z4ExtendedPerfectCode(s, m));
    distribution := LeeWeightDistribution(RM[m-1]);
    assert distribution eq extendedPerfectDistribution;
    
    assert RM[m] eq Z2Z4AdditiveEvenWeightCode(alphabeta[1], alphabeta[2]);
    
    assert RM[m+1] eq Z2Z4AdditiveUniverseCode(alphabeta[1], alphabeta[2]);

end procedure;

print "test 10: RM_s(r,m) (alpha=0) family with s=0, m=8";

s:=0;
m:=6;

CheckRMs(s, m);

print "test 11: RM_s(r,m) (alpha=0) family with s=3, m=7";

s:=3;
m:=7;

CheckRMs(s, m);