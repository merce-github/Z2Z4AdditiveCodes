/****************************************************************/
/*                                                              */
/* Package or project name: Z2Z4AdditiveCodes package           */
/* Test file name: Z2Z4Decoding_BB_test.m                       */
/*                                                              */
/* Comments: Black-box tests for the intrinsic functions        */
/*           InformationSpace, InformationSet,                  */
/*           IsInformationSet, SyndromeSpace,                   */
/*           Syndrome, CosetLeaders, CosetDecode                */
/*           and SyndromeDecode                                 */
/*           included in the Z2Z4Decode.m file                  */
/*                                                              */
/* Authors: J. Pujol and M. Villanueva                          */
/*                                                              */
/* Revision version and last date: v1.0, 2018/03/03             */
/*                                 v1.1, 2018/03/06             */
/*     User Defined type           v1.2, 01-02-2019             */
/*                                                              */
/****************************************************************/

//needs Z2Z4AdditiveCode file
//needs Z2Z4StandardForm file
//needs Z2Z4MinimumWeight file
//needs Z2Z4Decode file

SetAssertions(true);
Alarm(30*60);

/****************************************************************
	GLOBAL VARIABLES
*****************************************************************/	
Z4 := Integers(4);

/****************************************************************/
/*                                                              */
/* Function name: InformationSpace                              */
/* Parameters: C                                                */
/* Function description: Given a Z2Z4-additive code C of type   */
/*   (alpha,beta;gamma;delta;kappa), return the Z4-submodule of */
/*   Z4^(gamma+delta) isomorphic to Z2^gammaxZ4^delta such that */
/*   the first gamma coordinates are of order two, that is, the */
/*   space of information vectors for C. The function also      */
/*   returns the (gamma+2delta)-dimensional binary vector space,*/
/*   which is the space of information vectors for the corres-  */
/*   ponding binary code Cbin=Phi(C), where Phi is the Gray map.*/
/*   Finally, for the encoding process, it also returns the     */
/*   corresponding isomorphisms f and fbin from these spaces of */
/*   information vectors onto C and Cbin, respectively.         */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/* Output parameters description:                               */
/*   - Z4-submodule of length gamma+delta                       */
/*   - Binary vector space of dimension gamma+2delta            */
/*   - Map encoding information vectors over Z4 to C            */
/*   - Map encoding information binary vectors to Cbin          */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> ModTupRng, ModTupFld, Map, Map  */
/*                                                              */
/****************************************************************/
/****************************************************************/
/*                                                              */
/* Function name: InformationSet                                */
/* Parameters: C                                                */
/* Function description: Given a Z2Z4-additive code C of type   */
/*   (alpha,beta;gamma;delta;kappa), return an information set  */
/*   I=[i_1,...,i_(gamma+delta)] in [1,...,alpha+beta] for C    */
/*   such that [i_1,...,i_k] in [1,...,alpha] and the code C    */
/*   punctured on [1,...,alpha+beta] minus [i_(gamma+1),...,    */
/*   i_(gamma+delta)] is of type 4^delta, and the corresponding */
/*   information set Phi(I)=[i_1,...,i_k,2*i_(k+1)-1-alpha,..., */
/*   2*i_gamma-1-alpha,2*i_(gamma+1)-1-alpha,2*i_(gamma+1)-alpha*/
/*   ..., 2*i_(gamma+delta)-1-alpha, 2*i_(gamma+delta)-alpha]   */
/*   in [1,...,2*alpha+beta] for the binary code Cbin=Phi(C),   */
/*   where Phi is the Gray map. The information sets I and      */
/*   Phi(I) are returned as a sequence of gamma+delta and       */
/*   gamma+2*delta integers, giving the coordinate positions    */
/*   that correspond to the information set of C and Cbin,      */
/*   respectively.                                              */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/* Output parameters description:                               */
/*   - Sequence I with gamma+delta integers                     */
/*   - Sequence Ibin with gamma+2*delta integers                */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> SeqEnum[RngIntElt],             */
/*                              SeqEnum[RngIntElt]              */
/*                                                              */
/****************************************************************/
/****************************************************************/
/*                                                              */
/* Function name: IsInformationSet                              */
/* Parameters: C, I                                             */
/* Function description: Given a Z2Z4-additive code C of type   */
/*   (alpha,beta;gamma;delta;kappa) and a sequence I in [1,..., */
/*   alpha+beta] or I in [1,...,alpha+2*beta], return true if   */
/*   and only if I in [1,...,alpha+beta] is an information set  */
/*   for C. This function also returns another boolean, which   */
/*   is true if an only if I in [1,...,alpha+2*beta] is an      */
/*   information set for the corresponding binary code          */
/*   Cbin=Phi(C), where Phi is the Gray map.                    */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/*   - I : Sequence of integers in [1..alpha+beta] or           */
/*                                 [1..alpha+2*beta]            */
/* Output parameters description:                               */
/*   - Boolean, true if I is an information set for C           */
/*   - Boolean, true if I is an information set for Cbin        */
/*                                                              */
/* Signature: (<Z2Z4Code> C, <[RngIntElt]> I) -> BoolElt,       */
/*                                               BoolElt        */
/*                                                              */
/****************************************************************/
/****************************************************************/
/*                                                              */
/* Function name: SyndromeSpace                                 */
/* Parameters: C                                                */
/* Function description: Given a Z2Z4-additive code C of type   */
/*   (alpha,beta;gamma;delta;kappa), return the Z4-submodule of */
/*   Z4^(alpha+beta-gamma-kappa) isomorphic to Z2^(alpha+gamma  */
/*   -2*kappa) x Z4^(beta-gamma-delta+kappa) such that the first*/
/*   alpha+gamma-2*kappa coordinates are of order two, that is, */
/*   the space of syndrome vectors for C. The function also     */
/*   returns the (alpha+2*beta-gamma-2*delta)-dimensional binary*/
/*   vector space, which is the space of syndrome vectors for   */
/*   the corresponding binary code Cbin=Phi(C), where Phi is the*/
/*   Gray map. Note that these spaces are computed by using the */
/*   function InformationSpace(C) applied to the additive dual  */
/*   code of C, given by function Dual(C).                      */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/* Output parameters description:                               */
/*   - Z4-submodule of length alpha+beta-delta-kappa            */
/*   - Binary vector space of dimension alpha+2beta-gamma-2delta*/
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> ModTupRng, ModTupFld            */
/*                                                              */
/****************************************************************/
/****************************************************************/
/*                                                              */
/* Function name: Syndrome                                      */
/* Parameters: u, C                                             */
/* Function description: Given a Z2Z4-additive code C of type   */
/*   (alpha,beta;gamma;delta;kappa), and a vector u from the    */
/*   ambient space V=Z2^alpha x Z4^beta or V2=Z2^(alpha+2*beta),*/
/*   construct the syndrome of u relative to the code C, by     */
/*   using the parity check matrix given by the function        */
/*   MinRowsParityCheckMatrix(C). The vector u in V can be given*/
/*   either as a vector in Z4^(alpha+beta), where the ones in   */
/*   the first alpha coordinates are represented by twos; or as */
/*   a tuple in the cartesian product set Z2^alpha x Z4^beta.   */
/*   The syndrome is an element of the syndrome space of C,     */
/*   considered as the Z4-submodule of Z4^(alpha+beta-delta-    */
/*   kappa) isomorphic to Z2^(alpha+gamma-2*kappa) x Z4^(beta-  */
/*   gamma-delta+kappa) such that the first alpha+gamma-2*kappa */
/*   coordinates are of order two.                              */
/* Input parameters description:                                */
/*   - Vector over Z4 of length alpha+beta or                   */
/*     binary vector of length alpha+2beta                      */
/*   - C : A Z2Z4-additive code                                 */
/* Output parameters description:                               */
/*   - Vector over Z4 of length alpha+beta-delta-kappa          */
/*                                                              */
/* Signature: (<ModTupRngElt> u, <Z2Z4Code> C) -> ModTupRngElt  */
/* Signature: (<ModTupFldElt> u, <Z2Z4Code> C) -> ModTupRngElt  */
/* Signature: (<Tup> u, <Z2Z4Code> C) -> Tup                    */
/*                                                              */
/****************************************************************/
/****************************************************************/
/*                                                              */
/* Function name: CosetLeaders                                  */
/* Parameters: C                                                */
/* Function description: Given a Z2Z4-additive code C of type   */
/*   (alpha,beta;gamma;delta;kappa), with ambient space         */
/*   V=Z2^alpha x Z4^beta, return a set of coset leaders        */
/*   (vectors of minimal Lee weight in their cosets) for C in V */
/*   as an indexed set of vectors from V. The elements in       */
/*   V=Z2^alpha x Z4^beta are represented as elements in        */
/*   Z4^(alpha+beta) by changing the ones in the first alpha    */
/*   coordinates by twos. This function also returns a map from */
/*   the syndrome space of C onto the coset leaders (mapping a  */
/*   syndrome into its corresponding coset leader). Note that   */
/*   this function is only applicable when V and C are small.   */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/* Output parameters description:                               */
/*   - Indexed set of vectors from V                            */
/*   - Map from the syndrome space of C into the coset leaders  */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> SetIndx, Map                    */
/*                                                              */
/****************************************************************/
/****************************************************************/
/*                                                              */
/* Function name: CosetDecode                                   */
/* Parameters: C, u                                             */
/* Function description: Given a Z2Z4-additive code C of type   */
/*   (alpha, beta; gamma; delta; kappa), and a vector u from the*/
/*   ambient space V=Z2^alpha x Z4^beta or V2=Z2^(alpha+2*beta),*/
/*   attempt to decode u with respect to C. The elements in     */
/*   V=Z2^alpha x Z4^beta are represented as elements in        */
/*   Z4^(alpha+beta) by changing the ones in the first alpha    */
/*   coordinates by twos. If the decoding algorithm succeeds in */
/*   computing a vector u' in C as the decoded version of u in  */
/*   V, then the function returns true, u' and Phi(u'), where   */
/*   Phi is the Gray map. If the decoding algorithm does not    */
/*   succeed in decoding u, then the function returns false,    */
/*   the zero vector in V and the zero vector in V2.            */
/*   Instead of a vector u, we can also decode a sequence Q     */
/*   of vectors in V or V2.                                     */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/*   - u : Received vector to be decoded                        */
/* Output parameters description:                               */
/*   - Boolean, true if u is decoded and false otherwise        */
/*   - Decoded vector over Z4 of u, or the zero vector          */
/*   - Decoded binary vector of u, or the zero vector           */
/*                                                              */
/* case quaternary u                                            */
/* Signature: (<Z2Z4Code> C, <ModTupRngElt> u)                  */
/*               -> BoolElt, ModTupRngElt, ModTupFldElt         */
/* case binary u                                                */
/* Signature: (<Z2Z4Code> C, <ModTupFldElt> u)                  */
/*               -> BoolElt, ModTupRngElt, ModTupFldElt         */
/*                                                              */
/* case quaternary sequence Q                                   */
/* Signature: (<Z2Z4Code> C, <[ModTupRngElt]> Q)                */
/*               -> [BoolElt], [ModTupRngElt], [ModTupFldElt]   */
/* case binary sequence Q                                       */
/* Signature: (<Z2Z4Code> C, <[ModTupFldElt]> Q)                */
/*               -> [BoolElt], [ModTupRngElt], [ModTupFldElt]   */
/*                                                              */
/****************************************************************/
/****************************************************************/
/*                                                              */
/* Function name: SyndromeDecode                                */
/* Parameters: C, u                                             */
/* Function description: Given a Z2Z4-additive code C of type   */
/*   (alpha, beta; gamma; delta; kappa), and a vector u from the*/
/*   ambient space V=Z2^alpha x Z4^beta or V2=Z2^(alpha+2*beta),*/
/*   attempt to decode u with respect to C. The elements in     */
/*   V=Z2^alpha x Z4^beta are represented as elements in        */
/*   Z4^(alpha+beta) by changing the ones in the first alpha    */
/*   coordinates by twos. The decoding algorithm always succeeds*/
/*   in computing a vector u' in C as the decoded version of u  */
/*   in V, and the function returns true, u' and Phi(u'), where */
/*   Phi is the Gray map. Although the function never returns   */
/*   false, the first output parameter true is given to be      */
/*   consistent with the other decoding functions.              */
/*   Instead of a vector u, we can also decode a sequence Q     */
/*   of vectors in V or V2.                                     */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/*   - u : Received vector to be decoded                        */
/* Output parameters description:                               */
/*   - Boolean, true if u is decoded and false otherwise        */
/*   - Decoded vector over Z4 of u, or the zero vector          */
/*   - Decoded binary vector of u, or the zero vector           */
/*                                                              */
/* case quaternary u                                            */
/* Signature: (<Z2Z4Code> C, <ModTupRngElt> u)                  */
/*               -> BoolElt, ModTupRngElt, ModTupFldElt         */
/* case binary u                                                */
/* Signature: (<Z2Z4Code> C, <ModTupFldElt> u)                  */
/*               -> BoolElt, ModTupRngElt, ModTupFldElt         */
/*                                                              */
/* case quaternary sequence Q                                   */
/* Signature: (<Z2Z4Code> C, <[ModTupRngElt]> Q)                */
/*               -> [BoolElt], [ModTupRngElt], [ModTupFldElt]   */
/* case binary sequence Q                                       */
/* Signature: (<Z2Z4Code> C, <[ModTupFldElt]> Q)                */
/*               -> [BoolElt], [ModTupRngElt], [ModTupFldElt]   */
/*                                                              */
/****************************************************************/
print "INFORMATION SPACE AND INFORMATION SETS";
print "SYNDROME SPACE AND COSET LEADERS";
print "DECODE";

print "test 1: Trivial Z2Z4-additive zero code
               alpha = 2, beta = 4, gamma = 0, delta = 0, kappa = 0, length = 6, #C = 1";
C := Z2Z4AdditiveZeroCode(2, 4);
alpha := C`Alpha;
beta := Length(C)-alpha;
gamma := Z2Z4Type(C)[3];
delta := Z2Z4Type(C)[4];
kappa := Z2Z4Type(C)[5];
grayMap := GrayMap(Z2Z4AdditiveUniverseCode(alpha, beta));
R := RSpace(Z4, alpha+beta);
V := VectorSpace(GF(2), alpha+2*beta);
u := C`Code!0;                       // u in C
ubin := grayMap(u);                  // ubin in Cbin
v := R![2,0,1,3,3,3];                // v not in C
vbin := grayMap(v);                  // vbin not in Cbin

// Test InformationSpace function    C cannot be the zero code

// Test InformationSet function      C cannot be the zero code
 
// Test IsInformationSet function    C cannot be the zero code

// Test SyndromeSpace function
OutputRs, OutputVs := SyndromeSpace(C);
assert OutputRs eq RSpace(LinearCode(DiagonalMatrix(Z4,[2^^(alpha+gamma-2*kappa)] 
                                             cat [1^^(beta-gamma-delta+kappa)])));
assert OutputVs eq VectorSpace(GF(2), alpha+2*beta-gamma-2*delta);

// Test Syndrome function
OutputSyndrome_u := Syndrome(u, C);
OutputSyndrome_ubin := Syndrome(ubin, C);
OutputSyndrome_v := Syndrome(v, C);
OutputSyndrome_vbin := Syndrome(vbin, C);
assert OutputSyndrome_u eq OutputRs!0;
assert OutputSyndrome_ubin eq OutputRs!0;
assert OutputSyndrome_v eq OutputRs![0,2,3,3,3,1];       
assert OutputSyndrome_vbin eq OutputRs![0,2,3,3,3,1];
assert OutputSyndrome_u eq Z2Z4Mult(u, Transpose(MinRowsParityCheckMatrix(C)), C`Alpha);
assert OutputSyndrome_v eq Z2Z4Mult(v, Transpose(MinRowsParityCheckMatrix(C)), C`Alpha);

// Test CosetLeaders function
OutputCosetLeadersL, OutputCosetLeadersMap := CosetLeaders(C);     
assert #OutputCosetLeadersL eq #OutputRs;
assert OutputCosetLeadersL eq SetToIndexedSet(Set(OutputRs));
assert OutputCosetLeadersMap(OutputRs!0) eq u;
assert OutputCosetLeadersMap(OutputRs![0,2,3,3,3,1]) eq v;

// Test CosetDecode function
OutputIsDecoded_u, OutputDecoded_u, OutputDecoded_ubin := CosetDecode(C, u);
assert OutputIsDecoded_u eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;
OutputIsDecoded_ubin, OutputDecoded_u, OutputDecoded_ubin := CosetDecode(C, ubin);
assert OutputIsDecoded_ubin eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;
OutputIsDecoded_v, OutputDecoded_v, OutputDecoded_vbin := CosetDecode(C, v);
assert OutputIsDecoded_v eq true;
assert OutputDecoded_v eq C`Code!0;
assert OutputDecoded_vbin eq V!0;
OutputIsDecoded_vbin, OutputDecode_v, OutputDecoded_vbin := CosetDecode(C, vbin);
assert OutputIsDecoded_vbin eq true;
assert OutputDecoded_v eq C`Code!0;
assert OutputDecoded_vbin eq V!0;

// Test SyndromeDecode function
OutputIsDecoded_u, OutputDecoded_u, OutputDecoded_ubin := SyndromeDecode(C, u);
assert OutputIsDecoded_u eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;
OutputIsDecoded_ubin, OutputDecoded_u, OutputDecoded_ubin := SyndromeDecode(C, ubin);
assert OutputIsDecoded_ubin eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;
OutputIsDecoded_v, OutputDecoded_v, OutputDecoded_vbin := SyndromeDecode(C, v);
assert OutputIsDecoded_v eq true;
assert OutputDecoded_v eq C`Code!0;
assert OutputDecoded_vbin eq V!0;
OutputIsDecoded_vbin, OutputDecode_v, OutputDecoded_vbin := SyndromeDecode(C, vbin);
assert OutputIsDecoded_vbin eq true;
assert OutputDecoded_v eq C`Code!0;
assert OutputDecoded_vbin eq V!0;

/****************************************************************/
print "test 2: Trivial Z2Z4-additive universe code
               alpha = 2, beta = 4, gamma = 2, delta = 4, kappa = 2, length = 6, #C = 1024";
C := Z2Z4AdditiveUniverseCode(2, 4);
alpha := C`Alpha;
beta :=  Length(C)-alpha;
gamma := Z2Z4Type(C)[3];
delta := Z2Z4Type(C)[4];
kappa := Z2Z4Type(C)[5];
grayMap :=  GrayMap(Z2Z4AdditiveUniverseCode(alpha, beta));
R := RSpace(Z4, alpha+beta);
V := VectorSpace(GF(2), alpha+2*beta);
u := R![2,0,1,3,3,3];                // u in C
ubin := grayMap(u);                  // ubin in Cbin
Is := [1..alpha+beta];               // given information set for C
Isbin := [1..alpha+2*beta];          // given information set for Cbin
J := [1,2,3];                        // not an information set
Jbin := [1,2,3,4,5];                 // not an information set

// Test  InformationSpace function
OutputR, OutputV, Outputf, Outputfbin :=  InformationSpace(C);
assert OutputR eq RSpace(LinearCode(DiagonalMatrix(Z4,[2^^gamma] cat [1^^delta])));
assert OutputV eq VectorSpace(GF(2), gamma+2*delta);

// Test  InformationSet function
OutputInfoSet_C, OutputInfoSet_Cbin :=  InformationSet(C);
assert OutputInfoSet_C eq Is;
assert OutputInfoSet_Cbin eq Isbin;

// Test IsInformationSet function
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, Is);
assert OutputIsInfoSet_C eq true;
assert OutputIsInfoSet_Cbin eq false;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, Isbin);
assert OutputIsInfoSet_C eq false;
assert OutputIsInfoSet_Cbin eq true;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, J);
assert OutputIsInfoSet_C eq false;
assert OutputIsInfoSet_Cbin eq false;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, Jbin);
assert OutputIsInfoSet_C eq false;
assert OutputIsInfoSet_Cbin eq false;

// Test  SyndromeSpace function        C cannot be the universe code

// Test  Syndrome function
OutputSyndrome_u :=  Syndrome(u, C);
OutputSyndrome_ubin :=  Syndrome(ubin, C);
assert OutputSyndrome_u eq Matrix(Z4,1,0,[]);
assert OutputSyndrome_ubin eq Matrix(Z4,1,0,[]);
assert OutputSyndrome_u eq Z2Z4Mult(u, Transpose( MinRowsParityCheckMatrix(C)), C`Alpha);

// Test  CosetLeaders function
OutputCosetLeadersL, OutputCosetLeadersMap :=  CosetLeaders(C);
assert OutputCosetLeadersL eq {@ C`Code!0 @};
assert OutputCosetLeadersMap(Matrix(Z4,1,0,[])) eq C`Code!0;

// Test  CosetDecode function
OutputIsDecoded_u, OutputDecoded_u, OutputDecoded_ubin :=  CosetDecode(C, u);
assert OutputIsDecoded_u eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;
OutputIsDecoded_ubin, OutputDecoded_u, OutputDecoded_ubin :=  CosetDecode(C, ubin);
assert OutputIsDecoded_ubin eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;

// Test  SyndromeDecode function
OutputIsDecoded_u, OutputDecoded_u, OutputDecoded_ubin :=  SyndromeDecode(C, u);
assert OutputIsDecoded_u eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;
OutputIsDecoded_ubin, OutputDecoded_u, OutputDecoded_ubin :=  SyndromeDecode(C, ubin);
assert OutputIsDecoded_ubin eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;

/************************************************************/
print "test 3: Repetition Z2Z4-additive code
               alpha = 4, beta = 3, gamma = 1, delta = 0, kappa = 1, length = 7, #C = 2";
C := Z2Z4AdditiveRepetitionCode(4, 3);
alpha := C`Alpha;
beta :=  Length(C)-alpha;
gamma := Z2Z4Type(C)[3];
delta := Z2Z4Type(C)[4];
kappa := Z2Z4Type(C)[5];
grayMap :=  GrayMap(Z2Z4AdditiveUniverseCode(alpha, beta));
R := RSpace(Z4, alpha+beta);
V := VectorSpace(GF(2), alpha+2*beta);
u := R![2,2,2,2,2,2,2];         // u in C
ubin := grayMap(u);             // ubin in Cbin
v := R![2,0,2,2,3,1,2];         // v not in C, adding 4 binary errors
vbin := grayMap(v);             // vbin not in Cbin, adding 4 binary errors
Is := [1];                      // given information set for C
Isbin := [1];                   // given information set for Cbin
I := [2];                       // another information set
Ibin := [6];                    // another information set

// Test  InformationSpace function
OutputR, OutputV, Outputf, Outputfbin :=  InformationSpace(C);
assert OutputR eq RSpace(LinearCode(DiagonalMatrix(Z4,[2^^gamma] cat [1^^delta])));
assert OutputV eq VectorSpace(GF(2), gamma+2*delta);

// Test  InformationSet function
OutputInfoSet_C, OutputInfoSet_Cbin :=  InformationSet(C);
assert OutputInfoSet_C eq Is;
assert OutputInfoSet_Cbin eq Isbin;

// Test IsInformationSet function
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, Is);
assert OutputIsInfoSet_C eq true;
assert OutputIsInfoSet_Cbin eq true;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, Isbin);
assert OutputIsInfoSet_C eq true;
assert OutputIsInfoSet_Cbin eq true;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, I);
assert OutputIsInfoSet_C eq true;
assert OutputIsInfoSet_Cbin eq true;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, Ibin);
assert OutputIsInfoSet_C eq true;
assert OutputIsInfoSet_Cbin eq true;

// Test  SyndromeSpace function
OutputRs, OutputVs :=  SyndromeSpace(C);
assert OutputRs eq RSpace(LinearCode(DiagonalMatrix(Z4,[2^^(alpha+gamma-2*kappa)] 
                                             cat [1^^(beta-gamma-delta+kappa)])));
assert OutputVs eq VectorSpace(GF(2), alpha+2*beta-gamma-2*delta);

// Test  Syndrome function
OutputSyndrome_u :=  Syndrome(u, C);
OutputSyndrome_ubin :=  Syndrome(ubin, C);
OutputSyndrome_v :=  Syndrome(v, C);
OutputSyndrome_vbin :=  Syndrome(vbin, C);
assert OutputSyndrome_u eq OutputRs!0;
assert OutputSyndrome_ubin eq OutputRs!0;
assert OutputSyndrome_v eq OutputRs![2,0,0,1,3,0];
assert OutputSyndrome_vbin eq OutputRs![2,0,0,1,3,0];
assert OutputSyndrome_u eq Z2Z4Mult(u, Transpose( MinRowsParityCheckMatrix(C)), C`Alpha);
assert OutputSyndrome_v eq Z2Z4Mult(v, Transpose( MinRowsParityCheckMatrix(C)), C`Alpha);

// Test  CosetLeaders function
OutputCosetLeadersL, OutputCosetLeadersMap :=  CosetLeaders(C);
assert #OutputCosetLeadersL eq #OutputRs;
assert OutputCosetLeadersMap(OutputRs!0) eq C`Code!0;
assert OutputCosetLeadersMap(OutputRs![2,0,0,1,3,0]) eq v-u;

// Test  CosetDecode function
OutputIsDecoded_u, OutputDecoded_u, OutputDecoded_ubin :=  CosetDecode(C, u);
assert OutputIsDecoded_u eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;
OutputIsDecoded_ubin, OutputDecoded_u, OutputDecoded_ubin :=  CosetDecode(C, ubin);
assert OutputIsDecoded_ubin eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;
OutputIsDecoded_v, OutputDecoded_v, OutputDecoded_vbin :=  CosetDecode(C, v);
assert OutputIsDecoded_v eq true;
assert OutputDecoded_v eq u;
assert OutputDecoded_vbin eq ubin;
OutputIsDecoded_vbin, OutputDecode_v, OutputDecoded_vbin :=  CosetDecode(C, vbin);
assert OutputIsDecoded_vbin eq true;
assert OutputDecoded_v eq u;
assert OutputDecoded_vbin eq ubin;

// Test  SyndromeDecode function
OutputIsDecoded_u, OutputDecoded_u, OutputDecoded_ubin :=  SyndromeDecode(C, u);
assert OutputIsDecoded_u eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;
OutputIsDecoded_ubin, OutputDecoded_u, OutputDecoded_ubin :=  SyndromeDecode(C, ubin);
assert OutputIsDecoded_ubin eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;
OutputIsDecoded_v, OutputDecoded_v, OutputDecoded_vbin :=  SyndromeDecode(C, v);
assert OutputIsDecoded_v eq true;
assert OutputDecoded_v eq u;
assert OutputDecoded_vbin eq ubin;
OutputIsDecoded_vbin, OutputDecode_v, OutputDecoded_vbin :=  SyndromeDecode(C, vbin);
assert OutputIsDecoded_vbin eq true;
assert OutputDecoded_v eq u;
assert OutputDecoded_vbin eq ubin;

/************************************************************/
print "test 4: A Z2Z4-additive code with alpha=0
               alpha = 0, beta = 10, gamma = 5, delta = 2, kappa = 0, length = 10, #C = 512";
R := RSpace(Z4, 10);
L := [R![0,2,0,0,0,0,0,0,0,2],
      R![0,0,0,0,0,2,0,0,0,2],
      R![0,0,0,2,0,0,0,0,0,2],
      R![0,0,0,0,2,0,0,0,0,2],
      R![0,0,0,0,0,0,0,2,2,2],
      R![1,0,1,1,1,1,0,0,2,3],
      R![0,0,0,0,0,0,1,0,1,1]];
C := Z2Z4AdditiveCode(L, 0);
alpha := C`Alpha;
beta :=  Length(C)-alpha;
gamma := Z2Z4Type(C)[3];
delta := Z2Z4Type(C)[4];
kappa := Z2Z4Type(C)[5];
grayMap :=  GrayMap(Z2Z4AdditiveUniverseCode(alpha, beta));
R := RSpace(Z4, alpha+beta);
V := VectorSpace(GF(2), alpha+2*beta);
u := R![1,0,1,1,1,1,0,0,2,3];         // u in C
ubin := grayMap(u);                   // ubin in Cbin
v := R![1,0,2,1,1,1,0,0,2,3];         // v not in C, adding 1 binary error
vbin := grayMap(v);                   // vbin not in Cbin, adding 1 binary error
Is := [4,3,5,6,8,9,10];               // given information set for C
Isbin := [7,5,9,11,15,17,18,19,20];   // given information set for Cbin
I := [2,3,4,6,8,9,10];                // another information set
Ibin := [3,5,7,11,15,17,18,19,20];    // another information set
J := [1..7];                          // not an information set
Jbin := [1..9];                       // not an information set

// Test  InformationSpace function
OutputR, OutputV, Outputf, Outputfbin :=  InformationSpace(C);
assert OutputR eq RSpace(LinearCode(DiagonalMatrix(Z4,[2^^gamma] cat [1^^delta])));
assert OutputV eq VectorSpace(GF(2), gamma+2*delta);

// Test  InformationSet function
OutputInfoSet_C, OutputInfoSet_Cbin :=  InformationSet(C);
assert OutputInfoSet_C eq Is;
assert OutputInfoSet_Cbin eq Isbin;

// Test IsInformationSet function
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, Is);
assert OutputIsInfoSet_C eq true;
assert OutputIsInfoSet_Cbin eq false;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, Isbin);
assert OutputIsInfoSet_C eq false;
assert OutputIsInfoSet_Cbin eq true;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, I);
assert OutputIsInfoSet_C eq true;
assert OutputIsInfoSet_Cbin eq false;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, Ibin);
assert OutputIsInfoSet_C eq false;
assert OutputIsInfoSet_Cbin eq true;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, J);
assert OutputIsInfoSet_C eq false;
assert OutputIsInfoSet_Cbin eq false;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, Jbin);
assert OutputIsInfoSet_C eq false;
assert OutputIsInfoSet_Cbin eq false;

// Test  SyndromeSpace function
OutputRs, OutputVs :=  SyndromeSpace(C);
assert OutputRs eq RSpace(LinearCode(DiagonalMatrix(Z4,[2^^(alpha+gamma-2*kappa)] 
                                             cat [1^^(beta-gamma-delta+kappa)])));
assert OutputVs eq VectorSpace(GF(2), alpha+2*beta-gamma-2*delta);

// Test  Syndrome function
OutputSyndrome_u :=  Syndrome(u, C);
OutputSyndrome_ubin :=  Syndrome(ubin, C);
OutputSyndrome_v :=  Syndrome(v, C);
OutputSyndrome_vbin :=  Syndrome(vbin, C);
assert OutputSyndrome_u eq OutputRs!0;
assert OutputSyndrome_ubin eq OutputRs!0;
assert OutputSyndrome_v eq OutputRs![0,0,0,0,2,0,1,0];
assert OutputSyndrome_vbin eq OutputRs![0,0,0,0,2,0,1,0];
assert OutputSyndrome_u eq Z2Z4Mult(u, Transpose( MinRowsParityCheckMatrix(C)), C`Alpha);
assert OutputSyndrome_v eq Z2Z4Mult(v, Transpose( MinRowsParityCheckMatrix(C)), C`Alpha);

// Test  CosetLeaders function
OutputCosetLeadersL, OutputCosetLeadersMap :=  CosetLeaders(C);
assert #OutputCosetLeadersL eq #OutputRs;
assert OutputCosetLeadersMap(OutputRs!0) eq C`Code!0;
assert OutputCosetLeadersMap(OutputRs![0,0,0,0,2,0,1,0]) eq v-u;

// Test  CosetDecode function
OutputIsDecoded_u, OutputDecoded_u, OutputDecoded_ubin :=  CosetDecode(C, u);
assert OutputIsDecoded_u eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;
OutputIsDecoded_ubin, OutputDecoded_u, OutputDecoded_ubin :=  CosetDecode(C, ubin);
assert OutputIsDecoded_ubin eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;
OutputIsDecoded_v, OutputDecoded_v, OutputDecoded_vbin :=  CosetDecode(C, v);
assert OutputIsDecoded_v eq true;
assert OutputDecoded_v eq u;
assert OutputDecoded_vbin eq ubin;
OutputIsDecoded_vbin, OutputDecode_v, OutputDecoded_vbin :=  CosetDecode(C, vbin);
assert OutputIsDecoded_vbin eq true;
assert OutputDecoded_v eq u;
assert OutputDecoded_vbin eq ubin;

// Test  SyndromeDecode function
OutputIsDecoded_u, OutputDecoded_u, OutputDecoded_ubin :=  SyndromeDecode(C, u);
assert OutputIsDecoded_u eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;
OutputIsDecoded_ubin, OutputDecoded_u, OutputDecoded_ubin :=  SyndromeDecode(C, ubin);
assert OutputIsDecoded_ubin eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;
OutputIsDecoded_v, OutputDecoded_v, OutputDecoded_vbin :=  SyndromeDecode(C, v);
assert OutputIsDecoded_v eq true;
assert OutputDecoded_v eq u;
assert OutputDecoded_vbin eq ubin;
OutputIsDecoded_vbin, OutputDecode_v, OutputDecoded_vbin :=  SyndromeDecode(C, vbin);
assert OutputIsDecoded_vbin eq true;
assert OutputDecoded_v eq u;
assert OutputDecoded_vbin eq ubin;

/************************************************************/
print "test 5: A Z2Z4-additive code with beta=0, a Hamming code
               alpha = 15, beta = 0, gamma = 5, delta = 2, kappa = 0, length = 15, #C = 2048";
C := Z2Z4AdditiveCode(HammingCode(GF(2), 4));
alpha := C`Alpha;
beta :=  Length(C)-alpha;
gamma := Z2Z4Type(C)[3];
delta := Z2Z4Type(C)[4];
kappa := Z2Z4Type(C)[5];
grayMap :=  GrayMap(Z2Z4AdditiveUniverseCode(alpha, beta));
R := RSpace(Z4, alpha+beta);
V := VectorSpace(GF(2), alpha+2*beta);
u := R![2,2,0,2,0,0,0,2,2,2,2,2,2,0,2]; // u in C
ubin := grayMap(u);                     // ubin in Cbin
v := R![2,2,0,2,0,0,0,2,2,2,2,2,2,0,0]; // u not in C, adding 1 binary error
vbin := grayMap(v);                     // u not in Cbin, adding 1 binary error
w := R![2,0,2,2,0,0,0,2,2,2,2,2,2,0,2]; // u not in C, adding 2 binary errors
wbin := grayMap(w);                     // u not in Cbin, adding 2 binary errors
Is := [1..11];                          // given information set for C
Isbin := [1..11];                       // given information set for Cbin
I := [5..15];                           // another information set
Ibin := [4..14];                        // another information set
J := [1..12];                           // not an information set
Jbin := [6..15];                        // not an information set

// Test  InformationSpace function
OutputR, OutputV, Outputf, Outputfbin :=  InformationSpace(C);
assert OutputR eq RSpace(LinearCode(DiagonalMatrix(Z4,[2^^gamma] cat [1^^delta])));
assert OutputV eq VectorSpace(GF(2), gamma+2*delta);

// Test  InformationSet function
OutputInfoSet_C, OutputInfoSet_Cbin :=  InformationSet(C);
assert OutputInfoSet_C eq Is;
assert OutputInfoSet_Cbin eq Isbin;

// Test IsInformationSet function
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, Is);
assert OutputIsInfoSet_C eq true;
assert OutputIsInfoSet_Cbin eq true;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, Isbin);
assert OutputIsInfoSet_C eq true;
assert OutputIsInfoSet_Cbin eq true;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, I);
assert OutputIsInfoSet_C eq true;
assert OutputIsInfoSet_Cbin eq true;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, Ibin);
assert OutputIsInfoSet_C eq true;
assert OutputIsInfoSet_Cbin eq true;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, J);
assert OutputIsInfoSet_C eq false;
assert OutputIsInfoSet_Cbin eq false;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, Jbin);
assert OutputIsInfoSet_C eq false;
assert OutputIsInfoSet_Cbin eq false;

// Test  SyndromeSpace function
OutputRs, OutputVs :=  SyndromeSpace(C);
assert OutputRs eq RSpace(LinearCode(DiagonalMatrix(Z4,[2^^(alpha+gamma-2*kappa)] 
                                             cat [1^^(beta-gamma-delta+kappa)])));
assert OutputVs eq VectorSpace(GF(2), alpha+2*beta-gamma-2*delta);

// Test  Syndrome function
OutputSyndrome_u :=  Syndrome(u, C);
OutputSyndrome_ubin :=  Syndrome(ubin, C);
OutputSyndrome_v :=  Syndrome(v, C);
OutputSyndrome_vbin :=  Syndrome(vbin, C);
OutputSyndrome_w :=  Syndrome(w, C);
OutputSyndrome_wbin :=  Syndrome(wbin, C);
assert OutputSyndrome_u eq OutputRs!0;
assert OutputSyndrome_ubin eq OutputRs!0;
assert OutputSyndrome_v eq OutputRs![0,0,2,2];
assert OutputSyndrome_vbin eq OutputRs![0,0,2,2];
assert OutputSyndrome_w eq OutputRs![2,2,0,0];
assert OutputSyndrome_wbin eq OutputRs![2,2,0,0];
assert OutputSyndrome_u eq Z2Z4Mult(u, Transpose( MinRowsParityCheckMatrix(C)), C`Alpha);
assert OutputSyndrome_v eq Z2Z4Mult(v, Transpose( MinRowsParityCheckMatrix(C)), C`Alpha);
assert OutputSyndrome_w eq Z2Z4Mult(w, Transpose( MinRowsParityCheckMatrix(C)), C`Alpha);

// Test  CosetLeaders function
OutputCosetLeadersL, OutputCosetLeadersMap :=  CosetLeaders(C);
assert #OutputCosetLeadersL eq #OutputRs;
assert OutputCosetLeadersMap(OutputRs!0) eq C`Code!0;
assert OutputCosetLeadersMap(OutputRs![0,0,2,2]) eq v-u;
assert OutputCosetLeadersMap(OutputRs![2,2,0,0]) eq R![0,0,0,0,0,2,0,0,0,0,0,0,0,0,0];

// Test  CosetDecode function
OutputIsDecoded_u, OutputDecoded_u, OutputDecoded_ubin :=  CosetDecode(C, u);
assert OutputIsDecoded_u eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;
OutputIsDecoded_ubin, OutputDecoded_u, OutputDecoded_ubin :=  CosetDecode(C, ubin);
assert OutputIsDecoded_ubin eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;
OutputIsDecoded_v, OutputDecoded_v, OutputDecoded_vbin :=  CosetDecode(C, v);
assert OutputIsDecoded_v eq true;
assert OutputDecoded_v eq u;
assert OutputDecoded_vbin eq ubin;
OutputIsDecoded_vbin, OutputDecode_v, OutputDecoded_vbin :=  CosetDecode(C, vbin);
assert OutputIsDecoded_vbin eq true;
assert OutputDecoded_v eq u;
assert OutputDecoded_vbin eq ubin;
OutputIsDecoded_w, OutputDecoded_w, OutputDecoded_wbin :=  CosetDecode(C, w);
assert OutputIsDecoded_w eq true;
assert OutputDecoded_w eq R![2,0,2,2,0,2,0,2,2,2,2,2,2,0,2];
assert OutputDecoded_wbin eq V![1,0,1,1,0,1,0,1,1,1,1,1,1,0,1];
OutputIsDecoded_wbin, OutputDecode_w, OutputDecoded_wbin :=  CosetDecode(C, wbin);
assert OutputIsDecoded_wbin eq true;
assert OutputDecoded_w eq R![2,0,2,2,0,2,0,2,2,2,2,2,2,0,2];
assert OutputDecoded_wbin eq V![1,0,1,1,0,1,0,1,1,1,1,1,1,0,1];

// Test  SyndromeDecode function
OutputIsDecoded_u, OutputDecoded_u, OutputDecoded_ubin :=  SyndromeDecode(C, u);
assert OutputIsDecoded_u eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;
OutputIsDecoded_ubin, OutputDecoded_u, OutputDecoded_ubin :=  SyndromeDecode(C, ubin);
assert OutputIsDecoded_ubin eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;
OutputIsDecoded_v, OutputDecoded_v, OutputDecoded_vbin :=  SyndromeDecode(C, v);
assert OutputIsDecoded_v eq true;
assert OutputDecoded_v eq u;
assert OutputDecoded_vbin eq ubin;
OutputIsDecoded_vbin, OutputDecode_v, OutputDecoded_vbin :=  SyndromeDecode(C, vbin);
assert OutputIsDecoded_vbin eq true;
assert OutputDecoded_v eq u;
assert OutputDecoded_vbin eq ubin;
OutputIsDecoded_w, OutputDecoded_w, OutputDecoded_wbin :=  SyndromeDecode(C, w);
assert OutputIsDecoded_w eq true;
assert OutputDecoded_w eq R![2,0,2,2,0,2,0,2,2,2,2,2,2,0,2];
assert OutputDecoded_wbin eq V![1,0,1,1,0,1,0,1,1,1,1,1,1,0,1];
OutputIsDecoded_wbin, OutputDecode_w, OutputDecoded_wbin :=  SyndromeDecode(C, wbin);
assert OutputIsDecoded_wbin eq true;
assert OutputDecoded_w eq R![2,0,2,2,0,2,0,2,2,2,2,2,2,0,2];
assert OutputDecoded_wbin eq V![1,0,1,1,0,1,0,1,1,1,1,1,1,0,1];

/************************************************************/
print "test 6: A Z2Z4-additive code with alpha=0 and gamma=0, a Preparata code
               alpha = 0, beta = 8, gamma = 0, delta = 4, kappa = 0, length = 8, #C = 256";
C := Z2Z4AdditiveCode(PreparataCode(3));
alpha := C`Alpha;
beta :=  Length(C)-alpha;
gamma := Z2Z4Type(C)[3];
delta := Z2Z4Type(C)[4];
kappa := Z2Z4Type(C)[5];
grayMap :=  GrayMap(Z2Z4AdditiveUniverseCode(alpha, beta));
R := RSpace(Z4, alpha+beta);
V := VectorSpace(GF(2), alpha+2*beta);
d := 6;   // error correcting capability t=2
          // dK=8, mimimum weight of the kernel
u := R![1,1,3,3,1,3,1,3];               // u in C
ubin := grayMap(u);                     // ubin in Cbin
v := R![1,1,3,3,1,3,1,1];               // v not in C, with 2 errors
vbin := grayMap(v);                     // u not in Cbin, with 2 errors, 2 <=t
w := R![3,3,1,1,1,3,1,3];               // w not in C, with 8 errors
wbin := grayMap(w);                     // u not in Cbin, with 8 errors, 8 >= dK-1=7
Is := [5..8];                           // given information set for C
Isbin := [9..16];                       // given information set for Cbin
I := [1..4];                            // another information set
Ibin := [1..8];                         // another information set
J := [1..12];                           // not an information set
Jbin := [6..15];                        // not an information set

// Test  InformationSpace function
OutputR, OutputV, Outputf, Outputfbin :=  InformationSpace(C);
assert OutputR eq RSpace(LinearCode(DiagonalMatrix(Z4,[2^^gamma] cat [1^^delta])));
assert OutputV eq VectorSpace(GF(2), gamma+2*delta);

// Test  InformationSet function
OutputInfoSet_C, OutputInfoSet_Cbin :=  InformationSet(C);
assert OutputInfoSet_C eq Is;
assert OutputInfoSet_Cbin eq Isbin;

// Test IsInformationSet function
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, Is);
assert OutputIsInfoSet_C eq true;
assert OutputIsInfoSet_Cbin eq false;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, Isbin);
assert OutputIsInfoSet_C eq false;
assert OutputIsInfoSet_Cbin eq true;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, I);
assert OutputIsInfoSet_C eq true;
assert OutputIsInfoSet_Cbin eq false;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, Ibin);
assert OutputIsInfoSet_C eq false;
assert OutputIsInfoSet_Cbin eq true;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, J);
assert OutputIsInfoSet_C eq false;
assert OutputIsInfoSet_Cbin eq false;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, Jbin);
assert OutputIsInfoSet_C eq false;
assert OutputIsInfoSet_Cbin eq false;

// Test  SyndromeSpace function
OutputRs, OutputVs :=  SyndromeSpace(C);
assert OutputRs eq RSpace(LinearCode(DiagonalMatrix(Z4,[2^^(alpha+gamma-2*kappa)] 
                                             cat [1^^(beta-gamma-delta+kappa)])));
assert OutputVs eq VectorSpace(GF(2), alpha+2*beta-gamma-2*delta);

// Test  Syndrome function
OutputSyndrome_u :=  Syndrome(u, C);
OutputSyndrome_ubin :=  Syndrome(ubin, C);
OutputSyndrome_v :=  Syndrome(v, C);
OutputSyndrome_vbin :=  Syndrome(vbin, C);
OutputSyndrome_w :=  Syndrome(w, C);
OutputSyndrome_wbin :=  Syndrome(wbin, C);
assert OutputSyndrome_u eq OutputRs!0;
assert OutputSyndrome_ubin eq OutputRs!0;
assert OutputSyndrome_v eq OutputRs![0,2,2,2];
assert OutputSyndrome_vbin eq OutputRs![0,2,2,2];
assert OutputSyndrome_w eq OutputRs![2,2,2,2];
assert OutputSyndrome_wbin eq OutputRs![2,2,2,2];
assert OutputSyndrome_u eq Z2Z4Mult(u, Transpose( MinRowsParityCheckMatrix(C)), C`Alpha);
assert OutputSyndrome_v eq Z2Z4Mult(v, Transpose( MinRowsParityCheckMatrix(C)), C`Alpha);
assert OutputSyndrome_w eq Z2Z4Mult(w, Transpose( MinRowsParityCheckMatrix(C)), C`Alpha);

// Test  CosetLeaders function
OutputCosetLeadersL, OutputCosetLeadersMap :=  CosetLeaders(C);
assert #OutputCosetLeadersL eq #OutputRs;
assert OutputCosetLeadersMap(OutputRs!0) eq C`Code!0;
assert OutputCosetLeadersMap(OutputRs![0,2,2,2]) eq v-u;
assert OutputCosetLeadersMap(OutputRs![2,2,2,2]) eq R![0,1,3,3,0,0,1,0];

// Test  CosetDecode function
OutputIsDecoded_u, OutputDecoded_u, OutputDecoded_ubin :=  CosetDecode(C, u);
assert OutputIsDecoded_u eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;
OutputIsDecoded_ubin, OutputDecoded_u, OutputDecoded_ubin :=  CosetDecode(C, ubin);
assert OutputIsDecoded_ubin eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;
OutputIsDecoded_v, OutputDecoded_v, OutputDecoded_vbin :=  CosetDecode(C, v);
assert OutputIsDecoded_v eq true;
assert OutputDecoded_v eq u;
assert OutputDecoded_vbin eq ubin;
OutputIsDecoded_vbin, OutputDecode_v, OutputDecoded_vbin :=  CosetDecode(C, vbin);
assert OutputIsDecoded_vbin eq true;
assert OutputDecoded_v eq u;
assert OutputDecoded_vbin eq ubin;
OutputIsDecoded_w, OutputDecoded_w, OutputDecoded_wbin :=  CosetDecode(C, w);
assert OutputIsDecoded_w eq true;
assert OutputDecoded_w eq R![3,3,3,1,1,3,1,1];
assert OutputDecoded_wbin eq V![1,0,1,0,1,0,0,1,0,1,1,0,0,1,0,1];
OutputIsDecoded_wbin, OutputDecode_w, OutputDecoded_wbin :=  CosetDecode(C, wbin);
assert OutputIsDecoded_wbin eq true;
assert OutputDecoded_w eq R![3,3,3,1,1,3,1,1];
assert OutputDecoded_wbin eq V![1,0,1,0,1,0,0,1,0,1,1,0,0,1,0,1];

// Test  SyndromeDecode function
OutputIsDecoded_u, OutputDecoded_u, OutputDecoded_ubin :=  SyndromeDecode(C, u);
assert OutputIsDecoded_u eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;
OutputIsDecoded_ubin, OutputDecoded_u, OutputDecoded_ubin :=  SyndromeDecode(C, ubin);
assert OutputIsDecoded_ubin eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;
OutputIsDecoded_v, OutputDecoded_v, OutputDecoded_vbin :=  SyndromeDecode(C, v);
assert OutputIsDecoded_v eq true;
assert OutputDecoded_v eq u;
assert OutputDecoded_vbin eq ubin;
OutputIsDecoded_vbin, OutputDecode_v, OutputDecoded_vbin :=  SyndromeDecode(C, vbin);
assert OutputIsDecoded_vbin eq true;
assert OutputDecoded_v eq u;
assert OutputDecoded_vbin eq ubin;
OutputIsDecoded_w, OutputDecoded_w, OutputDecoded_wbin :=  SyndromeDecode(C, w);
assert OutputIsDecoded_w eq true;
assert OutputDecoded_w eq R![3,2,2,2,1,3,0,3];
assert OutputDecoded_wbin eq V![1,0,1,1,1,1,1,1,0,1,1,0,0,0,1,0];
OutputIsDecoded_wbin, OutputDecode_w, OutputDecoded_wbin :=  SyndromeDecode(C, wbin);
assert OutputIsDecoded_wbin eq true;
assert OutputDecoded_w eq R![3,2,2,2,1,3,0,3];
assert OutputDecoded_wbin eq V![1,0,1,1,1,1,1,1,0,1,1,0,0,0,1,0];

/****************************************************************/  
print "test 7: A Z2Z4-additive code with minimum distance 10
               alpha = 10, beta = 20, gamma = 6, delta = 2, kappa = 1, length = 30, #C = 1024";
C := Z2Z4AdditiveCode(Matrix(Integers(4),
        [[2,0,2,0,0,0,2,0,2,2,0,0,0,1,0,0,0,3,2,1,2,0,1,2,0,1,2,1,0,3],
         [0,2,2,0,0,2,0,0,2,2,1,1,1,1,0,1,0,3,2,0,0,2,2,0,1,0,3,2,1,1],
         [0,0,0,0,2,0,2,2,2,0,0,0,0,1,0,0,0,1,2,3,0,0,1,0,2,3,0,3,2,3],
         [0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,2,2,0,0,2,0,2,2,0,0,2,2],
         [0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2,0,2,2,2,2,2,0,0,2,2,2,0],
         [0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,2,0,2,0,0,0,2,2,0,2,2,2,2],
         [0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,0,2,0,0,2,0,0,2,0,2,0,2],
         [0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,2,0,2,0,0,2,2,0,2,0,0,2,2],
         [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,2,0,2,2,2,0,2,0,2,2,0,0],
         [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,0,0,0,0,0,0,2,2,0,2,2,2]]), 10);
alpha := C`Alpha;
beta :=  Length(C)-alpha;
gamma := Z2Z4Type(C)[3];
delta := Z2Z4Type(C)[4];
kappa := Z2Z4Type(C)[5];
grayMap :=  GrayMap(Z2Z4AdditiveUniverseCode(alpha, beta));
R := RSpace(Z4, alpha+beta);
V := VectorSpace(GF(2), alpha+2*beta);
d := 10;   // error correcting capability t=4
           // dK=10, mimimum weight of the kernel
u := R![0,0,0,0,2,0,2,2,2,0,0,0,0,1,0,0,0,1,2,3,0,0,1,0,2,3,0,3,2,3]; // u in C
ubin := grayMap(u);
v := u+R![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2]; // u not in C, 4 errors
vbin := grayMap(v);
w := u+R![2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]; // u not in C, 4 errors
wbin := grayMap(w);
a := u+R![0,0,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0]; // u not in C, 4 errors
abin := grayMap(a);
b := u+R![2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]; // u not in C, 4 errors
bbin := grayMap(b);
Is := [1,24,25,28,27,26,29,30];           // given information set for C
Isbin := [1,37,39,45,43,41,47,48,49,50];  // given information set for Cbin
I := [3] cat [24..30];                    // another information set
Ibin := [3,37,39,41,43,45] cat [47..50];  // another information set
J := [1..12];                             // not an information set
Jbin := [6..15];                          // not an information set
        
// Test  InformationSpace function
OutputR, OutputV, Outputf, Outputfbin :=  InformationSpace(C);
assert OutputR eq RSpace(LinearCode(DiagonalMatrix(Z4,[2^^gamma] cat [1^^delta])));
assert OutputV eq VectorSpace(GF(2), gamma+2*delta);

// Test  InformationSet function
OutputInfoSet_C, OutputInfoSet_Cbin :=  InformationSet(C);
assert OutputInfoSet_C eq Is;
assert OutputInfoSet_Cbin eq Isbin;

// Test IsInformationSet function
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, Is);
assert OutputIsInfoSet_C eq true;
assert OutputIsInfoSet_Cbin eq false;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, Isbin);
assert OutputIsInfoSet_C eq false;
assert OutputIsInfoSet_Cbin eq true;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, I);
assert OutputIsInfoSet_C eq true;
assert OutputIsInfoSet_Cbin eq false;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, Ibin);
assert OutputIsInfoSet_C eq false;
assert OutputIsInfoSet_Cbin eq true;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, J);
assert OutputIsInfoSet_C eq false;
assert OutputIsInfoSet_Cbin eq false;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, Jbin);
assert OutputIsInfoSet_C eq false;
assert OutputIsInfoSet_Cbin eq false;

// Test  SyndromeSpace function
OutputRs, OutputVs :=  SyndromeSpace(C);
assert OutputRs eq RSpace(LinearCode(DiagonalMatrix(Z4,[2^^(alpha+gamma-2*kappa)] 
                                             cat [1^^(beta-gamma-delta+kappa)])));
assert OutputVs eq VectorSpace(GF(2), alpha+2*beta-gamma-2*delta);

// Test  Syndrome function
OutputSyndrome_u :=  Syndrome(u, C);
OutputSyndrome_ubin :=  Syndrome(ubin, C);
OutputSyndrome_v :=  Syndrome(v, C);
OutputSyndrome_vbin :=  Syndrome(vbin, C);
OutputSyndrome_w :=  Syndrome(w, C);
OutputSyndrome_wbin :=  Syndrome(wbin, C);
OutputSyndrome_a :=  Syndrome(a, C);
OutputSyndrome_abin :=  Syndrome(abin, C);
OutputSyndrome_b :=  Syndrome(b, C);
OutputSyndrome_bbin :=  Syndrome(bbin, C);
assert OutputSyndrome_u eq OutputRs!0;
assert OutputSyndrome_ubin eq OutputRs!0;
assert OutputSyndrome_v eq OutputRs![0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,2,0,2,2,2,0,0,2,2,0,0];
assert OutputSyndrome_vbin eq OutputRs![0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,2,0,2,2,2,0,0,2,2,0,0];
assert OutputSyndrome_w eq OutputRs![0,0,0,2,0,2,0,0,0,2,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0];
assert OutputSyndrome_wbin eq OutputRs![0,0,0,2,0,2,0,0,0,2,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0];
assert OutputSyndrome_a eq OutputRs![0,0,0,2,0,0,0,0,0,2,0,0,0,0,0,2,0,2,2,2,2,2,0,0,0,0,2];
assert OutputSyndrome_abin eq OutputRs![0,0,0,2,0,0,0,0,0,2,0,0,0,0,0,2,0,2,2,2,2,2,0,0,0,0,2];
assert OutputSyndrome_b eq OutputRs![0,0,0,2,0,2,0,0,0,2,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0];
assert OutputSyndrome_bbin eq OutputRs![0,0,0,2,0,2,0,0,0,2,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0];
assert OutputSyndrome_u eq Z2Z4Mult(u, Transpose( MinRowsParityCheckMatrix(C)), C`Alpha);
assert OutputSyndrome_v eq Z2Z4Mult(v, Transpose( MinRowsParityCheckMatrix(C)), C`Alpha);
assert OutputSyndrome_w eq Z2Z4Mult(w, Transpose( MinRowsParityCheckMatrix(C)), C`Alpha);
assert OutputSyndrome_a eq Z2Z4Mult(a, Transpose( MinRowsParityCheckMatrix(C)), C`Alpha);
assert OutputSyndrome_b eq Z2Z4Mult(b, Transpose( MinRowsParityCheckMatrix(C)), C`Alpha);

// Test  CosetLeaders function     There are too many cosets

// Test  CosetDecode function
OutputIsDecoded_u, OutputDecoded_u, OutputDecoded_ubin :=  CosetDecode(C, u);
assert OutputIsDecoded_u eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;
OutputIsDecoded_ubin, OutputDecoded_u, OutputDecoded_ubin :=  CosetDecode(C, ubin);
assert OutputIsDecoded_ubin eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;
OutputIsDecoded_v, OutputDecoded_v, OutputDecoded_vbin :=  CosetDecode(C, v);
assert OutputIsDecoded_v eq true;
assert OutputDecoded_v eq u;
assert OutputDecoded_vbin eq ubin;
OutputIsDecoded_vbin, OutputDecode_v, OutputDecoded_vbin :=  CosetDecode(C, vbin);
assert OutputIsDecoded_vbin eq true;
assert OutputDecoded_v eq u;
assert OutputDecoded_vbin eq ubin;
OutputIsDecoded_w, OutputDecoded_w, OutputDecoded_wbin :=  CosetDecode(C, w);
assert OutputIsDecoded_w eq true;
assert OutputDecoded_w eq u;
assert OutputDecoded_wbin eq ubin;
OutputIsDecoded_wbin, OutputDecode_w, OutputDecoded_wbin :=  CosetDecode(C, wbin);
assert OutputIsDecoded_wbin eq true;
assert OutputDecoded_w eq u;
assert OutputDecoded_wbin eq ubin;
OutputIsDecoded_a, OutputDecoded_a, OutputDecoded_abin :=  CosetDecode(C, a);
assert OutputIsDecoded_a eq true;
assert OutputDecoded_a eq u;
assert OutputDecoded_abin eq ubin;
OutputIsDecoded_abin, OutputDecode_a, OutputDecoded_abin :=  CosetDecode(C, abin);
assert OutputIsDecoded_abin eq true;
assert OutputDecoded_a eq u;
assert OutputDecoded_abin eq ubin;
OutputIsDecoded_b, OutputDecoded_b, OutputDecoded_bbin :=  CosetDecode(C, b);
assert OutputIsDecoded_b eq true;
assert OutputDecoded_b eq u;
assert OutputDecoded_bbin eq ubin;
OutputIsDecoded_bbin, OutputDecode_b, OutputDecoded_bbin :=  CosetDecode(C, bbin);
assert OutputIsDecoded_bbin eq true;
assert OutputDecoded_b eq u;
assert OutputDecoded_bbin eq ubin;

// Test  SyndromeDecode function      There are too many cosets

/****************************************************************/  
print "test 8: A Hadamard Z2Z4-additive code
               alpha = 4, beta = 6, gamma = 1, delta = 2, kappa = 1, length = 10, #C = 32";
C := Z2Z4HadamardCode(2, 4);
alpha := C`Alpha;
beta :=  Length(C)-alpha;
gamma := Z2Z4Type(C)[3];
delta := Z2Z4Type(C)[4];
kappa := Z2Z4Type(C)[5];
grayMap :=  GrayMap(Z2Z4AdditiveUniverseCode(alpha, beta));
R := RSpace(Z4, alpha+beta);
V := VectorSpace(GF(2), alpha+2*beta);
d := 8;   // error correcting capability t=3
          // dK=8, mimimum weight of the kernel
u := R![2,2,0,0,3,3,0,1,2,3]; // u in C
ubin := grayMap(u);
v := u+R![0,0,0,0,0,0,0,0,1,2]; // u not in C, 3 errors
vbin := grayMap(v);
w := u+R![2,2,2,2,0,0,1,2,0,0]; // u not in C, 7 errors
wbin := grayMap(w);
a := u+R![0,0,2,2,2,2,2,2,0,0]; // u not in C, 10 errors
abin := grayMap(a);
Is := [1,9,10];            // given information set for C
Isbin := [1,13,14,15,16];  // given information set for Cbin
I := [2,8,9];              // another information set
Ibin := [2,11,12,13,14];   // another information set
J := [1..12];              // not an information set
Jbin := [6..15];           // not an information set

// Test  InformationSpace function
OutputR, OutputV, Outputf, Outputfbin :=  InformationSpace(C);
assert OutputR eq RSpace(LinearCode(DiagonalMatrix(Z4,[2^^gamma] cat [1^^delta])));
assert OutputV eq VectorSpace(GF(2), gamma+2*delta);

// Test  InformationSet function
OutputInfoSet_C, OutputInfoSet_Cbin :=  InformationSet(C);
assert OutputInfoSet_C eq Is;
assert OutputInfoSet_Cbin eq Isbin;

// Test IsInformationSet function
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, Is);
assert OutputIsInfoSet_C eq true;
assert OutputIsInfoSet_Cbin eq false;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, Isbin);
assert OutputIsInfoSet_C eq false;
assert OutputIsInfoSet_Cbin eq true;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, I);
assert OutputIsInfoSet_C eq true;
assert OutputIsInfoSet_Cbin eq false;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, Ibin);
assert OutputIsInfoSet_C eq false;
assert OutputIsInfoSet_Cbin eq true;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, J);
assert OutputIsInfoSet_C eq false;
assert OutputIsInfoSet_Cbin eq false;
OutputIsInfoSet_C, OutputIsInfoSet_Cbin := IsInformationSet(C, Jbin);
assert OutputIsInfoSet_C eq false;
assert OutputIsInfoSet_Cbin eq false;

// Test  SyndromeSpace function
OutputRs, OutputVs :=  SyndromeSpace(C);
assert OutputRs eq RSpace(LinearCode(DiagonalMatrix(Z4,[2^^(alpha+gamma-2*kappa)] 
                                             cat [1^^(beta-gamma-delta+kappa)])));
assert OutputVs eq VectorSpace(GF(2), alpha+2*beta-gamma-2*delta);

// Test  Syndrome function
OutputSyndrome_u :=  Syndrome(u, C);
OutputSyndrome_ubin :=  Syndrome(ubin, C);
OutputSyndrome_v :=  Syndrome(v, C);
OutputSyndrome_vbin :=  Syndrome(vbin, C);
OutputSyndrome_w :=  Syndrome(w, C);
OutputSyndrome_wbin :=  Syndrome(wbin, C);
OutputSyndrome_a :=  Syndrome(a, C);
OutputSyndrome_abin :=  Syndrome(abin, C);
assert OutputSyndrome_u eq OutputRs!0;
assert OutputSyndrome_ubin eq OutputRs!0;
assert OutputSyndrome_v eq OutputRs![0,2,2,0,1,3,1];
assert OutputSyndrome_vbin eq OutputRs![0,2,2,0,1,3,1];
assert OutputSyndrome_w eq OutputRs![0,0,0,2,1,2,2];
assert OutputSyndrome_wbin eq OutputRs![0,0,0,2,1,2,2];
assert OutputSyndrome_a eq OutputRs![2,0,2,2,2,0,0];
assert OutputSyndrome_abin eq OutputRs![2,0,2,2,2,0,0];
assert OutputSyndrome_u eq Z2Z4Mult(u, Transpose( MinRowsParityCheckMatrix(C)), C`Alpha);
assert OutputSyndrome_v eq Z2Z4Mult(v, Transpose( MinRowsParityCheckMatrix(C)), C`Alpha);
assert OutputSyndrome_w eq Z2Z4Mult(w, Transpose( MinRowsParityCheckMatrix(C)), C`Alpha);
assert OutputSyndrome_a eq Z2Z4Mult(a, Transpose( MinRowsParityCheckMatrix(C)), C`Alpha);

// Test  CosetLeaders function
OutputCosetLeadersL, OutputCosetLeadersMap :=  CosetLeaders(C);
assert #OutputCosetLeadersL eq #OutputRs;
assert OutputCosetLeadersMap(OutputRs!0) eq C`Code!0;
assert OutputCosetLeadersMap(OutputRs![0,2,2,0,1,3,1]) eq v-u;
assert OutputCosetLeadersMap(OutputRs![0,0,0,2,1,2,2]) eq R![0,0,0,0,0,0,1,0,0,2];
assert OutputCosetLeadersMap(OutputRs![2,0,2,2,2,0,0]) eq R![0,0,0,0,3,3,0,1,0,1];

// Test  CosetDecode function
OutputIsDecoded_u, OutputDecoded_u, OutputDecoded_ubin :=  CosetDecode(C, u);
assert OutputIsDecoded_u eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;
OutputIsDecoded_ubin, OutputDecoded_u, OutputDecoded_ubin :=  CosetDecode(C, ubin);
assert OutputIsDecoded_ubin eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;
OutputIsDecoded_v, OutputDecoded_v, OutputDecoded_vbin :=  CosetDecode(C, v);
assert OutputIsDecoded_v eq true;
assert OutputDecoded_v eq u;
assert OutputDecoded_vbin eq ubin;
OutputIsDecoded_vbin, OutputDecode_v, OutputDecoded_vbin :=  CosetDecode(C, vbin);
assert OutputIsDecoded_vbin eq true;
assert OutputDecoded_v eq u;
assert OutputDecoded_vbin eq ubin;
OutputIsDecoded_w, OutputDecoded_w, OutputDecoded_wbin :=  CosetDecode(C, w);
assert OutputIsDecoded_w eq true;
assert OutputDecoded_w eq R![0,0,2,2,3,3,0,3,2,1];
assert OutputDecoded_wbin eq grayMap(R![0,0,2,2,3,3,0,3,2,1]);
OutputIsDecoded_wbin, OutputDecode_w, OutputDecoded_wbin :=  CosetDecode(C, wbin);
assert OutputIsDecoded_wbin eq true;
assert OutputDecoded_w eq R![0,0,2,2,3,3,0,3,2,1];
assert OutputDecoded_wbin eq grayMap(R![0,0,2,2,3,3,0,3,2,1]);
OutputIsDecoded_a, OutputDecoded_a, OutputDecoded_abin :=  CosetDecode(C, a);
assert OutputIsDecoded_a eq true;
assert OutputDecoded_a eq R![2^^10];
assert OutputDecoded_abin eq grayMap(R![2^^10]);
OutputIsDecoded_abin, OutputDecode_a, OutputDecoded_abin :=  CosetDecode(C, abin);
assert OutputIsDecoded_abin eq true;
assert OutputDecoded_a eq R![2^^10];
assert OutputDecoded_abin eq grayMap(R![2^^10]);
OutputIsDecoded_a, OutputDecoded_a, OutputDecoded_abin := 
                                      CosetDecode(C, a : MinWeightCode := d);
assert OutputIsDecoded_a eq true;
assert OutputDecoded_a eq R![2^^10];
assert OutputDecoded_abin eq grayMap(R![2^^10]);
OutputIsDecoded_abin, OutputDecode_a, OutputDecoded_abin := 
                                      CosetDecode(C, abin : MinWeightCode := d);
assert OutputIsDecoded_abin eq true;
assert OutputDecoded_a eq R![2^^10];
assert OutputDecoded_abin eq grayMap(R![2^^10]);
OutputIsDecoded, OutputDecoded, OutputDecoded_bin :=  CosetDecode(C, [u,v,w,a]);
assert OutputIsDecoded eq [true, true, true, true];
assert OutputDecoded eq [u, u, R![0,0,2,2,3,3,0,3,2,1], R![2^^10]];
assert OutputDecoded_bin eq [ubin, ubin, grayMap(R![0,0,2,2,3,3,0,3,2,1]), 
                                         grayMap(R![2^^10])];
OutputIsDecoded, OutputDecoded, OutputDecoded_bin :=  CosetDecode(C, [ubin,vbin,wbin,abin]);
assert OutputIsDecoded eq [true, true, true, true];
assert OutputDecoded eq [u, u, R![0,0,2,2,3,3,0,3,2,1], R![2^^10]];
assert OutputDecoded_bin eq [ubin, ubin, grayMap(R![0,0,2,2,3,3,0,3,2,1]), 
                                         grayMap(R![2^^10])];

// Test  SyndromeDecode function
OutputIsDecoded_u, OutputDecoded_u, OutputDecoded_ubin :=  SyndromeDecode(C, u);
assert OutputIsDecoded_u eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;
OutputIsDecoded_ubin, OutputDecoded_u, OutputDecoded_ubin :=  SyndromeDecode(C, ubin);
assert OutputIsDecoded_ubin eq true;
assert OutputDecoded_u eq u;
assert OutputDecoded_ubin eq ubin;
OutputIsDecoded_v, OutputDecoded_v, OutputDecoded_vbin :=  SyndromeDecode(C, v);
assert OutputIsDecoded_v eq true;
assert OutputDecoded_v eq u;
assert OutputDecoded_vbin eq ubin;
OutputIsDecoded_vbin, OutputDecode_v, OutputDecoded_vbin :=  SyndromeDecode(C, vbin);
assert OutputIsDecoded_vbin eq true;
assert OutputDecoded_v eq u;
assert OutputDecoded_vbin eq ubin;
OutputIsDecoded_w, OutputDecoded_w, OutputDecoded_wbin :=  SyndromeDecode(C, w);
assert OutputIsDecoded_w eq true;
assert OutputDecoded_w eq w-R![0,0,0,0,0,0,1,0,0,2];
assert OutputDecoded_wbin eq grayMap(w-R![0,0,0,0,0,0,1,0,0,2]);
OutputIsDecoded_wbin, OutputDecoded_w, OutputDecoded_wbin :=  SyndromeDecode(C, wbin);
assert OutputIsDecoded_wbin eq true;
assert OutputDecoded_w eq w-R![0,0,0,0,0,0,1,0,0,2];
assert OutputDecoded_wbin eq grayMap(w-R![0,0,0,0,0,0,1,0,0,2]);
OutputIsDecoded_a, OutputDecode_a, OutputDecoded_abin :=  SyndromeDecode(C, a);
assert OutputIsDecoded_a eq true;
assert OutputDecoded_a eq a-R![0,0,0,0,3,3,0,1,0,1];
assert OutputDecoded_abin eq grayMap(a-R![0,0,0,0,3,3,0,1,0,1]);
OutputIsDecoded_abin, OutputDecode_a, OutputDecoded_abin :=  SyndromeDecode(C, abin);
assert OutputIsDecoded_abin eq true;
assert OutputDecoded_a eq a-R![0,0,0,0,3,3,0,1,0,1];
assert OutputDecoded_abin eq grayMap(a-R![0,0,0,0,3,3,0,1,0,1]);
OutputIsDecoded, OutputDecoded, OutputDecoded_bin :=  SyndromeDecode(C, [u,v,w,a]);
assert OutputIsDecoded eq [true, true, true, true];
assert OutputDecoded eq [u, u, w-R![0,0,0,0,0,0,1,0,0,2], a-R![0,0,0,0,3,3,0,1,0,1]];
assert OutputDecoded_bin eq [ubin, ubin, grayMap(w-R![0,0,0,0,0,0,1,0,0,2]), 
                                        grayMap(a-R![0,0,0,0,3,3,0,1,0,1])];
OutputIsDecoded, OutputDecoded, OutputDecoded_bin :=  SyndromeDecode(C, [ubin,vbin,wbin,abin]);
assert OutputIsDecoded eq [true, true, true, true];
assert OutputDecoded eq [u, u, w-R![0,0,0,0,0,0,1,0,0,2], a-R![0,0,0,0,3,3,0,1,0,1]];
assert OutputDecoded_bin eq [ubin, ubin, grayMap(w-R![0,0,0,0,0,0,1,0,0,2]), 
                                        grayMap(a-R![0,0,0,0,3,3,0,1,0,1])];

///****************************************************************/
//Add the following example in a future version
//print "test 3: Hadamard code
//               #alpha = 0, beta = 16, binary length = 32, #C = 64";
//C := Z2Z4AdditiveCode(HadamardCodeZ4(3, 5));
//V := VectorSpace(GF(2),32);
//R := RSpace(Integers(4),16);
//d := 16;  // error correcting capability t=7
//          // dK=16, mimimum weight of the kernel
//// u in C, u in Kernel
//u := R![2,0,2,0,0,2,0,2,2,0,2,0,0,2,0,2];
//ubin := V![1,1,0,0,1,1,0,0,0,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0,0,0,1,1,0,0,1,1]; 
//// v in C, v not in Kernel
//v := R![ 1, 0, 3, 2, 0, 3, 2,1, 3, 2,1, 0, 2,1, 0, 3 ];
//vbin := V![0,1,0,0,1,0,1,1,0,0,1,0,1,1,0,1,1,0,1,1,0,1,0,0,1,1,0,1,0,0,1,0]; 
//// w not in C,1 error (u in Kernel plus 1<=t error)
//w := R![ 1, 0, 2, 0, 0, 2, 0, 2, 2, 0, 2, 0, 0, 2, 0, 2 ];
//wbin := V![0,1,0,0,1,1,0,0,0,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0,0,0,1,1,0,0,1,1]; 
//// z not in C, 5 errors (u in Kernel plus 5<=t errors)
//z := R![ 1, 3,1, 3, 3, 2, 0, 2, 2, 0, 2, 0, 0, 2, 0, 2 ];
//zbin := V![0,1,1,0,0,1,1,0,1,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0,0,0,1,1,0,0,1,1]; 
//// a not in C,1 error (v not in Kernel plus 1<=t error)
//a := R![ 2, 0, 3, 2, 0, 3, 2,1, 3, 2,1, 0, 2,1, 0, 3 ];
//abin := V![1,1,0,0,1,0,1,1,0,0,1,0,1,1,0,1,1,0,1,1,0,1,0,0,1,1,0,1,0,0,1,0]; 
//// b not in C, 4 errors (v not in Kernel plus 4<=t errors)
//b := R![ 0, 3, 3, 2,1, 0, 2,1, 3, 2,1, 0, 2,1, 0, 3 ];
//bbin := V![0,0,1,0,1,0,1,1,0,1,0,0,1,1,0,1,1,0,1,1,0,1,0,0,1,1,0,1,0,0,1,0];
//// c not in C,10,errors (u plus 10,errors, t < 10,<= dK-1=15) 
//c := R![ 0, 2, 0, 2, 2, 2, 0, 2, 2, 0, 2, 0, 0, 2, 0, 2 ];
//cbin := V![0,0,1,1,0,0,1,1,1,1,1,1,0,0,1,1,1,1,0,0,1,1,0,0,0,0,1,1,0,0,1,1]; 

