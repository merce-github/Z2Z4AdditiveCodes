///////////////////////////////////////////////////////////////////////////////
/////////       Copyright 2014-2019 Roland Barrolleta, Jaume Pujol      ///////
/////////                     and Merce Villanueva                      ///////
/////////                                                               ///////
/////////       This program is distributed under the terms of GNU      ///////
/////////               General Public License                          ///////
/////////                                                               ///////
///////////////////////////////////////////////////////////////////////////////

/******************************************************************************
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
******************************************************************************/


/************************************************************/
/*                                                          */
/* Project name: Z2Z4-additive codes in MAGMA               */
/* File name: Z2Z4Decode.m                                  */
/*                                                          */
/* Comment: Package developed within the CCSG group         */
/*                                                          */
/* Authors: Roland Barrolleta, Jaume Pujol and              */
/*          Merce Villanueva                                */
/*                                                          */
/* Revision version and last date: v1.0   28-08-2014        */
/*                                 v1.1   07-02-2015        */ 
/*                                 v1.2   09-02-2015        */
/*                                 v2.0   26-02-2018        */
/*                                 v2.1   28-02-2018        */
/*    User defined type            v2.2   30-01-2019        */
/*                                                          */
/************************************************************/
//Uncomment freeze when package finished
freeze;

intrinsic Z2Z4Decode_version() -> SeqEnum
{Return the current version of this package.}

    version := [2,2];
    return version;

end intrinsic;

//needs Z2Z4AdditiveCode.m file
//needs Z2Z4StandardForm.m file

/******************************************************************
    GLOBAL VARIABLES
*******************************************************************/
//Maximum number of elements in a sequence for the current MAGMA distribution
//This number can be changed if your distribution allows longer sequences  
MAXSEQLENGTH := 2^28;

Z := IntegerRing();

///////////////////////////////////////////////////////////////////////////////
//////////// Functions we need from Z2Z4AdditiveCode.m package ////////////////
///////////////////////////////////////////////////////////////////////////////

import "Z2Z4AdditiveCodes.m" : Z4, Z2, F2, zero_Code, 
                               Z2Z4ChangeMatrixZ4toZ2,
                               FromVectorZ2Z4toZ4;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////                 THE INFORMATION SPACE FUNCTIONS                 ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/****************************************************************/
/*                                                              */  
/* Function name: GrayImageInverse                              */
/* Parameters: v, gamma, delta                                  */
/* Function description: Given a binary vector v of length      */
/*   gamma+delta with delta even, and two integer numbers gamma */
/*   and delta, return a vector of length gamma+delta/2 over Z, */
/*   where the first gamma coordinates are the same and the     */
/*   last delta coordinates are the inverse of the Gray map.    */
/*                                                              */
/****************************************************************/
GrayImageInverse := function(v, gamma, delta)
    Z4GrayInvSeq := [0, 3, 1, 2];
    return Vector([Z!v[i] : i in [1..gamma]] cat [Z4GrayInvSeq[Z!v[x - 1] + 
           2*Z!v[x] + 1 where x is 2*i+gamma] : i in [1..delta]]);
end function;

/****************************************************************/
/*                                                              */  
/* Function name: MultiplyByG                                   */
/* Parameters: v, gamma, G                                      */
/* Function description: Given a vector v over Z4, an integer   */
/*   number gamma and a generator matrix G of a code over Z4 of */
/*   type 2^gamma 4^delta with gamma+delta rows, where the      */
/*   first gamma are of order two and the last delta are of     */
/*   order four, return the product v*G.                        */
/*   The first gamma coordinates of v belong to {0,2}.          */
/*   Note that before multipling by G, it is necessary to       */
/*   replace {0,2} with {0,1} in the first gamma coordinates.   */
/*   This function also works when G is the generator matrix of */
/*   a Z2Z4-additive code of type (alpha,beta;gamma,delta;kappa)*/ 
/*                                                              */
/****************************************************************/
MultiplyByG := function(v, gamma, G)
    for i in [1..gamma] do
        v[i] div:= 2;
    end for;
    return v*G;
end function;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4Mult                                      */
/* Parameters: u, A, alpha                                      */
/* Function description: Given a vector u belonging to Z2^alpha */
/*   x Z4^beta represented as an element in Z4^(alpha+beta),    */
/*   where the ones in the first alpha coordinates are          */
/*   represented by twos; and a matrix A over Z4 having alpha+  */
/*   beta rows and the entries in the first alpha rows in {0,2};*/
/*   return the vector u*A.                                     */
/* Input parameters description:                                */
/*   - u : A vector over Z4                                     */
/*   - A : A matrix over Z4                                     */
/*   - alpha : a positive integer                               */
/* Output parameters description:                               */
/*   - A vector over Z4                                         */
/*                                                              */
/* Signature: (<ModTupRngElt> u, <ModMatRngElt> A,              */
/*             <RngIntElt> alpha) -> ModTupRngElt               */
/* Signature: (<ModTupRngElt> u, <AlgMatElt> A,                 */
/*             <RngIntElt> alpha) -> ModTupRngElt               */
/*                                                              */
/****************************************************************/
intrinsic Z2Z4Mult(u::ModTupRngElt, A::ModMatRngElt, alpha::RngIntElt) -> ModTupRngElt
{
Given a vector u belonging to Z2^alpha x Z4^beta represented as an element in 
Z4^(alpha+beta), where the ones in the first alpha coordinates are represented 
by twos; and a matrix A over Z4 having alpha+beta rows and the entries in the 
first alpha rows in [0,2]; return the vector u*A. 
}   
    require (BaseRing(u) eq Z4): "Argument 1 must be a vector over Z4";
    require (BaseRing(A) eq Z4): "Argument 2 must be a matrix over Z4";
    n := Degree(u);
    requirerange alpha, 0, n;
    require (Nrows(A) eq n): "Argument 2 must have", n,"rows";
    require (IsZero(2*ColumnSubmatrix(u, alpha))):
            "First", alpha, "coordinates of argument 1 must be in {0,2}";
    require (IsZero(2*RowSubmatrix(A, alpha))):
            "First", alpha, "rows of argument 2 must be in {0,2}";
    
    for i in [1..alpha] do
        u[i] div:= 2;
    end for;
    return u*A;
    
end intrinsic;

intrinsic Z2Z4Mult(u::ModTupRngElt, A::AlgMatElt, alpha::RngIntElt) -> ModTupRngElt
{
Given a vector u belonging to Z2^alpha x Z4^beta represented as an element in 
Z4^(alpha+beta), where the ones in the first alpha coordinates are represented 
by twos; and a matrix A over Z4 having alpha+beta rows and the entries in the 
first alpha rows in [0,2]; return the vector u*A. 
}   
    require (BaseRing(u) eq Z4): "Argument 1 must be a vector over Z4";
    require (BaseRing(A) eq Z4): "Argument 2 must be a matrix over Z4";
    n := Degree(u);
    requirerange alpha, 0, n;
    require (Nrows(A) eq n): "Argument 2 must have", n,"rows";
    require (IsZero(2*ColumnSubmatrix(u, alpha))):
            "First", alpha, "coordinates of argument 1 must be in {0,2}";
    require (IsZero(2*RowSubmatrix(A, alpha))):
            "First", alpha, "rows of argument 2 must be in {0,2}";
    
    for i in [1..alpha] do
        u[i] div:= 2;
    end for;
    return u*A;
    
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4Mult                                      */
/* Parameters: A, B, alpha                                      */
/* Function description: Given a m x (alpha+beta) matrix A where*/
/*   the rows belong to Z2^alpha x Z4^beta represented as       */
/*   Z4^(alpha+beta), that is, having the entries in the first  */
/*   alpha columns in {0,2}; and a (alpha+beta) x p matrix B    */
/*   over Z4 having the entries in the first alpha rows in {0,  */
/*   2}, return the m x p matrix A*B having in the i-th row the */
/*   vector Z2Z4Mult(u_i, B), where u_i is the i-th row of A.   */
/* Input parameters description:                                */
/*   - A : A matrix over Z4                                     */
/*   - B : A matrix over Z4                                     */
/*   - alpha : a positive integer                               */
/* Output parameters description:                               */
/*   - A matrix over Z4                                         */
/*                                                              */
/* Signature: (<ModMatRngElt> A, <ModMatRngElt> B,              */
/*             <RngIntElt> alpha) -> ModMatRngElt               */
/* Signature: (<AlgMatElt> A, <AlgMatElt> B,                    */
/*             <RngIntElt> alpha) -> AlgMatElt                  */
/*                                                              */
/****************************************************************/
intrinsic Z2Z4Mult(A::ModMatRngElt, B::ModMatRngElt, alpha::RngIntElt) -> ModMatRngElt
{
Given a m x (alpha+beta) matrix A where the rows belong to Z2^alpha x Z4^beta 
represented as Z4^(alpha+beta), that is, having the entries in the first alpha 
columns in [0,2]; and a (alpha+beta) x p matrix B over Z4 having the entries in 
the first alpha rows in [0,2], return the m x p matrix A*B having in the i-th 
row the vector Z2Z4Mult(u_i, B), where u_i is the i-th row of A. 
}   
    require (BaseRing(A) eq Z4): "Argument 1 must be a matrix over Z4";
    require (BaseRing(B) eq Z4): "Argument 2 must be a matrix over Z4";
    require Ncols(A) eq Nrows(B): "Arguments 1 and 2 are incompatible";
    requirerange alpha, 0, Ncols(A);
    require (IsZero(2*ColumnSubmatrix(A, alpha))):
            "First", alpha," columns of argument 1 must be in {0,2}";
    require (IsZero(2*RowSubmatrix(B, alpha))):
            "First", alpha," rows of argument 2 must be in {0,2}";
    
    for i in [1..Nrows(A)] do
        for j in [1..alpha] do
            A[i][j] div:= 2;
        end for;
    end for;
    return A*B;
    
end intrinsic;

intrinsic Z2Z4Mult(A::AlgMatElt, B::AlgMatElt, alpha::RngIntElt) -> AlgMatElt
{
Given a m x (alpha+beta) matrix A where the rows belong to Z2^alpha x Z4^beta 
represented as Z4^(alpha+beta), that is, having the entries in the first alpha 
columns in [0,2]; and a (alpha+beta) x p matrix B over Z4 having the entries in 
the first alpha rows in [0,2], return the m x p matrix A*B having in the i-th 
row the vector Z2Z4Mult(u_i, B), where u_i is the i-th row of A. 
}   
    require (BaseRing(A) eq Z4): "Argument 1 must be a matrix over Z4";
    require (BaseRing(B) eq Z4): "Argument 2 must be a matrix over Z4";
    require Ncols(A) eq Nrows(B): "Arguments 1 and 2 are incompatible";
    requirerange alpha, 0, Ncols(A);
    require (IsZero(2*ColumnSubmatrix(A, alpha))):
            "First", alpha, "columns of argument 1 must be in {0,2}";
    require (IsZero(2*RowSubmatrix(B, alpha))):
            "First", alpha, "rows of argument 2 must be in {0,2}";
    
    for i in [1..Nrows(A)] do
        for j in [1..alpha] do
            A[i][j] div:= 2;
        end for;
    end for;
    return A*B;
    
end intrinsic;

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
intrinsic InformationSpace(C::Z2Z4Code) -> ModTupRng, ModTupFld, Map, Map 
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma; delta; kappa), return the 
Z4-submodule of Z4^(gamma+delta) isomorphic to Z2^gamma x Z4^delta such that
the first gamma coordinates are of order two, that is, the space of information 
vectors for C. The function also returns the (gamma+2*delta)-dimensional binary 
vector space, which is the space of information vectors for the corresponding 
binary code Cbin=Phi(C), where Phi is the Gray map. Finally, for the encoding 
process, it also returns the corresponding isomorphisms f and fbin from these 
spaces of information vectors onto C and Cbin, respectively. 
}   
    require not(zero_Code(C)): "Code C cannot be the zero code";
    
    G, gamma, delta  := MinRowsGeneratorMatrix(C);
    type := Z2Z4Type(C);
    alpha := type[1];
    beta := type[2];
    kappa := type[5];

    diagonal := [2^^gamma] cat [1^^delta];
    InfCode := Z2Z4AdditiveCode(DiagonalMatrix(Z4, diagonal), kappa);
    InfRSpace := RSpace(InfCode`Code);  
    InfVSpace := VectorSpace(GF(2), gamma + 2*delta); 
    V := VectorSpace(GF(2), alpha + 2*beta);

    grayMap := GrayMap(C);
    grayMapInfCode := GrayMap(InfCode); 
    
    encodingZ4 := map<InfRSpace -> C`Code | x :-> MultiplyByG(x, gamma, G)>;
    encodingZ2 := map<InfVSpace -> V | 
        x :-> grayMap(GrayImageInverse(x, gamma, delta)*G) >;

    return InfRSpace, InfVSpace, encodingZ4, encodingZ2;

end intrinsic;

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
intrinsic InformationSet(C::Z2Z4Code) -> SeqEnum[RngIntElt], SeqEnum[RngIntElt]
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma; delta; kappa), return an 
information set I=[i_1,...,i_(gamma+delta)] in [1,...,alpha+beta] for C such that 
[i_1,...,i_k] in [1,...,alpha] and the code C punctured on [1,...,alpha+beta] 
minus [i_(gamma+1),..., i_(gamma+delta)] is of type 4^delta, and the corresponding 
information set Phi(I)=[i_1,...,i_k,2*i_(k+1)-1-alpha,...,2*i_gamma-1-alpha, 
2*i_(gamma+1)-1-alpha, 2*i_(gamma+1)-alpha,..., 2*i_(gamma+delta)-1-alpha, 
2*i_(gamma+delta)-alpha] in [1,...,2*alpha+beta] for the binary code Cbin=Phi(C), 
where Phi is the Gray map. The information sets I and Phi(I) are returned as a 
sequence of gamma+delta and gamma+2*delta integers, giving the coordinate positions 
that correspond to the information set of C and Cbin, respectively.

An information set I for C is an ordered set of gamma+delta coordinate positions 
such that |C^I|=2^gamma 4^delta, where C^I=[v^I : v in C] and v^I is the vector 
v restricted to the I coordinates. An information set J for Cbin is an ordered 
set of gamma+2*delta coordinate positions such that |C^J_bin|=2^(gamma+2*delta).
}
   require not(zero_Code(C)): "Code C cannot be the zero code";

   n := Length(C`Code); 
   type := Z2Z4Type(C);
   alpha := type[1];
   gamma := type[3];
   delta := type[4];
   kappa := type[5];
   _, _, _, perm := StandardForm(C);
   permInv := perm^(-1);

   InfSet_kappa := [i^permInv : i in [1..kappa]]; 
   InfSetZ4_gamma := [i^permInv : i in [n-(gamma-kappa+delta)+1..n-delta]];    
   InfSetZ4_delta := [i^permInv : i in [n-delta+1..n]];                   

   InfSetZ2_gamma := [2*i-1-alpha : i in InfSetZ4_gamma];
   InfSetZ2_delta := &cat[ [2*i-1-alpha, 2*i-alpha] : i in InfSetZ4_delta];

   return InfSet_kappa cat InfSetZ4_gamma cat InfSetZ4_delta,  
          InfSet_kappa cat InfSetZ2_gamma cat InfSetZ2_delta;

end intrinsic;

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
/* Signature: (<Z2Z4Code> C, <[RngIntElt]> I) -> BoolElt, BoolElt*/
/*                                                              */
/****************************************************************/
intrinsic IsInformationSet(C::Z2Z4Code, I::[RngIntElt]) -> BoolElt, BoolElt
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma; delta; kappa) and a sequence 
I in [1,...,alpha+beta] or I in [1,...,alpha+2*beta], return true if and only if 
I in [1,...,alpha+beta] is an information set for C. This function also returns 
another boolean, which is true if an only if I in [1,...,alpha+2*beta] is an 
information set for the corresponding binary code Cbin=Phi(C), where Phi is the 
Gray map.

An information set I for C is an ordered set of gamma+delta coordinate positions 
such that |C^I|=2^gamma 4^delta, where C^I=[v^I : v in C] and v^I is the vector 
v restricted to the I coordinates. An information set J for Cbin is an ordered 
set of gamma+2*delta coordinate positions such that |C^J_bin|=2^(gamma+2*delta).         
}
    require not(zero_Code(C)): "Code C cannot be the zero code";
    
    I := Set(I);  // Eliminate repeated coordinate positions in I
    k := #I;
    require (k ge 1): "Argument 2 cannot be an empty sequence";
    alpha := C`Alpha;
    beta := Length(C`Code) - alpha;
    n := alpha + beta;
    nbin := alpha + 2*beta;
    maxI := Max(I);
    require (Min(I) ge 1) and (maxI le nbin): 
          "Argument 2 should be a subset of", [1..alpha+beta], "or a subset of", [1..alpha+2*beta];

    type := Z2Z4Type(C);
    gamma := type[3];
    delta := type[4];
    kappa := type[5];

    // Check whether I is an information set for C or not
    if (maxI le n) and ((gamma + delta) eq k) then
        checkSet := {1..n} diff I;
        typeCp := Z2Z4Type(PunctureCode(C, checkSet));
        gammaCp := typeCp[3];
        deltaCp := typeCp[4];
        isInfSetZ4 := (gamma eq gammaCp) and (delta eq deltaCp);
    else
        isInfSetZ4 := false;
    end if;
    
    // Check whether I is an information set for Cbin or not
    if (gamma + 2*delta) eq k then 
        _, kernel := KernelZ2Code(C);
        checkSet := {1..nbin} diff I;
        punctureKernel := PunctureCode(kernel, checkSet);        
        if Dimension(kernel) eq Dimension(punctureKernel) then
            _, cosetRep := KernelCosetRepresentatives(C); 
            H := Transpose(ParityCheckMatrix(punctureKernel));  
            V := VectorSpace(GF(2), k);
            cosetSyndromes := [];
            isInfSetZ2 := true;
            Isorted := Sort(Setseq(I));
            for j := 1 to #cosetRep do
                punctureCosetRep := V![cosetRep[j][i] : i in Isorted]; 
                syndrome := punctureCosetRep * H; 
                if syndrome in cosetSyndromes then 
                    isInfSetZ2 := false;
                    break;
                else
                    Append(~cosetSyndromes, syndrome);
                end if;
            end for; 
        else
            isInfSetZ2 := false;
        end if; 
    else
        isInfSetZ2 := false;
    end if;

    return isInfSetZ4, isInfSetZ2;

end intrinsic;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////                 THE SYNDROME SPACE FUNCTIONS                    ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

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
intrinsic SyndromeSpace(C::Z2Z4Code) -> ModTupRng, ModTupFld
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma; delta; kappa), return the 
Z4-submodule of Z4^(alpha+beta-gamma-kappa) isomorphic to Z2^(alpha+gamma-2*kappa) x 
Z4^(beta-gamma-delta+kappa) such that the first alpha+gamma-2*kappa coordinates 
are of order two, that is, the space of syndrome vectors for C. The function also 
returns the (alpha+2*beta-gamma-2*delta)-dimensional binary vector space, which is 
the space of syndrome vectors for the corresponding binary code Cbin=Phi(C), where 
Phi is the Gray map. Note that these spaces are computed by using the function 
InformationSpace(C) applied to the additive dual code of C, given by function Dual(C).
}  
    require (#C`Code ne 2^C`Alpha*4^(Length(C`Code)-C`Alpha)) : 
                                           "Code C cannot be the universe code";

    InfRSpace, InfVSpace := InformationSpace(Dual(C));
    return InfRSpace, InfVSpace;

end intrinsic;

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
intrinsic Syndrome(u::ModTupRngElt, C::Z2Z4Code) -> ModTupRngElt
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma; delta; kappa), and a 
vector u from the ambient space V=Z2^alpha x Z4^beta, construct the syndrome 
of u relative to the code C, by using the parity check matrix given by the 
function MinRowsParityCheckMatrix(C). The vector u in V is given as a vector in 
Z4^(alpha+beta), where the ones in the first alpha coordinates are represented by twos.
The syndrome is an element of the syndrome space of C, considered as the 
Z4-submodule of Z4^(alpha+beta-delta-kappa) isomorphic to 
Z2^(alpha+gamma-2*kappa) x Z4^(beta-gamma-delta+kappa) such that the first 
alpha+gamma-2*kappa coordinates are of order two.
}
    require (BaseRing(u) eq Z4): "Argument 1 must be a vector over Z4";
    n := Length(C`Code);
    require (Degree(u) eq n): "Argument 1 must be of length", n;

    H := MinRowsParityCheckMatrix(C);
    
    return MultiplyByG(u, C`Alpha, Transpose(H));

end intrinsic;

/****************************************************************/
intrinsic Syndrome(u::ModTupFldElt, C::Z2Z4Code) -> ModTupRngElt
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma; delta; kappa), and a 
vector Phi(u) from the ambient space V2=Z2^(alpha+2*beta), construct the syndrome 
of u relative to the code C, by using the parity check matrix given by the 
function MinRowsParityCheckMatrix(C). The syndrome is an element of the syndrome space of C, 
considered as the Z4-submodule of Z4^(alpha+beta-delta-kappa) isomorphic to 
Z2^(alpha+gamma-2*kappa) x Z4^(beta-gamma-delta+kappa) such that the first 
alpha+gamma-2*kappa coordinates are of order two.
}
    require (BaseRing(u) eq GF(2)): "Argument 1 must be a vector over GF(2)";
    alpha := C`Alpha;
    beta := Length(C`Code) - alpha;
    nbin := alpha + 2*beta;
    require (Degree(u) eq nbin): "Argument 1 must be of length", nbin;

    grayMap := GrayMap(Z2Z4AdditiveUniverseCode(alpha, beta));
    H := MinRowsParityCheckMatrix(C);

    return MultiplyByG(u@@grayMap, alpha, Transpose(H));

end intrinsic;

/****************************************************************/
intrinsic Syndrome(u::Tup, C::Z2Z4Code) -> Tup
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma; delta; kappa), and a 
vector u from the ambient space V=Z2^alpha x Z4^beta, construct the syndrome 
of u relative to the code C, by using the parity check matrix given by the 
function MinRowsParityCheckMatrix(C). The vector u in V is given as a tuple in the cartesian
product set Z2^alpha x Z4^beta.
The syndrome is an element of the syndrome space of C, considered as the 
Z4-submodule of Z4^(alpha+beta-delta-kappa) isomorphic to 
Z2^(alpha+gamma-2*kappa) x Z4^(beta-gamma-delta+kappa) such that the first 
alpha+gamma-2*kappa coordinates are of order two.
}
    require #u eq 2: "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require Type(u[1]) cmpeq  ModTupRngElt and Type(u[2]) cmpeq ModTupRngElt: 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (BaseRing(u[2]) cmpeq Z4 and BaseRing(u[1]) cmpeq Z2): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    
    alpha := Ncols(u[1]);
    beta := Ncols(u[2]);
    n := Length(C`Code);
    require (alpha eq C`Alpha): "Argument 1 must have ", C`Alpha, " coordinates";
    require (alpha + beta eq n): "Argument 1 must be of length", n;
    u4 := FromVectorZ2Z4toZ4(u);

    H := MinRowsParityCheckMatrix(C);
    
    syn := Matrix(MultiplyByG(u4, C`Alpha, Transpose(H)));
    T := Z2Z4Type(C);
    synAlpha := alpha + T[3] - 2*T[5];
    Z2Z4ChangeMatrixZ4toZ2(~syn, synAlpha);

    return <Vector(Z2, [syn[1][j] : j in [1..synAlpha]]), 
              Vector(Z4, [syn[1][j] : j in [synAlpha+1..Ncols(syn)]])>;

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: CosetLeaders                                  */
/* Parameters: C                                                */
/* Function description: Given a Z2Z4-additive code C of type   */
/*   (alpha,beta;gamma;delta;kappa), with ambient space         */
/*   V=Z2^alpha x Z4^beta, return a set of coset leaders        */
/*   (vectors of minimal Lee weight in their cosets) for C in V */
/*   as an indexed set of vectors from V. The elements in       */
/*   V are represented as elements in                           */
/*   Z4^(alpha+beta) by replacing the ones in the first alpha   */
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
intrinsic CosetLeaders(C::Z2Z4Code) -> SetIndx, Map
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma; delta; kappa), with 
ambient space V=Z2^alpha x Z4^beta, return a set of coset leaders (vectors of 
minimal Lee weight in their cosets) for C in V as an indexed set of vectors from V. 
The elements in V are represented as elements in Z4^(alpha+beta) 
by replacing the ones in the first alpha coordinates by twos. This function also 
returns a map from the syndrome space of C onto the coset leaders (mapping a 
syndrome into its corresponding coset leader). Note that this function is only 
applicable when V and C are small.
}
    alpha := C`Alpha;
    beta := Length(C`Code) - alpha;
    totalNumSyndromes := (2^alpha*4^beta) / #C`Code;
    require (totalNumSyndromes le MAXSEQLENGTH): "Code C has too many cosets";    
    
    H := Transpose(MinRowsParityCheckMatrix(C));
    universeCode := Z2Z4AdditiveUniverseCode(alpha, beta);
    
    if zero_Code(C) then
        allSyndromes := [ MultiplyByG(x, alpha, H) : x in universeCode`Code ];
        allCosetLeaders := {@ x : x in universeCode`Code @};
        mapSyndromeCosets := map<allSyndromes -> allCosetLeaders | 
                                  x :-> allCosetLeaders[Position(allSyndromes, x)] >; 
        return allCosetLeaders, mapSyndromeCosets;
    end if;
    
    nbin := alpha + 2*beta; 
    grayMap := GrayMap(universeCode);     
    V := UniverseCode(GF(2), nbin);
  
    // Sequences that will have all syndromes and coset leaders 
    allSyndromes := [ C`Code!0 * H ]; 
    allCosetLeaders := {@ C`Code!0 @};

    errorCapability := Floor((Z2Z4MinimumLeeDistance(C)-1)/2);
    for i in [1..errorCapability] do

        // All binary vectors of weight i
        vectorsWeight_i := Setseq(Words(V, i));
        numVectorsWeight_i := #vectorsWeight_i;    
        
        // Compute all syndromes up to the error-correcting capability
        for j in [1..numVectorsWeight_i] do
            leaderZ4 := vectorsWeight_i[j]@@grayMap;
            s := MultiplyByG(leaderZ4, alpha, H);
            Append(~allSyndromes, s);
            Include(~allCosetLeaders, leaderZ4);
        end for;
        
    end for;

    i := errorCapability + 1;
    while ((#allSyndromes lt totalNumSyndromes) and (i le nbin)) do

        // All binary vectors of weight i
        vectorsWeight_i := Setseq(Words(V, i));
        numVectorsWeight_i := #vectorsWeight_i;

        // Compute the new sindromes from the leaders of weight i
        j := 1;
        while ((j le numVectorsWeight_i) and (#allSyndromes lt totalNumSyndromes)) do
            leaderZ4 := vectorsWeight_i[j]@@grayMap;
            s := MultiplyByG(leaderZ4, alpha, H);
            if (s notin allSyndromes) then
                Append(~allSyndromes, s);
                Include(~allCosetLeaders, leaderZ4);
            end if;
            j := j + 1;
        end while;
        
        i := i + 1;
    end while;

    mapSyndromeCosets := map<allSyndromes -> allCosetLeaders | 
                                  x :-> allCosetLeaders[Position(allSyndromes, x)] >; 

    return allCosetLeaders, mapSyndromeCosets;

end intrinsic;    

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////                 COSET DECODE FUNCTIONS                          ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/****************************************************************/
/*                                                              */
/* Function name: CosetDecodeAlg                                */
/* Parameters: kernel, cosetRep, d, u                           */
/* Function description: Given the binary kernel and the coset  */
/*   representatives of a binary code Cbin=Phi(C), where C is a */
/*   Z2Z4-additive code of type (alpha,beta;gamma,delta;kappa), */
/*   the minimum Lee distance of C, and a binary vector u from  */
/*   the ambient space V2=Z2^(alpha+2*beta),                    */
/*   attempt to decode u with respect to Cbin. If the decoding  */
/*   algorithm succeeds in computing a vector u' in Cbin as the */
/*   decoded version of u in V2, then the function returns true */
/*   and u'. If the decoding algorithm does not succeed in      */
/*   decoding u, then the function returns false and the zero   */
/*   vector in V2.                                              */
/*   See below the description of the intrinsic CosetDecode     */
/*   for more details.                                          */
/*                                                              */
/****************************************************************/
CosetDecodeAlg := function(kernel, cosetRep, d, u : minWeightKernel := -1)

    if IsEmpty(cosetRep) then
        // Decode binary linear codes 
        if u in kernel then
            return  true, u;
        end if;
        
        // Compute the minimum word of the extend coset Cu = kernel U (kernel+u)
        Cu := LinearCode(VerticalJoin(GeneratorMatrix(kernel), u));
        minCosetWord := MinimumWord(Cu);
        // Check whether the minCosetWord can be in kernel or not
        if Weight(minCosetWord) lt d then
            return true, u + minCosetWord;
        else  
            return false, kernel!0;
        end if;
    else                                
        // Decode binary nonlinear codes 
        minWord := u;
        minWeight := Weight(u);
    
        G := GeneratorMatrix(kernel);
        cosetRepresentatives := [kernel!0] cat cosetRep;
        t := Floor((d-1)/2); // error correcting capability 
       
        for v in cosetRepresentatives do
            // Check if u in C
            Cv := LinearCode(VerticalJoin(G, v)); 
            if u in Cv then
                return true, u;    
            end if;
                    
            // Compute the minimum word of the extend coset 
            // Cu = kernel U (kernel+u+v)
            Cu := LinearCode(VerticalJoin(G, u+v));
            minExCosetWord := MinimumWord(Cu);
            minExCosetWeight := Weight(minExCosetWord);
            if minExCosetWeight lt minWeight then
                minWord := minExCosetWord;
                minWeight := minExCosetWeight;
            end if;
       
            // Check whether the decoding is unique
            if minWeight le t then
                return true, u + minWord;  
            end if;
        end for;
                
        // Check whether the minWord can be in C`Kernel or not
        minWeightKernel := minWeightKernel eq -1 select 
                                   MinimumWeight(kernel) else minWeightKernel;
        if minWeight lt minWeightKernel then  
            return true, u + minWord;               
        else                           
            return false, kernel!0;   
        end if;
    end if;

end function;

/****************************************************************/
/*                                                              */
/* Function name: CosetDecode                                   */
/* Parameters: C, u                                             */
/* Function description: Given a Z2Z4-additive code C of type   */
/*   (alpha, beta; gamma; delta; kappa), and a vector u from the*/
/*   ambient space V=Z2^alpha x Z4^beta or V2=Z2^(alpha+2*beta),*/
/*   attempt to decode u with respect to C. The elements in     */
/*   V are represented as elements in                           */
/*   Z4^(alpha+beta) by replacing the ones in the first alpha   */
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
/* case tuple u                                                 */
/* Signature: (<Z2Z4Code> C, <Tup> u)                           */
/*               -> BoolElt, ModTupRngElt, ModTupFldElt         */
/* case binary u                                                */
/* Signature: (<Z2Z4Code> C, <ModTupFldElt> u)                  */
/*               -> BoolElt, ModTupRngElt, ModTupFldElt         */
/*                                                              */
/* case quaternary sequence Q                                   */
/* Signature: (<Z2Z4Code> C, <[ModTupRngElt]> Q)                */
/*               -> [BoolElt], [ModTupRngElt], [ModTupFldElt]   */
/* case binary sequence Q                                       */
/* Signature: (<Z2Z4Code C, <[ModTupFldElt]> Q)                 */
/*               -> [BoolElt], [ModTupRngElt], [ModTupFldElt]   */
/*                                                              */
/****************************************************************/
intrinsic CosetDecode(C::Z2Z4Code, u::ModTupRngElt : MinWeightCode:=-1,
                     MinWeightKernel:=-1) -> BoolElt, ModTupRngElt, ModTupFldElt
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma; delta; kappa), and a 
vector u from the ambient space V=Z2^alpha x Z4^beta, attempt to decode u with 
respect to C. The elements in V are represented as elements 
in Z4^(alpha+beta) by replacing the ones in the first alpha coordinates by twos. 
If the decoding algorithm succeeds in computing a vector u' in C as the decoded 
version of u in V, then the function returns true, u' and Phi(u'), where Phi is 
the Gray map. If the decoding algorithm does not succeed in decoding u, then the 
function returns false, the zero vector in V and the zero vector in V2.

The coset decoding algorithm considers the binary linear code Cu=Cbin U 
(Cbin+Phi(u)), when Cbin=Phi(C) is linear. On the other hand, when Cbin is 
nonlinear, we have that Cbin=U_(i=0)^t (Kbin + Phi(ci)), where Kbin=Phi(KC), KC 
is the kernel of C as a Z2Z4-additive subcode, [c0,c1,...,ct] are the coset 
representatives of C with respect to KC (not necessarily of minimal weight in 
their cosets) and c0 is the zero codeword. In this case, the algorithm considers 
the binary linear codes K0=Kbin U (Kbin+Phi(u)), K1=Kbin U (Kbin+Phi(c1)+Phi(u)), 
..., Kt=Kbin U (Kbin+Phi(ct)+Phi(u)).

If the parameter MinWeightCode is not assigned, then the minimum Lee weight of C, 
which coincides with the minimum weight of Cbin, denoted by d, is computed. 
Note that the minimum distance of Cbin coincides with its minimum weight.
If Cbin is linear and the minimum weight of Cu is less than d, then 
Phi(u')=Phi(u)+e, where e is a word of minimum weight of Cu; otherwise, the 
decoding algorithm returns false. On the other hand, if Cbin is nonlinear and 
the minimum weight of U_(i=0)^t Ki is less than the minimum weight of Kbin, then 
Phi(u')=Phi(u)+e, where e is a word of minimum weight of U_(i=0)^t Ki; otherwise, 
the decoding algorithm returns false. If the parameter MinWeightKernel is not 
assigned, then the minimum Hamming weight of Kbin is computed.
}
    require (BaseRing(u) eq Z4): "Argument 2 must be a vector over Z4";
    n := Length(C);
    require (Degree(u) eq n): "Argument 2 must be of length ", n;
    alpha := C`Alpha;
    require (IsZero(2*ColumnSubmatrix(u, alpha))):
            "First ", alpha, " coordinates must be in {0,2}";
    
    beta := n - alpha;
    nbin := alpha + 2*beta;
    grayMap := GrayMap(Z2Z4AdditiveUniverseCode(alpha, beta));
    ubin := grayMap(u);
    zero := C`Code!0;
            
    if zero_Code(C) then
        return true, zero, grayMap(zero);
    end if;
        
    d := (MinWeightCode eq -1) select Z2Z4MinimumLeeWeight(C) else MinWeightCode;

    _, kernel := KernelZ2Code(C);
    _, cosetRep := KernelCosetRepresentatives(C);
    isDecoded, ubinDecoded := CosetDecodeAlg(kernel, cosetRep, d, ubin :
                                                minWeightKernel := MinWeightKernel);
 
    if isDecoded then
        uDecoded := ubinDecoded@@grayMap;
        return isDecoded, uDecoded, ubinDecoded;
    else
        return isDecoded, zero, grayMap(zero);
    end if;

end intrinsic;

/****************************************************************/
intrinsic CosetDecode(C::Z2Z4Code, u::Tup : MinWeightCode:=-1,
                     MinWeightKernel:=-1) -> BoolElt, ModTupRngElt, ModTupFldElt
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma; delta; kappa), and a 
vector u from the ambient space V=Z2^alpha x Z4^beta, attempt to decode u with 
respect to C. The elements in V are represented as elements 
in Z4^(alpha+beta) by replacing the ones in the first alpha coordinates by twos. 
If the decoding algorithm succeeds in computing a vector u' in C as the decoded 
version of u in V, then the function returns true, u' and Phi(u'), where Phi is 
the Gray map. If the decoding algorithm does not succeed in decoding u, then the 
function returns false, the zero vector in V and the zero vector in V2.

The coset decoding algorithm considers the binary linear code Cu=Cbin U 
(Cbin+Phi(u)), when Cbin=Phi(C) is linear. On the other hand, when Cbin is 
nonlinear, we have that Cbin=U_(i=0)^t (Kbin + Phi(ci)), where Kbin=Phi(KC), KC 
is the kernel of C as a Z2Z4-additive subcode, [c0,c1,...,ct] are the coset 
representatives of C with respect to KC (not necessarily of minimal weight in 
their cosets) and c0 is the zero codeword. In this case, the algorithm considers 
the binary linear codes K0=Kbin U (Kbin+Phi(u)), K1=Kbin U (Kbin+Phi(c1)+Phi(u)), 
..., Kt=Kbin U (Kbin+Phi(ct)+Phi(u)).

If the parameter MinWeightCode is not assigned, then the minimum Lee weight of C, 
which coincides with the minimum weight of Cbin, denoted by d, is computed. 
Note that the minimum distance of Cbin coincides with its minimum weight.
If Cbin is linear and the minimum weight of Cu is less than d, then 
Phi(u')=Phi(u)+e, where e is a word of minimum weight of Cu; otherwise, the 
decoding algorithm returns false. On the other hand, if Cbin is nonlinear and 
the minimum weight of U_(i=0)^t Ki is less than the minimum weight of Kbin, then 
Phi(u')=Phi(u)+e, where e is a word of minimum weight of U_(i=0)^t Ki; otherwise, 
the decoding algorithm returns false. If the parameter MinWeightKernel is not 
assigned, then the minimum Hamming weight of Kbin is computed.
}
    require #u eq 2: "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require Type(u[1]) cmpeq  ModTupRngElt and Type(u[2]) cmpeq ModTupRngElt: 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (BaseRing(u[2]) cmpeq Z4 and BaseRing(u[1]) cmpeq Z2): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    
    alpha := Ncols(u[1]);
    beta := Ncols(u[2]);
    n := alpha + beta;
    require n gt 0 : "Argument 1 must be a vector of length greater than 0";
    u4 := FromVectorZ2Z4toZ4(u);
    require (Degree(u4) eq Length(C)): "Argument 2 must be of length ", Length(C);

    minWeightCode := MinWeightCode;
    minWeightKernel := MinWeightKernel;
    return CosetDecode(C, u4 : MinWeightCode := minWeightCode, MinWeightKernel := minWeightKernel);

end intrinsic;

/****************************************************************/
intrinsic CosetDecode(C::Z2Z4Code, u::ModTupFldElt : MinWeightCode:=-1,
                     MinWeightKernel:=-1) -> BoolElt, ModTupRngElt, ModTupFldElt
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma; delta; kappa), and a 
vector Phi(u) from the ambient space V2=Z2^(alpha+2*beta), attempt to decode u 
with respect to C. If the decoding algorithm succeeds in computing a vector u' 
in C as the decoded version of u in V, then the function returns true, u' and 
Phi(u'), where Phi is the Gray map. If the decoding algorithm does not succeed 
in decoding u, then the function returns false, the zero vector in V and the 
zero vector in V2.

The coset decoding algorithm considers the binary linear code Cu=Cbin U 
(Cbin+Phi(u)), when Cbin=Phi(C) is linear. On the other hand, when Cbin is 
nonlinear, we have that Cbin=U_(i=0)^t (Kbin + Phi(ci)), where Kbin=Phi(KC), KC 
is the kernel of C as a Z2Z4-additive subcode, [c0,c1,...,ct] are the coset 
representatives of C with respect to KC (not necessarily of minimal weight in 
their cosets) and c0 is the zero codeword. In this case, the algorithm considers 
the binary linear codes K0=Kbin U (Kbin+Phi(u)), K1=Kbin U (Kbin+Phi(c1)+Phi(u)), 
..., Kt=Kbin U (Kbin+Phi(ct)+Phi(u)).

If the parameter MinWeightCode is not assigned, then the minimum Lee weight of C, 
which coincides with the minimum weight of Cbin, denoted by d, is computed. 
Note that the minimum distance of Cbin coincides with its minimum weight.
If Cbin is linear and the minimum weight of Cu is less than d, then 
Phi(u')=Phi(u)+e, where e is a word of minimum weight of Cu; otherwise, the 
decoding algorithm returns false. On the other hand, if Cbin is nonlinear and 
the minimum weight of U_(i=0)^t Ki is less than the minimum weight of Kbin, then 
Phi(u')=Phi(u)+e, where e is a word of minimum weight of U_(i=0)^t Ki; otherwise, 
the decoding algorithm returns false. If the parameter MinWeightKernel is not 
assigned, then the minimum Hamming weight of Kbin is computed.
}
    require (BaseRing(u) eq GF(2)): "Argument 2 must be a vector over GF(2)";
    nbin := BinaryLength(C);
    require (Degree(u) eq nbin): "Argument 2 must be of length ", nbin;
    
    grayMap := GrayMap(Z2Z4AdditiveUniverseCode(C`Alpha, Length(C`Code)-C`Alpha));
    zero := C`Code!0;
    
    if zero_Code(C) then
        return true, zero, grayMap(zero);
    end if;
    
    d := (MinWeightCode eq -1) select Z2Z4MinimumLeeWeight(C) else MinWeightCode;

    _, kernel := KernelZ2Code(C);
    _, cosetRep := KernelCosetRepresentatives(C);
    isDecoded, ubinDecoded := CosetDecodeAlg(kernel, cosetRep, d, u :
                                                minWeightKernel := MinWeightKernel);
    if isDecoded then
        uDecoded := ubinDecoded@@grayMap;
        return isDecoded, uDecoded, ubinDecoded;
    else
        return isDecoded, zero, grayMap(zero);
    end if;

end intrinsic;

/****************************************************************/
intrinsic CosetDecode(C::Z2Z4Code, Q::[ModTupRngElt] : MinWeightCode:=-1,
               MinWeightKernel:=-1) -> SeqEnum[BoolElt], SeqEnum[ModTupRngElt], 
               SeqEnum[ModTupFldElt]
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma; delta; kappa), and a 
sequence Q of vectors from the ambient space V=Z2^alpha x Z4^beta, attempt to 
decode the vectors of Q with respect to C. This function is similar to the function 
CosetDecode(C, u) except that rather than decoding a single vector, it decodes a
sequence of vectors and returns a sequence of booleans and two sequences of decoded 
vectors corresponding to the given sequence. The algorithm used and effect of the 
parameters MinWeightCode and MinWeightKernel are as for the function CosetDecode(C, u).
}
    require not IsEmpty(Q): "Argument 2 cannot be an empty sequence";  
    require (BaseRing(Q[1]) eq Z4): "Argument 2 must contain vectors over Z4";
    n := Length(C);
    require (Degree(Q[1]) eq n): "Argument 2 must contain vectors of length ", n;
    alpha := C`Alpha;
    require (IsZero(2*ColumnSubmatrix(Matrix(Q), alpha))):
            "First ", alpha, " coordinates must be in {0,2}";
      
    beta := n - alpha;
    nbin := alpha + 2*beta;
    grayMap := GrayMap(Z2Z4AdditiveUniverseCode(alpha, beta));
    zero := C`Code!0;
    zerobin := grayMap(zero);
 
    if zero_Code(C) then
        return true, zero, zerobin;
    end if;
    
    d := (MinWeightCode eq -1) select Z2Z4MinimumLeeWeight(C) else MinWeightCode;

    _, kernel := KernelZ2Code(C);
    _, cosetRep := KernelCosetRepresentatives(C);
    
    isDecodedSeq := [];
    uDecodedSeq := [];
    ubinDecodedSeq := [];
    for u in Q do
        ubin := grayMap(u);
        isDecoded, ubinDecoded := CosetDecodeAlg(kernel, cosetRep, d, ubin :
                                                minWeightKernel := MinWeightKernel); 
        Append(~isDecodedSeq, isDecoded); 
        if isDecoded then          
            Append(~uDecodedSeq, ubinDecoded@@grayMap);
            Append(~ubinDecodedSeq, ubinDecoded);
        else
            Append(~uDecodedSeq, zero);
            Append(~ubinDecodedSeq, zerobin);        
        end if; 
    end for;
    
    return isDecodedSeq, uDecodedSeq, ubinDecodedSeq;

end intrinsic;

/****************************************************************/
intrinsic CosetDecode(C::Z2Z4Code, Q::[ModTupFldElt] : MinWeightCode:=-1,
          MinWeightKernel:=-1) -> SeqEnum[BoolElt], SeqEnum[ModTupRngElt], 
          SeqEnum[ModTupFldElt]
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma; delta; kappa), and a 
sequence Q of vectors from the ambient space V2=Z2^(alpha+2*beta), attempt to 
decode the vectors of Q with respect to C. This function is similar to the function 
CosetDecode(C, u) except that rather than decoding a single vector, it decodes a
sequence of vectors and returns a sequence of booleans and two sequences of decoded 
vectors corresponding to the given sequence. The algorithm used and effect of the 
parameters MinWeightCode and MinWeightKernel are as for the function CosetDecode(C, u).
}
    require not IsEmpty(Q): "Argument 2 cannot be an empty sequence";  
    require (BaseRing(Q[1]) eq GF(2)): "Argument 2 must contain vectors over GF(2)";
    nbin := BinaryLength(C);
    require (Degree(Q[1]) eq nbin): "Argument 2 must contain vectors of length ", nbin;
    
    grayMap := GrayMap(Z2Z4AdditiveUniverseCode(C`Alpha, Length(C`Code)-C`Alpha));
    zero := C`Code!0;
    zerobin := grayMap(zero);
    
    if zero_Code(C) then
        return true, zero, zerobin;
    end if;
    
    d := (MinWeightCode eq -1) select Z2Z4MinimumLeeWeight(C) else MinWeightCode;

    _, kernel := KernelZ2Code(C);
    _, cosetRep := KernelCosetRepresentatives(C);
    
    isDecodedSeq := [];
    uDecodedSeq := [];
    ubinDecodedSeq := [];
    for u in Q do
        isDecoded, ubinDecoded := CosetDecodeAlg(kernel, cosetRep, d, u :
                                                minWeightKernel := MinWeightKernel); 
        Append(~isDecodedSeq, isDecoded); 
        if isDecoded then          
            Append(~uDecodedSeq, ubinDecoded@@grayMap);
            Append(~ubinDecodedSeq, ubinDecoded);
        else
            Append(~uDecodedSeq, zero);
            Append(~ubinDecodedSeq, zerobin);        
        end if; 
    end for;
    
    return isDecodedSeq, uDecodedSeq, ubinDecodedSeq;

end intrinsic;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////                SYNDROME DECODE FUNCTIONS                            ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/****************************************************************/
/*                                                              */
/* Function name: SyndromeDecode                                */
/* Parameters: C, u                                             */
/* Function description: Given a Z2Z4-additive code C of type   */
/*   (alpha, beta; gamma; delta; kappa), and a vector u from the*/
/*   ambient space V=Z2^alpha x Z4^beta or V2=Z2^(alpha+2*beta),*/
/*   attempt to decode u with respect to C. The vector u in V   */
/*   can be given either as a vector in Z4^(\alpha+\beta), where*/
/*   the ones in the first alpha coordinates are represented by */
/*   twos; or as a tuple in the cartesian product set Z2^alpha x*/
/*   Z4^beta. The decoding algorithm always succeeds            */
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
/* case tuple u                                                 */
/* Signature: (<Z2Z4Code> C, <Tup> u)                           */
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
intrinsic SyndromeDecode(C::Z2Z4Code, u::ModTupRngElt) 
                                        -> BoolElt, ModTupRngElt, ModTupFldElt
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma; delta; kappa), and a 
vector u from the ambient space V=Z2^alpha x Z4^beta, attempt to decode u with 
respect to C. The vector u in V is given as a vector in Z4^(\alpha+\beta), where
the ones in the first alpha coordinates are represented by twos.
The decoding algorithm always succeeds in computing a vector u' in C as the 
decoded version of u in V, and the function returns true, u' and Phi(u'), where 
Phi is the Gray map. Although the function never returns false, the first output 
parameter true is given to be consistent with the other decoding functions.

The syndrome decoding algorithm consists of computing a table pairing each possible 
syndrome s with a vector of minimum Lee weight e_s, called coset leader, in the 
coset of C containing all vectors having syndrome s. After receiving a vector u, 
its syndrome s is computed using the parity check matrix. Then, u is decoded into the 
codeword c=u-e_s.
}
    require (BaseRing(u) eq Z4): "Argument 2 must be a vector over Z4";
    n := Length(C);
    require (Degree(u) eq n): "Argument 2 must be of length ", n;
    alpha := C`Alpha;
    require (IsZero(2*ColumnSubmatrix(u, alpha))):
            "First ", alpha, " coordinates must be in {0,2}";

    _, mapCosetLeaders := CosetLeaders(C);
    grayMap := GrayMap(C);
    
    c := u - mapCosetLeaders(Syndrome(u, C));

    return true, c, grayMap(c);

end intrinsic;

/****************************************************************/
intrinsic SyndromeDecode(C::Z2Z4Code, u::Tup) 
                                        -> BoolElt, ModTupRngElt, ModTupFldElt
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma; delta; kappa), and a 
vector u from the ambient space V=Z2^alpha x Z4^beta, attempt to decode u with 
respect to C. The vector u in V is given as a tuple in the cartesian product 
set Z2^alpha x Z4^beta. The decoding algorithm always succeeds in computing a 
vector u' in C as the decoded version of u in V, and the function returns true, 
u' and Phi(u'), where Phi is the Gray map. Although the function never returns 
false, the first output parameter true is given to be consistent with the other 
decoding functions.

The syndrome decoding algorithm consists of computing a table pairing each possible 
syndrome s with a vector of minimum Lee weight e_s, called coset leader, in the 
coset of C containing all vectors having syndrome s. After receiving a vector u, 
its syndrome s is computed using the parity check matrix. Then, u is decoded into the 
codeword c=u-e_s.
}
    require #u eq 2: "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require Type(u[1]) cmpeq  ModTupRngElt and Type(u[2]) cmpeq ModTupRngElt: 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (BaseRing(u[2]) cmpeq Z4 and BaseRing(u[1]) cmpeq Z2): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    
    alpha := Ncols(u[1]);
    beta := Ncols(u[2]);
    n := alpha + beta;
    require n gt 0 : "Argument 1 must be a vector of length greater than 0";
    u4 := FromVectorZ2Z4toZ4(u);
    require (Degree(u4) eq Length(C)): "Argument 2 must be of length ", Length(C);

    return SyndromeDecode(C, u4);

end intrinsic;

/****************************************************************/
intrinsic SyndromeDecode(C::Z2Z4Code, u::ModTupFldElt) 
                                        -> BoolElt, ModTupRngElt, ModTupFldElt
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma; delta; kappa), and a 
vector Phi(u) from the ambient space V=Z2^(alpha+2*beta), attempt to decode u with 
respect to C. The decoding algorithm always succeeds in computing a vector u' in C 
as the decoded version of u in V, and the function returns true, u' and Phi(u'), 
where Phi is the Gray map. Although the function never returns false, the first 
output parameter true is given to be consistent with the other decoding functions.

The syndrome decoding algorithm consists of computing a table pairing each possible 
syndrome s with a vector of minimum Lee weight e_s, called coset leader, in the 
coset of C containing all vectors having syndrome s. After receiving a vector u, 
its syndrome s is computed using the parity check matrix. Then, u is decoded into the 
codeword c=u-e_s.
}
    require (BaseRing(u) eq GF(2)): "Argument 2 must be a vector over GF(2)";
    nbin := BinaryLength(C);
    require (Degree(u) eq nbin): "Argument 2 must be of length ", nbin;
    
    _, mapCosetLeaders := CosetLeaders(C);
    grayMap := GrayMap(Z2Z4AdditiveUniverseCode(C`Alpha, Length(C`Code)-C`Alpha));

    c := u@@grayMap - mapCosetLeaders(Syndrome(u, C));

    return true, c, grayMap(c);

end intrinsic

/****************************************************************/
intrinsic SyndromeDecode(C::Z2Z4Code, Q::[ModTupRngElt]) 
               -> SeqEnum[BoolElt], SeqEnum[ModTupRngElt], SeqEnum[ModTupFldElt]
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma; delta; kappa), and a 
sequence Q of vectors from the ambient space V=Z2^alpha x Z4^beta, attempt to 
decode the vectors of Q with respect to C. This function is similar to the function 
SyndromeDecode(C, u) except that rather than decoding a single vector, it 
decodes a sequence of vectors and returns a sequence of booleans and two sequences 
of decoded vectors corresponding to the given sequence. The algorithm used is 
the same as that of function SyndromeDecode(C, u).
}
    require not IsEmpty(Q): "Argument 2 cannot be an empty sequence";  
    require (BaseRing(Q[1]) eq Z4): "Argument 2 must contain vectors over Z4";
    n := Length(C);
    require (Degree(Q[1]) eq n): "Argument 2 must contain vectors of length ", n;
    alpha := C`Alpha;
    require (IsZero(2*ColumnSubmatrix(Matrix(Q), alpha))):
            "First ", alpha, " coordinates must be in {0,2}";

    _, mapCosetLeaders := CosetLeaders(C);
    grayMap := GrayMap(C);
    
    uDecodedSeq := [];
    ubinDecodedSeq := [];
    for u in Q do
        c := u - mapCosetLeaders(Syndrome(u, C));           
        Append(~uDecodedSeq, c);
        Append(~ubinDecodedSeq, grayMap(c));
    end for;
    
    return [true : i in [1..#Q]], uDecodedSeq, ubinDecodedSeq;   

end intrinsic;

/****************************************************************/
intrinsic SyndromeDecode(C::Z2Z4Code, Q::[ModTupFldElt]) 
               -> SeqEnum[BoolElt], SeqEnum[ModTupRngElt], SeqEnum[ModTupFldElt]
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma; delta; kappa), and a 
sequence Q of vectors from the ambient space V=Z2^(alpha+2*beta), attempt to 
decode the vectors of Q with respect to C. This function is similar to the function 
SyndromeDecode(C, u) except that rather than decoding a single vector, it 
decodes a sequence of vectors and returns a sequence of booleans and two sequences 
of decoded vectors corresponding to the given sequence. The algorithm used is 
the same as that of function SyndromeDecode(C, u).
}
    require not IsEmpty(Q): "Argument 2 cannot be an empty sequence";  
    require (BaseRing(Q[1]) eq GF(2)): "Argument 2 must contain vectors over GF(2)";
    nbin := BinaryLength(C);
    require (Degree(Q[1]) eq nbin): "Argument 2 must contain vectors of length ", nbin;
    
    _, mapCosetLeaders := CosetLeaders(C);
    grayMap := GrayMap(Z2Z4AdditiveUniverseCode(C`Alpha, Length(C`Code)-C`Alpha));

    uDecodedSeq := [];
    ubinDecodedSeq := [];
    for u in Q do
        c := u@@grayMap - mapCosetLeaders(Syndrome(u, C));         
        Append(~uDecodedSeq, c);
        Append(~ubinDecodedSeq, grayMap(c));
    end for;
    
    return [true : i in [1..#Q]], uDecodedSeq, ubinDecodedSeq;  

end intrinsic;

/****************************************************************/ 
/*        Versions using the package for binary codes           */
/*        which are not included in Magma distribution          */             
/****************************************************************/
//Z2Z4AlphaFromZ4toGF2 := pmap< Integers(4) -> Integers(2) | [<0,0>,<2,1>]>;
//Z2Z4AlphaFromGF2toZ4 := pmap< Integers(2) -> Integers(4) | [<0,0>,<1,2>]>;
//Z2Z4GraySeq := [[0,0], [0,1], [1,1], [1,0]];
//Z2Z4GrayInvSeq := [0, 3, 1, 2];
//Z := IntegerRing();
//
//function Z2Z4GrayImageVector(v,alpha)
//    return Vector(GF(2),([Z2Z4AlphaFromZ4toGF2(v[i]) : i in [1..alpha]]) cat 
//           (&cat[Z2Z4GraySeq[Z!v[i]+1] : i in [alpha+1..Degree(v)]]));
//end function;
//
//function Z2Z4GrayImageVectorInverse(y,alpha,beta)
//    alphaseq := [Z2Z4AlphaFromGF2toZ4(y[i]) : i in [1..alpha]];
//    betaseq := [y[i] : i in [alpha+1..2*beta+alpha]];
//    return Vector(alphaseq cat [Z2Z4GrayInvSeq[Z!betaseq[2*i - 1] + 
//           2*Z!betaseq[2*i] + 1] : i in [1..beta]]);
//end function; 
//intrinsic Z2Z4CosetDecode(C::Rec, v::ModTupRngElt : minWeightCode:=-1, 
//                     minWeightKernel:=-1) -> BoolElt, ModTupRngElt, ModTupFldElt
//{
//}
//    require not(zero_Code(C)): "C cannot be the zero code";
//    require (BaseRing(v) eq Z4): "Argument 2 is not a vector over Z4";
//    n := Z2Z4Length(C);
//    require Degree(v) eq n: "Argument 2 must be of length ",n;
//    alpha := C`Alpha;
//    require (IsZero(2*ColumnSubmatrix(v,alpha))):
//            "First",alpha,"coordinates must be in {0,2}";
//
//    Cbin := BinaryCode(C);
//    nbin := BinaryLength(Cbin);
//    vbin := Z2Z4GrayImageVector(v, C`Alpha); 
//    
//    if minWeightCode eq -1 then
//        if alpha eq 0 then 
//            d := MinimumLeeWeight(C`Code);
//        else
//            // The minimum distance coincides with the minimum weight
//            d := BinaryMinimumWeight(Cbin);
//        end if;
//    else
//        d := minWeightCode;
//    end if;
//    
//    isDecoded, vbinDecoded := BinaryCosetDecode(Cbin, d, vbin : 
//                                                minWeightKernel := minWeightKernel);
//    if isDecoded then
//        vDecoded := Z2Z4GrayImageVectorInverse(vbinDecoded, alpha, n-alpha);
//        return isDecoded, vDecoded, vbinDecoded;
//    else
//        return isDecoded, C`Code!0, C`Kernel!0;
//    end if;
//    
//end intrinsic;    
//
///****************************************************************/ 
//intrinsic Z2Z4CosetDecode(C::Rec, v::ModTupFldElt : minWeightCode:=-1, 
//                     minWeightKernel:=-1) -> BoolElt, ModTupRngElt, ModTupFldElt 
//{
//}
//    require over_Z2Z4(C): "C is not a Z2Z4-additive code";
//    require not(zero_Code(C)): "C cannot be the zero code";
//    require (BaseRing(v) eq GF(2)): "Argument 2 is not a vector over GF(2)";
//    nbin := Z2Z4BinaryLength(C);
//    require Degree(v) eq nbin: "Argument 2 must be of length ", nbin;
//    
//    alpha := C`Alpha;
//    Cbin := BinaryCode(C);
//    //d := (minWeightCode eq -1) select BinaryMinimumDistance(Cbin) else minWeightCode;
//    if minWeightCode eq -1 then
//        if alpha eq 0 then 
//            d := MinimumLeeWeight(C`Code);
//        else
//            // The minimum distance coincides with the minimum weight
//            d := BinaryMinimumWeight(Cbin);
//        end if;
//    else
//        d := minWeightCode;
//    end if;
//    
//    isDecoded, vbinDecoded := BinaryCosetDecode(Cbin, d, v : 
//                                                minWeightKernel := minWeightKernel);
//    if isDecoded then
//        vDecoded := Z2Z4GrayImageVectorInverse(vbinDecoded, alpha, Z2Z4Length(C)-alpha); 
//        return isDecoded, vDecoded, vbinDecoded;
//    else
//        return isDecoded, C`Code!0, C`Kernel!0;
//    end if;
//    
//end intrinsic;    







