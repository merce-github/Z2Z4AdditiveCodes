///////////////////////////////////////////////////////////////////////////////
/////////       Copyright 2007-2019 Bernat Gaston, Jaume Pujol          ///////
/////////         and Merce Villanueva (with contributions of           ///////
/////////         Lorena Ronquillo, Jaume Pernas and Adrián Torres)     ///////
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
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
******************************************************************************/

 
/*************************************************************/
/*                                                           */
/* Project name: Z2Z4-additive codes in MAGMA                */
/* File name: Z2Z4CodeConstructions.m                        */
/*                                                           */
/* Comment: Package developed within the CCSG group          */
/*                                                           */
/* Authors: Bernat Gaston, Jaume Pujol and Merce Villanueva  */
/*          (with contributions of Lorena Ronquillo,         */
/*          Jaume Pernas and Adrian Torres)                  */
/*                                                           */
/* Revision version and last date: v1.0   14-03-2019         */
/*                                 v1.1   20-03-2019         */
/*              fix a bug in 'eq'  v1.2   20-05-2021         */
/*       add Reed-Mulller section  v1.3   31-05-2022         */
/*                                                           */
/*************************************************************/
//Uncomment freeze when package finished
freeze;

intrinsic Z2Z4CodeConstructions_version() -> SeqEnum
{Return the current version of this package.}
    version := [1, 3];
    return version;
end intrinsic;

////////////////////////////////////////////////////////////////////////////////
///////// Functions we need from Z2Z4AdditiveCode.m package ////////////////////
////////////////////////////////////////////////////////////////////////////////

import "Z2Z4AdditiveCodes.m": IsZ2Z4AlphaOverZ4,  
                              Z2Z4AlphaFromZ4toZ2,
                              Z2Z4AlphaFromZ2toZ4,
                              Z2Z4ChangeMatrixZ2toZ4,
                              Z2Z4ChangeMatrixZ4toZ2,
                              FromVectorZ2Z4toZ4,
                              FromSeqZ2Z4toZ4,
                              OrderTwoFourGenerators,
                              Z2Z4GrayVec,
                              UpdateMinimumLeeWeight,
                              UpdateMinimumLeeWeightLowerBound,
                              UpdateMinimumLeeWeightUpperBound;

/////        Type declaration

alphaMatrix := recformat<Alpha:RngIntElt, Matrix:Mtrx>;

/////        Global variables 
       
Z4 := Integers(4);
Z2 := Integers(2);
F2 := GF(2);

/*******************************************************************************
                            Checks
*******************************************************************************/

over_Z4 := func<C | Alphabet(C) cmpeq Z4>;
over_Z2 := func<C | Alphabet(C) cmpeq Z2>;
over_F2 := func<C | Alphabet(C) cmpeq F2>;

zero_Code := func<C | PseudoDimension(C) eq 0>;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////                SPAN AND KERNEL                                  ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/****************************************************************/
/*                                                              */
/* Procedure name: UpdateMinimumWeightLowerBound                */
/* Parameters: ~C, lowerBound                                   */
/* Procedure description: Given a binary linear code C and an   */
/*   integer with a lower bound for its minimum weight,         */
/*   update the attribute MinimumWeightLowerBound.              */  
/* Input parameters description:                                */
/*   - C: A binary linear code                                  */
/*   - lowerBound: An integer with a lower bound for the        */
/*                 minimum weight of C                          */
/*                                                              */
/****************************************************************/
UpdateMinimumWeightLowerBound := procedure(~C, lowerBound)
    if (lowerBound gt C`MinimumWeightLowerBound) and (lowerBound le C`MinimumWeightUpperBound) then
        C`MinimumWeightLowerBound := lowerBound;
    end if;
    if C`MinimumWeightLowerBound eq C`MinimumWeightUpperBound then
        C`MinimumWeight := C`MinimumWeightLowerBound;
    end if;
end procedure;

/****************************************************************/
/*                                                              */
/* Procedure name: UpdateMinimumWeightUpperBound                */
/* Parameters: ~C, upperBound                                   */
/* Procedure description: Given a binary linear code C and an   */
/*   integer with an upper bound for its minimum weight,        */
/*   update the attribute MinimumWeightUpperBound.              */  
/* Input parameters description:                                */
/*   - C: A binary linear code                                  */
/*   - upperBound: An integer with an upper bound for the       */
/*                 minimum weight of C                          */
/*                                                              */
/****************************************************************/
UpdateMinimumWeightUpperBound := procedure(~C, upperBound)
    if (upperBound lt C`MinimumWeightUpperBound) and 
       (upperBound ge C`MinimumWeightUpperBound) then
        C`MinimumWeightUpperBound := upperBound;
    end if;
    if C`MinimumWeightLowerBound eq C`MinimumWeightUpperBound then
        C`MinimumWeight := C`MinimumWeightLowerBound;
    end if;
end procedure;

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
intrinsic SpanZ2Code(C::Z2Z4Code) -> Z2Z4Code, CodeLinFld
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma, delta; kappa), 
return S_C=Phi^(-1)(Sbin) as a Z2Z4-additive code, and Sbin=<Cbin>, that is the 
linear span of Cbin, as a binary linear code of length alpha + 2*beta, 
where Cbin=Phi(C) and Phi is the Gray map.
}
    SZ4 := SpanZ2CodeZ4(C`Code);
    SC := Z2Z4AdditiveCode(SZ4, C`Alpha);
    _, SZ2 := HasLinearGrayMapImage(SC);
    
    if (C eq SC) then
        UpdateMinimumLeeWeightLowerBound(~SC, C`MinimumLeeWeightLowerBound);
        UpdateMinimumLeeWeightUpperBound(~SC, C`MinimumLeeWeightUpperBound);
        UpdateMinimumWeightLowerBound(~SZ2, C`MinimumLeeWeightLowerBound); 
        UpdateMinimumWeightUpperBound(~SZ2, C`MinimumLeeWeightUpperBound);
    else
        UpdateMinimumLeeWeightUpperBound(~SC, C`MinimumLeeWeightUpperBound);
        UpdateMinimumWeightUpperBound(~SZ2, C`MinimumLeeWeightUpperBound);
    end if;
    
    return SC, SZ2; 
end intrinsic;

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
intrinsic DimensionOfSpanZ2(C::Z2Z4Code) -> RngIntElt
{
Given a Z2Z4-additive code C, return the dimension of the linear span of Cbin,
that is, the dimension of <Cbin>, where Cbin=Phi(C) and Phi is the Gray map.
}
    return DimensionOfSpanZ2(C`Code);
end intrinsic;

/************************************************************/
/*                                                          */
/* Function name: RankZ2                                    */
/* Parameters: C                                            */
/* Function description: Given a Z2Z4-additive code C,      */
/*   return the dimension of the linear span of Cbin,       */
/*   that is, the dimension of <Cbin>, where Cbin=Phi(C)    */
/*   and Phi is the Gray map.                               */
/* Input parameters description:                            */
/*   - C: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - The dimension of the linear span of Cbin      */
/*                                                          */
/* Signature: (<Z2Z4Code> C) -> RngIntElt                   */
/*                                                          */
/************************************************************/
intrinsic RankZ2(C::Z2Z4Code) -> RngIntElt
{
Given a Z2Z4-additive code C, return the dimension of the linear span of Cbin,
that is, the dimension of <Cbin>, where Cbin=Phi(C) and Phi is the Gray map.
}
    return DimensionOfSpanZ2(C);
end intrinsic;

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
intrinsic KernelZ2Code(C::Z2Z4Code) -> Z2Z4Code, CodeLinFld
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma, delta; kappa),
return its kernel K_C as a Z2Z4-additive subcode of C, and Kbin=Phi(K_C) 
as a binary linear subcode of Cbin of length alpha + 2 beta, where Cbin=Phi(C) 
and Phi is the Gray map.

The kernel K_C contains the codewords v such that 2v*u in C for all u in C,
where * denotes the component-wise product. Equivalently, the kernel 
Kbin=Phi(K_C) contains the codewords c in Cbin such
that c+Cbin=Cbin, where Cbin=Phi(C) and Phi is the Gray map.
}
    KZ4 := KernelZ2CodeZ4(C`Code);
    KC := Z2Z4AdditiveCode(KZ4, C`Alpha);
    _, KZ2 := HasLinearGrayMapImage(KC);
    
    UpdateMinimumLeeWeightLowerBound(~KC, C`MinimumLeeWeightLowerBound);
    UpdateMinimumWeightLowerBound(~KZ2, C`MinimumLeeWeightLowerBound); 
  
    return KC, KZ2;  
end intrinsic;

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
intrinsic KernelCosetRepresentatives(C::Z2Z4Code) -> SeqEnum, SeqEnum
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma, delta; kappa),
return the coset representatives [c_1,..., c_t] as a sequence of codewords of C, 
such that C = KC U (U_(i=1)^t (KC + c_i)), where KC is the kernel of C as a 
Z2Z4-additive subcode. It also returns the coset representatives of the 
corresponding binary code Cbin = phi(C) as a sequence of binary codewords 
[phi(c_1),..., phi(c_t)], such that Cbin = Kbin U (U_(i=1)^t (Kbin + phi(c_i)))
where Kbin = phi(KC) and phi is the Gray map.
}
    grayMap := GrayMap(C);
    V := VectorSpace(GF(2), BinaryLength(C));
    leadersZ2Z4 := KernelCosetRepresentatives(C`Code);
    leadersZ2 := [V!grayMap(leader) : leader in leadersZ2Z4];
    
    return leadersZ2Z4, leadersZ2;     
end intrinsic;

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
intrinsic DimensionOfKernelZ2(C::Z2Z4Code) -> RngIntElt
{
Given a Z2Z4-additive code C, return the dimension of the Gray map image of its 
Z2Z4-additive kernel K_C, that is the dimension of Kbin = Phi(K_C), where Phi 
is the Gray map. Note that Kbin is always a binary linear code.
}
    _, KZ2 := KernelZ2CodeZ4(C`Code);

    return Dimension(KZ2);       
end intrinsic;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////        FROM LINEAR TO Z2Z4-ADDITIVE CODES                       ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4AdditiveFromBinaryLinearCode              */
/* Parameters: C                                                */
/* Function description: Given a binary linear code C of length */
/*   n, return the same code as a Z2Z4-additive code, so with   */ 
/*   alpha = n and beta = 0.                                    */
/* Input parameters description:                                */
/*   - C: A binary linear code C                                */
/* Output parameters description:                               */
/*   - Binary linear code C as a Z2Z-additive code              */
/*                                                              */
/* Signature: (<CodeLinFld> C) -> Z2Z4Code                      */
/*                                                              */
/****************************************************************/
/*intrinsic Z2Z4AdditiveFromBinaryLinearCode(C::CodeLinFld) -> Z2Z4Code
{
Given a binary linear code C of length n, return the same code as a
Z2Z4-additive code, so with alpha = n and beta = 0.
}
    require (over_Z2(C) or over_F2(C)): "Code C is not a code over GF(2) or Z2";

    return Z2Z4AdditiveCode(C);
end intrinsic;*/

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4AdditiveFromQuaternaryLinearCode          */
/* Parameters: C                                                */
/* Function description: Given a quaternary linear code C of    */
/*   length n, return the same code as a Z2Z4-additive code,    */
/*   so with alpha = 0 and beta = n..                           */
/* Input parameters description:                                */
/*   - C: A quaternary linear code C                            */
/* Output parameters description:                               */
/*   - Quaternary linear code C as a Z2Z-additive code          */
/*                                                              */
/* Signature: (<CodeLinRng> C) -> Z2Z4Code                      */
/*                                                              */
/****************************************************************/
/*intrinsic Z2Z4AdditiveFromQuaternaryLinearCode(C::CodeLinRng) -> Z2Z4Code
{
Given a quaternary linear code C of length n, return the same code as a
Z2Z4-additive code, so with alpha = 0 and beta = n.
}
    require over_Z4(C): "Code C is not a code over Z4";

    return Z2Z4AdditiveCode(C);
end intrinsic;*/

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////                SUBCODES Cx AND Cy                               ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/****************************************************************/
/*                                                              */
/* Function name: LinearBinaryCode                              */
/* Parameters: C                                                */
/* Function description: Given a Z2Z4-additive code C of type   */
/*   (alpha, beta; gamma, delta; kappa), return the binary      */
/*   linear code CX of length alpha which is the punctured code */
/*   of C by deleting the coordinates outside X, where X is the */
/*   set of the first alpha coordinates.                        */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - The binary linear code CX of length alpha                */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> CodeLinFld                      */
/*                                                              */
/****************************************************************/
intrinsic LinearBinaryCode(C::Z2Z4Code) -> CodeLinFld
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma, delta; kappa), 
return the binary linear code CX of length alpha which is the punctured code 
of C by deleting the coordinates outside X, where X is the set of the first 
alpha coordinates.
}
    require (C`Alpha gt 0): "Code C has not binary part";

    M := GeneratorMatrix(C`Code);
    Z2Z4ChangeMatrixZ4toZ2(~M, C`Alpha);
    CX := LinearCode(ChangeRing(ChangeRing(ColumnSubmatrix(M, 1, C`Alpha), Z2), F2));
    if C`Length eq C`Alpha then 
        UpdateMinimumWeightLowerBound(~CX, C`MinimumLeeWeightLowerBound); 
        UpdateMinimumWeightUpperBound(~CX, C`MinimumLeeWeightUpperBound);
    end if;
    
    return CX;
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: LinearQuaternaryCode                          */
/* Parameters: C                                                */
/* Function description: Given a Z2Z4-additive code C of type   */
/*   (alpha, beta; gamma, delta; kappa), return the quaternary  */
/*   linear code CY of length beta which is the punctured code  */
/*   of C by deleting the coordinates outside Y , where Y is    */
/*   the set of the last beta coordinates.                      */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - The quaternary linear code CY of length beta             */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> CodeLinRng                      */
/*                                                              */
/****************************************************************/
intrinsic LinearQuaternaryCode(C::Z2Z4Code) -> CodeLinRng
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma, delta; kappa), return 
the quaternary linear code CY of length beta which is the punctured code of C 
by deleting the coordinates outside Y , where Y is the set of the last beta 
coordinates.
}
   require ((C`Length-C`Alpha) gt 0): "Code C has not quaternary part";
            
    M := GeneratorMatrix(C`Code);
    CY := LinearCode(ColumnSubmatrix(M, C`Alpha+1, Ncols(M)-C`Alpha));
    if IsZero(C`Alpha) and (assigned C`MinimumLeeWeight) then 
        CY`MinimumLeeWeight := C`MinimumLeeWeight;
    end if;

    return CY;
end intrinsic;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////            CONSTRUCTION OF CODEWORD                             ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/****************************************************************/
/*                                                              */
/* Function name: Random                                        */
/* Parameters: C                                                */
/* Function description: A random codeword of the Z2Z4-additive */
/*   code C, which is a vector in Z4^(alpha+beta) where the     */
/*   ones in the first alpha coordinates are represented by twos*/
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - A random codeword of the Z2Z4-additive code C            */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> ModTupRngElt                    */
/*                                                              */
/****************************************************************/
intrinsic Random(C::Z2Z4Code) -> ModTupRngElt
{
A random codeword of the Z2Z4-additive code C, which is a vector in Z4^(alpha+beta) 
where the ones in the first alpha coordinates are represented by twos. 
}
    return Random(C`Code);
end intrinsic;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////            OPERATION ON CODEWORDS AND VECTORS                   ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/****************************************************************/
/*                                                              */
/* Function name: '+'                                           */
/* Parameters: u, v                                             */
/* Function description: Sum of the codewords u and v, where u  */
/*   and v belong to the same Z2Z4-additive code C of type      */
/*   (alpha, beta; gamma, delta; kappa). The vectors u and v    */
/*   are represented as tuples in the cartesian product set     */
/*   Z2^alpha x Z4^beta.                                        */
/* Input parameters description:                                */
/*   - u, v: Two tuples in Z2^alpha x Z4^beta                   */
/* Output parameters description:                               */
/*   - The sum u + v as a tuple in Z2^alpha x Z4^beta           */
/*                                                              */
/* Signature: (<Tup> u, <Tup> v) -> <Tup>                       */
/*                                                              */
/****************************************************************/
intrinsic '+'(u::Tup, v::Tup) -> Tup
{
Sum of the codewords u and v, where u and v belong to the same Z2Z4-additive code 
C of type (alpha, beta; gamma, delta; kappa). The vectors u and v are represented as 
tuples in the cartesian product set Z2^alpha x Z4^beta.
}
    require #u eq 2: "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require #v eq 2: "Argument 2 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (Type(u[1]) cmpeq ModTupRngElt) and (Type(u[2]) cmpeq ModTupRngElt): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (Type(v[1]) cmpeq ModTupRngElt) and (Type(v[2]) cmpeq ModTupRngElt): 
            "Argument 2 is not an element of the cartesian product Z2^alpha x Z4^beta";    
    require (BaseRing(u[2]) cmpeq Z4 and BaseRing(u[1]) cmpeq Z2): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (BaseRing(v[2]) cmpeq Z4 and BaseRing(v[1]) cmpeq Z2): 
            "Argument 2 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (Ncols(u[1]) eq Ncols(v[1]) and Ncols(u[2]) eq Ncols(v[2])): 
            "Arguments 1 and 2 must have the same length";

    alpha := Ncols(u[1]);   
    beta := Ncols(u[2]);
    R := CartesianProduct(RSpace(Z2, alpha), RSpace(Z4, beta));
    binSum := u[1] + v[1];
    quaSum := u[2] + v[2];

    return R ! <binSum, quaSum>;  
end intrinsic;
 
/****************************************************************/
/*                                                              */
/* Function name: '-'                                           */
/* Parameters: u                                                */
/* Function description: Additive inverse of the codeword u     */
/*   belonging to the Z2Z4-additive code C of type (alpha, beta;*/
/*   gamma, delta; kappa). The vector u is represented as a     */
/*   tuple in the cartesian product set Z2^alpha x Z4^beta.     */
/* Input parameters description:                                */
/*   - u: A tuple in Z2^alpha x Z4^beta                         */
/* Output parameters description:                               */
/*   - The additive inverse -u as a tuple in Z2^alpha x Z4^beta */
/*                                                              */
/* Signature: (<Tup> u) -> <Tup>                                */
/*                                                              */
/****************************************************************/
intrinsic '-'(u::Tup) -> Tup
{
Additive inverse of the codeword u belonging to the Z2Z4-additive code C of type 
(alpha, beta; gamma, delta; kappa). The vector u is represented as a tuple in 
the cartesian product set Z2^alpha x Z4^beta.
}
    require #u eq 2: "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (Type(u[1]) cmpeq ModTupRngElt) and (Type(u[2]) cmpeq ModTupRngElt): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";    
    require (BaseRing(u[2]) cmpeq Z4 and BaseRing(u[1]) cmpeq Z2): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";

    alpha := Ncols(u[1]);   
    beta := Ncols(u[2]);
    R := CartesianProduct(RSpace(Z2, alpha), RSpace(Z4, beta));
    binInverse := - u[1];
    quaInverse := - u[2];

    return R ! <binInverse, quaInverse>;  
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: '-'                                           */
/* Parameters: u, v                                             */
/* Function description: Difference of the codewords u and v,   */
/*   where u and v belong to the same Z2Z4-additive code C of   */
/*   type (alpha, beta; gamma, delta; kappa). The vectors u and */
/*   v are represented as tuples in the cartesian product set   */
/*   Z2^alpha x Z4^beta.                                        */
/* Input parameters description:                                */
/*   - u, v: Two tuples in Z2^alpha x Z4^beta                   */
/* Output parameters description:                               */
/*   - The difference u - v as a tuple in Z2^alpha x Z4^beta    */
/*                                                              */
/* Signature: (<Tup> u, <Tup> v) -> <Tup>                       */
/*                                                              */
/****************************************************************/
intrinsic '-'(u::Tup, v::Tup) -> Tup
{
Difference of the codewords u and v, where u and v belong to the same Z2Z4-additive 
code C of type (alpha, beta; gamma, delta; kappa). The vectors u and v are represented
as tuples in the cartesian product set Z2^alpha x Z4^beta.
}
    require #u eq 2: "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require #v eq 2: "Argument 2 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (Type(u[1]) cmpeq ModTupRngElt) and (Type(u[2]) cmpeq ModTupRngElt): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (Type(v[1]) cmpeq ModTupRngElt) and (Type(v[2]) cmpeq ModTupRngElt): 
            "Argument 2 is not an element of the cartesian product Z2^alpha x Z4^beta";    
    require (BaseRing(u[2]) cmpeq Z4 and BaseRing(u[1]) cmpeq Z2): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (BaseRing(v[2]) cmpeq Z4 and BaseRing(v[1]) cmpeq Z2): 
            "Argument 2 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (Ncols(u[1]) eq Ncols(v[1]) and Ncols(u[2]) eq Ncols(v[2])): 
            "Arguments 1 and 2 must have the same length";

    alpha := Ncols(u[1]);   
    beta := Ncols(u[2]);
    R := CartesianProduct(RSpace(Z2, alpha), RSpace(Z4, beta));
    binDif := u[1] - v[1];
    quaDif := u[2] - v[2];

    return R ! <binDif, quaDif>;  
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: '*'                                           */
/* Parameters: a, u                                             */
/* Function description: Given an element a belonging to Z4,    */
/*   and a codeword u belonging to the Z2Z4-additive code C of  */
/*   type (alpha, beta; gamma, delta; kappa) return the codeword*/
/*   a * u. The vector u is represented as a tuple in the       */
/*   cartesian product set Z2^alpha x Z4^beta.                  */
/* Input parameters description:                                */
/*   - a : a Z4 element                                         */
/*   - u : a tuple in Z2^alpha x Z4^beta                        */
/* Output parameters description:                               */
/*   - The product a * u as a tuple in Z2^alpha x Z4^beta       */
/*                                                              */
/* Signature: (<RngIntElt> a, <Tup> u) -> <Tup>                 */
/*                                                              */
/****************************************************************/
intrinsic '*'(a::., u::Tup) -> Tup
{
Given an element a belonging to Z4, and a codeword u belonging to the Z2Z4-additive 
code C of type (alpha, beta; gamma, delta; kappa) return the codeword a * u. The 
vector u is represented as a tuple in the cartesian product set Z2^alpha x Z4^beta.
}
    
    require #u eq 2: "Argument 2 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (Type(u[1]) cmpeq ModTupRngElt) and (Type(u[2]) cmpeq ModTupRngElt): 
            "Argument 2 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (BaseRing(u[2]) cmpeq Z4 and BaseRing(u[1]) cmpeq Z2): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    isCoerc, a := IsCoercible(Z4, a);
    require isCoerc: "Argument 1 is not an element of Z4";

    alpha := Ncols(u[1]);   
    beta := Ncols(u[2]);
    R := CartesianProduct(RSpace(Z2, alpha), RSpace(Z4, beta));
    binProd := a * u[1];
    quaProd := a * u[2];

    return R ! <binProd, quaProd>;  
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4InnerProduct                              */
/* Parameters: u, v, alpha                                      */
/* Function description: The inner product <u,v> in             */
/*   Z2^alpha x Z4^beta. The vectors u and v are represented as */
/*   elements in Z4^(alpha+beta) by replacing the ones in the   */
/*   first alpha coordinates by twos.                           */
/*space of the Z2Z4-additive code C.                            */
/* Input parameters description:                                */
/*   - u, v: Two vectors in Z4^n                                */
/*   - alpha: The length of the binary part                     */
/* Output parameters description:                               */
/*   - The inner product <u,v> in Z2^alpha x Z4^beta            */
/*                                                              */
/* Signature: (<ModTupRngElt> v, <ModTupRngElt> u,              */
/*                              <RngIntElt> alpha) -> RngIntElt */
/*                                                              */
/****************************************************************/
intrinsic Z2Z4InnerProduct(v::ModTupRngElt, u::ModTupRngElt, alpha::RngIntElt) 
                            -> RngIntElt
{
The inner product <u,v> in Z2^alpha x Z4^beta. The vectors u and v are represented 
as elements in Z4^(alpha+beta) by replacing the ones in the first alpha coordinates by twos.
}
    require (BaseRing(v) cmpeq Z4): "Argument 1 is not a vector over Z4";
    require (BaseRing(u) cmpeq Z4): "Argument 2 is not a vector over Z4";
    require (Ncols(v) eq Ncols(u)): 
            "Arguments 1 and 2 must have the same length";
    requirerange alpha, 0, Degree(u);
    require (IsZero(2*ColumnSubmatrix(v, alpha))) and
            (IsZero(2*ColumnSubmatrix(u, alpha))):
            "First", alpha, "coordinates must be in {0,2}";
    
    B := [Z2Z4AlphaFromZ4toZ2(v[i]) * Z2Z4AlphaFromZ4toZ2(u[i]) : i in [1..alpha]];
    Q := [v[i] * u[i] : i in [alpha+1..Ncols(v)]];
    sumB := (#B eq 0) select 0 else &+B;
    sumQ := (#Q eq 0) select 0 else &+Q;

    return 2 * sumB + sumQ;
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: '.'                                           */
/* Parameters: u, v                                             */
/* Function description: The inner product <u,v> in             */
/*   Z2^alpha x Z4^beta. The vectors u and v are represented as */
/*   tuples in the cartesian product set Z2^alpha x Z4^beta.    */
/* Input parameters description:                                */
/*   - u, v: Two tuples in Z2^alpha x Z4^beta                   */
/* Output parameters description:                               */
/*   - The inner product <u,v> in Z2^alpha x Z4^beta            */
/*                                                              */
/* Signature: (<Tup> v, <Tup> u) -> RngIntElt                   */
/*                                                              */
/****************************************************************/
intrinsic '.'(v::Tup, u::Tup) -> RngIntElt
{
The inner product <u,v> in Z2^alpha x Z4^beta. The vectors u and v are represented 
as tuples in the cartesian product set Z2^alpha x Z4^beta.
}
    require #v eq 2: "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require #u eq 2: "Argument 2 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (Type(v[1]) cmpeq ModTupRngElt) and (Type(v[2]) cmpeq ModTupRngElt): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (Type(u[1]) cmpeq ModTupRngElt) and (Type(u[2]) cmpeq ModTupRngElt): 
            "Argument 2 is not an element of the cartesian product Z2^alpha x Z4^beta";    
    require (BaseRing(v[2]) cmpeq Z4 and BaseRing(v[1]) cmpeq Z2): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (BaseRing(u[2]) cmpeq Z4 and BaseRing(u[1]) cmpeq Z2): 
            "Argument 2 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (Ncols(v[1]) eq Ncols(u[1]) and Ncols(v[2]) eq Ncols(u[2])): 
            "Arguments 1 and 2 must have the same length";
    
    alpha := Ncols(v[1]);   
    beta := Ncols(v[2]);
    binInnerProduct := v[1] * u[1];
    quaInnerProduct := v[2] * u[2];
    sumBin := &+[binInnerProduct[i] : i in [1..alpha]];
    sumQua := &+[quaInnerProduct[i] : i in [1..beta]];

    return 2 * (Z4!sumBin) + sumQua;
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4InnerProduct                              */
/* Parameters: u, v                                             */
/* Function description: The inner product <u,v> in             */
/*   Z2^alpha x Z4^beta. The vectors u and v are represented as */
/*   tuples in the cartesian product set Z2^alpha x Z4^beta.    */
/* Input parameters description:                                */
/*   - u, v: Two tuples in Z2^alpha x Z4^beta                   */
/* Output parameters description:                               */
/*   - The inner product <u,v> in Z2^alpha x Z4^beta            */
/*                                                              */
/* Signature: (<Tup> v, <Tup> u) -> RngIntElt                   */
/*                                                              */
/****************************************************************/
intrinsic Z2Z4InnerProduct(v::Tup, u::Tup) -> RngIntElt
{
The inner product <u,v> in Z2^alpha x Z4^beta. The vectors u and v are represented 
as tuples in the cartesian product set Z2^alpha x Z4^beta.
}
    require #v eq 2: "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require #u eq 2: "Argument 2 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (Type(v[1]) cmpeq ModTupRngElt) and (Type(v[2]) cmpeq ModTupRngElt): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (Type(u[1]) cmpeq ModTupRngElt) and (Type(u[2]) cmpeq ModTupRngElt): 
            "Argument 2 is not an element of the cartesian product Z2^alpha x Z4^beta";    
    require (BaseRing(v[2]) cmpeq Z4 and BaseRing(v[1]) cmpeq Z2): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (BaseRing(u[2]) cmpeq Z4 and BaseRing(u[1]) cmpeq Z2): 
            "Argument 2 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (Ncols(v[1]) eq Ncols(u[1]) and Ncols(v[2]) eq Ncols(u[2])): 
            "Arguments 1 and 2 must have the same length";

    return v . u;
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4StarProduct                               */
/* Parameters: u, v, alpha                                      */
/* Function description: The component-wise product of u and v  */
/*   in Z2^alpha x Z4^beta. The vectors u and v are represented */
/*   as elements in Z4^(alpha+beta) by replacing the ones in the*/
/*   first alpha coordinates by twos.                           */
/* Input parameters description:                                */
/*   - u, v: Two vectors in Z4^n                                */
/*   - alpha: The length of the binary part                     */
/* Output parameters description:                               */
/*   - The vector over Z4 having the component-wise product of  */
/*     vectors u and v in Z2^alpha x Z4^beta                    */
/*                                                              */
/* Signature: (<ModTupRngElt> v, <ModTupRngElt> u,              */
/*                           <RngIntElt> alpha) -> ModTupRngElt */
/*                                                              */
/****************************************************************/
intrinsic Z2Z4StarProduct(v::ModTupRngElt, u::ModTupRngElt, alpha::RngIntElt) 
                            -> ModTupRngElt
{
The component-wise product of u and v in Z2^alpha x Z4^beta. The vectors u and v are 
represented as elements in Z4^(alpha+beta) by replacing the ones in the first alpha 
coordinates by twos.
}
    require (BaseRing(v) cmpeq Z4): "Argument 1 is not a vector over Z4";
    require (BaseRing(u) cmpeq Z4): "Argument 2 is not a vector over Z4";
    require (Ncols(v) eq Ncols(u)): 
            "Arguments 1 and 2 must have the same length";
    requirerange alpha, 0, Degree(u);
    require (IsZero(2*ColumnSubmatrix(v, alpha))) and
            (IsZero(2*ColumnSubmatrix(u, alpha))):
            "First", alpha, "coordinates must be in {0,2}";
    
    S := [Z2Z4AlphaFromZ2toZ4(Z2Z4AlphaFromZ4toZ2(v[i])*Z2Z4AlphaFromZ4toZ2(u[i])) 
                                                                : i in [1..alpha]];
    S cat:= [v[i]*u[i] : i in [alpha+1..Ncols(v)]];

    return Vector(Z4, S);
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: '*'                                           */
/* Parameters: u, v                                             */
/* Function description: The component-wise product of u and v  */
/*   in Z2^alpha x Z4^beta. The vectors u and v are represented */
/*   as tuples in the cartesian product set Z2^alpha x Z4^beta. */
/* Input parameters description:                                */
/*   - u, v: Two tuples in Z2^alpha x Z4^beta                   */
/* Output parameters description:                               */
/*   - The tuple having the component-wise product of vectors   */
/*     u and v in Z2^alpha x Z4^beta                            */
/*                                                              */
/* Signature: (<Tup> v, <Tup> u) -> Tup                         */
/*                                                              */
/****************************************************************/
intrinsic '*'(v::Tup, u::Tup) -> Tup 
{
The component-wise product of u and v in Z2^alpha x Z4^beta. The vectors u and v are
represented as tuples in the cartesian product set Z2^alpha x Z4^beta.  
}
    require #v eq 2: "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require #u eq 2: "Argument 2 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (Type(v[1]) cmpeq ModTupRngElt) and (Type(v[2]) cmpeq ModTupRngElt): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (Type(u[1]) cmpeq ModTupRngElt) and (Type(u[2]) cmpeq ModTupRngElt): 
            "Argument 2 is not an element of the cartesian product Z2^alpha x Z4^beta";    
    require (BaseRing(v[2]) cmpeq Z4 and BaseRing(v[1]) cmpeq Z2): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (BaseRing(u[2]) cmpeq Z4 and BaseRing(u[1]) cmpeq Z2): 
            "Argument 2 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (Ncols(v[1]) eq Ncols(u[1]) and Ncols(v[2]) eq Ncols(u[2])): 
            "Arguments 1 and 2 must have the same length";
    
    alpha := Ncols(v[1]);   
    beta := Ncols(v[2]);
    R := CartesianProduct(RSpace(Z2, alpha), RSpace(Z4, beta));
    binStarProduct := v[1] * u[1];
    quaStarProduct := v[2] * u[2];
    
    return R ! <binStarProduct, quaStarProduct>;  
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4StarProduct                               */
/* Parameters: u, v                                             */
/* Function description: The component-wise product of u and v  */
/*   in Z2^alpha x Z4^beta. The vectors u and v are represented */
/*   as tuples in the cartesian product set Z2^alpha x Z4^beta. */
/* Input parameters description:                                */
/*   - u, v: Two tuples in Z2^alpha x Z4^beta                   */
/* Output parameters description:                               */
/*   - The tuple having the component-wise product of vectors   */
/*     u and v in Z2^alpha x Z4^beta                            */
/*                                                              */
/* Signature: (<Tup> v, <Tup> u) -> Tup                         */
/*                                                              */
/****************************************************************/
intrinsic Z2Z4StarProduct(v::Tup, u::Tup) -> Tup
{
The component-wise product of u and v in Z2^alpha x Z4^beta. The vectors u and v are
represented as tuples in the cartesian product set Z2^alpha x Z4^beta.  
}
    require #v eq 2: "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require #u eq 2: "Argument 2 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (Type(v[1]) cmpeq ModTupRngElt) and (Type(v[2]) cmpeq ModTupRngElt): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (Type(u[1]) cmpeq ModTupRngElt) and (Type(u[2]) cmpeq ModTupRngElt): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";    
    require (BaseRing(v[2]) cmpeq Z4 and BaseRing(v[1]) cmpeq Z2): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (BaseRing(u[2]) cmpeq Z4 and BaseRing(u[1]) cmpeq Z2): 
            "Argument 2 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (Ncols(v[1]) eq Ncols(u[1]) and Ncols(v[2]) eq Ncols(u[2])): 
            "Arguments 1 and 2 must have the same length";
    
    return v * u;
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Normalize                                     */
/* Parameters: u                                                */
/* Function description: Given a codeword u belonging to the    */
/*   Z2Z4-additive code C, which is represented as a subspace of*/
/*   Z4^alpha+beta, return the normalization of u, which        */
/*   is the unique vector v such that v = a�u for some scalar a */
/*   in Z4 such that the first non-zero entry of v is the       */
/*   canonical associate in Z4 of the first non-zero entry of   */
/*   u (v is zero if u is zero). The codeword u is given as a   */
/*   tuple in the cartesian product set Z2^alpha x Z^beta.      */
/* Input parameters description:                                */
/*   - u: A vector in Z2^alpha x Z4^beta                        */
/* Output parameters description:                               */
/*   - The normalization of u                                   */
/*                                                              */
/* Signature: (<Tup> u)-> ModTupRngElt                          */
/*                                                              */
/****************************************************************/
intrinsic Normalize(u::Tup) -> ModTupRngElt
{
Given a codeword u belonging to the Z2Z4-additive code C, which is
represented as a codeword of Z2^alpha x Z4^beta, return the normalization of u, 
which is the unique vector v such that v = a�u for some scalar a in Z4 such
that the first non-zero entry of v is the canonical associate in Z4 of the
first non-zero entry of u (v is zero if u is zero). The codeword u is given 
as a tuple in the cartesian product set Z2^alpha x Z^beta.
}
    require #u eq 2: "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require Type(u[1]) cmpeq ModTupRngElt and Type(u[2]) cmpeq ModTupRngElt: 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (BaseRing(u[2]) cmpeq Z4 and BaseRing(u[1]) cmpeq Z2): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    
    alpha := Ncols(u[1]);
    beta := Ncols(u[2]);
    n := alpha + beta;
    require n gt 0 : "Argument 1 must be a vector of length greater than 0";
    u4 := FromVectorZ2Z4toZ4(u);

    return Normalize(u4);

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Support                                       */
/* Parameters: u                                                */
/* Function description: Given a codeword u belonging to the    */
/*   Z2Z4-additive code C, which is represented as a subspace   */
/*   of Z4^alpha+beta, return its support as a subset of the    */
/*   integer set {1,...,alpha+beta}. The support of u consists  */
/*   of the coordinates at which u has non-zero entries. The    */
/*   codeword u is given as a tuple in the cartesian product    */
/*   set Z2^alpha x Z^beta.                                     */
/* Input parameters description:                                */
/*   - u: A vector in Z2^alpha x Z4^beta                        */
/* Output parameters description:                               */
/*   - The support of u                                         */
/*                                                              */
/* Signature: (<Tup> u)-> Set                                   */
/*                                                              */
/****************************************************************/
intrinsic Support(u::Tup) -> Set
{
Given a codeword u belonging to the Z2Z4-additive code C, which is represented 
as a subspace of Z4^alpha+beta, return its support as a subset of the 
integer set (1,...,alpha+beta). The support of u consists of the coordinates 
at which u has non-zero entries. The codeword u is given as a tuple in the 
cartesian product set Z2^alpha x Z^beta.
}
    require #u eq 2: "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require Type(u[1]) cmpeq ModTupRngElt and Type(u[2]) cmpeq ModTupRngElt: 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (BaseRing(u[2]) cmpeq Z4 and BaseRing(u[1]) cmpeq Z2): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    
    alpha := Ncols(u[1]);
    beta := Ncols(u[2]);
    n := alpha + beta;
    require n gt 0 : "Argument 1 must be a vector of length greater than 0";
    u4 := FromVectorZ2Z4toZ4(u);

    return Support(u4);

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Coordinates                                   */
/* Parameters: C, u                                             */
/* Function description: Given a Z2Z4-additive code C and a     */
/*   codeword u of C, return the coordinates of u with respect  */
/*   to C. The coordinates of u are returned as a sequence      */
/*   Q = [a1,..., ak] of elements from Z4 so that               */
/*   u = a1*C.1 +...+ ak*C.k. The codeword u is given as a      */
/*   tuple in the cartesian product set Z2^alpha x Z^beta.      */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/*   - u: A vector in Z4^n                                      */
/* Output parameters description:                               */
/*   - The coordinates of u with respect to C                   */
/*                                                              */
/* Signature: (<Z2Z4Code> C, <ModTupRngElt> u )-> SeqEnum        */
/*                                                              */
/****************************************************************/
intrinsic Coordinates(C::Z2Z4Code, u::ModTupRngElt) -> SeqEnum
{
Given a Z2Z4-additive code C and a codeword u of C, return the coordinates of 
u with respect to C. The coordinates of u are returned as a sequence 
Q = [a1,..., ak] of elements from Z4 so that u = a1*C.1 +...+ ak*C.k. 
The codeword u is given as a tuple in the cartesian product set Z2^alpha x Z^beta.
}
    require (BaseRing(u) cmpeq Z4): "Argument 2 is not a vector over Z4";
    require u in C`Code: "Argument 2 is not a codeword in C";

    return Coordinates(C`Code, u);
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Coordinates                                   */
/* Parameters: C, u                                             */
/* Function description: Given a Z2Z4-additive code C and a     */
/*   codeword u of C, return the coordinates of u with respect  */
/*   to C. The coordinates of u are returned as a sequence      */
/*   Q = [a1,..., ak] of elements from Z4 so that               */
/*   u = a1*C.1 +...+ ak*C.k.                                   */
/*   The codeword u is given as a tuple in the cartesian        */
/*   product set Z2^alpha x Z^beta                              */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/*   - u: A vector in Z2^alpha x Z4^beta                        */
/* Output parameters description:                               */
/*   - The coordinates of u with respect to C                   */
/*                                                              */
/* Signature: (<Z2Z4Code> C, <Tup> u )-> SeqEnum                 */
/*                                                              */
/****************************************************************/
intrinsic Coordinates(C::Z2Z4Code, u::Tup) -> SeqEnum
{
Given a Z2Z4-additive code C and a codeword u of C, return the coordinates of 
u with respect to C. The coordinates of u are returned as a sequence 
Q = [a1,..., ak] of elements from Z4 so that u = a1*C.1 +...+ ak*C.k.
The codeword u is given as a tuple in the cartesian product set Z2^alpha x Z^beta.
}
    require #u eq 2: "Argument 2 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require Type(u[1]) cmpeq ModTupRngElt and Type(u[2]) cmpeq ModTupRngElt: 
            "Argument 2 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (BaseRing(u[2]) cmpeq Z4 and BaseRing(u[1]) cmpeq Z2): 
            "Argument 2 is not an element of the cartesian product Z2^alpha x Z4^beta";
    
    alpha := Ncols(u[1]);
    beta := Ncols(u[2]);
    n := alpha + beta;
    require n gt 0 : "Argument 2 must be a vector of length greater than 0";

    u4 := FromVectorZ2Z4toZ4(u);

    require u4 in C`Code: "Argument 2 is not a codeword in C";

    return Coordinates(C`Code, u4);
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Rotate                                        */
/* Parameters: u, k                                             */
/* Function description: Given a vector u belonging to the      */
/*   cartesian product set Z2^alpha x Z4^beta, return the vector*/
/*   obtained from u by cyclically shifting its first alpha     */
/*   components and last beta components, separetely, to the    */
/*   right by k coordinate positions.                           */
/* Input parameters description:                                */
/*   - u: A vector in Z2^alpha x Z4^beta                        */
/*   - k: an integer                                            */
/* Output parameters description:                               */
/*   - the vector u cyclically shifting k components            */
/*                                                              */
/* Signature: (<Tup> u, <RngIntEly > k)-> <Tup>                 */
/*                                                              */
/****************************************************************/
intrinsic Rotate(u::Tup, k:RngIntElt) -> Tup
{
Given a vector u belonging to the cartesian product set Z2^alpha x Z4^beta, 
return the vector obtained from u by cyclically shifting its first alpha 
components and last beta components, separetely, to the right by k 
coordinate positions.
}
    require #u eq 2: "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require Type(u[1]) cmpeq ModTupRngElt and Type(u[2]) cmpeq ModTupRngElt: 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (BaseRing(u[2]) cmpeq Z4 and BaseRing(u[1]) cmpeq Z2): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";

    u1 := u[1];
    u2 := u[2];
    Rotate(~u1, k);
    Rotate(~u2, k);

    return <u1, u2>;

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Rotate                                        */
/* Parameters: u, k                                             */
/* Function description: Given a vector u belonging to the      */
/*   cartesian product set Z2^alpha x Z4^beta, return the vector*/
/*   obtained from u by cyclically shifting its first alpha     */
/*   components and last beta components, separetely, to the    */
/*   right by k coordinate positions.                           */
/* Input parameters description:                                */
/*   - u: A vector in Z2^alpha x Z4^beta                        */
/*   - k: an integer                                            */
/* Output parameters description:                               */
/*   - the vector u rotated k coordinate positions              */
/*                                                              */
/* Signature: (<Tup> ~u, <RngIntEly > k)                        */
/*                                                              */
/****************************************************************/
intrinsic Rotate(~u::Tup, k:RngIntElt)
{
Given a vector u belonging to the cartesian product set Z2^alpha x Z4^beta, 
return the vector obtained from u by cyclically shifting its first alpha 
components and last beta components, separetely, to the right by k 
coordinate positions.
}
    require #u eq 2: "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require Type(u[1]) cmpeq ModTupRngElt and Type(u[2]) cmpeq ModTupRngElt: 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (BaseRing(u[2]) cmpeq Z4 and BaseRing(u[1]) cmpeq Z2): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";

    u1 := u[1];
    u2 := u[2];
    Rotate(~u1, k);
    Rotate(~u2, k);

    u[1] := u1;
    u[2] := u2;

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Parent                                        */
/* Parameters: u                                                */
/* Function description: Given a vector u belonging to the      */
/*   cartesian product set Z2^alpha x Z4^beta, return this      */
/*   cartesian product set.                                     */
/* Input parameters description:                                */
/*   - u: A vector in Z2^alpha x Z4^beta                        */
/* Output parameters description:                               */
/*   - The parent of u, that is the cartesian product set       */
/*                                                              */
/* Signature: (<Tup> u)-> ModTupRng                             */
/*                                                              */
/****************************************************************/
intrinsic Parent(u::Tup) -> ModTupRng
{
Given a vector u belonging to the cartesian product set Z2^alpha x Z4^beta, 
return this cartesian product set.
}
    require #u eq 2: "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require Type(u[1]) cmpeq ModTupRngElt and Type(u[2]) cmpeq ModTupRngElt: 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (BaseRing(u[2]) cmpeq Z4 and BaseRing(u[1]) cmpeq Z2): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    
    alpha := Ncols(u[1]);
    beta := Ncols(u[2]);
    n := alpha + beta;
    require n gt 0 : "Argument 1 must be a vector of length greater than 0";
    
    return CartesianProduct(RSpace(Z2, alpha), RSpace(Z4, beta));

end intrinsic;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////                DISTANCE AND WEIGHT                              ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/****************************************************************/
/*                                                              */
/* Function name: LeeWeight                                     */
/* Parameters: u                                                */
/* Function description: The Lee weight of the codeword u, i.e.,*/
/*   the Hamming weight of the binary part (that is, the first  */
/*   alpha coordinates) of u plus the Lee weight of the         */
/*   quaternary part (that is, the rest of the coordinates) of  */
/*   u. Equivalently, it corresponds to the number of non-zero  */
/*   components of Phi(u), where Phi is the Gray map. The       */
/*   codeword u is represented as a tuple in the cartesian      */
/*   product set Z2^alpha x Z4^beta.                            */
/* Input parameters description:                                */
/*   - u: A codeword in Z2^alpha x Z4^beta                      */
/* Output parameters description:                               */
/*   - The Lee weight of u                                      */
/*                                                              */
/* Signature: (<Tup> u) -> RngIntElt                            */
/*                                                              */
/****************************************************************/
intrinsic LeeWeight(u::Tup) -> RngIntElt
{
The Lee weight of the codeword u, i.e., the Hamming weight of the binary part 
(that is, the first alpha coordinates) of u plus the Lee weight of
the quaternary part (that is, the rest of the coordinates) of u.
Equivalently, it corresponds to the number of non-zero components of
Phi(u), where Phi is the Gray map. The codeword u is represented as a tuple in 
the cartesian product set Z2^alpha x Z4^beta.
}
    require #u eq 2: "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require Type(u[1]) cmpeq ModTupRngElt and Type(u[2]) cmpeq ModTupRngElt: 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (BaseRing(u[2]) cmpeq Z4 and BaseRing(u[1]) cmpeq Z2): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    
    alpha := Ncols(u[1]);
    beta := Ncols(u[2]);
    n := alpha + beta;
    require n gt 0 : "Argument 1 must be a vector of length greater than 0";
    u4 := FromVectorZ2Z4toZ4(u);
    
    return Weight(Z2Z4GrayVec(u4, alpha));
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: LeeWeight                                     */
/* Parameters: v, alpha                                         */
/* Function description: The Lee weight of the codeword v, i.e.,*/
/*   the Hamming weight of the binary part (that is, the first  */
/*   alpha coordinates) of v plus the Lee weight of the         */
/*   quaternary part (that is, the rest of the coordinates) of  */
/*   v. Equivalently, it corresponds to the number of non-zero  */
/*   components of Phi(v), where Phi is the Gray map. The code- */
/*   word v is represented as an element in Z4^(alpha+beta) by  */
/*   replacing the ones in the first alpha coordinates by twos. */
/* Input parameters description:                                */
/*   - v: A codeword                                            */
/*   - alpha: The length of the binary part                     */
/* Output parameters description:                               */
/*   - The Lee weight of v                                      */
/*                                                              */
/* Signature: (<ModTupRngElt> v, <RngIntElt> alpha) -> RngIntElt*/
/*                                                              */
/****************************************************************/
intrinsic LeeWeight(v::ModTupRngElt, alpha::RngIntElt) -> RngIntElt
{
The Lee weight of the codeword v, i.e., the Hamming weight of the binary part 
(that is, the first alpha coordinates) of v plus the Lee weight of
the quaternary part (that is, the rest of the coordinates) of v.
Equivalently, it corresponds to the number of non-zero components of
Phi(v), where Phi is the Gray map. The codeword v is represented as an element
in Z4^(alpha+beta) by replacing the ones in the first alpha coordinates by twos.
}
    require (BaseRing(v) eq Z4): "Argument 1 is not a vector over Z4";
    requirerange alpha, 0, Degree(v);
    require (IsZero(2 * ColumnSubmatrix(v, alpha))):
           "First", alpha, "coordinates must be in {0,2}";

    return Weight(Z2Z4GrayVec(v, alpha));
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Weight                                        */
/* Parameters: u                                                */
/* Function description: The Hamming weight of the codeword u,  */
/*    i.e., the number of non-zero components of u. The codeword*/
/*    u is given as a tuple in the cartesian product            */
/*    set Z2^alpha x Z4^beta.                                   */
/* Input parameters description:                                */
/*   - u: A codeword in Z2^alpha x Z4^beta                      */
/* Output parameters description:                               */
/*   - The Hamming weight of u                                  */
/*                                                              */
/* Signature: (<Tup> u) -> RngIntElt                            */
/*                                                              */
/****************************************************************/
intrinsic Weight(u::Tup) -> RngIntElt
{
The Hamming weight of the codeword u, i.e., the number of non-zero
components of u. The codeword u is given as a tuple in the cartesian product
set Z2^alpha x Z4^beta.
}
    require #u eq 2: "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require Type(u[1]) cmpeq ModTupRngElt and Type(u[2]) cmpeq ModTupRngElt: 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (BaseRing(u[2]) cmpeq Z4 and BaseRing(u[1]) cmpeq Z2): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    
    alpha := Ncols(u[1]);
    beta := Ncols(u[2]);
    n := alpha + beta;
    require n gt 0 : "Argument 1 must be a vector of length greater than 0";
    u4 := FromVectorZ2Z4toZ4(u);

    return Weight(u4);
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: LeeDistance                                   */
/* Parameters: u, v, alpha                                      */
/* Function description: The Lee distance between the codewords */
/*   u and v, where u and v belong to the same Z2Z4-additive    */
/*   code C. This is defined to be the Lee weight of u-v. The   */
/*   codewords u and v are represented as elements in Z4^(alpha */
/*   +beta) by replacing the ones in the first alpha coordinates*/
/*   by twos.                                                   */
/* Input parameters description:                                */
/*   - u, v: Two codewords                                      */
/*   - alpha: The length of the binary part                     */
/* Output parameters description:                               */
/*   - The Lee distance between u and v                         */
/*                                                              */
/* Signature: (<ModTupRngElt> u, <ModTupRngElt> v,              */
/*                             <RngIntElt> alpha) -> RngIntElt  */
/*                                                              */
/****************************************************************/
intrinsic LeeDistance(u::ModTupRngElt, v::ModTupRngElt, alpha::RngIntElt) 
                           -> RngIntElt
{
The Lee distance between the codewords u and v, where u and v belong
to the same Z2Z4-additive code C. This is defined to be the Lee weight
of u-v. The codewords u and v are represented as elements in Z4^(alpha+beta) 
by replacing the ones in the first alpha coordinates by twos.
}
    require (BaseRing(u) cmpeq Z4): "Argument 1 is not a vector over Z4";
    require (BaseRing(v) cmpeq Z4): "Argument 2 is not a vector over Z4";
    require (Ncols(v) eq Ncols(u)): 
            "Arguments 1 and 2 must have the same length";
    requirerange alpha, 0, Degree(u);
    require (IsZero(2*ColumnSubmatrix(v, alpha))) and
            (IsZero(2*ColumnSubmatrix(u, alpha))):
            "First", alpha, "coordinates must be in {0,2}";

    return Weight(Z2Z4GrayVec(u-v, alpha));
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: LeeDistance                                   */
/* Parameters: u, v                                             */
/* Function description: The Lee distance between the codewords */
/*   u and v, where u and v belong to the same Z2Z4-additive    */
/*   code C. This is defined to be the Lee weight of u-v. The   */
/*   codewords u and v are given as tuples in the cartesian     */
/*   product set Z2^alpha x Z4^beta.                            */
/* Input parameters description:                                */
/*   - u, v: Two codewords in Z2^alpha x Z4^beta                */
/* Output parameters description:                               */
/*   - The Lee distance between u and v                         */
/*                                                              */
/* Signature: (<Tup> u, <Tup> v) -> RngIntElt                   */
/*                                                              */
/****************************************************************/
intrinsic LeeDistance(u::Tup, v::Tup) -> RngIntElt
{
The Lee distance between the codewords u and v, where u and v belong
to the same Z2Z4-additive code C. This is defined to be the Lee weight
of u-v. The codewords u and v are given as tuples in the cartesian product
set Z2^alpha x Z4^beta.                       
}
    require #u eq 2: "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require Type(u[1]) cmpeq ModTupRngElt and Type(u[2]) cmpeq ModTupRngElt: 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (BaseRing(u[2]) cmpeq Z4 and BaseRing(u[1]) cmpeq Z2): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
			
    require #v eq 2: "Argument 2 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require Type(v[1]) cmpeq ModTupRngElt and Type(v[2]) cmpeq ModTupRngElt: 
            "Argument 2 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (BaseRing(v[2]) cmpeq Z4 and BaseRing(v[1]) cmpeq Z2): 
            "Argument 2 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (Ncols(u[1]) eq Ncols(v[1])) and (Ncols(u[2]) eq Ncols(v[2])):
            "Argument 1 and 2 must have de same length";

    alpha := Ncols(u[1]);
    beta := Ncols(u[2]);
    n := alpha + beta;
    require n gt 0 : "Argument 1 and 2 must be vectors of length greater than 0";
    u4 := FromVectorZ2Z4toZ4(u);
    v4 := FromVectorZ2Z4toZ4(v);

    return Weight(Z2Z4GrayVec(u4-v4, alpha));
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Distance                                      */
/* Parameters: u, v                                             */
/* Function description: The Hamming distance between the       */
/*    codewords u and v, where u and v belong to the same       */
/*    Z2Z4-additive code C. This is defined to be the Hamming   */
/*    weight of u-v. The codewords u and v are given as tuples  */
/*    in the cartesian product set Z2^alpha x Z4^beta.          */
/* Input parameters description:                                */
/*   - u, v: Two codewords in Z2^alpha x Z4^beta                */
/* Output parameters description:                               */
/*   - The Hamming distance between u and v                     */
/*                                                              */
/* Signature: (<Tup> u, <Tup> v) -> RngIntElt                   */
/*                                                              */
/****************************************************************/
intrinsic Distance(u::Tup, v::Tup) -> RngIntElt
{
The Hamming distance between the codewords u and v, where u and
v belong to the same Z2Z4-additive code C. This is defined to be the Hamming weight
of u-v. The codewords u and v are given as tuples in the cartesian product
set Z2^alpha x Z4^beta. 
}
    require #u eq 2: "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require Type(u[1]) cmpeq ModTupRngElt and Type(u[2]) cmpeq ModTupRngElt: 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (BaseRing(u[2]) cmpeq Z4 and BaseRing(u[1]) cmpeq Z2): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
			
    require #v eq 2: "Argument 2 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require Type(v[1]) cmpeq ModTupRngElt and Type(v[2]) cmpeq ModTupRngElt: 
            "Argument 2 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (BaseRing(v[2]) cmpeq Z4 and BaseRing(v[1]) cmpeq Z2): 
            "Argument 2 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (Ncols(u[1]) eq Ncols(v[1])) and (Ncols(u[2]) eq Ncols(v[2])):
            "Argument 1 and 2 must have de same length";

    alpha := Ncols(u[1]);
    beta := Ncols(u[2]);
    n := alpha + beta;
    require n gt 0 : "Argument 1 and 2 must be vectors of length greater than 0";
    u4 := FromVectorZ2Z4toZ4(u);
    v4 := FromVectorZ2Z4toZ4(v);

    return Weight(u4-v4);
end intrinsic;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////                BOOLEAN PREDICATES                               ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/****************************************************************/
/*                                                              */
/* Function name: 'in'                                          */
/* Parameters: u, C                                             */
/* Function description: Return true if and only if the vector  */
/*   u belongs to the Z2Z4-additive code C of type (alpha, beta;*/
/*   gamma, delta; kappa). The vector u is given as a tuple in  */
/*   the cartesian product set Z2^alpha x Z4^beta.              */
/* Input parameters description:                                */
/*   - u: An object of type Tup                                 */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - True if u is in C, and false otherwise                   */
/*                                                              */
/* Signature: (<Tup> u, <Z2Z4Code> C) -> BoolElt                */
/*                                                              */
/****************************************************************/
intrinsic 'in'(u::Tup, C::Z2Z4Code) -> BoolElt
{
Return true if and only if the vector u belongs to the Z2Z4-additive code C of 
type (alpha, beta; gamma, delta; kappa). The vector u is given as a tuple in the 
cartesian product set Z2^alpha x Z4^beta. 
}
    require #u eq 2: "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require Type(u[1]) cmpeq ModTupRngElt and Type(u[2]) cmpeq ModTupRngElt: 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (BaseRing(u[2]) cmpeq Z4 and BaseRing(u[1]) cmpeq Z2): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    
    alpha := Ncols(u[1]);
    beta := Ncols(u[2]);
    n := alpha + beta;
    require n gt 0 : "Argument 1 must be a vector of length greater than 0";
    
    u4 := FromVectorZ2Z4toZ4(u);
    coercible, _ := IsCoercible(C, u4);
    
    return coercible;     
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: 'in'                                          */
/* Parameters: u, C                                             */
/* Function description: Return true if and only if the vector  */
/*   u belongs to the Z2Z4-additive code C of type (alpha, beta;*/
/*   gamma, delta; kappa). The vector u is given as a vector in */ 
/*   Z4^(alpha+beta), where the ones in the first alpha         */
/*   coordinates are represented by twos.                       */
/* Input parameters description:                                */
/*   - u: An object of any type                                 */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - True if u is in C, and false otherwise                   */
/*                                                              */
/* Signature: (<.> u, <Z2Z4Code> C) -> BoolElt                  */
/*                                                              */
/****************************************************************/
intrinsic 'in'(u::., C::Z2Z4Code) -> BoolElt
{
Return true if and only if the vector u belongs to the Z2Z4-additive code C of 
type (alpha, beta; gamma, delta; kappa). The vector u is given as a vector in 
Z4^(alpha+beta), where the ones in the first alpha coordinates are represented by twos.
}
    coercible, _ := IsCoercible(C, u);
    
    return coercible;     
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: 'notin'                                       */
/* Parameters: u, C                                             */
/* Function description: Return true if and only if the vector  */
/*   u does not belong to the Z2Z4-additive code C of type      */
/*   (alpha, beta; gamma, delta; kappa).                        */
/* Input parameters description:                                */
/*   - u: An object of any type                                 */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - True if u is not in C, and false otherwise               */
/*                                                              */
/* Signature: (<.> u, <Z2Z4Code> C) -> BoolElt                  */
/*                                                              */
/****************************************************************/
intrinsic 'notin'(u::., C::Z2Z4Code) -> BoolElt
{
Return true if and only if the vector u does not belong to the Z2Z4-additive 
code C of type (alpha, beta; gamma, delta; kappa). The vector u can be given
either as a vector in Z4^(alpha+beta), where the ones in the first alpha coordinates 
are represented by twos; or as a tuple in the cartesian product set Z2^alpha x Z4^beta. 
}
    return not (u in C);   
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: subset                                        */
/* Parameters: C, D                                             */
/* Function description: Return true if and only if the Z2Z4-   */
/*   additive code C is a subcode of the Z2Z4-additive code D.  */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/*   - D: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - True if C is a subcode of D, and false otherwise         */
/*                                                              */
/* Signature: (<Z2Z4Code> C, <Z2Z4Code> D) -> BoolElt           */
/*                                                              */
/****************************************************************/
intrinsic 'subset'(C::Z2Z4Code, D::Z2Z4Code) -> BoolElt
{
Return true if and only if the Z2Z4-additive code C is a subcode of the
Z2Z4-additive code D.
}
    if (C`Alpha eq D`Alpha) then
        return C`Code subset D`Code;
    else
        return false;
    end if;
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: notsubset                                     */
/* Parameters: C, D                                             */
/* Function description: Return true if and only if the         */
/*   Z2Z4-additive code C is not a subcode of the Z2Z4-additive */
/*   code D.                                                    */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/*   - D: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - True if C is not a subcode of D, and false otherwise     */
/*                                                              */
/* Signature: (<Z2Z4Code> C, <Z2Z4Code> D) -> BoolElt           */
/*                                                              */
/****************************************************************/
intrinsic 'notsubset'(C::Z2Z4Code, D::Z2Z4Code) -> BoolElt
{
Return true if and only if the Z2Z4-additive code C is not a subcode of
the Z2Z4-additive code D.
}
    return not (C subset D);
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: eq                                            */
/* Parameters: C, D                                             */
/* Function description: Return true if and only if the         */
/*   Z2Z4-additive codes C and D are equal.                     */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/*   - D: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - True if C is equal to D, and false otherwise             */
/*                                                              */
/* Signature: (<Z2Z4Code> C, <Z2Z4Code> D) -> BoolElt           */
/*                                                              */
/****************************************************************/
intrinsic 'eq'(C::Z2Z4Code, D::Z2Z4Code) -> BoolElt
{
Return true if and only if the Z2Z4-additive codes C and D are equal.
}
    if (C`Alpha eq D`Alpha) and (C`Length eq D`Length) then
        return C`Code eq D`Code;
    else
        return false;
    end if;
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: ne                                            */
/* Parameters: C, D                                             */
/* Function description: Return true if and only if the         */
/*   Z2Z4-additive codes C and D are not equal.                 */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/*   - D: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - True if C is not equal to D, and false otherwise         */
/*                                                              */
/* Signature: (<Z2Z4Code> C, <Z2Z4Code> D) -> BoolElt           */
/*                                                              */
/****************************************************************/
intrinsic 'ne'(C::Z2Z4Code, D::Z2Z4Code) -> BoolElt
{
Return true if and only if the Z2Z4-additive codes C and D are not equal.
}
    return not(C eq D);
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: IsZero                                        */
/* Parameters: u                                                */
/* Function description: Return true if and only if the         */
/*    codeword u is the zero vector. The codeword u is given as */
/*    a tuple in the cartesian product set Z2^alpha x Z4^beta.  */
/* Input parameters description:                                */
/*   - u: A codeword in Z2^alpha x Z4^beta                      */
/* Output parameters description:                               */
/*   - True if u is the zero vector                             */
/*                                                              */
/* Signature: (<Tup> u) -> BoolElt                              */
/*                                                              */
/****************************************************************/
intrinsic IsZero(u::Tup) -> BoolElt
{
Return true if and only if the codeword u is the zero vector. 
The codeword u i given as a tuple in the cartesian product
set Z2^alpha x Z4^beta.
}
    require #u eq 2: "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require Type(u[1]) cmpeq ModTupRngElt and Type(u[2]) cmpeq ModTupRngElt: 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (BaseRing(u[2]) cmpeq Z4 and BaseRing(u[1]) cmpeq Z2): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    
    alpha := Ncols(u[1]);
    beta := Ncols(u[2]);
    n := alpha + beta;
    require n gt 0 : "Argument 1 must be a vector of length greater than 0";
    u4 := FromVectorZ2Z4toZ4(u);

    return IsZero(u4);

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: IsSelfDual                                    */
/* Parameters: C                                                */
/* Function description: Return true if and only if the         */
/*   Z2Z4-additive code C is additive self-dual, i.e., C equals */
/*   the additive dual code of C.                               */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - True if C is additive self-dual, and false otherwise     */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> BoolElt                         */
/*                                                              */
/****************************************************************/
intrinsic IsSelfDual(C::Z2Z4Code) -> BoolElt
{
Return true if and only if the Z2Z4-additive code C is additive self-dual,
i.e., C equals the additive dual code of C.
}
    return (C eq Dual(C));
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: IsSelfOrthogonal                              */
/* Parameters: C                                                */
/* Function description: Return true if and only if the         */
/*   Z2Z4-additive code C is additive self-orthogonal, that is, */
/*   return whether C is contained in the additive dual code of */
/*   C.                                                         */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - True if C is additive self-orthogonal, and false otherwise*/
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> BoolElt                         */
/*                                                              */
/****************************************************************/
intrinsic IsSelfOrthogonal(C::Z2Z4Code) -> BoolElt
{
Return true if and only if the Z2Z4-additive code C is additive self-
orthogonal, that is, return whether C is contained in the additive dual code of
C.
}
    return (C subset Dual(C));
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: IsZ2Z4AdditiveCode                            */
/* Parameters: C                                                */
/* Function description: Return true if and only if C is a      */
/*   Z2Z4-additive code.                                        */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - True if C is a Z2Z4-additive code, and false otherwise   */
/*                                                              */
/* Signature: (<.> C) -> BoolElt                                */
/*                                                              */
/****************************************************************/
intrinsic IsZ2Z4AdditiveCode(C::.) -> BoolElt
{
Return true if and only if C is a Z2Z4-additive code.
}
    return (Type(C) eq Z2Z4Code);
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: IsSeparable                                   */
/* Parameters: C                                                */
/* Function description: Return true if and only if C is a      */
/*   separable Z2Z4-additive code. A Z2Z4-additive code C of    */
/*   type (alpha, beta; gamma, delta; kappa) is separable if    */
/*   C=CX x CY, where CX and CY are the punctured codes of C by */
/*   deleting the last beta and first alpha coordinates,        */
/*   respectively.                                              */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - True if and only if C is a separable Z2Z4-additive code  */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> BoolElt                         */
/*                                                              */
/****************************************************************/
intrinsic IsSeparable(C::Z2Z4Code) -> BoolElt
{
Return true if and only if C is a separable Z2Z4-additive code.
A Z2Z4-additive code C of type (alpha, beta; gamma, delta; kappa) is separable if
C=CX x CY, where CX and CY are the punctured codes of C by deleting the last 
beta and first alpha coordinates, respectively.
}
    if IsZero(C`Alpha) or IsZero(C`Length-C`Alpha) then
        return true;
    else 
        CX := LinearBinaryCode(C);
        CY := LinearQuaternaryCode(C);
        Cnew := DirectSum(Z2Z4AdditiveCode(CX), Z2Z4AdditiveCode(CY));
        return (C eq Cnew); 
    end if; 
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: IsAntipodal                                   */
/* Parameters: C                                                */
/* Function description: Return true if and only if C is an     */
/*   antipodal Z2Z4-additive code. A Z2Z4-additive code C of    */
/*   type (alpha, beta; gamma, delta; kappa) is antipodal if    */
/*   Cbin = Phi(C) is antipodal, o equivalently if (1..1|2..2)  */
/*   belongs to C.                                              */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - True if and only if C is an antipodal Z2Z4-additive code */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> BoolElt                         */
/*                                                              */
/****************************************************************/
intrinsic IsAntipodal(C::Z2Z4Code) -> BoolElt
{
Return true if and only if C is an antipodal Z2Z4-additive code. A Z2Z4-additive 
code C of type (alpha, beta; gamma, delta; kappa) is antipodal if Cbin = Phi(C) 
is antipodal, o equivalently if (1..1|2..2) belongs to C.
}
    n := Length(C);
    return RSpace(Z4,n)![2^^n] in C;
end intrinsic;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////            CONSTRUCTION NEW CODES FROM OLD                      ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/****************************************************************/
/*                                                              */
/* Function name: Subcode                                       */
/* Parameters: C, t1, t2                                        */
/* Function description: Given a Z2Z4-additive code C of type   */
/*   (alpha, beta; gamma, delta; kappa) and two integers,       */
/*   t1 and t2, such that 1 <= t1 <= gamma and 1 <= t2 <= delta,*/
/*   return a Z2Z4-additive subcode of C of type (alpha, beta;  */
/*   t1, t2; kappa'), where kappa' <= kappa. This Z2Z4-additive */
/*   subcode is generated by the first t1 rows of order two and */
/*   the first t2 rows of order four in the generator matrix    */
/*   given by the function MinRowsGeneratorMatrix(C).           */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/*   - t1: An integer, 1 <= t1 <= gamma                         */
/*   - t2: An integer, 1 <= t2 <= delta                         */
/* Output parameters description:                               */
/*   - The Z2Z4-additive subcode of C                           */
/*                                                              */
/* Signature: (<Z2Z4Code> C, <RngIntElt> t1, <RngIntElt> t2)    */
/*                                                  -> Z2Z4Code */
/*                                                              */
/****************************************************************/
intrinsic Subcode(C::Z2Z4Code, t1::RngIntElt, t2::RngIntElt) -> Z2Z4Code
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma, delta; kappa) and 
two integers, t1 and t2, such that 1 <= t1 <= gamma and 1 <= t2 <= delta, 
return a Z2Z4-additive subcode of C of type (alpha, beta; t1, t2; kappa'), 
where kappa' <= kappa. This Z2Z4-additive subcode is generated by the first t1
rows of order two and the first t2 rows of order four in the generator matrix 
given by the function MinRowsGeneratorMatrix(C).
}
    Order2, Order4 := OrderTwoFourGenerators(C);
    requirerange t1, 0, #Order2;
    requirerange t2, 0, #Order4;
    
    if (t1 eq 0) and (t2 eq 0) then
        return Z2Z4AdditiveZeroCode(C`Alpha, C`Length-C`Alpha);
    end if;
    subC := Z2Z4AdditiveCode(Matrix(Setseq(Order2)[1..t1] cat 
                             Setseq(Order4)[1..t2]), C`Alpha);
    UpdateMinimumLeeWeightLowerBound(~subC, C`MinimumLeeWeightLowerBound);
	
    return subC;
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Subcode                                       */
/* Parameters: C, S1, S2                                        */
/* Function description: Given a Z2Z4-additive code C of type   */
/*   (alpha, beta; gamma, delta; kappa) and two sets of         */
/*   integers, S1 and S2, such that each of their elements      */
/*   lies in the range [1,gamma] and [1,delta], respectively,   */
/*   return a Z2Z4-additive subcode of C of type (alpha, beta;  */
/*   |S1|, |S2|; kappa'), where kappa' <= kappa. This           */
/*   Z2Z4-additive subcode is generated by the rows of order    */
/*   two and four whose positions appear in S1 and S2,          */
/*   respectively, in the generator matrix given by the         */ 
/*   function MinRowsGeneratorMatrix(C).                        */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/*   - S1: A set of integers in the range [1,gamma]             */
/*   - S2: A set of integers in the range [1,delta]             */
/* Output parameters description:                               */
/*   - A Z2Z4-additive subcode of C                             */
/*                                                              */
/* Signature: (<Z2Z4Code> C, <SetEnum> S1, <SetEnum> S2)        */
/*                                                -> Z2Z4Code   */
/*                                                              */
/****************************************************************/
intrinsic Subcode(C::Z2Z4Code, S1::SetEnum, S2::SetEnum) -> Z2Z4Code
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma, delta; kappa) and two
sets of integers, S1 and S2, such that each of their elements lies in the range
[1,gamma] and [1,delta], respectively, return a Z2Z4-additive subcode of C of
type (alpha, beta; |S1|, |S2|; kappa'), where kappa' <= kappa. This 
Z2Z4-additive subcode is generated by the rows of order two and four whose 
positions appear in S1 and S2, respectively, in the generator matrix given by 
the function MinRowsGeneratorMatrix(C).
}
    Order2, Order4 := OrderTwoFourGenerators(C);
    require S1 subset {1..#Order2}: "Argument 2 must be a subset of ",{1..#Order2};
    require S2 subset {1..#Order4}: "Argument 3 must be a subset of ",{1..#Order4};

    if (#S1 eq 0) and (#S2 eq 0) then
        return Z2Z4AdditiveZeroCode(C`Alpha, C`Length-C`Alpha);
    end if;
    subC := Z2Z4AdditiveCode(Matrix( [ Setseq(Order2)[i] : i in S1] cat 
                            [Setseq(Order4)[i] : i in S2]), C`Alpha);
    UpdateMinimumLeeWeightLowerBound(~subC, C`MinimumLeeWeightLowerBound);

    return subC;
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: OrderTwoSubcode                               */
/* Parameters: C                                                */
/* Function description: Given a Z2Z4-additive code C of type   */
/*   (alpha, beta; gamma, delta; kappa) return the Z2Z4-additive*/
/*   subcode Cb which contains all order two codewords of C.    */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - The Z2Z4-additive subcode which contains all order two   */
/*     codewords of C                                           */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> Z2Z4Code                        */
/*                                                              */
/****************************************************************/
intrinsic OrderTwoSubcode(C::Z2Z4Code) -> Z2Z4Code
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma, delta; kappa) return
the Z2Z4-additive subcode Cb which contains all order two codewords of C.
}
    Order2 := OrderTwoSubcodeGenerators(C);
    
    if IsZero(#Order2) then
        return Z2Z4AdditiveZeroCode(C`Alpha, C`Length-C`Alpha);
    end if; 
    subC := Z2Z4AdditiveCode(Setseq(Order2), C`Alpha);

    // update the lower bound from the lower bound of C
    // the new lower bound must be even   
    parityLowerBound := C`MinimumLeeWeightLowerBound mod 2;
    UpdateMinimumLeeWeightLowerBound(~subC, C`MinimumLeeWeightLowerBound + parityLowerBound);
 
    // update the current upper bound to have an even upper bound
    parityUpperBound := subC`MinimumLeeWeightUpperBound mod 2;
    if not IsZero(parityUpperBound) then
        UpdateMinimumLeeWeightUpperBound(~subC, subC`MinimumLeeWeightUpperBound - 1);  
    end if; 

    return subC;
end intrinsic;

/************************************************************/
/*                                                          */
/* Function name: '+'                                       */
/* Parameters: C, D                                         */
/* Function description: The sum of the Z2Z4-additive codes */
/*   C and D, where C and D have the same parameters alpha  */
/*   and beta, hence also the same length.                  */
/* Input parameters description                             */
/*   - C: A Z2Z4-additive code                              */
/*   - D: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - The sum of the Z2Z4-additive codes C and D           */ 
/*                                                          */
/* Signature: (<Z2Z4Code> C, <Z2Z4Code> D) -> Z2Z4Code      */
/*                                                          */
/************************************************************/
intrinsic '+'(C::Z2Z4Code, D::Z2Z4Code) -> Z2Z4Code
{
The sum of the Z2Z4-additive codes C and D, where C
and D have the same parameters alpha and beta, hence also the same length.
}
    require (C`Alpha eq D`Alpha) and ((C`Length-C`Alpha) eq (D`Length-D`Alpha)): 
        "Arguments 1 and 2 must have same parameters alpha and beta";

    newC := Z2Z4AdditiveCode(C`Code + D`Code, C`Alpha);
    newUpperBound := Minimum(C`MinimumLeeWeightUpperBound, D`MinimumLeeWeightUpperBound);
    UpdateMinimumLeeWeightUpperBound(~newC, newUpperBound);

    return newC;
end intrinsic;

/************************************************************/
/*                                                          */
/* Function name: 'meet'                                    */
/* Parameters: C, D                                         */
/* Function description: The intersection of the            */
/*   Z2Z4-additive codes C and D, where C and D have the    */
/*   same parameters alpha and beta, hence also the same    */
/*   length.                                                */
/* Input parameters description                             */
/*   - C: A Z2Z4-additive code                              */
/*   - D: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - The intersection of the Z2Z4-additive codes C and D  */ 
/*                                                          */
/* Signature: (<Z2Z4Code> C, <Z2Z4Code> D) -> Z2Z4Code      */
/*                                                          */
/************************************************************/
intrinsic 'meet'(C::Z2Z4Code, D::Z2Z4Code) -> Z2Z4Code
{
The intersection of the Z2Z4-additive codes C and D, where C and D
have the same parameters alpha and beta, hence also the same length.
}
    require (C`Alpha eq D`Alpha) and ((C`Length-C`Alpha) eq (D`Length-D`Alpha)): 
        "Arguments 1 and 2 must have same parameters alpha and beta";

    newC := Z2Z4AdditiveCode(C`Code meet D`Code, C`Alpha);
    newLowerBound := Maximum(C`MinimumLeeWeightLowerBound, D`MinimumLeeWeightLowerBound);
    UpdateMinimumLeeWeightLowerBound(~newC, newLowerBound);

    return newC; 
end intrinsic;

/************************************************************/
/*                                                          */
/* Function name: DirectSum                                 */
/* Parameters: C, D                                         */
/* Function description: Given Z2Z4-additive codes C and D, */
/*   construct the direct sum of C and D. The direct sum is */
/*   a Z2Z4-additive code that consists of all vectors of   */
/*   the form (u1, v1|u2, v2), where (u1|u2) in C and       */
/*   (v1|v2) in D.                                          */
/* Input parameters description                             */
/*   - C: A Z2Z4-additive code                              */
/*   - D: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - The Z2Z4-additive code direct sum of C and D         */ 
/*                                                          */
/* Signature: (<Z2Z4Code> C, <Z2Z4Code> D) -> Z2Z4Code      */
/*                                                          */
/************************************************************/
intrinsic DirectSum(C::Z2Z4Code, D::Z2Z4Code) -> Z2Z4Code
{
Given Z2Z4-additive codes C and D, construct the direct sum of C and
D. The direct sum is a Z2Z4-additive code that consists of all vectors of
the form (u1, v1|u2, v2), where (u1|u2) in C and (v1|v2) in D.
}
    nc := C`Length;
    nd := D`Length;
    Gc := GeneratorMatrix(C);
    Gd := GeneratorMatrix(D);
    Gc_alpha := ColumnSubmatrix(Gc, 1, C`Alpha);
    Gd_alpha := ColumnSubmatrix(Gd, 1, D`Alpha);
    Gc_beta := ColumnSubmatrix(Gc, C`Alpha+1, nc-C`Alpha);
    Gd_beta := ColumnSubmatrix(Gd, D`Alpha+1, nd-D`Alpha);
    
    M1 := VerticalJoin(Gc_alpha, ZeroMatrix(Z4, Nrows(Gd), C`Alpha));
    M2 := VerticalJoin(ZeroMatrix(Z4, Nrows(Gc), D`Alpha), Gd_alpha);
    M3 := VerticalJoin(Gc_beta, ZeroMatrix(Z4, Nrows(Gd), nc-C`Alpha));
    M4 := VerticalJoin(ZeroMatrix(Z4, Nrows(Gc), nd-D`Alpha), Gd_beta);
    L := HorizontalJoin(<M1, M2, M3, M4>);

    newC := Z2Z4AdditiveCode(L, C`Alpha + D`Alpha);
    if (assigned C`MinimumLeeWeight) and (assigned D`MinimumLeeWeight) then
        newMinWeight := Minimum(C`MinimumLeeWeight, D`MinimumLeeWeight);
        UpdateMinimumLeeWeight(~newC, newMinWeight);
    else 
        newLowerBound := Minimum(C`MinimumLeeWeightLowerBound, D`MinimumLeeWeightLowerBound);
        newUpperBound := Minimum(C`MinimumLeeWeightUpperBound, D`MinimumLeeWeightUpperBound);
        UpdateMinimumLeeWeightLowerBound(~newC, newLowerBound);
        UpdateMinimumLeeWeightUpperBound(~newC, newUpperBound);
    end if;

    return newC;
end intrinsic;

/************************************************************/
/*                                                          */
/* Function name: Concatenation                             */
/* Parameters: C, D                                         */
/* Function description: Given Z2Z4-additive codes C and D, */
/*   return the concatenation of C and D. If                */
/*   Gc = ( A_alpha | A_beta ) and Gd = ( B_alpha | B_beta )*/
/*   are the generator matrices of C and D, respectively,   */
/*   the concatenation of C and D is the Z2Z4-additive code */
/*   with generator matrix whose rows consist of each row   */
/*   (a_alpha | a_beta) of Gc concatenated with each row    */
/*   (b_alpha | b_beta) of Gd in the following way:         */ 
/*   (a_alpha, b_alpha | a_beta, b_beta).                   */
/* Input parameters description                             */
/*   - C: A Z2Z4-additive code                              */
/*   - D: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - The Z2Z4-additive code concatenation of C and D      */ 
/*                                                          */
/* Signature: (<Z2Z4Code> C, <Z2Z4Code> D) -> Z2Z4Code      */
/*                                                          */
/************************************************************/
intrinsic Concatenation(C::Z2Z4Code, D::Z2Z4Code) -> Z2Z4Code
{
Given Z2Z4-additive codes C and D, return the concatenation of C and
D. If Gc = ( A_alpha | A_beta ) and Gd = ( B_alpha | B_beta ) are the generator 
matrices of C and D, respectively, the concatenation of C and D is the 
Z2Z4-additive code with generator matrix whose rows consist of each row 
(a_alpha | a_beta) of Gc concatenated with each row (b_alpha | b_beta) of Gd  
in the following way: (a_alpha, b_alpha | a_beta, b_beta).
}
    Gc := GeneratorMatrix(C);
    Gd := GeneratorMatrix(D);
    mc := Nrows(Gc);
    md := Nrows(Gd);
    nc := Ncols(Gc);
    nd := Ncols(Gd);

    // when C or D is the zero code, the generator matrix has a zero row,
    // and in this case the concatenation coincides with the direct sum.
    if (mc eq 0) then
        Gc := ZeroMatrix(Z4, 1, C`Length);
        mc := 1;
    end if;
    if (md eq 0) then
        Gd := ZeroMatrix(Z4, 1, D`Length);
        md := 1;
    end if;
    
    L := [];
    for z := 1 to mc do
        v := Eltseq(Gc[z]);
        for j := 1 to md do
            w := Eltseq(Gd[j]);
            Ve := v[1..C`Alpha] cat w[1..D`Alpha] cat
                  v[C`Alpha+1..nc] cat w[D`Alpha+1..nd];
            Append(~L, Ve);
        end for;
    end for;

    newC := Z2Z4AdditiveCode(Matrix(Z4,L), C`Alpha+D`Alpha);
    if (assigned C`MinimumLeeWeight) and (assigned D`MinimumLeeWeight) then
        newMinWeight := C`MinimumLeeWeight + D`MinimumLeeWeight;
        UpdateMinimumLeeWeight(~newC, newMinWeight);
    else 
        newLowerBound := C`MinimumLeeWeightLowerBound + D`MinimumLeeWeightLowerBound;
        newUpperBound := C`MinimumLeeWeightUpperBound + D`MinimumLeeWeightUpperBound;
        UpdateMinimumLeeWeightLowerBound(~newC, newLowerBound);
        UpdateMinimumLeeWeightUpperBound(~newC, newUpperBound);
    end if;

    return newC;
end intrinsic;

/************************************************************/
/*                                                          */
/* Function name: PlotkinSum                                */
/* Parameters: C, D                                         */
/* Function description: Given Z2Z4-additive codes C and D  */
/*   both with the same parameters alpha and beta, hence    */
/*   also the same length, construct the Plotkin sum of C   */
/*   and D. The Plotkin sum is a Z2Z4-additive code that    */
/*   consists of all vectors of the form                    */
/*   (u1, u1 + v1 | u2, u2 + v2), where (u1|u2) in C and    */
/*   (v1|v2) in D.                                          */
/* Input parameters description                             */
/*   - C: A Z2Z4-additive code                              */
/*   - D: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - The Z2Z4-additive Plotkin sum code of C and D        */ 
/*                                                          */
/* Signature: (<Z2Z4Code> C, <Z2Z4Code> D) -> Z2Z4Code      */
/*                                                          */
/************************************************************/
intrinsic PlotkinSum(C::Z2Z4Code, D::Z2Z4Code) -> Z2Z4Code
{
Given Z2Z4-additive codes C and D both with the same parameters alpha
and beta, hence also the same length, construct the Plotkin sum of C and
D. The Plotkin sum is a Z2Z4-additive code that consists of all vectors
of the form (u1, u1 + v1 | u2, u2 + v2), where (u1|u2) in C and (v1|v2) in D.
}
    require (C`Alpha eq D`Alpha) and ((C`Length-C`Alpha) eq (D`Length-D`Alpha)):
        "Arguments 1 and 2 must have same parameters alpha and beta";

    Gc := GeneratorMatrix(C);
    Gd := GeneratorMatrix(D);
    Ca := ColumnSubmatrix(Gc, C`Alpha);
    Da := ColumnSubmatrix(Gd, D`Alpha);
    Cb := ColumnSubmatrix(Gc, C`Alpha+1, C`Length-C`Alpha);
    Db := ColumnSubmatrix(Gd, D`Alpha+1, D`Length-D`Alpha);
    M1 := HorizontalJoin(<Ca, Ca, Cb, Cb>);
    Zero1 := ZeroMatrix(Z4, Nrows(Da), Ncols(Ca));
    Zero2 := ZeroMatrix(Z4, Nrows(Db), Ncols(Cb));
    M2 := HorizontalJoin(<Zero1, Da, Zero2, Db>);
    M := VerticalJoin(M1, M2);

    newC := Z2Z4AdditiveCode(M, 2*C`Alpha);
    if (assigned C`MinimumLeeWeight) and (assigned D`MinimumLeeWeight) then
        newMinWeight := Minimum(2*C`MinimumLeeWeight, D`MinimumLeeWeight);
        UpdateMinimumLeeWeight(~newC, newMinWeight);
    else 
        newLowerBound := Minimum(2*C`MinimumLeeWeightLowerBound, D`MinimumLeeWeightLowerBound);
        newUpperBound := Minimum(2*C`MinimumLeeWeightUpperBound, D`MinimumLeeWeightUpperBound);
        UpdateMinimumLeeWeightLowerBound(~newC, newLowerBound);
        UpdateMinimumLeeWeightUpperBound(~newC, newUpperBound);
    end if;

    return newC;
end intrinsic;

/************************************************************/
/*                                                          */
/* Function name: PunctureCode                              */
/* Parameters: C, i                                         */
/* Function description: Given a Z2Z4-additive code C and an*/
/*    integer i, 1 <= i <= alpha + beta, construct a new    */
/*    Z2Z4-additive code C' by deleting the i-th coordinate */
/*    from each codeword of C.                              */
/* Input parameters description                             */
/*   - C: A Z2Z4-additive code                              */
/*   - i: An integer, 1 <= i <= alpha + beta                */
/* Output parameters description:                           */
/*   - The Z2Z4-additive puncture code of C at i position   */ 
/*                                                          */
/* Signature: (<Z2Z4Code> C, <RngIntElt> i) -> Z2Z4Code     */
/*                                                          */
/************************************************************/
intrinsic PunctureCode(C::Z2Z4Code, i::RngIntElt) -> Z2Z4Code
{
Given a Z2Z4-additive code C and an integer i, 1 <= i <= alpha + beta, 
construct a new Z2Z4-additive code C' by deleting the i-th coordinate from each
codeword of C.
}
    requirerange i, 0, C`Length;

    alpha := C`Alpha;
    if (i le C`Alpha) then
        alpha := alpha - 1;
    end if;

    newC := Z2Z4AdditiveCode(PunctureCode(C`Code, i), alpha);
    if zero_Code(newC) then
        return newC;
    end if;
    if assigned C`MinimumLeeWeight then
        d := C`MinimumLeeWeight;
        if d eq 1 then
            oneCodeword := C!0;
            if i le C`Alpha then 
                oneCodeword[i] := 2;
            else
                oneCodeword[i] := 1;
            end if;
            if (not oneCodeword in C) then
                UpdateMinimumLeeWeight(~newC, 1);
            end if; 
        elif d eq 2 then
            UpdateMinimumLeeWeightLowerBound(~newC, d-1);
            UpdateMinimumLeeWeightUpperBound(~newC, d);  
        else 
            UpdateMinimumLeeWeightLowerBound(~newC, d-2);
            UpdateMinimumLeeWeightUpperBound(~newC, d);     
        end if;
    else
        UpdateMinimumLeeWeightLowerBound(~newC, Maximum(1, C`MinimumLeeWeightLowerBound-2));
        UpdateMinimumLeeWeightUpperBound(~newC, C`MinimumLeeWeightUpperBound);     
    end if; 

    return newC;
end intrinsic;

/************************************************************/
/*                                                          */
/* Function name: PunctureCode                              */
/* Parameters: C, S                                         */
/* Function description: Given a Z2Z4-additive code C and a */
/*   set S of distinct integers set(i_1,..., i_r) each of   */
/*   which lies in the range [1, alpha + beta], construct a */
/*   new Z2Z4-additive code C' by deleting the components   */
/*   i_1,..., i_r from each codeword of C.                  */
/* Input parameters description                             */
/*   - C: A Z2Z4-additive code                              */
/*   - S: A set of distinct integers in the range           */
/*        [1, alpha + beta]                                 */
/* Output parameters description:                           */
/*   - The Z2Z4-additive puncture code of C at S positions  */ 
/*                                                          */
/* Signature: (<Z2Z4Code> C, <SetEnum> S) -> Z2Z4Code       */
/*                                                          */
/************************************************************/
intrinsic PunctureCode(C::Z2Z4Code, S::SetEnum) -> Z2Z4Code
{
Given a Z2Z4-additive code C and a set S of distinct integers set(i_1,..., i_r) 
each of which lies in the range [1, alpha + beta], construct a new 
Z2Z4-additive code C' by deleting the components i_1,..., i_r from each 
codeword of C.2 must be a subset of C.
}
    require S subset {1..Length(C`Code)}: 
           "Argument 2 must be a subset of ",{1..Length(C`Code)};
    alpha := C`Alpha;
    for i in S do
        if (i le C`Alpha) then
            alpha := alpha - 1;
        end if;
    end for;

    return Z2Z4AdditiveCode(PunctureCode(C`Code, S), alpha);
end intrinsic;

/************************************************************/
/*                                                          */
/* Function name: ShortenCode                               */
/* Parameters: C, i                                         */
/* Function description: Given a Z2Z4-additive code C and   */
/*   an integer i, 1 <= i <= alpha + beta, construct a new  */
/*   Z2Z4-additive code from C by selecting only those      */
/*   codewords of C having a zero as their i-th component   */
/*   and deleting the i-th component from these codewords.  */
/* Input parameters description                             */
/*   - C: A Z2Z4-additive code                              */
/*   - i: An integer, 1 <= i <= alpha + beta                */
/* Output parameters description:                           */
/*   - The Z2Z4-additive shorten code at position i         */ 
/*                                                          */
/* Signature: (<Z2Z4Code> C, <RngIntElt> i) -> Z2Z4Code     */
/*                                                          */
/************************************************************/
intrinsic ShortenCode(C::Z2Z4Code, i::RngIntElt) -> Z2Z4Code
{
Given a Z2Z4-additive code C and an integer i, 1 <= i <= alpha + beta, 
construct a new Z2Z4-additive code from C by selecting only those codewords of 
C having a zero as their i-th component and deleting the i-th component
from these codewords.
}
    requirerange i, 0, C`Length;

    alpha := C`Alpha;
    if (i le C`Alpha) then
        alpha := alpha - 1;
    end if;

    newC := Z2Z4AdditiveCode(ShortenCode(C`Code, i), alpha);
    UpdateMinimumLeeWeightLowerBound(~newC, C`MinimumLeeWeightLowerBound);

    return newC;
end intrinsic;

/************************************************************/
/*                                                          */
/* Function name: ShortenCode                               */
/* Parameters: C, S                                         */
/* Function description: Given a Z2Z4-additive code C and a */
/*  set S of distinct integers set(i_1,..., i_r) each of    */
/*   which lies in the range [1, alpha + beta], construct a */
/*   new Z2Z4-additive code from C by selecting only those  */
/*   codewords of C having zeros in each of the coordinate  */
/*   positions i_1,..., i_r, and deleting these components. */
/* Input parameters description                             */
/*   - C: A Z2Z4-additive code                              */
/*   - S: A set of distinct integers in the range           */
/*        [1, alpha + beta]                                 */
/* Output parameters description:                           */
/*   - The Z2Z4-additive shorten code at positions S        */  
/*                                                          */
/* Signature: (<Z2Z4Code> C, <SetEnum> S) -> Z2Z4Code       */
/*                                                          */
/************************************************************/
intrinsic ShortenCode(C::Z2Z4Code, S::SetEnum) -> Z2Z4Code
{
Given a Z2Z4-additive code C and a set S of distinct integers set(i_1,..., i_r)
each of which lies in the range [1, alpha + beta], construct a new 
Z2Z4-additive code from C by selecting only those codewords of C having zeros 
in each of the coordinate positions i_1,..., i_r, and deleting these components.
}
    require S subset {1..C`Length}: 
           "Argument 2 must be a subset of ",{1..C`Length};

    alpha := C`Alpha;
    for i in S do
        if (i le C`Alpha) then
            alpha := alpha - 1;
        end if;
    end for;

    newC := Z2Z4AdditiveCode(ShortenCode(C`Code, S), alpha);
    UpdateMinimumLeeWeightLowerBound(~newC, C`MinimumLeeWeightLowerBound);

    return newC;
end intrinsic;

/************************************************************/
/*                                                          */
/* Function name: ExtendCode                                */
/* Parameters: C                                            */
/* Function description: Given a Z2Z4-additive code C form  */
/*   a new Z2Z4-additive code C' from C by adding the       */
/*   appropriate extra binary coordinate to each codeword v */
/*   of such that Phi(v) has even Hamming weight, where Phi */
/*   is the Gray map.                                       */
/* Input parameters description                             */
/*   - C: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - The Z2Z4-additive extend code of C                   */ 
/*                                                          */
/* Signature: (<Z2Z4Code> C) -> Z2Z4Code                    */
/*                                                          */
/************************************************************/
intrinsic ExtendCode(C::Z2Z4Code) -> Z2Z4Code
{
Given a Z2Z4-additive code C form a new Z2Z4-additive code C' from C
by adding the appropriate extra binary coordinate to each codeword v of
C such that Phi(v) has even Hamming weight, where Phi is the Gray map.
}
    beta := C`Length - C`Alpha;
    H := ParityCheckMatrix(C);

    //adds the all-zero column in the binary part of the parity check matrix of C     
    allZeroColumn := Transpose(Matrix(Z4, [[0 : i in [1..Nrows(H)]]] ) );
    newH := HorizontalJoin(allZeroColumn, H);

    //adds the all-twos row (first alpha+1 twos representing binary ones)
    allTwosRow := Matrix(Z4, [[2 : i in [1..C`Alpha+1+beta]]]);
    newH := VerticalJoin(newH, allTwosRow);
  
    newC := Dual(Z2Z4AdditiveCode(newH, C`Alpha+1));
    if assigned C`MinimumLeeWeight then
        UpdateMinimumLeeWeight(~newC, C`MinimumLeeWeight + (C`MinimumLeeWeight mod 2));
    else
        UpdateMinimumLeeWeightLowerBound(~newC, C`MinimumLeeWeightLowerBound + 
                                                C`MinimumLeeWeightLowerBound mod 2);
        UpdateMinimumLeeWeightUpperBound(~newC, C`MinimumLeeWeightUpperBound + 
                                                C`MinimumLeeWeightUpperBound mod 2);
    end if;

    return newC;
end intrinsic;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////            COSET REPRESENTATIVES                                ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/****************************************************************/
/*                                                              */
/* Function name: CosetRepresentatives                          */
/* Parameters: C, S                                             */
/* Function description: Given a Z2Z4-additive code C of type   */
/*   (alpha, beta; gamma, delta; kappa), and a Z2Z4-additive    */
/*   subcode S of C, return a set of coset representatives      */
/*   (not necessarily of minimal weight in their cosets) for S  */
/*   in C as an indexed set of codewords from C. The codewords  */
/*   in C subset Z2^alpha x Z4^beta are represented as elements */
/*   in Z4^(alpha+beta) by replacing the ones in the first alpha*/
/*   coordinates by twos. The set of coset representatives      */
/*   \{c_0, c_1,..., c_t \} satisfies the conditions that       */
/*   c_0 is the zero codeword, and C = U_(i=0)^t (S + c_i).     */
/*   Note that this function is only applicable when S and C    */
/*   are small.                                                 */  
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/*   - S : A subcode of C                                       */
/* Output parameters description:                               */
/*   - Set of quaternary coset representatives for S in C       */
/*   - Set of binary coset representatives for S in C           */
/*                                                              */
/* Remark: This function returns also another parameter not     */
/*   included in the description nor in the manual. This second */
/*   parameter may be used in other funcions.                   */
/*                                                              */
/* Signature: (Z2Z4Code C, Z2Z4Code S) -> SeqEnum, SeqEnum      */
/*                                                              */
/****************************************************************/ 
intrinsic CosetRepresentatives(C::Z2Z4Code, S::Z2Z4Code) -> SeqEnum, SeqEnum
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma, delta; kappa), 
and a Z2Z4-additive subcode S of C, return a set of coset representatives 
(not necessarily of minimal weight in their cosets) for S in C as an indexed set 
of codewords from C. The codewords in C subset Z2^alpha x Z4^beta are represented 
as elements in Z4^(alpha+beta) by replacing the ones in the first alpha 
coordinates by twos. The set of coset representatives \{c_0, c_1,..., c_t \} 
satisfies the conditions that c_0 is the zero codeword, and 
C = U_(i=0)^t (S + c_i). Note that this function is only applicable when S and C 
are small.
}    
    require (S subset C): "Argument 2 must be a subcode of argument 1";

    //short version which seems slower than long version    
    //leadersZ2Z4 := CosetRepresentatives(C`Code, S`Code);
    //V := VectorSpace(GF(2), Z2Z4BinaryLength(C));
    //grayMap := Z2Z4GrayMap(C);
    //leadersZ2 := {@ V!grayMap(v) : v in leadersZ2Z4 @};

    grayMap := GrayMap(C);
    V := VectorSpace(GF(2), BinaryLength(C));
    Q, f := quo<RSpace(C`Code) | RSpace(S`Code)>;
    degreeQ := Degree(Q);
    if degreeQ gt 0 then
        R := RSpace(Z4, degreeQ);
        leadersZ2Z4 := {@ (Q!x)@@f : x in R @};
        leadersZ2 := {@ V!grayMap(v) : v in leadersZ2Z4 @};
    else 
        leadersZ2Z4 := {@ C`Code!0 @};
        leadersZ2 := {@ V!0 @};
    end if;

    return leadersZ2Z4, leadersZ2;
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: CosetRepresentatives                          */
/* Parameters: C                                                */
/* Function description: Given a Z2Z4-additive code C of type   */
/*   (alpha, beta; gamma, delta; kappa), with ambient space     */
/*   V = Z2^alpha x Z4^beta, return a set of coset representa-  */
/*   tives (not necessarily of minimal weight in their cosets)  */
/*   for C in V as an indexed set of vectors from V. The        */
/*   elements in V are represented as elements in               */
/*   Z4^(alpha+beta) by replacing the ones in the first alpha   */
/*   coordinates by twos. The set of coset representatives      */
/*   \{c_0, c_1,..., c_t \} satisfies the conditions that       */
/*   c_0 is the zero codeword, and V = U_(i=0)^t (C + c_i).     */
/*   Note that this function is only applicable when S and C    */
/*   are small.                                                 */  
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/* Output parameters description:                               */
/*   - Set of quaternary coset representatives for C in         */
/*     the ambient space                                        */
/*   - Set of binary coset representatives for C in             */
/*     the ambient space                                        */
/*                                                              */
/* Remark: This function returns also another parameter not     */
/*   included in the description nor in the manual. This second */
/*   parameter may be used in other funcions.                   */
/*                                                              */
/* Signature: (Z2Z4Code C) -> SeqEnum, SeqEnum                  */
/*                                                              */
/****************************************************************/ 
intrinsic CosetRepresentatives(C::Z2Z4Code) -> SeqEnum, SeqEnum
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma, delta; kappa), 
with ambient space V = Z2^alpha x Z4^beta, return a set of coset representatives 
(not necessarily of minimal weight in their cosets) for C in V as an indexed set 
of vectors from V. The elements in V are represented 
as elements in Z4^(alpha+beta) by replacing the ones in the first alpha 
coordinates by twos. The set of coset representatives \{c_0, c_1,..., c_t\}
satisfies the conditions that c_0 is the zero codeword, and 
V = U_(i=0)^t (C + c_i). Note that this function is only applicable when S and C 
are small.
}    
    U := Z2Z4AdditiveUniverseCode(C`Alpha, C`Length-C`Alpha);
    repZ2Z4, repZ2 := CosetRepresentatives(U, C);

    return repZ2Z4, repZ2;
end intrinsic;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////            Z2Z4-ADDITIVE HADAMARD CODES                         ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/****************************************************************/
/* Functions in this section developed by J. Pernas             */
/****************************************************************/

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4HadamardConstr1                           */
/* Parameters: alpha, G                                         */
/* Function description: Given an integer alpha and a generator */
/*   matrix G of a Z2Z4-additive Hadamard code C of binary      */
/*   length 2^m, return an integer alpha' and generator matrix  */
/*   G' for the Z2Z4-additive Hadamard code C' of binary length */
/*   2^(m+1).                                                   */   
/* Input parameters description:                                */
/*   - alpha : integer with the number of binary coordinates    */
/*   - G : generator matrix over Z4 of a Z2Z4-additive          */
/*              Hadamard code                                   */
/* Output parameters description:                               */
/*   - The integer with the number of binary coordinates of C'  */
/*   - A generator matrix over Z4 of C'                         */
/*                                                              */
/****************************************************************/ 
function Z2Z4HadamardConstr1(alpha, G)
    // split generator matrix G in alpha and beta coordinates
    beta := Ncols(G) - alpha; 
    Galpha := ColumnSubmatrix(G, alpha);
    Gbeta := ColumnSubmatrix(G, alpha+1, beta);

    newG := HorizontalJoin(<Galpha, Galpha, Gbeta, Gbeta>); 
    newRow := Vector(Z4, [0^^alpha] cat [2^^alpha] cat [0^^beta] cat [2^^beta]); 
    
    // New row goes first to have all rows of order two in the first rows
    return 2*alpha, Matrix(Rows(newRow) cat Rows(newG));
end function;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4HadamardConstr2                           */
/* Parameters: alpha, G                                         */
/* Function description: Given an integer alpha and a generator */
/*   matrix G of a Z2Z4-additive Hadamard code C of binary      */
/*   length 2^m, return an integer alpha' and generator matrix  */
/*   G' for the Z2Z4-additive Hadamard code C' of binary length */
/*   2^(m+2).                                                   */   
/* Input parameters description:                                */
/*   - alpha : integer with the number of binary coordinates    */
/*   - G : generator matrix over Z4 of a Z2Z4-additive          */
/*              Hadamard code                                   */
/* Output parameters description:                               */
/*   - The integer with the number of binary coordinates of C'  */
/*   - A generator matrix over Z4 of C'                         */
/*                                                              */
/****************************************************************/ 
function Z2Z4HadamardConstr2(alpha, G)
    // split generator matrix in alpha and beta coordinates
    beta := Ncols(G) - alpha;
    Galpha := ColumnSubmatrix(G, alpha);
    Gbeta := ColumnSubmatrix(G, alpha+1, beta);

    newG := HorizontalJoin(<Galpha, Galpha, Galpha, Gbeta, Gbeta, Gbeta, Gbeta>);
    newRow := Vector(Z4, [0^^alpha] cat [2^^alpha] cat [1^^alpha] cat 
                         [0^^beta] cat [1^^beta] cat [2^^beta] cat [3^^beta]);
    
    // New row goes at the end to have all rows of order 4 in the last rows
    return 2*alpha, Matrix(Rows(newG) cat Rows(newRow)); 
end function;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4HadamardMatrix                            */
/* Parameters: delta, m                                         */
/* Function description: Given two integers delta and m, return */
/*   an integer with the number of binary coordinates and a     */
/*   generator matrix over Z4, for a Z2Z4-additive Hadamard     */
/*   code C of type (alpha, beta; gamma, delta; gamma) where    */
/*   alpha = 2^(m-delta), beta = 2^(m-1)-2^(m-delta-1) and      */
/*   gamma = m+1-2*delta.                                       */
/* Input parameters description:                                */
/*   - delta : A positive integer                               */
/*   - m : A positive integer                                   */
/* Output parameters description:                               */
/*   - The integer with the number of binary coordinates of C   */
/*   - A generator matrix over Z4 of C                          */
/*                                                              */
/****************************************************************/ 
function Z2Z4HadamardMatrix(delta, m)    
    case m:
        when 0:
            // base case, type (1,0;1,0)
            G := Matrix(Z4, [[2]]);
            return 1, G;
        else
            if (delta in [0..((m-1) div 2)]) then
                alpha, G := Z2Z4HadamardMatrix(delta, m-1);
                return Z2Z4HadamardConstr1(alpha, G);
            else
                alpha, G := Z2Z4HadamardMatrix(delta-1, m-2);
                return Z2Z4HadamardConstr2(alpha, G);
            end if;
    end case;
end function;

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
intrinsic Z2Z4HadamardCode(delta::RngIntElt, m::RngIntElt : OverZ4 := false) -> Z2Z4Code
{
The function returns a Z2Z4-additive Hadamard code of type (alpha, beta; gamma, 
delta; kappa). The parameter OverZ4 specifies whether the code is over Z4, that is 
alpha=0, or, otherwise, alpha<>0. The default value is false. When OverZ4 is true, 
given an integer m >= 1 and an integer delta such that 1 <= delta <= Floor((m+1)/2), 
return a Z2Z4-additive Hadamard code of type (0, beta; gamma, delta; 0), where 
beta=2^(m-1) and gamma=m+1-2*delta. When OverZ4 is false, given an integer m >= 1 
and an integer delta such that 0 <= delta <= Floor(m/2), return a Z2Z4-additive 
Hadamard code of type (alpha, beta; gamma, delta; gamma), where alpha=2^(m-delta), 
beta=2^(m-1)-2^(m-delta-1) and gamma=m+1-2*delta. Moreover, return a generator 
matrix with gamma+delta rows constructed in a recursive way, where the ones in 
the first alpha coordinates are represented by twos.

A Z2Z4-additive Hadamard code of type (alpha, beta; gamma, delta; kappa) is a 
Z2Z4-additive code such that, after the Gray map, give a binary (not necessarily 
linear) code with the same parameters as the binary Hadamard code of length 2^m.
}
    requirege m, 1;
    if OverZ4 then
        requirerange delta, 1, (m+1) div 2;    
    else
        requirerange delta, 0, m div 2;
    end if;
    
    if OverZ4 then
        CZ4, G := HadamardCodeZ4(delta, m);
        C := Z2Z4AdditiveCode(CZ4);
    else
        alpha, G := Z2Z4HadamardMatrix(delta, m);
        C := Z2Z4AdditiveCode(LinearCode(G), alpha);
    end if;
    UpdateMinimumLeeWeight(~C, 2^(m-1));
    
    return C, G;
end intrinsic;

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
intrinsic Z2Z4ExtendedPerfectCode(delta::RngIntElt, m::RngIntElt : OverZ4 := false) -> Z2Z4Code
{
The function returns a Z2Z4-additive extended perfect code of type alpha, beta; gamma, 
delta; kappa). The parameter OverZ4 specifies whether the code is over Z4, that is 
alpha=0, or, otherwise, alpha<>0. The default value is false. When OverZ4 is true, 
given an integer m >= 1 and an integer delta such that 1 <= delta <= Floor((m+1)/2), 
return a Z2Z4-additive extended perfect code, such that its additive dual code is of 
type of type (0, beta; gamma, delta; 0), where beta=2^(m-1) and gamma=m+1-2*delta. 
When OverZ4 is false, given an integer m >= 1 and an integer delta such that 
0 <= delta <= Floor(m/2), return a Z2Z4-additive extended perfect code, such that 
its additive dual code is of type (alpha, beta; gamma, delta; gamma), where alpha=2^(m-delta), 
beta=2^(m-1)-2^(m-delta-1) and gamma=m+1-2*delta. Moreover, return a parity check 
matrix with gamma+delta rows constructed in a recursive way, where the ones in the 
first alpha coordinates are represented by twos.

A Z2Z4-additive extended perfect code of type (alpha, beta; gamma, delta; kappa) is 
a Z2Z4-additive code such that, after the Gray map, give a binary (not necessarily 
linear) code with the same parameters as the binary extended perfect code of length 2^m.
}
    requirege m, 1;
    if OverZ4 then
        requirerange delta, 1, (m+1) div 2;    
    else
        requirerange delta, 0, m div 2;
    end if;
    
    if OverZ4 then
        CZ4, G := HadamardCodeZ4(delta, m);
        C := Dual(Z2Z4AdditiveCode(CZ4));
        return C, G;
    else
        alpha, G := Z2Z4HadamardMatrix(delta, m);
        C := Dual(Z2Z4AdditiveCode(LinearCode(G), alpha));
        return C, G;
    end if;
    UpdateMinimumLeeWeight(~C, 3);
    
    return C, G;
end intrinsic;

intrinsic KlemmCodeZ4(m::RngIntElt) -> Z2Z4Code
{

}
    requirege m, 1;
    
    n := 4*m;
    firstRow := Matrix([[Z4!1^^n]]);
    lastRows := HorizontalJoin(<ZeroMatrix(Z4, n-2, 1), 
                DiagonalMatrix([Z4!2^^(n-2)]), Matrix(1, [Z4!2^^(n-2)])>);
    G := VerticalJoin(firstRow, lastRows);
    
    return Z2Z4AdditiveCode(G, 0), G;
end intrinsic;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////            Z2Z4-ADDITIVE REED-MULLER CODES                      ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/****************************************************************/
/* Functions in this section first developed by J. Pernas and   */
/*   L. Ronquillo and updated by A. Torres                      */
/****************************************************************/

/****************************************************************/
/*                                                              */
/* Function name: GetQuaternaryCoordinates                      */
/* Parameters: G, alpha                                         */
/* Function description: Given a matrix G over Z4, representing */
/*   a matrix over Z2Z4, return the submatrix of G consisting   */
/*   of all the columns after the first alpha columns. This is, */
/*   the last beta columns, corresponding to quaternary entries.*/
/* Input parameters description:                                */
/*   - G: A matrix over Z2Z4 represented over Z4                */
/*   - alpha: A nonnegative integer.                            */
/* Output parameters description:                               */
/*   - A matrix over Z2Z4 represented over Z4                   */
/*                                                              */
/****************************************************************/
function GetQuaternaryCoordinates(G, alpha)
	return ColumnSubmatrix(G, alpha+1, Ncols(G)-alpha);
end function;

/****************************************************************/
/*                                                              */
/* Function name: GetBinaryCoordinates                          */
/* Parameters: G, alpha                                         */
/* Function description: Given a matrix G over Z4, representing */
/*   a matrix over Z2Z4, return the submatrix of G consisting of*/
/*   the first alpha columns.                                   */
/* Input parameters description:                                */
/*   - G: A matrix over Z2Z4 represented over Z4                */
/*   - alpha: A nonnegative integer.                            */
/* Output parameters description:                               */
/*   - A matrix over Z2Z4 represented over Z4                   */
/*                                                              */
/****************************************************************/
function GetBinaryCoordinates(G,alpha)
	return ColumnSubmatrix(G, 1, alpha);
end function;

/****************************************************************/
/*                                                              */
/* Function name: GetOrder2And4Rows                             */
/* Parameters: G                                                */
/* Function description: Given a matrix G, return two sequences:*/
/*   one with the rows of order two of G and the other with the */
/*   ones of order four.                                        */
/* Input parameters description:                                */
/*   - G: A matrix over Z2Z4 represented over Z4                */
/* Output parameters description:                               */
/*   - A sequence with the rows of order two                    */
/*   - A sequence with the rows of order four                   */
/*                                                              */
/****************************************************************/
function GetOrder2And4Rows(G)
	listOrder2 := [];
	listOrder4 := [];
	n := Ncols(G);
	m := Nrows(G);
	zero := ZeroMatrix(Z4, 1, n);
	
	for i := 1 to m do
		if IsZero(G[i] + G[i]) then 
			Append(~listOrder2, G[i]);
		else
			Append(~listOrder4, G[i]);
		end if;
	end for;
	
	if IsEmpty(listOrder2) then
		matOrder2 := zero;
	else 
		matOrder2 := Matrix(listOrder2);
	end if;
	
	if IsEmpty(listOrder4) then
		matOrder4 := zero;
	else
		matOrder4 := Matrix(listOrder4);	
	end if;
	
	return matOrder2, matOrder4; 
end function;

/****************************************************************/
/*                                                              */
/* Function name: GetBinaryOrder2Coordinates                    */
/* Parameters: G, alpha                                         */
/* Function description: Given a matrix G over Z4, representing */
/*   a matrix over Z2Z4, return the submatrix of G consisting of*/
/*   the first alpha columns of the rows of order two.          */
/* Input parameters description:                                */
/*   - G: A matrix over Z2Z4 represented over Z4                */
/*   - alpha: A nonnegative integer                             */
/* Output parameters description:                               */
/*   - A matrix over Z2Z4 represented over Z4                   */
/*                                                              */
/****************************************************************/
function GetBinaryOrder2Coordinates(G, alpha)
	order2Rows, order4Rows := GetOrder2And4Rows(G);
	return GetBinaryCoordinates(order2Rows, alpha);
end function;

/****************************************************************/
/*                                                              */
/* Function name: GetBinaryOrder4Coordinates                    */
/* Parameters: G, alpha                                         */
/* Function description: Given a matrix G over Z4, representing */
/*   a matrix over Z2Z4, return the submatrix of G consisting of*/
/*   the first alpha columns of the rows of order four.         */
/* Input parameters description:                                */
/*   - G: A matrix over Z2Z4 represented over Z4                */
/*   - alpha: A nonnegative integer                             */
/* Output parameters description:                               */
/*   - A matrix over Z2Z4 represented over Z4                   */
/*                                                              */
/****************************************************************/
function GetBinaryOrder4Coordinates(G, alpha)
	order2Rows, order4Rows := GetOrder2And4Rows(G);
	return GetBinaryCoordinates(order4Rows, alpha);
end function;

/****************************************************************/
/*                                                              */
/* Function name: GetQuaternaryOrder2Coordinates                */
/* Parameters: G, alpha                                         */
/* Function description: Given a matrix G over Z4, representing */
/*   a matrix over Z2Z4, return the rows of G of order two,     */
/*   restricted to the columns from the first alpha columns.    */
/*   This is, the last beta columns, corresponding to quaternary*/
/*   entries.                                                   */
/* Input parameters description:                                */
/*   - G: A matrix over Z2Z4 represented over Z4                */
/*   - alpha: A nonnegative integer                             */
/* Output parameters description:                               */
/*   - A matrix over Z2Z4 represented over Z4                   */
/*                                                              */
/****************************************************************/
function GetQuaternaryOrder2Coordinates(G, alpha)
	order2Rows, order4Rows := GetOrder2And4Rows(G);
	return GetQuaternaryCoordinates(order2Rows, alpha);
end function;

/****************************************************************/
/*                                                              */
/* Function name: GetQuaternaryOrder2Coordinates                */
/* Parameters: G, alpha                                         */
/* Function description: Given a matrix G over Z4, representing */
/*   a matrix over Z2Z4, return the rows of G of order four,    */
/*   restricted to the columns from the first alpha columns.    */
/*   This is, the last beta columns, corresponding to quaternary*/
/*   entries.                                                   */
/* Input parameters description:                                */
/*   - G: A matrix over Z2Z4 represented over Z4                */
/*   - alpha: A nonnegative integer                             */
/* Output parameters description:                               */
/*   - A matrix over Z2Z4 represented over Z4                   */
/*                                                              */
/****************************************************************/
function GetQuaternaryOrder4Coordinates(G, alpha)
	order2Rows, order4Rows := GetOrder2And4Rows(G);
	return GetQuaternaryCoordinates(order4Rows, alpha);
end function;

/****************************************************************/
/*                                                              */
/* Function name: Change2By1                                    */
/* Parameters: G                                                */
/* Function description: Given a matrix G, change all 2s by 1s. */
/* Input parameters description:                                */
/*   - G: A matrix over Z2Z4 represented over Z4                */
/* Output parameters description:                               */
/*   - A matrix over Z4                                         */
/*                                                              */
/****************************************************************/
function Change2By1(G)
	Z2Z4ChangeMatrixZ4toZ2(~G, Ncols(G));
	return G;
end function;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4BAPlotkinConstruction                     */
/* Parameters: A, B, C                                          */
/* Function description: Given three matrices A, C and D over   */
/*   Z2Z4, with the same alpha and beta, construct the          */
/*   BA-Plotkin sum of A, B and C.                              */
/* Input parameters description:                                */
/*   - A: An object of type alphaMatrix                         */
/*   - B: An object of type alphaMatrix                         */
/*   - C: An object of type alphaMatrix                         */
/* Output parameters description:                               */
/*   - An object of type alphaMatrix                            */
/*                                                              */
/****************************************************************/
function Z2Z4BAPlotkinConstruction(A, B, C)
	alpha := A`Alpha;
	beta := Ncols(A`Matrix)-alpha;
	
	Abin := GetBinaryCoordinates(A`Matrix, alpha);
	Aquat := GetQuaternaryCoordinates(A`Matrix, alpha);
	Bbin2 := GetBinaryOrder2Coordinates(B`Matrix, alpha);
	Bbin4 := GetBinaryOrder4Coordinates(B`Matrix, alpha);
	Bquat2 := GetQuaternaryOrder2Coordinates(B`Matrix, alpha);
    Bquat2 := Change2By1(Bquat2);
	Bquat4 := GetQuaternaryOrder4Coordinates(B`Matrix, alpha);
	Cbin := GetBinaryCoordinates(C`Matrix, alpha);
	Cquat := GetQuaternaryCoordinates(C`Matrix, alpha);

	zeroAlpha := ZeroMatrix(Z4, Nrows(B`Matrix), alpha);
	zeroAlpha2 := ZeroMatrix(Z4, Nrows(Bbin2), alpha);
	zeroAlpha4 := ZeroMatrix(Z4, Nrows(Bbin4), alpha);
	zeroAlphac := ZeroMatrix(Z4, Nrows(C`Matrix), alpha);
	zeroBeta2 := ZeroMatrix(Z4, Nrows(Bbin2), beta);
	zeroBeta4 := ZeroMatrix(Z4, Nrows(Bbin4), beta);
	zeroBetac := ZeroMatrix(Z4, Nrows(C`Matrix), beta);
	
	f1 := HorizontalJoin(<Abin, Abin,		Abin, Aquat, Aquat, Aquat, Aquat>);
	f2 := HorizontalJoin(<zeroAlpha2, Bbin2, Change2By1(Bbin2), zeroBeta2, 2*Bquat2, Bquat2, 3*Bquat2>);
	f3 := HorizontalJoin(<zeroAlpha4, Bbin4, Change2By1(Bbin4), zeroBeta4, Bquat4, 2*Bquat4, 3*Bquat4>);
	f4 := HorizontalJoin(<Bbin4, Bbin4, 	zeroAlpha4, zeroBeta4, zeroBeta4, Bquat4, Bquat4>);
	f5 := HorizontalJoin(<zeroAlphac, Cbin, zeroAlphac, zeroBetac, zeroBetac, zeroBetac, Cquat>);
    
	Matrix := VerticalJoin(<f1, f2, f3, f4, f5>);
	
	return rec<alphaMatrix | Alpha := 2*Ncols(Abin), Matrix := Matrix>;
end function;

/****************************************************************/
/*                                                              */
/* Function name: RemoveZeroRows                                */
/* Parameters: G                                                */
/* Function description: Given a matrix G, remove all zero rows.*/
/* Input parameters description:                                */
/*   - G: A matrix over Z2Z4 represented over Z4                */
/* Output parameters description:                               */
/*   - A matrix over Z2Z4 represented over Z4                   */
/*                                                              */
/****************************************************************/
function RemoveZeroRows(G)
	ring := BaseRing(G);
	zero :=	ZeroMatrix(ring, 1, Ncols(G));
	list := [G[i] : i in [1..Nrows(G)] | G[i] ne zero];
	if IsEmpty(list) then
		return Matrix(zero);
	else
		return Matrix(list);
	end if;
end function;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4PlotkinConstruction                       */
/* Parameters: C, D                                             */
/* Function description: Given two matrices C and D over Z2Z4,  */
/*   with the same alpha and beta, construct the Plotkin sum of */
/*   C and D.                                                   */
/* Input parameters description:                                */
/*   - C: An object of type alphaMatrix                         */
/*   - D: An object of type alphaMatrix                         */
/* Output parameters description:                               */
/*   - An object of type alphaMatrix                            */
/*                                                              */
/****************************************************************/
function Z2Z4PlotkinConstruction(C, D)
	G1 := C`Matrix;
	G2 := D`Matrix;
    Ca := ColumnSubmatrix(G1, C`Alpha);
    Da := ColumnSubmatrix(G2, D`Alpha);
    Cb := ColumnSubmatrix(G1, C`Alpha+1, Ncols(G1)-C`Alpha);
    Db := ColumnSubmatrix(G2, D`Alpha+1, Ncols(G2)-D`Alpha);
    M1 := HorizontalJoin(<Ca, Ca, Cb, Cb>);
    Zero1 := ZeroMatrix(Z4, Nrows(Da), Ncols(Ca));
    Zero2 := ZeroMatrix(Z4, Nrows(Db), Ncols(Cb));
    M2 := HorizontalJoin(<Zero1, Da, Zero2, Db>);
    return rec<alphaMatrix | Alpha := (C`Alpha+D`Alpha), 
                                Matrix := RemoveZeroRows(VerticalJoin(M1, M2))>;
end function;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4ReedMullerMatrix                          */
/* Parameters: s, r, m                                          */
/* Function description: Given three integers s, r and m, return*/
/*   a generator matrix for a Z2Z4-additive Reed-Muller code    */
/*   ARM_s(r,m) of type (alpha, beta; gamma, delta), where      */
/*   alpha = 2^(m-s) and beta = 2^(m-1)-2^(m-s-1).              */
/* Input parameters description:                                */
/*   - s : An integer between 0 and Floor(m/2)                  */
/*   - r : An integer between 0 and m                           */
/*   - m : A positive integer                                   */
/* Output parameters description:                               */
/*   - A generator matrix over Z2Z4 of ARM_s(r,m)               */
/*                                                              */
/****************************************************************/
function Z2Z4ReedMullerMatrix(s, r, m)
	case m:
		when 1: // m=1
			case s:
				when 0: // m=1, s=0
					case r:
						when 0:
							//repetition code
							return rec<alphaMatrix | Alpha := 2^m,
								   //Matrix := Matrix(Z4, [[2 : i in [1..2^m]]])>;
                                   Matrix := Matrix(Z4, [[2^^(2^m)]])>;
					    else if r ge 1 then
							return rec<alphaMatrix | Alpha := 2^m,
							       Matrix := Matrix(Z4, [[2,2],[0,2]])>;
						elif r lt 0 then 
						 	return rec<alphaMatrix | Alpha := 2^m,
						           Matrix := ZeroMatrix(Z4, 1, 2^m)>;
                        end if;
					end case; // end case r
			end case; // end case s
		
		when 2: // m=2
			case s: 
				when 1: // m=2, s=1 starting BA family 
					case r:
						when 0: //repetition code
							return rec<alphaMatrix | Alpha := 2^(m-1),
								Matrix := Matrix(Z4, [[2: i in [1..2^(m-1)]] 
								         cat [2 : j in [1..2^(m-1)-2^(m-2)]]])>; 
						when 1: //hadamard
							return rec<alphaMatrix | Alpha := 2^(m-1),
								Matrix := Matrix(Z4, [[2,2,2],[0,2,1]])>;
						else if r ge 2 then
							return rec<alphaMatrix | Alpha := 2^(m-1),
								Matrix := Matrix(Z4, [[2,2,2],[0,2,0],[0,2,1]])>;
					    elif r lt 0 then 
						    return rec<alphaMatrix | Alpha := 2^(m-1),
						     	Matrix := ZeroMatrix(Z4, 1, 2^(m-1)+2^(m-1)-2^(m-2))>;
						end if;
					end case; // end case r
				else 
					return Z2Z4PlotkinConstruction(
						Z2Z4ReedMullerMatrix(s, r, m-1), 
						Z2Z4ReedMullerMatrix(s, r-1, m-1));	
				end case; // end case s
		else // else case m
			case s:
				when m/2:
					A := Z2Z4ReedMullerMatrix(s-1, r, m-2);
					B := Z2Z4ReedMullerMatrix(s-1, r-1, m-2);
					C := Z2Z4ReedMullerMatrix(s-1, r-2, m-2);
					return Z2Z4BAPlotkinConstruction(A, B, C); 
				else
					// Plotkin construction
					return Z2Z4PlotkinConstruction(
                        Z2Z4ReedMullerMatrix(s, r, m-1), 
					    Z2Z4ReedMullerMatrix(s, r-1, m-1));
			end case; //end case s
	end case; // end case m
end function;

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
intrinsic Z2Z4ReedMullerCode(s::RngIntElt, r::RngIntElt, m::RngIntElt : 
                             OverZ4 := false) -> Z2Z4Code, ModMatRngElt
{
The function returns a Z2Z4-additive Reed-Muller code of type 
(alpha, beta; gamma, delta;kappa). The parameter OverZ4 specifies whether the
code is over Z4, that is alpha=0, or, otherwise, alpha<>0. The default value 
is false. When OverZ4 is true, given an integer m>=1, an integer s such that 
0 <= s <= Floor((m-1)/2) and an integer r such that 0 <= r <= m, return an 
r-th order Z2Z4-Additive Reed-Muller Code RM_s(r,m), with alpha = 0 and 
beta = 2^(m-1). In this case, the code coincides with the one obtained by 
using ReedMullerCodeRMZ4(s,r,m). When OverZ4 is false, given m>=1, 
0 <= s <= Floor(m/2) and 0 <= r <= m, return an r-th order Z2Z4-additive 
Reed-MullerCode ARM_s(r,m), with alpha = 2^(m-s) and beta = 2^(m-1)-2^(m-s-1). 
Moreover, in both cases, it returns a generator matrix with gamma+delta rows
constructed in a recursive way, where the ones in the first alpha coordinates 
are represented by twos.

An r-th order Z2Z4-additive Reed-Muller code is a Z2Z4-additive code such that, 
after the Gray map, gives a binary (not necessarily linear) code with the same 
parameters as the binary linear r-th order Reed-Muller code of length 2^m. 
For RM_s(r,m) codes, the inclusion and duality properties are also satisfied, 
that is, RM_s(r-1,m) is a subcode of RM_s(r,m), r>0, and RM_s(r,m) is the Kronecker 
dual code of RM_s(m-r-1,m), r<m. For ARM_s(r,m) codes, only the inclusion property 
is satisfied, that is, ARM_s(r-1,m) is a subcode of ARM_s(r,m), r>0.
}
    requirege m, 1;
    if OverZ4 then
        requirerange s, 0, (m-1) div 2;
        requirerange r, 0, m;    
    else
        requirerange s, 0, m div 2;
        requirerange r, 0, m;
    end if;
    
    if OverZ4 then
        CZ4, G := ReedMullerCodeRMZ4(s, r, m);
        C := Z2Z4AdditiveCode(CZ4);
    else
        alphaG := Z2Z4ReedMullerMatrix(s, r, m);
        G := alphaG`Matrix;
        alpha := alphaG`Alpha;
        C := Z2Z4AdditiveCode(LinearCode(G), alpha);
    end if;
    UpdateMinimumLeeWeight(~C, 2^(m-r));
    
    return C, G;
end intrinsic;

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
intrinsic Z2Z4ReedMullerCodes(s::RngIntElt, m::RngIntElt : 
                              OverZ4 := false) -> Z2Z4Code, ModMatRngElt
{
The function returns a sequence containing a family of Z2Z4-additive 
Reed-Muller codes. The parameter OverZ4 specifies whether the codes are over
Z4, that is alpha=0, or, otherwise, alpha<>0. The default value is false.
When OverZ4 is true, return the codes RM_s(r,m), given any m>=1 and 0<=s<=Floor((m-1)/2).
When OverZ4 is false, return the codes ARM_s(r,m), given any m>=1 and 0<=s<=Floor(m/2).

The binary image of these codes under the Gray map gives a family of binary
(not necessarily linear) codes with the same parameters as the binary linear 
Reed-Muller family of codes of length 2^m. Note that RM_s(0,m) subset RM_s(1,m) 
subset ... subset RM_s(m,m) and ARM_s(0,m) subset ARM_s(1,m) subset ... subset 
ARM_s(m,m).
}
    requirege m, 1;
    if OverZ4 then
        requirerange s, 0, (m-1) div 2;  
    else
        requirerange s, 0, m div 2;
    end if;
    
    return [Z2Z4ReedMullerCode(s, r, m : OverZ4 := OverZ4) : r in [0..m]];
end intrinsic;
