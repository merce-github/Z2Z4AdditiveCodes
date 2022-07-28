///////////////////////////////////////////////////////////////////////////////
/////////       Copyright 2007-2019 Bernat Gaston, Jaume Pujol          ///////
/////////         and Merce Villanueva (with contributions of           ///////
/////////         Lorena Ronquillo and Jaume Pernas)                    ///////
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

 
/*************************************************************/
/*                                                           */
/* Project name: Z2Z4-additive codes in MAGMA                */
/* File name: Z2Z4AdditiveCodes.m                            */
/*                                                           */
/* Comment: Package developed within the CCSG group          */
/*                                                           */
/* Authors: Bernat Gaston, Jaume Pujol and Merce Villanueva  */
/*          (with contributions of Lorena Ronquillo,         */
/*          Jaume Pernas, and Adrian Torres-Martin)          */
/*                                                           */
/* Revision version and last date: v1.7   14-10-2008         */
/*                                 v1.8   09-03-2009         */
/*                                 v1.9   18-03-2010         */    
/*                                 v3.1   14-01-2011         */ 
/*                                 v3.3   15-02-2011         */
/*                                 v3.4   19-06-2011  Split  */
/*                                 v3.5   06-11-2014         */
/*                                 v3.6   13-02-2015         */
/*                                 v3.7   30-05-2016         */
/*                                 v3.8   16-08-2016         */
/*                                 v4.0   10-02-2018         */
/*                                 v4.1   23-02-2018         */
/*                                 v4.2   08-10-2018         */
/*                                 v4.3   26-01-2019         */
/*     User Defined type           v4.4   01-02-2019         */
/*                                 v4.5   14-03-2019  Split  */
/*     Two new attributes added    v4.6   03-06-2022         */
/*                                                           */
/*************************************************************/
//Uncomment freeze when package finished
freeze;

intrinsic Z2Z4AdditiveCodes_version() -> SeqEnum
{Return the current version of this package.}
    version := [4, 6];
    return version;
end intrinsic;

/////        Type declaration 

declare type Z2Z4Code;

declare attributes Z2Z4Code:
    Alpha,                //Length of the binary part
    Length,               //the code length
    Code,                 //the Z4 representation of the code
    MinimumLeeWeight,              //minimum Lee weight of the code
    MinimumLeeWeightLowerBound,    //minimum Lee weight lower bound
    MinimumLeeWeightUpperBound,    //minimum Lee weight upper bound
    MinimumLeeWeightWord,          //codeword of minimum Lee weight
    LeeWeightDistribution,         //Lee weight distribution of the code
    CoveringRadius;                //covering radius  

/////        Global variables 
       
Z4 := Integers(4);
Z2 := Integers(2);
F2 := GF(2);

/////        Constants        

NUMBER_OF_PARAMETERS := 5;    // Parameters which define a Z2Z4-Additive code

/*******************************************************************************
                            Checks
*******************************************************************************/

over_Z4 := func<C | Alphabet(C) cmpeq Z4>;
over_Z2 := func<C | Alphabet(C) cmpeq Z2>;
over_F2 := func<C | Alphabet(C) cmpeq F2>;

zero_Code := func<C | PseudoDimension(C) eq 0>;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////                        CONVERSION FUNCTIONS                     ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/******************************************************************************
                Conversion functions from/to Z4 to/from Z2 and
******************************************************************************/

// M is a matrix over Z4, where the first alpha coordinates are multiplied by 2
procedure Z2Z4ChangeMatrixZ2toZ4(~M, alpha);
    for i in [1..alpha] do
        MultiplyColumn(~M, 2, i);
    end for;
end procedure;

// M is a matrix over Z4, where the first alpha coordinates, which are in {0,2},
// are replaced by {0,1}, respectevely. 
procedure Z2Z4ChangeMatrixZ4toZ2(~M, alpha);
    for i in [1..Nrows(M)] do
        for j in [1..alpha] do
            if M[i][j] eq 2 then 
                M[i][j] := 1;
            end if;
        end for;
    end for;
end procedure;

// M is a matrix over Z4 to check whether the first alpha coordinates are 
// in {0,2} or not
function IsZ2Z4AlphaOverZ4(M, alpha)
    Mt := Transpose(M);
    for i in [1..alpha] do
        if not IsZero(2 * Mt[i]) then
            return false;
        end if;
    end for;
    return true;
end function;

Z2Z4AlphaFromZ4toZ2 := pmap<Z4 -> Z4 | [<0,0>,<2,1>]>;
Z2Z4AlphaFromZ2toZ4 := pmap<Z4 -> Z4 | [<0,0>,<1,2>]>;

function OrderTwoFourGenerators(C)
    Rs, f := StandardForm(C`Code);
    delta, gamma := Z4Type(C`Code);

    return {Rs.i@@f : i in [delta+1..delta+gamma]}, {Rs.i@@f : i in [1..delta]};
end function;

//from Z2Z4 to Z4
//u is a vector in Z2^alpha x Z4^beta
//return a vector in Z4^alpha+beta
function FromVectorZ2Z4toZ4(u)
    u4 := ChangeRing(u[1], Z4);
    U := HorizontalJoin(2*u4, u[2]);

    return U[1];
end function;

//T is a sequence of vectors in Z2^alpha x Z4^beta
//return a matrix in Z4^alpha+beta
function FromSeqZ2Z4toZ4(T)
    alpha := Ncols(T[1][1]);
    Mbin := Matrix([T[i][1] : i in [1..#T]]);
    MbinZ4 := ChangeRing(Mbin, Z4);
    Z2Z4ChangeMatrixZ2toZ4(~MbinZ4, alpha);
    Mqua := Matrix([T[i][2] : i in [1..#T]]);
    MZ4 := HorizontalJoin(MbinZ4, Mqua);

    return MZ4;
end function;

/******************************************************************************
                FromZ4toZ2Z4 function
******************************************************************************/

/****************************************************************/
/*                                                              */
/* Function name: FromZ4toZ2Z4                                  */
/* Parameters: u, alpha                                         */
/* Function description: Given a vector over                    */
/*   Z4^(alpha+beta) and an integer alpha, return               */
/*   the conversion of this vector to a tuple                   */
/*   in the cartesian product set Z2^alpha x Z4^beta, replacing */
/*   the twos over Z4 in the first alpha coordinates            */
/*   to ones over Z2. It is checked whether the elements in the */
/*   first alpha coordinates are in set {0,2}.                  */
/* Input parameters description:                                */
/*   - u: a vector over Z4^(alpha+beta)                         */
/*   - alpha: the length of the binary part of the code         */
/*                                                              */
/* Signature: (<ModTupRngElt> u, <RngIntElt> alpha) -> SeqEnum  */
/*                                                              */
/****************************************************************/
intrinsic FromZ4toZ2Z4(u::ModTupRngElt, alpha::RngIntElt) -> SeqEnum
{
Given a vector over Z4^(alpha+beta) and an integer alpha, return the conversion 
of this vector to a tuple in the cartesian product set Z2^alpha x Z4^beta, 
replacing the twos over Z4 in the first alpha coordinates to ones over Z2. 
It is checked whether the elements in the first alpha coordinates are in set (0,2).
}
    require (Type(BaseRing(u)) cmpeq RngIntRes) and (BaseRing(u) cmpeq Z4): 
            "Argument 1 does not contain a vector over Z4";
    n := Degree(u);
    requirerange alpha, 0, n;
    M := Matrix([u]);
    require IsZ2Z4AlphaOverZ4(M, alpha): 
            "First", alpha, "coordinates must be in {0,2}"; 
            
    Z2Z4ChangeMatrixZ4toZ2(~M, alpha);
    return <Vector(Z2, [M[1][j] : j in [1..alpha]]), 
              Vector(Z4, [M[1][j] : j in [alpha+1..n]])>;

end intrinsic;


/****************************************************************/
/*                                                              */
/* Function name: FromZ4toZ2Z4                                  */
/* Parameters: L, alpha                                         */
/* Function description: Given a sequence L of vectors over     */
/*   Z4^(alpha+beta) and an integer alpha, return               */
/*   the conversion of these vectors to a sequence of tuples    */
/*   in the cartesian product set Z2^alpha x Z4^beta, replacing */
/*   the twos over Z4 in the first alpha coordinates            */
/*   to ones over Z2. It is checked whether the elements in the */
/*   first alpha coordinates are in set {0,2}.                  */
/* Input parameters description:                                */
/*   - L: a sequence of vectors over Z4^(alpha+beta)            */
/*   - alpha: the length of the binary part of the code         */
/*                                                              */
/* Signature: (<SeqEnum> L, <RngIntElt> alpha) -> SeqEnum       */
/*                                                              */
/****************************************************************/
intrinsic FromZ4toZ2Z4(L::SeqEnum, alpha::RngIntElt) -> SeqEnum
{
Given a sequence L of vectors over Z4^(alpha+beta) and an integer alpha, return 
the conversion of these vectors to a sequence of tuples in the cartesian product
set Z2^alpha x Z4^beta, replacing the twos over Z4 in the first alpha coordinates
to ones over Z2. It is checked whether the elements in the first alpha
coordinates are in set (0,2).
}
    require not(IsEmpty(L)): "Argument 1 cannot be empty";
    require ElementType(L) cmpeq ModTupRngElt: "Argument 1 does not contain vectors over Z4";
    require (Type(BaseRing(L[1])) cmpeq RngIntRes) and (BaseRing(L[1]) cmpeq Z4): 
            "Argument 1 does not contain vectors over Z4";
    n := Degree(L[1]);
    requirerange alpha, 0, n;
    M := Matrix(L);
    require IsZ2Z4AlphaOverZ4(M,alpha): 
            "First", alpha, "coordinates must be in {0,2}"; 

    Z2Z4ChangeMatrixZ4toZ2(~M, alpha);
    return [ <Vector(Z2, [M[i][j] : j in [1..alpha]]), 
              Vector(Z4, [M[i][j] : j in [alpha+1..n]])> : i in [1..Nrows(M)] ];

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: FromZ4toZ2Z4                                  */
/* Parameters: S, alpha                                         */
/* Function description: Given a set S of vectors over          */
/*   Z4^(alpha+beta) and an integer alpha, return               */
/*   the conversion of these vectors to a sequence of tuples    */
/*   in the cartesian product set Z2^alpha x Z4^beta, replacing */
/*   the twos over Z4 in the first alpha coordinates            */
/*   to ones over Z2. It is checked whether the elements in the */
/*   first alpha coordinates are in set {0,2}.                  */
/* Input parameters description:                                */
/*   - S: a set of vectors over Z4^(alpha+beta)                 */
/*   - alpha: the length of the binary part of the code         */
/*                                                              */
/* Signature: (<SetEnum> S, <RngIntElt> alpha) -> SetEnum       */
/*                                                              */
/****************************************************************/
intrinsic FromZ4toZ2Z4(S::SetEnum, alpha::RngIntElt) -> SetEnum
{
Given a set S of vectors over Z4^(alpha+beta) and an integer alpha, return 
the conversion of these vectors to a set of tuples in the cartesian product 
set Z2^alpha x Z4^beta, replacing the twos over Z4 in the first alpha coordinates 
to ones over Z2. It is checked whether the elements in the first alpha 
coordinates are in set (0,2).
}
    require not(IsEmpty(S)): "Argument 1 cannot be empty";
    require ElementType(S) cmpeq ModTupRngElt: "Argument 1 does not contain vectors over Z4";
    v := Random(S);
    require (Type(BaseRing(v)) cmpeq RngIntRes) and (BaseRing(v) cmpeq Z4): 
            "Argument 1 does not contain vectors over Z4";
    n := Degree(v);
    requirerange alpha, 0, n;
    M := Matrix(Setseq(S));
    require IsZ2Z4AlphaOverZ4(M,alpha): 
            "First", alpha, "coordinates must be in {0,2}"; 

    Z2Z4ChangeMatrixZ4toZ2(~M, alpha);
    return { <Vector(Z2, [M[i][j] : j in [1..alpha]]), 
              Vector(Z4, [M[i][j] : j in [alpha+1..n]])> : i in [1..Nrows(M)] };

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: FromZ4toZ2Z4                                  */
/* Parameters: M, alpha                                         */
/* Function description: Given matrix M over Z4 with alpha+beta */
/*   columns and an integer alpha, return                       */
/*   the conversion of these vectors to a sequence of tuples    */
/*   in the cartesian product set Z2^alpha x Z4^beta, replacing */
/*   the twos over Z4 in the first alpha coordinates            */
/*   to ones over Z2. It is checked whether the elements in the */
/*   first alpha coordinates are in set {0,2}.                  */
/* Input parameters description:                                */
/*   - M: a matrix M over Z4 with alpha+beta columns            */
/*   - alpha: the length of the binary part of the code         */
/*                                                              */
/* Signature: (<ModMatRngElt> M, <RngIntElt> alpha) -> SeqEnum  */
/*                                                              */
/****************************************************************/
intrinsic FromZ4toZ2Z4(M::ModMatRngElt, alpha::RngIntElt) -> SeqEnum
{
Given a matrix M over Z4 with alpha+beta columns and an integer alpha, return 
the conversion from M to a sequence of tuples in the cartesian product set 
Z2^alpha x Z4^beta, replacing the twos over Z4 in the first alpha coordinates 
to ones over Z2. It is checked whether the elements in the first alpha 
coordinates are in set (0,2).
}
    n := Ncols(M);
    require (n gt 0): "Argument 1 cannot be a matrix with 0 columns";
    requirerange alpha, 0, n;
    require (Type(BaseRing(M)) cmpeq RngIntRes) and (BaseRing(M) cmpeq Z4): 
            "Argument 1 is not a matrix over Z4";
    require IsZ2Z4AlphaOverZ4(M, alpha): 
            "First", alpha, "coordinates must be in {0,2}"; 
        
    Z2Z4ChangeMatrixZ4toZ2(~M, alpha);
    return [<Vector(Z2, [M[i][j] : j in [1..alpha]]), 
             Vector(Z4, [M[i][j]: j in [alpha+1..n]])> : i in [1..Nrows(M)]];

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: FromZ4toZ2Z4                                  */
/* Parameters: M, alpha                                         */
/* Function description: Given a matrix M over Z4 with          */
/*   alpha+beta columns and an integer alpha, return            */
/*   the conversion of these vectors to a sequence of tuples    */
/*   in the cartesian product set Z2^alpha x Z4^beta, replacing */
/*   the twos over Z4 in the first alpha coordinates            */
/*   to ones over Z2. It is checked whether the elements in the */
/*   first alpha coordinates are in set {0,2}.                  */
/* Input parameters description:                                */
/*   - M: a matrix M over Z4 with alpha+beta columns            */
/*   - alpha: the length of the binary part of the code         */
/*                                                              */
/* Signature: (<AlgMatElt> M, <RngIntElt> alpha) -> SeqEnum     */
/*                                                              */
/****************************************************************/
intrinsic FromZ4toZ2Z4(M::AlgMatElt, alpha::RngIntElt) -> SeqEnum
{
Given a matrix M over Z4 with alpha+beta columns and an integer alpha, return 
the conversion from M to a sequence of tuples in the cartesian product set 
Z2^alpha x Z4^beta, replacing the twos over Z4 in the first alpha coordinates 
to ones over Z2. It is checked whether the elements in the first alpha 
coordinates are in set(0,2).
}
    n := Ncols(M);
    require (n gt 0): "Argument 1 cannot be a matrix with 0 columns";
    requirerange alpha, 0, n;
    require (Type(BaseRing(M)) cmpeq RngIntRes) and (BaseRing(M) cmpeq Z4): 
            "Argument 1 is not a matrix over Z4";
    require IsZ2Z4AlphaOverZ4(M,alpha): 
            "First", alpha, "coordinates must be in {0,2}"; 
            
    Z2Z4ChangeMatrixZ4toZ2(~M, alpha);
    return [<Vector(Z2, [M[i][j] : j in [1..alpha]]), 
             Vector(Z4, [M[i][j]: j in [alpha+1..n]])> : i in [1..Nrows(M)]];

end intrinsic;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////            USER DEFINED TYPE: PRINTING                          ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/****************************************************************/
/*                                                              */
/* Procedure name: Print                                        */
/* Parameters: C, L                                             */
/* Procedure description: Given a Z2Z4-additive code C, print   */
/*   the information associated to the Z2Z4-additive code C. The*/
/*   procedure is called automatically by Magma whenever the    */
/*   object C is to be printed.                                 */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/*   - L: "Default", "Minimal", "Maximal", "Magma"              */
/*                                                              */
/* Signature: (<Z2Z4Code> C, <MonStgElt> L)                     */
/*                                                              */
/****************************************************************/
intrinsic Print(C::Z2Z4Code, L::MonStgElt)
{
Print the Z2Z4-additive code C information at level L.
}    
    cardinalCode := #C;
    generatorMatrix := GeneratorMatrix(C`Code);
    typeC := Z2Z4Type(C);
    
    case L:
        when "Default" :
            //N, M, d Z2Z4Code 
            if assigned C`MinimumLeeWeight then
                printf "(%o, %o, %o) Z2Z4-additive code of type (%o, %o; %o, %o; %o)\n", 
                        C`Length, cardinalCode, C`MinimumLeeWeight, 
                        typeC[1], typeC[2], typeC[3], typeC[4], typeC[5];                
            else
                printf "(%o, %o) Z2Z4-additive code of type (%o, %o; %o, %o; %o)\n",
                        C`Length, cardinalCode,
                        typeC[1], typeC[2], typeC[3], typeC[4], typeC[5];
            end if;
            printf "Generator matrix:\n%o", GeneratorMatrix(C`Code);
        when "Minimal" :
            //N, M, d Z2Z4Code 
            if assigned C`MinimumLeeWeight then
                printf "(%o, %o, %o) Z2Z4-additive code of type (%o, %o; %o, %o; %o)", 
                        C`Length, cardinalCode, C`MinimumLeeWeight,
                        typeC[1], typeC[2], typeC[3], typeC[4], typeC[5];              
            else
                printf "(%o, %o) Z2Z4-additive code of type (%o, %o; %o, %o; %o)",
                        C`Length, cardinalCode,
                        typeC[1], typeC[2], typeC[3], typeC[4], typeC[5];
            end if;
        when "Maximal" :
            //N, M, d Z2Z4Code 
            if assigned C`MinimumLeeWeight then
                printf "(%o, %o, %o) Z2Z4-additive code of type (%o, %o; %o, %o; %o)\n",  
                        C`Length, cardinalCode, C`MinimumLeeWeight, 
                        typeC[1], typeC[2], typeC[3], typeC[4], typeC[5];    
                if assigned C`CoveringRadius then
                    printf "having a (%o, %o, %o) binary code as its Gray map image, 
                            with covering radius %o\n",
                        BinaryLength(C), cardinalCode, C`MinimumLeeWeight, C`CoveringRadius; 
                else
                    printf "having a (%o, %o, %o) binary code as its Gray map image \n",
                        BinaryLength(C), cardinalCode, C`MinimumLeeWeight;  
                end if;          
            else
                printf "(%o, %o) Z2Z4-additive code of type (%o, %o; %o, %o; %o)\n",  
                        C`Length, cardinalCode, 
                        typeC[1], typeC[2], typeC[3], typeC[4], typeC[5];
                printf "having a (%o, %o) binary code as its Gray map image \n",
                        BinaryLength(C), cardinalCode;
            end if;
            printf "Generator matrix:\n%o", generatorMatrix;
        when "Magma" : 
            printf "Z2Z4AdditiveCode(%O, %O)", generatorMatrix, "Magma", 
                                          C`Alpha, "Magma";
    end case;

end intrinsic;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////            USER DEFINED TYPE: PARENT                            ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//non needed
//intrinsic Parent(C::Z2Z4Code) -> .
//end intrinsic;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////            USER DEFINED TYPE: COERCION                          ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/****************************************************************/
/*                                                              */
/* Function name: IsMemberOf                                    */
/* Parameters: u, C                                             */
/* Function description: Given a vector u over Z4 such that the */
/*   first alpha coordinates are equal to 1 or 2,               */
/*   return true if u is in the code C and false otherwise.     */
/* Input parameters description:                                */
/*   - u: A vector over Z4                                      */
/*   - C: A q-ary code                                          */
/* Output parameters description:                               */
/*   - true if u is in C, and false otherwise                   */
/*                                                              */
/****************************************************************/
function IsMemberOf(u, C) 
    n := C`Length;    
    alpha := C`Alpha;
    R := RSpace(Z4, n);
    if Parent(u) cmpeq R then
        Mu := Matrix(u);
        if not IsZ2Z4AlphaOverZ4(Mu, alpha) then 
            Z2Z4ChangeMatrixZ2toZ4(~Mu, alpha);
        end if;
        return Mu[1] in C`Code;
    else
        //case u is not coercible to R
        return false;
    end if;
end function;

/****************************************************************/
/*                                                              */
/* Function name: IsCoercible                                   */
/* Parameters: C, c                                             */
/* Function description: Given a Z2Z4-additive code C, and an   */
/*   object c of any type, return whether c is coercible into   */
/*   C and the coerced element if c is coercible to C.          */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/*   - c: An object of any type                                 */
/* Output parameters description:                               */
/*   - true if c is coercible to C, and false otherwise         */
/*   - the object c coerced to C, if c is coercible to C        */
/*                                                              */
/* Signature: (<Z2Z4Code> C, <.> c) -> BoolElt, .               */
/*                                                              */
/****************************************************************/
intrinsic IsCoercible(C::Z2Z4Code, c::.) -> BoolElt, .
{
Return whether c is coercible into C and the result if so.
}
    R4 := RSpace(Z4, C`Length);
    isCoercible, v := IsCoercible(R4, c);
    if isCoercible then
        if IsMemberOf(v, C) then
            return true, v;
        else        
            return false, "Illegal coercion";
        end if;
    else
        return false, "Illegal coercion";
    end if; 
    
end intrinsic;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////        CONSTRUCTION OF GENERAL Z2Z4-ADDITIVE CODES              ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/****************************************************************/
/*                                                              */
/* Procedure name: UpdateMinimumLeeWeight                       */
/* Parameters: ~C, minLeeWeight                                 */
/* Procedure description: Given a Z2Z4-additive code C and an   */
/*   integer with its minimum Lee weight, update the attributes */
/*   of C related to the minimum Lee weight and lower/upper     */
/*   bounds.                                                    */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/*   - minLeeWeight: An integer with the minimum Lee weight of C*/
/*                                                              */
/****************************************************************/
UpdateMinimumLeeWeight := procedure(~C, minLeeWeight)
    C`MinimumLeeWeight := minLeeWeight;
    C`MinimumLeeWeightLowerBound := minLeeWeight;
    C`MinimumLeeWeightUpperBound := minLeeWeight;
end procedure;

/****************************************************************/
/*                                                              */
/* Procedure name: UpdateMinimumLeeWeightLowerBound             */
/* Parameters: ~C, lowerBound                                   */
/* Procedure description: Given a Z2Z4-additive code C and an   */
/*   integer with a lower bound for its minimum Lee weight,     */
/*   update the attribute MinimumLeeWeightLowerBound.           */  
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/*   - lowerBound: An integer with a lower bound for the        */
/*                 minimum Lee weight of C                      */
/*                                                              */
/****************************************************************/
UpdateMinimumLeeWeightLowerBound := procedure(~C, lowerBound)
    if lowerBound gt C`MinimumLeeWeightLowerBound then
        C`MinimumLeeWeightLowerBound := lowerBound;
    end if;
    if C`MinimumLeeWeightLowerBound eq C`MinimumLeeWeightUpperBound then
        C`MinimumLeeWeight := C`MinimumLeeWeightLowerBound;
    end if;
end procedure;

/****************************************************************/
/*                                                              */
/* Procedure name: UpdateMinimumLeeWeightUpperBound             */
/* Parameters: ~C, upperBound                                   */
/* Procedure description: Given a Z2Z4-additive code C and an   */
/*   integer with an upper bound for its minimum Lee weight,    */
/*   update the attribute MinimumLeeWeightUpperBound.           */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/*   - upperBound: An integer with an upper bound for the       */
/*                 minimum Lee weight of C                      */
/*                                                              */
/****************************************************************/
UpdateMinimumLeeWeightUpperBound := procedure(~C, upperBound)
    if upperBound lt C`MinimumLeeWeightUpperBound then
        C`MinimumLeeWeightUpperBound := upperBound;
    end if;
    if C`MinimumLeeWeightLowerBound eq C`MinimumLeeWeightUpperBound then
        C`MinimumLeeWeight := C`MinimumLeeWeightLowerBound;
    end if;
end procedure;

/****************************************************************/
/*                                                              */
/* Procedure name: UpdateMinimumLeeWeightWord                   */
/* Parameters: ~C, word                                         */
/* Procedure description: Given a Z2Z4-additive code C and a    */
/*   codeword with minimum Lee weight, update the attribute     */ 
/*   MinimumLeeWeightWord. At the same time, it also updates    */
/*   the attribute MinimumLeeWeight.                            */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/*   - word: A codeword of C with minimum Lee weight            */
/*                                                              */
/****************************************************************/
UpdateMinimumLeeWeightWord := procedure(~C, word)
    C`MinimumLeeWeightWord := word;
    UpdateMinimumLeeWeight(~C, LeeWeight(word, C`Alpha));
end procedure;

/****************************************************************/
/*                                                              */
/* Procedure name: UpdateLeeWeightDistribution                  */
/* Parameters: ~C, distribution                                 */
/* Procedure description: Given a Z2Z4-additive code C and a    */
/*   sequence of tuples, where the i-th tuple contains the i-th */
/*   weight, wi say, and the number of codewords having weight  */
/*   wi, update the attribute LeeWeightDistribution. At the     */
/*   time, it also updates the attribute LeeWeightDistribution. */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/*   - distribution: Sequence of tuples <Lee weight,            */
/*                                       number of codewords>   */
/*                                                              */
/****************************************************************/
UpdateLeeWeightDistribution := procedure(~C, distribution)
    C`LeeWeightDistribution := distribution;
    if #distribution gt 1 then
        UpdateMinimumLeeWeight(~C, distribution[2][1]);
    end if;
end procedure;

/****************************************************************/
/*                                                              */
/* Function name: NewCodeZ2Z4                                   */
/* Parameters: code, alpha                                      */
/* Function description: Create a new Z2Z4-additive code        */
/*   return the associated Z2Z4-additive code.                  */
/* Input parameters description:                                */
/*   - code: The quaternary linear code equal to the            */
/*           Z2Z4-additive code, where the ones in the first    */
/*           alpha coordinates are represented by twos          */
/*   - alpha: The length of the binary part of the Z2Z4-additive*/
/*            code                                              */
/* Output parameters description:                               */
/*   - C: The new Z2Z4-additive code                            */
/*                                                              */
/****************************************************************/
NewCodeZ2Z4 := function(code, alpha)
    C := New(Z2Z4Code);
    C`Alpha := alpha;
    C`Length := Length(code);
    C`Code := code;
    
    if (#code eq 1) then
        C`MinimumLeeWeightLowerBound := 0;
        C`MinimumLeeWeightUpperBound := 0;
        C`MinimumLeeWeight := 0;
        C`MinimumLeeWeightWord := C`Code!0;
        C`LeeWeightDistribution := [<0,1>];
        return C;
    end if;
    
    // Trivial lower bound
    C`MinimumLeeWeightLowerBound := 1;
    // Singleton upper bound
    singletonBound := 2*C`Length - C`Alpha - Ceiling(Log(2, #C)) +1; 
    C`MinimumLeeWeightUpperBound := singletonBound; 
    if C`MinimumLeeWeightLowerBound eq C`MinimumLeeWeightUpperBound then
        C`MinimumLeeWeight := C`MinimumLeeWeightLowerBound;
    end if;   

    return C;
end function;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4AdditiveCode                              */
/* Parameters: M, alpha                                         */
/* Function description: Create a Z2Z4-additive code C in       */
/*   Z2^alpha x Z4^beta given as a linear code over Z4 after    */
/*   replacing the ones in the first alpha coordinates by twos. */
/*   The corresponding linear code over Z4 of C is obtained by  */
/*   using the function LinearCode() as a subspace of Z4^(alpha */
/*   +beta) generated by a m x n matrix M over the ring Z4.     */
/* Input parameters description:                                */
/*   - M: A matrix M over Z4 with alpha+beta columns            */
/*   - alpha: The length of the binary part of the code         */
/* Output parameters description:                               */
/*   - The Z2Z4-additive code                                   */
/*                                                              */
/* Signature: (<ModMatRngElt> M, <RngIntElt> alpha) -> Z2Z4Code */
/*                                                              */
/****************************************************************/
intrinsic Z2Z4AdditiveCode(M::ModMatRngElt, alpha::RngIntElt) -> Z2Z4Code
{
Create a Z2Z4-additive code C in Z2^alpha x Z4^beta given as a linear code over 
Z4 after replacing the ones in the first alpha coordinates by twos.
The corresponding linear code over Z4 of C is obtained by using the function 
LinearCode() as a subspace of Z4^(alpha+beta) generated by a m x n matrix M over 
the ring Z4.
}
    requirerange alpha, 0, Ncols(M); 
    require (Type(BaseRing(M)) cmpeq RngIntRes) and (BaseRing(M) cmpeq Z4): 
             "Argument 1 is not a matrix over Z4";
     
    if not IsZ2Z4AlphaOverZ4(M, alpha) then
        Z2Z4ChangeMatrixZ2toZ4(~M, alpha);
    end if;

    C := NewCodeZ2Z4(LinearCode(M), alpha);
    return C;

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4AdditiveCode                              */
/* Parameters: M, alpha                                         */
/* Function description: Create a Z2Z4-additive code C in       */
/*   Z2^alpha x Z4^beta given as a linear code over Z4 after    */
/*   replacing the ones in the first alpha coordinates by twos. */
/*   The corresponding linear code over Z4 of C is obtained by  */
/*   using the function LinearCode() as a subspace of Z4^(alpha */
/*   +beta) generated by a m x n matrix M over the ring Z4.     */
/* Input parameters description:                                */
/*   - M: A matrix M over Z4 with alpha+beta columns            */
/*   - alpha: The length of the binary part of the code         */
/* Output parameters description:                               */
/*   - The Z2Z4-additive code                                   */
/*                                                              */
/* Signature: (<AlgMatElt> M, <RngIntElt> alpha) -> Z2Z4Code    */
/*                                                              */
/****************************************************************/
intrinsic Z2Z4AdditiveCode(M::AlgMatElt, alpha::RngIntElt) -> Z2Z4Code
{
Create a Z2Z4-additive code C in Z2^alpha x Z4^beta given as a linear code over 
Z4 after replacing the ones in the first alpha coordinates by twos.
The corresponding linear code over Z4 of C is obtained by using the function 
LinearCode() as a subspace of Z4^(alpha+beta) generated by a m x n matrix M over 
the ring Z4.
}
    requirerange alpha, 0, Ncols(M); 
    require (Type(BaseRing(M)) cmpeq RngIntRes) and (BaseRing(M) cmpeq Z4): 
             "Argument 1 is not a matrix over Z4";
              
    if not IsZ2Z4AlphaOverZ4(M, alpha) then
        Z2Z4ChangeMatrixZ2toZ4(~M, alpha);
    end if;

    C := NewCodeZ2Z4(LinearCode(M), alpha);
    return C;

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4AdditiveCode                              */
/* Parameters: L, alpha                                         */
/* Function description: Create a Z2Z4-additive code C in       */
/*   Z2^alpha x Z4^beta given as a linear code over Z4 after    */
/*   replacing the ones in the first alpha coordinates by twos. */
/*   The corresponding linear code over Z4 of C is obtained by  */
/*   using the function LinearCode() as a subspace of Z4^(alpha */
/*   +beta) generated by the sequence L.                        */
/* Input parameters description:                                */
/*   - L: A sequence of vectors over Z4^(alpha+beta)            */
/*   - alpha: The length of the binary part of the code         */
/* Output parameters description:                               */
/*   - The Z2Z4-additive code                                   */
/*                                                              */
/* Signature: (<SeqEnum> L, <RngIntElt> alpha) -> Z2Z4Code      */
/*                                                              */
/****************************************************************/
intrinsic Z2Z4AdditiveCode(L::SeqEnum, alpha::RngIntElt) -> Z2Z4Code
{
Create a Z2Z4-additive code C in Z2^alpha x Z4^beta given as a linear code over 
Z4 after replacing the ones in the first alpha coordinates by twos.
The corresponding linear code over Z4 of C is obtained by using the function 
LinearCode() as a subspace of Z4^(alpha+beta) generated by the sequence L.
}
    require not(IsEmpty(L)): "Argument 1 cannot be empty";
    require ElementType(L) cmpeq ModTupRngElt: "Argument 1 does not contain vectors over Z4";
    require (Type(BaseRing(L[1])) cmpeq RngIntRes) and (BaseRing(L[1]) cmpeq Z4): 
            "Argument 1 does not contain vectors over Z4";
    requirerange alpha, 0, Degree(L[1]);
    
    M := Matrix(L);
    if not IsZ2Z4AlphaOverZ4(M, alpha) then
        Z2Z4ChangeMatrixZ2toZ4(~M, alpha);
    end if;

    C := NewCodeZ2Z4(LinearCode(M), alpha);
    return C;

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4AdditiveCode                              */
/* Parameters: V, alpha                                         */
/* Function description: Create a Z2Z4-additive code C in       */
/*   Z2^alpha x Z4^beta given as a linear code over Z4 after    */
/*   replacing the ones in the first alpha coordinates by twos. */
/*   The corresponding linear code over Z4 of C is obtained by  */
/*   using the function LinearCode() as a subspace of Z4^(alpha */
/*   +beta) generated by the subspece V of Z4^(alpha+beta).     */
/* Input parameters description:                                */
/*   - V: A subspace V of Z4^n                                  */
/*   - alpha: The length of the binary part of the code         */
/* Output parameters description:                               */
/*   - The Z2Z4-additive code                                   */
/*                                                              */
/* Signature: (<ModTupRng> V, <RngIntElt> alpha) -> Z2Z4Code    */
/*                                                              */
/****************************************************************/
intrinsic Z2Z4AdditiveCode(V::ModTupRng, alpha::RngIntElt) -> Z2Z4Code
{
Create a Z2Z4-additive code C in Z2^alpha x Z4^beta given as a linear code over 
Z4 after replacing the ones in the first alpha coordinates by twos.
The corresponding linear code over Z4 of C is obtained by using the function 
LinearCode() as a subspace of Z4^(alpha+beta) generated by the subspece V of 
Z4^(alpha+beta).
}
    requirerange alpha, 0, Degree(V); 
    require (Type(BaseRing(V)) cmpeq RngIntRes) and (BaseRing(V) cmpeq Z4): 
             "Argument 1 is not a subspace over Z4";

    M := Matrix(Basis(V));
    if not IsZ2Z4AlphaOverZ4(M, alpha) then
        Z2Z4ChangeMatrixZ2toZ4(~M, alpha);
    end if;

    C := NewCodeZ2Z4(LinearCode(M), alpha);
    return C;

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4AdditiveCode                              */
/* Parameters: C                                                */
/* Function description: Given a quaternary linear code C of    */
/*   length beta and type 2^gamma 4^delta, return the Z2Z4-     */
/*   additive code of type (0, beta; gamma, delta; 0)           */
/*   corresponding to C.                                        */
/* Input parameters description:                                */
/*   - C: A quaternary linear code C                            */
/* Output parameters description:                               */
/*   - The Z2Z4-additive code                                   */
/*                                                              */
/* Signature: (<CodeLinRng> C) -> Z2Z4Code                      */
/*                                                              */
/****************************************************************/
intrinsic Z2Z4AdditiveCode(C::CodeLinRng) -> Z2Z4Code
{
Given a quaternary linear code C of length beta and type 2^gamma 4^delta, return 
the Z2Z4-additive code of type (0, beta; gamma, delta; 0) corresponding to C.
}   
    require (over_Z4(C)): "Argument 1 is not a code over Z4"; 

    alpha := 0;
    if zero_Code(C) then
        return Z2Z4AdditiveZeroCode(alpha, Length(C));
    end if;
    
    M := GeneratorMatrix(C);
    if not IsZ2Z4AlphaOverZ4(M, alpha) then
        Z2Z4ChangeMatrixZ2toZ4(~M, alpha);
    end if;

    D := NewCodeZ2Z4(LinearCode(M), alpha);
    
    //assign the minimum Lee weight or update the lower and upper bounds 
    if assigned(C`MinimumLeeWeight) then
        UpdateMinimumLeeWeight(~D, C`MinimumLeeWeight);
    else
        UpdateMinimumLeeWeightLowerBound(~D, C`MinimumWeightLowerBound);
        UpdateMinimumLeeWeightUpperBound(~D, 2*C`MinimumWeightUpperBound); 
    end if;

    return D;

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4AdditiveCode                              */
/* Parameters: C, alpha                                         */
/* Function description: Create a Z2Z4-additive code in Z2^alpha*/
/*   x Z4^beta given as a linear code over Z4 after replacing   */
/*   the ones in the first alpha coordinates by twos. The       */
/*   corresponding linear code over Z4 is obtained directly form*/
/*   the linear code over Z4 C.                                 */
/* Input parameters description:                                */
/*   - C: A quaternary linear code C                            */
/*   - alpha: The length of the binary part of the code         */
/* Output parameters description:                               */
/*   - The Z2Z4-additive code                                   */
/*                                                              */
/* Signature: (<CodeLinRng> C, <RngIntElt> alpha) -> Z2Z4Code   */
/*                                                              */
/****************************************************************/
intrinsic Z2Z4AdditiveCode(C::CodeLinRng, alpha::RngIntElt) -> Z2Z4Code
{
Create a Z2Z4-additive code in Z2^alpha x Z4^beta given as a linear code over 
Z4 after replacing the ones in the first alpha coordinates by twos.
The corresponding linear code over Z4 is obtained directly form the linear 
code over Z4 C.
}
    requirerange alpha, 0, Length(C);
    require (over_Z4(C)): "Argument 1 is not a code over Z4"; 

    if zero_Code(C) then
        return Z2Z4AdditiveZeroCode(alpha, Length(C)-alpha);
    end if;
      
    M := GeneratorMatrix(C);
    if not IsZ2Z4AlphaOverZ4(M, alpha) then
        Z2Z4ChangeMatrixZ2toZ4(~M, alpha);
    end if;

    D := NewCodeZ2Z4(LinearCode(M), alpha);
    
    //assign the minimum Lee weight or update the lower and upper bounds
    if assigned(C`MinimumLeeWeight) then
        UpdateMinimumLeeWeightLowerBound(~D, C`MinimumLeeWeight - alpha);
        UpdateMinimumLeeWeightUpperBound(~D, C`MinimumLeeWeight);
    else
        UpdateMinimumLeeWeightLowerBound(~D, C`MinimumWeightLowerBound);
        UpdateMinimumLeeWeightUpperBound(~D, 2*C`MinimumWeightUpperBound); 
    end if;
    
    return D;

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4AdditiveCode                              */
/* Parameters: C                                                */
/* Function description: Given a binary linear code C of length */
/*   alpha and dimension k, return the Z2Z4-additive code of    */
/*   type (alpha,0;k,0;k) corresponding to C. The code is       */
/*   represented as a linear code over Z4 after replacing the   */
/*   ones by twos.                                              */
/* Input parameters description:                                */
/*   - C: A binary linear code C                                */
/* Output parameters description:                               */
/*   - The Z2Z4-additive code                                   */
/*                                                              */
/* Signature: (<CodeLinFld> C) -> Z2Z4Code                      */
/*                                                              */
/****************************************************************/
intrinsic Z2Z4AdditiveCode(C::CodeLinFld) -> Z2Z4Code
{
Given a binary linear code C of length alpha and dimension k, return the 
Z2Z4-additive code of type (alpha,0;k,0;k) corresponding to C. The code is 
represented as a linear code over Z4 after replacing the ones by twos.
}
    require (over_Z2(C) or over_F2(C)): 
             "Argument 1 is not a code over GF(2) or Z2"; 

    if IsZero(Dimension(C)) then
        return Z2Z4AdditiveZeroCode(Length(C), 0);
    end if;
    
    alpha := Length(C);
    M := 2 * ChangeRing(ChangeRing(GeneratorMatrix(C), Z2), Z4);

    D := NewCodeZ2Z4(LinearCode(M), alpha);
    
    //assign the minimum Lee weight or update the lower and upper bounds
    if assigned(C`MinimumWeight) then
        UpdateMinimumLeeWeight(~D, C`MinimumWeight);
    else  
        UpdateMinimumLeeWeightLowerBound(~D, C`MinimumWeightLowerBound);
        UpdateMinimumLeeWeightUpperBound(~D, C`MinimumWeightUpperBound); 
    end if;
    if assigned(C`MinimumWeightWord) then
        UpdateMinimumLeeWeightWord(~D, D![Z2Z4AlphaFromZ2toZ4(Z2!y) : 
                                          y in Eltseq(C`MinimumWeightWord)]);
    end if;
    if assigned(C`WeightDistribution) then
        UpdateLeeWeightDistribution(~D, C`WeightDistribution);
    end if;

    return D;

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4AdditiveCode                              */
/* Parameters: T                                                */
/* Function description: Given a sequence T of elements of the  */
/*   cartesian product set Z2^alpha x Z4^beta, return the Z2Z4- */
/*   additive code of type (alpha, beta; gamma, delta; kappa)   */
/*   generated by T. The code is represented as a linear code   */
/*   over Z4 after replacing the ones in the first alpha        */
/*   coordinates by twos.                                       */
/* Input parameters description:                                */
/*   - L: A sequence of tuples over Z2^alpha x Z4^beta          */
/* Output parameters description:                               */
/*   - The Z2Z4-additive code                                   */
/*                                                              */
/* Signature: (<SeqEnum> T) -> Z2Z4Code                         */
/*                                                              */
/****************************************************************/
intrinsic Z2Z4AdditiveCode(T::SeqEnum) -> Z2Z4Code
{
Given a sequence T of elements of the cartesian product set Z2^alpha x Z4^beta,
return the Z2Z4-additive code of type (alpha, beta; gamma, delta; kappa) 
generated by T. The code is represented as a linear code over Z4 after replacing 
the ones in the first alpha coordinates by twos.
}
    require not(IsEmpty(T)): "Argument 1 cannot be empty";
    require ElementType(T) cmpeq Tup: "Argument 1 does not contain tuples over Z2 x Z4";  
    require #T[1] eq 2: "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require Type(T[1][1]) cmpeq  ModTupRngElt and Type(T[1][2]) cmpeq ModTupRngElt: 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    require (BaseRing(T[1][2]) cmpeq Z4 and BaseRing(T[1][1]) cmpeq Z2): 
            "Argument 1 is not an element of the cartesian product Z2^alpha x Z4^beta";
    
    alpha := Ncols(T[1][1]);
    M := FromSeqZ2Z4toZ4(T);

    C := NewCodeZ2Z4(LinearCode(M), alpha);
    return C;

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4AdditiveCode                              */
/* Parameters: M, N                                             */
/* Function description: Given a matrix M over Z2 with alpha    */
/*   columns and a matrix N over Z4 with beta columns, return   */
/*   the Z2Z4-additive code of type (alpha, beta; gamma, delta; */
/*   kappa) generated by the rows of the matrix (M | N). When it*/
/*   is necessary, zero rows are added at the end of matrices M */
/*   and N to assure that both have the same number of rows.    */
/*   The code is represented as a linear code over Z4 after     */
/*   replacing the ones in the first alpha coordinates by twos. */
/* Input parameters description:                                */
/*   - M : A matrix over Z2 with alpha columns                  */
/*   - N : A matrix over Z4 with beta columns                   */ 
/* Output parameters description:                               */
/*   - The Z2Z4-additive code                                   */
/*                                                              */
/* Signature: (<ModMatRngElt> M, <ModMatRngElt> N) -> Z2Z4Code  */
/*                                                              */
/****************************************************************/
intrinsic Z2Z4AdditiveCode(M::ModMatRngElt, N::ModMatRngElt) -> Z2Z4Code
{
Given a matrix M over Z2 with alpha columns and a matrix N over Z4 with beta 
columns, return the Z2Z4-additive code of type (alpha, beta; gamma, delta; kappa) 
generated by the rows of the matrix (M | N). When it is necessary, zero rows are 
added at the end of matrices M and N to assure that both have the same number of rows.
The code is represented as a linear code over Z4 after replacing the ones in the 
first alpha coordinates by twos.
} 
    require (Type(BaseRing(M)) cmpeq RngIntRes or Type(BaseRing(M)) cmpeq FldFin) and 
            (BaseRing(M) cmpeq Z2 or BaseRing(M) cmpeq GF(2)): 
             "Argument 1 is not a matrix over Z2 or GF(2)";
    require (Type(BaseRing(N)) cmpeq RngIntRes) and (BaseRing(N) cmpeq Z4): 
             "Argument 2 is not a matrix over Z4";
    
    alpha := Ncols(M);
    beta := Ncols(N);
    difRows := Nrows(N) - Nrows(M);
    if difRows gt 0 then
        zeroM := ZeroMatrix(Z2, difRows, alpha);
        M := ChangeRing(M, Z2);
        M := VerticalJoin(M, zeroM);
    elif difRows lt 0 then
        zeroN := ZeroMatrix(Z4, -difRows, beta);
        N := VerticalJoin(N, zeroN);
    end if;
    MZ4 := ChangeRing(ChangeRing(M, Z2), Z4);
    Z2Z4ChangeMatrixZ2toZ4(~MZ4, alpha);
    MN := HorizontalJoin(MZ4, N);
    
    C := NewCodeZ2Z4(LinearCode(MN), alpha);
    return C;

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4AdditiveCode                              */
/* Parameters: M, N                                             */
/* Function description: Given a matrix M over Z2 with alpha    */
/*   columns and a matrix N over Z4 with beta columns, return   */
/*   the Z2Z4-additive code of type (alpha, beta; gamma, delta; */
/*   kappa) generated by the rows of the matrix (M | N). When it*/
/*   is necessary, zero rows are added at the end of matrices M */
/*   and N to assure that both have the same number of rows.    */
/*   The code is represented as a linear code over Z4 after     */
/*   replacing the ones in the first alpha coordinates by twos. */
/* Input parameters description:                                */
/*   - M : A matrix over Z2 with alpha columns                  */
/*   - N : A matrix over Z4 with beta columns                   */ 
/* Output parameters description:                               */
/*   - The Z2Z4-additive code                                   */
/*                                                              */
/* Signature: (<AlgMatElt> M, <AlgMatElt> N) -> Z2Z4Code        */
/*                                                              */
/****************************************************************/
intrinsic Z2Z4AdditiveCode(M::AlgMatElt, N::AlgMatElt) -> Z2Z4Code
{
Given a matrix M over Z2 with alpha columns and a matrix N over Z4 with beta 
columns, return the Z2Z4-additive code of type (alpha, beta; gamma, delta; kappa) 
generated by the rows of the matrix (M | N). When it is necessary, zero rows are 
added at the end of matrices M and N to assure that both have the same number of rows.
The code is represented as a linear code over Z4 after replacing the ones in the 
first alpha coordinates by twos.
} 
    require (Type(BaseRing(M)) cmpeq RngIntRes or Type(BaseRing(M)) cmpeq FldFin) and 
            (BaseRing(M) cmpeq Z2 or BaseRing(M) cmpeq GF(2)): 
             "Argument 1 is not a matrix over Z2 or GF(2)";
    require (Type(BaseRing(N)) cmpeq RngIntRes) and (BaseRing(N) cmpeq Z4): 
             "Argument 2 is not a matrix over Z4";
    
    alpha := Ncols(M);
    beta := Ncols(N);
    diffRows := Nrows(N) - Nrows(M);
    if diffRows gt 0 then
        zeroM := ZeroMatrix(Z2, diffRows, alpha);
        M := ChangeRing(M, Z2);
        M := VerticalJoin(M, zeroM);
    elif diffRows lt 0 then
        //afegir files a N
        zeroN := ZeroMatrix(Z4, -diffRows, beta);
        N := VerticalJoin(N, zeroN);
    end if;
    MZ4 := ChangeRing(ChangeRing(M, Z2), Z4);
    Z2Z4ChangeMatrixZ2toZ4(~MZ4, alpha);
    MN := HorizontalJoin(MZ4, N);
    
    C := NewCodeZ2Z4(LinearCode(MN), alpha);
    return C;

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4AdditiveCode                              */
/* Parameters: M, N                                             */
/* Function description: Given a matrix M over Z2 with alpha    */
/*   columns and a matrix N over Z4 with beta columns, return   */
/*   the Z2Z4-additive code of type (alpha, beta; gamma, delta; */
/*   kappa) generated by the rows of the matrix (M | N). When it*/
/*   is necessary, zero rows are added at the end of matrices M */
/*   and N to assure that both have the same number of rows.    */
/*   The code is represented as a linear code over Z4 after     */
/*   replacing the ones in the first alpha coordinates by twos. */
/* Input parameters description:                                */
/*   - M : A matrix over Z2 with alpha columns                  */
/*   - N : A matrix over Z4 with beta columns                   */ 
/* Output parameters description:                               */
/*   - The Z2Z4-additive code                                   */
/*                                                              */
/* Signature: (<.> M, <.> N) -> Z2Z4Code                        */
/*                                                              */
/****************************************************************/
intrinsic Z2Z4AdditiveCode(M::., N::.) -> Z2Z4Code
{
Given a matrix M over Z2 with alpha columns and a matrix N over Z4 with beta 
columns, return the Z2Z4-additive code of type (alpha, beta; gamma, delta; kappa) 
generated by the rows of the matrix (M | N). When it is necessary, zero rows are 
added at the end of matrices M and N to assure that both have the same number of rows.
The code is represented as a linear code over Z4 after replacing the ones in the 
first alpha coordinates by twos.
} 
    require (Type(M) cmpeq AlgMatElt) or (Type(M) cmpeq ModMatRngElt) :
             "Argument 1 is not a matrix over Z2 or GF(2)";
    require (Type(N) cmpeq AlgMatElt) or (Type(N) cmpeq ModMatRngElt) :
             "Argument 2 is not a matrix over Z4";
    require (Type(BaseRing(M)) cmpeq RngIntRes or Type(BaseRing(M)) cmpeq FldFin) and 
            (BaseRing(M) cmpeq Z2 or BaseRing(M) cmpeq GF(2)): 
             "Argument 1 is not a matrix over Z2 or GF(2)";
    require (Type(BaseRing(N)) cmpeq RngIntRes) and (BaseRing(N) cmpeq Z4): 
             "Argument 2 is not a matrix over Z4";
    
    alpha := Ncols(M);
    beta := Ncols(N);
    diffRows := Nrows(N) - Nrows(M);
    if diffRows gt 0 then
        zeroM := ZeroMatrix(Z2, diffRows, alpha);
        M := ChangeRing(M, Z2);
        M := VerticalJoin(M, zeroM);
    elif diffRows lt 0 then
        //afegir files a N
        zeroN := ZeroMatrix(Z4, -diffRows, beta);
        N := VerticalJoin(N, zeroN);
    end if;
    MZ4 := ChangeRing(ChangeRing(M, Z2), Z4);
    Z2Z4ChangeMatrixZ2toZ4(~MZ4, alpha);
    MN := HorizontalJoin(MZ4, N);
    
    C := NewCodeZ2Z4(LinearCode(MN), alpha);
    return C;

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4AdditiveCode                              */
/* Parameters: D, E                                             */
/* Function description: Given a binary linear code D of length */
/*   alpha and a quaternary linear code E of length beta, return*/
/*   the Z2Z4-additive code of type (alpha, beta; gamma, delta; */
/*   kappa) generated by the vectors of the form (u|0) and      */
/*   (0|v), where u in D and v in E. The code is represented as */
/*   a linear code over Z4 after replacing the ones in the first*/
/*   alpha coordinates by twos. Note that the Z2Z4-additive code*/
/*   is separable.                                              */
/* Input parameters description:                                */
/*   - D : A linear code over GF(2) of length alpha             */
/*   - E : A linear code over Z4 of length beta                 */ 
/* Output parameters description:                               */
/*   - The Z2Z4-additive code                                   */
/*                                                              */
/* Signature: (<CodeLinFld> D, <CodeLinRng> E) -> Z2Z4Code      */
/*                                                              */
/****************************************************************/
intrinsic Z2Z4AdditiveCode(D::CodeLinFld, E::CodeLinRng) -> Z2Z4Code
{
Given a binary linear code D of length alpha and a quaternary linear code E of 
length beta, return the Z2Z4-additive code of type (alpha, beta; gamma, delta; kappa) 
generated by the vectors of the form (u|0) and (0|v), where u in D and v in E. 
The code is represented as a linear code over Z4 after replacing the ones in the 
first alpha coordinates by twos. Note that the Z2Z4-additive code is separable.
} 
    require (over_Z2(D) or over_F2(D)): "Argument 1 is not a code over GF(2) or Z2"; 
    require (over_Z4(E)): "Argument 2 is not a code over Z4";

    alpha := Length(D);
    GD := GeneratorMatrix(D);
    GE := GeneratorMatrix(E);
    GDZ4 := ChangeRing(ChangeRing(GD, Z2), Z4);
    for i in [1..alpha] do
        MultiplyColumn(~GDZ4, 2, i);
    end for;
    zeroQ := ZeroMatrix(Z4, Nrows(GD), Ncols(GE));
    zeroB := ZeroMatrix(Z4, Nrows(GE), Ncols(GD));
    M := HorizontalJoin(VerticalJoin(GDZ4, zeroB), VerticalJoin(zeroQ, GE));
    
    C := NewCodeZ2Z4(LinearCode(M), alpha);
    if (assigned D`MinimumWeight) and (assigned E`MinimumLeeWeight) then
        newMinWeight := Minimum(D`MinimumWeight, E`MinimumLeeWeight);
        UpdateMinimumLeeWeight(~C, newMinWeight);
    else 
        newLowerBound := Minimum(D`MinimumWeightLowerBound, E`MinimumWeightLowerBound);
        newUpperBound := Minimum(D`MinimumWeightUpperBound, E`MinimumWeightUpperBound);
        UpdateMinimumLeeWeightLowerBound(~C, newLowerBound);
        UpdateMinimumLeeWeightUpperBound(~C, newUpperBound);
    end if;

    return C;

end intrinsic;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////            SOME TRIVIAL Z2Z4-ADDITIVE CODES                     ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4AdditiveUniverseCode                      */
/* Parameters: alpha, beta                                      */
/* Function description: Given two non-negative integers alpha  */
/*   and beta, return the generic Z2Z4-additive code of type    */
/*   (alpha, beta; alpha, beta; alpha) consisting of all        */
/*   possible codewords.                                        */
/* Input parameters description:                                */
/*   - alpha: An integer with the binary length of the code     */
/*   - beta: An integer with the quaternary length of the code  */
/* Output parameters description:                               */
/*   - The Z2Z4-additive universe code                          */
/*                                                              */
/* Signature: (<RngIntElt> alpha, <RngIntElt> beta) -> Z2Z4Code */
/*                                                              */
/****************************************************************/
intrinsic Z2Z4AdditiveUniverseCode(alpha::RngIntElt, beta::RngIntElt) -> Z2Z4Code
{
Given two non-negative integers alpha and beta, return the generic Z2Z4-additive
code of type (alpha, beta; alpha, beta; alpha) consisting of all possible 
codewords.
}
    requirege alpha, 0;
    requirege beta, 0;
    require (alpha + beta) gt 0: 
             "Argument 1 plus argument 2 must be greater than 0";

    C := Z2Z4AdditiveCode(DiagonalMatrix(Z4,([2: i in [1..alpha]] 
                                 cat [1: i in [1..beta]])), alpha);
    //lower bound (1) coincides with singleton upper bound (1), 
    //so the minimum Lee weight is already assigned
    if beta gt 0 then
        UpdateMinimumLeeWeightWord(~C, C`Code!([0^^(alpha+beta-1)] cat [1]));
    else
        UpdateMinimumLeeWeightWord(~C, C`Code!([2] cat [0^^(alpha-1)]));
    end if;                             

    return C;

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4AdditiveZeroCode                          */
/* Parameters: alpha, beta                                      */
/* Function description: Given two non-negative integers alpha  */
/*   and beta, return the generic Z2Z4-additive code of type    */
/*   (alpha, beta; 0, 0; 0) consisting of only the zero codeword*/
/* Input parameters description:                                */
/*   - alpha: An integer with the binary length of the code     */
/*   - beta: An integer with the quaternary length of the code  */
/* Output parameters description:                               */
/*   - The Z2Z4-additive zero code                              */
/*                                                              */
/* Signature: (<RngIntElt> alpha, <RngIntElt> beta) -> Z2Z4Code */
/*                                                              */
/****************************************************************/
intrinsic Z2Z4AdditiveZeroCode(alpha::RngIntElt, beta::RngIntElt) -> Z2Z4Code
{
Given two non-negative integers alpha and beta, return the Z2Z4-additive code of
type (alpha, beta; 0, 0; 0) consisting of only the zero codeword.
}
    requirege alpha, 0;
    requirege beta, 0;
    require (alpha + beta) gt 0: 
             "Argument 1 plus argument 2 must be greater than 0";

    C := Z2Z4AdditiveCode(Matrix(Z4,[[0: i in [1..alpha+beta]]]), alpha);
    //the minimum weight and bounds are assigned in NewCodeZ2Z4 function
    
    return C;   

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4AdditiveRepetitionCode                    */
/* Parameters: alpha, beta                                      */
/* Function description: Given two non-negative integers alpha  */
/*   and beta, return the Z2Z4-additive code of type            */
/*   (alpha, beta; 1, 0; 1) generated by the vector             */
/*   (1,...,1|2,...,2).                                         */
/* Input parameters description:                                */
/*   - alpha: An integer with the binary length of the code     */
/*   - beta: An integer with the quaternary length of the code  */
/* Output parameters description:                               */
/*   - The Z2Z4-additive repetition code                        */
/*                                                              */
/* Signature: (<RngIntElt> alpha, <RngIntElt> beta) -> Z2Z4Code */
/*                                                              */
/****************************************************************/
intrinsic Z2Z4AdditiveRepetitionCode(alpha::RngIntElt, beta::RngIntElt) -> Z2Z4Code
{
Given two non-negative integers alpha and beta, return the Z2Z4-additive code of
type (alpha, beta; 1, 0; 1) generated by the vector (1,...,1|2,...,2).
}
    requirege alpha, 0;
    requirege beta, 0;
    require (alpha + beta) gt 0: 
             "Argument 1 plus argument 2 must be greater than 0";

    C := Z2Z4AdditiveCode(Matrix(Z4,[[2: i in [1..alpha + beta]]]), alpha);
    UpdateMinimumLeeWeight(~C, alpha + 2*beta);
    UpdateMinimumLeeWeightWord(~C, C`Code![2^^(alpha+beta)]);
    UpdateLeeWeightDistribution(~C, [<0,1>, <alpha+2*beta,1>]);

    return C;

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4AdditiveEvenWeightCode                    */
/* Parameters:  alpha, beta                                     */
/* Function description:  Given two non-negative integers alpha */
/*   and beta, return the Z2Z4-additive code of type            */
/*   (alpha, beta; alpha-1,beta; alpha-1) such that all vectors */
/*   have even weight. If beta=0, then alpha must be greater    */
/*    than or equal to 2.                                       */   
/* Input parameters description:                                */
/*   - alpha: An integer with the binary length of the code     */
/*   - beta: An integer with the quaternary length of the code  */
/* Output parameters description:                               */
/*   - A Z2Z4-additive even weight code                         */
/*                                                              */
/* Signature: (<RngIntElt> alpha, <RngIntElt> beta)-> Z2Z4Code  */
/*                                                              */
/****************************************************************/ 
intrinsic Z2Z4AdditiveEvenWeightCode(alpha::RngIntElt, beta::RngIntElt) -> Z2Z4Code
{
Given two non-negative integers alpha and beta return the Z2Z4-additive code of 
type (alpha, beta; alpha-1,beta; alpha-1) such that all vectors have even weight. 
If beta=0, then alpha must be greater than or equal to 2.
}   
    requirege alpha, 0;
    requirege beta, 0;
    require (alpha + beta) gt 0: 
        "Argument 1 plus argument 2 must be greater than 0";

    if IsZero(beta) then 
        requirege alpha, 2;
        C := Z2Z4AdditiveCode(EvenWeightCode(alpha));
        UpdateMinimumLeeWeightWord(~C, C`Code!([2,2] cat [0^^(alpha-2)]));
        return C; 
    else 
        Order2Matrix := HorizontalJoin(<DiagonalMatrix([Z4!2^^alpha]), 
                     ZeroMatrix(Z4, alpha, beta-1), Matrix(1, [Z4!1^^(alpha)])>);
        Order4Matrix := HorizontalJoin(<ZeroMatrix(Z4, beta-1, alpha), 
                     DiagonalMatrix([Z4!1^^(beta-1)]), Matrix(1, [Z4!1^^(beta-1)])>);
        Order2Extra := HorizontalJoin(<ZeroMatrix(Z4, 1, alpha+beta-1), Matrix(Z4,[[2]])>);
        GeneratorMatrixEvenCode := VerticalJoin(<Order2Matrix, Order4Matrix, Order2Extra>);
    
        C := Z2Z4AdditiveCode(GeneratorMatrixEvenCode, alpha);
        UpdateMinimumLeeWeight(~C, 2);
        UpdateMinimumLeeWeightWord(~C, C`Code!([0^^(alpha+beta-1)] cat [2]));
    end if;

    return C;    

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: RandomZ2Z4AdditiveCode                        */
/* Parameters: alpha, beta                                      */
/* Function description: Given two integers alpha and beta,     */
/*   return a random Z2Z4-additive code of type                 */
/*   (alpha, beta; gamma, delta; kappa), where gamma, delta     */
/*   and kappa are computed randomly. Note that there exists    */
/*   a Z2Z4-additive code of type                               */
/*   (alpha, beta; gamma, delta; kappa) if and only if          */
/*   alpha, beta, gamma, delta, kappa >= 0, (alpha+beta)>0,     */
/*   k <= min(alpha,beta) and 0 < delta + gamma <= beta + kappa */
/* Input parameters description:                                */
/*   - alpha: An integer with the binary length of the code     */
/*   - beta: An integer with the quaternary length of the code  */
/* Output parameters description:                               */
/*   - A random Z2Z4-additive code                              */
/*                                                              */
/* Signature: (<RngIntElt> alpha, <RngIntElt> beta)-> Z2Z4Code  */
/*                                                              */
/****************************************************************/
intrinsic RandomZ2Z4AdditiveCode(alpha::RngIntElt, beta::RngIntElt) -> Z2Z4Code
{
Given two integers alpha and beta, return a random Z2Z4-additive code of type 
(alpha, beta; gamma, delta; kappa), where gamma, delta and kappa are computed 
randomly. Note that there exists a Z2Z4-additive code of type (alpha, beta; 
gamma, delta; kappa) if and only if alpha, beta, gamma, delta, kappa >= 0, 
(alpha+beta)>0, k <= min(alpha,beta) and 0 < delta + gamma <= beta + kappa.
}
    requirege alpha,0;
    requirege beta,0;
    require (alpha + beta) gt 0: 
             "Argument 1 plus argument 2 must be greater than 0";

    delta := Random(beta);
    gamma := Random(alpha + beta - delta);
    kappa := Random(Max(0, delta + gamma - beta), Min(alpha, gamma));

    if IsZero(delta) and IsZero(gamma) and IsZero(kappa) then
        return Z2Z4AdditiveZeroCode(alpha, beta);
    end if;
    repeat 

        Mt := ZeroMatrix(Z4, 0, alpha + beta);
        if (kappa gt 0) then
            Ik := ScalarMatrix(Z4, kappa, 2);
            Tb := 2 * ChangeRing(RandomMatrix(Z2, kappa, alpha - kappa), Z4);
            T2_2 := 2 * ChangeRing(RandomMatrix
                    (Z2, kappa, beta - gamma + kappa - delta), Z4);
            Zero1 := ZeroMatrix(Z4, kappa, gamma - kappa + delta);
            Mt := HorizontalJoin(<Ik, Tb, T2_2, Zero1>);
        end if;
        if (gamma - kappa gt 0) then
            Zero1 := ZeroMatrix(Z4, gamma - kappa, alpha);
            Igk_2 := ScalarMatrix(Z4, gamma - kappa, 2);
            T1_2 := 2 * ChangeRing(RandomMatrix
                    (Z2, gamma - kappa, beta - gamma + kappa - delta), Z4);
            Zero2 := ZeroMatrix(Z4, gamma - kappa, delta);
            Mgk := HorizontalJoin(<Zero1, T1_2, Igk_2, Zero2>);
            Mt := VerticalJoin(Mt, Mgk);
        end if;
        if (delta gt 0) then
            Zero1 := ZeroMatrix(Z4, delta, kappa);
            Sb := 2 * ChangeRing(RandomMatrix(Z2, delta, alpha - kappa), Z4);
            Sq := RandomMatrix(Z4, delta, beta - gamma + kappa - delta);
            R := ChangeRing(RandomMatrix(Z2, delta, gamma - kappa), Z4);
            Id := ScalarMatrix(Z4, delta,1);
            Md := HorizontalJoin(<Zero1, Sb, Sq, R, Id>);
            Mt := VerticalJoin(Mt, Md);
        end if;
    until (Mt ne ZeroMatrix(Z4,gamma + delta, alpha + beta));
    
    C := Z2Z4AdditiveCode(Mt, alpha);
    return C; 

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: RandomZ2Z4AdditiveCode                        */
/* Parameters: alpha, beta, gamma                               */
/* Function description: Given three integers alpha, beta and   */
/*   gamma, return a random Z2Z4-additive code of type          */
/*   (alpha, beta; gamma, delta; kappa), where delta and kappa  */
/*   are computed randomly. Note that there exists a            */
/*   Z2Z4-additive code of type (alpha, beta; gamma, delta;     */
/*   kappa) if and only if alpha, beta, gamma, delta,           */
/*   kappa >= 0, (alpha+beta)>0, k <= min(alpha,beta) and       */
/*   0 < delta + gamma <= beta + kappa.                         */
/* Input parameters description:                                */
/*   - alpha: An integer with the binary length of the code     */
/*   - beta: An integer with the quaternary length of the code  */
/*   - gamma: integer with the order two generators of the code */
/* Output parameters description:                               */
/*   - A random Z2Z4-additive code                              */
/*                                                              */
/* Signature: (<RngIntElt> alpha, <RngIntElt> beta,             */
/*                                <RngIntElt> gamma)-> Z2Z4Code */
/*                                                              */
/****************************************************************/
intrinsic RandomZ2Z4AdditiveCode(alpha::RngIntElt, beta::RngIntElt,
                                 gamma::RngIntElt) -> Z2Z4Code
{
Given three integers alpha, beta and gamma, return a random Z2Z4-additive code 
of type (alpha, beta; gamma, delta; kappa), where delta and kappa are computed 
randomly. Note that there exists a Z2Z4-additive code of type (alpha, beta; 
gamma, delta; kappa) if and only if alpha, beta, gamma, delta, kappa >= 0, 
(alpha+beta)>0, k <= min(alpha,beta) and 0 < delta + gamma <= beta + kappa.
}
    requirege alpha, 0;
    requirege beta, 0;
    requirege gamma, 0;
    require (alpha + beta) gt 0: 
        "Argument 1 plus argument must be greater than 0";
    require (alpha + beta) ge gamma: 
        "Argument 1 plus argument 2 must be greater than or equal to argument 3";

    kappa := Random(Max(0, gamma - beta), Min(alpha, gamma));
    delta := Random(beta + kappa - gamma);

    if IsZero(delta) and IsZero(gamma) and IsZero(kappa) then
        return Z2Z4AdditiveZeroCode(alpha, beta);
    end if;
    repeat 
        Mt := ZeroMatrix(Z4, 0, alpha + beta);
        if (kappa gt 0) then
            Ik := ScalarMatrix(Z4, kappa, 2);
            Tb := 2 * ChangeRing(RandomMatrix(Z2, kappa, alpha - kappa), Z4);
            T2_2 := 2 * ChangeRing(RandomMatrix
                   (Z2, kappa, beta - gamma + kappa - delta), Z4);
            Zero1 := ZeroMatrix(Z4, kappa, gamma - kappa + delta);
            Mt := HorizontalJoin(<Ik, Tb, T2_2, Zero1>);
        end if;
        if (gamma - kappa gt 0) then
            Zero1 := ZeroMatrix(Z4, gamma - kappa, alpha);
            Igk_2 := ScalarMatrix(Z4, gamma - kappa, 2);
            T1_2 := 2 * ChangeRing(RandomMatrix
                   (Z2, gamma - kappa, beta - gamma + kappa - delta), Z4);
            Zero2 := ZeroMatrix(Z4, gamma - kappa, delta);
            Mgk := HorizontalJoin(<Zero1, T1_2, Igk_2, Zero2>);
            Mt := VerticalJoin(Mt, Mgk);
        end if;
        if (delta gt 0) then
            Zero1 := ZeroMatrix(Z4, delta, kappa);
            Sb := 2 * ChangeRing(RandomMatrix(Z2, delta, alpha - kappa), Z4);
            Sq := RandomMatrix(Z4, delta, beta - gamma + kappa - delta);
            R := ChangeRing(RandomMatrix(Z2, delta, gamma - kappa), Z4);
            Id := ScalarMatrix(Z4, delta, 1);
            Md := HorizontalJoin(<Zero1, Sb, Sq, R, Id>);
            Mt := VerticalJoin(Mt, Md);
        end if;
    until (Mt ne ZeroMatrix(Z4, gamma + delta, alpha + beta));
    
    C := Z2Z4AdditiveCode(Mt, alpha);
    return C;

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: RandomZ2Z4AdditiveCode                        */
/* Parameters: alpha, beta, gamma, delta                        */
/* Function description: Given four integers alpha, beta, gamma */
/*   and delta, return a random Z2Z4-additive code of type      */
/*   (alpha, beta; gamma, delta; kappa), where kappa is         */
/*   computed randomly. Note that there exists a Z2Z4-additive  */
/*   code of type (alpha, beta; gamma, delta; kappa) if and     */
/*   only if alpha, beta, gamma, delta, kappa >= 0,             */
/*   (alpha+beta)>0, k <= min(alpha,beta) and                   */
/*   0 < delta + gamma <= beta + kappa.                         */
/* Input parameters description:                                */
/*   - alpha: integer with the binary length of the code        */
/*   - beta: integer with the quaternary length of the code     */
/*   - gamma: integer with the order two generators of the code */
/*   - delta: integer with the order four generators of the code*/
/* Output parameters description:                               */
/*   - A random Z2Z4-additive code                              */
/*                                                              */
/* Signature: (<RngIntElt> alpha, <RngIntElt> beta,             */
/*             <RngIntElt> gamma, <RngIntElt> delta)-> Z2Z4Code */
/*                                                              */
/****************************************************************/
intrinsic RandomZ2Z4AdditiveCode(alpha::RngIntElt, beta::RngIntElt,
                                 gamma::RngIntElt, delta::RngIntElt) -> Z2Z4Code
{
Given four integers alpha, beta, gamma and delta, return a random Z2Z4-additive 
code of type (alpha, beta; gamma, delta; kappa), where kappa is computed 
randomly. Note that there exists a Z2Z4-additive code of type (alpha, beta; 
gamma, delta; kappa) if and only if alpha, beta, gamma, delta, kappa >= 0, 
(alpha+beta)>0, k <= min(alpha,beta) and 0 < delta + gamma <= beta + kappa.
}
    requirege alpha, 0;
    requirege beta, 0;
    requirege gamma, 0;
    requirege delta, 0;
    require (alpha + beta) gt 0: 
        "Argument 1 plus argument 2 must be greater than 0";
    require (gamma + delta) gt 0: 
        "Argument 3 plus argument 4 must be greater than 0";
    require beta ge delta: 
        "Argument 2 must be greater than or equal to argument 4";
    require (alpha + beta) ge (gamma + delta): 
        "Argument 1 plus argument 2 must be greater than or equal to 
         argument 3 plus argument 4";
    
    kappa := Random(Max(0, delta + gamma - beta), Min(alpha, gamma));

    if IsZero(delta) and IsZero(gamma) and IsZero(kappa) then
        return Z2Z4AdditiveZeroCode(alpha, beta);
    end if;
    repeat 
        Mt := ZeroMatrix(Z4, 0, alpha + beta);
        if (kappa gt 0) then
            Ik := ScalarMatrix(Z4, kappa, 2);
            Tb := 2 * ChangeRing(RandomMatrix(Z2, kappa, alpha - kappa), Z4);
            T2_2 := 2 * ChangeRing(RandomMatrix
                   (Z2, kappa, beta - gamma + kappa - delta), Z4);
            Zero1 := ZeroMatrix(Z4, kappa, gamma - kappa + delta);
            Mt := HorizontalJoin(<Ik, Tb, T2_2, Zero1>);
        end if;
        if (gamma - kappa gt 0) then
            Zero1 := ZeroMatrix(Z4, gamma - kappa, alpha);
            Igk_2 := ScalarMatrix(Z4, gamma - kappa, 2);
            T1_2 := 2 * ChangeRing(RandomMatrix
                   (Z2, gamma - kappa, beta - gamma + kappa - delta), Z4);
            Zero2 := ZeroMatrix(Z4, gamma - kappa, delta);
            Mgk := HorizontalJoin(<Zero1, T1_2, Igk_2, Zero2>);
            Mt := VerticalJoin(Mt, Mgk);
        end if;
        if (delta gt 0) then
            Zero1 := ZeroMatrix(Z4, delta, kappa);
            Sb := 2 * ChangeRing(RandomMatrix(Z2, delta, alpha - kappa), Z4);
            Sq := RandomMatrix(Z4, delta, beta - gamma + kappa - delta);
            R := ChangeRing(RandomMatrix(Z2, delta, gamma - kappa), Z4);
            Id := ScalarMatrix(Z4, delta, 1);
            Md := HorizontalJoin(<Zero1, Sb, Sq, R, Id>);
            Mt := VerticalJoin(Mt, Md);
        end if;
    until (Mt ne ZeroMatrix(Z4, gamma + delta, alpha + beta));
    
    C := Z2Z4AdditiveCode(Mt, alpha);
    return C;
    
end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: RandomZ2Z4AdditiveCode                        */
/* Parameters: alpha, beta, gamma, delta, kappa                 */
/* Function description: Given five integers alpha, beta, gamma,*/
/*   delta and kappa, return a random Z2Z4-additive code of type*/
/*   (alpha, beta; gamma, delta; kappa). Note that there exists */
/*   a Z2Z4-additive code of type (alpha, beta; gamma, delta;   */
/*   kappa) if and only if alpha, beta, gamma, delta, kappa >= 0,*/ 
/*   (alpha+beta)>0, k <= min(alpha,beta) and                   */
/*   0 < delta + gamma <= beta + kappa.                         */
/* Input parameters description:                                */
/*   - alpha: integer with the binary length of the code        */
/*   - beta: integer with the quaternary length of the code     */
/*   - gamma: integer with the order two generators of the code */
/*   - delta: integer with the order four generators of the code*/
/*   - kappa: integer with the dimension the binary part of the */
/*            code                                              */
/* Output parameters description:                               */
/*   - A random Z2Z4-additive code                              */
/*                                                              */
/* Signature: (<RngIntElt> alpha, <RngIntElt> beta,             */
/*             <RngIntElt> gamma, <RngIntElt> delta,            */
/*                                <RngIntElt> kappa)-> Z2Z4Code */
/*                                                              */
/****************************************************************/
intrinsic RandomZ2Z4AdditiveCode(alpha::RngIntElt, beta::RngIntElt,
           gamma::RngIntElt, delta::RngIntElt, kappa::RngIntElt) -> Z2Z4Code
{
Given five integers alpha, beta, gamma, delta and kappa, return a random 
Z2Z4-additive code of type (alpha, beta; gamma, delta; kappa). 
Note that there exists a Z2Z4-additive code of type (alpha, beta; 
gamma, delta; kappa) if and only if alpha, beta, gamma, delta, kappa >= 0, 
(alpha+beta)>0, k <= min(alpha,beta) and 0 < delta + gamma <= beta + kappa.
}
    requirege alpha, 0;
    requirege beta, 0;
    requirege gamma, 0;
    requirege delta, 0;
    requirege kappa, 0;
    require (alpha + beta) gt 0: 
        "Argument 1 plus argument 2 must be greater than 0";
    require (gamma + delta) gt 0: 
        "Argument 3 plus argument 4 must be greater than 0";
    require beta ge delta: 
        "Argument 2 must be greater than or equal to argument 4";
    require (alpha + beta) ge (gamma + delta): 
        "Argument 1 plus argument 2 must be greater than or equal to 
         argument 3 plus argument 4";
    require kappa le Min(gamma, alpha): 
        "Argument 5 must be lesser than or equal to Min(argument1,argument3)";
    require (beta + kappa) ge (delta + gamma): 
        "Argument 2 plus argument 5 must be greater than or equal to 
         argument 3 plus argument 4";

    if IsZero(delta) and IsZero(gamma) and IsZero(kappa) then
        return Z2Z4AdditiveZeroCode(alpha, beta);
    end if;
    repeat 
        Mt := ZeroMatrix(Z4, 0, alpha + beta);
        if (kappa gt 0) then
            Ik := ScalarMatrix(Z4, kappa, 2);
            Tb := 2 * ChangeRing(RandomMatrix(Z2, kappa, alpha - kappa), Z4);
            T2_2 := 2 * ChangeRing(RandomMatrix
                    (Z2, kappa, beta - gamma + kappa - delta), Z4);
            Zero1 := ZeroMatrix(Z4, kappa, gamma - kappa + delta);
            Mt := HorizontalJoin(<Ik, Tb, T2_2, Zero1>);
        end if;
        if (gamma - kappa gt 0) then
            Zero1 := ZeroMatrix(Z4, gamma - kappa, alpha);
            Igk_2 := ScalarMatrix(Z4, gamma - kappa, 2);
            T1_2 := 2 * ChangeRing(RandomMatrix
                   (Z2, gamma - kappa, beta - gamma + kappa - delta), Z4);
            Zero2 := ZeroMatrix(Z4, gamma - kappa, delta);
            Mgk := HorizontalJoin(<Zero1, T1_2, Igk_2, Zero2>);
            Mt := VerticalJoin(Mt, Mgk);
        end if;
        if (delta gt 0) then
            Zero1 := ZeroMatrix(Z4, delta, kappa);
            Sb := 2 * ChangeRing(RandomMatrix(Z2, delta, alpha - kappa), Z4);
            Sq := RandomMatrix(Z4, delta, beta - gamma + kappa - delta);
            R := ChangeRing(RandomMatrix(Z2, delta, gamma - kappa), Z4);
            Id := ScalarMatrix(Z4, delta, 1);
            Md := HorizontalJoin(<Zero1, Sb, Sq, R, Id>);
            Mt := VerticalJoin(Mt, Md);
        end if;
    until (Mt ne ZeroMatrix(Z4, gamma + delta, alpha + beta));
    
    C:= Z2Z4AdditiveCode(Mt, alpha);
    return C;

end intrinsic;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////            BASIC NUMERICAL INVARIANTS                           ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/****************************************************************/
/*                                                              */
/* Function name: Length                                        */
/* Parameters: C                                                */
/* Function description: Given a Z2Z4-additive code C of type   */
/*   (alpha, beta; gamma, delta; kappa), return the length      */
/*   n = alpha + beta, and the sequence [alpha, beta] with the  */
/*   number of coordinates over Z2 and the number of            */
/*   coordinates over Z4, respectively.                         */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - The length of the Z2Z4-additive code C                   */
/*   - The sequence [alpha, beta]                               */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> RngIntElt, SeqEnum              */
/*                                                              */
/****************************************************************/
intrinsic Length(C::Z2Z4Code) -> RngIntElt, SeqEnum
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma, delta; kappa), 
return the length n = alpha + beta, and the sequence [alpha, beta] with the number 
of coordinates over Z2 and the number of coordinates over Z4, respectively.
}
    return C`Length, [C`Alpha, C`Length-C`Alpha];

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: BinaryLength                                  */
/* Parameters: C                                                */
/* Function description: Given a Z2Z4-additive code C of type   */
/*   (alpha, beta; gamma, delta; kappa), return the binary      */
/*   length nbin = alpha + 2 * beta of C, which corresponds to  */
/*   the length of the binary code Cbin = Phi(C), where Phi is  */
/*   the Gray map.                                              */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - The length of the the binary code Cbin = Phi(C)          */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> RngIntElt                       */
/*                                                              */
/****************************************************************/
intrinsic BinaryLength(C::Z2Z4Code) -> RngIntElt
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma, delta; kappa), 
return the binary length nbin = alpha + 2 * beta of C, which corresponds to 
the length of the binary code Cbin = Phi(C), where Phi is the Gray map.
}  
    return 2 * C`Length - C`Alpha;

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: PseudoDimension                               */
/* Parameters: C                                                */
/* Function description: The number of generators (which equals */
/*   the pseudo-dimension k) of the Z2Z4-additive code C as a   */
/*   quaternary linear code.                                    */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - The pseudodimension of the code C                        */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> RngIntElt                       */
/*                                                              */
/****************************************************************/
intrinsic PseudoDimension(C::Z2Z4Code) -> RngIntElt
{
The number of generators (which equals the pseudo-dimension k) of the
Z2Z4-additive code C as a quaternary linear code.
}
    return PseudoDimension(C`Code);

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: NumberOfGenerators                            */
/* Parameters: C                                                */
/* Function description: The number of generators (which equals */
/*   the pseudo-dimension k) of the                             */
/*   Z2Z4-additive code C as a quaternary linear code.          */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - The number of the generators of the code C               */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> RngIntElt                       */
/*                                                              */
/****************************************************************/
intrinsic NumberOfGenerators(C::Z2Z4Code) -> RngIntElt
{
The number of generators (which equals the pseudo-dimension k) of the
Z2Z4-additive code C as a quaternary linear code.
}
    return NumberOfGenerators(C`Code);

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Ngens                                         */
/* Parameters: C                                                */
/* Function description: The number of generators (which equals */
/*   the pseudo-dimension k) of the Z2Z4-additive code C as a   */
/*   quaternary linear code.                                    */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - The number of the generators of the code C               */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> RngIntElt                       */
/*                                                              */
/****************************************************************/
intrinsic Ngens(C::Z2Z4Code) -> RngIntElt
{
The number of generators (which equals the pseudo-dimension k) of the
Z2Z4-additive code C as a quaternary linear code.
}
    return Ngens(C`Code);

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4Type                                      */
/* Parameters: C                                                */
/* Function description: Given a Z2Z4-additive code C, return a */
/*   sequence with the parameters [alpha,beta,gamma,delta,kappa],*/
/*   that is, the type of the code.                             */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - The parameters [alpha, beta, gamma, delta, kappa]        */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> RngIntElt                       */
/*                                                              */
/****************************************************************/
intrinsic Z2Z4Type(C::Z2Z4Code) -> SeqEnum
{
Given a Z2Z4-additive code C, return a sequence with the parameters
[alpha,beta,gamma,delta,kappa], that is, the type of the code.
}
    Rs, f := StandardForm(C`Code);
    delta, gamma := Z4Type(C`Code);
    gens := [ Eltseq(Rs.i@@f)[1..C`Alpha] : i in [delta + 1..gamma + delta] ];
    kappa := Dimension(sub<RSpace(Z4, C`Alpha) | gens>);

    return [C`Alpha, C`Length-C`Alpha, gamma, delta, kappa];

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: #                                             */
/* Parameters: C                                                */
/* Function description: Given a Z2Z4-additive code C of type   */
/*   (alpha, beta; gamma, delta; kappa), return the number of   */
/*   codewords belonging to C, that is 2^gamma 4^delta.         */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - The number of codewords of C                             */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> RngIntElt                       */
/*                                                              */
/****************************************************************/
intrinsic '#'(C::Z2Z4Code) -> RngIntElt
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma, delta; kappa), return 
the number of codewords belonging to C, that is 2^gamma 4^delta.
}
    return #C`Code;

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: InformationRate                               */
/* Parameters: C                                                */
/* Function description: Given a Z2Z4-additive code C of type   */
/*   (alpha, beta; gamma, delta; kappa), return the information */
/*   rate of C, that is the ratio                               */
/*   (gamma + 2*delta)/(alpha + 2*beta).                        */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - The information rate of C                                */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> FldRatElt                       */
/*                                                              */
/****************************************************************/
intrinsic InformationRate(C::Z2Z4Code) -> FldRatElt
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma, delta; kappa), 
return the information rate of C, that is the ratio (gamma + 2*delta)/
(alpha + 2*beta).
}
    return Log(2, #C`Code) / BinaryLength(C);

end intrinsic;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////                THE CODE SPACE                                   ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/****************************************************************/
/*                                                              */
/* Function name: Alphabet                                      */
/* Parameters: C                                                */
/* Function description: The underlying finite ring             */
/*   (or alphabet) R of the Z2Z4-additive code C.               */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - The finite ring Z4                                       */
/*                                                              */
/* Signature: (<CodeFld> C) -> RngFin                           */
/*                                                              */
/****************************************************************/
/* intrinsic Alphabet(C::Z2Z4Code) -> RngFin
/*{
/*The underlying finite ring (or alphabet) Z4 of the Z2Z4-additive code C.
/*}
/*    return Z4;
/*end intrinsic;*/

/****************************************************************/
/*                                                              */
/* Function name: AmbientSpace                                  */
/* Parameters: C                                                */
/* Function description: Given a Z2Z4-additive code C of type   */
/*   (alpha,beta;gamma,delta;kappa), return the ambient space   */
/*   of C, i.e., the Z4-submodule of Z4^(alpha+beta) isomorphic */
/*   to Z2^alpha x Z4^beta, where C is contained after replacing*/
/*   the ones in the first alpha coordinates by twos. It also   */
/*   returns the generic space Z2^alpha x Z4^beta as a          */
/*   cartesian product.                                         */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - The cartesian product Z2^alpha x Z4^beta                 */
/*   - The Z4-space isomorphic to Z2^alpha x Z4^beta            */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> SetCart, ModTupRng              */
/*                                                              */
/****************************************************************/
intrinsic AmbientSpace(C::Z2Z4Code) -> ModTupRng, SetCart
{
Given a Z2Z4-additive code C of type (alpha,beta;gamma,delta;kappa), return the 
ambient space of C, i.e., the Z4-submodule of Z4^(alpha+beta) isomorphic to 
Z2^alpha x Z4^beta, where C is contained after replacing the ones in the first 
alpha coordinates by twos. It also returns the generic space Z2^alpha x Z4^beta 
as a cartesian product.
}
    alpha := C`Alpha;
    beta := C`Length - alpha;
    V2 := RSpace(Z2, alpha);
    V4 := RSpace(Z4, beta);
    V := CartesianProduct(V2, V4);
    U := Z2Z4AdditiveUniverseCode(alpha, beta);
    R := RSpace(U`Code);
    
    return R, V;

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Generic                                       */
/* Parameters: C                                                */
/* Function description: Given a Z2Z4-additive code C of type   */
/*   (alpha, beta; gamma, delta; kappa), return the generic     */
/*   Z2Z4-additive code of type (alpha, beta; alpha, beta; alpha)*/
/*   in which C is contained.                                   */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - The generic (n, 4^n, 1) Z2Z4-additive code               */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> Z2Z4Code                        */
/*                                                              */
/****************************************************************/
intrinsic Generic(C::Z2Z4Code) -> Z2Z4Code
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma, delta; kappa), 
return the generic Z2Z4-additive code of type (alpha, beta; alpha, beta; alpha)
in which C is contained.
}
    return Z2Z4AdditiveUniverseCode(C`Alpha, C`Length-C`Alpha);

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: Name                                          */
/* Parameters: C, i                                             */
/* Function description: Given a Z2Z4-Additive code C of type   */
/*   (alpha, beta; gamma, delta; kappa), and a positive integer */
/*   i, return the i-th generator of C as a quaternary linear   */
/*   code, where the ones in the first alpha coordinates are    */
/*   represented by twos.                                       */
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/*   - i: A positive integer                                    */
/* Output parameters description:                               */
/*   - The i-th generator of C                                  */
/*                                                              */
/* Signature: (<Z2Z4Code> C, <RngIntElt> i) -> ModTupRngElt     */
/*                                                              */
/****************************************************************/
intrinsic Name(C::Z2Z4Code, i::RngIntElt) -> ModTupRngElt
{
Given a Z2Z4-Additive code C of type (alpha, beta; gamma, delta; kappa), and 
a positive integer i, return the i-th generator of C as a quaternary linear code, 
where the ones in the first alpha coordinates are represented by twos.
}
    requirerange i, 1, Ngens(C);

    return C`Code.i;

end intrinsic;

intrinsic '.'(C::Z2Z4Code, i::RngIntElt) -> ModTupRngElt
{
Given a Z2Z4-Additive code C of type (alpha, beta; gamma, delta; kappa), and 
a positive integer i, return the i-th generator of C as a quaternary linear code, 
where the ones in the first alpha coordinates are represented by twos.
}
    requirerange i, 1, Ngens(C);

    return C`Code.i;

end intrinsic;

/************************************************************/
/*                                                          */
/* Function name: Set                                       */
/* Parameters: C                                            */
/* Function description: Given a Z2Z4-additive code C of    */
/*   type (alpha, beta; gamma, delta; kappa), return the    */
/*   set containing all its codewords, where the ones in    */
/*   the first alpha coordinates are represented by twos.   */
/* Input parameters description:                            */
/*   - C: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - The set with all codewords of the code C             */
/*                                                          */
/* Signature: (<Z2Z4Code> C) -> SetEnum                     */
/*                                                          */
/************************************************************/
intrinsic Set(C::Z2Z4Code) -> SetEnum
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma, delta; kappa), 
return the set containing all its codewords, where the ones in the first alpha 
coordinates are represented by twos.
}
    return Set(C`Code);

end intrinsic;

/************************************************************/
/*                                                          */
/* Function name: Z2Z4Set                                   */
/* Parameters: C                                            */
/* Function description: Given a Z2Z4-additive code C of    */
/*   type (alpha, beta; gamma, delta; kappa), return the    */
/*   set containing all its codewords as tuples in the      */
/*   cartesian product set Z2^alpha x Z4^beta.              */
/* Input parameters description:                            */
/*   - C: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - The set with all codewords of the code C             */
/*                                                          */
/* Signature: (<Z2Z4Code> C) -> SetEnum                     */
/*                                                          */
/************************************************************/
intrinsic Z2Z4Set(C::Z2Z4Code) -> SetEnum
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma, delta; kappa), 
return the set containing all its codewords as tuples in the cartesian product 
set  Z2^alpha x Z4^beta.
}
    return FromZ4toZ2Z4(Set(C`Code), C`Alpha);

end intrinsic;

/************************************************************/
/*                                                          */
/* Function name: Basis                                     */
/* Parameters: C                                            */
/* Function description: The basis for the Z2Z4-additive    */
/*   code C of type (alpha, beta; gamma, delta; kappa), as  */
/*   a quaternary linear code, returned as a sequence of    */
/*   elements of C, where the ones in the first alpha       */
/*   coordinates are represented by twos.                   */
/* Input parameters description:                            */
/*   - C: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - The basis of C as a sequence of vectors              */
/*                                                          */
/* Signature: (<Z2Z4Code> C) -> SeqEnum                     */
/*                                                          */
/************************************************************/
intrinsic Basis(C::Z2Z4Code) -> SeqEnum
{
The basis for the Z2Z4-additive code C of type (alpha, beta; gamma, delta; 
kappa), as a quaternary linear code, returned as a sequence of elements of 
C, where the ones in the first alpha coordinates are represented by twos.
}
    return Basis(C`Code);

end intrinsic;

/************************************************************/
/*                                                          */
/* Function name: Generators                                */
/* Parameters: C                                            */
/* Function description: The generators for the             */
/*   Z2Z4-additive code C of type (alpha, beta; gamma, delta;*/
/*   kappa), as a quaternary linear code, returned as a set */
/*   of elements of C, where the ones in the first          */
/*   alpha coordinates are represented by twos.             */
/* Input parameters description:                            */
/*   - C: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - The generators of C as a set of vectors              */
/*                                                          */
/* Signature: (<Z2Z4Code> C) -> SeqEnum                     */
/*                                                          */
/************************************************************/
intrinsic Generators(C::Z2Z4Code) -> SetEnum
{
The generators for the Z2Z4-additive code C of type (alpha, beta; gamma, delta; 
kappa), as a quaternary linear code, returned as a set of elements of C, 
where the ones in the first alpha coordinates are represented by twos.
}
    return Generators(C`Code);

end intrinsic;

/************************************************************/
/*                                                          */
/* Function name: GeneratorMatrix                           */
/* Parameters: C                                            */
/* Function description: The unique generator matrix for    */
/*   the Z2Z4-additive code C of type (alpha, beta; gamma,  */
/*   delta; kappa), as a quaternary linear code,            */
/*   corresponding to the Howell form. The ones in the      */
/*   first alpha coordinates are represented by twos.       */
/* Input parameters description:                            */
/*   - C: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - The unique generator matrix of C in Howell form      */
/*                                                          */
/* Signature: (<Z2Z4Code> C) -> ModMatRngElt                */
/*                                                          */
/************************************************************/
intrinsic GeneratorMatrix(C::Z2Z4Code) -> ModMatRngElt
{
The unique generator matrix for the Z2Z4-additive code C of type (alpha, beta; 
gamma, delta; kappa), as a quaternary linear code, corresponding to the Howell 
form. The ones in the first alpha coordinates are represented by twos.
}
    return GeneratorMatrix(C`Code);

end intrinsic;

/************************************************************/
/*                                                          */
/* Function name: OrderTwoGenerators                        */
/* Parameters: C                                            */
/* Function description: The gamma generators of order two  */
/*   of the Z2Z4-additive code C of type (alpha, beta;      */
/*   gamma, delta; kappa), returned as a set of elements of */
/*   C, where the ones in the first alpha coordinates are   */
/*   represented by twos.                                   */
/* Input parameters description:                            */
/*   - C: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - The gamma generators of order two of the code C      */
/*                                                          */
/* Signature: (<Z2Z4Code> C) -> SetEnum                     */
/*                                                          */
/************************************************************/
intrinsic OrderTwoGenerators(C::Z2Z4Code) -> SetEnum
{
The gamma generators of order two of the Z2Z4-additive code C of type (alpha, 
beta; gamma, delta; kappa), returned as a set of elements of C, where the ones in 
the first alpha coordinates are represented by twos.
}
    Order2, _ := OrderTwoFourGenerators(C);
    return Order2;

end intrinsic;

/************************************************************/
/*                                                          */
/* Function name: OrderFourGenerators                       */
/* Parameters: C                                            */
/* Function description: The delta generators of order four */
/*   of the Z2Z4-additive code C of type (alpha, beta;      */
/*   gamma, delta; kappa), returned as a set of elements of */
/*   C, where the ones in the first alpha coordinates are   */
/*   represented by twos.                                   */
/* Input parameters description:                            */
/*   - C: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - The delta generators of order four of the code C     */
/*                                                          */
/* Signature: (<Z2Z4Code> C) -> SetEnum                     */
/*                                                          */
/************************************************************/
intrinsic OrderFourGenerators(C::Z2Z4Code) -> SetEnum
{
The delta generators of order four of the Z2Z4-additive code C of type 
(alpha, beta; gamma, delta; kappa), returned as a set of elements of C, where 
the ones in the first alpha coordinates are represented by twos.
}
    _, Order4 := OrderTwoFourGenerators(C);
    return Order4; 

end intrinsic;

/************************************************************/
/*                                                          */
/* Function name: OrderTwoSubcodeGenerators                 */
/* Parameters: C                                            */
/* Function description: The gamma + delta generators of    */
/*   the Z2Z4-additive subcode Cb which contains all order  */
/*   two codewords of the Z2Z4-additive code C of type      */
/*   (alpha, beta; gamma, delta; kappa), returned as a set  */
/*   of elements of Cb, where the ones in the first alpha   */
/*   coordinates are represented by twos.                   */
/* Input parameters description:                            */
/*   - C: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - The gamma + delta generators of the code C           */
/*                                                          */
/* Signature: (<Z2Z4Code> C) -> SetEnum                     */
/*                                                          */
/************************************************************/
intrinsic OrderTwoSubcodeGenerators(C::Z2Z4Code) -> SetEnum
{
The gamma + delta generators of the Z2Z4-additive subcode Cb which contains all
order two codewords of the Z2Z4-additive code C of type (alpha, beta; gamma, 
delta; kappa), returned as a set of elements of Cb, where the ones in the 
first alpha coordinates are represented by twos.
}
    Order2, Order4 := OrderTwoFourGenerators(C);
    return Order2 join {2*v : v in Order4};

end intrinsic;

/************************************************************/
/*                                                          */
/* Function name: MinRowsGeneratorMatrix                    */
/* Parameters: C                                            */
/* Function description: A generator matrix for the         */
/*   Z2Z4-additive code C of type (alpha, beta; gamma,      */
/*   delta; kappa), with the minimum number of rows, that   */
/*   is with gamma + delta rows: gamma rows of order two    */
/*   and delta rows of order four. It also returns the      */
/*   parameters gamma and delta. The ones in the first      */
/*   alpha coordinates are represented by twos.             */
/* Input parameters description:                            */
/*   - C: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - A generator matrix for C, with the minimum number    */
/*     of rows                                              */
/*   - The parametres gamma and delta                       */
/*                                                          */
/* Signature: (<Z2Z4Code> C) -> ModMatRngElt, RngIntElt,    */
/*                                            RngIntElt     */
/*                                                          */
/************************************************************/
intrinsic MinRowsGeneratorMatrix(C::Z2Z4Code) -> ModMatRngElt, RngIntElt, RngIntElt
{
A generator matrix for the Z2Z4-additive code C of type
(alpha, beta; gamma, delta; kappa), with the minimum number of rows, 
that is with gamma + delta rows: gamma rows of order two and delta rows of 
order four. It also returns the parameters gamma and delta. The ones in the 
first alpha coordinates are represented by twos.
}
    Order2, Order4 := OrderTwoFourGenerators(C);
    gamma := #Order2;
    delta := #Order4;

    if (gamma eq 0) and (delta eq 0) then
        return GeneratorMatrix(C`Code), 0, 0;
    else
        return Matrix(Setseq(Order2) cat Setseq(Order4)), gamma, delta;
    end if;

end intrinsic;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////                            THE DUAL SPACE                       ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/************************************************************/
/*                                                          */
/* Function name: Z2Z4DualType                              */
/* Parameters: C                                            */
/* Function description: Given a Z2Z4-additive code C,      */
/*   return a sequence with the parameters [alpha, beta,    */
/*   gamma', delta', kappa'], that is the type of the       */
/*   additive dual code of C. If C is of type (alpha, beta; */
/*   gamma, delta; kappa), then its additive dual code is   */
/*   of type (alpha, beta; alpha + gamma - (2*kappa),       */
/*   beta - gamma - delta + kappa; alpha - kappa).          */
/* Input parameters description:                            */
/*   - C: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - A sequence with the type of the additive dual code   */
/*     of C                                                 */
/*                                                          */
/* Signature: (<Z2Z4Code> C) -> SeqEnum                     */
/*                                                          */
/************************************************************/
intrinsic Z2Z4DualType(C::Z2Z4Code) -> SeqEnum
{
Given a Z2Z4-additive code C, return a sequence with the parameters
[alpha, beta, gamma', delta', kappa'], that is the type of the additive dual 
code of C. If C is of type (alpha, beta; gamma, delta; kappa), then its 
additive dual code is of type (alpha, beta; alpha + gamma - (2*kappa), 
beta - gamma - delta + kappa; alpha - kappa).
}
    type := Z2Z4Type(C);
    return [type[1], type[2], type[1]+type[3]-2*type[5],
            type[2]-type[3]-type[4]+type[5], type[1]-type[5]];

end intrinsic;

/************************************************************/
/*                                                          */
/* Function name: Dual                                      */
/* Parameters: C                                            */
/* Function description: Given a Z2Z4-additive code C,      */
/*   return the additive dual code of C.                    */
/* Input parameters description:                            */
/*   - C: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - The additive dual code of C                          */
/*                                                          */
/* Signature: (<Z2Z4Code> C) -> Z2Z4Code                    */
/*                                                          */
/************************************************************/
intrinsic Dual(C::Z2Z4Code) -> Z2Z4Code
{
Given a Z2Z4-additive code C, return the additive dual code of C.
}
    if (zero_Code(C)) then
        return Z2Z4AdditiveUniverseCode(C`Alpha, C`Length - C`Alpha);
    elif Length(C) eq C`Alpha then
        return Z2Z4AdditiveCode(Dual(LinearBinaryCode(C)));
    elif IsZero(C`Alpha) then
        return Z2Z4AdditiveCode(Dual(C`Code));
    end if;

    M := GeneratorMatrix(C`Code);
    Z2Z4ChangeMatrixZ4toZ2(~M, C`Alpha);
    M := VerticalJoin(HorizontalJoin(ScalarMatrix(Z4, C`Alpha, 2),
    ZeroMatrix(Z4, C`Alpha, Ncols(M)-C`Alpha)), M);
    D := Z2Z4AdditiveCode(Dual(LinearCode(M)), C`Alpha);
    
    return D;

end intrinsic;

/************************************************************/
/*                                                          */
/* Function name: ParityCheckMatrix                         */
/* Parameters: C                                            */
/* Function description: The unique parity-check matrix for */
/*   the Z2Z4-additive code C of type (alpha, beta; gamma,  */
/*   delta; kappa), that is the unique generator matrix for */
/*   the additive dual code of C as a quaternary linear code,*/
/*   corresponding to the Howell form. The ones in the first*/
/*   alpha coordinates are represented by twos.             */
/* Input parameters description:                            */
/*   - C: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - The unique generator matrix for the additive dual    */
/*     code of C in Howell form                             */
/*                                                          */
/* Signature: (<Z2Z4Code> C) -> ModMatRngElt                */
/*                                                          */
/************************************************************/
intrinsic ParityCheckMatrix(C::Z2Z4Code) -> ModMatRngElt
{
The unique parity-check matrix for the Z2Z4-additive code C of type (alpha, 
beta; gamma, delta; kappa), that is the unique generator matrix for the additive 
dual code of C as a quaternary linear code, corresponding to the Howell form. 
The ones in the first alpha coordinates are represented by twos.
}   
    M := GeneratorMatrix(Dual(C));  
    return M;

end intrinsic;

/************************************************************/
/*                                                          */
/* Function name: MinRowsParityCheckMatrix                  */
/* Parameters: C                                            */
/* Function description: A parity-check matrix for the      */
/*   Z2Z4-additive code C of type (alpha, beta; gamma,      */
/*   delta; kappa), that is a generator matrix for the      */
/*   additive dual code of C, with the minimum number of    */
/*   rows, that is with gamma + delta rows: gamma rows of   */
/*   order two and delta rows of order four. The ones in    */
/*   the first alpha coordinates are represented by twos.   */
/* Input parameters description:                            */
/*   - C: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - A generator matrix for the additive dual code of C,  */
/*     with the minimum number of rows                      */
/*                                                          */
/* Signature: (<Z2Z4Code> C) -> ModMatRngElt                */
/*                                                          */
/************************************************************/
intrinsic MinRowsParityCheckMatrix(C::Z2Z4Code) -> ModMatRngElt
{
A parity-check matrix for the Z2Z4-additive code C of type (alpha, beta; gamma, 
delta; kappa), that is a generator matrix for the additive dual code of C, with 
the minimum number of rows, that is with gamma + delta rows: gamma rows of order 
two and delta rows of order four. The ones in the first alpha coordinates are 
represented by twos.
}
    return MinRowsGeneratorMatrix(Dual(C));

end intrinsic;

/************************************************************/
/*                                                          */
/* Function name: Hull                                      */
/* Parameters: C                                            */
/* Function description: Given a Z2Z4-additive code C,      */
/*   return the Hull of C, which is the intersection        */
/*   between itself and its additive dual.                  */ 
/* Input parameters description:                            */
/*   - C: A Z2Z4-additive code                              */
/* Output parameters description:                           */
/*   - The Hull of C                                        */
/*                                                          */
/* Signature: (<Z2Z4Code> C) -> Z2Z4Code                    */
/*                                                          */
/************************************************************/
intrinsic Hull(C::Z2Z4Code) -> Z2Z4Code
{
Given a Z2Z4-additive code C, return the Hull of C, which is the intersection 
between itself and its additive dual.
}
    D := C meet Dual(C);
    return D;

end intrinsic;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////                THE GRAY MAP                                     ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/******************************************************************************
            GrayMap function    
******************************************************************************/

Z2Z4GraySeq := [[0,0], [0,1], [1,1], [1,0]];
Z2Z4GrayInvSeq := [0, 3, 1, 2];
Z := IntegerRing();

function Z2Z4GrayVec(v, alpha)
    return Vector(Z2,([Z2Z4AlphaFromZ4toZ2(v[i]) : i in [1..alpha]]) cat 
           (&cat[Z2Z4GraySeq[Z!v[i] + 1] : i in [alpha + 1..Degree(v)]]));
end function;

function Z2Z4GrayVecInv(y, alpha, beta)
    alphaseq := [Z2Z4AlphaFromZ2toZ4(Z2!y[i]) : i in [1..alpha]];
    betaseq := [y[i] : i in [alpha + 1..2 * beta + alpha]];
    return alphaseq cat [Z2Z4GrayInvSeq[Z!betaseq[2*i - 1] + 
           2*Z!betaseq[2*i] + 1] : i in [1..beta]];
end function;

function Z2Z4GrayMapCodomain(R, V)
    a, len := Length(R);
    f := map<R`Code -> V |
        x :-> Z2Z4GrayVec(x,R`Alpha),
        y :-> R`Code!Z2Z4GrayVecInv(y,len[1],len[2])
    >;
    return f;
end function;

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
intrinsic GrayMap(C::Z2Z4Code) -> Map
{
Given a Z2Z4-additive code C, return the Gray map for C. This is the
map Phi from C to Phi(C), as defined above.
}
    return Z2Z4GrayMapCodomain(C, VectorSpace(GF(2), BinaryLength(C)));

end intrinsic;

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
intrinsic GrayMapImage(C::Z2Z4Code) -> SeqEnum
{
Given a Z2Z4-additive code C, return the image of C under the Gray
map as a sequence of vectors in Z2^(alpha+2*beta). As the resulting image may 
not be a binary linear code, a sequence of vectors is returned rather than a
code.
}    
    f := GrayMap(C); 
    return [f(x): x in C`Code];

end intrinsic;

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
intrinsic HasLinearGrayMapImage(C::Z2Z4Code) -> BoolElt, CodeLinFld, Map
{
Given a Z2Z4-additive code C, return true if and only if the image of
C under the Gray map is a binary linear code. If so, the function also
returns the image B as a binary linear code, together with the bijection
Phi: C -> B.
}
    Order2, Order4 := OrderTwoFourGenerators(C);
    binaryGen :=  [ Z2Z4GrayVec(v, C`Alpha) : v in Order2 ];
    if (#Order4 ne 0) then
        Order4seq := Setseq(Order4);
        delta := #Order4seq;
        V := AmbientSpace(C`Code);
        n := C`Length;
        for j in [1..delta-1] do
            v := Order4seq[j];
            for k in [j+1..delta] do
                w := Order4seq[k];
                if 2*(v*w) notin C`Code then
                    return false, _, _;
                end if;
            end for;
        end for;
        binaryGen := binaryGen cat &cat[ [Z2Z4GrayVec(v, C`Alpha),
                                Z2Z4GrayVec(-v, C`Alpha)] : v in Order4 ];    
    end if;
    
    if (#binaryGen eq 0) then 
        B := ZeroCode(F2, BinaryLength(C));       
        return true, B, Z2Z4GrayMapCodomain(C, B);
    else
        B := LinearCode(ChangeRing(Matrix(binaryGen), F2));
        return true, B, Z2Z4GrayMapCodomain(C, B);
    end if;

end intrinsic;

