///////////////////////////////////////////////////////////////////////////////
/////////       Copyright 2007-2022 Bernat Gaston, Jaume Pujol          ///////
/////////         and Merce Villanueva (with contributions of           ///////
/////////         Adrian Torres)                                        ///////
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
/* File name: Z2Z4StandardForm.m                             */
/*                                                           */
/* Comment: Package developed within the CCSG group          */
/*                                                           */
/* Authors: Bernat Gaston, Jaume Pujol and Merce Villanueva  */
/*          (with contributions of Adrian Torres)            */
/*                                                           */
/* Revision version and last date: v1.0   29-06-2012         */
/*                                 v1.1   06-11-2014         */
/*                                 v1.2   10-02-2018         */
/*                                 v2.0   01-02-2019         */
/*                                 v2.1   30-05-2022         */
/* History: Improvements of Z2Z4StandardForm function.       */
/*          Split from Z2Z4AdditiveCodes.m file.             */
/*          v1.1 Version function added                      */
/*          v2.0 User defined user type by J. Pujol          */
/*          v2.1 Some functions updated and added to use them*/
/*            for the minimum weight computation by A. Torres*/ 
/*                                                           */
/*************************************************************/
//Uncomment freeze when package finished
freeze;

intrinsic Z2Z4StandardForm_version() -> SeqEnum
{Return the current version of this package.}
    version := [2,1];
    return version;
end intrinsic;

//needs Z2Z4AdditiveCode.m file

import "Z2Z4AdditiveCodes.m" : Z4, Z2, F2, over_Z4, over_Z2, over_F2,
                               zero_Code, NUMBER_OF_PARAMETERS;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////                THE STANDARD FORM                                ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


/******************************************************************************
                DiagonalizeKappaMatrix function
******************************************************************************/
/* Function first developed by B. Gaston and updated by A. Torres            */

function DiagonalizeKappaMatrix(M, gamma, kappa)
    ncols := Ncols(M);
    mrows := Nrows(M);
    colsSym := Sym(ncols); 
    P := colsSym!1;
    Q := IdentityMatrix(Z4, mrows);

    for j := 1 to kappa do
        found := false;
        columnPivot := j;
        while not found do
            rowPivot := j;                        
            //search for non zero element in columnPivot column
            while (rowPivot lt gamma) and (M[rowPivot][columnPivot] eq 0) do
                 rowPivot := rowPivot + 1;
            end while;
             
            if (M[rowPivot][columnPivot] ne 0) then
                // swap rows and columns if non zero element found
                SwapRows(~M, j, rowPivot);
                SwapRows(~Q, j, rowPivot);
                if (j ne columnPivot) then                    
                    SwapColumns(~M, j, columnPivot);                        
                    P := P*colsSym!(j, columnPivot);                    
                end if;
                found := true;
            else
                // search for non zero element in next column
                columnPivot := columnPivot + 1;                                        
            end if;                        
        end while;
        //diagonalize j column
        for i := j-1 to 1 by -1 do
            if (M[i][j] ne 0) then
                   M[i] := M[i] - M[j];
                   Q[i] := Q[i] - Q[j];
            end if;
        end for;
        for i := j+1 to mrows do
            if (M[i][j] ne 0) then
                M[i] := M[i] - M[j];
                Q[i] := Q[i] - Q[j];
            end if;
        end for;
    end for;
    return M, P, Q;
end function;



/******************************************************************************
                        DiagonalizeDeltaMatrix function
******************************************************************************/
/* Function first developed by B. Gaston and updated by A. Torres            */

function DiagonalizeDeltaMatrix(M, alpha, delta)
    ncols := Ncols(M);
    mrows := Nrows(M);
    nc := ncols;
    mc := mrows;
    colsSym := Sym(ncols);
    P := colsSym!1;
    Q := IdentityMatrix(Z4, mrows);

    for j := ncols to ncols-delta+1 by -1 do
        if (M[mc][nc] eq 0) or (M[mc][nc] eq 2) then
            m2 := mc;
            while (m2 gt mrows-delta) and ((M[m2][nc] eq 0) or (M[m2][nc] eq 2)) do
                 m2 := m2-1;
            end while;
            if (m2 gt mrows-delta) then
                SwapRows(~M, mc, m2);
                SwapRows(~Q, mc, m2);
            else
                n2 := nc;
                while ((M[mc][n2] eq 0) or (M[mc][n2] eq 2)) and (n2 gt alpha) do
                    n2 := n2-1;
                end while;
                if (n2 gt alpha) and (n2 ne nc) then
                    SwapColumns(~M, nc, n2);
                    P := P*colsSym!(nc, n2);
                end if;
            end if;
        end if;

        if (M[mc][nc] eq 3) then
            M[mc] := M[mc]*(-1);
            Q[mc] := Q[mc]*(-1);
        end if;
        for i := mc-1 to 1 by -1 do
            if (M[i][nc] ne 0) then
                   Q[i] := Q[i]-Q[mc]*M[i][nc];
                   M[i] := M[i]-M[mc]*M[i][nc];
            end if;
        end for;
        for i := mc+1 to mrows do
            if (M[i][nc] ne 0) then
                Q[i] := Q[i]-Q[mc]*M[i][nc];
                M[i] := M[i]-M[mc]*M[i][nc];
            end if;
        end for;

        mc := mc-1;
        nc := nc-1;
    end for;
    
    return M, P, Q;
end function;



/******************************************************************************
                        DiagonalizeGammaKappaMatrix function
******************************************************************************/
/* Function first developed by B. Gaston and updated by A. Torres            */

function DiagonalizeGammaKappaMatrix(M, alpha, gamma, delta, kappa)
    ncols := Ncols(M);
    mrows := Nrows(M);
    nc := ncols-delta;
    mc := gamma;
    colsSym := Sym(ncols);
    P := colsSym!1;
    Q := IdentityMatrix(Z4, mrows);
    for j := nc to nc-gamma+kappa+1 by -1 do

        if (M[mc][nc] eq 0) then
            m2 := mc;
            while (m2 gt 0) and (m2 gt nc-gamma+kappa) and (M[m2][nc] eq 0) do
                 m2 := m2-1;
            end while;
            if (m2 gt nc-gamma+kappa) and (m2 ne 0) then
                SwapRows(~M, mc, m2);
                SwapRows(~Q, mc, m2);
            else
                n2 := nc;
                while (M[mc][n2] eq 0) and (n2 gt alpha) do
                    n2 := n2-1;
                end while;
                if (n2 gt alpha) and (n2 ne nc) then
                    SwapColumns(~M, nc, n2);
                    P := P*colsSym!(nc, n2);
                end if;
            end if;
        end if;

        for i := mc-1 to 1 by -1 do
            if (M[i][nc] ne 0) then
                   M[i] := M[i]-M[mc];
                   Q[i] := Q[i]-Q[mc];
            end if;
        end for;
        for i := mc+1 to gamma do
            if (M[i][nc] ne 0) then
                M[i] := M[i]-M[mc];
                Q[i] := Q[i]-Q[mc];
            end if;
        end for;
        for i := gamma+1 to mrows do
            if (M[i][nc] ne 0) and (M[i][nc] ne 1) then
                M[i] := M[i]-M[mc];
                Q[i] := Q[i]-Q[mc];
            end if;
        end for;

        mc := mc-1;
        nc := nc-1;        
    end for;
    
    return M, P, Q;
end function;



/******************************************************************************
            Z2Z4StandardForm function
******************************************************************************/
function Z2Z4StandardFormInfo(R)
    M := MinRowsGeneratorMatrix(R);
    type := Z2Z4Type(R);
    n := Ncols(M);
    Mk, Pk := DiagonalizeKappaMatrix(M, type[3], type[5]);  
    Md, Pd := DiagonalizeDeltaMatrix(Mk, type[1], type[4]);  
    Mgk, Pgk := DiagonalizeGammaKappaMatrix(Md, type[1], type[3], type[4], type[5]);
    P := Pk*Pd*Pgk; 
    return Mgk, P;
end function;

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
/* Signature:(<Z2Z4Code> C) -> Z2Z4Code, Map, ModMatRngElt, */
/*                             GrpPermElt                   */
/*                                                          */
/************************************************************/
intrinsic StandardForm(C::Z2Z4Code) -> Z2Z4Code, Map, ModMatRngElt, GrpPermElt
{
Given any Z2Z4-additive code C, return a permutation-equivalent Z2Z4-
additive code Csf in standard form, together with the corresponding
isomorphism from C to Csf, the generator matrix in standard form, and
the coordinate permutation used to define the isomorphism.
}
    require not(zero_Code(C)): "Code C can not be the zero code";

    S, p := Z2Z4StandardFormInfo(C);
    RC := Z2Z4AdditiveCode(S, C`Alpha);
    f := map<C`Code -> RC`Code | v :-> v^p, w :-> w^(p^-1)>;

    return RC, f, S, p;
end intrinsic;



/******************************************************************************
            IsStandardFormMatrix function
******************************************************************************/


// M is a matrix over Z4 to check whether all entrances are in {0,1} or not
function IsBinaryMatrix(M)
    n := Ncols(M);
    m := Nrows(M);
    for i := 1 to m do
        for j := 1 to n do
             if (M[i][j] gt 1) then
                return false;
               end if;
        end for;
    end for;
    return true;
end function;

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
/* Signature:(<AlgMatElt> M, <SeqEnum> L) -> BoolElt        */
/*                                                          */
/************************************************************/
intrinsic IsStandardFormMatrix(M::AlgMatElt, L::SeqEnum) -> BoolElt
{
Return true if and only if the matrix M over Z4, where the ones in the first alpha 
coordinates are represented by twos, is a generator matrix of a Z2Z4-additive 
code of type (alpha, beta; gamma, delta; kappa) in standard form.
}
    require #L eq NUMBER_OF_PARAMETERS: "Argument 2 must contain 5 parameters";
    require ElementType(L) eq RngIntElt: "Argument 2 must contain integer numbers"; 
    require (Type(BaseRing(M)) eq RngIntRes) and (BaseRing(M) eq Z4): 
            "Argument 1 is not a matrix over Z4";

    alpha := L[1];
    beta := L[2];
    gamma := L[3];
    delta := L[4];
    kappa := L[5];
    n := Ncols(M);
    m := Nrows(M);
    
    cond := (alpha ge 0) and (beta ge 0) and (gamma ge 0) and (delta ge 0) and
            (kappa ge 0) and (alpha+beta gt 0) and (alpha+beta eq n) and 
            (gamma+delta gt 0) and (gamma+delta eq m) and (beta ge delta) and
            (alpha+beta ge gamma+delta) and (kappa le Min(gamma,alpha)) and
            (beta+kappa ge delta+gamma);
    if (not cond) then 
        return false;
    end if;

    Mk := Submatrix(M, 1, 1, kappa, kappa); 
    
    Mkzero1 := Submatrix(M, kappa+1,1, gamma-kappa,kappa);
    Mkzero2 := Submatrix(M, gamma+1, 1, delta, kappa);
    boolk := IsOne(Mk-1) and IsZero(Mkzero1) and IsZero(Mkzero2);

    MTb := Submatrix(M, 1, kappa+1, kappa, alpha-kappa);
    Mbzero := Submatrix(M, kappa+1, kappa+1, gamma-kappa, alpha-kappa);
    MSb := Submatrix(M, gamma+1, kappa+1, delta, alpha-kappa);
    boolb :=  IsZero(2*MTb) and IsZero(Mbzero) and IsZero(2*MSb);

    MT2 := Submatrix(M, 1, alpha+1, kappa, beta-delta-gamma+kappa);
    MT1 := Submatrix(M, kappa+1, alpha+1, gamma-kappa, beta-delta-gamma+kappa);
    MSq := Submatrix(M, gamma+1, alpha+1, delta, beta-delta-gamma+kappa);
    boolq :=  IsZero(2*MT2) and IsZero(2*MT1);

    Mgkzero := Submatrix(M,1,alpha+beta-delta-gamma+kappa+1,kappa,gamma-kappa);
    Mgk := Submatrix(M,kappa+1,alpha+beta-delta-gamma+kappa+1,
                    gamma-kappa,gamma-kappa);
    MR := Submatrix(M,gamma+1, alpha+beta-delta-gamma+kappa+1,
                    delta,gamma-kappa);
    boolgk :=  IsZero(Mgkzero) and IsOne(Mgk-1) and IsBinaryMatrix(MR);

    Mdzero1 := Submatrix(M, 1, alpha+beta-delta+1, kappa,delta);
    Mdzero2 := Submatrix(M, kappa+1, alpha+beta-delta+1, gamma-kappa,delta);
    Md := Submatrix(M, gamma+1, alpha+beta-delta+1, delta,delta);
    boold :=  IsZero(Mdzero1) and IsZero(Mdzero2) and IsOne(Md);

    return  boolk and boolb and boolq and boolgk and boold;
end intrinsic;

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
/*                                                          */
/************************************************************/
intrinsic IsStandardFormMatrix(M::ModMatRngElt, L::SeqEnum) -> BoolElt
{
Return true if and only if the matrix M over Z4, where the ones in the first alpha 
coordinates are represented by twos, is a generator matrix of a Z2Z4-additive 
code of type (alpha, beta; gamma, delta; kappa) in standard form.
}
    require #L eq NUMBER_OF_PARAMETERS: "Argument 2 must contain 5 parameters";
    require ElementType(L) eq RngIntElt: "Argument 2 must contain integer numbers"; 
    require (Type(BaseRing(M)) eq RngIntRes) and (BaseRing(M) eq Z4): 
            "Argument 1 is not a matrix over Z4";

    alpha := L[1];
    beta := L[2];
    gamma := L[3];
    delta := L[4];
    kappa := L[5];
    n := Ncols(M);
    m := Nrows(M);
    
    cond := (alpha ge 0) and (beta ge 0) and (gamma ge 0) and (delta ge 0) and
            (kappa ge 0) and (alpha+beta gt 0) and (alpha+beta eq n) and 
            (gamma+delta gt 0) and (gamma+delta eq m) and (beta ge delta) and
            (alpha+beta ge gamma+delta) and (kappa le Min(gamma,alpha)) and
            (beta+kappa ge delta+gamma);
    if (not cond) then 
        return false;
    end if;

    Mk := Submatrix(M, 1, 1, kappa, kappa); 
    
    Mkzero1 := Submatrix(M, kappa+1,1, gamma-kappa,kappa);
    Mkzero2 := Submatrix(M, gamma+1, 1, delta, kappa);
    boolk := IsOne(Mk-1) and IsZero(Mkzero1) and IsZero(Mkzero2);

    MTb := Submatrix(M, 1, kappa+1, kappa, alpha-kappa);
    Mbzero := Submatrix(M, kappa+1, kappa+1, gamma-kappa, alpha-kappa);
    MSb := Submatrix(M, gamma+1, kappa+1, delta, alpha-kappa);
    boolb :=  IsZero(2*MTb) and IsZero(Mbzero) and IsZero(2*MSb);

    MT2 := Submatrix(M, 1, alpha+1, kappa, beta-delta-gamma+kappa);
    MT1 := Submatrix(M, kappa+1, alpha+1, gamma-kappa, beta-delta-gamma+kappa);
    MSq := Submatrix(M, gamma+1, alpha+1, delta, beta-delta-gamma+kappa);
    boolq :=  IsZero(2*MT2) and IsZero(2*MT1);

    Mgkzero := Submatrix(M,1,alpha+beta-delta-gamma+kappa+1,kappa,gamma-kappa);
    Mgk := Submatrix(M,kappa+1,alpha+beta-delta-gamma+kappa+1,
                    gamma-kappa,gamma-kappa);
    MR := Submatrix(M,gamma+1, alpha+beta-delta-gamma+kappa+1,
                    delta,gamma-kappa);
    boolgk :=  IsZero(Mgkzero) and IsOne(Mgk-1) and IsBinaryMatrix(MR);

    Mdzero1 := Submatrix(M, 1, alpha+beta-delta+1, kappa,delta);
    Mdzero2 := Submatrix(M, kappa+1, alpha+beta-delta+1, gamma-kappa,delta);
    Md := Submatrix(M, gamma+1, alpha+beta-delta+1, delta,delta);
    boold :=  IsZero(Mdzero1) and IsZero(Mdzero2) and IsOne(Md);

    return  boolk and boolb and boolq and boolgk and boold;
end intrinsic;

