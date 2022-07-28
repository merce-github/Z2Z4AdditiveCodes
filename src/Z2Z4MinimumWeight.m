///////////////////////////////////////////////////////////////////////////////
/////////       Copyright 2015-2022 Cristina Diéguez-Martí, Cristina    ///////           
/////////        Fernández-Córdoba, Jaume Pujol, Adrián Torres-Martín   ///////
/////////        and Mercè Villanueva                                   ///////
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
/* Project name: Z2Z4-additive codes package                 */
/* File name: Z2Z4MinimumWeight.m                            */
/*                                                           */
/* Comments: Package developed within the CCSG group         */
/*                                                           */
/* Authors: Cristina Diéguez, Cristina Fernández-Córdoba,    */
/*          Jaume Pujol, Marta Pujol, Adrián Torres-Martín   */
/*          and Mercè Villanueva                             */
/*                                                           */
/* Revision version and last date: v1.2   19-07-2016         */
/*                                 v1.3   14-08-2016         */ 
/*                                 v2.0   12-02-2018         */
/*    User defined type            v2.1   30-01-2019         */
/*                                 v2.2   05-04-2019         */
/*    Added B-Z method             v2.3   28-04-2022         */
/*                                 v2.4   14-07-2022         */
/*                                                           */
/*************************************************************/
//Uncomment freeze when package finished
freeze;

intrinsic Z2Z4MinimumWeight_version() -> SeqEnum
{Return the current version of this package.}
    version := [2, 4];
    return version;
end intrinsic;

//needs Z2Z4AdditiveCode.m file
//needs Z2Z4StandardForm.m file

///////////////////////////////////////////////////////////////////////////////
declare verbose IgnoreWeightAttributes, 1;
///////////////////////////////////////////////////////////////////////////////

/******************************************************************
    GLOBAL VARIABLES
*******************************************************************/
Z4 := Integers(4);
MAX_CARDINAL_DISTRIBUTION := 32768; // Maximum code cardinal for brute force
MAX_DELTA_COSETS := 6; // Maximum parameter delta for the KernelCosets method
MAX_TORSION := 3; // Maximum value gamma-kappa for the Brouwer-Zimmermann method
MIN_GAUSSIAN_APPROX := 80; // Minimum length n to use the gaussian approximation
                           // instead of the binomial distribution
MAX_EXPECTED_WEIGHT := 2; // Maximimum value of the expected minimum weight 
                          // to use the Brouwer-Zimmermann method
///////////////////////////////////////////////////////////////////////////////
//////////// Functions we need from other files of the package ////////////////
///////////////////////////////////////////////////////////////////////////////

import "Z2Z4AdditiveCodes.m": zero_Code,
                              UpdateMinimumLeeWeight,
                              UpdateMinimumLeeWeightLowerBound,
                              UpdateMinimumLeeWeightUpperBound,
                              UpdateMinimumLeeWeightWord,
                              UpdateLeeWeightDistribution;
import "Z2Z4StandardForm.m" : DiagonalizeDeltaMatrix,
                              DiagonalizeKappaMatrix;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////                      THE MINIMUM WEIGHT                        ////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/****************************************************************/
/*                                                              */
/* Function name: KernelCosetRepresentativesZ2                  */
/* Parameters:  C                                               */
/* Function description: Given a Z2Z4-additive code, return the */
/*   kernel and coset representatives of the corresponding      */
/*   Gray map image of C.                                       */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/* Output parameters description:                               */
/*   - The kernel of C as a binary linear code                  */
/*   - Coset representatives as a sequence of codewords of C    */
/*                                                              */
/****************************************************************/ 
function KernelCosetRepresentativesZ2(C)
    codeZ4 := C`Code;
    kernelZ4 := KernelZ2CodeZ4(codeZ4);
    kernelZ2Z4 := Z2Z4AdditiveCode(kernelZ4, C`Alpha);
      _, kernelZ2 := HasLinearGrayMapImage(kernelZ2Z4);
  
    grayMap := GrayMap(C);
    V := VectorSpace(GF(2), BinaryLength(C));
    Q, f := quo<RSpace(codeZ4) | RSpace(kernelZ4)>;
    degreeQ := Degree(Q);

    leadersZ2Z4 := [];
    leadersZ2 := [];
    if degreeQ gt 0 then
        R := RSpace(Integers(2), degreeQ);
        leadersZ2Z4 := [(Q!x)@@f : x in R | not IsZero(x)];
        leadersZ2 := [V!grayMap(leader) : leader in leadersZ2Z4];
    end if;

    return kernelZ2, leadersZ2; 
end function;

/****************************************************************/
/*                                                              */
/* Function name: MinimumLeeWeightBruteForce                    */
/* Parameters:  C                                               */
/* Function description: Given a Z2Z4-additive code C, return   */
/*   the minimum Lee weight of the codewords belonging to the   */
/*   code C, which is also the minimum Lee distance between any */
/*   two codewords. The method used is by exhaustive search.    */                  
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/* Output parameters description:                               */
/*   - The minimum Lee weight of the code                       */
/*                                                              */
/****************************************************************/ 
function MinimumLeeWeightBruteForce(C) 
    alpha := C`Alpha;
    if not IsVerbose("IgnoreWeightAttributes") then
        assignedLowerBound := C`MinimumLeeWeightLowerBound;
        minimumWeight := C`MinimumLeeWeightUpperBound;
    else
        assignedLowerBound := 0;
        minimumWeight := BinaryLength(C);
    end if;

    for c in C`Code do
        if minimumWeight eq assignedLowerBound then
            return minimumWeight;
        end if;
        weightLeeCodeword := LeeWeight(c, alpha);
        if ((minimumWeight gt weightLeeCodeword) and 
                                           (not IsZero(weightLeeCodeword))) then 
            minimumWeight := weightLeeCodeword;
        end if; 
    end for; 

    return minimumWeight;
end function;

/****************************************************************/
/*                                                              */
/* Function name: MinimumLeeWeightKernelCosets                  */
/* Parameters:  C                                               */
/* Function description: Given a Z2Z4-additive code C, return   */
/*   the minimum Lee weight of the codewords belonging to the   */
/*   code C, which is also the minimum Lee distance between any */
/*   two codewords. The method used is by dividing the code     */
/*   into cosets of the kernel and applying the same function   */
/*   for linear codes.                                          */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/* Output parameters description:                               */
/*   - The minimum Lee weight of the code                       */
/*                                                              */
/****************************************************************/ 
function MinimumLeeWeightKernelCosets(C)  
    kernel, cosetRepresentatives := KernelCosetRepresentativesZ2(C);        

    minimumWeight := MinimumWeight(kernel);
    numCosets := #cosetRepresentatives; 

    if not IsVerbose("IgnoreWeightAttributes") then
        assignedLowerBound := C`MinimumLeeWeightLowerBound;
    else
        assignedLowerBound := 0;
    end if;

    if (numCosets gt 0) then
        kernelBasis := Basis(kernel);
        n := BinaryLength(C);     
        for i := 1 to numCosets do
            if minimumWeight eq assignedLowerBound then
                return minimumWeight;
            end if;
            coset := LinearCode<GF(2), n | kernelBasis cat [cosetRepresentatives[i]]>;
            minimumWeightCoset := MinimumWeight(coset);
            
            if (minimumWeight gt minimumWeightCoset) then
                minimumWeight := minimumWeightCoset;
            end if;
        end for;
    end if;

    return minimumWeight; 
end function;

/****************************************************************/
/*                                                              */
/* Function name: MinimumWordBruteForce                         */
/* Parameters:  C                                               */
/* Function description: Given a Z2Z4-additive code C, return   */
/*   one codeword of the code C having minimum Lee weight.      */
/*   The method used is by exhaustive search.                   */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/* Output parameters description:                               */
/*   - A codeword with minimum Lee weight                       */
/*                                                              */
/****************************************************************/ 
function MinimumWordBruteForce(C)   
    alpha := C`Alpha;
    minimumWord := C`Code.1;
    minimumWeight := LeeWeight(minimumWord, alpha); 

    for c in C`Code do
        weightLeeCodeword := LeeWeight(c, alpha);
        if ((minimumWeight gt weightLeeCodeword) and 
                                           (not IsZero(weightLeeCodeword))) then 
            minimumWeight := weightLeeCodeword;
            minimumWord := c;    
        end if; 
    end for; 

    return minimumWord;
end function;

/****************************************************************/
/*                                                              */
/* Function name: MinimumWordKernelCosets                       */
/* Parameters:  C                                               */
/* Function description: Given a Z2Z4-additive code C, return   */
/*   one codeword of the code C having minimum Lee weight.      */
/*   The method used is by dividing the code into cosets of the */
/*   kernel and applying the same function for linear codes.    */    
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/* Output parameters description:                               */
/*   - A codeword with minimum Lee weight                       */
/*                                                              */
/****************************************************************/ 
function MinimumWordKernelCosets(C)  
    kernel, cosetRepresentatives := KernelCosetRepresentativesZ2(C);
    
    minimumWord := MinimumWord(kernel);
    minimumWeight := Weight(minimumWord);
    numCosets := #cosetRepresentatives; 

    if (numCosets gt 0) then
        kernelBasis := Basis(kernel);        
        n := BinaryLength(C);
        for i := 1 to numCosets do
            coset := LinearCode<GF(2), n | kernelBasis cat [cosetRepresentatives[i]]>;
            minimumWordCoset := MinimumWord(coset);
            minimumWeightCoset := Weight(minimumWordCoset);
            if (minimumWeight gt minimumWeightCoset) then
                minimumWeight := minimumWeightCoset;
                minimumWord := minimumWordCoset;
            end if;
        end for;
    end if;
    
    return minimumWord@@GrayMap(C); 
end function;

/****************************************************************/
/*                                                              */
/* Function name: MinimumWordsBruteForce                        */
/* Parameters:  C                                               */
/* Function description: Given a Z2Z4-additive code C and an    */
/*   integer number numWords, return the set of all codewords   */
/*   of C having minimum Lee weight if numWords is negative.    */
/*   If numWords is a non-negative integer, then the function   */
/*   return at least that total of codewords of minimum weight. */
/*   The method used is by exhaustive search.                   */  
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/* Output parameters description:                               */
/*   - All codewords with minimum Lee weight                    */
/*                                                              */
/****************************************************************/ 
function MinimumWordsBruteForce(C, numWords)
    alpha := C`Alpha;     
    minimumWeight := BinaryLength(C); 
    minimumWords := [ ];
    contMinWords := 0;

    for c in C`Code do
        weightLeeCodeword := LeeWeight(c, alpha);
        if ((minimumWeight gt weightLeeCodeword) and 
                                           (not IsZero(weightLeeCodeword))) then 
            minimumWeight := weightLeeCodeword;
            minimumWords := [c];
            contMinWords := 1;
        else
            if (minimumWeight eq weightLeeCodeword) then
                Append(~minimumWords, c);
                contMinWords := contMinWords + 1;
            end if;  
        end if;
        if contMinWords eq numWords then 
            return Set(minimumWords);
        end if;
    end for; 

    return Set(minimumWords);     
end function;

/****************************************************************/
/*                                                              */
/* Function name: MinimumWordsKernelCosets                      */
/* Parameters:  C                                               */
/* Function description: Given a Z2Z4-additive code C and an    */
/*   integer number numWords, return the set of all codewords   */
/*   of C having minimum Lee weight if numWords is negative.    */
/*   If numWords is a non-negative integer, then the function   */
/*   return at least that total of codewords of minimum weight. */
/*   The method used is by dividing the code into cosets of the */
/*   kernel and applying the same function for linear codes.    */        
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/* Output parameters description:                               */
/*   - All codewords with minimum Lee weight                    */
/*                                                              */
/****************************************************************/ 
function MinimumWordsKernelCosets(C, numWords)
    kernel, cosetRepresentatives := KernelCosetRepresentativesZ2(C); 

    if numWords eq -1 then
        numWords := #C;
    end if;
    minimumWords := MinimumWords(kernel);
    if #minimumWords ge numWords then
        return minimumWords@@GrayMap(C);
    end if;
    minimumWeight := Weight(Setseq(minimumWords)[1]);
    numCosets := #cosetRepresentatives; 

    if (numCosets gt 0) then  
        kernelBasis := Basis(kernel);        
        n := BinaryLength(C);   
        for i := 1 to numCosets do
            coset := LinearCode<GF(2), n | kernelBasis cat [cosetRepresentatives[i]]>;
            minimumWeightCoset := MinimumWeight(coset);
            if ((minimumWeight gt minimumWeightCoset) and 
                                          (not IsZero(minimumWeightCoset))) then
                minimumWeight := minimumWeightCoset;
                minimumWords := MinimumWords(coset);
                contMinWords := #minimumWords;
            else
                if (minimumWeight eq (MinimumWeight(coset))) then
                    minimumWords join:= (MinimumWords(coset));
                    contMinWords := #minimumWords;
                end if;
            end if; 
            if contMinWords ge numWords then 
                return minimumWords@@GrayMap(C);
            end if;       
        end for;
    end if;
    
    return minimumWords@@GrayMap(C);
end function;

/****************************************************************/
/*                                                              */
/* Function name: InitializeComb                                */
/* Parameters: n, t                                             */
/* Function description: Initialize a t+1 sequence in           */
/*   increasing order. The t+1 position is the sentinel         */
/* Input parameters description:                                */
/*   - n : list length                                          */ 
/*   - t : combinations number                                  */
/*                                                              */
/* Function developed by J. Pujol                               */
/*                                                              */
/* Comments: Revolving-door combinations                        */
/*   Compute t-combinations of {0,1,...,n-1}, 1 < t < n, using  */
/*   an algorithm of W.H. Payne.                                */
/*   See D.E. Knuth, Combinatorial Algorithms, p. 363           */
/*                                                              */
/****************************************************************/
InitializeComb := function (n, t)
   L := [j-1 : j in [1..t]] cat [n];
   return L;
end function;

/****************************************************************/
/*                                                              */
/* Function name: NextCombPair                                  */
/* Parameters: L, n, t, R3                                      */
/* Function description: Compute the next combination and return*/
/*   a sequence of two numbers representing which number has    */ 
/*   left the combination and which one has entered.            */
/* Input parameters description:                                */
/*   - L : current combination                                  */
/*   - R3 : algorithm status                                    */
/*        : (false) first step                                  */
/*        : (true) next step                                    */
/*   - n : list length                                          */ 
/*   - t : combinations number                                  */
/*                                                              */
/* Function first developed by J. Pujol and updated by A. Torres*/
/*                                                              */
/* Comments: Revolving-door combinations                        */
/*   Compute t-combinations of {0,1,...,n-1}, 1 < t < n, using  */
/*   an algorithm of W.H. Payne.                                */
/*   See D.E. Knuth, Combinatorial Algorithms, p. 363           */
/*                                                              */
/****************************************************************/
function NextCombPair(L, n, t, R3)
    R4 := false; //R4 step
    R5 := false; //R5 step
    R6 := false; //recursion control
	if R3 then
	    if IsOdd(t) then
	        if L[1]+1 lt L[2] then
                pair := [L[1], L[1]+1];
                L[1] := L[1] + 1;			    
                R4 := false;
                R5 := false;
                R6 := false;
            else
                j := 2;
                R4 := true;
                R5 := false;
                R6 := true;
            end if;
        else
            if L[1] gt 0 then
                pair := [L[1], L[1]-1];
                L[1] := L[1]-1;				
                R4 := false;
                R5 := false;
                R6 := false;
            else
                j := 2;
                R4 := false;
                R5 := true;
                R6 := true;
            end if;
        end if;
	end if;
	while R6 do
	    if R4 then
            if L[j] ge j then
                pair := [L[j], j-2];
                L[j] := L[j-1];
                L[j-1] := j-2;
                R4 := false;			
                R5 := false;			        
                R6 := false;
            else
                j := j+1;
                R4 := false;
                R5 := true;			        
                R6 := false;
            end if;
	    end if;
	    if R5 then
            if L[j]+1 lt L[j+1] then
                pair := [L[j-1], L[j]+1];
                L[j-1] := L[j];
                L[j] := L[j]+1;	
                R4 := false;
                R5 := false;
                R6 := false;		
	        else
		        j := j+1;
                if j le t then	
                    R4 := true;	
                    R5 := true;	          
                    R6 := true;
                else
                    R4 := false;
                    R5 := false;
                    R6 := false;
	            end if;	          
            end if;
        end if;
	end while;
	return L, pair;
end function;

/****************************************************************/
/*                                                              */
/* Function name: SupportInverseGrayMap                         */
/* Parameters: supportZ2Vector, gamma                           */
/* Function description: Given a set of integers, representing  */
/*   the support of a binary vector, and an integer gamma,      */
/*   return a sequence of tuples representing the support of    */
/*   the inverse of the binary vector under the Gray map taking */
/*   the first gamma coordinates as the identity.               */
/* Input parameters description:                                */
/*   - supportZ2Vector : Subset of integers                     */
/*   - gamma : Integer                                          */
/* Output parameters description:                               */
/*   - A sequence of tuples representing the support of a vector*/
/*     over Z4. Tuple <i,x> means that in position i there is x.*/
/*                                                              */
/* Function first developed by M. Pujol and updated by A. Torres*/
/*                                                              */
/****************************************************************/
SupportInverseGrayMap := function(supportZ2Vector, gamma)  
    // old version 
    // supportZ4Vector := [<position, 1> : position in 
    //                                     (supportZ2Vector meet {1..gamma})];
    // supportZ2VectorLast := supportZ2Vector diff {1..gamma};
    supportZ4Vector := [];
    supportZ2VectorLast := {};
    for position in supportZ2Vector do
        if position in {1..gamma} then
            Append(~supportZ4Vector, <position, 1>);
        else
            Include(~supportZ2VectorLast, position);
        end if;
    end for;

    for position in supportZ2VectorLast do
        q, r := Quotrem(position-gamma, 2);
        if r eq 1 then
            // when the subset contains "11" corresponds to 2
            if position+1 in supportZ2VectorLast then 
                Append(~supportZ4Vector, <q+gamma+1, 2>);
            // when the subset contains "10" corresponds to 3 
            else 
                Append(~supportZ4Vector, <q+gamma+1, 3>);
            end if;
        else
            // when the subset contains "01" corresponds to 1
            if position-1 notin supportZ2VectorLast then 
                Append(~supportZ4Vector, <q+gamma, 1>);
            end if;
        end if;
    end for;
    
    return supportZ4Vector;
end function;

/****************************************************************/
/*                                                              */
/* Function name: BitsSwapInverseGrayMap                        */
/* Parameters: supportZ2Vector, gamma, pair                     */
/* Function description: Given a set of integers representing   */
/*   the support of a binary vector, an integer gamma, and a    */
/*   pair of integers representing a coordinate position x      */
/*   leaving the support and a coordinate position y entering it*/
/*   (from the last step in the Revolving Door algorithm),      */
/*   return a pair of tuples representing the support of two    */
/*   quaternary vectors (and their sign) that need to be added  */
/*   to obtain the inverse of the new binary vector (after      */
/*   changing the bits in x, y positions) under the Gray map    */
/*   taking the first gamma coordinates as the identity.        */
/* Input parameters description:                                */
/*   - supportZ2Vector : Subset of integers                     */
/*   - gamma : Integer                                          */
/*   - pair : Sequence of two integers                          */
/* Output parameters description:                               */
/*   - A pair of tuples, each one representing a vector over Z4 */
/*     of weight one. Tuple <i,x> represents the vector having  */
/*     x (+1/-1) in position i.                                 */  
/*     These two vectors also represent the rows of the gene-   */
/*     rator matrix in standard form that need to be added to   */
/*     obtain the next codeword from a codeword such it has     */
/*     the same weight in the information coordinates.          */
/*                                                              */
/* Function developed by A. Torres                              */
/*                                                              */
/****************************************************************/
BitsSwapInverseGrayMap := function(supportZ2Vector, gamma, pair)
    supportZ2VectorLast := supportZ2Vector diff {1..gamma};

    // if it is in the first gamma coordinates, substract 1
    if pair[1] le gamma then
        supportZ4PairVectors := [<pair[1], -1>];
    else
        q, r := Quotrem(pair[1]-gamma, 2);
        if r eq 1 then
            if pair[1]+1 in supportZ2VectorLast then 
                supportZ4PairVectors := [<q+gamma+1, -1>]; 
            else 
                supportZ4PairVectors := [<q+gamma+1, 1>];
            end if;
        else
            if pair[1]-1 in supportZ2VectorLast then 
                supportZ4PairVectors := [<q+gamma, 1>];
            else 
                supportZ4PairVectors := [<q+gamma, -1>];
            end if;
        end if;
        Exclude(~supportZ2VectorLast, pair[1]);
    end if;

    // if it is in the first gamma coordinates, add 1
    if pair[2] le gamma then
        Append(~supportZ4PairVectors, <pair[2], 1>);
    else
        q, r := Quotrem(pair[2]-gamma, 2);
        if r eq 1 then
            if pair[2]+1 in supportZ2VectorLast then 
                Append(~supportZ4PairVectors, <q+gamma+1, 1>);
            else 
                Append(~supportZ4PairVectors, <q+gamma+1, -1>);
            end if;
        else
            if pair[2]-1 in supportZ2VectorLast then 
                Append(~supportZ4PairVectors, <q+gamma, -1>);
            else 
                Append(~supportZ4PairVectors, <q+gamma, 1>);
            end if;
        end if;
    end if;
    
    // return both tuples, ordered by the first element in tuple
    if supportZ4PairVectors[1][1] le supportZ4PairVectors[2][1] then
        return supportZ4PairVectors;
    else
        return [supportZ4PairVectors[2], supportZ4PairVectors[1]];
    end if;
end function;

/****************************************************************/
/*                                                              */
/* Function name: MatricesStandardForm_Brouwer                  */
/* Parameters: C, type                                          */
/* Function description: Given a Z2Z4-additive code and its type*/
/*   as a sequence [alpha, beta, gamma, delta, kappa], return a */
/*   sequence of generator matrices of the same code with       */
/*   disjoint information sets for the kappa and delta columns. */
/*   It also returns a sequence of permutations needed to       */
/*   transform the generator matrices into standard form.       */
/*   Finally, it returns a sequence with the coordinates that   */
/*   are not in the free information set for each generator     */
/*   matrix.                                                    */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/*   - type : Sequence [alpha, beta, gamma, delta, kappa]       */
/* Output parameters description:                               */
/*   - A sequence of generator matrices in standard form with   */
/*     dijoint free information sets                            */
/*   - A sequence of permutations to transform the generator    */
/*     matrices in standard form                                */
/*   - A sequence of sets with the coordinate position outside  */
/*     the free information set                                 */
/*                                                              */
/* Function first developed by M. Pujol and updated by A. Torres*/
/*                                                              */
/****************************************************************/
MatricesStandardForm_Brouwer := function(C, type)
    _, _, G_sf, p := StandardForm(C);
    nRows := NumberOfRows(G_sf);
    nCols := NumberOfColumns(G_sf);

    alpha := type[1];
    kappa := type[5];
    delta := type[4];

    GList := [ G_sf ];
    pList := [ p ];
    infoCoord := {1..kappa} join {(nCols-delta+1)..(nCols)};
    iList := [{1..nCols} diff infoCoord];

    kappaCols := kappa;
    deltaCols := delta;

    subG := ColumnSubmatrixRange(G_sf, kappa + 1, nCols - delta);
    if Ncols(subG) eq 0 then
        return GList, pList, iList; 
    end if;
    Ik := ColumnSubmatrixRange(G_sf, 1, kappa);
    Id := ColumnSubmatrixRange(G_sf, nCols - delta + 1, nCols);

    subGalpha := alpha - kappa;
    subC := Z2Z4AdditiveCode(subG, subGalpha);
    subC_type := Z2Z4Type(subC);  

    while type[3..5] eq subC_type[3..5] do
        _, f, subG_sf, subp := StandardForm(subC);
        // Transform subp over Sym(Ncols(subG)) as a permutation over Sym(nCols)
        pSeq := [1..kappaCols] cat [i+kappaCols : i in Eltseq(subp)] cat
                 [(nCols-deltaCols+1)..nCols];
        p := Sym(nCols)!pSeq;
        // M gives the linear combinations from subG to subG_sf 
        M := Matrix(Solution(subG, Matrix([subG_sf[i]@@f : i in [1..nRows]])));
        G_sf := HorizontalJoin(< M * Ik, subG_sf, M * Id>); 

        Append(~GList, G_sf);
        Append(~pList, pList[#pList]*p); 
        infoCoord := {(kappaCols+1)..(kappaCols+kappa)} join 
                     {(nCols-deltaCols-delta+1)..(nCols-deltaCols)};
        Append(~iList, {1..nCols} diff infoCoord);
        
        kappaCols := kappaCols + kappa;
        deltaCols := deltaCols + delta;
        
        subG := ColumnSubmatrixRange(G_sf, kappaCols + 1, nCols - deltaCols);
        if Ncols(subG) eq 0 then
            return GList, pList, iList; 
        end if;
        Ik := ColumnSubmatrixRange(G_sf, 1, kappaCols);
        Id := ColumnSubmatrixRange(G_sf, nCols - deltaCols + 1, nCols);

        subGalpha := subGalpha - kappa;
        subC := Z2Z4AdditiveCode(subG, subGalpha);
        subC_type := Z2Z4Type(subC); 
    end while;

    return GList, pList, iList;
end function;

/****************************************************************/
/*                                                              */
/* Function name: MinWeightVectorCombinations                   */
/* Parameters: G, subG, seqDifferenceTwoRowsG, seqSumTwoRowsG,  */
/*             r, type, distanceLower                           */
/* Function description: Given a generator matrix G of a Z2Z4-  */
/*   additive code, a submatrix of G (explained below), two     */
/*   sequences with precomuted vectors to speed the calculations*/
/*   (explained below), an integer r, the type of the code as a */
/*   sequence [alpha, beta, gamma, delta, kappa], and the       */
/*   current lower bound for the minimum weight, compute all    */
/*   codewords from the information vectors of Hamming weight r */
/*   and return one of minimum Lee weight and the new lower	    */
/*   bound for the next step.                                   */          
/* Input parameters description:                                */
/*   - G : Generator matrix of a Z2Z4-additive code             */
/*   - subG : G restricted to the columns corresponding to the  */
/*     positions outside the free information positions         */
/*   - seqDifferenceTwoRowsG : Sequence of sequences containing */
/*     in the i-th sequence the difference between the i-th row */
/*     and the j-th row of subG                                 */
/*   - seqSumTwoRowsG : Sequence of sequences containing in the */
/*     i-th sequence the sum of the i-th row with the j-th row  */
/*     of subG                                                  */
/*   - r : Integer						                        */
/*   - type : Sequence [alpha, beta, gamma, delta, kappa]       */
/*   - distanceLower : Current lower bound of the minimum weight*/
/* Output parameters description:                               */
/*   - The minimum weight, say w, of all codewords generated    */
/*     from G and the information vectors of weight r           */
/*   - A codeword of weight w                                   */
/*                                                              */
/* Function first developed by M. Pujol and updated by A. Torres*/
/*                                                              */
/****************************************************************/
MinWeightVectorCombinations := function(G, subG, seqDifferenceTwoRowsG, seqSumTwoRowsG, 
                                        r, type, distanceLower)
    lenBinInfoVector := type[3] + 2 * type[4];
    steps := Binomial(lenBinInfoVector, r);
    positionsOrder2 := {(type[5]+1)..type[3]};
    // an upper bound is the number of columns of G
    minWeight := (type[1]+type[5]) + 2 * (type[2]+type[4]);
    
    for k in [1..steps] do
        if k eq 1 then
            L := InitializeComb(lenBinInfoVector, r);
            minL := L;
            suppZ2InfoVector := { j + 1 : j in L[1..r]};
            suppZ4InfoVector := SupportInverseGrayMap(suppZ2InfoVector, type[3]);
            codewordZ2Z4 := &+[suppZ4InfoVector[i][2]*subG[suppZ4InfoVector[i][1]] : 
                               i in [1..#suppZ4InfoVector]];
        else
            L_old := L;
            L, pair := NextCombPair(L, lenBinInfoVector, r, true);
            pair := [pair[1]+1, pair[2]+1];
            suppZ2InfoVector := { j + 1 : j in L_old[1..r]};
            suppZ4InfoVector := BitsSwapInverseGrayMap(suppZ2InfoVector, type[3], pair);
            i := suppZ4InfoVector[1][1];
            j := suppZ4InfoVector[2][1];
            if suppZ4InfoVector[1][2]*suppZ4InfoVector[2][2] eq 1 then
                if suppZ4InfoVector[1][2] eq 1 then
                    codewordZ2Z4 +:= seqSumTwoRowsG[i][j-i+1];
                else
                    codewordZ2Z4 -:= seqSumTwoRowsG[i][j-i+1];
                end if;
            else
                if suppZ4InfoVector[1][2] eq 1 then
                    codewordZ2Z4 +:= seqDifferenceTwoRowsG[i][j-i+1];
                else
                    codewordZ2Z4 -:= seqDifferenceTwoRowsG[i][j-i+1];
                end if;
            end if;
        end if;

        vWeight := LeeWeight(codewordZ2Z4, type[1]);
        rFree := r-#({j+1 : j in L[1..r]} meet positionsOrder2);
        fullWeight := vWeight + rFree;
        if (fullWeight lt minWeight) then
	        minL := L;
            minWeight := fullWeight;
            
            if minWeight le distanceLower then
                suppZ2InfoVector := { j + 1 : j in minL[1..r]};
                suppZ4InfoVector := SupportInverseGrayMap(suppZ2InfoVector, type[3]);
                minVector := &+[suppZ4InfoVector[i][2]*G[suppZ4InfoVector[i][1]] : 
                                i in [1..#suppZ4InfoVector]];
                return minWeight, minVector;
            end if;
        end if;
    end for;

    suppZ2InfoVector := { j + 1 : j in minL[1..r]};
    suppZ4InfoVector := SupportInverseGrayMap(suppZ2InfoVector, type[3]);
    minVector := &+[suppZ4InfoVector[i][2]*G[suppZ4InfoVector[i][1]] : 
                    i in [1..#suppZ4InfoVector]];
    return minWeight, minVector;
end function;

/****************************************************************/
/*                                                              */
/* Function name: DifferenceAndSumOfTwoRows                     */
/* Parameters: G                                                */
/* Function description: Given a matrix G, return a sequence of */
/*   sequences, where in the i-th sequence there is the         */
/*   difference between the i-th row of G with the j-th row,    */
/*   for j>=1. It also returns a sequence of sequences, where   */
/*   in the i-th sequence there is the sum of the i-th row of G */
/*   with the j-th row, for j>=1.                               */
/* Input parameters description:                                */
/*   - G : Matrix                                               */
/* Output parameters description:                               */
/*   - A sequence of sequences containing in the i-th sequence  */
/*     the difference between the i-th row and the j-th row     */
/*     of G                                                     */
/*   - A sequence of sequences containing in the i-th sequence  */
/*     the sum of the i-th row and the j-th row of G            */
/*                                                              */
/* Function developed by A. Torres                              */
/*                                                              */
/****************************************************************/
DifferenceAndSumOfTwoRows := function(G)
    nRows := Nrows(G);
    seqDifferenceTwoRowsG := [];
    seqSumTwoRowsG := [];
    for i in [1..nRows] do
        DifferenceTwoRowsGi := [];
        SumTwoRowsGi := [];
        for j in [i..nRows] do
            Append(~DifferenceTwoRowsGi, G[i]-G[j]);
            Append(~SumTwoRowsGi, G[i]+G[j]);
        end for;
        Append(~seqDifferenceTwoRowsG, DifferenceTwoRowsGi);
        Append(~seqSumTwoRowsG, SumTwoRowsGi);
    end for;

    return seqDifferenceTwoRowsG, seqSumTwoRowsG;
end function;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4MinimumLeeWeight_ForBinaryQuaternaryLinear*/
/* Parameters: C                                                */
/* Function description: Given a Z2Z4-additive code C, determine*/
/*   if C is either a binary code, a quaternary code or its Gray*/
/*   map image is a binary linear code. If C is any of these    */
/*   cases, return true, the minimum weight of C and a codeword */
/*   of minimum weight. If the code is not binary, quaternary or*/
/*   linear under the Gray map, then the function returns false,*/
/*   a minimum weight of 0 and the codeword 0.                  */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/* Output parameters description:                               */
/*   - A boolean indicating whether C is binary, quaternary or  */
/*     its image by the Gray map is linear.                     */
/*   - Integer with the minimum weight of C (if the boolean is  */
/*     true)                                                    */
/*   - A codeword of minimum weight (if the boolean is true)    */
/*                                                              */
/* Function first developed by M. Pujol and updated by A. Torres*/
/*                                                              */
/****************************************************************/
function Z2Z4MinimumLeeWeight_ForBinaryQuaternaryLinear(C)
    type := Z2Z4Type(C);

    // C is a binary linear code
    if type[2] eq 0 then
        LBC := LinearBinaryCode(C);
        grayMap := GrayMap(C);
        minWord := MinimumWord(LBC);
        return true, Weight(minWord), minWord@@grayMap;
    end if;

    // C is a quaternary linear code
    if type[1] eq 0 then
        LQC := LinearQuaternaryCode(C);
        minWeight := MinimumLeeDistance(LQC);
        minWords := WordsOfLeeWeight(LQC, minWeight : NumWords := 1);
        return true, minWeight, Rep(minWords);
    end if;

    // the Gray map image is a binary linear code
    isBinaryLinear, grayMapImageCode := HasLinearGrayMapImage(C);
    if isBinaryLinear then
        grayMap := GrayMap(C);
        minWord := MinimumWord(grayMapImageCode);
        return true, Weight(minWord), minWord@@grayMap;
    end if;

    return false, _, _;
end function;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4MinimumLeeWeight_Brouwer                  */
/* Parameters: C                                                */
/* Function description: Given a Z2Z4-additive code, determine  */
/*   the minimum weight of the words belonging to the code C    */
/*   using the Brouwer's algorithm, which is also the minimum   */
/*   distance between any two codewords. The function also      */
/*   returns a codeword of minimum weight.                      */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/* Output parameters description:                               */
/*   - Integer with the minimum weight of C                     */
/*   - A codeword of minimum weight                             */
/*                                                              */
/* Function first developed by M. Pujol and updated by A. Torres*/
/*                                                              */
/****************************************************************/
function Z2Z4MinimumLeeWeight_Brouwer(C)

    type := Z2Z4Type(C);

    GList, pList, iList := MatricesStandardForm_Brouwer(C, type);
    subGList := [];
    difList := [];
    sumList := [];
    h := #GList;
    for i := 1 to h do
        subG := ColumnSubmatrix(GList[i], Sort(Setseq(iList[i])));
        Append(~subGList, subG);
        seqDifferenceTwoRowsG, seqSumTwoRowsG := DifferenceAndSumOfTwoRows(subG);
        Append(~difList, seqDifferenceTwoRowsG);
        Append(~sumList, seqSumTwoRowsG);
    end for;

    newType := type;
    newType[1] -:= newType[5];
    newType[2] -:= newType[4];

    // inizialice the variables distanceUpper and word 
    if type[1] ge type[3] then
        distanceUpper := type[1] + 2*type[2] - type[3] - 2*type[4] + 1;
    else
        distanceUpper := 2*type[1] + 2*type[2] - 2*type[3] - 2*type[4] + 2;
    end if;
    word := C!0;

    lenBinInfoVector := type[3] + 2 * type[4];
    for r := 1 to lenBinInfoVector do
        distanceLower := Maximum(0, h*(r-(type[3]-type[5])));
        for i := 1 to h do
            minWeight, minVector := MinWeightVectorCombinations(GList[i], subGList[i], 
                                    difList[i], sumList[i], r, newType, distanceLower);
            if minWeight lt distanceUpper then
                distanceUpper := minWeight;
                word := minVector^(pList[i]^(-1));
            end if;
            if r+1-(type[3]-type[5]) gt 0 then
                distanceLower +:= 1;
            end if;
        end for;

        if distanceLower ge distanceUpper then
            return distanceUpper, word;
        end if;
    end for;

    return distanceUpper, word;
end function;

/****************************************************************/
/*                                                              */
/* Function name: RearrangeRows                                 */
/* Parameters: M, type                                          */
/* Function description: Given a matrix M and a sequence        */
/*   type=[alpha,beta,gamma,delta,kappa], return a matrix       */
/*   rearranging the m rows in M the following order:           */
/*   [delta+1,...,m,1,...,delta]. The function also returns a   */
/*   matrix representing the rows operations.                   */
/* Input parameters description:                                */
/*   - M : Matrix over Z4                                       */
/*   - type : Sequence containing the type of a code            */
/* Output parameters description:                               */
/*   - Matrix M after rearranging the rows                      */
/*   - A matrix representing the row transformation             */ 
/*                                                              */
/* Function developed by A. Torres                              */
/*                                                              */
/****************************************************************/
function RearrangeRows(M, type)
    nRows := Nrows(M);
    I := IdentityMatrix(Z4, nRows);
    listRows := [M[i] : i in [(type[4]+1)..nRows]] cat
                [M[i] : i in [1..type[4]]];
    listRowsI := [I[i] : i in [(type[4]+1)..nRows]] cat
                 [I[i] : i in [1..type[4]]];
    E := Matrix(listRows);
    Q := Matrix(listRowsI);
    return E, Q;
end function;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4StandardFormDeficit                       */
/* Parameters: subG, alpha, kappa_Order2, type                  */
/* Function description: Given a matrix subG with alpha         */
/*   binary coordinates, an integer kappa_Order2 and a          */
/*   sequence type, return a matrix with an identity submatrix  */
/*   in the first kappa_Order2 rows and columns (and 0s in the  */
/*   rest of each column) and an identity submatrix in the last */
/*   delta (where delta is the number of order 4 rows in subG)  */
/*   rows and columns (and 0s in the rest of each column), after*/
/*   performing elementary rows operations and column           */
/*   permutations. The function also returns the full column    */
/*   permutation and a matrix representing the row operations.  */
/* Input parameters description:                                */
/*   - subG : Matrix over Z4                                    */
/*   - alpha : Number of quaternary coordinates of subG in subG */
/*   - kappa_Order2 : Dimension of the code generated by the    */
/*     quaternary coordinates of the order 2 rows of subG       */
/*   - type : Sequence containing the type of the original code */
/* Output parameters description:                               */
/*   - Matrix subG having the identity in the first kapp_Order2 */
/*     columns and in the last delta columns                    */
/*   - Permutation that represents the column operations        */
/*   - Matrix that represents the row operations                */ 
/*                                                              */
/* Function developed by A. Torres                              */
/*                                                              */
/****************************************************************/
function Z2Z4StandardFormDeficit(subG, alpha, kappa_Order2, type)
    subC := Z2Z4AdditiveCode(subG, alpha);
    subC_type := Z2Z4Type(subC);
    if type[4] gt 0 then
        subG4 := RowSubmatrixRange(subG, Nrows(subG)-type[4]+1, Nrows(subG));
        subC4 := Z2Z4AdditiveCode(subG4, alpha);
        _, _, _, p4 := StandardFormInfo(subC4`Code);
        subG4_sf, Q1 := EchelonForm(subG4^p4);
        subG4_rr, Q2 := RearrangeRows(subG4_sf^(p4^(-1)), Z2Z4Type(subC4));
        Q12 := DiagonalJoin(IdentityMatrix(Z4, Nrows(subG) - Nrows(subG4)), Q2*Q1);
        subG := VerticalJoin(RowSubmatrix(subG, Nrows(subG)-type[4]), subG4_rr);
    else
        Q12 := IdentityMatrix(Z4, Nrows(subG));
    end if;
    Ed, pd, Q3 := DiagonalizeDeltaMatrix(subG, subC_type[1], subC_type[4]);
    Ek, pk, Q4 := DiagonalizeKappaMatrix(Ed, type[3], kappa_Order2);
    P := pd*pk;
    Q := Q4*Q3*Q12;
    return Ek, P, Q;
end function;

/****************************************************************/
/*                                                              */
/* Function name: DiagonalizeDeltaDeficit                       */
/* Parameters: M, delta, def                                    */
/* Function description: Given a matrix M of size m x n, having */
/*   the identity matrix of size def in the last columns and    */
/*   between the rows m-delta-def+1 and m-delta, return M with  */
/*   zeros in the rest of positions of these columns, after     */
/*   performing elementary row operations. The function also    */
/*   returns a matrix representing the row operations.          */
/* Input parameters description:                                */
/*   - M : Matrix over Z4                                       */
/*   - delta : Integer with type of the last restricted matrix  */
/*   - def : Difference between delta, and the original delta   */
/* Output parameters description:                               */
/*   - Matrix M diagonalized in the last def columns            */
/*   - Matrix that represents the row operations                */
/*                                                              */
/* Function developed by A. Torres                              */
/*                                                              */
/****************************************************************/
function DiagonalizeDeltaDeficit(M, delta, def)
    ncols := Ncols(M);
    mrows := Nrows(M);

    mc := mrows-delta;
    
    // matrix to store the row linear combinations performed
    Q := IdentityMatrix(Z4, mrows);

    for j := ncols to ncols-def+1 by -1 do
        for i := mrows to mrows-delta+1 by -1 do
            if M[i][j] ne 0 then
                M[i] := M[mc][j]*M[i]-M[i][j]*M[mc];
                Q[i] := M[mc][j]*Q[i]-M[i][j]*Q[mc];
            end if;
        end for;
        for i := 1 to mrows-delta-def-1 do
            if M[i][j] ne 0 then
                M[i] := M[mc][j]*M[i]-M[i][j]*M[mc];
                Q[i] := M[mc][j]*Q[i]-M[i][j]*Q[mc];
            end if;
        end for;
        mc -:= 1;
    end for;
    
    return M, Q;
end function;

/****************************************************************/
/*                                                              */
/* Function name: DiagonalizeKappaDeficit                       */
/* Parameters: M, kappa, def                                    */
/* Function description: Given a matrix M of size m x n, having */
/*   the identity matrix of size def in the first columns and   */
/*   between the rows kappa+1 and kappa+def, return M with      */
/*   zeros in the rest of positions of these columns, after     */
/*   performing elementary row operations. The function also    */
/*   returns a matrix representing the row operations.          */
/* Input parameters description:                                */
/*   - M : Matrix over Z4                                       */
/*   - kappa : Integer with type of the last restricted matrix. */
/*   - def : Difference between kappa, and the original kappa.  */
/* Output parameters description:                               */
/*   - Matrix M diagonalized in the first def columns           */
/*   - Matrix that represents the row operations                */
/*                                                              */
/* Function developed by A. Torres                              */
/*                                                              */
/****************************************************************/
function DiagonalizeKappaDeficit(M, kappa, def)
    ncols := Ncols(M);
    mrows := Nrows(M);

    mc := kappa+1;

    // matrix to store the row linear combinations performed
    Q := IdentityMatrix(Z4, mrows);

    for j := 1 to def do
        for i := 1 to kappa do
            if M[i][j] ne 0 then
                M[i] := M[i]-M[mc];
                Q[i] := Q[i]-Q[mc];
            end if;
        end for;
        for i := kappa+def+1 to mrows do
            if M[i][j] ne 0 then
                M[i] := M[i]-M[mc];
                Q[i] := Q[i]-Q[mc];
            end if;
        end for;
        mc +:= 1;
    end for;
    
    return M, Q;
end function;

/****************************************************************/
/*                                                              */
/* Function name: MatricesStandardForm_Zimmermann               */
/* Parameters: C, type                                          */
/* Function description: Given a Z2Z4-additive code and its type*/
/*   as a sequence [alpha, beta, gamma, delta, kappa], return a */
/*   sequence of generator matrices of the same code with       */
/*   non-disjoint information sets for the kappa and delta      */
/*   columns. It also returns a sequence with the rank deficits */
/*   of each generator matrix and a sequence of permutations    */
/*   to transform the generator matrices into standard form.    */
/*   Finally, it returns a sequence with the coordinates that   */
/*   are not in the free information set for each generator     */
/*   matrix.                                                    */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/*   - type : Sequence [alpha, beta, gamma, delta, kappa]       */
/* Output parameters description:                               */
/*   - A sequence of generator matrices in standard form without*/
/*     dijoint free information sets                            */
/*   - A sequence of permutations to transform the generator    */
/*     matrices in standard form                                */
/*   - A sequence of sets with the coordinate position outside  */
/*     the free information set                                 */
/*                                                              */
/* Function developed by A. Torres                              */
/*                                                              */
/****************************************************************/
MatricesStandardForm_Zimmermann := function(C, type)
    _, _, G_sf, p := StandardForm(C);
    nRows := NumberOfRows(G_sf);
    nCols := NumberOfColumns(G_sf);

    alpha := type[1];    
    kappa := type[5];
    delta := type[4];
    gamma := type[3];

    GList := [ G_sf ];
    kList := [gamma-kappa]; 
    pList := [p];
    infoCoord := {1..kappa} join {(nCols-delta+1)..nCols};
    iList := [{1..nCols} diff infoCoord];

    kappaCols := kappa;
    deltaCols := delta;
    
    subG := ColumnSubmatrixRange(G_sf, kappa + 1, nCols - delta);
    if Ncols(subG) eq 0 then
        return GList, kList, pList, iList; 
    end if;
    Ik := ColumnSubmatrix(G_sf, kappa);
    Id := ColumnSubmatrixRange(G_sf, nCols - delta + 1, nCols);
    
    subGalpha := alpha - kappa;
    subC := Z2Z4AdditiveCode(subG, subGalpha);
    subC_type := Z2Z4Type(subC);

    // The column submatrix subG may generate a code with more rows of order 2
    // than the original code (some of the order 4 rows transform into order 2 rows).
    // However, we cannot use any of their coordinate positions as binary information
    // coordinates (kappa) once we append the rest of the matrix. Instead, we create 
    // another code, only with the first gamma rows in order to obtain the real
    // number of order 2 generators.
    subG_order2 := RowSubmatrix(subG, gamma);
    subC_order2 := Z2Z4AdditiveCode(subG_order2, subGalpha);
    subC_order2_type := Z2Z4Type(subC_order2);

    while subC_type[4] + subC_order2_type[5] gt 0 do
        defG, subP, M := Z2Z4StandardFormDeficit(subG, subC_type[1], subC_order2_type[5], type);
        pSeq := [1..kappaCols] cat [i+kappaCols : i in Eltseq(subP)] cat
                 [(nCols-deltaCols+1)..nCols];
        p := Sym(nCols)!pSeq;
        G := HorizontalJoin(<M*Ik, defG, M*Id>);
        midRowsG := RowSubmatrixRange(G, subC_order2_type[5]+1, nRows-subC_type[4]);
        if Nrows(midRowsG) eq 0 then
            H := G;
            infoDef := {};
        else
            midRowsC := Z2Z4AdditiveCode(midRowsG, alpha);
            _, _, midRowsG_sf, p2 := StandardForm(midRowsC);

            H := VerticalJoin(<RowSubmatrix(G, subC_order2_type[5])^p2, midRowsG_sf,
                               RowSubmatrixRange(G, nRows-subC_type[4]+1, nRows)^p2>);
            H := DiagonalizeKappaDeficit(H, subC_order2_type[5], kappa-subC_order2_type[5]);
            H := DiagonalizeDeltaDeficit(H, subC_type[4], delta-subC_type[4])^(p2^(-1));

            infoDef := {1..(kappa-subC_order2_type[5])} join 
                       {(nCols-delta+subC_type[4]+1)..nCols};
            infoDef := infoDef^(p2^(-1));
        end if;
        kdef := gamma-kappa + kappa-subC_order2_type[5] + 2*(delta - subC_type[4]);

        Append(~GList, H);
        Append(~kList, kdef);
        Append(~pList, pList[#pList]*p);
        infoCoord := {(kappaCols+1)..(kappaCols+subC_order2_type[5])} join 
                     {(nCols-deltaCols-subC_type[4]+1)..(nCols-deltaCols)};
        Append(~iList, {1..nCols} diff (infoCoord join infoDef));

        kappaCols +:= subC_order2_type[5];
        deltaCols +:= subC_type[4];
        subGalpha -:= subC_order2_type[5]; 

        subG := ColumnSubmatrixRange(H, kappaCols + 1, nCols - deltaCols);
        if Ncols(subG) eq 0 then
            return GList, kList, pList, iList;
        end if;
        Ik := ColumnSubmatrix(H, kappaCols);
        Id := ColumnSubmatrixRange(H, nCols - deltaCols + 1, nCols);

        subC := Z2Z4AdditiveCode(subG, subGalpha);
        subC_type := Z2Z4Type(subC);
        subG_order2 := RowSubmatrix(subG, gamma);
        subC_order2 := Z2Z4AdditiveCode(subG_order2, subGalpha);
        subC_order2_type := Z2Z4Type(subC_order2);
    end while;

    return GList, kList, pList, iList;
end function;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4MinimumLeeWeight_Zimmermann               */
/* Parameters: C                                                */
/* Function description: Given a Z2Z4-additive code, determine  */
/*   the minimum weight of the words belonging to the code C    */
/*   using the Brouwer-Zimmermann algorithm, which is also the  */
/*   minimum distance between any two codewords. The function   */
/*   also returns a codeword of minimum weight.                 */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/* Output parameters description:                               */
/*   - Integer with the minimum weight of C                     */
/*   - A codeword of minimum weight                             */
/*                                                              */
/* Function developed by A. Torres                              */
/*                                                              */
/****************************************************************/
function Z2Z4MinimumLeeWeight_Zimmermann(C)

    type := Z2Z4Type(C);

    GList, kList, pList, iList := MatricesStandardForm_Zimmermann(C, type);
    subGList := [];
    difList := [];
    sumList := [];
    h := #GList;
    for i := 1 to h do
        subG := ColumnSubmatrix(GList[i], Sort(Setseq(iList[i])));
        Append(~subGList, subG);
        seqDifferenceTwoRowsG, seqSumTwoRowsG := DifferenceAndSumOfTwoRows(subG);
        Append(~difList, seqDifferenceTwoRowsG);
        Append(~sumList, seqSumTwoRowsG);
    end for;

    newType := type;
    newType[1] -:= newType[5];
    newType[2] -:= newType[4];
    
    // inizialice the variables distanceUpper and word
    if type[1] ge type[3] then
        distanceUpper := type[1] + 2*type[2] - type[3] - 2*type[4] + 1;
    else
        distanceUpper := 2*type[1] + 2*type[2] - 2*type[3] - 2*type[4] + 2;
    end if;
    word := C!0;

    lenBinInfoVector := type[3] + 2 * type[4];
    for r := 1 to lenBinInfoVector do
        lb := 0;
        w0 := 1;
        while (lb lt distanceUpper) and (w0 lt lenBinInfoVector) do
            w0 +:= 1;
            lb := &+[Maximum(0, w0+1-x) : x in kList];
        end while;

        h := Max([i : i in [1..h] | kList[i] lt w0+1]);
        
        distanceLower := &+[Maximum(0, r-kList[i]) : i in [1..h]];
        for i := 1 to h do
            minWeight, minVector := MinWeightVectorCombinations(GList[i], subGList[i], 
                                    difList[i], sumList[i], r, newType, distanceLower);
            if minWeight lt distanceUpper then
                distanceUpper := minWeight;
                word := minVector^(pList[i]^(-1));
            end if;
            if r+1-kList[i] gt 0 then
                distanceLower +:= 1;
            end if;
        end for;

        if distanceLower ge distanceUpper then
                return distanceUpper, word;
        end if;
    end for;

    return distanceUpper, word;
end function;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4MinimumLeeWeight_Quaternary               */
/* Parameters: C                                                */
/* Function description: Given a Z2Z4-additive code, determine  */
/*   the minimum weight of the words belonging to the code C    */
/*   using the quaternary method, which is also the minimum     */
/*   distance between any two codewords.                        */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/* Output parameters description:                               */
/*   - Integer with the minimum weight of C                     */
/*                                                              */
/* Function developed by A. Torres                              */
/*                                                              */
/* Comments: If a Z2Z4-additive code C has a generator matrix G */
/*   of the form G=(B|Q), where B consists of the first alpha   */
/*   columns of G and Q of the last beta columns. Then the      */
/*   quaternary method computes the minimum distance of the     */
/*   quaternary code generated by (2B|Q|Q), which is twice the  */
/*   minimum distance of C.                                     */
/*                                                              */
/****************************************************************/
function Z2Z4MinimumLeeWeight_Quaternary(C)
    type := Z2Z4Type(C);
    G := GeneratorMatrix(C);
    binaryG := ColumnSubmatrix(G, type[1]);
    quaternaryG := ColumnSubmatrixRange(G, type[1]+1, type[1]+type[2]);
    doubleG := HorizontalJoin(<binaryG, quaternaryG, quaternaryG>);
    LQC := LinearCode(doubleG);
    return MinimumLeeDistance(LQC)/2;
end function;

/****************************************************************/
/*                                                              */
/* Function name: MinWeightVectorCombinationsWords              */
/* Parameters: G, seqDifferenceTwoRowsG, seqSumTwoRowsG, r,     */
/*             type, minWords, minWeight, numWords              */
/* Function description: Given a generator matrix G of a Z2Z4-  */
/*   additive code, two sequences containing the difference and */
/*   sum between each pair of rows of G (seqDifferenceTwoRowsG  */
/*   and seqSumTwoRowsG, respectively), an integer r, and the   */
/*   type of the code as a sequence [alpha, beta, gamma, delta, */
/*   kappa], compute all codewords from the information vectors */
/*   of Hamming weight r and return the minimum weight. It also */
/*   returns a set of codewords with minimum weight, expanding  */
/*   the one given in minWords.		                            */
/* Input parameters description:                                */
/*   - G : Generator matrix of a Z2Z4-additive code             */
/*   - seqDifferenceTwoRowsG : Sequence of sequences containing */
/*     in the i-th sequence the difference between the i-th row */
/*     and the j-th row of G                                    */
/*   - seqSumTwoRowsG : Sequence of sequences containing in the */
/*     i-th sequence the sum of the i-th row with the j-th row  */
/*     of G                                                     */
/*   - r : Integer						                        */
/*   - type : Sequence [alpha, beta, gamma, delta, kappa]       */
/*   - minWords : Set with current minimum weight codewords     */
/*   - minWeight : Current minimum weight                       */
/*   - numWords : Maximum number of minimum weight codewords    */
/* Output parameters description:                               */
/*   - The minimum weight, say w, of all codewords generated    */
/*     from G and the information vectors of weight r           */
/*   - All codewords of weight w                                */
/*                                                              */
/* Function developed by A. Torres                              */
/*                                                              */
/****************************************************************/
MinWeightVectorCombinationsWords := function(G, seqDifferenceTwoRowsG, seqSumTwoRowsG, 
                                             r, type, minWords, minWeight, numWords)
    lenBinInfoVector := type[3] + 2 * type[4];
    steps := Binomial(lenBinInfoVector, r);

    for k in [1..steps] do
        if k eq 1 then
            L := InitializeComb(lenBinInfoVector, r);
            minL := L;
            suppZ2InfoVector := { j + 1 : j in L[1..r]};
            suppZ4InfoVector := SupportInverseGrayMap(suppZ2InfoVector, type[3]);
            codewordZ2Z4 := &+[suppZ4InfoVector[i][2]*G[suppZ4InfoVector[i][1]] : 
                               i in [1..#suppZ4InfoVector]];
        else
            L_old := L;
            L, pair := NextCombPair(L, lenBinInfoVector, r, true);
            pair := [pair[1]+1, pair[2]+1];
            suppZ2InfoVector := { j + 1 : j in L_old[1..r]};
            suppZ4InfoVector := BitsSwapInverseGrayMap(suppZ2InfoVector, type[3], pair);
            i := suppZ4InfoVector[1][1];
            j := suppZ4InfoVector[2][1];
            if suppZ4InfoVector[1][2]*suppZ4InfoVector[2][2] eq 1 then
                if suppZ4InfoVector[1][2] eq 1 then
                    codewordZ2Z4 +:= seqSumTwoRowsG[i][j-i+1];
                else
                    codewordZ2Z4 -:= seqSumTwoRowsG[i][j-i+1];
                end if;
            else
                if suppZ4InfoVector[1][2] eq 1 then
                    codewordZ2Z4 +:= seqDifferenceTwoRowsG[i][j-i+1];
                else
                    codewordZ2Z4 -:= seqDifferenceTwoRowsG[i][j-i+1];
                end if;
            end if;
        end if;

        vWeight := LeeWeight(codewordZ2Z4, type[1]);
        if (vWeight lt minWeight) then
            minWords := {codewordZ2Z4};
            minWeight := vWeight;
        elif (vWeight eq minWeight) then
            Include(~minWords, codewordZ2Z4);
        end if;

        if #minWords eq numWords then
            return minWeight, minWords;
        end if;
    end for;

    return minWeight, minWords;
end function;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4MinimumLeeWords_ForBinaryQuaternaryLinear */
/* Parameters: C, numWords                                      */
/* Function description: Given a Z2Z4-additive code C, determine*/
/*   if C is either a binary code, a quaternary code or its Gray*/
/*   map image is a binary linear code. If C is any of these    */
/*   cases, return true and a set of numWords codewords of      */
/*   minimum weight. If numWords is not a positive integer, it  */
/*   returns all codewords of minimum weight. If the code is not*/
/*   binary, quaternary or linear under the Gray map, then the  */
/*   function returns false and an empty set.                   */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/*   - numWords : Integer with the maximum number of returned   */
/*     codewords                                                */
/* Output parameters description:                               */
/*   - A boolean indicating whether C is binary, quaternary or  */
/*     its image by the Gray map is linear.                     */
/*   - A set of numWords codewords of minimum weight (if the    */
/*     boolean is true)                                         */
/*                                                              */
/* Function developed by A. Torres                              */
/*                                                              */
/****************************************************************/
function Z2Z4MinimumLeeWords_ForBinaryQuaternaryLinear(C, numWords)
    type := Z2Z4Type(C);

    // C is a binary linear code
    if type[2] eq 0 then
        LBC := LinearBinaryCode(C);
        grayMap := GrayMap(C);
        if numWords lt 0 then
            return true, {w@@grayMap : w in MinimumWords(LBC)};   
        else
            return true, {w@@grayMap : w in MinimumWords(LBC : NumWords := numWords)};
        end if;
    end if;

    // C is a quaternary linear code
    if type[1] eq 0 then
        LQC := LinearQuaternaryCode(C);
        minWeight := MinimumLeeDistance(LQC);
        if numWords lt 0 then
            return true, WordsOfLeeWeight(LQC, minWeight);
        else
            return true, WordsOfLeeWeight(LQC, minWeight : NumWords := numWords);
        end if;
    end if;

    // the Gray map image is a binary linear code
    isBinaryLinear, grayMapImageCode := HasLinearGrayMapImage(C);
    if isBinaryLinear then
        grayMap := GrayMap(C);
        return true, {w@@grayMap : w in MinimumWords(grayMapImageCode)};
    end if;

    return false, _;
end function;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4MinimumLeeWords_Brouwer                   */
/* Parameters: C, numWords                                      */
/* Function description: Given a Z2Z4-additive code, returns a  */
/*   set of numWords codewords of minimum weight using the      */
/*   Brouwer's algorithm. If numWords is not a positive integer,*/
/*   it returns all codewords of minimum weight.                */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/*   - numWords : Integer with the maximum number of returned   */
/*     codewords                                                */
/* Output parameters description:                               */
/*   - A set of numWords codewords of minimum weight            */
/*                                                              */
/* Function developed by A. Torres                              */
/*                                                              */
/****************************************************************/
function Z2Z4MinimumLeeWords_Brouwer(C, numWords)

    type := Z2Z4Type(C);

    GList, pList, iList := MatricesStandardForm_Brouwer(C, type);
    difList := [];
    sumList := [];
    h := #GList;
    for i := 1 to h do
        seqDifferenceTwoRowsG, seqSumTwoRowsG := DifferenceAndSumOfTwoRows(GList[i]);
        Append(~difList, seqDifferenceTwoRowsG);
        Append(~sumList, seqSumTwoRowsG);
    end for;

    // inizialice the variables distanceUpper and words 
    if type[1] ge type[3] then
        distanceUpper := type[1] + 2*type[2] - type[3] - 2*type[4] + 1;
    else
        distanceUpper := 2*type[1] + 2*type[2] - 2*type[3] - 2*type[4] + 2;
    end if;
    words := {};

    lenBinInfoVector := type[3] + 2 * type[4];
    for r := 1 to lenBinInfoVector do
        for i := 1 to h do
            distanceUpper, words := MinWeightVectorCombinationsWords(GList[i]^(pList[i]^(-1)), 
                                    difList[i]^(pList[i]^(-1)), sumList[i]^(pList[i]^(-1)), 
                                    r, type, words, distanceUpper, numWords);
        end for;
        distanceLower := Maximum(0, h*(r+1-(type[3]-type[5])));
        
        if (distanceLower gt distanceUpper) then
            if (numWords gt 0) and (#words ge numWords) then 
                return words;
            elif numWords le 0 then
                return words;
            end if;
        end if;
    end for;

    return words;
end function;

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4MinimumLeeWords_Zimmermann                */
/* Parameters: C, numWords                                      */
/* Function description: Given a Z2Z4-additive code, returns a  */
/*   set of numWords codewords of minimum weight using the      */
/*   Brouwer-Zimmermann algorithm. If numWords is not a positive*/
/*   integer, it returns all codewords of minimum weight.       */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/*   - numWords : Integer with the maximum number of returned   */
/*     codewords                                                */
/* Output parameters description:                               */
/*   - A set of numWords codewords of minimum weight            */
/*                                                              */
/* Function developed by A. Torres                              */
/*                                                              */
/****************************************************************/
function Z2Z4MinimumLeeWords_Zimmermann(C, numWords)

    type := Z2Z4Type(C);
    
    GList, kList, pList, iList := MatricesStandardForm_Zimmermann(C, type);
    difList := [];
    sumList := [];
    h := #GList; 
    for i := 1 to h do
        seqDifferenceTwoRowsG, seqSumTwoRowsG := DifferenceAndSumOfTwoRows(GList[i]);
        Append(~difList, seqDifferenceTwoRowsG);
        Append(~sumList, seqSumTwoRowsG);
    end for;

    // inizialice the variables distanceUpper and word
    if type[1] ge type[3] then
        distanceUpper := type[1] + 2*type[2] - type[3] - 2*type[4] + 1;
    else
        distanceUpper := 2*type[1] + 2*type[2] - 2*type[3] - 2*type[4] + 2;
    end if;
    words := {};

    lenBinInfoVector := type[3] + 2 * type[4];
    for r := 1 to lenBinInfoVector do
        for i := 1 to h do
            distanceUpper, words := MinWeightVectorCombinationsWords(GList[i]^(pList[i]^(-1)), 
                                    difList[i]^(pList[i]^(-1)), sumList[i]^(pList[i]^(-1)), 
                                    r, type, words, distanceUpper, numWords);
        end for;
        distanceLower := &+[Maximum(0, r+1-kList[i]) : i in [1..h]];

        if (distanceLower gt distanceUpper) then
            if (numWords gt 0) and (#words ge numWords) then 
                return words;
            elif numWords le 0 then
                return words;
            end if;
        end if;
    end for;

    return words;
end function; 

/****************************************************************/
/*                                                              */
/* Function name: Z2Z4MinimumLeeWords_Quaternary                */
/* Parameters: C, numWords                                      */
/* Function description: Given a Z2Z4-additive code, returns a  */
/*   set of numWords codewords of minimum weight using the      */
/*   quaternary method. If numWords is not a positive integer,  */
/*   it returns all codewords of minimum weight.                */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/*   - numWords : Integer with the maximum number of returned   */
/*     codewords                                                */
/* Output parameters description:                               */
/*   - A set of numWords codewords of minimum weight            */
/*                                                              */
/* Function developed by A. Torres                              */
/*                                                              */
/* Comments: If a Z2Z4-additive code C has a generator matrix G */
/*   of the form G=(B|Q), where B consists of the first alpha   */
/*   columns of G and Q of the last beta columns. Then the      */
/*   quaternary method computes the minimum distance of the     */
/*   quaternary code generated by (2B|Q|Q), which is twice the  */
/*   minimum distance of C.                                     */
/*                                                              */
/****************************************************************/
function Z2Z4MinimumLeeWeightWords_Quaternary(C, numWords)
    type := Z2Z4Type(C);
    G := MinRowsGeneratorMatrix(C);
    binaryG := ColumnSubmatrix(G, type[1]);
    quaternaryG := ColumnSubmatrixRange(G, type[1]+1, type[1]+type[2]);
    doubleG := HorizontalJoin(<binaryG, quaternaryG, quaternaryG>);
    LQC := LinearCode(doubleG);
    minWeightDouble := MinimumLeeDistance(LQC);

    if numWords lt 0 then
        minWordsDouble := WordsOfLeeWeight(LQC, minWeightDouble);
    else
        minWordsDouble := WordsOfLeeWeight(LQC, minWeightDouble : NumWords := numWords);
    end if;

    n := Length(C`Code);
    return {C!(Eltseq(word)[1..n]) : word in minWordsDouble};
end function;

/****************************************************************/
/*                                                              */
/* Function name: ExpectedMinWeight                             */
/* Parameters: type                                             */
/* Function description: Given the type of a Z2Z4-additive code,*/
/*   return the expected minimum weight.                        */
/* Input parameters description:                                */
/*   - type : The type of a Z2Z4-additive code                  */
/* Output parameters description:                               */
/*   - An integer representing the expected minimum weight      */
/*                                                              */
/* Comments: Probabilistic argument assuming a random code.     */
/*                                                              */
/* Function developed by A. Torres                              */
/*                                                              */
/****************************************************************/
ExpectedMinWeight := function(type)
    alpha := type[1];
    beta := type[2];
    gamma := type[3];
    delta := type[4];
    kappa := type[5];

    n := alpha+beta;
    nBin := alpha+2*beta;

    M1 := 2^(gamma-kappa+delta);
    M2 := 2^(gamma+delta)-M1;
    M3 := 2^(2*delta+gamma)-M2-M1;

    sum := 0;
    for y in [0..(nBin-1)] do
        sum +:= 2^(
            ((beta gt 0) select M1*(-beta+Log(2, 2^(beta)-&+[Binomial(beta, r)
                : r in [0..Min(Floor(y/2), beta-1)]]))
                else 0)
            + ((alpha gt 0) select M2*(-n+Log(2, 2^n-&+[&+[
                Binomial(alpha, r-2*rp)*Binomial(beta, rp)
                : rp in [Max(Ceiling((r-alpha)/2), 0)..Min(Floor(r/2), beta)]]
                : r in [0..y]]))
                else 0)
            + M3*(-nBin+Log(2, 2^nBin-&+[Binomial(nBin, r)
                : r in [0..y]])));
    end for;

    return sum;
end function;

/****************************************************************/
/*                                                              */
/* Function name: GaussianCDF                                   */
/* Parameters: x, mu, std                                       */
/* Function description: Value at point x of the cumulative     */
/*   distribution function of a normal distribution with mean   */
/*   mu and standard deviation std.                             */
/*                                                              */
/* Function developed by A. Torres                              */
/*                                                              */
/****************************************************************/
GaussianCDF := func< x, mu, var | 1/2*(1+(Erf((x-mu)/(Sqrt(2*var))))) >;

/****************************************************************/
/*                                                              */
/* Function name: ExpectedMinWeightGaussian                     */
/* Parameters: type                                             */
/* Function description: Given the type of a Z2Z4-additive code,*/
/*   return the expected minimum weight approximating the       */
/*   binomial distribution by a normal distribution.            */
/* Input parameters description:                                */
/*   - type : The type of a Z2Z4-additive code                  */
/* Output parameters description:                               */
/*   - An integer representing expected minimum weight          */
/*                                                              */
/* Comments: Probabilistic argument assuming a random code.     */
/*                                                              */
/* Function developed by A. Torres                              */
/*                                                              */
/****************************************************************/
ExpectedMinWeightGaussian := function(type)
    alpha := type[1];
    beta := type[2];
    gamma := type[3];
    delta := type[4];
    kappa := type[5];

    M1 := 2^(gamma-kappa+delta);
    M2 := 2^(gamma+delta)-M1;
    M3 := 2^(2*delta+gamma)-M2-M1;
    sum := 1;
    for y in [1..(alpha+2*beta-1)] do
        f1 := (beta gt 0) select GaussianCDF(Floor(y/2), beta/2, beta/4) else 0;
        f2 := GaussianCDF(y, alpha/2+beta, alpha/4+beta);
        f3 := GaussianCDF(y, alpha/2+beta, alpha/4+beta/2);
        if (f1 lt 1) and (f2 lt 1) and (f3 lt 1) then
            sum +:= 2^(M1*Log(2, 1-f1) + M2*Log(2, 1-f2) + M3*Log(2, 1-f3));
        end if;
    end for;

    return sum;
end function;

/****************************************************************/
/*                                                              */
/* Function name: WorkFactorBruteForce                          */
/* Parameters: C                                                */
/* Function description: Given a Z2Z4-additive code, return the */
/*   work factor corresponding to the Brute Force method.       */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/* Output parameters description:                               */
/*   - An integer having the work factor                        */
/*                                                              */
/* Function developed by A. Torres                              */
/*                                                              */
/****************************************************************/
WorkFactorBruteForce := function(C)
    return #C*Length(C);
end function;

/****************************************************************/
/*                                                              */
/* Function name: WorkFactorKernelCosets                        */
/* Parameters: C, minWeight                                     */
/* Function description: Given a Z2Z4-additive code C and an    */
/*   integer minWeight, return the work factor corresponding to */
/*   the KernelCosets method, assuming that minWeight is the    */
/*   minimum weight of C.                                       */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/*   - minWeight : an integer representing the expected minimum */
/*                 weight of C                                  */
/* Output parameters description:                               */
/*   - An integer having the work factor                        */
/*                                                              */
/* Remark: This function is based on the formulas for the work  */
/*   factors given in p. 59 of F. Zeng, "Nonlinear codes:       */
/*   representation, constructions, minimum distance            */
/*   computation and decoding", PhD Thesis, UAB, 2014.          */
/*                                                              */
/* Function developed by A. Torres                              */
/*                                                              */
/****************************************************************/
WorkFactorKernelCosets := function(C, minWeight)
    type := Z2Z4Type(C);

    // Length of the binary codewords in the Gray map image of C
    n := type[1]+2*type[2];
    // Dimension of the extended cosets with a kernel of minimum 
    // dimension, that is, having only order two codewords
    k := type[3]+type[4]+1;
    // Maximum number of cosets, when the kernel has minimum dimension
    t := 2^type[4];

    w0 := 0;
    repeat
        w0 +:= 1;
        lowerBound := Floor(n/k)*(w0+1) + Max(0, w0+1-(k-(n mod (k))));
    until lowerBound ge minWeight;

    return t*(n-k-1)*Ceiling(n/k)*&+[Binomial(k, r) : r in [0..w0]];
end function;

/****************************************************************/
/*                                                              */
/* Function name: WorkFactorQuaternary                          */
/* Parameters: C, minWeight                                     */
/* Function description: Given a Z2Z4-additive code C and an    */
/*   integer minWeight, return the work factor corresponding to */
/*   the Quaternary method, assuming that minWeight is the      */
/*   minimum weight of C.                                       */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/*   - minWeight : an integer representing the expected minimum */
/*                 weight of C                                  */
/* Output parameters description:                               */
/*   - An integer having the work factor                        */
/*                                                              */
/* Remark: This function is based on the formulas for the work  */
/*   factors given in p. 123 of G. White, "Enumeration-based    */
/*   algorithms in Coding Theory", PhD Thesis, University of    */
/*   Sidney, 2006.                                              */
/*                                                              */
/* Function developed by A. Torres                              */
/*                                                              */
/****************************************************************/
WorkFactorQuaternary := function(C, minWeight)
    type := Z2Z4Type(C);
    minWeight *:= 2; 

    // Length of the binary codewords in the Gray map image of C
    n := type[1]+2*type[2];
    beta := type[2];
    delta := type[4];
    gamma := type[3];
    lenBinInfoVector := gamma+2*delta;

    if delta eq 0 then
        return WorkFactorBruteForce(C);
    end if;

    rBrouwer := Ceiling(gamma-1+minWeight/Floor(beta/delta));
    rZimmermann := Ceiling(gamma-1+(minWeight+2*(delta-(beta mod delta)))/Ceiling(beta/delta));

    wBrouwer := 2*(n-delta)*Floor(beta/delta)
                    *&+[Binomial(lenBinInfoVector, r) : r in [1..rBrouwer]];
    wZimmermann := 2*(n-delta)*Ceiling(beta/delta)
                    *&+[Binomial(lenBinInfoVector, r) : r in [1..rZimmermann]];

    return Minimum(wBrouwer, wZimmermann);
end function;

/****************************************************************/
/*                                                              */
/* Function name: WorkFactorBrouwer                             */
/* Parameters: C, minWeight                                     */
/* Function description: Given a Z2Z4-additive code C and an    */
/*   integer minWeight, return the work factor corresponding to */
/*   the Brouwer method, assuming that minWeight is the minimum */
/*   weight of C.                                               */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/*   - minWeight : an integer representing the expected minimum */
/*                 weight of C                                  */
/* Output parameters description:                               */
/*   - An integer having the work factor                        */
/*                                                              */
/* Remark: This function is based (adapted to Z2Z4) on the      */
/*   formulas for the work factors given in p. 123 of G. White, */
/*   "Enumeration-based algorithms in Coding Theory", PhD       */
/*   Thesis, University of Sidney, 2006.                        */
/*                                                              */
/* Function developed by A. Torres                              */
/*                                                              */
/****************************************************************/
WorkFactorBrouwer := function(C, minWeight)
    type := Z2Z4Type(C);

    h := #MatricesStandardForm_Brouwer(C, type);
    lenBinInfoVector := type[3]+2*type[4];

    w0 := 0;
    repeat
        w0 +:= 1;
        lowerBound := h*(w0+1);
    until lowerBound ge minWeight;

    return (type[1]-type[5]+2*(type[2]-type[4]))*h
            *&+[Binomial(lenBinInfoVector, r) : r in [0..w0]];
end function;

/****************************************************************/
/*                                                              */
/* Function name: WorkFactorZimmermann                          */
/* Parameters: C, minWeight                                     */
/* Function description: Given a Z2Z4-additive code C and an    */
/*   integer minWeight, return the work factor corresponding to */
/*   the Brouwer-Zimmermann method, assuming that minWeight is  */
/*   the minimum weight of C.                                   */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/*   - minWeight : an integer representing the expected minimum */
/*                 weight of C                                  */
/* Output parameters description:                               */
/*   - An integer having the work factor                        */
/*                                                              */
/* Remark: This function is based (adapted to Z2Z4) on the      */
/*   formulas for the work factors given in p. 123 of G. White, */
/*   "Enumeration-based algorithms in Coding Theory", PhD       */
/*   Thesis, University of Sidney, 2006.                        */
/*                                                              */
/* Function developed by A. Torres                              */
/*                                                              */
/****************************************************************/
WorkFactorZimmermann := function(C, minWeight)
    type := Z2Z4Type(C);

    _, kList := MatricesStandardForm_Zimmermann(C, type);
    lenBinInfoVector := type[3]+2*type[4];

    w0 := 0;
    repeat
        w0 +:= 1;
        lowerBound := &+[Maximum(0, w0+1-x) : x in kList];
    until lowerBound ge minWeight;

    h := Max([1] cat [i : i in [1..#kList] | kList[i] lt w0]);
    return (type[1]-type[5]+2*(type[2]-type[4]))*h
            *&+[Binomial(lenBinInfoVector, r) : r in [0..w0]];
end function;

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
intrinsic MinimumWords(C::Z2Z4Code : NumWords := -1, Method := "Auto") -> Set 
{
Given a Z2Z4-additive code C, return the set of all codewords of C having minimum 
Lee weight. If NumWords is set to a non-negative integer, then the algorithm will 
terminate after at least that total of codewords have been found. 

Depending on the parameters of the code C, some methods  to collect the codewords
of minimum Lee weight may be faster than others. For example, sometimes, a brute 
force calculation of the entire Lee weight distribution can be a faster way for
small codes. When the parameter Method is set to the default "Auto", then the 
method is internally chosen. The user can specify which method they want to use,
setting the parameter Method to "Distribution", "KernelCosets", "Brouwer", 
"Zimmermann", or "Quaternary". 
}
    require not(zero_Code(C)): "Code C cannot be the zero code";
    
    if NumWords eq 0 then 
        return {};
    end if;
    
    type := Z2Z4Type(C);

    case Method :
        when "Distribution": 
            // By using an exhaustive search
            words := MinimumWordsBruteForce(C, NumWords);
        when "KernelCosets": 
            // By dividing C into cosets of the kernel 
            words := MinimumWordsKernelCosets(C, NumWords);
        when "Brouwer":
            isBinaryQuaternaryLinear, minWords := 
                      Z2Z4MinimumLeeWords_ForBinaryQuaternaryLinear(C, NumWords);
            if isBinaryQuaternaryLinear then
                words := minWords;
            else
                words := Z2Z4MinimumLeeWords_Brouwer(C, NumWords);
            end if;
        when "Zimmermann":
            isBinaryQuaternaryLinear, minWords := 
                      Z2Z4MinimumLeeWords_ForBinaryQuaternaryLinear(C, NumWords);
            if isBinaryQuaternaryLinear then
                words := minWords;
            else
                words := Z2Z4MinimumLeeWords_Zimmermann(C, NumWords);
            end if;
        when "Quaternary":
            words := Z2Z4MinimumLeeWeightWords_Quaternary(C, NumWords);
        else 
            if (type[1]+2*type[2]) lt MIN_GAUSSIAN_APPROX then
                expected := Ceiling(ExpectedMinWeight(type));
            else 
                expected := Ceiling(ExpectedMinWeightGaussian(type));
            end if;

            WFBF := WorkFactorBruteForce(C);
            WFK := WorkFactorKernelCosets(C, expected);
            WFB := WorkFactorBrouwer(C, expected);
            WFZ := WorkFactorZimmermann(C, expected);
            WFQ := WorkFactorQuaternary(C, expected);

            Work := [WFBF, WFK, WFB, WFZ, WFQ];
            Ind := ["BF", "K", "B", "Z", "Q"];
            ParallelSort(~Work, ~Ind);
            case Ind[1]:
                when "BF":
                    words := MinimumWordsBruteForce(C, NumWords);
                when "K":
                    words := MinimumWordsKernelCosets(C, NumWords);
                when "Q":
                    words := Z2Z4MinimumLeeWeightWords_Quaternary(C, NumWords);
                when "B":
                    isBinaryQuaternaryLinear, minWords := 
                            Z2Z4MinimumLeeWords_ForBinaryQuaternaryLinear(C, NumWords);
                    if isBinaryQuaternaryLinear then
                        words := minWords;
                    else
                        if (type[4] lt MAX_DELTA_COSETS) then
                            words := MinimumWordsKernelCosets(C, NumWords);
                        elif (type[3]-type[5] gt MAX_TORSION) or
                        expected gt MAX_EXPECTED_WEIGHT then
                            words := Z2Z4MinimumLeeWeightWords_Quaternary(C, NumWords);
                        else
                            words := Z2Z4MinimumLeeWords_Brouwer(C, NumWords);
                        end if;
                    end if;
                when "Z":
                    isBinaryQuaternaryLinear, minWords := 
                            Z2Z4MinimumLeeWords_ForBinaryQuaternaryLinear(C, NumWords);
                    if isBinaryQuaternaryLinear then
                        words := minWords;
                    else
                        if (type[4] lt MAX_DELTA_COSETS) then
                            words := MinimumWordsKernelCosets(C, NumWords);
                        elif (type[3]-type[5] gt MAX_TORSION) or
                        expected gt MAX_EXPECTED_WEIGHT then
                            words := Z2Z4MinimumLeeWeightWords_Quaternary(C, NumWords);
                        else
                            words := Z2Z4MinimumLeeWords_Zimmermann(C, NumWords);
                        end if;
                    end if;
            end case;               
    end case;

    if #words gt 0 then
        UpdateMinimumLeeWeightWord(~C, Random(words));
    end if;

    return words;

end intrinsic;

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
intrinsic MinimumWord(C::Z2Z4Code : Method := "Auto") -> ModTupFldElt  
{
Given a Z2Z4-additive code C, return one codeword of the code C having minimum 
Lee weight. 

Depending on the parameters of the code C, some methods to obtain one codeword
of minimum Lee weight may be faster than others. For example, sometimes, a brute
force calculation of the entire Lee weight distribution can be a faster way for
small codes. When the parameter Method is set to the default "Auto", then the
method is internally chosen. The user can specify which method they want to use, 
setting the parameter Method to "Distribution", "KernelCosets", "Brouwer", 
"Zimmermann", or "Quaternary".
}
    require not(zero_Code(C)): "Code C cannot be the zero code";
    
    if not IsVerbose("IgnoreWeightAttributes") then
        if assigned C`MinimumLeeWeightWord then
            return C`MinimumLeeWeightWord;
        end if;
    end if;

    type := Z2Z4Type(C);

    case Method :
        when "Distribution": 
            // By using an exhaustive search
            word := MinimumWordBruteForce(C);
        when "KernelCosets": 
            // By dividing C into cosets of the kernel 
            word := MinimumWordKernelCosets(C);
        when "Brouwer":
            isBinaryQuaternaryLinear, _, minWord := 
                             Z2Z4MinimumLeeWeight_ForBinaryQuaternaryLinear(C);
            if isBinaryQuaternaryLinear then
                word := minWord;
            else
                _, word := Z2Z4MinimumLeeWeight_Brouwer(C);
            end if;
        when "Zimmermann":
            isBinaryQuaternaryLinear, _, minWord := 
                             Z2Z4MinimumLeeWeight_ForBinaryQuaternaryLinear(C);
            if isBinaryQuaternaryLinear then
                word := minWord;
            else
                _, word := Z2Z4MinimumLeeWeight_Zimmermann(C);
            end if;
        when "Quaternary":
            word := Rep(Z2Z4MinimumLeeWeightWords_Quaternary(C, 1));
        else   
            if (type[1]+2*type[2]) lt MIN_GAUSSIAN_APPROX then
                expected := Ceiling(ExpectedMinWeight(type));
            else 
                expected := Ceiling(ExpectedMinWeightGaussian(type));
            end if;

            WFBF := WorkFactorBruteForce(C);
            WFK := WorkFactorKernelCosets(C, expected);
            WFB := WorkFactorBrouwer(C, expected);
            WFZ := WorkFactorZimmermann(C, expected);
            WFQ := WorkFactorQuaternary(C, expected);

            Work := [WFBF, WFK, WFB, WFZ, WFQ];
            Ind := ["BF", "K", "B", "Z", "Q"];
            ParallelSort(~Work, ~Ind);
            case Ind[1]:
                when "BF":
                    word := MinimumWordBruteForce(C);
                when "K":
                    word := MinimumWordKernelCosets(C);
                when "Q":
                    word := Rep(Z2Z4MinimumLeeWeightWords_Quaternary(C, 1));
                when "B":
                    isBinaryQuaternaryLinear, _, minWord := 
                               Z2Z4MinimumLeeWeight_ForBinaryQuaternaryLinear(C);
                    if isBinaryQuaternaryLinear then
                        word := minWord;
                    else
                        if (type[4] lt MAX_DELTA_COSETS) then
                            word := MinimumWordKernelCosets(C);
                        elif (type[3]-type[5] gt MAX_TORSION) or
                        expected gt MAX_EXPECTED_WEIGHT then
                            word := Rep(Z2Z4MinimumLeeWeightWords_Quaternary(C, 1));
                        else
                            _, word := Z2Z4MinimumLeeWeight_Brouwer(C);
                        end if;
                    end if;
                when "Z":
                    isBinaryQuaternaryLinear, _, minWord := 
                               Z2Z4MinimumLeeWeight_ForBinaryQuaternaryLinear(C);
                    if isBinaryQuaternaryLinear then
                        word := minWord;
                    else
                        if (type[4] lt MAX_DELTA_COSETS) then
                            word := MinimumWordKernelCosets(C);
                        elif (type[3]-type[5] gt MAX_TORSION) or
                        expected gt MAX_EXPECTED_WEIGHT then
                            word := Rep(Z2Z4MinimumLeeWeightWords_Quaternary(C, 1));
                        else
                            _, word := Z2Z4MinimumLeeWeight_Zimmermann(C);
                        end if;
                    end if;
            end case;          
    end case;

    UpdateMinimumLeeWeightWord(~C, word);
    
    return word;

end intrinsic;

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
intrinsic Z2Z4MinimumLeeWeight(C::Z2Z4Code : Method := "Auto") -> RngIntElt 
{
Given a Z2Z4-additive code C, return the minimum Lee weight of the codewords 
belonging to the code C, which is also the minimum Lee distance between any two 
codewords.

Depending on the parameters of the code C, some methods to obtain the minimum
Lee weight may be faster than others. For example, sometimes, a brute force 
calculation of the entire Lee weight distribution can be a faster way for 
small codes. When the parameter Method is set to the default "Auto", then the 
method is internally chosen. The user can specify which method they want to use, 
setting the parameter Method to "Distribution", "KernelCosets", "Brouwer", 
"Zimmermann", or "Quaternary".
}
    require not(zero_Code(C)): "Code C cannot be the zero code";

    if not IsVerbose("IgnoreWeightAttributes") then
        if assigned C`MinimumLeeWeight then
            return C`MinimumLeeWeight;
        end if;
    end if;
    
    type := Z2Z4Type(C);
    
    case Method :
        when "Distribution":
            // By using an exhaustive search
            d := MinimumLeeWeightBruteForce(C);
        when "KernelCosets": 
            // By dividing C into cosets of the kernel
            d := MinimumLeeWeightKernelCosets(C);
        when "Brouwer":
            isBinaryQuaternaryLinear, minWeight := 
                        Z2Z4MinimumLeeWeight_ForBinaryQuaternaryLinear(C);
            if isBinaryQuaternaryLinear then
                d := minWeight;
            else
                d := Z2Z4MinimumLeeWeight_Brouwer(C);
            end if;
        when "Zimmermann":
            isBinaryQuaternaryLinear, minWeight := 
                       Z2Z4MinimumLeeWeight_ForBinaryQuaternaryLinear(C);
            if isBinaryQuaternaryLinear then
                d := minWeight;
            else
                d := Z2Z4MinimumLeeWeight_Zimmermann(C);
            end if;
        when "Quaternary":
            d := Z2Z4MinimumLeeWeight_Quaternary(C);
        else
            if (type[1]+2*type[2]) lt MIN_GAUSSIAN_APPROX then
                expected := Ceiling(ExpectedMinWeight(type));
            else 
                expected := Ceiling(ExpectedMinWeightGaussian(type));
            end if;

            WFBF := WorkFactorBruteForce(C);
            WFK := WorkFactorKernelCosets(C, expected);
            WFB := WorkFactorBrouwer(C, expected);
            WFZ := WorkFactorZimmermann(C, expected);
            WFQ := WorkFactorQuaternary(C, expected);

            Work := [WFBF, WFK, WFB, WFZ, WFQ];
            Ind := ["BF", "K", "B", "Z", "Q"];
            ParallelSort(~Work, ~Ind);
            case Ind[1]:
                when "BF":
                    d := MinimumLeeWeightBruteForce(C);
                when "K":
                    d := MinimumLeeWeightKernelCosets(C);
                when "Q":
                    d := Z2Z4MinimumLeeWeight_Quaternary(C);
                when "B":
                    isBinaryQuaternaryLinear, minWeight := 
                              Z2Z4MinimumLeeWeight_ForBinaryQuaternaryLinear(C);
                    if isBinaryQuaternaryLinear then
                        d := minWeight;
                    else     
                        if (type[4] lt MAX_DELTA_COSETS) then
                            d := MinimumLeeWeightKernelCosets(C);
                        elif (type[3]-type[5] gt MAX_TORSION) or
                        expected gt MAX_EXPECTED_WEIGHT then
                            d := Z2Z4MinimumLeeWeight_Quaternary(C);
                        else
                            d := Z2Z4MinimumLeeWeight_Brouwer(C);
                        end if;
                    end if;
                when "Z":
                    isBinaryQuaternaryLinear, minWeight := 
                              Z2Z4MinimumLeeWeight_ForBinaryQuaternaryLinear(C);
                    if isBinaryQuaternaryLinear then
                        d := minWeight;
                    else 
                        if (type[4] lt MAX_DELTA_COSETS) then
                            d := MinimumLeeWeightKernelCosets(C);
                        elif (type[3]-type[5] gt MAX_TORSION) or
                        expected gt MAX_EXPECTED_WEIGHT then
                            d := Z2Z4MinimumLeeWeight_Quaternary(C);
                        else
                            d := Z2Z4MinimumLeeWeight_Zimmermann(C);
                        end if;
                    end if;
            end case;          
    end case;

    UpdateMinimumLeeWeight(~C, d);
    
    return d;

end intrinsic;

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
intrinsic Z2Z4MinimumLeeDistance(C::Z2Z4Code : Method := "Auto") -> RngIntElt 
{
Given a Z2Z4-additive code C, return the minimum Lee weight of the codewords 
belonging to the code C, which is also the minimum Lee distance between any two 
codewords.

Depending on the parameters of the code C, some methods to obtain the minimum
Lee weight may be faster than others. For example, sometimes, a brute force 
calculation of the entire Lee weight distribution can be a faster way for 
small codes. When the parameter Method is set to the default "Auto", then the 
method is internally chosen. The user can specify which method they want to use, 
setting the parameter Method to "Distribution", "KernelCosets", "Brouwer", 
"Zimmermann", or "Quaternary".
}
    return Z2Z4MinimumLeeWeight(C : Method := Method);

end intrinsic;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////                      THE WEIGHT DISTRIBUTION                   ////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/****************************************************************/
/*                                                              */
/* Function name: WeightDistributionBruteForce                  */
/* Parameters:  C                                               */
/* Function description: Determine the Lee weight distribution  */
/*   for the Z2Z4-additive code C. The distribution is returned */
/*   in the form of a sequence of tuples, where the i-th tuple  */
/*   contains the i-th weight, wi say, and the number of        */
/*   codewords having weight wi.                                */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/* Output parameters description:                               */
/*   - Sequence of tuples <Lee weight, number of codewords>     */
/*                                                              */
/****************************************************************/ 
function WeightDistributionBruteForce(C) 
    if (zero_Code(C)) then 
        return [<0, 1>];
    else
        alpha := C`Alpha;
        n := BinaryLength(C) + 1; 
        weightLeeSeq := [0^^n];

        for c in C`Code do 
            weightLeeCodeword := LeeWeight(c, alpha); 
            weightLeeSeq[weightLeeCodeword + 1] +:= 1;
        end for; 

        return [<i-1, weightLeeSeq[i]> : i in [1..n] | not IsZero(weightLeeSeq[i])];
    end if;
end function;

/****************************************************************/
/*                                                              */
/* Function name: AddWeightDistribution                         */
/* Parameters:  weightSeq, weightDistribution                   */
/* Function description: Given a sequence which has the number  */
/*   of codewords of each Lee weight, update the sequence by    */
/*   adding the values given  by weightDistribution             */  
/* Input parameters description:                                */
/*   - weightSeq : Sequence of integers                         */
/*   - weightDistribution : Sequence of tuples.                 */   
/*                                                              */
/****************************************************************/ 
procedure AddWeightDistribution(~weightSeq, weightDistribution)  
    for tuple in weightDistribution do 
        weightSeq[tuple[1] + 1] := weightSeq[tuple[1] + 1] + tuple[2]; 
    end for; 
end procedure; 

/****************************************************************/
/*                                                              */
/* Function name: WeightDistributionKernelCosets                */
/* Parameters:  C                                               */
/* Function description: Determine the Lee weight distribution  */
/*   for the Z2Z4-additive code C. The distribution is returned */
/*   in the form of a sequence of tuples, where the i-th tuple  */
/*   contains the i-th weight, wi say, and the number of        */
/*   codewords having weight wi.                                */
/* Input parameters description:                                */
/*   - C : A Z2Z4-additive code                                 */
/* Output parameters description:                               */
/*   - Sequence of tuples <Lee weight, number of codewords>     */
/*                                                              */
/****************************************************************/ 
function WeightDistributionKernelCosets(C)
    if (zero_Code(C)) then 
        weightDistribution := [ <0, 1> ];
    else
        kernel, cosetRepresentatives := KernelCosetRepresentativesZ2(C);  
        
        numCosets := #cosetRepresentatives; 
        weightDistribution := WeightDistribution(kernel);
    
        if (numCosets gt 0) then
            n := BinaryLength(C) + 1;  
            weightSeq := [0^^n];
            AddWeightDistribution(~weightSeq, weightDistribution);

            for i := 1 to numCosets do  
                weightDistributionCoset := WeightDistribution(kernel, 
                                                        cosetRepresentatives[i]); 
                AddWeightDistribution(~weightSeq, weightDistributionCoset);
            end for;
            weightDistribution := [<i-1, weightSeq[i]> : i in [1..n] | 
                                                        not IsZero(weightSeq[i])];
        end if;
    end if;

    return weightDistribution;
end function;

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
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - Sequence of tuples <Lee weight, number of codewords>     */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> SeqEnum                         */
/*                                                              */
/****************************************************************/
intrinsic LeeWeightDistribution(C::Z2Z4Code : Method := "Auto") -> SeqEnum
{
Determine the Lee weight distribution for the Z2Z4-additive code C. The 
distribution is returned in the form of a sequence of tuples, where the i-th tuple 
contains the i-th weight, wi say, and the number of codewords having weight wi. 
}
    if not IsVerbose("IgnoreWeightAttributes") then
        if assigned C`LeeWeightDistribution then
            return C`LeeWeightDistribution;
        end if;
    end if;
    // when Method is either "Distribution" or "KernelCosets"
    if Method eq "Distribution" then 
        distribution := WeightDistributionBruteForce(C);
    elif Method eq "KernelCosets" then
        distribution := WeightDistributionKernelCosets(C);
    end if;
    
    alpha := C`Alpha;
    beta := Length(C) - alpha;
    isBinaryLinear, grayMapImageCode := HasLinearGrayMapImage(C);
    
    // when Method is equal to "Auto"
    if (beta eq 0) then 
        // C is a binary linear code
        distribution := WeightDistribution(LinearBinaryCode(C));
    elif (alpha eq 0) then 
        // C is a quaternary linear code
        distribution := LeeWeightDistribution(LinearQuaternaryCode(C));
    elif (isBinaryLinear) then
        // Gray map image of C is a binary linear code
        distribution := WeightDistribution(grayMapImageCode);
    else
        // General case: C has alpha <> 0 and beta <> 0
        if (InformationRate(C) gt 0.5) then
            // Compute Lee weigth distribution by Macwilliams transform
            dualCode := Dual(C);
            dualCodeType := Z2Z4Type(dualCode);
            dimension := dualCodeType[3] + 2*dualCodeType[4];
            binaryLength := alpha + 2*beta;                         

            W := LeeWeightDistribution(dualCode);
            distribution := MacWilliamsTransform(binaryLength, dimension, 2, W);
        else
            if (#C lt MAX_CARDINAL_DISTRIBUTION) then
                // Compute Lee weight distribution by exhaustive search
                distribution := WeightDistributionBruteForce(C);
            else
                // Compute Lee weigth distribution by dividing C into kernel cosets 
                distribution := WeightDistributionKernelCosets(C);     
            end if;
        end if;
    end if;
    UpdateLeeWeightDistribution(~C, distribution);
    return distribution;
end intrinsic;

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
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - Enumerated sequence of tuples with the weight            */
/*     distribution                                             */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> SeqEnum                         */
/*                                                              */
/****************************************************************/
intrinsic DualLeeWeightDistribution(C::Z2Z4Code : Method := "Auto") -> SeqEnum
{
Determine the Lee weight distribution for the additive dual code of C. The 
distribution is returned in the form of a sequence of tuples, where the i-th tuple 
contains the i-th weight, wi say, and the number of codewords having weight wi. 
}
    if (InformationRate(C) gt 0.5) then
        return LeeWeightDistribution(Dual(C));
    else
        codeType := Z2Z4Type(C);
        dimension := codeType[3] + 2*codeType[4];
        binaryLength := BinaryLength(C);                             
        W := LeeWeightDistribution(C);
        return MacWilliamsTransform(binaryLength, dimension, 2, W);
    end if;

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: ExternalDistance                              */
/* Parameters: C                                                */
/* Function description: Determine the external distance for    */
/*   the Z2Z4-additive code C. The external distance of a code  */
/*   is the number of different nonzero weights of the dual     */
/*   code of C.                                                 */      
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - Integer with the external distance                       */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> RngIntElt                       */
/*                                                              */
/****************************************************************/
intrinsic ExternalDistance(C::Z2Z4Code : Method := "Auto") -> RngIntElt
{
Determine the external distance for the Z2Z4-additive code C. The external distance 
of a code is the number of different nonzero weights of the dual code of C.
}
    return #DualLeeWeightDistribution(C : Method := Method) - 1;

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: LeeWeightEnumerator                           */
/* Parameters: C                                                */
/* Function description: Determine the Lee weight enumerator    */
/*   for the Z2Z4-additive code C. The Lee weight enumerator of */
/*   C is the polynomial sum_{v in C}(X^(alpha+2*beta-w_L(v))*  */
/*   Y^(w_L(v))), where w_L(v) is the Lee weight of v. The      */
/*   result will lie in a global multivariate polynomial ring   */
/*   over Z with two variables. The angle-bracket notation may  */
/*   be used to assign names to the indeterminates.             */      
/* Input parameters description:                                */
/*   - C: A Z2Z4-additive code                                  */
/* Output parameters description:                               */
/*   - Multivariate polynomial                                  */
/*                                                              */
/* Signature: (<Z2Z4Code> C) -> RngIntElt                       */
/*                                                              */
/****************************************************************/
intrinsic LeeWeightEnumerator(C::Z2Z4Code : Method := "Auto") -> RngMPolElt
{
Determine the Lee weight enumerator for the Z2Z4-additive code C. The Lee weight 
enumerator of C is the polynomial sum_(v in C)(X^(alpha+2*beta - w_L(v))* Y^(w_L(v))), 
where w_L(v) is the Lee weight of v. The result will lie in a global multivariate 
polynomial ring over Z with two variables. The angle-bracket notation may be used 
to assign names to the indeterminates. 
}
    weightDistribution := LeeWeightDistribution(C : Method := Method);
    P := PolynomialRing(Integers(), 2);
    n := BinaryLength(C);
        
    return &+[tuple[2]*P.1^(n-tuple[1])*P.2^tuple[1] : tuple in weightDistribution];

end intrinsic;