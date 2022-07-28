///////////////////////////////////////////////////////////////////////////////
/////////    Copyright 2010-2019 D. Molinero, J. Pujol, J. Rifà         ///////
/////////                   and M. Villanueva                           ///////
/////////                                                               ///////
/////////    This program is distributed under the terms of GNU         ///////
/////////               General Public License                          ///////
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
/* File name: Z2Z4CoveringRadius.m                           */
/*                                                           */
/* Comment: Package developed within the CCSG group          */
/*                                                           */
/* Authors: Daniel Molinero, Jaume Pujol, Josep Rifà and     */
/*          Mercè Villanueva                                 */
/*                                                           */
/* Revision version and last date: 1.0, 2010/03/29           */
/*                                 1.1, 2010/08/15           */
/*                                 1.2, 2010/09/24           */
/*                                 1.3, 2010/11/04           */
/*                                 1.4  2019/03/18           */
/*                                 1.5  2022/07/24           */
/*                                                           */
/*************************************************************/
//Uncomment freeze when package finished
freeze;

intrinsic Z2Z4CoveringRadius_version() -> SeqEnum
{Return the current version of this package.}
    version := [1, 5];
    return version;
end intrinsic;

/******************************************************************
    GLOBAL VARIABLES
*******************************************************************/
Z4 := Integers(4);

//Maximum number of elements in a sequence for the current MAGMA distribution
//This number can be changed if your distribution allows longer sequences  
MAXSEQLENGTH := 2^28;

import "Z2Z4AdditiveCodes.m": Z2Z4ChangeMatrixZ4toZ2;
import "Z2Z4CodeConstructions.m": zero_Code;
import "Z2Z4Decode.m": MultiplyByG; 

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////                COVERING RADIUS                                  ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

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
intrinsic CoveringRadiusBounds(C::Z2Z4Code) -> RngIntElt, RngIntElt
{
Return a lower and upper bounds on the covering radius of the Z2Z4-additive code C.
The lower bound is given by the maximum between the error correcting capability of
C and the covering radius of the linear span of Cbin, where Cbin=Phi(C) and Phi is 
the Gray map. The upper bound is given by the minimum between the external distance 
of C and the covering radius of Kbin=Phi(K_C), where K_C is the kernel of C. 
Note that this function is only applicable when C is small.
}
    errorCapability := Floor((Z2Z4MinimumLeeWeight(C)-1)/2);
    _, spanZ2Code := SpanZ2Code(C);
    _, kernelZ2Code := KernelZ2Code(C);
    lowerBound := Max(errorCapability, CoveringRadius(spanZ2Code)); 
    upperBound := Min(ExternalDistance(C), CoveringRadius(kernelZ2Code));
    
    return lowerBound, upperBound;  

end intrinsic;

intrinsic CoveringRadiusBruteForce(C::Z2Z4Code : MaximumTime := 0) -> RngIntElt
{}
	cw_code := SetToSequence(Set(C));
	code_cardinal := #C;
	alpha := C`Alpha;
	beta := Length(C) - alpha;
	
	universe := Z2Z4AdditiveUniverseCode(alpha, beta);
	cov_radius := 0;

    for cw_universe in universe`Code do
        minDist := LeeDistance(cw_universe, Random(C), alpha);
        for cw_code in C`Code do
            dist := LeeDistance(cw_universe, cw_code, alpha);
            if dist lt minDist then
                minDist := dist;
            end if;
        end for;
        if minDist gt cov_radius then
            cov_radius := minDist;
        end if;
    end for;
    
	return cov_radius;

end intrinsic;

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
intrinsic CoveringRadius(C::Z2Z4Code : MaximumTime := 0) -> RngIntElt
{
The covering radius of the Z2Z4-additive code C, which  is the smallest 
radius rho (considering the Lee distance) such that the spheres of radius
rho centered at all codewords of C cover the ambient space V=Z2^alpha x
Z4^beta. Note that this function is only applicable when C is small.

The parameter MaximumTime sets a time limit (in seconds of "user time") 
after which the calculation is aborted. The default value is infinity, 
when there is no restriction on time.
}
    if assigned(C`CoveringRadius) then
        return C`CoveringRadius;
    end if;

    if zero_Code(C) then
        return BinaryLength(C);
    end if;

    isLinearCbin, Cbin := HasLinearGrayMapImage(C);
    if isLinearCbin then 
        rho := CoveringRadius(Cbin);
        C`CoveringRadius := rho;
        return rho;
    end if;
    
    requirege MaximumTime, 0;  
    iniTime := Cputime();
	limitTime := (MaximumTime ne 0);
    
    alpha := C`Alpha;
    beta := Length(C`Code) - alpha;
    nbin := alpha + 2*beta; 
    totalNumSyndromes := 2^nbin / #C;
    require (totalNumSyndromes le MAXSEQLENGTH): "Code C has too many cosets";    
    
    H := Transpose(MinRowsParityCheckMatrix(C));
    universeCode := Z2Z4AdditiveUniverseCode(alpha, beta);
    grayMap := GrayMap(universeCode);     
    V := UniverseCode(GF(2), nbin);
  
    // Sequence that will have all syndromes 
    allSyndromes := [ C!0 * H ]; 

    if assigned C`MinimumLeeWeight then
        errorCapability := Floor((C`MinimumLeeWeight-1)/2);
    else
        errorCapability := Floor((Z2Z4MinimumLeeWeight(C)-1)/2);
    end if;
    
    rho := errorCapability;
    for i in [1..errorCapability] do
        if ((limitTime) and (Cputime(iniTime) gt MaximumTime)) then
		    return rho;
		end if;
    
        // All binary vectors of weight i
        vectorsWeight_i := Setseq(Words(V, i));
        numVectorsWeight_i := #vectorsWeight_i;    
        
        // Compute all syndromes up to the error-correcting capability
        for j in [1..numVectorsWeight_i] do
            leaderZ4 := vectorsWeight_i[j]@@grayMap;
            s := MultiplyByG(leaderZ4, alpha, H);
            Append(~allSyndromes, s);
        end for;        
    end for;
    
    i := errorCapability + 1;
    while ((#allSyndromes lt totalNumSyndromes) and (i le nbin)) do
        if ((limitTime) and (Cputime(iniTime) gt MaximumTime)) then
		    return rho;
		end if;
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
                rho := i;
            end if;
            j := j + 1;
        end while;
        
        i := i + 1;
    end while;

    C`CoveringRadius := rho;
    return rho;

end intrinsic;

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
intrinsic IsCompletelyRegular(C::Z2Z4Code : MaximumCosetWeight := 0,
                              MaximumTime := 0) -> BoolElt, SeqEnum,
                                                   SeqEnum, RngIntElt
{
Given a Z2Z4-additive code C of type (alpha, beta; gamma, delta; kappa),
return true if and only if Cbin is completely regular, where Cbin=Phi(C) 
and Phi is the Gray map. If return true, then it also returns the 
intersection array, a sequence with the number of coset leaders of Cbin 
of each weight, and the covering radius. Note that, if no optional parameters 
are specified, then this function is only applicable when C is small. 

A code Cbin subset of Z2^(alpha+2beta) is called completely
regular if, for all i>=0, every vector x in C(i) has the same number c(i) of
neighbours in C(i-1) and the same number b(i) of neighbours in C(i+1), where
C(i)=[x in Z2^(alpha+2beta) | d(x,Cbin)=i]. The intersection array is given as
a sequence with the parameters [ [b(0),...,b(rho-1)], [c(1),...,c(rho)] ], where
rho is the covering radius of Cbin.

The parameter MaximumCosetWeight sets the maximum weight of the coset leaders
which are used in the calculation to check whether the code is completely
regular or not. In this case, the function returns true if and only if for all
0<=i<=MaximumCosetWeight, every vector x in C(i) has the same number c(i) of
neighbours in C(i-1) and the same number b(i) of neighbours in C(i+1). The
intersection array is given as a sequence with the parameters [ [b(0),...,b(m-1)],
[c(1),...,c(m)] ], where m=MaximumCosetWeight. The sequence with the number of coset
leaders of Cbin of each weight is given until MaximumCosetWeight. The covering
radius is not given unless MaximumCosetWeight>=rho. The default value is rho.

The parameter MaximumTime sets a time limit (in seconds of "user time") after
which the calculation is aborted. According to this parameter, there is a
MaximumCosetWeight depending on the maximum weight of the computed coset
leaders. In this case, the function returns the same as if this parameter
MaximumCosetWeight had been specified. The default value is infinity, when
there is no restriction on time.
}
	requirege MaximumCosetWeight, 0;
    requirege MaximumTime, 0;
	
    iniTime := Cputime();
	limitTime := (MaximumTime ne 0);

	intersectionArray := [ [], [] ];
    alpha := C`Alpha;
    beta := Length(C`Code) - alpha;
    f := GrayMap(Z2Z4AdditiveUniverseCode(alpha, beta));
    binaryLength := BinaryLength(C);
	maxCosetWeight := ( (MaximumCosetWeight eq 0) or
                        (MaximumCosetWeight gt binaryLength))
                       select binaryLength else MaximumCosetWeight;
    V := UniverseCode(GF(2), binaryLength);

    // sequence with the sindromes of binary vectors of weight one
    parityCheckMatrixOriginal := ParityCheckMatrix(C);
    parityCheckMatrix := Transpose(parityCheckMatrixOriginal);
    eVectorsSyndromes := [];
    for i := 1 to alpha do
        Append(~eVectorsSyndromes, parityCheckMatrix[i]);
    end for;
    for j := alpha + 1 to alpha + beta do
        Append(~eVectorsSyndromes, parityCheckMatrix[j]);
        Append(~eVectorsSyndromes, 3*parityCheckMatrix[j]);
    end for;
	
    totalNumSyndromes := (2^binaryLength) / #C;
	cosetLeadersWeightDistribution := [];

	//errorCapability := Floor((MinimumLeeDistance(C)-1)/2);
    errorCapability := 0;

    // Parity check matrix replacing twos in the alpha coordinates by ones
    Z2Z4ChangeMatrixZ4toZ2(~parityCheckMatrixOriginal, alpha);
    parityCheckMatrix := Transpose(parityCheckMatrixOriginal);

    // sequence that will have all syndromes of the binary code
    allSyndromes := [ C!0 * parityCheckMatrix ];

    // syndromes of the previous subconstituent C(i-1)
    if (errorCapability eq 0) then
        prevSyndromes := [ C!0 * parityCheckMatrix ];
    else
        prevSyndromes := [];
    end if;

    // syndromes of the current subconstituent C(i)
    classSyndromes := [];

	i := 1;
	computing := (i le maxCosetWeight);
	while ((i le errorCapability) and (computing)) do

        // all binary vectors of weight i
        vectorsWeight_i := Setseq(Words(V,i));
        numVectorsWeight_i := #vectorsWeight_i;

        // the parameter c that should have all coset leaders in C(i)
        intersectionArray[2][i]:= i;	
		
        // compute all coset leaders in C(i), which have the same parameter c
        for j in [1..numVectorsWeight_i] do
            s := vectorsWeight_i[j]@@f * parityCheckMatrix;
            Append(~allSyndromes, s);
            //Append(~classSyndromes, s);
            if (i eq errorCapability) then
                Append(~prevSyndromes,s);
            end if;
        end for;

        // the parameter b that should have all coset leaders in C(i-1)
        intersectionArray[1][i]:= binaryLength - i + 1;
		
		Append(~cosetLeadersWeightDistribution, numVectorsWeight_i);
        //Append(~cosetLeadersWeightDistribution, #classSyndromes);
        //prevSyndromes := classSyndromes;
        //classSyndromes := [];
		i := i+1;
		computing := (i le maxCosetWeight);
    end while;

    i := errorCapability +1;
    computing := (i le maxCosetWeight);	
    while ((#allSyndromes lt totalNumSyndromes) and (computing)) do

        // all binary vectors of weight i
        vectorsWeight_i := Setseq(Words(V,i));
        numVectorsWeight_i := #vectorsWeight_i;

        // find the first new syndrome in the current subconstituent C(i)
        j := 1;
        foundNewSyndrome := false;
        while (not foundNewSyndrome) do
            s := vectorsWeight_i[j]@@f * parityCheckMatrix;
            if (s notin allSyndromes) then
                foundNewSyndrome := true;
                Append(~allSyndromes, s);
                Append(~classSyndromes, s);
            end if;
            j := j + 1;
        end while;

        // the parameter c that should have all coset leaders in C(i)
        Append(~intersectionArray[2], 0);
	
        for v in eVectorsSyndromes do
			// time control
		    if ((limitTime) and (Cputime(iniTime) gt MaximumTime)) then
			     intersectionArray := [intersectionArray[1][1..i-1],
                                       intersectionArray[2][1..i-1]];
				 return true, intersectionArray, cosetLeadersWeightDistribution, _;
			end if;
            if ((s + v) in prevSyndromes) then
                intersectionArray[2][i] := intersectionArray[2][i] + 1;
            end if;
        end for;

        // check whether all coset leaders in C(i) have the same parameter c
        while ((j le numVectorsWeight_i) and
               (#allSyndromes lt totalNumSyndromes)) do
            s := vectorsWeight_i[j]@@f * parityCheckMatrix;
            if (s notin allSyndromes) then
                Append(~allSyndromes, s);
                Append(~classSyndromes, s);

                intersectionParameters := 0;
					
                for v in eVectorsSyndromes do
				    // time control
					if ((limitTime) and (Cputime(iniTime) gt MaximumTime)) then
						intersectionArray := [intersectionArray[1][1..i-1],
                                              intersectionArray[2][1..i-1]];
						return true, intersectionArray, cosetLeadersWeightDistribution, _;
					end if;
                    if ((s + v) in prevSyndromes) then
                        intersectionParameters := intersectionParameters + 1;
                    end if;
                end for;

                if (intersectionParameters ne intersectionArray[2][i]) then
                    return false, _, _, _;
                end if;
            end if;

            j := j + 1;
        end while;

        // the parameter b that should have all coset leaders in C(i-1)
        Append(~intersectionArray[1], 0);
        for v in eVectorsSyndromes do
			// time control
			if ((limitTime) and (Cputime(iniTime) gt MaximumTime)) then
			    intersectionArray := [intersectionArray[1][1..i-1],
                                      intersectionArray[2][1..i-1]];
				return true, intersectionArray, cosetLeadersWeightDistribution, _;
			end if;
            if ((prevSyndromes[1] + v) in classSyndromes) then
                intersectionArray[1][i] := intersectionArray[1][i] + 1;
            end if;
        end for;

        // check whether all coset leaders in C(i) have the same parameter b
        for j in [2..#prevSyndromes] do
            intersectionParameters := 0;
            for v in eVectorsSyndromes do
				// time control
			    if ((limitTime) and (Cputime(iniTime) gt MaximumTime)) then
					intersectionArray := [intersectionArray[1][1..i-1],
                                          intersectionArray[2][1..i-1]];
					return true, intersectionArray, cosetLeadersWeightDistribution, _;
			    end if;
                if ((prevSyndromes[j] + v) in classSyndromes) then
                    intersectionParameters := intersectionParameters + 1;
                end if;
            end for;

            if (intersectionParameters ne intersectionArray[1][i]) then
                return false, _, _, _;
            end if;
        end for;

        Append(~cosetLeadersWeightDistribution, #classSyndromes);

        prevSyndromes := classSyndromes;
        classSyndromes := [];
        i := i + 1;
		computing := (i le maxCosetWeight);
    end while;

    if (#allSyndromes eq totalNumSyndromes) then
        return true, intersectionArray,
               cosetLeadersWeightDistribution, #intersectionArray[1];
    else
        return true, intersectionArray, cosetLeadersWeightDistribution, _;
    end if;

end intrinsic;
