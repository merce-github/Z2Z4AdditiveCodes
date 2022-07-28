Z4 := Integers(4);
SetVerbose("IgnoreWeightAttributes", 1);

// /*************************************************************/
// print "test 1: A Z2Z4-additive code with alpha, beta even, p=Id.
//                alpha = 10, beta = 20, gamma = 6, delta = 2, kappa = 1, length = 30, #C = 1024";
// M := Matrix(Z4,[[2,0,2,0,0,0,2,0,2,2,0,0,0,1,0,0,0,3,2,1,2,0,1,2,0,1,2,1,0,3],
//                 [0,2,2,0,0,2,0,0,2,2,1,1,1,1,0,1,0,3,2,0,0,2,2,0,1,0,3,2,1,1],
//                 [0,0,0,0,2,0,2,2,2,0,0,0,0,1,0,0,0,1,2,3,0,0,1,0,2,3,0,3,2,3],
//                 [0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,2,2,0,0,2,0,2,2,0,0,2,2],
//                 [0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2,0,2,2,2,2,2,0,0,2,2,2,0],
//                 [0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,2,0,2,0,0,0,2,2,0,2,2,2,2],
//                 [0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,0,2,0,0,2,0,0,2,0,2,0,2],
//                 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,2,0,2,0,0,2,2,0,2,0,0,2,2],
//                 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,2,0,2,2,2,0,2,0,2,2,0,0],
//                 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,0,0,0,0,0,0,2,2,0,2,2,2]]);
// alpha := 10;
// beta := 20;
// gamma := 6;
// delta := 2;
// C := Z2Z4AdditiveCode(M, alpha);

// expectedOutputMinimumLeeDistance :=  10;

// t1 := Cputime();
// outputBruteForce := Z2Z4MinimumLeeWeight(C : Method := "Distribution");
// assert outputBruteForce eq expectedOutputMinimumLeeDistance;
// print "Brute Force:", Cputime(t1);

// t2 := Cputime();
// outputBrouwer := Z2Z4MinimumLeeWeight(C : Method := "Brouwer");
// assert outputBrouwer eq expectedOutputMinimumLeeDistance;
// print "Brouwer:", Cputime(t2);

// t3 := Cputime();
// outputZimmermann := Z2Z4MinimumLeeWeight(C : Method := "Zimmermann");
// assert outputZimmermann eq expectedOutputMinimumLeeDistance;
// print "Zimmermann:", Cputime(t3);

// /*************************************************************/
// print "test 2: A Z2Z4-additive code with alpha, beta even, p=Id.
//                alpha = 10, beta = 20, gamma = 2, delta = 6, kappa = 1, length = 30, #C = 1024";
// C := RandomZ2Z4AdditiveCode(1,30,2,8,1);

// t1 := Cputime();
// outputBruteForce := Z2Z4MinimumLeeWeight(C: Method := "Distribution"); 
// print "Brute Force:", Cputime(t1);

// t2 := Cputime();
// outputBrouwer := Z2Z4MinimumLeeWeight(C : Method := "Brouwer");
// assert outputBrouwer eq outputBruteForce;
// print "Brouwer:", Cputime(t2);

// t3 := Cputime();
// outputZimmermann := Z2Z4MinimumLeeWeight(C : Method := "Zimmermann");
// assert outputZimmermann eq outputBruteForce;
// print "Zimmermann:", Cputime(t3);

// M:=Matrix(Integers(4),[[2,0,2,2,0,2,0,0,0],[0,2,0,0,2,0,2,0,0],[0,0,0,0,2,2,2,0,0],[0,0,2,0,1,3,1,1,0],[0,0,2,0,2,1,0,0,1]]);

// G:=Matrix(Integers(4),[[2,2,2,2,2,0,0,0,0,0,2,2,0,2,2,0,2,2,0,0,0,2,2,0,0,0,0,2,0,0],
// [0,0,0,0,0,0,0,0,0,0,2,2,0,0,2,0,0,0,0,0,2,0,0,2,0,0,0,0,2,0],
// [0,0,0,0,0,0,0,0,0,0,2,0,0,2,2,0,0,2,0,0,0,0,2,0,2,0,0,0,0,0],
// [0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,2,2,2,0,0,2,0,0,2,0,0,2,0,0,0],
// [0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,2,0,2,0,0,0,2,2,2,0,2,0,0,0,0],
// [0,0,0,0,0,0,0,0,0,0,2,2,0,2,0,0,0,0,0,0,2,2,2,0,0,0,0,2,0,2],
// [2,0,0,0,0,2,2,0,0,2,1,3,0,0,0,0,0,0,1,0,1,3,2,1,1,3,3,3,3,2],
// [0,0,2,2,0,0,2,0,2,0,2,0,0,0,0,0,0,0,0,1,0,3,1,1,0,3,0,1,0,1]]);


// WorkFactorKernelCosets := function(n, k, t, w)
//     lb := 0;
//     w0 := 0;
//     while lb lt w do
//         w0+:=1;
//         lb := Floor(n/(k+1))*(w0+1) + Max(0, w0+1-(k+1-( n mod (k+1))));
//     end while;
//     print w0, lb, w;
//     return t*(n-k-1)*Ceiling(n/(k+1))*&+[ Binomial(k+1, r) : r in [1..w0]];
// end function;
/****************************************************************/
print "Move t1";
timeBF := [];
timeK := [];
timeB := [];
timeZ := [];
timeBO := [];
timeQ := [];
timeAuto := [];
kernelDim := [];
cosetLead := [];
linear := [];
weight := [];
WFKernel := [];
WFZimmermann := [];
WFQuaternary := [];

p := 3;
s := 2;
n := 20;

range := &cat[[i^^5] : i in [31..35]];
// range := &cat[[i^^5] : i in [0..15]];
// range := [1..1];
// range := &cat[[i^^5] : i in [2..10]];
// range := &cat[[i^^5] : i in [5..25]];
// range := &cat[[i^^5] : i in [10..100 by 10]];

// import "src/Z2Z4MinimumWeight.m": KernelCosetRepresentativesZ2;
// import "src/Z2Z4MinimumWeight.m": Z2Z4MinimumLeeWeight_Zimmermann_Old;
// import "src/Z2Z4MinimumWeight.m": WorkFactorKernelCosets;
// import "src/Z2Z4MinimumWeight.m": WorkFactorZimmermann;

for d in range do
    // C := RandomZ2Z4AdditiveCode(0,30,0,d,0);
    // C := RandomZ2Z4AdditiveCode(30,25,2,d,2);
    // C := RandomZ2Z4AdditiveCode(10,20,6,5,2);
    // C := RandomZ2Z4AdditiveCode(50,20,d,13,d);
    // C := RandomZ2Z4AdditiveCode(2,d,2,10,2);
    C := RandomZ2Z4AdditiveCode(d,25,2,13,2);
    // C := RandomZ2Z4AdditiveCode(4,25, 2,14, 2);
    // C := RandomZ2Z4AdditiveCode(20,d,4,10,4);
    // C := RandomZ2Z4AdditiveCode(20,100,4,d,4);
    // C := Z2Z4HadamardCode(5,10);

    // t1 := Cputime();
    // outputBruteForce := Z2Z4MinimumLeeWeight(C: Method := "Distribution");
    // Append(~timeBF, Cputime(t1));

    // t2 := Cputime();
    // outputKernel := Z2Z4MinimumLeeWeight(C: Method := "KernelCosets");
    // // assert outputKernel eq outputBruteForce;
    // Append(~timeK, Cputime(t2));

    // t3 := Cputime();
    // outputBrouwer := Z2Z4MinimumLeeWeight(C: Method := "Brouwer");
    // assert outputBrouwer eq outputKernel;
    // Append(~timeB, Cputime(t3));

    t4 := Cputime();
    outputZimmermann := Z2Z4MinimumLeeWeight(C: Method := "Zimmermann");
    // assert outputZimmermann eq outputKernel;
    Append(~timeZ, Cputime(t4));

    // t5 := Cputime();
    // outputBrouwerOld := Z2Z4MinimumLeeWeight_Brouwer(C);
    // assert outputBrouwerOld eq outputKernel;
    // Append(~timeBO, Cputime(t5));

    t5 := Cputime();
    outputZimmermannO := Z2Z4MinimumLeeWeight(C: Method := "Quaternary");
    // assert outputZimmermannO eq outputZimmermann;
    Append(~timeQ, Cputime(t5));

    // t6 := Cputime();
    // outputAuto := Z2Z4MinimumLeeWeight(C);
    // assert outputAuto eq outputZimmermann;
    // Append(~timeAuto, Cputime(t6));

    K, leaders := KernelCosetRepresentativesZ2(C);
    // type := Z2Z4Type(C);
    // print WorkFactorKernelCosets(type[1]+2*type[2], Dimension(K), #leaders,outputKernel);
    Append(~kernelDim, Dimension(K));
    Append(~cosetLead, #leaders);
    Append(~linear, HasLinearGrayMapImage(C));
    Append(~weight, outputZimmermannO);
    Append(~WFKernel, WorkFactorKernelCosets(C, outputZimmermann));
    Append(~WFZimmermann, WorkFactorZimmermann(C, outputZimmermann));
    Append(~WFQuaternary, WorkFactorQuaternary(C, 2*outputZimmermann));
    print d;

end for;

for i in [1..#range] do
    print range[i], timeZ[i], timeQ[i], WFZimmermann[i], WFQuaternary[i], weight[i];
end for;

/****************************************************************/
// print "Move beta";
// timeBF := [];
// timeK := [];
// timeB := [];
// timeZ := [];
// timeZO := [];
// kernelDim := [];
// cosetLead := [];
// linear := [];

// p := 3;
// s := 2;
// n := 20;

// range := &cat[[i,i] : i in [15..25 by 5]];

// // import "src/Z2Z4MinimumWeight.m": KernelCosetRepresentativesZ2;
// // import "src/Z2Z4MinimumWeight.m": Z2Z4MinimumLeeWeight_Zimmermann_Old;

// for d in range do
//     C := RandomZ2Z4AdditiveCode(3,d,2,10,1);
//     // C := Z2Z4HadamardCode(10,25);

//     // t1 := Cputime();
//     // outputBruteForce := Z2Z4MinimumLeeWeight(C: Method := "Distribution");
//     // Append(~timeBF, Cputime(t1));

//     t2 := Cputime();
//     outputKernel := Z2Z4MinimumLeeWeight(C: Method := "KernelCosets");
//     // assert outputKernel eq outputBruteForce;
//     Append(~timeK, Cputime(t2));

//     t3 := Cputime();
//     outputBrouwer := Z2Z4MinimumLeeWeight(C: Method := "Brouwer");
//     assert outputBrouwer eq outputKernel;
//     Append(~timeB, Cputime(t3));

//     t4 := Cputime();
//     outputZimmermann := Z2Z4MinimumLeeWeight(C: Method := "Zimmermann");
//     assert outputZimmermann eq outputKernel;
//     Append(~timeZ, Cputime(t4));

//     // t5 := Cputime();
//     // outputZimmermannO := Z2Z4MinimumLeeWeight_Zimmermann_Old(C);
//     // assert outputZimmermannO eq outputZimmermann;
//     // Append(~timeZO, Cputime(t5));

//     K, leaders := KernelCosetRepresentativesZ2(C);

//     Append(~kernelDim, Dimension(K));
//     Append(~cosetLead, #leaders);
//     Append(~linear, HasLinearGrayMapImage(C));

//     print d, Cputime(t2), Cputime(t3), Cputime(t4), Dimension(K), #leaders;
// end for;

// for i in [1..#range] do
//     print range[i], timeK[i], timeB[i], timeZ[i], kernelDim[i], cosetLead[i];
// end for;

SetVerbose("IgnoreWeightAttributes", 0);
