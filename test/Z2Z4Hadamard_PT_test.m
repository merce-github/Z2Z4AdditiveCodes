// import "src/Z2Z4CodeConstructions.m": Z4HadamardMatrix;
// import "src/Z2Z4CodeConstructions.m": HadamardCodeZ4Matrix;
// alphaMatrix := recformat<Alpha:RngIntElt, Matrix:Mtrx>;

timeRec := [];
timeDir := [];

for m in [0..23] do
    time1 := 0.0;
    time2 := 0.0;
    for delta in [1..Floor((m+1)/2)] do
        t1 := Cputime();
        G1 := HadamardCodeZ4Matrix(delta-1, m);
        time1 +:= Cputime(t1);

        t2 := Cputime();
        G2 := Z4HadamardMatrix(delta, m);
        time2 +:= Cputime(t2);

        // assert G1 eq G2;
    end for;
    Append(~timeRec, time1);
    Append(~timeDir, time2);
end for;

for i in [1..24] do
    print i-1, timeRec[i], timeDir[i];
end for;