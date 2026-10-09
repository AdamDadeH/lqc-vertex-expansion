# TestingExactFRW.nb

_cells: {'Input': 13, 'Output': 4}_

```mathematica
In this notebook we carry out a quick check of the exact solution? Does it actually solve the difference equation?  C=(p^2-Θ)Aexact[n,x,m]=0.
```
```mathematica
ThK[n_,m_]:= -Sqrt[n*m]((n+m))/(2);
ThD[n_]:=2 (n)^(2)
Aexact[n_,m_]:= Simplify[2 Sqrt[n*m] SeriesCoefficient[FId[(E)^(x)(1+s)/(1-s)(1-t)/(1+t)],{s,0,2n},{t,0,2m}]];
FId[x_]:=-2(Log[1+x]+Log[Gamma[1/2+I (Log[x])/(2π)]])
```
```mathematica
Psi=Table[1,{n,1,50}];
```
```mathematica
Psi[[1]]=Aexact[1,1]
Psi[[0]]=Aexact[0,1]
```
(* Out: -(16 (E)^(x) (1-4 (E)^(x)+(E)^(2 x)) (π)^(4)+((1+(E)^(x)))^(4) PolyGamma[3,(π+I Log[(E)^(x)])/(2 π)])/(((1+(E)^(x)))^(4) (π)^(4)) *)
(* Out: 0 *)
```mathematica
Psi[[2]]=Aexact[2,1];
Psi[[3]]=Aexact[3,1];
```
```mathematica
Psi[[4]]=Aexact[4,1];
```
```mathematica
Psi[[5]]=Aexact[5,1];
```
```mathematica

```
```mathematica
Diff[n_]:=ThK[n+1,n]Simplify[Aexact[n+1,1]]+ThD[n] Simplify[Aexact[n,1]] + ThK[n-1,n]Simplify[Aexact[n-1,1]];
```
```mathematica
Diff[n_]:=ThK[n+1,n]Simplify[Psi[[n+1]]]+ThD[n] Simplify[Psi[[n]]] + ThK[n-1,n]Simplify[Psi[[n-1]]]+ D[Psi[[n]],{x,2}];
```
```mathematica
Diff[3]
```
(* Out: -(1)/(4 Sqrt[3] ((1+(E)^(x)))^(6) (π)^(6))5 (-192 (E)^(x) ((-1+(E)^(x)))^(2) (1-8 (E)^(x)+(E)^(2 x)) (π)^(6)-8 ((1+(E)^(x)))^(6) (π)^(2) PolyGamma[3,(π+I Log[(E)^(x)])/(2 π)]+((1+(E)^(x)))^(6) PolyGamma[5,(π+I Log[(E)^(x)])/(2 π)])+(1)/(20 ((1+(E)^(x)))^(8) (π)^(8))Sqrt[3] (-5760 (E)^(x) (π)^(8)+921 ...[7089 chars] *)
```mathematica
FullSimplify[Diff[3]]+FullSimplify[D[Aexact[3,1],{x,2}]]
```
(* Out: (1)/(480 Sqrt[3])((1)/((π)^(10))(184 (π)^(4) PolyGamma[5,(π+I Log[(E)^(x)])/(2 π)]-40 (π)^(2) PolyGamma[7,(π+I Log[(E)^(x)])/(2 π)]+PolyGamma[9,(π+I Log[(E)^(x)])/(2 π)])-180 (-1044+887 Cosh[x]-84 Cosh[2 x]+Cosh[3 x]) (Sech[(x)/(2)])^(8) (Tanh[(x)/(2)])^(2))+(1)/(480 Sqrt[3])(-(1)/((π)^(10))(184 (π) ...[518 chars] *)
```mathematica
Diff[n_]:=ThK[n+1,n]Simplify[Aexact[n+1,1]]+ThD[n] Simplify[Aexact[n,1]] + ThK[n-1,n]Simplify[Aexact[n-1,1]];
```
