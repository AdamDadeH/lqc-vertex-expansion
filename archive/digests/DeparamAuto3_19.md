# DeparamAuto3_19.nb

_cells: {'Subsubtitle': 8, 'Input': 12, 'Output': 5, 'Text': 3}_

```mathematica
Fsqth[x_]:=2I (x)/(1+x) - (1)/(π)PolyGamma[1/2+I (Log[x])/(2π)]
FId[x_]:=-2(Log[1+x]+Log[Gamma[1/2+I (Log[x])/(2π)]])
FTheta[x_]:=(4 (π)^(2) x-((1+x))^(2) PolyGamma[1,(π+I Log[x])/(2 π)])/(2 (π)^(2) ((1+x))^(2))
```
```mathematica
Simplify[-I x D[2I (x)/(1+x) - (1)/(π)PolyGamma[1/2+I (Log[x])/(2π)],x]]
```
(* Out: (4 (π)^(2) x-((1+x))^(2) PolyGamma[1,(π+I Log[x])/(2 π)])/(2 (π)^(2) ((1+x))^(2)) *)
```mathematica
Sqth[n_,m_]:= 2 Sqrt[n*m] SeriesCoefficient[Fsqth[(1+s)/(1-s)(1-t)/(1+t)],{s,0,2n},{t,0,2m}];
Id[n_,m_]:= 2 Sqrt[n*m] SeriesCoefficient[FId[(1+s)/(1-s)(1-t)/(1+t)],{s,0,2n},{t,0,2m}];

Theta[n_,m_]:=2 Sqrt[n*m] SeriesCoefficient[FTheta[(1+s)/(1-s)(1-t)/(1+t)],{s,0,2n},{t,0,2m}];
Aexact[n_,m_]:= 2 Sqrt[n*m] SeriesCoefficient[FId[(E)^(x)(1+s)/(1-s)(1-t)/(1+t)],{s,0,2n},{t,0,2m}];
```
```mathematica
Clear[x,y]
x=5
y=6
(-2FullSimplify[Theta[x,y]])/((x+y)Sqrt[x*y])
```
(* Out: 5 *)
(* Out: 6 *)
(* Out: 1 *)
```mathematica
degen[v_]:=Table[Length[Split[Sort[v]][[i]]],{i,1,Length[Split[Sort[v]]]}];
w[v_]:=Table[Split[Sort[v]][[i]][[1]],{i,1,Length[Split[Sort[v]]]}];
c[w_]:=Table[N[Sqth[w[[i]],w[[i]]]],{i,1,Length[w]}];
```
> The following produces the diagonal matrix elements corresponding to the distinct values of volume.

```mathematica

```
```mathematica
IPartBase[w_,deg_]:=Product[(1)/((deg[[i]]-1)!),{i,1,Length[w]}]Sum[(E)^(I*d[[i]]*x)Product[(1)/(d[[i]]-d[[j]]),{j,1,i-1}]* Product[(1)/(d[[i]]-d[[j]]),{j,i+1,Length[w]}],{i,1,Length[w]}];
IPartDer[w_,deg_]:=Fold[D[#1,#2]&,IPartBase[w,deg],Table[{d[[i]],deg[[i]]-1},{i,1,Length[w]}]];
IPart[w_,deg_,c_]:=Evaluate[IPartDer[w,deg]]/.d->c
ProdPart[v_]:=Product[N[Sqth[v[[i]],v[[i+1]]]],{i,1,Length[v]-1}];
Amplitude[v_]:=( wv=w[v]; deg=degen[v]; diagl=c[wv];  ProdPart[v]*IPart[wv, deg, diagl]);
```
> Changes generates a set of ordered changes given a unorded set called base.
Paths then generates the actual sets of paths here starting from volume = 4
ManyAmp then computes the amplitude for a whole set of paths
ManyChange compute the amplitude corresponding to many different choice of base.

```mathematica
Changes[base_]:=(
changes=Permutations[base];
)
Paths[changes_]:=(
numpaths=Length[changes];
m=Length[changes[[1]]];
paths=Table[1,{j,1,numpaths}];
path=Table[1,{i,1,m+1}];
k=1;
For[j=1,j<=numpaths,j++,
valid=1;
For[i=1,i<=m,i++,
path[[i+1]]=path[[i]]+changes[[j]][[i]];
If[path[[i+1]]<1,valid=0];];
If[valid==1,
paths[[k]]=path;
k++;];
];
numpaths=k-1;);

ManyAmp[paths_]:=(
Many=0;
For[i2=1,i2<=numpaths,i2++,
Many=Many+Amplitude[paths[[i2]]];
];)
ManyChange[bases_]:=(
ManyC=0;
numchange=Length[bases];
For[i1=1,i1<=numchange,i1++,
Changes[bases[[i1]]];
Paths[changes];
ManyAmp[paths];
ManyC=ManyC+Many;];
)
```
> Finally Generate base gives every possible set of unordered changes that go a distance m from the initial point (for same volume same volume transitions).

```mathematica
GenerateBase[m_,n_]:=(
part=IntegerPartitions[m,n];
len1=Length[part];
len2=Length[part]^2;
bases=Table[1,{i,1,len2}];
For[i=1,i<=len1,i++,
For[j=1,j<=len1,j++,
bases[[(j-1)*len1+i]]=Join[part[[i]],-part[[j]]];
];
];
)
```
```mathematica
Generate base with M changes where greatest deviation is D
```
(* Out: base changes D deviation Generate greatest is M where with *)
```mathematica
Generate2[M_,D_]:=(
basetot={};
For[d=0,d<=D,d++,
For[i=1,i<=M-1,i++,
part1=IntegerPartitions[d,{i}];
part2=IntegerPartitions[d,{M-i}];
len1=Length[part1];
len2=Length[part2];
len12=len1*len2;
bases=Table[1,{j,1,len12}];
For[k=1,k<=len1,k++,
For[l=1,l<=len2,l++,
bases[[(k-1)*len2+l]]=Join[part1[[k]],-part2[[l]]];
];
];
basetot=Join[basetot,bases];
];
];)
```
```mathematica
Generate3[M_,d_]:=(
basetot={};
For[i=1,i<=M-1,i++,
part1=IntegerPartitions[d,{i}];
part2=IntegerPartitions[d,{M-i}];
len1=Length[part1];
len2=Length[part2];
len12=len1*len2;
bases=Table[1,{j,1,len12}];
For[k=1,k<=len1,k++,
For[l=1,l<=len2,l++,
bases[[(k-1)*len2+l]]=Join[part1[[k]],-part2[[l]]];
];
];
basetot=Join[basetot,bases];
];
)
```
