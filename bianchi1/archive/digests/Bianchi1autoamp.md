# Bianchi1autoamp.nb

_cells: {'Input': 10, 'Output': 7, 'Print': 1}_

```mathematica
OffD[v1_,v2_,m1_,m2_]:=Piecewise[{{-Sqrt[v1 v2]((v1+v2))/(2)((E)^(I m1 Log[ (2 v1)/(v1+v2)])(E)^(I m2 Log[ (v1+v2)/(2 v2)])+(E)^(I m1 Log[ (2 v1)/(v1+v2)])+(E)^(I m1 Log[ (v1+v2)/(2 v2)])+(E)^(I m2 Log[ (2 v1)/(v1+v2)])(E)^(I m1 Log[ (v1+v2)/(2 v2)])+(E)^(I m2 Log[ (2 v1)/(v1+v2)])+(E)^(I m2 Log[ (v1+v2)/(2 v2)])),v1!=0&&v2!=0},{0,v1==0},{0,v2==0}},0]
Diag[v1_,m1_,m2_]:=Piecewise[{{v1 (v1+2)((E)^(I m1 Log[ (v1)/(v1+2)])(E)^(I m2 Log[ (v1+2)/(v1)])+(E)^(I m1 Log[ (v1)/(v1+2)])+(E)^(I m1 Log[ (v1+2)/(v1)])+(E)^(I m2 Log[ (v1)/(v1+2)])(E)^(I m1 Log[ (v1+2)/(v1)])+(E)^(I m2 Log[ (v1)/(v1+2)])+(E)^(I m2 Log[ (v1+2)/(v1)]))+v1 (v1-2)((E)^(I m1 Log[ (v1)/(v1-2)])(E)^(I m2 Log[ (v1-2)/(v1)])+(E)^(I m1 Log[ (v1)/(v1-2)])+(E)^(I m1 Log[ (v1-2)/(v1)])+(E)^(I m2 Log[ (v1)/(v1-2)])(E)^(I m1 Log[ (v1-2)/(v1)])+(E)^(I m2 Log[ (v1)/(v1-2)])+(E)^(I m2 Log[ (v1-2)/(v1)])),v1!=0},{1,v1==0}},1]
```
```mathematica
Initial Volume;
vinit =4;
Final Volume;
vfin = 20;
Difference between them;
diff = (vfin-vinit)/4;
Anisotropy Parameters;
m1=5;
m2=5;
Number of Steps;
```
```mathematica
amplitude=0;
δ=0.001
m=diff;
base=Join[Table[4,{i,1,m/2+diff/2}],Table[-4,{i,1,m/2-diff/2}]];
paths=Permutations[base];
numpaths=Length[paths];

Evaluate Amplitude for Single m;
For[ii=1,ii<=numpaths,ii++,
steps=paths[[ii]];
v=Table[vinit,{i,1,m+1}];
For[i=2,i<=m+1,i++,v[[i]]=v[[i-1]]+steps[[i-1]]];
K = Product[OffD[v[[i]],v[[i+1]],m1,m2],{i,1,m}];
Denom = Product[(Diag[v[[i]],m1,m2]+I δ),{i,1,m+1}];
Denom2= Product[(Diag[v[[i]],m1,m2]-I δ),{i,1,m+1}];
amplitude=amplitude+K/Denom - K/Denom2;
]
amplitude = (I )/(2 π)((-1))^(m) amplitude;
Print[N[amplitude]]
```
(* Out: 0.001` *)
(* Print: 4.7669860236151725`*^-8+7.468163193550472`*^-8 I *)
```mathematica
0
```
```mathematica
0.00001059945079786934`+2.3978345423244067`*^-6 I-(0.000010599089750617887`+2.399429967130521`*^-6 I)
```
(* Out: 3.6104725145260697`*^-10-1.595424806114094`*^-9 I *)
```mathematica
0.000012425762788965483`+2.810183680691274`*^-6 I-(0.000012424993464154583`+2.8135832354097267`*^-6 I)
```
(* Out: 7.693248108995367`*^-10-3.3995547184526187`*^-9 I *)
```mathematica
Psi(4) = 1
Psi(20)=-0.176845-0.277053 I
```
```mathematica
(4.250496315080809`*^-6+0.` I)/(4.7669860236151725`*^-8+7.468163193550472`*^-8 I)
```
(* Out: 25.812320476937824`-40.438679864180386` I *)
```mathematica
N[(1)/(-0.176845-0.277053 I)]
```
(* Out: -1.6369608142871825`+2.5645333737493665` I *)
```mathematica
25.8123/1.63696
40.4387/2.56453
```
(* Out: 15.768436614211709` *)
(* Out: 15.768464396985022` *)
