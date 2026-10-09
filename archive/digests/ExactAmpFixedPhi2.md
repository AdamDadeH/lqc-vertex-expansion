# ExactAmpFixedPhi2.nb

_cells: {'Input': 17, 'Output': 14}_

```mathematica
a[n_Integer]=FullSimplify[(1)/((n))(2I*k*a[n-1]+(n-2)a[n-2])];
a[0]=0;
a[1]=I*k;
```
```mathematica
b1=Table[a[n],{n,2,5}];
a[5]=b1[[4]];
b2=Table[a[n],{n,6,10}];
a[10]=b2[[5]];
b3=Table[a[n],{n,11,15}];
a[15]=b3[[5]];
```
```mathematica
b4=Table[a[n],{n,16,20}];
a[20]=b4[[5]];
b5=Table[a[n],{n,21,25}];
a[25]=b5[[5]];
b6=Table[a[n],{n,26,30}];
a[30]=b6[[5]];
```
```mathematica
b7=Table[a[n],{n,31,35}];
a[35]=b7[[5]];
```
```mathematica
b=Join[b1,b2,b3,b4,b5,b6,b7];
```
```mathematica
For[i=0,i<30,i++,b[[i]]=FullSimplify[b[[i]]]]
```
```mathematica
b[[5]]
```
(* Out: (1)/(45) (-23 (k)^(2)+20 (k)^(4)-2 (k)^(6)) *)
```mathematica
ImPart[n_,m_,x_]:=-(1)/(π)Integrate[Sin[(m-n)/2*k - x/2*Sqrt[Sqrt[n*m]*(n+m)/2]*Sin[k]],{k,0,π}];
RePart[n_,m_,x_]:=BesselJ[(m-n)/2,x/2*Sqrt[Sqrt[n*m](n+m)/2]];

ImExact[n_,m_,x_]:=2Sqrt[n*m]NIntegrate[Im[(E)^(I *x* k)(
1)/(k Sinh[k*π])b[[n/2]]*b[[m/2]]],{k,0,100}];
ReExact[n_,m_,x_]:=2Sqrt[n*m]NIntegrate[Re[(E)^(I *x* k)(
1)/(k Sinh[k*π])b[[n/2]]*b[[m/2]]],{k,0,100}]
```
```mathematica
x=1
```
```mathematica
n=1
a41=Table[{4*n,4*m,ReExact[4*n,4*m,1]},{m,1,10}];
n=2
a81=Table[{4*n,4*m,ReExact[4*n,4*m,1]},{m,1,10}];
n=3
a121=Table[{4*n,4*m,ReExact[4*n,4*m,1]},{m,1,10}];
n=4
a161=Table[{4*n,4*m,ReExact[4*n,4*m,1]},{m,1,10}];
n=5
a201=Table[{4*n,4*m,ReExact[4*n,4*m,1]},{m,1,10}];
n=6
a241=Table[{4*n,4*m,ReExact[4*n,4*m,1]},{m,1,10}];
n=7
a281=Table[{4*n,4*m,ReExact[4*n,4*m,1]},{m,1,10}];
n=8
a321=Table[{4*n,4*m,ReExact[4*n,4*m,1]},{m,1,10}];
n=9
a361=Table[{4*n,4*m,ReExact[4*n,4*m,1]},{m,1,10}];
n=10
a401=Table[{4*n,4*m,ReExact[4*n,4*m,1]},{m,1,10}];
```
```mathematica
n=1
a42=Table[{4*n,4*m,ReExact[4*n,4*m,1]},{m,11,15}];
n=2
a82=Table[{4*n,4*m,ReExact[4*n,4*m,1]},{m,11,15}];
n=3
a122=Table[{4*n,4*m,ReExact[4*n,4*m,1]},{m,11,15}];
n=4
a162=Table[{4*n,4*m,ReExact[4*n,4*m,1]},{m,11,15}];
n=5
a202=Table[{4*n,4*m,ReExact[4*n,4*m,1]},{m,11,15}];
n=6
a242=Table[{4*n,4*m,ReExact[4*n,4*m,1]},{m,11,15}];
n=7
a282=Table[{4*n,4*m,ReExact[4*n,4*m,1]},{m,11,15}];
n=8
a322=Table[{4*n,4*m,ReExact[4*n,4*m,1]},{m,11,15}];
n=9
a362=Table[{4*n,4*m,ReExact[4*n,4*m,1]},{m,11,15}];
n=10
a402=Table[{4*n,4*m,ReExact[4*n,4*m,1]},{m,11,15}];
```
(* Out: 1 *)
(* Out: 2 *)
(* Out: 3 *)
(* Out: 4 *)
(* Out: 5 *)
(* Out: 6 *)
(* Out: 7 *)
(* Out: 8 *)
(* Out: 9 *)
(* Out: 10 *)
```mathematica

```
```mathematica

```
```mathematica
Data=Join[a41,a81,a121,a161,a201,a241,a281,a321,a361,a401,a42,a82,a122,a162,a202,a242,a282,a322,a362,a402];
```
```mathematica
Export["out.dat",Data]
```
(* Out: "out.dat" *)
```mathematica
ListPlot3D[Join[a41,a81,a121,a161,a201,a241,a281,a321,a361,a401,a42,a82,a122,a162,a202,a242,a282,a322,a362,a402],InterpolationOrder->10]
```
(* Out: <<Graphics3D>> *)
```mathematica
a361
```
(* Out: {{36,4,-0.0005830486108207827`},{36,8,-0.02131704646100334`},{36,12,-0.16960277250402`},{36,16,-0.19089031416366756`},{36,20,0.20678025123624766`},{36,24,-0.19091408925695552`},{36,28,0.17412706066549732`},{36,32,-0.10921566420191318`},{36,36,-0.030222765972962418`},{36,40,0.15724994768591913`},{36, ...[513 chars] *)
