# AutoAmplitude.nb

_cells: {'Input': 9, 'Output': 5, 'Print': 2}_

```mathematica
3/19/2011
This program generates the exact amplitude and the partial amplitudes for the timeless framework.  I have not tested this system, so it may need to be reworked!!
```
```mathematica
Amplitude for a generic path is a product of two terms.   One that does not depend on the order of the path and one that does.
```
```mathematica
Dependent on order of path.  v is the sequence of volumes describing the path.
```
```mathematica
Independent of Order of path.  w is the set of distinct volumes.  n is the set of multiplicities
```
```mathematica
a[n_Integer]=FullSimplify[(1)/((n))(2I*k*a[n-1]+(n-2)a[n-2])];
a[0]=0;
a[1]=I*k;
```
```mathematica
ImExact[n_Integer,m_Integer]:=2Sqrt[n*m]NIntegrate[Im[(E)^(I *x* k)(
1)/(k Sinh[k*π])a[n/2]*a[m/2]],{k,0,100}];
ReExact[n_Integer,m_Integer]:=2Sqrt[n*m]NIntegrate[Re[(E)^(I *x* k)(
1)/(k Sinh[k*π])a[n/2]*a[m/2]],{k,0,100}];.10.10
```
(* Out: (.10)^(2) *)
```mathematica
m=16
base=Join[Table[4,{i,1,(m-4)/2+4}],Table[-4,{i,1,(m-4)/2}]];
paths=Permutations[base];
numpaths=Length[paths]
Print[numpaths]
```
(* Out: 16 *)
(* Out: 8008 *)
(* Print: 8008 *)
```mathematica
amplitude=0
For[ii=1,ii<=numpaths,ii++,
steps=paths[[ii]];
v=Table[20,{i,1,m+1}];
For[i=2,i<=m+1,i++,v[[i]]=v[[i-1]]+steps[[i-1]]];

K = Product[(-1)/(4)Sqrt[v[[i]] v[[i+1]]](v[[i]]+v[[i+1]]),{i,1,m}];
vsort=Sort[v];
vsplit=Split[vsort];
q=Length[vsplit];
w=Table[vsplit[[i]][[1]],{i,1,q}];
degen=Table[Length[vsplit[[i]]],{i,1,q}];
cat =If[K==0,0, Product[(1)/((degen[[i]]-1)!),{i,1,q}]
Sum[(E)^(I*Sqrt[d[i]]*x/Sqrt[8])Product[(1)/(d[i]-d[j]),{j,1,i-1}]* Product[(1)/(d[i]-d[j]),{j,i+1,q}],{i,1,q}]];
For[i=1,i<=q,i++, cat=D[cat,{d[i],degen[[i]]-1}]];
For[i=1,i<=q,i++,d[i]=(w[[i]])^(2)];
amplitude=amplitude+Expand[K*cat];
Clear[d];
]
Print[amplitude]
```
(* Out: 0 *)
(* Print: (6861518937 Sqrt[5] (E)^(I Sqrt[2] x))/(140737488355328)+(8389611129 Sqrt[5] (E)^(2 I Sqrt[2] x))/(4398046511104)+(2527409816883 Sqrt[5] (E)^(3 I Sqrt[2] x))/(281474976710656)-(646104087 Sqrt[5] (E)^( *)
```mathematica
amplitude
```
(* Out: amplitude *)
