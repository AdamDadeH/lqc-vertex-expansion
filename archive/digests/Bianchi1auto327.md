# Bianchi1auto327.nb

_cells: {'Title': 1, 'Section': 3, 'Input': 6, 'Output': 1}_


# Automatic Generator of Amplitudes in the vacuum model.


## Definition of the matrix elements of the Bianchi I constraint operator

```mathematica
The definition here corresponds to Theta_{v2,v1}.  v1 is the second element.
```
```mathematica
OffD[v1_,v2_,m1_,m2_]:=Piecewise[{{-Sqrt[v1 v2]((v1+v2))/(2)((E)^(I m1 Log[ (2 v1)/(v1+v2)])(E)^(I m2 Log[ (v1+v2)/(2 v2)])+(E)^(I m1 Log[ (2 v1)/(v1+v2)])+(E)^(I m1 Log[ (v1+v2)/(2 v2)])+(E)^(I m2 Log[ (2 v1)/(v1+v2)])(E)^(I m1 Log[ (v1+v2)/(2 v2)])+(E)^(I m2 Log[ (2 v1)/(v1+v2)])+(E)^(I m2 Log[ (v1+v2)/(2 v2)])),v1!=0&&v2!=0},{0,v1==0},{0,v2==0}},0];
Diag[v1_,m1_,m2_]:=Piecewise[{{v1 (v1+2)((E)^(I m1 Log[ (v1)/(v1+2)])(E)^(I m2 Log[ (v1+2)/(v1)])+(E)^(I m1 Log[ (v1)/(v1+2)])+(E)^(I m1 Log[ (v1+2)/(v1)])+(E)^(I m2 Log[ (v1)/(v1+2)])(E)^(I m1 Log[ (v1+2)/(v1)])+(E)^(I m2 Log[ (v1)/(v1+2)])+(E)^(I m2 Log[ (v1+2)/(v1)]))+v1 (v1-2)((E)^(I m1 Log[ (v1)/(v1-2)])(E)^(I m2 Log[ (v1-2)/(v1)])+(E)^(I m1 Log[ (v1)/(v1-2)])+(E)^(I m1 Log[ (v1-2)/(v1)])+(E)^(I m2 Log[ (v1)/(v1-2)])(E)^(I m1 Log[ (v1-2)/(v1)])+(E)^(I m2 Log[ (v1)/(v1-2)])+(E)^(I m2 Log[ (v1-2)/(v1)])),v1!=0},{1,v1==0}},1];
```

## Definition of the parameters of the expansion

```mathematica
m1=5;
m2=5;
δ=0.01;
```

## Functions Used In Evaluating Amplitudes

```mathematica
Amp generates the amplitude corresponding to a generic discrete history.
Changes and paths together generate the set of all possible paths between the initial and final volume with m transitions.
ManyAmp them takes that set of paths and computes the overall amplitude.
```
```mathematica
Nume[v_]:=(I ((-1))^(Length[v]-1))/(2π)Product[OffD[v[[i]],v[[i+1]],m1,m2],{i,1,Length[v]-1}];
Denompos[v_]:=Product[(Diag[v[[i]],m1,m2]+I δ),{i,1,Length[v]}];
Denomneg[v_]:=Product[(Diag[v[[i]],m1,m2]-I δ),{i,1,Length[v]}];
Amp[v_]:=Nume[v]((1)/(Denompos[v])-(1)/(Denomneg[v]));
Changes[m_,vi_,vf_]:=(
diff=(vf-vi)/4;
base=Join[Table[4,{i,1,m/2+diff/2}],Table[-4,{i,1,m/2-diff/2}]];
changes=Permutations[base];
numpaths=Length[changes];)
Paths[m_,vi_,vf_]:=(
paths=Table[1,{j,1,numpaths}];
path=Table[vi,{i,1,m+1}];
For[j=1,j<=numpaths,j++,
For[i=1,i<=m,i++,
path[[i+1]]=path[[i]]+changes[[j]][[i]];];
paths[[j]]=path;
];);
NullPath[v_]:=Product[v[[i]],{i,1,Length[v]}];
ManyAmp[m_,vi_,vf_]:=(
Changes[m,vi,vf];
Paths[m,vi,vf];
Many=0;
For[i=1,i<=Length[paths],i++,
If[NullPath[paths[[i]]]==0,Many=Many,Many=Many+Amp[paths[[i]]]];
];)
```
```mathematica
ManyAmp[14,4,4]
Many
```
(* Out: -18.55597294432733`-4.0689106531141147`*^-16 I *)
