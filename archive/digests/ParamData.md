# ParamData.nb

_cells: {'Input': 22, 'Output': 26}_

```mathematica
Clear[x]
rmin=-5
rmax=5
```
(* Out: -5 *)
(* Out: 5 *)
```mathematica
Amp0 = ManyAmp[0,1,1]
Plot[{Re[Amp0],Im[Amp0]},{x,0,10}]
```
(* Out: (E)^(I Sqrt[2] x) *)
(* Out: <<Graphics>> *)
```mathematica
pre0 = Plot[{Re[Amp0],Im[Amp0]},{x,0,10},AxesLabel->{"\!\(\*
StyleBox[SqrtBox[
RowBox[{"12", "π", " ", "G"}]], "Text"]\)\!\(\*
StyleBox["φ", "Text"]\)","\!\(\*
StyleBox[SubscriptBox["A", "0"], "Text"]\)\!\(\*
StyleBox["(", "Text"]\)\!\(\*
StyleBox["4", "Text"]\)\!\(\*
StyleBox["λ", "Text"]\)\!\(\*
StyleBox[",", "Text"]\)\!\(\*
StyleBox[" ", "Text"]\)\!\(\*
StyleBox["4", "Text"]\)\!\(\*
StyleBox["λ", "Text"]\)\!\(\*
StyleBox[")", "Text"]\)"},PlotRange->{rmin,rmax},LabelStyle->Medium,PlotStyle->{Thickness[0.004], {Thickness[0.004],Dashed}}]
```
(* Out: <<Graphics>> *)
```mathematica
Amp2=ManyAmp[2,1,1]
Plot[{Re[Amp2],Im[Amp2]},{x,0,10}]
```
(* Out: (9)/(2) (-(1)/(36) (E)^(I Sqrt[2] x)+(1)/(36) (E)^(2 I Sqrt[2] x)-(I (E)^(I Sqrt[2] x) x)/(12 Sqrt[2])) *)
(* Out: <<Graphics>> *)
```mathematica
pre0 = Plot[{Re[Amp2],Im[Amp2]},{x,0,10},AxesLabel->{"\!\(\*
StyleBox[SqrtBox[
RowBox[{"12", "π", " ", "G"}]], "Text"]\)\!\(\*
StyleBox["φ", "Text"]\)","\!\(\*
StyleBox[SubscriptBox["A", "2"], "Text"]\)\!\(\*
StyleBox["(", "Text"]\)\!\(\*
StyleBox["4", "Text"]\)\!\(\*
StyleBox["λ", "Text"]\)\!\(\*
StyleBox[",", "Text"]\)\!\(\*
StyleBox[" ", "Text"]\)\!\(\*
StyleBox["4", "Text"]\)\!\(\*
StyleBox["λ", "Text"]\)\!\(\*
StyleBox[")", "Text"]\)"},PlotRange->{rmin,rmax},LabelStyle->Medium,PlotStyle->{Thickness[0.004], {Thickness[0.004],Dashed}}]
```
(* Out: <<Graphics>> *)
```mathematica
Amp4=ManyAmp[4,1,1];
Plot[{Re[Amp4],Im[Amp4]},{x,0,10}]
```
(* Out: <<Graphics>> *)
```mathematica
pre0 = Plot[{Re[Amp4],Im[Amp4]},{x,0,10},AxesLabel->{"\!\(\*
StyleBox[SqrtBox[
RowBox[{"12", "π", " ", "G"}]], "Text"]\)\!\(\*
StyleBox["φ", "Text"]\)","\!\(\*
StyleBox[SubscriptBox["A", "2"], "Text"]\)\!\(\*
StyleBox["(", "Text"]\)\!\(\*
StyleBox["4", "Text"]\)\!\(\*
StyleBox["λ", "Text"]\)\!\(\*
StyleBox[",", "Text"]\)\!\(\*
StyleBox[" ", "Text"]\)\!\(\*
StyleBox["4", "Text"]\)\!\(\*
StyleBox["λ", "Text"]\)\!\(\*
StyleBox[")", "Text"]\)"},PlotRange->{rmin,rmax},LabelStyle->Medium,PlotStyle->{Thickness[0.004], {Thickness[0.004],Dashed}}]
```
(* Out: <<Graphics>> *)
```mathematica
Amp6=ManyAmp[6,1,1];
Plot[{Re[Amp6],Im[Amp6]},{x,0,10}]
```
(* Out: <<Graphics>> *)
```mathematica
pre0 = Plot[{Re[Amp6],Im[Amp6]},{x,0,10},AxesLabel->{"\!\(\*
StyleBox[SqrtBox[
RowBox[{"12", "π", " ", "G"}]], "Text"]\)\!\(\*
StyleBox["φ", "Text"]\)","\!\(\*
StyleBox[SubscriptBox["A", "6"], "Text"]\)\!\(\*
StyleBox["(", "Text"]\)\!\(\*
StyleBox["4", "Text"]\)\!\(\*
StyleBox["λ", "Text"]\)\!\(\*
StyleBox[",", "Text"]\)\!\(\*
StyleBox[" ", "Text"]\)\!\(\*
StyleBox["4", "Text"]\)\!\(\*
StyleBox["λ", "Text"]\)\!\(\*
StyleBox[")", "Text"]\)"},PlotRange->{rmin,rmax},LabelStyle->Medium,PlotStyle->{Thickness[0.004], {Thickness[0.004],Dashed}}]
```
(* Out: <<Graphics>> *)
```mathematica
Amp8=ManyAmp[8,1,1];
Plot[{Re[Amp8],Im[Amp8]},{x,0,10}]
```
(* Out: <<Graphics>> *)
```mathematica
pre0 = Plot[{Re[Amp8],Im[Amp8]},{x,0,10},AxesLabel->{"\!\(\*
StyleBox[SqrtBox[
RowBox[{"12", "π", " ", "G"}]], "Text"]\)\!\(\*
StyleBox["φ", "Text"]\)","\!\(\*
StyleBox[SubscriptBox["A", "2"], "Text"]\)\!\(\*
StyleBox["(", "Text"]\)\!\(\*
StyleBox["4", "Text"]\)\!\(\*
StyleBox["λ", "Text"]\)\!\(\*
StyleBox[",", "Text"]\)\!\(\*
StyleBox[" ", "Text"]\)\!\(\*
StyleBox["4", "Text"]\)\!\(\*
StyleBox["λ", "Text"]\)\!\(\*
StyleBox[")", "Text"]\)"},PlotRange->{rmin,rmax},LabelStyle->Medium,PlotStyle->{Thickness[0.004], {Thickness[0.004],Dashed}}]
```
(* Out: <<Graphics>> *)
```mathematica
Amp10=ManyAmp[10,1,1];
Plot[{Re[Amp10],Im[Amp10]},{x,0,10}]
```
(* Out: <<Graphics>> *)
```mathematica
Amp12=ManyAmp[12,1,1];
Plot[{Re[Amp12],Im[Amp12]},{x,0,10}]
```
(* Out: <<Graphics>> *)
```mathematica
Amp14=ManyAmp[14,1,1];
Plot[{Re[Amp14],Im[Amp14]},{x,0,10}]
```
(* Out: <<Graphics>> *)
```mathematica
Amp16=ManyAmp[16,1,1];
Plot[{Re[Amp16],Im[Amp16]},{x,0,10}]
```
(* Out: <<Graphics>> *)
```mathematica
Amp18=ManyAmp[18,1,1];
Plot[{Re[Amp18],Im[Amp18]},{x,0,10}]
```
(* Out: <<Graphics>> *)
```mathematica
TotalAmp=Amp0+Amp2+Amp4+Amp6+Amp8+Amp10+Amp12+Amp14+Amp16+Amp18;
```
```mathematica
aex = Aexact[1,1];
```
```mathematica
Plot[{Re[Amp0],Re[Amp0+Amp2+Amp4],Re[Amp0+Amp2+Amp4+Amp6+Amp8+Amp10],Re[TotalAmp],Re[aex]},{x,0,5}]
Plot[{Im[Amp0],Im[Amp0+Amp2+Amp4],Im[Amp0+Amp2+Amp4+Amp6+Amp8+Amp10],Im[TotalAmp],Im[aex]},{x,0,5}]
```
(* Out: <<Graphics>> *)
(* Out: <<Graphics>> *)
```mathematica
Plot[{Re[aex],Re[Amp0],Re[Amp0+Amp2+Amp4],Re[Amp0+Amp2+Amp4+Amp6+Amp8+Amp10],Re[TotalAmp]},{x,0,5},AxesLabel->{"\!\(\*
StyleBox[SqrtBox[
RowBox[{"12", "π", " ", "G"}]], "Text"]\)\!\(\*
StyleBox["Δφ", "Text"]\)",},PlotRange->{-1.5,1.5},LabelStyle->Medium,PlotStyle->{{Thickness[0.001]}, {Thickness[0.003],Dashing[{Large}]},{Thickness[0.003],Dashing[Medium]},{Thickness[0.003],DotDashed},{Thickness[0.003],Dotted}}]
```
(* Out: <<Graphics>> *)
```mathematica
x=4.1
aaa=N[{{0,Amp0},{2,Amp2},{4,Amp4},{6,Amp6},{8,Amp8},{10,Amp10},{12,Amp12},{14,Amp14},{16,Amp16},{18,Amp18}}]
ListPlot[Re[aaa],PlotStyle->PointSize[Large],AxesLabel->{"M",}]
```
(* Out: 4.1` *)
(* Out: {{0.`,0.8847170434406316`-0.4661284726828726` I},{2.`,-0.5466725722464592`-1.0066756265225956` I},{4.`,-0.8908140104361512`-0.12953749378889112` I},{6.`,-0.65673988315602`+0.35571259184487924` I},{8.`,-0.3096610558505122`+0.5461083112505509` I},{10.`,-0.012763505851162094`+0.5162864840170958` I},{12 ...[500 chars] *)
(* Out: <<Graphics>> *)
```mathematica
ListPlot[Re[aaa]]
```
(* Out: <<Graphics>> *)
