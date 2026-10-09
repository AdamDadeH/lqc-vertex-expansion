# GaussExpansion2.nb

_cells: {'Input': 9, 'Output': 11}_

```mathematica
Sum[(((-a))^(n))/(n!)(2n)! (1)/(((-I b + δ))^(2 n+1)),{n,0,Infinity},Regularization->"Borel"]
```
(* Out: ((E)^(-(((b+I δ))^(2))/(4 a)) (-I b+δ) Gamma[(1)/(2),-(((b+I δ))^(2))/(4 a)])/(2 a Sqrt[-(((b+I δ))^(2))/(a)]) *)
```mathematica
Series[((E)^(-(((b+I δ))^(2))/(4 a)) (-I b+δ) Gamma[(1)/(2),-(((b+I δ))^(2))/(4 a)])/(2 a Sqrt[-(((b+I δ))^(2))/(a)]),{δ,0,2}]
```
(* Out: (I Sqrt[-((b)^(2))/(a)] (E)^(-((b)^(2))/(4 a)) Gamma[(1)/(2),-((b)^(2))/(4 a)])/(2 b)-(((E)^(-((b)^(2))/(4 a)) (2 (E)^(((b)^(2))/(4 a))-Sqrt[-((b)^(2))/(a)] Gamma[(1)/(2),-((b)^(2))/(4 a)])) δ)/(4 a)+(I (E)^(-((b)^(2))/(4 a)) (2 (b)^(2) (E)^(((b)^(2))/(4 a))+2 a Sqrt[-((b)^(2))/(a)] Gamma[(1)/(2),-( ...[412 chars] *)
```mathematica
Subscript[∫, 0]^(Infinity)(E)^(-a (x)^(2))(E)^(I b x)dx
```
(* Out: If[Re[a]>0,((E)^(-((b)^(2))/(4 a)) Sqrt[π] (1+I Erfi[(b)/(2 Sqrt[a])]))/(2 Sqrt[a]),Integrate[(E)^(x (I b-a x)),{x,0,Infinity},Assumptions->Re[a]<=0]] *)
```mathematica
Subscript[∫, 0]^(Infinity)(E)^(-a (x)^(2))(E)^(I b x-δ x)dx
```
(* Out: If[Re[a]>0,((E)^(-(((b+I δ))^(2))/(4 a)) Sqrt[π] Erfc[(-I b+δ)/(2 Sqrt[a])])/(2 Sqrt[a]),Integrate[(E)^(-x (-I b+a x+δ)),{x,0,Infinity},Assumptions->Re[a]<=0]] *)
```mathematica
FunctionExpand[Erfi[-I x]]
FunctionExpand[Gamma[1/2,x]]
FunctionExpand[Erfc[a]]
```
(* Out: -I Erf[x] *)
(* Out: Sqrt[π] (1-Erf[Sqrt[x]]) *)
(* Out: 1-Erf[a] *)
```mathematica
Sum[(((a))^(n))/(n!)(2n)! (1)/(((I b + δ))^(2 n+1)),{n,0,Infinity},Regularization->"Borel"]
```
(* Out: -(I (E)^((((b-I δ))^(2))/(4 a)) Sqrt[(((b-I δ))^(2))/(a)] Gamma[(1)/(2),(((b-I δ))^(2))/(4 a)])/(2 (b-I δ)) *)
```mathematica
Series[-(I (E)^((((b-I δ))^(2))/(4 a)) Sqrt[(((b-I δ))^(2))/(a)] Gamma[(1)/(2),(((b-I δ))^(2))/(4 a)])/(2 (b-I δ)),{δ,0,1}]
```
(* Out: -(I Sqrt[((b)^(2))/(a)] (E)^(((b)^(2))/(4 a)) Gamma[(1)/(2),((b)^(2))/(4 a)])/(2 b)-((-2+Sqrt[((b)^(2))/(a)] (E)^(((b)^(2))/(4 a)) Gamma[(1)/(2),((b)^(2))/(4 a)]) δ)/(4 a)+(O[δ])^(2) *)
```mathematica
FullSimplify[FunctionExpand[(I (E)^((((b+I δ))^(2))/(4 a)) Sqrt[(((b+I δ))^(2))/(a)] Gamma[(1)/(2),(((b+I δ))^(2))/(4 a)])/(2 (b+I δ))+-(I (E)^((((b-I δ))^(2))/(4 a)) Sqrt[(((b-I δ))^(2))/(a)] Gamma[(1)/(2),(((b-I δ))^(2))/(4 a)])/(2 (b-I δ))]]
```
(* Out: -(1)/(2) I (E)^((((b-I δ))^(2))/(4 a)) Sqrt[π] ((Sqrt[(((b-I δ))^(2))/(a)])/(b-I δ)-((E)^((I b δ)/(a)) Sqrt[(((b+I δ))^(2))/(a)])/(b+I δ)+(-Erf[(b-I δ)/(2 Sqrt[a])]+(E)^((I b δ)/(a)) Erf[(b+I δ)/(2 Sqrt[a])])/(Sqrt[a])) *)
```mathematica
Series[-(1)/(2) I (E)^((((b-I δ))^(2))/(4 a)) Sqrt[π] ((Sqrt[(((b-I δ))^(2))/(a)])/(b-I δ)-((E)^((I b δ)/(a)) Sqrt[(((b+I δ))^(2))/(a)])/(b+I δ)+(-Erf[(b-I δ)/(2 Sqrt[a])]+(E)^((I b δ)/(a)) Erf[(b+I δ)/(2 Sqrt[a])])/(Sqrt[a])),{δ,0,1}]
```
(* Out: ((1)/(a)-(Sqrt[((b)^(2))/(a)] (E)^(((b)^(2))/(4 a)) Sqrt[π])/(2 a)+(b (E)^(((b)^(2))/(4 a)) Sqrt[π] Erf[(b)/(2 Sqrt[a])])/(2 (a)^(3/2))) δ+(O[δ])^(2) *)
