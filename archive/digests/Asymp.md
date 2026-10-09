# Asymp.nb

_cells: {'Input': 4, 'Output': 7}_

```mathematica
v=1/y
Series[Sqrt[v](v+2)Sqrt[v+4]((E)^(I m1 Log[(v+4)/(v+2)])(E)^(I m2 Log[(v+2)/(v)])+(E)^( I m2 Log[(v+4)/(v+2)])(E)^( I m1 Log[(v+2)/(v)])+(E)^(I m1 Log[(v+4)/(v+2)])+(E)^(I m2 Log[(v+4)/(v+2)])+(E)^(I m2 Log[(v+2)/(v)])+(E)^(I m1 Log[(v+2)/(v)])),{y,0,0}]
```
(* Out: (1)/(y) *)
(* Out: (6)/((y)^(2))+(8 I (-3 I+m1+m2))/(y)+(12+16 I m1-8 (m1)^(2)+16 I m2-8 m1 m2-8 (m2)^(2))+Sqrt[O[y]] *)
```mathematica
v=1/y
Series[Sqrt[v](v)Sqrt[v]((E)^(I m1 Log[(v-4)/(v-2)])(E)^(I m2 Log[(v-2)/(v)])+(E)^( I m2 Log[(v-4)/(v-2)])(E)^( I m1 Log[(v-2)/(v)])+(E)^(I m1 Log[(v-4)/(v-2)])+(E)^(I m2 Log[(v-4)/(v-2)])+(E)^(I m2 Log[(v-2)/(v)])+(E)^(I m1 Log[(v-2)/(v)])),{y,0,0}]
```
(* Out: (1)/(y) *)
(* Out: (6)/((y)^(2))-(8 I (m1+m2))/(y)-8 (2 I m1+(m1)^(2)+2 I m2+m1 m2+(m2)^(2))+(O[y])^(1) *)
```mathematica
v=1/y
Series[Sqrt[v](v+2)Sqrt[v]((E)^(I m1 Log[(v)/(v+2)])(E)^(I m2 Log[(v+2)/(v)])+(E)^( I m2 Log[(v)/(v+2)])(E)^( I m1 Log[(v+2)/(v)])+(E)^(I m1 Log[(v)/(v+2)])+(E)^(I m2 Log[(v)/(v+2)])+(E)^(I m2 Log[(v+2)/(v)])+(E)^(I m1 Log[(v+2)/(v)])),{y,0,0}]
```
(* Out: (1)/(y) *)
(* Out: (6)/((y)^(2))+(12)/(y)-8 ((m1)^(2)-m1 m2+(m2)^(2))+(O[y])^(1) *)
```mathematica
Series[Sqrt[v](v)Sqrt[v]((E)^(I m1 Log[(v)/(v-2)])(E)^(I m2 Log[(v-2)/(v)])+(E)^( I m2 Log[(v)/(v-2)])(E)^( I m1 Log[(v-2)/(v)])+(E)^(I m1 Log[(v)/(v-2)])+(E)^(I m2 Log[(v)/(v-2)])+(E)^(I m2 Log[(v-2)/(v)])+(E)^(I m1 Log[(v-2)/(v)])),{y,0,0}]
```
(* Out: (6)/((y)^(2))-8 ((m1)^(2)-m1 m2+(m2)^(2))+(O[y])^(1) *)
