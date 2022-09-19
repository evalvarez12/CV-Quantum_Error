#!/usr/bin/env wolframscript
(* ::Package:: *)

Clear["Global`*"];
Import["/home/eduardo/quantum/CV-Quantum_Error_Correction/erasure/ToMatlab.m", "Package"];
(*Import["\\Users\\z5239621\\Documents\\Repos\\CV-Quantum_Error_Correction\\erasure\\ToMatlab.m", "Package"];*)
(*Import["/home/eduardo/QuantumComms/CV-Quantum_Error_Correction/chain/cm.wls", "Package"];*)

V=.;
V2=Sqrt[V^2-1];
T1r=.;
T2r=.;
T3r=.;
epsilon=.;
T1p=Sqrt[epsilon + 1 - T1r]; 
T2p=Sqrt[epsilon + 1 - T2r];
T3p=Sqrt[epsilon + 1 - T3r];
T1=Sqrt[T1r];
T2=Sqrt[T2r];
T3=Sqrt[T3r];
gx=.;
gp=.;
ra=.;
ri=.;
r=.;


CFtmsv[xa_, pa_, xb_, pb_] = Exp[-V/4*(xa^2 + pa^2 + xb^2 + pb^2) + V2/2*(xa*xb - pa*pb)];
CF0[x_, p_] = Exp[-1/4*(x^2 + p^2)];
CFalpha[x_, p_] = Exp[-1/4*(x^2 + p^2) + Sqrt[2]*I*(p*ra - x*ri)];
(*CFalpha[x_, p_] = Exp[-1/2*(x^2 + p^2) + 2*I*(p*ra)];*)

CFsb[xa_, pa_, xb_, pb_] = 1/Sqrt[Cos[d]^2 + Sin[d]^2]*Exp[-V/2*(xa^2 + pa^2 + xb^2 + pb^2) + V2*(xa*xb - pa*pb)](Cos[d]^2 - 2*Cos[d]*Sin[d]*(V*(xa*xb-pa*pb) - V2/2*(xa^2 + pa^2 + xb^2 + pb^2)) +
							Sin[d]^2*(1 - (xa^2 + pa^2 + xb^2 + pb^2)*V + 2*V2*(xa*xb - pa*pb)+ (V^2 + V2^2/2)*((xa*xb)^2 + (pa*pb)^2) + (V^2 - V2^2/2)*((xa*pb)^2 + (pa*xb)^2) + V2^2/4*(xa^4 + xb^4 + pa^4 + pb^4) +
							V*V2*(-xa^3*xb + xa^2*pa*pb - xa*xb^3 + xb^2*pa*pb - xa*xb*pa^2 - xa*xb*pb^2 + pa*pb^3 + pa^3*pb) + V2^2/2*(xa^2*pa^2 + xb^2*pb^2) - 2*V2^2*xa*xb*pa*pb));

CFsbred[xa_, pa_, xb_, pb_] = Exp[-V/2*(xa^2 + pa^2 + xb^2 + pb^2) + V2*(xa*xb - pa*pb)]*(1 - (xa^2 + pa^2 + xb^2 + pb^2)*V + 2*V2*(xa*xb - pa*pb)+ (V^2 + V2^2/2)*((xa*xb)^2 + (pa*pb)^2) + (V^2 - V2^2/2)*((xa*pb)^2 + (pa*xb)^2) + V2^2/4*(xa^4 + xb^4 + pa^4 + pb^4) +
							V*V2*(-xa^3*xb + xa^2*pa*pb - xa*xb^3 + xb^2*pa*pb - xa*xb*pa^2 - xa*xb*pb^2 + pa*pb^3 + pa^3*pb) + V2^2/2*(xa^2*pa^2 + xb^2*pb^2) - 2*V2^2*xa*xb*pa*pb);

(*----------------------------------------*)
CFs[x_, p_] = CFalpha[x, p];
CFab[xa_, pa_, xb_, pb_] = CFtmsv[xa, pa, xb, pb];
CF0[a_, a_c] = Exp[(a*ac)/2];

(*------- Output CF --------------*)
CF[x_, p_] = CFs[x/2*((T1+T2) + gx*(T1-T2)), p/2*((T1+T2) + gp*(T1-T2))]*CFab[x/2*((T1-T2) + gx*(T1+T2)),p/2*((T1-T2) + gp*(T1+T2)), T3*gx*x, -T3*gp*p]*
			 CF0[(1+gx)*x*T1p/Sqrt[2],(1+gp)*p*T1p/Sqrt[2]]*CF0[(1-gx)*x*T2p/Sqrt[2],(1-gp)*p*T2p/Sqrt[2]]*CF0[x*T3p*gx,-p*T3p*gp];




CF[x,p]


Fid[x_, p_] = CF[x, p] * CFs[-x, -p]/Pi;
FidPolar[r_, theta_] = Simplify[r*Fid[r*Cos[theta], r*Sin[theta]]];
Print["-----------------  Fid"];
ToMatlab[FidPolar[r, theta]]


(*------- g = 0 --------------*)
T0pp = Sqrt[(T1p^2 + T2p^2)/2];
CFg0[x_, p_] = CFs[x*(T1+T2)/2, p*(T1+T2)/2]*CFab[x*(T1-T2)/2, p*(T1-T2)/2, 0, 0]*
			 CF0[x*T0pp,p*T0pp];

(*------- g = 1 --------------*)
T1pp = Sqrt[2*T1p^2 + T3p^2];
CFg1[x_, p_] = CFs[x*T1, p*T1]*CFab[x*T1, p*T1, x*T3, -p*T3]*
			 CF0[x*T1pp,p*T1pp];
			 
(*------- g = -1 --------------*)
T2pp = Sqrt[2*T2p^2 + T3p^2];
CFg2[x_, p_] = CFs[x*T2, p*T2]*CFab[-x*T2, -p*T2, -x*T3, p*T3]*
			 CF0[x*T2pp,p*T2pp];
			 


Fidg0[x_, p_] = CFg0[x, p] * CFs[-x, -p]/Pi;
FidPolarg0[r_, theta_] = Simplify[r*Fidg0[r*Cos[theta], r*Sin[theta]]];
Print["-----------------  Fid g = 0"];
ToMatlab[FidPolarg0[r, theta]]



Fidg1[x_, p_] = CFg1[x, p] * CFs[-x, -p]/Pi;
FidPolarg1[r_, theta_] = Simplify[r*Fidg1[r*Cos[theta], r*Sin[theta]]];
Print["-----------------  Fid g = 1"];
ToMatlab[FidPolarg1[r, theta]]



Fidg2[x_, p_] = CFg2[x, p] * CFs[-x, -p]/Pi;
FidPolarg2[r_, theta_] = Simplify[r*Fidg2[r*Cos[theta], r*Sin[theta]]];
Print["-----------------  Fid g = -1"];
ToMatlab[FidPolarg2[r, theta]]



CFdir[x_, p_] = CFs[x*T1, p*T1]*CF0[x*T1p, p*T1p];
Fiddir[x_, p_] = CFdir[x, p] * CFs[-x, -p]/Pi;
FidPolardir[r_, theta_] = Simplify[r*Fiddir[r*Cos[theta], r*Sin[theta]]];
Print["-----------------  Fid dir"];
ToMatlab[FidPolardir[r, theta]]






