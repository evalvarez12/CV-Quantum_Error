#!/usr/bin/env wolframscript
(* ::Package:: *)

Clear["Global`*"];
(*Import["/home/eduardo/quantum/CV-Quantum_Error_Correction/erasure/ToMatlab.m", "Package"];*)
Import["\\Users\\z5239621\\Documents\\Repos\\CV-Quantum_Error_Correction\\erasure\\ToMatlab.m", "Package"];
(*Import["/home/eduardo/QuantumComms/CV-Quantum_Error_Correction/chain/cm.wls", "Package"];*)

V=.;
gx=.;
gp=.;
ar=.;
ai=.;


CFtmsv[xa_, pa_, xb_, pb_] = Exp[-V/2*(xa^2 + pa^2 + xb^2 + pb^2) + Sqrt[V^2-1]*(xa*xb - pa*pb)];
CF0[x_, p_] = Exp[-1/2*(x^2 + p^2)];
(*CFalpha[x_, p_] = Exp[-1/2*(x^2 + p^2) + 2*I*(p*ar - x*ai)];*)
CFalpha[x_, p_] = Exp[-1/2*(x^2 + p^2) + 2*I*(p*ar)];

(*----------------------------------------*)
CFs[x_, p_] = CFalpha[x, p];
CFab[xa_, pa_, xb_, pb_] = CFtmsv[xa, pa, xb, pb];
CF0[a_, a_c] = Exp[(a*ac)/2];

(*------- Mode 1 lost --------------*)
CF1[x_, p_] = CFs[(1-gx)*x/2, (1-gp)*p/2]*CFab[-(1-gx)*x/2, -(1-gp)*p/2, gx*x, -gp*p]*CF0[(1+gx)*x/Sqrt[2], (1+gp)*p/Sqrt[2]];

(*------- Mode 2 lost --------------*)
CF2[x_, p_] = CFs[(1+gx)*x/2, (1+gp)*p/2]*CFab[(1+gx)*x/2, (1+gp)*p/2, gx*x, -gp*p]*CF0[(1-gx)*x/Sqrt[2], (1-gp)*p/Sqrt[2]];

(*------- Mode 3 lost --------------*)
CF3[x_, p_] = CFs[x, p]*CFab[gx*x, gp*p, 0, 0]*CF0[gx*x, -gp*p];

(*------- Mode 1,2 lost --------------*)
CF12[x_, p_] = CFab[0, 0, gx*x, -gp*p]*CF0[(1+gx)*x/Sqrt[2], (1+gp)*p/Sqrt[2]]*CF0[(1-gx)*x/Sqrt[2], (1-gp)*p/Sqrt[2]];

(*------- Mode 1,3 lost --------------*)
CF13[x_, p_] = CFs[(1-gx)*x/2, (1-gp)*p/2]*CFab[-(1-gx)*x/2, -(1-gp)*p/2, 0, 0]*CF0[(1+gx)*x/2, (1+gp)*p/2]*CF0[(1+gx)*x/2, (1-gp)*p/2];

(*------- Mode 2,3 lost --------------*)
CF23[x_, p_] = CFs[(1+gx)*x/2, (1+gp)*p/2]*CFab[(1+gx)*x/2, (1+gp)*p/2, 0, 0]*CF0[(1-gx)*x/Sqrt[2], (1-gp)*p/Sqrt[2]]*CF0[gx*x, -gp*p];



Fid1[x_, p_] = CF1[x, p] * CFs[-x, -p]/Pi;
Fid1Polar[r_, theta_] = Simplify[r*Fid1[r*Cos[theta], r*Sin[theta]]];
Print["-----------------  Fid1"];
ToMatlab[Fid1Polar[r, theta]]

Fid2[x_, p_] = CF2[x, p] * CFs[-x, -p]/Pi;
Fid2Polar[r_, theta_] = Simplify[r*Fid2[r*Cos[theta], r*Sin[theta]]];
Print["-----------------  Fid2"];
ToMatlab[Fid2Polar[r, theta]]

Fid3[x_, p_] = CF3[x, p] * CFs[-x, -p]/Pi;
Fid3Polar[r_, theta_] = Simplify[r*Fid3[r*Cos[theta], r*Sin[theta]]];
Print["-----------------  Fid3"];
ToMatlab[Fid3Polar[r, theta]]

Fid12[x_, p_] = CF12[x, p] * CFs[-x, -p]/Pi;
Fid12Polar[r_, theta_] = Simplify[r*Fid12[r*Cos[theta], r*Sin[theta]]];
Print["-----------------  Fid12"];
ToMatlab[Fid12Polar[r, theta]]

Fid13[x_, p_] = CF13[x, p] * CFs[-x, -p]/Pi;
Fid13Polar[r_, theta_] = Simplify[r*Fid13[r*Cos[theta], r*Sin[theta]]];
Print["-----------------  Fid13"];
ToMatlab[Fid13Polar[r, theta]]

Fid23[x_, p_] = CF23[x, p] * CFs[-x, -p]/Pi;
Fid23Polar[r_, theta_] = Simplify[r*Fid23[r*Cos[theta], r*Sin[theta]]];
Print["-----------------  Fid23"];
ToMatlab[Fid23Polar[r, theta]]

(*ToMatlab[CF0[((1+gx) x)/Sqrt[2],((1+gp) p)/Sqrt[2]] CFalpha[1/2 (1-gx) x,1/2 (1-gp) p] CFtmsv[1/2 (-1+gx) x,1/2 (-1+gp) p,gx x,-gp p]]*)


CFs[0,0]
