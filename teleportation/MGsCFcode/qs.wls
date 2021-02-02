#!/usr/bin/env wolframscript
(* ::Package:: *)

Clear["Global`*"];
(*Import["\\Users\\z5239621\\Documents\\Repos\\CV-Quantum_Error_Correction\\teleportation\\MGsCFcode\\ToMatlab.m", "Package"]*)

t=.;
g=.;
e=.;
T=.;
R=1-T;
l=.;
ra=.
ta=.

(*----------------EPR-----------------*)
Fepr[a_,ac_,b_,bc_]=Exp[-((1+l^2)/(2(1-l^2)))(a*ac+b*bc)+l/(1-l^2)*(a*b+ac*bc)];
(*Feprt[a_,ac_,b_,bc_]=Exp[-(1+e)(1-t)(a*ac)/2]Fepr[a*Sqrt[t],ac*Sqrt[t],b,bc];*)
Feprt[a_,ac_,b_,bc_]=Exp[-(1+e)(1-t)(b*bc)/2]Fepr[a,ac,b*Sqrt[t],bc*Sqrt[t]];

Fcoh[a_,ac_]=Exp[-((a*ac)/2)]Exp[a*ra-ac*ra];
Feprcoh[a_,ac_]=Feprt[-a,-ac,-g ac,-g a]*Fcoh[-g a,-g ac]*Fcoh[a,ac]/Pi;
FeprcohPolar[r_,theta_]=Simplify[r*Feprcoh[r*Exp[theta*I],r*Exp[-theta*I]]];
ToMatlab[FeprcohPolar[r,theta]]

(*Integrate[FeprcohPolar[r,theta],{theta,0,2*Pi},GenerateConditions->False]*)

b1=(1+l^2)/(1-l^2);
b2=t(b1+(1-t)/t+e);
b3=Sqrt[t]((2l)/(1-l^2));

d1=((1+T)b1+R(b1*b2-b3^2))/(1+T+R*b2);
d2=((1+T)b2+R)/(1+T+R*b2);
d3=(2Sqrt[T]*b3)/(1+T+R*b2);

Fzpc[a_,ac_,b_,bc_]=Exp[-(d1/2)a*ac-d2/2 b*bc+d3/2 (a*b+ac*bc)];

(*-----------------------PS-----------------------*)
Fsa[a_,ac_,b_,bc_]=\
-((1-T)/T)Exp[-(b*bc)/2] * D[Exp[(b*bc)/2] * Fzpc[a,ac,b,bc],b,bc];

Ps=Simplify[Fsa[0,0,0,0]];

Fs[a_,ac_,b_,bc_]=Fsa[a,ac,b,bc]/Ps;

Fscoh[a_,ac_]=Fs[-a,-ac,-g ac,-g a]*Fcoh[-g a,-g ac]*Fcoh[a,ac]/Pi;

FscohPolar[r_,theta_]=Simplify[r*Fscoh[r*Exp[theta*I],r*Exp[-theta*I]]];
ToMatlab[FscohPolar[r,theta]]

