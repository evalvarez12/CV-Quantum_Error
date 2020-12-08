#!/usr/bin/env wolframscript
(* ::Package:: *)

Clear["Global`*"];
Import["\\Users\\z5239621\\Documents\\Repos\\MGsCFcode\\ToMatlab.m", "Package"];

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

(*-----------------------PA-----------------------*)
Faa[a_,ac_,b_,bc_]=\
-(1-T)Exp[(b*bc)/2] * D[Exp[-(b*bc)/2] * Fzpc[a,ac,b,bc],b,bc];

Pa=Simplify[Faa[0,0,0,0]];

Fa[a_,ac_,b_,bc_]=Faa[a,ac,b,bc]/Pa;

Facoh[a_,ac_]=Fa[-a,-ac,-g ac,-g a]*Fcoh[-g a,-g ac]*Fcoh[a,ac]/Pi;

FacohPolar[r_,theta_]=Simplify[r*Facoh[r*Exp[theta*I],r*Exp[-theta*I]]];
ToMatlab[FacohPolar[r,theta]]

(*-----------------------PC-----------------------*)
Tp=(T-1)/T;
Fca[a_,ac_,b_,bc_]=\
Tp^2*Exp[(b*bc)/2] * D[Exp[-(b*bc)] * D[Exp[(b*bc)/2] * Fzpc[a,ac,b,bc],b,bc],b,bc]+\
-Tp*Exp[(b*bc)/2] * D[Exp[-(b*bc)] * D[Exp[(b*bc)/2] * Fzpc[a,ac,b,bc],b],bc]+\
-Tp*Exp[(b*bc)/2] * D[Exp[-(b*bc)] * D[Exp[(b*bc)/2] * Fzpc[a,ac,b,bc],bc],b]+\
Fzpc[a,ac,b,bc];

Pc=T*Fca[0,0,0,0];

Fc[a_,ac_,b_,bc_]=T*Fca[a,ac,b,bc]/Pc;

Fccoh[a_,ac_]=Fc[-a,-ac,-g ac,-g a]*Fcoh[-g a,-g ac]*Fcoh[a,ac]/Pi;

FccohPolar[r_,theta_]=Simplify[r*Fccoh[r*Exp[theta*I],r*Exp[-theta*I]]];
ToMatlab[FccohPolar[r,theta]]


Clear["Global`*"];
Import["\\Users\\z5239621\\Documents\\Repos\\MGsCFcode\\ToMatlab.m", "Package"];

g=.;
t=.;
e=.;
T=.;
a0=.;
R=1-T^2;
l=.;
ra=.;

b1=(1+l^2)/(1-l^2);
b2=t(b1+(1-t)/t+e);
b3=Sqrt[t]((2l)/(1-l^2));

d1=((1+T^2)b1+R(b1*b2-b3^2))/(1+T^2+R*b2);
d2=((1+T^2)b2+R)/(1+T^2+R*b2);
d3=(2T*b3)/(1+T^2+R*b2);

Fzpc[a_,ac_,b_,bc_]=Exp[-(d1/2)a*ac-d2/2 b*bc+d3/2 (a*b+ac*bc)];

(*-----------------------PS-PA (a+a)-----------------------*)

Tp=(T-1)/T;
Fsa0[a_,ac_,b_,bc_]=\
Tp^2*Exp[(b*bc)/2] * D[Exp[-(b*bc)] * D[Exp[(b*bc)/2] * Fzpc[a,ac,b,bc],b,bc],b,bc];

Psa=Fsa0[0,0,0,0];

Fsa[a_,ac_,b_,bc_]=Fsa0[a,ac,b,bc]/Psa;

Fcoh[a_,ac_]=Exp[-((a*ac)/2)]Exp[a*ra-ac*ra];
Fsacoh[a_,ac_]=Fsa[-a,-ac,-g ac,-g a]*Fcoh[-g a,-g ac]*Fcoh[a,ac]/Pi;

FsacohPolar[r_,theta_]=Simplify[r*Fsacoh[r*Exp[theta*I],r*Exp[-theta*I]]];
ToMatlab[FsacohPolar[r,theta]]

(*-----------------------PA-PS (aa+)-----------------------*)

Tp=(T-1);
Fas0[a_,ac_,b_,bc_]=\
Tp^2*Exp[-(b*bc)/2] * D[Exp[(b*bc)] * D[Exp[-(b*bc)/2] * Fzpc[a,ac,b,bc],b,bc],b,bc];

Pas=Fas0[0,0,0,0];

Fas[a_,ac_,b_,bc_]=Fas0[a,ac,b,bc]/Pas;

Fcoh[a_,ac_]=Exp[-((a*ac)/2)]Exp[a*ra-ac*ra];
Fascoh[a_,ac_]=Fas[-a,-ac,-g ac,-g a]*Fcoh[-g a,-g ac]*Fcoh[a,ac]/Pi;

FascohPolar[r_,theta_]=Simplify[r*Fascoh[r*Exp[theta*I],r*Exp[-theta*I]]];
ToMatlab[FascohPolar[r,theta]]






