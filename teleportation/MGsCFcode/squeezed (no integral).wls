#!/usr/bin/env wolframscript
(* ::Package:: *)

Clear["Global`*"];
(*Import["\\Users\\z5239621\\Documents\\Repos\\CV-Quantum_Error_Correction\\teleportation\\MGsCFcode\\ToMatlab.m", "Package"];*)
Import["ToMatlab.m", "Package"];

t=.;
g=.;
eta=.;
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

Fsq[a_,ac_]=Exp[-((a*ac)*(1 + 2*Sinh[ra]^2) + (Sinh[2*ra]/2)*(a*a + ac*ac))/2];
Feprsq[a_,ac_]=Exp[-g^2*(1-eta^2)*((a*ac)/2)]*Feprt[-a,-ac,-g eta ac,-g eta a]*Fsq[-g eta a,-g eta ac]*Fsq[a,ac]/Pi;
FeprsqPolar[r_,theta_]=Simplify[r*Feprsq[r*Exp[theta*I],r*Exp[-theta*I]]];
Print["TMSV"]
ToMatlab[FeprsqPolar[r,theta]]

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

Fssq[a_,ac_]=Exp[-g^2*(1-eta^2)*((a*ac)/2)]*Fs[-a,-ac,-g eta ac,-g eta a]*Fsq[-g eta a,-g eta ac]*Fsq[a,ac]/Pi;

FssqPolar[r_,theta_]=Simplify[r*Fssq[r*Exp[theta*I],r*Exp[-theta*I]]];
Print["PS"]
ToMatlab[FssqPolar[r,theta]]

(*-----------------------PA-----------------------*)
Faa[a_,ac_,b_,bc_]=\
-(1-T)Exp[(b*bc)/2] * D[Exp[-(b*bc)/2] * Fzpc[a,ac,b,bc],b,bc];

Pa=Simplify[Faa[0,0,0,0]];

Fa[a_,ac_,b_,bc_]=Faa[a,ac,b,bc]/Pa;

Fasq[a_,ac_]=Exp[-g^2*(1-eta^2)*((a*ac)/2)]*Fa[-a,-ac,-g eta ac,-g eta a]*Fsq[-g eta a,-g eta ac]*Fsq[a,ac]/Pi;

FasqPolar[r_,theta_]=Simplify[r*Fasq[r*Exp[theta*I],r*Exp[-theta*I]]];
Print["PA"]
ToMatlab[FasqPolar[r,theta]]

(*-----------------------PC-----------------------*)
Tp=(T-1)/T;
Fca[a_,ac_,b_,bc_]=\
Tp^2*Exp[(b*bc)/2] * D[Exp[-(b*bc)] * D[Exp[(b*bc)/2] * Fzpc[a,ac,b,bc],b,bc],b,bc]+\
-Tp*Exp[(b*bc)/2] * D[Exp[-(b*bc)] * D[Exp[(b*bc)/2] * Fzpc[a,ac,b,bc],b],bc]+\
-Tp*Exp[(b*bc)/2] * D[Exp[-(b*bc)] * D[Exp[(b*bc)/2] * Fzpc[a,ac,b,bc],bc],b]+\
Fzpc[a,ac,b,bc];

Pc=T*Fca[0,0,0,0];

Fc[a_,ac_,b_,bc_]=T*Fca[a,ac,b,bc]/Pc;

Fcsq[a_,ac_]=Exp[-g^2*(1-eta^2)*((a*ac)/2)]*Fc[-a,-ac,-g eta ac,-g eta a]*Fsq[-g eta a,-g eta ac]*Fsq[a,ac]/Pi;

FcsqPolar[r_,theta_]=Simplify[r*Fcsq[r*Exp[theta*I],r*Exp[-theta*I]]];
Print["PC"]
ToMatlab[FcsqPolar[r,theta]]


Clear["Global`*"];
(*Import["\\Users\\z5239621\\Documents\\Repos\\CV-Quantum_Error_Correction\\teleportation\\MGsCFcode\\ToMatlab.m", "Package"];*)
Import["ToMatlab.m", "Package"];

g=.;
eta=.;
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

Fsq[a_,ac_]=Exp[-Norm[a*Cosh[ra] + ac*Sinh[ra]]^2/2];
Fsasq[a_,ac_]=Exp[-g^2*(1-eta^2)*((a*ac)/2)]*Fsa[-a,-ac,-g eta ac,-g eta a]*Fsq[-g eta a,-g eta ac]*Fsq[a,ac]/Pi;

FsasqPolar[r_,theta_]=Simplify[r*Fsasq[r*Exp[theta*I],r*Exp[-theta*I]]];
Print["PS-PA"]
ToMatlab[FsasqPolar[r,theta]]

(*-----------------------PA-PS (aa+)-----------------------*)

Tp=(T-1);
Fas0[a_,ac_,b_,bc_]=\
Tp^2*Exp[-(b*bc)/2] * D[Exp[(b*bc)] * D[Exp[-(b*bc)/2] * Fzpc[a,ac,b,bc],b,bc],b,bc];

Pas=Fas0[0,0,0,0];

Fas[a_,ac_,b_,bc_]=Fas0[a,ac,b,bc]/Pas;

Fsq[a_,ac_]=Exp[-Norm[a*Cosh[ra] + ac*Sinh[ra]]^2/2];
Fassq[a_,ac_]=Exp[-g^2*(1-eta^2)*((a*ac)/2)]*Fas[-a,-ac,-g eta ac,-g eta a]*Fsq[-g eta a,-g eta ac]*Fsq[a,ac]/Pi;

FassqPolar[r_,theta_]=Simplify[r*Fassq[r*Exp[theta*I],r*Exp[-theta*I]]];
Print["PA-PS"]
ToMatlab[FassqPolar[r,theta]]















