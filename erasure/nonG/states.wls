#!/usr/bin/env wolframscript
(* ::Package:: *)

Clear["Global`*"];
Import["/home/eduardo/quantum/CV-Quantum_Error_Correction/erasure/nonG/ToMatlab.m", "Package"];


(*----------------TMSV-----------------*)
CFtmsv[a_,ac_,b_,bc_]=Exp[-1/2*(V*(a*ac+b*bc)- Sqrt[V^2-1]*(a*b+ac*bc))];


(*------------------Ancilla definition of f function --------------------*)
d1 = ((1+T)*V+R)/(1+T+R*V);
d2 = ((1+T)V+R)/(1+T+R*V);
d3 = (2*Sqrt[T*(V^2-1)])/(1+T+R*V);
CFa[a_,ac_,b_,bc_]=Exp[-(d1/2)a*ac-d2/2 b*bc+d3/2 (a*b+ac*bc)];

(*-----------------------PS-----------------------*)
Fsa[a_,ac_,b_,bc_] = \
-((1-T)/T)Exp[-(b*bc)/2] * D[Exp[(b*bc)/2] * CFa[a,ac,b,bc],b,bc];
Pps = Simplify[Fsa[0,0,0,0]];
CFps = Fsa[a_,ac_,b_,bc_]/Pps;

(*-----------------------PA-----------------------*)
Faa[a_,ac_,b_,bc_]=\
-(1-T)Exp[(b*bc)/2] * D[Exp[-(b*bc)/2] * CFa[a,ac,b,bc],b,bc];

Ppa=Simplify[Faa[0,0,0,0]];

CFpa[a_,ac_,b_,bc_]=Faa[a,ac,b,bc]/Ppa;

(*-----------------------PC-----------------------*)
Tp=(T-1)/T;
Fca[a_,ac_,b_,bc_]=\
Tp^2*Exp[(b*bc)/2] * D[Exp[-(b*bc)] * D[Exp[(b*bc)/2] * CFa[a,ac,b,bc],b,bc],b,bc]+\
-Tp*Exp[(b*bc)/2] * D[Exp[-(b*bc)] * D[Exp[(b*bc)/2] * CFa[a,ac,b,bc],b],bc]+\
-Tp*Exp[(b*bc)/2] * D[Exp[-(b*bc)] * D[Exp[(b*bc)/2] * CFa[a,ac,b,bc],bc],b]+\
CFa[a,ac,b,bc];

Ppc=T*Fca[0,0,0,0];

CFpc[a_,ac_,b_,bc_]=T*Fca[a,ac,b,bc]/Ppc;

(*-----------------------PS-PA (a+a)-----------------------*)
Tp=(T-1)/T;
Fsa0[a_,ac_,b_,bc_]=\
Tp^2*Exp[(b*bc)/2] * D[Exp[-(b*bc)] * D[Exp[(b*bc)/2] * CFa[a,ac,b,bc],b,bc],b,bc];

Ppspa=Fsa0[0,0,0,0];

CFpspa[a_,ac_,b_,bc_]=Fsa0[a,ac,b,bc]/Ppspa;

(*-----------------------PA-PS (aa+)-----------------------*)
Tp=(T-1);
Fas0[a_,ac_,b_,bc_]=\
Tp^2*Exp[-(b*bc)/2] * D[Exp[(b*bc)] * D[Exp[-(b*bc)/2] * CFa[a,ac,b,bc],b,bc],b,bc];

Ppaps=Fas0[0,0,0,0];

CFpaps[a_,ac_,b_,bc_]=Fas0[a,ac,b,bc]/Ppaps;




