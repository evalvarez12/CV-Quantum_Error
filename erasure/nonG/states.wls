#!/usr/bin/env wolframscript
(* ::Package:: *)

Clear["Global`*"];
Import["/home/eduardo/quantum/CV-Quantum_Error_Correction/erasure/nonG/ToMatlab.m", "Package"];

(*---------------- Variables -----------------*)
T1 = .;
T2 = .;
T3 = .;
T = .;
V = .;
R = 1 - T;
g = .;
ra = .;
ri = .;

(*---------------- TMSV -----------------*)
CFtmsv[a_,ac_,b_,bc_]=Exp[-1/2*(V*(a*ac+b*bc)- Sqrt[V^2-1]*(a*b+ac*bc))];

(*---------------- Coherent state -----------------*)
CFcoh[a_,ac_] = Exp[-(a*ac)/2+(a*(ra - I*ri)-ac*(ra + I*ri))];

(*---------------- Vacuum  state -----------------*)
CF0[a_,ac_] = Exp[-(a*ac)/2];

(*------------------Ancilla definition of f function --------------------*)
d1 = ((1+T)*V+R)/(1+T+R*V);
d2 = ((1+T)V+R)/(1+T+R*V);
d3 = (2*Sqrt[T*(V^2-1)])/(1+T+R*V);
CFa[a_,ac_,b_,bc_]=Exp[-(d1/2)a*ac-d2/2 b*bc+d3/2 (a*b+ac*bc)];

(*-----------------------PS-----------------------*)
Fsa[a_,ac_,b_,bc_] = \
-((1-T)/T)Exp[-(b*bc)/2] * D[Exp[(b*bc)/2] * CFa[a,ac,b,bc],b,bc];
Pps = Simplify[Fsa[0,0,0,0]];
CFps[a_,ac_,b_,bc_] = Fsa[a,ac,b,bc]/Pps;

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
		
				
(*------- Output CF of code --------------*)
CFcode[a_,ac_] = CFs[a*((T1+T2)/2 + g*(T1-T2)/2), ac*((T1+T2)/2 + g*(T1-T2)/2)]*CFab[a*((T1-T2)/2 + g*(T1+T2)/2), ac*((T1-T2)/2 + g*(T1+T2)/2), ac*T3*g, a*T3*g]*
				 CF0[a*Sqrt[(1-T1^2)/2*(1+g)^2 + (1-T2^2)/2*(1-g)^2 + (1-T3^2)*g^2], ac*Sqrt[(1-T1^2)/2*(1+g)^2 + (1-T2^2)/2*(1-g)^2 + (1-T3^2)*g^2]];								
				


(*------- Fidelity TMSV --------------*)
CFs[a_, ac_] = CFcoh[a, ac];
CFab[a_, ac_, b_, bc_] = CFtmsv[a, ac, b, bc];

				 
Ftmsv[a_, ac_] = CFcode[a, ac]*CFs[-a, -ac]/Pi;
FidPolarTmsv[r_, theta_] = Simplify[r*Ftmsv[r*Exp[theta*I], r*Exp[-theta*I]]];
Print["-----------------  Fid tmsv"];

ToMatlab[FidPolarTmsv[r, theta]]


(*------- Fidelity Ps --------------*)
CFs[a_, ac_] = CFcoh[a, ac];
CFab[a_, ac_, b_, bc_] = Simplify[CFps[a, ac, b, bc]];	
				 
Fps[a_, ac_] = CFcode[a, ac]*CFs[-a, -ac]/Pi;
FidPolarPs[r_, theta_] = Simplify[r*Fps[r*Exp[theta*I], r*Exp[-theta*I]]];
Print["-----------------  Fid ps"];

ToMatlab[FidPolarPs[r, theta]]


(*------- Fidelity Pa --------------*)
CFs[a_, ac_] = CFcoh[a, ac];
CFab[a_, ac_, b_, bc_] = Simplify[CFpa[a, ac, b, bc]];	
				 
Fpa[a_, ac_] = CFcode[a, ac]*CFs[-a, -ac]/Pi;
FidPolarPa[r_, theta_] = Simplify[r*Fpa[r*Exp[theta*I], r*Exp[-theta*I]]];
Print["-----------------  Fid pa"];

ToMatlab[FidPolarPa[r, theta]]


(*------- Fidelity Pc --------------*)
CFs[a_, ac_] = CFcoh[a, ac];
CFab[a_, ac_, b_, bc_] = Simplify[CFpc[a, ac, b, bc]];	
				 
Fpc[a_, ac_] = CFcode[a, ac]*CFs[-a, -ac]/Pi;
FidPolarPc[r_, theta_] = Simplify[r*Fpc[r*Exp[theta*I], r*Exp[-theta*I]]];
Print["-----------------  Fid pc"];

ToMatlab[FidPolarPc[r, theta]]


(*------- Fidelity Paps --------------*)
CFs[a_, ac_] = CFcoh[a, ac];
CFab[a_, ac_, b_, bc_] = Simplify[CFpaps[a, ac, b, bc]];	
				 
Fpaps[a_, ac_] = CFcode[a, ac]*CFs[-a, -ac]/Pi;
FidPolarPaps[r_, theta_] = Simplify[r*Fpaps[r*Exp[theta*I], r*Exp[-theta*I]]];
Print["-----------------  Fid paps"];

ToMatlab[FidPolarPaps[r, theta]]


(*------- Fidelity Pspa --------------*)
CFs[a_, ac_] = CFcoh[a, ac];
CFab[a_, ac_, b_, bc_] = Simplify[CFpspa[a, ac, b, bc]];	
				 
Fpspa[a_, ac_] = CFcode[a, ac]*CFs[-a, -ac]/Pi;
FidPolarPspa[r_, theta_] = Simplify[r*Fpspa[r*Exp[theta*I], r*Exp[-theta*I]]];
Print["-----------------  Fid pspa"];

ToMatlab[FidPolarPspa[r, theta]]
