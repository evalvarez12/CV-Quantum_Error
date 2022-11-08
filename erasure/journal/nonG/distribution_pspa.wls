#!/usr/bin/env wolframscript
(* ::Package:: *)

Clear["Global`*"];
Import["/home/eduardo/quantum/CV-Quantum_Error_Correction/erasure/journal/nonG/ToMatlab.m", "Package"];

(*---------------- Variables -----------------*)
T1 = .;
T2 = .;
T3 = .;
rs = .;
d = .;
g = .;
ra = .;
ri = .;
da = .;
di = .;
xbar = .;
pbar = .;

(*---------------- Coherent state -----------------*)
CFcoh[a_,ac_] = Exp[-(a*ac)/2+(a*(ra - I*ri)-ac*(ra + I*ri))];

(*---------------- Vacuum  state -----------------*)
CF0[a_,ac_] = Exp[-(a*ac)/2];

(*---------------- TMSV -----------------*)
CFtmsv[a_,ac_,b_,bc_]=Exp[-1/2*(V*(a*ac+b*bc)- Sqrt[V^2-1]*(a*b+ac*bc))];



d1 = ((1+T)*V+R)/(1+T+R*V);
d2 = ((1+T)V+R)/(1+T+R*V);
d3 = (2*Sqrt[T*(V^2-1)])/(1+T+R*V);
CFa[a_,ac_,b_,bc_]=Exp[-(d1/2)a*ac-d2/2 b*bc+d3/2 (a*b+ac*bc)];

Tp=(T-1)/T;
Fsa0[a_,ac_,b_,bc_]=\
Tp^2*Exp[(b*bc)/2] * D[Exp[-(b*bc)] * D[Exp[(b*bc)/2] * CFa[a,ac,b,bc],b,bc],b,bc];

Ppspa=Fsa0[0,0,0,0];

CFpspa[a_,ac_,b_,bc_]=Fsa0[a,ac,b,bc]/Ppspa;



(*------- PROBABILITY DISTRIBUTION --------------*)
CFs[a_, ac_] = CFcoh[a, ac];
CFab[a_, ac_, b_, bc_] = Simplify[CFpspa[a, ac, b, bc]];	

(*------- Displacement CF --------------*)
DCF[a_,ac_] = Exp[(da - I*di)*a - (da + I*di)*ac];

(*------- Distribtion --------------*)
P[a_, ac_] = CFs[a*((T1-T2)/(Sqrt[2]*2)),ac*((T1-T2)/(Sqrt[2]*2))]*CFab[a*((T1+T2)/(Sqrt[2]*2)), ac*((T1+T2)/(Sqrt[2]*2)), ac*T3/Sqrt[2], a*T3/Sqrt[2]]*
					CF0[a*Sqrt[1-T1^2]/2, ac*Sqrt[1-T1^2]/2]*CF0[a*Sqrt[1-T2^2]/2, ac*Sqrt[1-T2^2]/2]*CF0[a*Sqrt[(1-T3^2)/2], ac*Sqrt[(1-T3^2)/2]]*
					DCF[a*((T1+T2)/(Sqrt[2]*2)), ac*((T1+T2)/(Sqrt[2]*2))]*Exp[(xbar - I*pbar)*a - (xbar + I*pbar)*ac];

PPolar[r_, theta_] = Simplify[P[r*Exp[theta*I], r*Exp[-theta*I]]];

ToMatlab[PPolar[r, theta]]



(*------- PROBABILITY DISTRIBUTION with modes inverted --------------*)
CFs[a_, ac_] = CFcoh[a, ac];
CFab[a_, ac_, b_, bc_] = Simplify[CFpspa[b, bc, a, ac]];	

(*------- Displacement CF --------------*)
DCF[a_,ac_] = Exp[(da - I*di)*a - (da + I*di)*ac];

(*------- Distribtion --------------*)
P[a_, ac_] = CFs[a*((T1-T2)/(Sqrt[2]*2)),ac*((T1-T2)/(Sqrt[2]*2))]*CFab[a*((T1+T2)/(Sqrt[2]*2)), ac*((T1+T2)/(Sqrt[2]*2)), ac*T3/Sqrt[2], a*T3/Sqrt[2]]*
					CF0[a*Sqrt[1-T1^2]/2, ac*Sqrt[1-T1^2]/2]*CF0[a*Sqrt[1-T2^2]/2, ac*Sqrt[1-T2^2]/2]*CF0[a*Sqrt[(1-T3^2)/2], ac*Sqrt[(1-T3^2)/2]]*
					DCF[a*((T1+T2)/(Sqrt[2]*2)), ac*((T1+T2)/(Sqrt[2]*2))]*Exp[(xbar - I*pbar)*a - (xbar + I*pbar)*ac];

PPolar[r_, theta_] = Simplify[P[r*Exp[theta*I], r*Exp[-theta*I]]];

ToMatlab[PPolar[r, theta]]

