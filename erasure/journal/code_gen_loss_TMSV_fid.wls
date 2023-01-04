#!/usr/bin/env wolframscript
(* ::Package:: *)

Clear["Global`*"];
Import["/home/eduardo/quantum/CV-Quantum_Error_Correction/erasure/journal/nonG/ToMatlab.m", "Package"];

(*---------------- Variables -----------------*)
T1 = .;
T2 = .;
T3 = .;
V = .;
B = .;
g = .;


(*---------------- TMSV -----------------*)
CFtmsv[a_,ac_,b_,bc_]=Exp[-1/2*(V*(a*ac+b*bc)- Sqrt[V^2-1]*(a*b+ac*bc))];
CFtmsv2[a_,ac_,b_,bc_]=Exp[-1/2*(B*(a*ac+b*bc)- Sqrt[B^2-1]*(a*b+ac*bc))];

(*---------------- Vacuum  state -----------------*)
CF0[a_,ac_] = Exp[-(a*ac)/2];


(*------- Output CF of code --------------*)
CFcode[a1_, b1_, a_, ac_] = CFtmsv2[a1, b1, a*((T1+T2)/2 + g*(T1-T2)/2), ac*((T1+T2)/2 + g*(T1-T2)/2)]*CFtmsv[a*((T1-T2)/2 + g*(T1+T2)/2), ac*((T1-T2)/2 + g*(T1+T2)/2), ac*T3*g, a*T3*g]*
				 CF0[a*Sqrt[(1+eps)*(2*g^2+1) - (T1^2)/2*(1+g)^2 - (T2^2)/2*(1-g)^2 - (T3^2)*g^2], ac*Sqrt[(1+eps)*(2*g^2+1) - (T1^2)/2*(1+g)^2 - (T2^2)/2*(1-g)^2 - (T3^2)*g^2]];								
				


(*------- Fidelity --------------*)
F[a_, ac_, b_, bc_] = CFcode[a, ac, b, bc]*CFtmsv2[-a, -ac, -b, -bc]/Pi;
FidPolar[r1_, r2_, theta1_, theta2_] = Simplify[r1*r2*F[r1*Exp[theta1*I], r1*Exp[-theta1*I], r2*Exp[theta2*I], r2*Exp[-theta2*I]]];
Print["-----------------  Fid"];

ToMatlab[FidPolar[r1, r2, theta1, theta2]]


*------- Direct transmission --------------*)
CFcode[a1_, b1_, a_, ac_] = CFtmsv2[a1, b1, a*((T1+T2)/2 + g*(T1-T2)/2), ac*((T1+T2)/2 + g*(T1-T2)/2)]*CFtmsv[a*((T1-T2)/2 + g*(T1+T2)/2), ac*((T1-T2)/2 + g*(T1+T2)/2), ac*T3*g, a*T3*g]*
				 CF0[a*Sqrt[(1+eps)*(2*g^2+1) - (T1^2)/2*(1+g)^2 - (T2^2)/2*(1-g)^2 - (T3^2)*g^2], ac*Sqrt[(1+eps)*(2*g^2+1) - (T1^2)/2*(1+g)^2 - (T2^2)/2*(1-g)^2 - (T3^2)*g^2]];								
				
