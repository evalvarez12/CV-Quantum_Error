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


(*-----------------------SB-----------------------*)		
(*CFsb[xa_, pa_, xb_, pb_] = 1/Sqrt[Cos[d]^2 + Sin[d]^2]*Exp[-V/2*(xa^2 + pa^2 + xb^2 + pb^2) + V2*(xa*xb - pa*pb)](Cos[d]^2 - 2*Cos[d]*Sin[d]*(V*(xa*xb-pa*pb) - V2/2*(xa^2 + pa^2 + xb^2 + pb^2)) +
							Sin[d]^2*(1 - (xa^2 + pa^2 + xb^2 + pb^2)*V + 2*V2*(xa*xb - pa*pb)+ (V^2 + V2^2/2)*((xa*xb)^2 + (pa*pb)^2) + (V^2 - V2^2/2)*((xa*pb)^2 + (pa*xb)^2) + V2^2/4*(xa^4 + xb^4 + pa^4 + pb^4) +
							V*V2*(-xa^3*xb + xa^2*pa*pb - xa*xb^3 + xb^2*pa*pb - xa*xb*pa^2 - xa*xb*pb^2 + pa*pb^3 + pa^3*pb) + V2^2/2*(xa^2*pa^2 + xb^2*pb^2) - 2*V2^2*xa*xb*pa*pb));		
*)
ap[a_, ac_, b_ ,bc_] = Cosh[rs]*a - Sinh[rs]*bc;
bp[a_, ac_, b_ ,bc_] = Cosh[rs]*b - Sinh[rs]*ac;
acp[a_, ac_, b_ ,bc_] = Cosh[rs]*ac - Sinh[rs]*b;
bcp[a_, ac_, b_ ,bc_] = Cosh[rs]*bc - Sinh[rs]*a;

CFsb[a_, ac_, b_, bc_] = 1/Sqrt[Cos[d]^2 + Sin[d]^2]*Exp[-1/2*(ap[a, ac, b ,bc]*acp[a, ac, b ,bc] + bp[a, ac, b ,bc]*bcp[a, ac, b ,bc])]*
						(Cos[d]^2 + Cos[d]*Sin[d]*(ap[a, ac, b ,bc]*bp[a, ac, b ,bc] + acp[a, ac, b ,bc]*bcp[a, ac, b ,bc]) + 
						Sin[d]^2*(1-ap[a, ac, b ,bc]*acp[a, ac, b ,bc])*(1-bp[a, ac, b ,bc]*bcp[a, ac, b ,bc]));

			
(*CFsb[a_, ac_, b_, bc_] = 1/Sqrt[Cos[d]^2 + Sin[d]^2]*Exp[-V/2*(a*ac + b*bc) +V2*(a*b + ac*bc)](Cos[d]^2 - 2*Cos[d]*Sin[d]*(V*(a*b+ac*bc) -V2/2*(a*ac + b*bc)) +
						Sin[d]^2*(1 - V*(a*ac + b*bc) + 2*V2*(a*b + ac*bc) + (V^2 + V2^2/2)*(()*)
						

						
																		
(*------- Output CF of code --------------*)
CFcode[a_,ac_] = CFs[a*((T1+T2)/2 + g*(T1-T2)/2), ac*((T1+T2)/2 + g*(T1-T2)/2)]*CFab[a*((T1-T2)/2 + g*(T1+T2)/2), ac*((T1-T2)/2 + g*(T1+T2)/2), ac*T3*g, a*T3*g]*
				 CF0[a*Sqrt[(1-T1^2)/2*(1+g)^2 + (1-T2^2)/2*(1-g)^2 + (1-T3^2)*g^2], ac*Sqrt[(1-T1^2)/2*(1+g)^2 + (1-T2^2)/2*(1-g)^2 + (1-T3^2)*g^2]];


(*------- PROBABILITY DISTRIBUTION --------------*)
CFs[a_, ac_] = CFcoh[a, ac];
CFab[a_, ac_, b_, bc_] = Simplify[CFsb[a, ac, b, bc]];	

(*------- Displacement CF --------------*)
DCF[a_,ac_] = Exp[(da - I*di)*a - (da + I*di)*ac];

(*------- Distribtion --------------*)
P[a_, ac_] = CFs[a*((T1-T2)/(Sqrt[2]*2)),ac*((T1-T2)/(Sqrt[2]*2))]*CFab[a*((T1+T2)/(Sqrt[2]*2)), ac*((T1+T2)/(Sqrt[2]*2)), ac*T3/Sqrt[2], a*T3/Sqrt[2]]*
					CF0[a*Sqrt[1-T1^2]/2, ac*Sqrt[1-T1^2]/2]*CF0[a*Sqrt[1-T2^2]/2, ac*Sqrt[1-T2^2]/2]*CF0[a*Sqrt[(1-T3^2)/2], ac*Sqrt[(1-T3^2)/2]]*
					DCF[a*((T1+T2)/(Sqrt[2]*2)), ac*((T1+T2)/(Sqrt[2]*2))]*Exp[(xbar - I*pbar)*a - (xbar + I*pbar)*ac];

PPolar[r_, theta_] = Simplify[P[r*Exp[theta*I], r*Exp[-theta*I]]];

ToMatlab[PPolar[r, theta]]

