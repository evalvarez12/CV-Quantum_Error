#!/usr/bin/env wolframscript
(* ::Package:: *)

Clear["Global`*"];
Import["/home/eduardo/QuantumComms/CV-Quantum_Error_Correction/chain/ToMatlab.m", "Package"];
Import["/home/eduardo/QuantumComms/CV-Quantum_Error_Correction/chain/cm.wls", "Package"];

p1=.;
p2=.;
p3=.;
eta=.;
g1=.;
g2=.;
T=.;
R = 1-T;
(*----------------Initial Gaussian CF -----------------*)
CFgauss[a_, ac_, b_, bc_] = Exp[p1*(a*ac) + p2*(b*bc) + p3*(a*b + ac*bc)];

(*---------------------------------*)
d1=((1+T)*p1+R(p1*p2-p3^2))/(1+T+R*p2);
d2=((1+T)*p2+R)/(1+T+R*p2);
d3=(2Sqrt[T]*p3)/(1+T+R*p2);

Fz[a_,ac_,b_,bc_]=Exp[-(d1/2)a*ac-d2/2 b*bc+d3/2 (a*b+ac*bc)];
(*-----------------------PS-----------------------*)
CFps[a_,ac_,b_,bc_]=\
-((1-T)/T)Exp[-(b*bc)/2] * D[Exp[(b*bc)/2] * Fz[a,ac,b,bc],b,bc];

Pps=Simplify[CFps[0,0,0,0]];
CFps[a_,ac_,b_,bc_]=CFps[a,ac,b,bc]/Pps;

(*-----------------------PA-----------------------*)
CFpa[a_,ac_,b_,bc_]=\
-(1-T)Exp[(b*bc)/2] * D[Exp[-(b*bc)/2] * Fz[a,ac,b,bc],b,bc];

Ppa=Simplify[CFpa[0,0,0,0]];
CFpa[a_,ac_,b_,bc_]=Faa[a,ac,b,bc]/Ppa;

(*-----------------------PC-----------------------*)
Tp=(T-1)/T;
CFpc[a_,ac_,b_,bc_]=\
Tp^2*Exp[(b*bc)/2] * D[Exp[-(b*bc)] * D[Exp[(b*bc)/2] * Fz[a,ac,b,bc],b,bc],b,bc]+\
-Tp*Exp[(b*bc)/2] * D[Exp[-(b*bc)] * D[Exp[(b*bc)/2] * Fz[a,ac,b,bc],b],bc]+\
-Tp*Exp[(b*bc)/2] * D[Exp[-(b*bc)] * D[Exp[(b*bc)/2] * Fz[a,ac,b,bc],bc],b]+\
Fz[a,ac,b,bc];

Ppc=CFpc[0,0,0,0];

CFpc[a_,ac_,b_,bc_]=CFpc[a,ac,b,bc]/Ppc;

(*-----------------------PS-PA (a+a)-----------------------*)
Tp=(T-1)/T;
CFsa[a_,ac_,b_,bc_]=\
Tp^2*Exp[(b*bc)/2] * D[Exp[-(b*bc)] * D[Exp[(b*bc)/2] * Fz[a,ac,b,bc],b,bc],b,bc];

Psa=Simplify[CFsa[0,0,0,0]];
CFsa[a_,ac_,b_,bc_]=CFsa[a,ac,b,bc]/Psa;
			  						  					  
(*-----------------------PA-PS (aa+)-----------------------*)
Tp=(T-1);
CFas[a_,ac_,b_,bc_]=\
Tp^2*Exp[-(b*bc)/2] * D[Exp[(b*bc)] * D[Exp[-(b*bc)/2] * Fz[a,ac,b,bc],b,bc],b,bc];

Pas=Simplify[CFas[0,0,0,0]];
CFas[a_,ac_,b_,bc_]=CFas[a,ac,b,bc]/Pas;


(*------------------ Swap CF -----------------------------*)

CFin[a_,ac_,b_,bc_] = CFgauss[a,ac,b,bc];
CFout[a_, ac_, b_, bc_] = CFin[a, ac, eta*(a*g1 + b*g2), eta*(-ac*g1 + bc*g2)] 
						  CFin[eta*(a*g1 + b*g2), eta*(ac*g1-bc*g2),b, bc]
						  Exp[-(1/2)*(g1^2*(a^2 +ac^2) + g2^2*(b^2+bc^2))*(1-eta^2)];

(*------- Ideal swap ---------------*)
(*CFout[a_, ac_, b_, bc_] = CFin[a, ac, b, bc]CFin[b, -bc,b, bc];*)

(*Expected values of a and adag*)
CF = CFout;
DCFa[a_, ac_, b_, bc_] = -D[CF[a,ac,b,bc],ac];
Print["-------------------- a --------------------"];
ToMatlab[Simplify[DCFa[0,0,0,0]]]
				
DCFadag[a_, ac_, b_, bc_] = -D[CF[a,ac,b,bc],a];
Print["-------------------- a^dag --------------------"];
ToMatlab[Simplify[DCFadag[0,0,0,0]]]
				
DCFa2[a_, ac_, b_, bc_] = -D[CF[a,ac,b,bc],{ac,2}];
Print["-------------------- a^2 --------------------"];
ToMatlab[Simplify[DCFa2[0,0,0,0]]]
			
DCFadag2[a_, ac_, b_, bc_] = D[CF[a,ac,b,bc],{a,2}];
Print["-------------------- a^dag^2 --------------------"];
ToMatlab[Simplify[DCFadag2[0,0,0,0]]]
				
DCFadaga[a_, ac_, b_, bc_] = -D[CF[a,ac,b,bc],a,ac];
Print["-------------------- a^dag a --------------------"];
ToMatlab[Simplify[DCFadaga[0,0,0,0]]]
					
DCFb[a_, ac_, b_, bc_] = -D[CF[a,ac,b,bc],bc];
Print["-------------------- b --------------------"];
ToMatlab[Simplify[DCFb[0,0,0,0]]]
				
DCFbdag[a_, ac_, b_, bc_] = -D[CF[a,ac,b,bc],b];
Print["-------------------- b^dag --------------------"];
ToMatlab[Simplify[DCFbdag[0,0,0,0]]]
										
DCFb2[a_, ac_, b_, bc_] = -D[CF[a,ac,b,bc],{bc,2}];
Print["-------------------- b^2 --------------------"];
ToMatlab[Simplify[DCFb2[0,0,0,0]]]
				
DCFbdag2[a_, ac_, b_, bc_] = D[CF[a,ac,b,bc],{b,2}];
Print["-------------------- b^dag^2 --------------------"];
ToMatlab[Simplify[DCFbdag2[0,0,0,0]]]
				
DCFbdagb[a_, ac_, b_, bc_] = -D[CF[a,ac,b,bc],b,bc];
Print["-------------------- b^dag b --------------------"];
ToMatlab[Simplify[DCFbdagb[0,0,0,0]]]

DCFab[a_, ac_, b_, bc_] = D[CF[a,ac,b,bc],ac,bc];
Print["-------------------- a b --------------------"];
ToMatlab[Simplify[DCFab[0,0,0,0]]]

DCFadagbdag[a_, ac_, b_, bc_] = D[CF[a,ac,b,bc],a,b];
Print["-------------------- a^dag b^dag --------------------"];
ToMatlab[Simplify[DCFadagbdag[0,0,0,0]]]

DCFadagb[a_, ac_, b_, bc_] = -D[CF[a,ac,b,bc],a,bc];
Print["-------------------- a^dag b --------------------"];
ToMatlab[Simplify[DCFadagb[0,0,0,0]]]

DCFabdag[a_, ac_, b_, bc_] = -D[CF[a,ac,b,bc],ac,b];
Print["-------------------- a b^dag --------------------"];
ToMatlab[Simplify[DCFabdag[0,0,0,0]]]










