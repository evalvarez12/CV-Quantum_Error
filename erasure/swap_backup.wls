#!/usr/bin/env wolframscript
(* ::Package:: *)

Clear["Global`*"];
Import["/home/eduardo/QuantumComms/CV-Quantum_Error_Correction/chain/ToMatlab.m", "Package"];
Import["/home/eduardo/QuantumComms/CV-Quantum_Error_Correction/chain/cm.wls", "Package"];

t=.;
g1=.;
g2=.;
e=.;
T=.;
V=.;
to=.;

(*----------------EPR-----------------*)
CFepr[a_, ac_, b_, bc_] = Exp[-(1/2)*((V (a*ac + b*bc)) - Sqrt(V^2 - 1) * (a*b + ac*bc))];
CFch[a_, ac_, b_, bc_] = Exp[-(1+e)*(1-T)/2]CFepr[a, ac, b*Sqrt[T], bc*Sqrt[T]];

(*Swap CF*)
CFout[a_, ac_, b_, bc_] = CFch[a, ac, eta*(a*g1 + b*g2), eta*(-ac*g1 + bc*g2)] CFch[eta*(a*g1 + b*g2), eta*(ac*g1-bc*g2),b, bc]Exp[-(1/2)*(g1^2*(a^2 +ac^2) + g2^2*(b^2+bc^2))*(1-eta^2)];
	
CM[CFout]
						  					  						  					  
(*(*Expected values of a and adag*)
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
ToMatlab[Simplify[DCFbdagb[0,0,0,0]]]*)



