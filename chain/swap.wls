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
CFout[a_, ac_, b_, bc_] = CFch[a, ac, eta*(a*g1 + b*g2), eta*(-ac*g1 + bc*g2)] CFch[eta*(a*g1 + b*g2), eta*(ac*g1-bc*g2), eta*(a*g1 + b*g2),b, bc]Exp[-(1/2)*(g1^2*(a^2 +ac^2) + g2^2*(b^2+bc^2))*(1-eta^2)];
						  
(*Expected values of a and adag*)
CF = CFout;
DCFa = -D[CF[a,ac,b,bc],ac];
valDCFa = DCFa[0,0,0,0];
Print["-------------------- a --------------------"];
ToMatlab[valDCFa];
				
DCFadag = -D[CF[a,ac,b,bc],a];
valDCFadag = DCFadag[0,0,0,0];
Print["-------------------- a^dag --------------------"];
ToMatlab[valDCFadag];
				
DCFa2 = -D[CF[a,ac,b,bc],{ac,2}];
valDCFa2 = DCFa2[0,0,0,0];
Print["-------------------- a^2 --------------------"];
ToMatlab[valDCFa2];
			
DCFadag2 = D[CF[a,ac,b,bc],{a,2}];
valDCFadag2 = DCFadag2[0,0,0,0];
Print["-------------------- a^dag^2 --------------------"];
ToMatlab[valDCFadag2];
				
DCFadaga = -D[CF[a,ac,b,bc],a,ac];
valDCFadaga = DCFadaga[0,0,0,0];
Print["-------------------- a^dag a --------------------"];
ToMatlab[valDCFadaga];
					
CFb = -D[CF[a,ac,b,bc],bc];
valDCFb = DCFb[0,0,0,0];
Print["-------------------- b --------------------"];
ToMatlab[valDCFb];
				
DCFbdag = -D[CF[a,ac,b,bc],b];
valDCFbdag = DCFbdag[0,0,0,0];
Print["-------------------- b^dag --------------------"];
ToMatlab[valDCFbdag];
										
DCFb2 = -D[CF[a,ac,b,bc],{bc,2}];
valDCFb2 = DCFb2[0,0,0,0];
Print["-------------------- b^2 --------------------"];
ToMatlab[valDCFb2];
				
DCFbdag2 = D[CF[a,ac,b,bc],{b,2}];
valDCFbdag2 = DCFadag2[0,0,0,0];
Print["-------------------- b^dag^2 --------------------"];
ToMatlab[valDCFbdag2];
				
DCFbdagb = -D[CF[a,ac,b,bc],b,bc];
valDCFbdagb = DCFbdagb[0,0,0,0];
Print["-------------------- b^dag b --------------------"];
ToMatlab[valDCFbdagb];


Simplify[DCFa[0,0,0,0]]
