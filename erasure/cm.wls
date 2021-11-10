#!/usr/bin/env wolframscript
(* ::Package:: *)

BeginPackage["cm`"];
Begin["`Private`"];


CM[CF0_] := Module[{CF=CF0},
(*Expected values of a and adag*)

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
]
				
End[];
EndPackage[];
