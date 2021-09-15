#!/usr/bin/env wolframscript
(* ::Package:: *)

BeginPackage["cm`"];
Begin["`Private`"];


CM[CF0_] := Module[{CF=CF0},
				CFa = -D[CF[a,ac,b,bc],ac];
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
				ToMatlab[valDCFbdagb];]
				
				
End[];
EndPackage[];
