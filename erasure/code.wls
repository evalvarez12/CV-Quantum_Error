#!/usr/bin/env wolframscript
(* ::Package:: *)

Clear["Global`*"];
Import["/home/eduardo/QuantumComms/CV-Quantum_Error_Correction/erasure/ToMatlab.m", "Package"];
Import["/home/eduardo/QuantumComms/CV-Quantum_Error_Correction/chain/cm.wls", "Package"];

V=.;
gx=.;
gp=.;

CFtmsv[xa_, pa_, xb_, pb_] = Exp[-V/2*(xa^2 + pa^2 + xb^2 + pb^2) + Sqrt[V^2-1]*(xa*xb - pa*pb)];
CF0[x_, p_] = Exp[-1/2*(x^2 + p^2)];
CFalpha[x_, p_] = Exp[-1/2*(x^2 + p^2) + 2*i*(p*ar - x*ai)];


CFs[a_, ac_] = ...
CFab[a_, ac_, b_, bc_] = ...
CF0[a_, a_c] = Exp[(a*ac)/2];

(*------- Mode 1 lost --------------*)
CF1[x_, p_] = CFs[(1-gx)*x/2, (1+gp)*p/2]*CFab[-(1-gx)*x/2, -(1+gp)*p/2, -gx*x, -gp*p]*CF0[(1+gx)*x/Sqrt[2], (1-gp)*p/Sqrt[2]];

(*------- Mode 2 lost --------------*)
CF2[x_, p_] = CFs[(1+gx)*x/2, (1-gp)*p/2]*CFab[(1+gx)*x/2, (1-gp)*p/2, -gx*x, -gp*p]*CF0[(1-gx)*x/Sqrt[2], (1+gp)*p/Sqrt[2]];

(*------- Mode 3 lost --------------*)
CF3[x_, p_] = CFs[x, p]*CFab[gx*x, -gp*p, 0, 0]*CF0[-gx*x, -gp*p];

(*------- Mode 1,2 lost --------------*)
CF12[x_, p_] = CF0[x, p];

(*------- Mode 1,3 lost --------------*)
CF13[x_, p_] = CFs[(1-gx)*x/2, (1+gp)*p/2]*CFab[-(1-gx)*x/2, -(1+gp)*p/2, 0, 0]*CF0[(1+gx)*x/2, (1-gp)*p/2]*CF0[(1-gx)*x/2, (1+gp)*p/2];

(*------- Mode 2,3 lost --------------*)
CF23[x_, p_] = CFs[(1+gx)*x/2, (1-gp)*p/2]*CFab[(1+gx)*x/2, (1-gp)*p/2, 0, 0]*CF0[(1-gx)*x/Sqrt[2], (1+gp)*p/Sqrt[2]]*CF0[-gx*x, -gp*p];









