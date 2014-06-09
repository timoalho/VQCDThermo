(* ::Package:: *)

(***********************************************
VQCDCore.m

Core code for solving the IHQCD + tachyon model background differential equations in Mathematica

Written by Timo Alho, based on arxiv:1112.1261, arxiv:1210.4516 and arxiv:1312.5199.

Package created on 10.11.2011
************************************************)
BeginPackage["VQCDThermo`VQCDCore`", 
{"VQCDThermo`VQCDInitialConditions`", "VQCDThermo`NumericsExtra`","VQCDThermo`Verbosity`", "DifferentialEquations`InterpolatingFunctionAnatomy`"}]

VQCDEquationsOfMotion::usage = "VQCDEquationsOfMotion[{Vf, Vg, \[Kappa], \[Omega]}]\n\n
 returns the equations of motion with the given potentials."

VQCDBHInitialConditions::usage = "VQCDBHInitialConditions[\[Tau]h, \[Lambda]h, bh, fdh, {Vf, Vg, \[Kappa], \[Omega]}, \[Epsilon]]\n
returns the initial conditions at \[Epsilon]0 for a black hole solution where the horizon is at 0. The initial conditions are generated to assure the regularity of the solution at the horizon.

Options and their defaults:
EpsilonMeaning -> \"fd\", possible values \"fd\" or \"Absolute\". With fd, \[Epsilon]0 is defined such that \[Epsilon] is the maximum relative change in the fields or their derivatives during the first step.
With Absolute, \[Epsilon]0 = \[Epsilon].
"

VQCDExtremalBHInitialConditions::usage = "Computes the initial conditions with an extremal horizon. NOTE: only EpsilonValue -> \"Absolute\" works currently."

SolveVQCDBH::usage = "SolveVQCDBlackHole[\[Lambda]h, \[Tau]h, nt, {Vg, Vf, \[Kappa], \[Omega]}] solves the EoM's with the given potentials and \[Lambda]h, \[Tau]h, nt for a black hole with horizon at \[Lambda]h. The other boundary conditions are given in the options.

Returns a list {sols, success}, where sols is the solution in the format returned by NDSolve, and success == True if the solution was found and it extends all the way to the UV AdS -singularity, false otherwise (the part of the solution that was reached is still returned in sols).)

Options and their defaults:
HorizonEpsilon \[Rule] 10^(-5), The distance from the horizon at which to set the initial conditions\[IndentingNewLine]fpInitialValue \[Rule] 1, Initial value for df/dA at the horizon\[IndentingNewLine]ARange \[Rule] 50, The range of A-values to solve for. Should be large enough to guarantee that the UV asymptote is reached, but too large values may cause a too large initial step
In addition any options for NDSolve or VQCDBHInitialConditions may be given.
"
SolveAndScaleVQCDBH::usage = "SolveAndScaleVQCDBlackHole[\[Lambda]h, \[Tau]h, nt, {Vg, Vf, \[Kappa], \[Omega]}] is the same as SolveVQCDBlackHole[\[Lambda]h, \[Tau]h, nt, {Vg, Vf, \[Kappa], \[Omega]}], except that the output is automatically scaled by ScaleSolutions. In addition the output is returned as a list of functions instead of replacement rules.

Note that any imaginary parts are stripped from the solutions. There will be warning messages if they are large either in absolute value or relative to the real part.

Returns a list {q, f, \[Lambda], \[Tau]d, \[Tau], Amin, Amax, \[CapitalLambda]scale, fscale, {b0, b1}, ell}, where the solutions are scaled such that \[CapitalLambda]scale is the scaling factor needed to bring the UV -scale to \[CapitalLambda] = 1 and b0, b1 are the \[Beta]-function expansion parameters. success is True if the solution extends to the requested A-range, false otherwise.
It is recommended to use the ...FromSols -functions to extract observables from this list, in order to proof against future changes.

For options see SolveVQCDBlackHole.
 "
SolveVQCDExtremalBH::usage = "SolveVQCDExtremalBH[C\[Lambda]_?NumericQ, {Vg, Vf, \[Kappa], \[Omega]}] finds an extremal black hole solution, with the size of the non-analytic terms being determined by C\[Lambda]"

SolveAndScaleExtremalTachyons::usage = "SolveAndScaleVQCDExtremalBH[C\[Lambda]_?NumericQ, {Vg, Vf, \[Kappa], \[Omega]}] finds an extremal black hole solution and scales it to the correct units, with the size of the non-analytic terms being determined by C\[Lambda]"

SolverCore::usage = "The core of the black hole EoM -solver, directly calling this should only be necessary for debugging."

ScaleSolution::usage = "ScaleSolution[q_, f_, \[Lambda]_, \[Tau]d_, \[Tau]_, Amin_, Amax_, betacoefs_List, ell_, \[CapitalLambda]_] scales the functions such that f(Amax) = 1 and \[Lambda](z) asymptotes to the two-loop IHQCD -formula with the scale \[CapitalLambda].

Returns list {q, f, \[Lambda], \[Tau]d, \[Tau], Amin, Amax, \[CapitalLambda]scale, fscale}, where \[CapitalLambda]scale is the scaling factor needed to bring the UV -scale to \[CapitalLambda] = 1 and fscale is the factor needed to bring f(Amax) = 1. Note that fscale has already been applied to the returned functions, whereas the \[CapitalLambda] -scaling has not.

ScaleSolution[{sols, finished}, \[CapitalLambda]] does the same when  finished is true, sols being a list of replacement rules in the format given by NDSolve.

Options and their defaults:
lambdaFitRMSLimit -> 10^(-7), a limit for the RMS difference between the actual \[Lambda](z) and the fitted form, in the area defined by FitRange. If the limit is exceeded, a warning is issued but the calculation continues without interruption.
WorkingPrecision -> Automatic, not used yet. To be supported since solving the EoM's to higher precision requires also scaling to a higher precision.
"

DefineOldPotI::usage = "DefineOldPotI[xf] gives potential I, old style W"
DefineJKPotI::usage = "DefineJKPotI[xf] gives potential I"
DefineJKPotIMod::usage = "DefineOJKPotIMod[xf] gives modified potential I"
DefineJKPotII::usage = "DefineJKPotII[xf] gives potential II"
DefineJKPotIIMod::usage = "DefineJKPotIIMod[xf] gives potential II"
DefineJKPotIII::usage = "DefineJKPotIII[xf] gives potential III"
DefineJKPotIKappaMod::usage = "DefineJKPotIKappaMod[xf] gives potential I with modified Kappa"
DefineAIJKOptPot::usage = "DefineAIJKOptPot[xf] defines the potential given as optimal in 1309.2286, around eq. (4.26)"


VEffective::usage = "Returns the effective potential, given the full potentials as input."
ntCritical::usage = "Returns the limiting nt given the potentials"

PotFromAnsatzAndBeta::usage = "Given a form for the potential functions and a beta function, produces the potentials which give a matching beta power series"

bCoefsFromPotential::usage = "Computes the beta coefficients given a potential V[\[Lambda]]"

GammaFromKappa::usage = "Computes \[Gamma] given \[Kappa] and b0"

QuarkMass::usage = "Computes mq."

\[Tau]hFromQuarkMass::usage = "Computes \[Tau]h given mq, \[Lambda]h."

qFromSols::usage = "A helper to extract q given a scaled solution to the EoMs"
\[Lambda]FromSols::usage = "A helper to extract \[Lambda] given a scaled solution to the EoMs"
bFromSols::usage = "A helper to extract b = Exp(A) given a scaled solution to the EoMs"
fFromSols::usage = "A helper to extract f given a scaled solution to the EoMs"
\[Tau]FromSols::usage = "A helper to extract \[Tau] given a scaled solution to the EoMs"
zFromSols::usage = "A helper to extract z(A) given a scaled solution to the EoMs"

\[CapitalLambda]ScaleFromSols::usage = "A helper to extract \[CapitalLambda] given a scaled solution to the EoMs"
fScaleFromSols::usage = "A helper to extract the f-scaling factor 1/Sqrt[f_UV] give a scaled solution to the EoMs"
qhFromSols::usage = "A helper to extract the value of q(A) at horizon, given a scaled solution to the EoMs"
qhFromThermoData::usage = "A helper to extract qh, given the relevant boundary scaling parameters"
TemperatureFromThermoData::usage = "A helper to extract the temperature, given the boundary scaling parameters"
APrimeFromSols::usage = "A helper to extract the A-derivative of the gauge field, given a scaled solution to the EoMs"
AFromSols::usage = "A helper to extract the gauge field, given a scaled solution to the EoMs"
AAndMuFromSols::usage = "A helper to extract the gauge field and , given a scaled solution to the EoMs"
TemperatureFromSols::usage = "A helper to extract the temperature, given a scaled solution to the EoMs"
s4G5FromSols::usage = "A helper to extract the entropy density, given a scaled solution to the EoMs"
ntildeFromSols::usage = "A helper to extract the effective ntilde, given a scaled solution to the EoMs"

ARangeFromSols::usage = "Returns the range in A on which the solutions are defined"
AMinFromSols::usage = "Returns the lower limit of the A-range on which the solutions are defined"
AMaxFromSols::usage = "Returns the upper limit of the A-range on which the solutions are defined"

(*Export the symbols appearing in the EoM's.
Probably bad form, but without these the output consists of terms such as FiniteTTachyons`Private`\[Lambda][FiniteTTachyons`Private`z]

Any suggestions to deal with the problem more elegantly are welcome.
*)
\[Lambda];
q;
f;
\[Tau];
A;


Begin["`Private`"]
(*Returns the equations of motion given Vf(\[Lambda], \[Tau]), Vg(\[Lambda]), \[Kappa](\[Lambda]) and w(\[Lambda])*)
VQCDEquationsOfMotion[{Vg_Symbol, Vf_Symbol, \[Kappa]_Symbol, w_Symbol}, nt_] := {Derivative[2][\[Tau]][A] == 
  ((-3*q[A]^2*\[Lambda][A]^2*(3*E^(3*A)*w[\[Lambda][A]]*(nt^2 + E^(6*A)*Vf[\[Lambda][A], \[Tau][A]]^2*
          w[\[Lambda][A]]^2)*\[Kappa][\[Lambda][A]]*Derivative[1][f][A]*Derivative[1][\[Tau]][A] + 
       q[A]*(nt^2 + E^(6*A)*Vf[\[Lambda][A], \[Tau][A]]^2*w[\[Lambda][A]]^2)^(3/2)*\[Kappa][\[Lambda][A]]*
        Derivative[1][\[Tau]][A]*Sqrt[q[A]^2 + f[A]*\[Kappa][\[Lambda][A]]*Derivative[1][\[Tau]][A]^
            2] + E^(3*A)*q[A]^2*w[\[Lambda][A]]*(nt^2*Vg[\[Lambda][A]]*\[Kappa][\[Lambda][A]]*
          Derivative[1][\[Tau]][A] + E^(6*A)*Vf[\[Lambda][A], \[Tau][A]]*w[\[Lambda][A]]^2*
          (Vf[\[Lambda][A], \[Tau][A]]*Vg[\[Lambda][A]]*\[Kappa][\[Lambda][A]]*Derivative[1][\[Tau]][A] - 
           6*Derivative[0, 1][Vf][\[Lambda][A], \[Tau][A]]))))/E^(3*A) - 
    9*f[A]^2*\[Kappa][\[Lambda][A]]*\[Lambda][A]^2*Derivative[1][\[Tau]][A]^3*
     (nt^2*(-2*\[Kappa][\[Lambda][A]]*Derivative[1][w][\[Lambda][A]]*Derivative[1][\[Lambda]][A] + 
        w[\[Lambda][A]]*(2*\[Kappa][\[Lambda][A]] + Derivative[1][\[Kappa]][\[Lambda][A]]*Derivative[1][\[Lambda]][
            A])) + E^(6*A)*Vf[\[Lambda][A], \[Tau][A]]*w[\[Lambda][A]]^3*
       (Vf[\[Lambda][A], \[Tau][A]]*(8*\[Kappa][\[Lambda][A]] + Derivative[1][\[Kappa]][\[Lambda][A]]*
           Derivative[1][\[Lambda]][A]) + 2*\[Kappa][\[Lambda][A]]*Derivative[1][\[Lambda]][A]*
         Derivative[1, 0][Vf][\[Lambda][A], \[Tau][A]])) + f[A]*Derivative[1][\[Tau]][A]*
     (-9*w[\[Lambda][A]]*(nt^2 + E^(6*A)*Vf[\[Lambda][A], \[Tau][A]]^2*w[\[Lambda][A]]^2)*\[Kappa][\[Lambda][A]]^2*
       \[Lambda][A]^2*Derivative[1][f][A]*Derivative[1][\[Tau]][A]^2 + 
      2*q[A]^2*(nt^2*(9*w[\[Lambda][A]]*\[Kappa][\[Lambda][A]]*\[Lambda][A]^2 + 
          9*\[Lambda][A]^2*(\[Kappa][\[Lambda][A]]*Derivative[1][w][\[Lambda][A]] - 
            w[\[Lambda][A]]*Derivative[1][\[Kappa]][\[Lambda][A]])*Derivative[1][\[Lambda]][A] + 
          2*w[\[Lambda][A]]*\[Kappa][\[Lambda][A]]*Derivative[1][\[Lambda]][A]^2) - 
        E^(6*A)*Vf[\[Lambda][A], \[Tau][A]]*w[\[Lambda][A]]^3*(Vf[\[Lambda][A], \[Tau][A]]*
           (9*\[Lambda][A]^2*Derivative[1][\[Kappa]][\[Lambda][A]]*Derivative[1][\[Lambda]][A] + 
            2*\[Kappa][\[Lambda][A]]*(9*\[Lambda][A]^2 - Derivative[1][\[Lambda]][A]^2)) + 
          9*\[Kappa][\[Lambda][A]]*\[Lambda][A]^2*(-(Derivative[1][\[Tau]][A]*Derivative[0, 1][Vf][\[Lambda][
                A], \[Tau][A]]) + Derivative[1][\[Lambda]][A]*Derivative[1, 0][Vf][\[Lambda][A], 
              \[Tau][A]])))))/(18*f[A]*q[A]^2*w[\[Lambda][A]]*
    (nt^2 + E^(6*A)*Vf[\[Lambda][A], \[Tau][A]]^2*w[\[Lambda][A]]^2)*\[Kappa][\[Lambda][A]]*\[Lambda][A]^2), 
 Derivative[1][q][A] == (q[A]*(12 + (4*Derivative[1][\[Lambda]][A]^2)/(3*\[Lambda][A]^2) + 
     (3*Derivative[1][f][A] + q[A]^2*(-Vg[\[Lambda][A]] + Vf[\[Lambda][A], \[Tau][A]]*
          Sqrt[1 + nt^2/(E^(6*A)*Vf[\[Lambda][A], \[Tau][A]]^2*w[\[Lambda][A]]^2)]*
          Sqrt[1 + (f[A]*\[Kappa][\[Lambda][A]]*Derivative[1][\[Tau]][A]^2)/q[A]^2]))/f[A]))/6, 
 Derivative[2][\[Lambda]][A] == (-6*E^(3*A)*q[A]^4*w[\[Lambda][A]]^2*\[Lambda][A]^2*
     (9*\[Lambda][A]^2*Derivative[1][Vg][\[Lambda][A]] + 4*Vg[\[Lambda][A]]*Derivative[1][\[Lambda]][A]) + 
    8*E^(3*A)*f[A]*w[\[Lambda][A]]^2*\[Kappa][\[Lambda][A]]*Derivative[1][\[Lambda]][A]*
     (-9*\[Lambda][A]^2*Derivative[1][f][A] - 2*f[A]*(3*\[Lambda][A] - 2*Derivative[1][\[Lambda]][A])*
       (6*\[Lambda][A] + Derivative[1][\[Lambda]][A]))*Derivative[1][\[Tau]][A]^2 - 
    2*E^(3*A)*q[A]^2*w[\[Lambda][A]]^2*(4*Derivative[1][\[Lambda]][A]*
       (9*\[Lambda][A]^2*Derivative[1][f][A] + 2*f[A]*(3*\[Lambda][A] - 2*Derivative[1][\[Lambda]][A])*
         (6*\[Lambda][A] + Derivative[1][\[Lambda]][A])) + 3*f[A]*\[Kappa][\[Lambda][A]]*\[Lambda][A]^2*
       (9*\[Lambda][A]^2*Derivative[1][Vg][\[Lambda][A]] + 4*Vg[\[Lambda][A]]*Derivative[1][\[Lambda]][A])*
       Derivative[1][\[Tau]][A]^2) + 6*q[A]^3*\[Lambda][A]^2*
     (-4*w[\[Lambda][A]]*Sqrt[nt^2 + E^(6*A)*Vf[\[Lambda][A], \[Tau][A]]^2*w[\[Lambda][A]]^2]*
       Derivative[1][\[Lambda]][A]*Sqrt[q[A]^2 + f[A]*\[Kappa][\[Lambda][A]]*Derivative[1][\[Tau]][A]^
           2] + 9*nt^2*\[Lambda][A]^2*Derivative[1][w][\[Lambda][A]]*
       Sqrt[(q[A]^2 + f[A]*\[Kappa][\[Lambda][A]]*Derivative[1][\[Tau]][A]^2)/
         (nt^2 + E^(6*A)*Vf[\[Lambda][A], \[Tau][A]]^2*w[\[Lambda][A]]^2)] - 
      9*E^(6*A)*Vf[\[Lambda][A], \[Tau][A]]*w[\[Lambda][A]]^3*\[Lambda][A]^2*
       Sqrt[(q[A]^2 + f[A]*\[Kappa][\[Lambda][A]]*Derivative[1][\[Tau]][A]^2)/
         (nt^2 + E^(6*A)*Vf[\[Lambda][A], \[Tau][A]]^2*w[\[Lambda][A]]^2)]*
       Derivative[1, 0][Vf][\[Lambda][A], \[Tau][A]]) - 3*f[A]*q[A]*\[Lambda][A]^2*
     Derivative[1][\[Tau]][A]^2*(-18*nt^2*\[Kappa][\[Lambda][A]]*\[Lambda][A]^2*Derivative[1][w][\[Lambda][A]]*
       Sqrt[(q[A]^2 + f[A]*\[Kappa][\[Lambda][A]]*Derivative[1][\[Tau]][A]^2)/
         (nt^2 + E^(6*A)*Vf[\[Lambda][A], \[Tau][A]]^2*w[\[Lambda][A]]^2)] + 
      w[\[Lambda][A]]*(8*Sqrt[nt^2 + E^(6*A)*Vf[\[Lambda][A], \[Tau][A]]^2*w[\[Lambda][A]]^2]*\[Kappa][\[Lambda][A]]*
         Derivative[1][\[Lambda]][A]*Sqrt[q[A]^2 + f[A]*\[Kappa][\[Lambda][A]]*Derivative[1][\[Tau]][A]^
             2] + 9*nt^2*\[Lambda][A]^2*Derivative[1][\[Kappa]][\[Lambda][A]]*
         Sqrt[(q[A]^2 + f[A]*\[Kappa][\[Lambda][A]]*Derivative[1][\[Tau]][A]^2)/
           (nt^2 + E^(6*A)*Vf[\[Lambda][A], \[Tau][A]]^2*w[\[Lambda][A]]^2)]) + 
      9*E^(6*A)*Vf[\[Lambda][A], \[Tau][A]]*w[\[Lambda][A]]^3*\[Lambda][A]^2*
       Sqrt[(q[A]^2 + f[A]*\[Kappa][\[Lambda][A]]*Derivative[1][\[Tau]][A]^2)/
         (nt^2 + E^(6*A)*Vf[\[Lambda][A], \[Tau][A]]^2*w[\[Lambda][A]]^2)]*
       (Vf[\[Lambda][A], \[Tau][A]]*Derivative[1][\[Kappa]][\[Lambda][A]] + 2*\[Kappa][\[Lambda][A]]*
         Derivative[1, 0][Vf][\[Lambda][A], \[Tau][A]])))/(144*E^(3*A)*f[A]*w[\[Lambda][A]]^2*
    \[Lambda][A]^2*(q[A]^2 + f[A]*\[Kappa][\[Lambda][A]]*Derivative[1][\[Tau]][A]^2)), 
 Derivative[2][f][A] == Derivative[1][f][A]^2/(2*f[A]) - 
   (nt^2*q[A]*Sqrt[q[A]^2 + f[A]*\[Kappa][\[Lambda][A]]*Derivative[1][\[Tau]][A]^2])/
    (E^(3*A)*w[\[Lambda][A]]*Sqrt[nt^2 + E^(6*A)*Vf[\[Lambda][A], \[Tau][A]]^2*w[\[Lambda][A]]^2]) + 
   Derivative[1][f][A]*(-2 + (2*Derivative[1][\[Lambda]][A]^2)/(9*\[Lambda][A]^2) + 
     (q[A]^2*(-Vg[\[Lambda][A]] + Vf[\[Lambda][A], \[Tau][A]]*
         Sqrt[1 + nt^2/(E^(6*A)*Vf[\[Lambda][A], \[Tau][A]]^2*w[\[Lambda][A]]^2)]*
         Sqrt[1 + (f[A]*\[Kappa][\[Lambda][A]]*Derivative[1][\[Tau]][A]^2)/q[A]^2]))/(6*f[A]))};


(*Constructs the full set of initial conditions at distance \[Epsilon] from the horizon by requiring regularity of the solution at the horizon*)
Options[VQCDExtremalBHInitialConditions] = {EpsilonMeaning->"fd" (*Sets the meaning of epsilon. fd means that it is the relative change in f', absolute means it is the actual offset from horizon in z*)
											, EpsilonLimit -> Infinity
											, nt1rootrange -> {0.5, 1.5} (*The area where to search for the critical points of the potential*)
}
VQCDExtremalBHInitialConditions[bdh_?NumericQ, C\[Lambda]_?NumericQ, {Vg_, Vf_, \[Kappa]_, w_}, \[Epsilon]_?NumericQ, opts:OptionsPattern[]] := Module[
 {\[Lambda]hcrit, nt1crit, nt1critfun, \[Lambda]h, \[Lambda]dh, \[Lambda]ddh, \[Lambda]3h, \[Lambda]3s1, \[Lambda]3s2, \[Lambda]3s3, bddh, bdds1, bdds2, f3s, \[Tau]h, \[Tau]dh, \[Tau]ddh, fddh, f4h, f4s1, f4s2, \[Lambda]dds1, \[Lambda]dds2, b3h, b3s1, b3s2, b3s3, \[Epsilon]eff, \[Alpha], \[Beta], f3h, a, bh, bz, bdz, \[Lambda]z, \[Lambda]dz, fz, fdz, qz, Veff, A\[Epsilon]},
 Block[{$vcontext = "VQCDExtremalBHInitialConditions"},

(*First find the critical point that admits regular solutions with fph = 0*)
 nt1critfun[\[Lambda]h_] := w[\[Lambda]h ]Sqrt[Vg[\[Lambda]h]^2 - Vf[\[Lambda]h, 0]^2];
 {\[Lambda]hcrit, nt1crit} = {\[Lambda]h, nt1critfun[\[Lambda]h]} /. FindRoot[nt1critfun'[\[Lambda]h] == 0, Evaluate[Join[{\[Lambda]h}, OptionValue[nt1rootrange]]], Method -> "Brent"];
 PrintV[StringForm["Found an extremal point at \[Lambda]h = `1`, nt = `2`", \[Lambda]hcrit, nt1crit], "Checks"];

 \[Lambda]h = \[Lambda]hcrit;

 \[Tau]h = 0;

 {a, fddh, f3s, f4s2, f3h, f4s1, f4h, \[Lambda]dds2, \[Lambda]3s3, \[Lambda]dh, \[Lambda]dds1, \[Lambda]3s2, \[Lambda]ddh, \[Lambda]3s1, \[Lambda]3h, bh, bdds2, b3s3, bdds1, b3s2, bddh, b3s1, b3h, \[Tau]dh, \[Tau]ddh}
 = ExtremalInitialConditions[\[Lambda]h, Vg, Vf, \[Kappa], w, C\[Lambda], bdh];

 (*Veff[\[Lambda]h_] := Vg[\[Lambda]h] - Vf[\[Lambda]h, \[Tau]h]*Sqrt[1 + nt1crit^2/(Vf[\[Lambda]h, \[Tau]h]^2*\[Kappa][\[Lambda]h]^2)];*)

 (*a = 1/2(-3 + Sqrt[1 - (3 \[Lambda]h^2 Vg[\[Lambda]h])/(-Vf[\[Lambda]h, \[Tau]h]^2+Vg[\[Lambda]h]^2) Veff''[\[Lambda]h]]);*)

 PrintV[StringForm["Initial conditions are {fddh = `1`, f3s = `2`, f4s2 = `3`, f3h = `4`, f4s1 = `5`, f4h = `6`, \[Lambda]dds2 = `7`, \[Lambda]3s3 = `8`, \[Lambda]dh = `9`, \[Lambda]dds1 = `10`, \[Lambda]3s2 = `11`, \[Lambda]ddh = `12`, \[Lambda]3s1 = `13`\.1d, \[Lambda]3h = `14`, bh = `15`, bdds2 = `16` , b3s3 = `17`, bdh = `18`, bdds1 = `19`\.1d, b3s2 = `20`, bddh = `21`, b3s1 = `22`, b3h = `23`, \[Tau]dh = `24`, \[Tau]ddh = `25`}, a = `26`", 
	fddh, f3s, f4s2, f3h, f4s1, f4h, \[Lambda]dds2, \[Lambda]3s3, \[Lambda]dh, \[Lambda]dds1, \[Lambda]3s2, \[Lambda]ddh, \[Lambda]3s1, \[Lambda]3h, bh, bdds2, b3s3, bdh, bdds1, b3s2, bddh, b3s1, b3h, \[Tau]dh, \[Tau]ddh,
	a], "Debug"];


 If[OptionValue[EpsilonMeaning] == "Absolute",
{{
 (*epsilon is directly the derivative itself*)

 (*Convert to A -coordinates. Note that in z, \[Epsilon] is slightly negative, which inverts
  the sign of odd integer power terms*)
 (*Print[StringForm[]]*)
 bz = bh + bdds2/2 \[Epsilon]^(2 + 2 a) - 1/6 b3s3 \[Epsilon]^(3 + 3a) - bdh \[Epsilon] + bdds1/2 \[Epsilon]^(2 + a) - 1/6 b3s2 \[Epsilon]^(3 + 2a) + bddh/2 \[Epsilon]^2 - 1/6 b3s1 \[Epsilon]^(3 + a)- 1/6 b3h \[Epsilon]^3;
 qz = bz^2/(-(1 + a)bdds2 \[Epsilon]^(1 + 2a) + (3 + 3a)/6 b3s3 \[Epsilon]^(2 + 3a) + bdh + (3 + 2a)/6 b3s2 \[Epsilon]^(2 + 2a) - (1 + a/2)bdds1 \[Epsilon]^(1 + a) - \[Epsilon] bddh + (3 + a)/6 b3s1 \[Epsilon]^(2 + a) + 1/2 b3h \[Epsilon]^2);


 fz = fddh/2 \[Epsilon]^2 - f3s/6 \[Epsilon]^(3 + a) + f4s2/24 \[Epsilon]^(4 + 2 a) - f3h/6 \[Epsilon]^3 + f4s1/24 \[Epsilon]^(4 + a) + f4h/24 \[Epsilon]^4;
 fdz = -fddh \[Epsilon] + (3 + a)f3s/6 \[Epsilon]^(2 + a) - (4+ 2a)f4s2/24 \[Epsilon]^(3 + 2a) + f3h/2 \[Epsilon]^2 - (4 + a)f4s1/24 \[Epsilon]^(3 + a) - f4h/6 \[Epsilon]^3;
 \[Lambda]z = \[Lambda]h - C\[Lambda] \[Epsilon]^(1 + a) + \[Lambda]dds2/2 \[Epsilon]^(2 + 2a) - \[Lambda]3s3/6  \[Epsilon]^(3 + 3a) - \[Lambda]dh \[Epsilon] + \[Lambda]dds1/2 \[Epsilon]^(2 + a) - \[Lambda]3s2/6 \[Epsilon]^(3 + 2a) + \[Lambda]ddh/2 \[Epsilon]^2 - \[Lambda]3s1/6 \[Epsilon]^(3 + a) - \[Lambda]3h/6 \[Epsilon]^3;
 \[Lambda]dz = (1 + a)C\[Lambda] \[Epsilon]^a - (2 + 2a)/2 \[Lambda]dds2 \[Epsilon]^(1 + 2a) + (3 + 3a)\[Lambda]3s3/6 \[Epsilon]^(2 + 3a) + \[Lambda]dh - (2 + a)\[Lambda]dds1/2 \[Epsilon]^(1 + a) + (3 + 2a)\[Lambda]3s2/6\[Epsilon]^(2 + 2a) - \[Lambda]ddh \[Epsilon] + (3 + a)\[Lambda]3s1 / 6 \[Epsilon]^(2 + a) + \[Lambda]3h/2 \[Epsilon]^2;

 A\[Epsilon] = Log[bz];

 PrintV[StringForm["\[Epsilon]_A = `1`", A\[Epsilon]], "Debug"];
 
 PrintV[StringForm["z boundary conds: b = `1`, b' = `2`, q = `3`, f = `4`, f' = `5`, \[Lambda] = `6`, \[Lambda]' = `7`\.1d", bz, bz^2/qz, qz, fz, fdz, \[Lambda]z, \[Lambda]dz], "Debug"];

 {\[Lambda][A\[Epsilon]] == \[Lambda]z,
  \[Lambda]'[A\[Epsilon]] == \[Lambda]dz qz/bz,
 q[A\[Epsilon]] == qz,
 f[A\[Epsilon]] ==  fz,
 f'[A\[Epsilon]] ==  fdz qz/bz,
 \[Tau][A\[Epsilon]] == 0,
 \[Tau]'[A\[Epsilon]] == 0
 },
 PrintV[StringForm["Initial conditions: \[Lambda] = `1`, \[Lambda]dh = `2`, q = `3`, f = `4`, fdh = `5`", \[Lambda]z, \[Lambda]dz qz / bz, qz, fz, fdz qz/bz], "Debug"];
 (fddh > 0) && (qz < 0) (*Return whether the initial conditions are viable*)
},
 nt1crit,
 Log[bz] (*\[Epsilon] in A*)
}
,
{
 (*\[Epsilon] is the maximum allowed relative change in any of the fields during the first step.*)


 (*The mess is related to making Min cope with the case \[Tau]ph/\[Tau]pph = 0/0*)
 \[Epsilon]eff = Min[OptionValue[EpsilonLimit], Block[{Indeterminate = Infinity},

 Min[1/2 Min[\[Epsilon] \[Lambda]h/Abs[\[Lambda]ph], (\[Epsilon] \[Lambda]h /Abs[C\[Lambda]])^(1/(1 + a))],
Quiet[1/3 Min[\[Epsilon] Abs[fpph]/Abs[f3h], (6 \[Epsilon] Abs[fpph]/((3 + a)(2 + a)Abs[f3s1]))^(1/(1 + a)),
			  (6 \[Epsilon] Abs[fpph]/((3 + 2 a)(2 + 2 a)Abs[f3s2]))^(1/(1 + 2a))]],
		   1/3 Min[\[Epsilon] Abs[qh]/Abs[qph], (\[Epsilon] Abs[qh]/Abs[qphs1])^(1/(1 + a)), (\[Epsilon] Abs[qh]/Abs[qphs2])^(1/(1 + 2 a))],
		Quiet[Abs[\[Tau]h/\[Tau]ph]], Quiet[Abs[\[Tau]ph/\[Tau]pph]]]]];

 
 PrintV[StringForm["Expansion parameters: a = `9`, f3s1 = `10`, f3s2 = `11`, fpph = `1`, fppph = `8`, \[Tau]ph = `2`, \[Tau]pph = `3`, qh = `7`,  qph = `4`, \[Lambda]ph = `5`, effective \[Epsilon] = `6`", fpph, \[Tau]ph, \[Tau]pph, qph, \[Lambda]ph, \[Epsilon]eff, qh, f3h, a, f3s1, f3s2], "Debug"];

 {{
 f'[\[Epsilon]eff] ==  \[Epsilon]eff fpph + 1/2 \[Epsilon]eff^2 (f3h + f3s1 \[Epsilon]eff^a + f3s2 \[Epsilon]eff^(2a)) + 1/6 (f3s1 a \[Epsilon]eff^(a + 2) + 2 f3s2 a \[Epsilon]eff^(2 a + 2)),
 \[Lambda][\[Epsilon]eff] == \[Lambda]h + \[Epsilon]eff (\[Lambda]ph + C\[Lambda] \[Epsilon]eff^a),
\[Lambda]'[\[Epsilon]eff] == C\[Lambda]*\[Epsilon]eff^a + a*C\[Lambda]*\[Epsilon]eff^a + \[Lambda]ph,
 q[\[Epsilon]eff] == qh + \[Epsilon]eff (qph + qphs1 \[Epsilon]eff^a + qphs2 \[Epsilon]eff^(2 a)),
 f[\[Epsilon]eff] ==  \[Epsilon]eff^2 fpph/2 + 1/6 \[Epsilon]eff^3(f3h + f3s1 \[Epsilon]eff^a + f3s2 \[Epsilon]eff^(2a)),
 \[Tau][\[Epsilon]eff] == \[Tau]h + \[Epsilon]eff \[Tau]ph + 1/2 \[Epsilon]eff^2 \[Tau]pph,
 \[Tau]'[\[Epsilon]eff] == \[Tau]ph + \[Epsilon]eff \[Tau]pph},
 (fpph > 0) && (qh + \[Epsilon]eff (qph + qphs1 \[Epsilon]eff^a + qphs2 \[Epsilon]eff^(2 a)) < 0) (*Return whether the initial conditions are viable*)
 },
 nt1crit,
 \[Epsilon]eff (*\[Epsilon] in A*)
}
]

](*vcontext*)
];


(*Constructs the full set of initial conditions at distance \[Epsilon] from the horizon by requiring regularity of the solution at the horizon*)
Options[VQCDBHInitialConditions] = {EpsilonMeaning->"fd" (*Sets the meaning of epsilon. fd means that it is the relative change in f', absolute means it is the actual offset from horizon in z*)
											,EpsilonLimit -> Infinity
}
VQCDBHInitialConditions[\[Tau]h_?NumericQ, \[Lambda]h_?NumericQ, nt_?NumericQ, fph_?NumericQ, {Vg_, Vf_, \[Kappa]_, w_}, \[Epsilon]_?NumericQ, opts:OptionsPattern[]] := 
 Module[{qh, \[Tau]ph, \[Tau]pph, fpph, qph, \[Lambda]ph, \[Epsilon]eff, qpph, \[Lambda]pph, fppph},

 {qh, \[Tau]ph, \[Tau]pph, fpph, qph, \[Lambda]ph, qpph, \[Lambda]pph, fppph} = FiniteTInitialConditions[\[Lambda]h, nt, fph, \[Tau]h, Vg, Vf, \[Kappa], w];

 
If[OptionValue[EpsilonMeaning] == "Absolute",
{{
 (*epsilon is directly the derivative itself*)
 {\[Lambda][\[Epsilon]] == \[Lambda]h + \[Epsilon] \[Lambda]ph + 1/2 \[Epsilon]^2 \[Lambda]pph,
 q[\[Epsilon]] == qh + \[Epsilon] qph + 1/2 \[Epsilon]^2qpph,
 \[Lambda]'[\[Epsilon]] == \[Lambda]ph + \[Lambda]pph \[Epsilon], 
 f[\[Epsilon]] ==  \[Epsilon] fph + 1/2 \[Epsilon]^2 fpph + 1/6 \[Epsilon]^3 fppph,
 f'[\[Epsilon]] ==  fph + \[Epsilon] fpph + 1/2 \[Epsilon]^2 fppph,
 \[Tau][\[Epsilon]] == \[Tau]h + \[Epsilon] \[Tau]ph + 1/2 \[Epsilon]^2 \[Tau]pph,
 \[Tau]'[\[Epsilon]] == \[Tau]ph + \[Epsilon] \[Tau]pph
 },
 (\[Lambda]ph < 0) && (fph > 0) && (qh + \[Epsilon]eff qph < 0) (*Return whether the initial conditions are viable*)
},
 \[Epsilon] (*This is returned here for consistency with the non-absolute version, even though it's obviously precisely what was input*)
}
,
{
 (*\[Epsilon] is the maximum allowed relative change in any of the fields during the first step.*)

 (*since fpph may become zero, we use the third derivative to estimate it's change also. Without this the automatic algorithm
  would increase \[Epsilon]eff to very large values near the point where fpph changes sign*)

 


 (*The mess is related to making Min cope with the case \[Tau]ph/\[Tau]pph = 0/0*)
 \[Epsilon]eff = Min[OptionValue[EpsilonLimit], \[Epsilon] Block[{Indeterminate = Infinity},Min[Abs[fph/fpph], Abs[\[Lambda]h/\[Lambda]ph], Abs[qh/qph], Abs[fpph/fppph], Quiet[Abs[\[Tau]h/\[Tau]ph]], Quiet[Abs[\[Tau]ph/\[Tau]pph]]]]];

 (*Debug printing*)
 Block[{$vcontext = "VQCDBHInitialConditions"},
 PrintV[StringForm["Expansion parameters: fpph = `1`, fppph = `7`, \[Tau]ph = `2`, \[Tau]pph = `3`, qh = `7`,  qph = `4`, \[Lambda]ph = `5`, effective \[Epsilon] = `6`", fpph, \[Tau]ph, \[Tau]pph, qph, \[Lambda]ph, \[Epsilon]eff, qh, fppph], "Debug"];
 ];

 {{
 f'[\[Epsilon]eff] ==  fph + \[Epsilon]eff fpph + 1/2 \[Epsilon]eff^2 fppph,
 \[Lambda][\[Epsilon]eff] == \[Lambda]h + \[Epsilon]eff \[Lambda]ph + 1/2 \[Epsilon]eff^2 \[Lambda]pph,
\[Lambda]'[\[Epsilon]eff] == \[Lambda]ph + \[Epsilon]eff \[Lambda]pph,
 q[\[Epsilon]eff] == qh + \[Epsilon]eff qph + 1/2 \[Epsilon]eff^2 qpph,
 f[\[Epsilon]eff] ==  \[Epsilon]eff fph + 1/2 \[Epsilon]eff^2fpph + 1/6 \[Epsilon]eff^3 fppph,
 \[Tau][\[Epsilon]eff] == \[Tau]h + \[Epsilon]eff \[Tau]ph + 1/2 \[Epsilon]eff^2 \[Tau]pph,
 \[Tau]'[\[Epsilon]eff] == \[Tau]ph + \[Epsilon]eff \[Tau]pph},
 (fph > 0) && (qh + \[Epsilon]eff qph < 0) (*Return whether the initial conditions are viable*)
 },
 \[Epsilon]eff (*Return the final epsilon*)
 }
]
]


(*ScaleSolution: Scales the solution UV asymptotics*)
ScaleSolution::lazhNotReal = "\[Lambda][zh] is not real, value is `1`";
ScaleSolution::badb0b1 = "Invalid value `1` for betacoefs, must a list of numeric values {b0, b1, ...}";
ScaleSolution::\[CapitalLambda]0errorlarge = "The estimated relative error `1` in \[CapitalLambda]0 is larger than the the set limit `2`";
Options[ScaleSolution] = {
	A0Improvement -> Differential, (*Options are Differential, {FiniteDifference, howmuchdifference} and None.*)
	AUV -> Automatic, \[CapitalLambda]relerrLimit-> 10^-3, WorkingPrecision -> Automatic, \[CapitalLambda]errorEstimate -> 9/10};
ScaleSolution[q_, f_, \[Lambda]_, \[Tau]d_, \[Tau]_, Amin_, Amax_, betacoefs_List, ell_, \[CapitalLambda]_, opts : OptionsPattern[]] := Module[{A0, q1, \[CapitalLambda]0relerror, b0, b1, f1, \[Lambda]1, \[Tau]d1, \[Tau]1, q2, f2, \[Lambda]2, \[Tau]d2, \[Tau]2, f3, \[Lambda]3,\[Tau]d3, \[Tau]3, fscale, \[CapitalLambda]0, zUV, bExp, \[CapitalLambda]scale, zh, betafitRMS, \[Lambda]fitRMS, zshift, \[Lambda]f1, \[Lambda]f2, A0improved, AUV0},
Block[{$vcontext = "ScaleSolution"},
(*The b0, b1 factors are analytically given*)
(*TODO: use more coefs from beta if they are available*)
PrintV["Starting to scale", "Debug"];
If[OptionValue[AUV] === Automatic, AUV0 = Amax, AUV0 = OptionValue[AUV]];
{b0, b1} = betacoefs;
If[!(NumericQ[b0] && NumericQ[b1]), Message[ScaleSolution::badb0b1, betacoefs]; Throw[Null]];
A0[AUV_] := AUV - Log[ell] - 1/\[Lambda][AUV]/b0 - b1/b0^2 Log[b0 \[Lambda][AUV]];

(*Eliminate the O(\[Lambda]) corrections, with the chosen method.*)
If[(Length[OptionValue[A0Improvement]] > 0) && OptionValue[A0Improvement][[1]] === FiniteDifference,
 A0improved[AUV_] = Module[{AUV2 = (1 - OptionValue[A0Improvement][[2]] AUV)  + OptionValue[A0Improvement][[2]]Amin}, (A0[AUV2]\[Lambda][AUV] - A0[AUV]\[Lambda][AUV2])/(\[Lambda][AUV] - \[Lambda][AUV2])],
 If[OptionValue[A0Improvement] === Differential, A0improved[AUV_] = A0[AUV] - A0'[AUV]\[Lambda][AUV]/ \[Lambda]'[AUV],
	A0improved = A0];
];

\[CapitalLambda]0 = Exp[A0improved[AUV0]];
\[CapitalLambda]0relerror = Exp[A0improved[AUV0] - A0improved[((1-OptionValue[\[CapitalLambda]errorEstimate])Amin+ OptionValue[\[CapitalLambda]errorEstimate]AUV0)]] - 1; (*Estimate the error in A0 by taking the exponential of the difference at AUV and 3/4 AUV*)
If[Abs[\[CapitalLambda]0relerror] > OptionValue[\[CapitalLambda]relerrLimit], Message[ScaleSolution::\[CapitalLambda]0errorlarge, \[CapitalLambda]0relerror, OptionValue[\[CapitalLambda]relerrLimit]]];
PrintV[StringForm["\[CapitalLambda]0 estimated relative error is `1`", \[CapitalLambda]0relerror], "Checks"];
PrintV[LogLinearPlot[{A0[Auv], A0[Amax], A0[((1-OptionValue[\[CapitalLambda]errorEstimate])Amin+ OptionValue[\[CapitalLambda]errorEstimate]Amax)]}, {Auv, Amax/2, Amax}, PlotStyle -> {Blue, {Red, Dashed}, {Red, Dashed}}, AxesLabel-> {"A", "A0"}], "Checks"];
(*scale f[0] to 1*)
fscale = 1/Sqrt[f[Amax]];
q1[A_] := q[A] * fscale;
f1[A_] := f[A]*fscale^2;
\[Lambda]1[A_] := \[Lambda][A];
\[Tau]d1[A_] := \[Tau]d[A];
\[Tau]1[A_] := \[Tau][A];
(*Scale \[CapitalLambda]*)
\[CapitalLambda]scale = \[CapitalLambda]0/\[CapitalLambda];

(*return the scaled functions and the scaling factors*)
{q1, f1, \[Lambda]1, \[Tau]d1, \[Tau]1, Amin, Amax, \[CapitalLambda]scale,
	fscale Sqrt[f'[Amin]] (*Multiplying by the Sqrt[f'\.1d] -term adjusts for possible
							non-zero settings of fpInitialValue.*)
	}
]
];
(*A helper function to get the scaled solutions directly from ndsolve output*)
ScaleSolution[{sols_, finished_}, betacoefs_, ell_, \[CapitalLambda]_, opts : OptionsPattern[]] := Module[{q0, f0, \[Tau]0, W0, \[Lambda]0, \[Tau]d0, Amin, Amax},
Clear[f, \[Tau], \[Lambda]];
Block[{$vcontext = "ScaleSolution"},
PrintV["Scaling q0...", "Debug"];
q0 = RemoveImaginaries[q /. sols];
PrintV["Scaling f0...", "Debug"];
f0 = RemoveImaginaries[f /. sols];
PrintV["Scaling \[Tau]0...", "Debug"];
\[Tau]0 = RemoveImaginaries[\[Tau] /. sols];
PrintV["Scaling \[Lambda]0...", "Debug"];
\[Lambda]0 = RemoveImaginaries[\[Lambda] /. sols];
PrintV["Scaling \[Tau]d0...", "Debug"];
\[Tau]d0 = RemoveImaginaries[\[Tau]' /. sols];
If[finished == False, 
(
Module[{q1, f1, \[Tau]1, W1, \[Lambda]1, \[Tau]d1},
q1 = RemoveImaginaries[q /. sols];
f1 = RemoveImaginaries[f /. sols];
\[Tau]1 = RemoveImaginaries[\[Tau] /. sols];
\[Lambda]1 = RemoveImaginaries[\[Lambda] /. sols];
\[Tau]d1 = RemoveImaginaries[\[Tau]' /. sols];
{Amin, Amax} = If[Head[q1] === InterpolatingFunction,
Re[First[InterpolatingFunctionDomain[q1]]],
{1, 2}
];
Return[{q1, f1, \[Lambda]1, \[Tau]d1, \[Tau]1, Amin, Amax, 1, 1, False}]
]
);
];
{Amin, Amax} = Re[First[InterpolatingFunctionDomain[q0]]];
Append[ScaleSolution[q0, f0, \[Lambda]0, \[Tau]d0, \[Tau]0, Amin, Amax, betacoefs, ell, \[CapitalLambda], opts], True]
]
]


Options[bCoefsFromPotential] = {bOrder -> 2};
bCoefsFromPotential[V_, opts : OptionsPattern[]] := Module[{bseries, Vseries, bcoef, sols, ford, beta, ellsqrd},

Block[{$vcontext = "bCoefsFromPotential"},

PrintV[StringForm["Potential: `1`", V[\[Lambda]]], "Debug"];

Vseries = Series[V[\[Lambda]], {\[Lambda], 0, OptionValue[bOrder] }, Assumptions -> {\[Lambda] > 0}];
ellsqrd = 12/Coefficient[Vseries, \[Lambda], 0];

beta[\[Lambda]_] := Sum[-bcoef[n-2] \[Lambda]^n, {n, 2, 2 + OptionValue[bOrder]}];
bseries = (12 / ellsqrd Series[Exp[-8/9 Integrate[beta[\[Lambda]]/\[Lambda]^2, {\[Lambda], 0, \[Lambda]}]](1-(beta[\[Lambda]] )^2/9/\[Lambda]^2), {\[Lambda], 0, OptionValue[bOrder]+ 2}]);

PrintV[StringForm["\.1dnumel is `1`", Sqrt[12/Coefficient[Vseries, \[Lambda], 0]]], "Debug"];

PrintV[StringForm["Computing beta coefs, bseries = `1`, Vseries = `2`", bseries, Vseries], "Debug"];

(*A function which gives the next coefficient*)
NextCoef[previousOrders_, order_] := Module[{eq},
PrintV[StringForm["b-series coefficient at order `1` is `2`", order, Coefficient[bseries, \[Lambda]^(order + 1)]], "Debug"];
eq = (Coefficient[bseries, \[Lambda]^(order+1)]/.previousOrders) == (Coefficient[Vseries, \[Lambda]^(order+1)]/.previousOrders);
PrintV[StringForm["Solving `1`, solution is `2`", eq, Solve[eq, bcoef[order]]], "Debug"];
{Join[previousOrders, First[Solve[eq, bcoef[order]]]], order + 1}
];

{sols, ford} = Nest[NextCoef[#[[1]], #[[2]]]&, {{}, 0}, OptionValue[bOrder]];

PrintV[StringForm["Final solution list is `1`, replacing in `2`", sols, Table[bcoef[order], {n, 0, OptionValue[bOrder]}]], "Debug"];

{Table[bcoef[n], {n, 0,  OptionValue[bOrder]-1}]/.sols, Sqrt[ellsqrd]}
]

]

(*A shortcut for using the potential list*)
bCoefsFromPotential[{Vg_, Vf_, \[Kappa]_, \[Omega]_}, opts : OptionsPattern[]] := bCoefsFromPotential[Vg[#] - Vf[#, 0]&, opts];


GammaFromKappa[\[Kappa]_, b0_, Vf_, ell_] := Module[{h1},
(*Assumes that Vf is of the form Vf=exp(-a(\[Lambda])\[Tau]^2)Vf0(\[Lambda])*)
h1 = 3/2/ell^2SeriesCoefficient[\[Kappa][\[Lambda]]/-SeriesCoefficient[Log[Vf[\[Lambda], \[Tau]]], {\[Tau], 0, 2}], {\[Lambda], 0, 1}];
PrintV[StringForm["h1 = `1`", h1], "Debug"];
(3(8/9 + h1/b0)/2)
];

GammaFromKappa[{Vg_, Vf_, \[Kappa]_, \[Omega]_}, b0_, ell_] := GammaFromKappa[\[Kappa], b0, Vf, ell];


Options[SolverCore] = Join[{ARange -> 160 (*The range of A-values to solve for. Should be large enough to guarantee that the boundary asymptotic is reached, but too large values may cause a too large initial step*)}
					, Options[NDSolve]];
SolverCore[{initconds_, initvalid_}, nt_, Ah_, pots_List, opts : OptionsPattern[]] := Module[{invalidsol, sols, WentImaginary = False},
Block[{$vcontext = "SolverCore"},

(*Build the invalid return list, in case we need it*)
invalidsol := {{f -> 1, \[Lambda] -> Function[A, \[Lambda]h], \[Tau] -> Function[A, \[Tau]h], \[Tau]'-> 1, q -> 1},False};

PrintV[StringForm["Initial conditions are `1`\.1d", initconds], "Debug"];

(*Check the initial conditions*)
(*Real*)
If[!FreeQ[initconds, _Complex],
(PrintV["The initial conditions are not real", "Checks"];
Return[invalidsol];)];
(*positive \[Lambda]ph*)
If[!initvalid,
(
 PrintV[StringForm["The initial conditions `1` are not valid at \[Lambda]h = `2`, \[Tau]h = `3`, ntilde = `4`\.08", initconds, \[Lambda]h, \[Tau]h, nt], "Checks"];
 Return[invalidsol];)
];


PrintV[StringForm["Initial conditions are `1`", initconds], "Debug"];
Check[ sols = Quiet[NDSolve[Evaluate[Join[VQCDEquationsOfMotion[pots, nt], initconds]], {f, \[Lambda], \[Tau], q, \[Tau]'}, {A, Ah, OptionValue[ARange]},
 Evaluate[FilterRules[{opts}, Options[NDSolve]]], AccuracyGoal -> 60, Compiled-> True,
Method -> {"EventLocator"(*, Method -> "Adams"*), "Event"-> ((Im[\[Lambda][A]] != 0) || (Im[\[Tau][A]] != 0) || (Im[q[A]] != 0) || (Im[f[A]] != 0)), 
  "EventAction":> (WentImaginary = True;
	PrintV[StringForm["\[Lambda](A) = `2` is complex at A = `1`, \[Lambda]h = `3`, \[Tau]h = `4`", A, \[Lambda][A], \[Lambda]h, \[Tau]h], "Checks"]; Throw[Null, "StopIntegration"]), Sequence @@ If[$VersionNumber >= 9., {"EquationSimplification" -> {Automatic, "TimeConstraint" -> 120}}, {}]}], NDSolve::precw],
  Return[invalidsol](*, NDSolCheckve::mxst*)];


{Last[sols], !WentImaginary(*Return whether the computation finished or not*)}
]
] 


(*SolveFiniteTTachyons: Constructs the initial value problem and numerically solves it.*)
Options[SolveVQCDBH] = Join[{
HorizonEpsilon -> 10^(-8), (*The distance from the horizon at which to set the initial conditions*)
fpInitialValue -> 10^-3, (*Initial value for f derivative at the horizon*)
\[CapitalLambda] -> 1 (*Lambda to scale to*)
}, Options[SolverCore], Options[VQCDBHInitialConditions]]; (*Allow any options for NDSolve or TachyonInitialConditionsAtHorizon to be passed on*)
SolveVQCDBH[\[Lambda]h_?NumericQ, \[Tau]h_?NumericQ, nt_?NumericQ, {Vg_, Vf_, \[Kappa]_, \[Omega]_}, opts : OptionsPattern[]] := Module[{initconds, initvalid, sols, Vgh, Vfh, invalidsol, \[Epsilon]},

Block[{$vcontext = "SolveVQCDBH"},

(*Build the invalid return list, in case we need it*)
invalidsol := {{f -> 1, \[Lambda] -> Function[A, \[Lambda]h], \[Tau] -> Function[A, \[Tau]h], \[Tau]'-> 1, q -> 1},False};

(*Check that this is a reasonable starting point at all*)
Vgh = Vg[\[Lambda]h];
Vfh = Vf[\[Lambda]h, \[Tau]h];
If[!Element[Vgh, Reals] || !Element[Vfh, Reals] || TrueQ[Vg[\[Lambda]h] < Vf[\[Lambda]h, \[Tau]h]],
(PrintV[StringForm["The potential is not real."],"Checks"];
Return[invalidsol];)
];

(*Solve*)
{{initconds, initvalid}, \[Epsilon]} = VQCDBHInitialConditions[\[Tau]h, \[Lambda]h, nt, OptionValue[fpInitialValue], {Vg, Vf, \[Kappa], \[Omega]}, OptionValue[HorizonEpsilon], Evaluate[FilterRules[{opts}, Options[VQCDBHInitialConditions]]]];

SolverCore[{initconds, initvalid}, nt, \[Epsilon], {Vg, Vf, \[Kappa], \[Omega]}, Evaluate[FilterRules[{opts}, Options[SolverCore]]]]

]
]

(*SolveAndScaleVQCDBH: A convenience wrapper for solving the equations and scaling them*)
Options[SolveAndScaleVQCDBH] = Join[Options[SolveVQCDBH], Options[ScaleSolution], Options[bCoefsFromPotential]];
SolveAndScaleVQCDBH[\[Lambda]h_?NumericQ, \[Tau]h_?NumericQ, nt_?NumericQ, pots_List, opts : OptionsPattern[]] := Module[{scaledsol, bcoefs},
Block[{$vcontext = "SolveAndScaleFiniteTTachyons"},
(*Find the b-expansion coefficients*)
bcoefs = bCoefsFromPotential[pots, Evaluate[FilterRules[{opts}, Options[bCoefsFromPotential]]]];
(*Print[StringForm["\.1dSolving FiniteTTachyons at la = `1`, \[Tau]h = `2`", \[Lambda]h, \[Tau]h]];*)
(*Compute the solution*)
scaledsol = ScaleSolution[
	SolveVQCDBH[\[Lambda]h, \[Tau]h, nt, pots, Evaluate[FilterRules[{opts}, Options[SolveVQCDBH]]]],
	First[bcoefs], (*The first element is the list of b0, b1, ...*)
	Last[bcoefs], (*The last element is ell*)
	OptionValue[\[CapitalLambda]], Evaluate[FilterRules[{opts},
	Options[ScaleSolution]]]];
	
	(*return the scaled solution followed by the b-coefficients*)
	Join[scaledsol, bcoefs]
]
];


Options[SolveVQCDExtremalBH] = Join[{
 HorizonEpsilon -> 10^(-8), qHorizon -> -1},
 Options[SolverCore],
 Options[VQCDExtremalBHInitialConditions],
 Options[ScaleSolution]
]
SolveVQCDExtremalBH[C\[Lambda]_?NumericQ, pots_List, opts : OptionsPattern[]] := Module[{initconds, initvalid, nt1crit, Ah},
Block[{$cvontext = "SolveExtremalTachyons"},
 
(*Solve*)
{{initconds, initvalid}, nt1crit, Ah} = VQCDExtremalBHInitialConditions[OptionValue[qHorizon], C\[Lambda], pots, OptionValue[HorizonEpsilon], Evaluate[FilterRules[{opts}, Options[VQCDExtremalBHInitialConditions]]]];

 SolverCore[{initconds, initvalid}, nt1crit, Ah, pots, Evaluate[FilterRules[{opts}, Options[SolverCore]]]]
]
];

(*SolveAndScaleFiniteTTachyons: A convenience wrapper for solving the equations and scaling them*)
Options[SolveAndScaleVQCDExtremalBH] = Join[Options[SolveVQCDExtremalBH], Options[ScaleSolution], Options[bCoefsFromPotential]];
SolveAndScaleVQCDExtremalBH[C\[Lambda]_?NumericQ, pots_List, opts : OptionsPattern[]] := Module[{scaledsol, bcoefs},
Block[{$vcontext = "SolveAndScaleFiniteTTachyons"},
(*Find the b-expansion coefficients*)
bcoefs = bCoefsFromPotential[pots, Evaluate[FilterRules[{opts}, Options[bCoefsFromPotential]]]];
(*Print[StringForm["\.1dSolving FiniteTTachyons at la = `1`, \[Tau]h = `2`", \[Lambda]h, \[Tau]h]];*)
(*Compute the solution*)
scaledsol = ScaleSolution[
	SolveVQCDExtremalBH[C\[Lambda], pots, Evaluate[FilterRules[{opts}, Options[SolveVQCDExtremalBH]]]],
	First[bcoefs], (*The first element is the list of b0, b1, ...*)
	Last[bcoefs], (*The last element is ell*)
	1, Evaluate[FilterRules[{opts},
	Options[ScaleSolution]]]];
	
	(*return the scaled solution followed by the b-coefficients*)
	Join[scaledsol, bcoefs]
]
];


VEffective[{Vg_, Vf_, \[Kappa]_, \[Omega]_}] := Function[{\[Lambda]h, nt, \[Tau]h}, Vg[\[Lambda]h] - Vf[\[Lambda]h, \[Tau]h]Sqrt[1 + nt^2/(\[Omega][\[Lambda]h] Vf[\[Lambda]h, \[Tau]h])^2]]


ntCritical[pots_List, \[Tau]h_ : 0] := Module[{Vg, Vf, \[Kappa], \[Omega], ntfun}, 
	ntfun[\[Lambda]h_] = (nt /. Last[Solve[VEffective[{Vg, Vf, \[Kappa], \[Omega]}][\[Lambda]h, nt, \[Tau]h] == 0, nt]]) 
					/. {Vg -> pots[[1]], Vf -> pots[[2]], \[Kappa] -> pots[[3]], \[Omega] -> pots[[4]]};
	ntfun
]


DefineOldPotI[xf_, \[Lambda]0v_ : 1,
\[Lambda]s_ :  1(*8 Pi^2*) (*Choice of units*),
(*beta coefs*)
b0_ : Function[xf, 1/3 (11-2 xf)], b1_ : Function[xf, 1/6 (34-13 xf)]] := Module[{pots, Vg, Vf, \[Kappa], Vf0, Vg0, V0, W0,\[Kappa]0, xfd},

Vg0[\[Lambda]_] := Vuv[0]+Vuv[1]\[Lambda] +Vuv[2] (\[Lambda])^2/(1+\[Lambda]/\[Lambda]0)^(2/3) Sqrt[1+Log[\[Lambda]/\[Lambda]0+1]]/. \[Lambda]0 -> \[Lambda]0v;
Vf0[\[Lambda]_] := Wuv[0]+Wuv[1]\[Lambda]+Wuv[2]\[Lambda]^2/. \[Lambda]0 -> \[Lambda]0v;
\[Kappa]0[\[Lambda]_] := 1/(1+\[Gamma]fix \[Lambda])^(4/3);
V0 = 12;
W0 = 12/11;
pots = (PotFromAnsatzAndBeta[Vg0, Vf0, \[Kappa]0, \[Kappa]0, 1&, V0, W0, xfd, b0[xfd], b1[xfd], \[Lambda]s]) /. xfd -> xf;
 Vg[\[Lambda]v_] = pots[[1]] /. \[Lambda] -> \[Lambda]v;
 Vf[\[Lambda]v_, \[Tau]v_] = pots[[2]] /. {\[Lambda] -> \[Lambda]v, \[Tau] -> \[Tau]v};
 \[Kappa][\[Lambda]v_] = pots[[3]] /. \[Lambda] -> \[Lambda]v;
 {Vg, Vf, \[Kappa], \[Kappa]}
]


Options[DefineJKPotIMod] = {W0 -> (12/#(1 - 1/(1 + 7 # /4)^(2/3))&),
V0 -> (12&),
b0 -> (1/3 (11-2 #)&),
b1 ->  (1/6 (34-13 #)&),
\[Lambda]0v -> 1,
\[Lambda]s -> 1};
DefineJKPotIMod[xf_(*, \[Lambda]0v_ : 1(*(8 Pi^2)*),
\[Lambda]s_ :  1(*8 Pi^2*) (*Choice of units*),
(*beta coefs*)
b0_ : Function[xf, 1/3 (11-2 xf)], b1_ : Function[xf, 1/6 (34-13 xf)]*), opts : OptionsPattern[]] := Module[{pots, Vg, Vf, \[Kappa], Vf0, Vg0, V0f, W0f,\[Kappa]0, xfd},
Vg0[\[Lambda]_] := Vuv[0]+Vuv[1]\[Lambda] +Vuv[2] (\[Lambda])^2/(1+\[Lambda]/\[Lambda]0)^(2/3) Sqrt[1+Log[\[Lambda]/\[Lambda]0+1]]/. \[Lambda]0 -> OptionValue[\[Lambda]0v];
Vf0[\[Lambda]_] := Wuv[0]+Wuv[1]\[Lambda]+Wuv[2] (\[Lambda])^2/(1+\[Lambda]/\[Lambda]0)^(2/3) Sqrt[1+Log[\[Lambda]/\[Lambda]0+1]] /. \[Lambda]0 -> OptionValue[\[Lambda]0v];
\[Kappa]0[\[Lambda]_] := 1/(1+\[Gamma]fix \[Lambda])^(4/3);
V0f = OptionValue[V0][xfd];
W0f = OptionValue[W0][xfd];
pots = (PotFromAnsatzAndBeta[Vg0, Vf0, \[Kappa]0, \[Kappa]0, 1&, V0f, W0f, xfd, OptionValue[b0][xfd], OptionValue[b1][xfd], OptionValue[\[Lambda]s]]) /. xfd -> xf;
 Vg[\[Lambda]v_] = pots[[1]] /. \[Lambda] -> \[Lambda]v;
 Vf[\[Lambda]v_, \[Tau]v_] = pots[[2]] /. {\[Lambda] -> \[Lambda]v, \[Tau] -> \[Tau]v};
 \[Kappa][\[Lambda]v_] = pots[[3]] /. \[Lambda] -> \[Lambda]v;
 {Vg, Vf, \[Kappa], \[Kappa]}
]


Options[DefineJKPotI] = {W0 -> (12/#(1 - 1/(1 + 7 # /4)^(2/3))&),
V0 -> (12&),
b0 -> (1/3 (11-2 #)&),
b1 ->  (1/6 (34-13 #)&),
\[Lambda]0v -> 1,
\[Lambda]s -> 1};
DefineJKPotI[xf_(*, \[Lambda]0v_ : 1(*(8 Pi^2)*),
\[Lambda]s_ :  1(*8 Pi^2*) (*Choice of units*),
(*beta coefs*)
b0_ : Function[xf, 1/3 (11-2 xf)], b1_ : Function[xf, 1/6 (34-13 xf)]*), opts : OptionsPattern[]] := Module[{pots, Vg, Vf, \[Kappa], Vf0, Vg0, V0f, W0f,\[Kappa]0, xfd},
Vg0[\[Lambda]_] := Vuv[0]+Vuv[1]\[Lambda] +Vuv[2] (\[Lambda])^2/(1+\[Lambda]/\[Lambda]0)^(2/3) Sqrt[1+Log[\[Lambda]/\[Lambda]0+1]]/. \[Lambda]0 -> OptionValue[\[Lambda]0v];
Vf0[\[Lambda]_] := Wuv[0]+Wuv[1]\[Lambda]+Wuv[2]\[Lambda]^2/. \[Lambda]0 -> OptionValue[\[Lambda]0v];
\[Kappa]0[\[Lambda]_] := 1/(1+\[Gamma]fix \[Lambda])^(4/3);
V0f = OptionValue[V0][xfd];
W0f = OptionValue[W0][xfd];
pots = (PotFromAnsatzAndBeta[Vg0, Vf0, \[Kappa]0, \[Kappa]0, 1&, V0f, W0f, xfd, OptionValue[b0][xfd], OptionValue[b1][xfd], OptionValue[\[Lambda]s]]) /. xfd -> xf;
Vg[\[Lambda]v_] = pots[[1]] /. \[Lambda] -> \[Lambda]v;
Vf[\[Lambda]v_, \[Tau]v_] = pots[[2]] /. {\[Lambda] -> \[Lambda]v, \[Tau] -> \[Tau]v};
\[Kappa][\[Lambda]v_] = pots[[3]] /. \[Lambda] -> \[Lambda]v;
{Vg, Vf, \[Kappa], \[Kappa]}
]


Options[DefineJKPotIIMod] = {W0 -> (12/#(1 - 1/(1 + 7 # /4)^(2/3))&),
V0 -> (12&),
b0 -> (1/3 (11-2 #)&),
b1 ->  (1/6 (34-13 #)&),
\[Lambda]0v -> 1,
\[Lambda]s -> 1};
DefineJKPotIIMod[xf_(*,
W0 : Function[xf, (12/xf(1 - 1/(1 + 7 xf /4)^(2/3)))],
\[Lambda]0v_ : 1(*(8 Pi^2)*),
\[Lambda]s_ :  1(*8 Pi^2*) (*Choice of units*), 
(*beta coefs*)
b0_ : Function[xf, 1/3 (11-2 xf)], b1_ : Function[xf, 1/6 (34-13 xf)]*), opts:OptionsPattern[]] := Module[{pots, Vg, Vf, \[Kappa], Vf0, Vg0, V0f, W0f,\[Kappa]0, xfd, a0},

Vg0[\[Lambda]_] := Vuv[0]+Vuv[1]\[Lambda] +Vuv[2] (\[Lambda])^2/(1+\[Lambda]/\[Lambda]0)^(2/3) Sqrt[1+Log[\[Lambda]/\[Lambda]0+1]]/. \[Lambda]0 -> OptionValue[\[Lambda]0v];
Vf0[\[Lambda]_] := Wuv[0]+Wuv[1]\[Lambda]+Wuv[2](\[Lambda])^2/(1+\[Lambda]/\[Lambda]0)^(2/3) Sqrt[1+Log[\[Lambda]/\[Lambda]0+1]]/. \[Lambda]0 -> OptionValue[\[Lambda]0v];
\[Kappa]0[\[Lambda]_] := 1/(1 + \[Lambda]/\[Lambda]0)^(4/3) /. \[Lambda]0 -> OptionValue[\[Lambda]0v];
V0f = OptionValue[V0][xfd];
W0f = OptionValue[W0][xfd];
a0[\[Lambda]_] := (1-\[Gamma]fix \[Lambda] + \[Lambda]^2/\[Lambda]0^2)/(1+ \[Lambda]/\[Lambda]0)^(4/3) /. \[Lambda]0 -> OptionValue[\[Lambda]0v];
pots = (PotFromAnsatzAndBeta[Vg0, Vf0, \[Kappa]0, \[Kappa]0, a0, V0f, W0f, xfd, OptionValue[b0][xfd], OptionValue[b1][xfd], OptionValue[\[Lambda]s]]) /. xfd -> xf;
 Vg[\[Lambda]v_] = pots[[1]] /. \[Lambda] -> \[Lambda]v;
 Vf[\[Lambda]v_, \[Tau]v_] = pots[[2]] /. {\[Lambda] -> \[Lambda]v, \[Tau] -> \[Tau]v};
 \[Kappa][\[Lambda]v_] = pots[[3]] /. \[Lambda] -> \[Lambda]v;
 {Vg, Vf, \[Kappa], \[Kappa]}
]


Options[DefineJKPotII] = {W0 -> (12/#(1 - 1/(1 + 7 # /4)^(2/3))&),
V0 -> (12&),
b0 -> (1/3 (11-2 #)&),
b1 ->  (1/6 (34-13 #)&),
\[Lambda]0v -> 1,
\[Lambda]s -> 1};
DefineJKPotII[xf_(*,
W0 : Function[xf, (12/xf(1 - 1/(1 + 7 xf /4)^(2/3)))],
\[Lambda]0v_ : 1(*(8 Pi^2)*),
\[Lambda]s_ :  1(*8 Pi^2*) (*Choice of units*), 
(*beta coefs*)
b0_ : Function[xf, 1/3 (11-2 xf)], b1_ : Function[xf, 1/6 (34-13 xf)]*), opts:OptionsPattern[]] := Module[{pots, Vg, Vf, \[Kappa], Vf0, Vg0, V0f, W0f,\[Kappa]0, xfd, a0},

Vg0[\[Lambda]_] := Vuv[0]+Vuv[1]\[Lambda] +Vuv[2] (\[Lambda])^2/(1+\[Lambda]/\[Lambda]0)^(2/3) Sqrt[1+Log[\[Lambda]/\[Lambda]0+1]]/. \[Lambda]0 -> OptionValue[\[Lambda]0v];
Vf0[\[Lambda]_] := Wuv[0]+Wuv[1]\[Lambda]+Wuv[2]\[Lambda]^2/. \[Lambda]0 -> OptionValue[\[Lambda]0v];
\[Kappa]0[\[Lambda]_] := 1/(1 + \[Lambda]/\[Lambda]0)^(4/3) /. \[Lambda]0 -> OptionValue[\[Lambda]0v];
V0f = OptionValue[V0][xfd];
W0f = OptionValue[W0][xfd];
a0[\[Lambda]_] := (1-\[Gamma]fix \[Lambda] + \[Lambda]^2/\[Lambda]0^2)/(1+ \[Lambda]/\[Lambda]0)^(4/3) /. \[Lambda]0 -> OptionValue[\[Lambda]0v];
pots = (PotFromAnsatzAndBeta[Vg0, Vf0, \[Kappa]0, \[Kappa]0, a0, V0f, W0f, xfd, OptionValue[b0][xfd], OptionValue[b1][xfd], OptionValue[\[Lambda]s]]) /. xfd -> xf;
Vg[\[Lambda]v_] = pots[[1]] /. \[Lambda] -> \[Lambda]v;
Vf[\[Lambda]v_, \[Tau]v_] = pots[[2]] /. {\[Lambda] -> \[Lambda]v, \[Tau] -> \[Tau]v};
\[Kappa][\[Lambda]v_] = pots[[3]] /. \[Lambda] -> \[Lambda]v;
{Vg, Vf, \[Kappa], \[Kappa]}
]


DefineJKPotIII[xf_,  \[Lambda]0v_ : 1(*(8 Pi^2)*),
\[Lambda]s_ :  1(*8 Pi^2*) (*Choice of units*),
(*beta coefs*)
b0_ : Function[xf, 1/3 (11-2 xf)], b1_ : Function[xf, 1/6 (34-13 xf)]
] := Module[{pots, Vg, Vf, \[Kappa], Vf0, Vg0, V0, W0,\[Kappa]0, xfd},

Vg0[\[Lambda]_] := Vuv[0]+Vuv[1]\[Lambda] +Vuv[2] (\[Lambda])^2/(1+\[Lambda]/\[Lambda]0)^(2/3) Sqrt[1+Log[\[Lambda]/\[Lambda]0+1]]/. \[Lambda]0 -> \[Lambda]0v;
Vf0[\[Lambda]_] := Wuv[0]+Wuv[1]\[Lambda]+Wuv[2] (\[Lambda])^2/(1+\[Lambda]/\[Lambda]0)^(2/3) Sqrt[1+Log[\[Lambda]/\[Lambda]0+1]] /. \[Lambda]0 -> \[Lambda]0v;
\[Kappa]0[\[Lambda]_] := 1/(1+\[Gamma]fix \[Lambda])^(4/3);
V0 = 12;
W0 = 12/xfd(1 - 1/(1 + 7 xfd /4)^(2/3));
pots = (PotFromAnsatzAndBeta[Vg0, Vf0, \[Kappa]0, \[Kappa]0, 1&, V0, W0, xfd, b0[xfd], b1[xfd], \[Lambda]s]) /. xfd -> xf;
 Vg[\[Lambda]v_] = pots[[1]] /. \[Lambda] -> \[Lambda]v;
 Vf[\[Lambda]v_, \[Tau]v_] = pots[[2]] /. {\[Lambda] -> \[Lambda]v, \[Tau] -> \[Tau]v};
 \[Kappa][\[Lambda]v_] = pots[[3]] /. \[Lambda] -> \[Lambda]v;
{Vg, Vf, \[Kappa], \[Kappa]}
]


Options[DefineJKPotIKappaMod] = {W0 -> (12/#(1 - 1/(1 + 7 # /4)^(2/3))&),
V0 -> (12&),
b0 -> (1/3 (11-2 #)&),
b1 ->  (1/6 (34-13 #)&),
\[Lambda]0v -> 1,
\[Lambda]s -> 1,
\[Mu] -> -1/2};
DefineJKPotIKappaMod[xf_(*, \[Lambda]0v_ : 1(*(8 Pi^2)*),
\[Lambda]s_ :  1(*8 Pi^2*) (*Choice of units*),
(*beta coefs*)
b0_ : Function[xf, 1/3 (11-2 xf)], b1_ : Function[xf, 1/6 (34-13 xf)]*), opts : OptionsPattern[]] := Module[{Vg, Vf, \[Kappa], pots, Vf0, Vg0, V0f, W0f,\[Kappa]0, xfd},
Vg0[\[Lambda]_] := Vuv[0]+Vuv[1]\[Lambda] +Vuv[2] (\[Lambda])^2/(1+\[Lambda]/\[Lambda]0)^(2/3) Sqrt[1+Log[\[Lambda]/\[Lambda]0+1]]/. \[Lambda]0 -> OptionValue[\[Lambda]0v];
Vf0[\[Lambda]_] := Wuv[0]+Wuv[1]\[Lambda]+Wuv[2]\[Lambda]^2/. \[Lambda]0 -> OptionValue[\[Lambda]0v];
\[Kappa]0[\[Lambda]_] := (1+ Log[1 + \[Lambda]])^OptionValue[\[Mu]]/(1+\[Gamma]fix \[Lambda])^(4/3);
V0f = OptionValue[V0][xfd];
W0f = OptionValue[W0][xfd];
pots = (PotFromAnsatzAndBeta[Vg0, Vf0, \[Kappa]0, \[Kappa]0, 1&, V0f, W0f, xfd, OptionValue[b0][xfd], OptionValue[b1][xfd], OptionValue[\[Lambda]s]]) /. xfd -> xf;
Vg[\[Lambda]v_] = pots[[1]] /. \[Lambda] -> \[Lambda]v;
Vf[\[Lambda]v_, \[Tau]v_] = pots[[2]] /. {\[Lambda] -> \[Lambda]v, \[Tau] -> \[Tau]v};
\[Kappa][\[Lambda]v_] = pots[[3]] /. \[Lambda] -> \[Lambda]v;
 {Vg, Vf, \[Kappa], \[Kappa]}
]


Options[DefineAIJKOptPot] = {W0 -> (3/11&)(*(12/#(1 - 1/(1 + 7 # /4)^(2/3))&)*),
V0 -> (12&),
b0 -> (1/3 (11-2 #)&),
b1 ->  (1/6 (34-13 #)&),
\[Lambda]0v -> 1, (*Choice of units, 1 is consistent with the beta given here*)
\[Lambda]s -> 1,
\[Mu] -> 1/2,
\[Mu]\[Omega] -> 1,
D0 -> 200};
DefineAIJKOptPot[xf_(*, \[Lambda]0v_ : 1(*(8 Pi^2)*),
\[Lambda]s_ :  1(*8 Pi^2*) (*Choice of units*),
(*beta coefs*)
b0_ : Function[xf, 1/3 (11-2 xf)], b1_ : Function[xf, 1/6 (34-13 xf)]*), opts : OptionsPattern[]] := Module[{Vg, Vf, \[Kappa], \[Omega], pots, Vf0, Vg0, V0f, W0f,\[Kappa]0, \[Omega]0, xfd},
Vg0[\[Lambda]_] := Vuv[0]+Vuv[1]\[Lambda] +Vuv[2] (\[Lambda])^2/(1+\[Lambda]/\[Lambda]0)^(2/3) Sqrt[1+Log[\[Lambda]/\[Lambda]0+1]]/. \[Lambda]0 -> OptionValue[\[Lambda]0v];
Vf0[\[Lambda]_] := Wuv[0]+Wuv[1]\[Lambda]+Wuv[2]\[Lambda]^2/. \[Lambda]0 -> OptionValue[\[Lambda]0v];
\[Kappa]0[\[Lambda]_] := (1+ 1/OptionValue[D0]Log[1 + (\[Lambda]/\[Lambda]0)^2])^OptionValue[\[Mu]]/(1+\[Gamma]fix \[Lambda])^(4/3) /. \[Lambda]0 -> OptionValue[\[Lambda]0v];
\[Omega]0[\[Lambda]_] := (1+ 1/OptionValue[D0]Log[1 + (\[Lambda]/\[Lambda]0)^2])^OptionValue[\[Mu]\[Omega]]/(1+\[Gamma]fix \[Lambda])^(4/3) /. \[Lambda]0 -> OptionValue[\[Lambda]0v];
V0f = OptionValue[V0][xfd];
W0f = OptionValue[W0][xfd];
pots = (PotFromAnsatzAndBeta[Vg0, Vf0, \[Kappa]0, \[Omega]0, 1&, V0f, W0f, xfd, OptionValue[b0][xfd], OptionValue[b1][xfd], OptionValue[\[Lambda]s]]) /. xfd -> xf;
Vg[\[Lambda]v_] = pots[[1]] /. \[Lambda] -> \[Lambda]v;
Vf[\[Lambda]v_, \[Tau]v_] = pots[[2]] /. {\[Lambda] -> \[Lambda]v, \[Tau] -> \[Tau]v};
\[Kappa][\[Lambda]v_] = pots[[3]] /. \[Lambda] -> \[Lambda]v;
\[Omega][\[Lambda]v_] = pots[[4]] /. \[Lambda] -> \[Lambda]v;
 {Vg, Vf, \[Kappa], \[Omega]}
]


PotFromAnsatzAndBeta[Vg0_, Vf0_, \[Kappa]0_, \[Omega]0_, a0_, V0_, W0_, xf_, b0_, b1_, \[Lambda]s_] := Module[{ Vf, Vg, \[Kappa], \[Omega], \[Kappa]p0v, coefs, a, \[Gamma]0,\[Gamma]1,h1sol, h1, \[Gamma]fixv},

(*a = 3/2 (V0 - xf W0)/12 a0[\[Lambda]];*)
\[Gamma]0 = 3/2 / \[Lambda]s;
\[Gamma]1 = (203 - 20 xf)/(12 (2 \[Lambda]s)^2);
coefs = PotCoefsFromBeta[Vg0, Vf0, W0, V0, xf, \[Lambda]s, b0, b1];
PrintV[StringForm["coefs are `1`", coefs], "Debug"];

(*Compute h1, the O(\[Lambda]) coef of \[Kappa]/a*)
h1 = Coefficient[Series[\[Kappa]0[\[Lambda]]/a0[\[Lambda]], {\[Lambda], 0, 1}], \[Lambda]];
PrintV[StringForm["h1 is `1`", h1]];
(*\[Kappa]sol = Solve[3 ( 8/9 + \[Kappa]0'[0]/ (b0/\[Lambda]s))/2 == -\[Gamma]0/(b0/\[Lambda]s), \[Kappa]p0];*)
h1sol = Solve[3 ( 8/9 + h1/ (b0/\[Lambda]s))/2 == -\[Gamma]0/(b0/\[Lambda]s), \[Gamma]fix];
\[Gamma]fixv = \[Gamma]fix/. First[h1sol];
PrintV[StringForm["gamma fixing factor is `1`", \[Gamma]fixv], "Debug"];

{\[Kappa][\[Lambda]_], a[\[Lambda]_], \[Omega][\[Lambda]_]} = Evaluate[{\[Kappa]0[\[Lambda]], a0[\[Lambda]] 3/2 (V0 - xf W0)/12, \[Omega]0[\[Lambda]]} /. \[Gamma]fix -> \[Gamma]fixv];
Vg[\[Lambda]_] = Evaluate[Vg0[\[Lambda]]/.coefs];
Vf[\[Lambda]_, \[Tau]_] = Exp[- a[\[Lambda]] \[Tau]^2] Simplify[Evaluate[xf Vf0[\[Lambda]] /. coefs]];
{Vg[\[Lambda]], Vf[\[Lambda], \[Tau]], \[Kappa][\[Lambda]], \[Omega][\[Lambda]]}
]


PotCoefsFromBeta[Vg0_, Vf0_, W0_, V0_, xf_, \[Lambda]s_, b0_, b1_]:=Module[{Vseries, ellsqrd, bseries, \[Beta]YM, sols, ford},

\[Beta]YM[\[Lambda]_] :=\[Lambda]s(- b0(\[Lambda]/\[Lambda]s)^2 - b1 (\[Lambda]/\[Lambda]s)^3);

Vseries = Series[(Vg0[\[Lambda]]- xf Vf0[\[Lambda]]), {\[Lambda], 0,2 }];
ellsqrd = 12/Coefficient[Vseries, \[Lambda], 0];
PrintV[StringForm["ell^2 is `1`", ellsqrd], "Debug"];
bseries = (12 / ellsqrd Series[Exp[-8/9 Integrate[\[Beta]YM[\[Lambda]]/\[Lambda]^2, {\[Lambda], 0, \[Lambda]}]](1-(\[Beta]YM[\[Lambda]] )^2/9/\[Lambda]^2), {\[Lambda], 0, 2}]);

(*A function which gives the next pair of coefficients*)
NextCoef[previousOrders_, order_] := Module[{eq, sol1, Vn, Wn, sol2},
PrintV[StringForm["b-series coefficient at order `1` is `2`", order, Coefficient[bseries, \[Lambda]^(order + 1)]], "Debug"];
eq = (Coefficient[bseries, \[Lambda]^(order+1)]/.previousOrders) == (Coefficient[Vseries, \[Lambda]^(order+1)]/.previousOrders);
sol1 = Solve[eq/.xf-> 0, Vuv[order+1]];
PrintV[StringForm["Solving `1`, solution is `2`", eq, sol1], "Debug"];
Vn = First[Vuv[order+1]/.sol1];
PrintV[StringForm["V`2` is `1`", Vn, order+1], "Debug"];
sol2 = Solve[(eq/. Vuv[order+1] -> Vn), Wuv[order+1]];
Wn = First[Wuv[order + 1] /. sol2];
PrintV[StringForm["W`2` is `1`", Wn, order+1], "Debug"];

{Join[previousOrders, {Vuv[order + 1] -> Vn, Wuv[order + 1]-> Wn}], order + 1}
];

{sols, ford} = Nest[NextCoef[#[[1]], #[[2]]]&, {{Vuv[0]-> V0, Wuv[0] -> W0}, 0}, 2];
sols
]


Options[QuarkMass] = {MassAUV -> 120 (*At which A to extract mq. Larger gives smaller errors in the expansion but larger errors numerically.*)
, mqErrorLimit -> 10^(-2),
  mqMethod -> {FiniteDifference, 10^-3}};
QuarkMass::mqErrorLarge = "The estimated error `1` in mq is larger than the set limit `2`";
QuarkMass::SolutionNotValid = "The `1`solution was not valid at \[Lambda]h = `2`, \[Tau]h = `3`";

QuarkMass[{q_, f_, \[Lambda]_, \[Tau]d_, \[Tau]_, Amin_, Amax_, \[CapitalLambda]scale_, fscale_, success_, bcoefs_, ell_}, {_, Vf_, \[Kappa]_?((!ListQ[#] && Head[#] =!= Rule)&), \[Omega]_}, opts:OptionsPattern[]]:= 
	If[success == True, QuarkMass[\[Tau][A],\[Lambda][A],  Amin, Amax, \[CapitalLambda]scale, ell, GammaFromKappa[\[Kappa], bcoefs[[1]], Vf, ell], opts], Message[QuarkMass::SolutionNotValid, "scaled ", \[Lambda][Amin], \[Tau][Amin]]; Undefined];

QuarkMass[{q_, f_, \[Lambda]_, \[Tau]d_, \[Tau]_, Amin_, Amax_, \[CapitalLambda]scale_, fscale_, success_, bcoefs_, ell_}, \[Gamma]_?NumericQ, opts:OptionsPattern[]]:= 
	If[success == True, QuarkMass[\[Tau][A],\[Lambda][A],  Amin, Amax, \[CapitalLambda]scale, ell, \[Gamma], opts], Message[QuarkMass::SolutionNotValid, "scaled ", \[Lambda][Amin], \[Tau][Amin]];Undefined];

QuarkMass[\[Tau]_InterpolatingFunction[arg_], \[Lambda]_InterpolatingFunction[arg_], Amin_, Amax_, \[CapitalLambda]_,ell_, \[Gamma]_,  opts : OptionsPattern[]]:=  Module[{\[Tau]AUV, AUV,mqbest, mq, deltamq, mqfun},
	Block[{$vcontext = "QuarkMass"},
	mq[AUV_]:=Exp[AUV]\[Tau][AUV]/ell^2/(AUV - Log[ell])^\[Gamma]/\[CapitalLambda];
	
	If[(Length[OptionValue[mqMethod]] > 1) && First[OptionValue[mqMethod]] === FiniteDifference, 
    (   
	  mqfun[AUV_] = Module[{AUV2 = AUV*(1 - OptionValue[mqMethod][[2]]) + OptionValue[mqMethod][[2]] Amin},
		 -(mq[AUV2]AUV2 - mq[AUV]AUV)/(AUV - AUV2)
	  ]
    ),
	  mqfun[AUV_] = mq[AUV] + AUV mq'[AUV];
    ];

	mqbest = mqfun[OptionValue[MassAUV]];

	(*Error control*)
	deltamq = mqbest - mqfun[(Amin +OptionValue[MassAUV])/2];
	PrintV[StringForm["mq = `1`, estimated error `2`", mqbest, deltamq], "Checks"];
	If[deltamq > OptionValue[mqErrorLimit], Message[QuarkMass::mqErrorLarge, deltamq, OptionValue[mqErrorLimit]]];

	PrintV[Plot[mqfun[Auv], {Auv, Amin, OptionValue[MassAUV]}, AxesLabel-> {"A", "mq"}], "Debug"];
	
	mqbest
	]
]

QuarkMass[{sols_, success_},ell_?NumericQ, \[Gamma]_?NumericQ, opts: OptionsPattern[]] := Module[{\[Tau]0, \[Tau]1, bcoefs, Amin, Amax, \[Lambda]0},
	(*If the solution does not reach the boundary, there's nothing we can do*)
	If[success == False, Message[QuarkMass::SolutionNotValid, "", (\[Lambda] /. sols)[0], (\[Tau] /. sols)[0]]; Return[Undefined]];

	\[Tau]0 = RemoveImaginaries[\[Tau] /. sols];
	\[Lambda]0 = RemoveImaginaries[\[Lambda] /. sols];
	{Amin, Amax} = First[InterpolatingFunctionDomain[\[Tau]0]];
	QuarkMass[\[Tau]0[A], \[Lambda]0[A], Amin, Amax, 1, ell, \[Gamma], opts]
]


\[Tau]hFromQuarkMass::largeerror = "The relative error `2` in mq for the final solution exceeds the set limit `3`. Actual mq = `1`";
\[Tau]hFromQuarkMass::largeabserror = "The actual value of mq = `1` exceeds the set limit `2`.";
Options[\[Tau]hFromQuarkMass] := Join[{MaxRecursion -> 80,(*How deep to recurse at maximum*)
NotFound -> Undefined,(*What to return if no result is found*)
InitialGuess -> 1, (*Initial guess for \[Tau]h. The algorithm seems quite insensitive to this, so there should be no need to change it*)
NodeCountSubdivision -> 10000, (*Once an interval where a solution exists is found, to how many parts should the interval be subdivided to approximate the number of nodes*)
NodeNumber -> 0, (*Which node to look for. UNTESTED for > 0*)
NodeCountHeuristic -> 20, (*A parameter for a heuristic speed up: once the value of mq exceeds the largest negative value encountered by this factor, the search is stopped.*)
NodeCountMinSteps -> 20, (*The minimum number of points to evaluate*)
NodeCountSteps -> 5,
NodeStretch -> 1,
ScaleAlways -> False, (*Forces the algorithm to always scale the result, affects some bordercases. Slightly slower.*)
MqErrorLimit :> If[mqv != 0, 10^(-2), 2*10^(-5)] (*An error is generated if the relative error in value of mq(\[Tau]h) is larger than this limit*)
}, Options[SolveAndScaleVQCDBH], Options[QuarkMass], Options[bCoefsFromPotential], Options[GammaFromKappa]]
\[Tau]hFromQuarkMass[mq_?NumericQ,\[Lambda]h_?NumericQ, nt_?NumericQ, pots_List,  opts: OptionsPattern[]] := Module[{expression, \[Tau]hinitial, mqat\[Tau]hinitial, \[Tau]hmin, \[Tau]hmax, \[Tau]hv, bcoefs, ell, \[Gamma], \[Tau]hsol, mqact},
Block[{$vcontext = "\[Tau]hFromQuarkMass"},

(*First check that there is any chance: is the effective potential real at \[Tau]h -> Infty*)
(*This is necessary in Mathematica 9, which will attempt to evaluate bignums when \[Tau]h is ridiculously large, and choke on it.*)
Veff[\[Tau]h_] = VEffective[pots][\[Lambda]h, nt, \[Tau]h];
If[!(Limit[Veff[\[Tau]h], \[Tau]h -> Infinity] > 0),
 (
	PrintV[StringForm["Veff < 0 at \[Lambda]h = `1`, nt = `2`", \[Lambda]h, nt], "Debug"];
	Return[OptionValue[NotFound]];
 ),
 PrintV["Veff > 0", "Debug"];,
 PrintV["Sign[Veff] Undetermined", "Debug"];
];

PrintV[StringForm["\[Tau]hFromQuarkMassOptions: `1`", opts], "Debug"];

(*b0, kappa and \[Gamma] are needed repeatedly, so solve them now*)
{bcoefs, ell} = bCoefsFromPotential[pots, Evaluate[FilterRules[{opts}, Options[bCoefsFromPotential]]]];
\[Gamma] = GammaFromKappa[pots, First[bcoefs], ell(*, Evaluate[FilterRules[{opts}, Options[GammaFromKappa]]]*)];
(*First choose whether we need to bother with scaled output or not*)
If[mq == 0 && !OptionValue[ScaleAlways],
(*mq == 0 is special, since then we don't have to scale, which makes the solution considerably faster*)
(expression := QuarkMass[SolveVQCDBH[\[Lambda]h, \[Tau]hv, nt, pots,Evaluate[FilterRules[{opts}, Options[SolveVQCDBH]]]], ell, \[Gamma], Evaluate[FilterRules[{opts}, Options[QuarkMass]]]](* // (Sign[Re[#]]Abs[#])&*);),
(expression :=  QuarkMass[SolveAndScaleVQCDBH[\[Lambda]h, \[Tau]hv, nt, pots,Evaluate[FilterRules[{opts}, Options[SolveAndScaleVQCDBH]]]], \[Gamma], Evaluate[FilterRules[{opts}, Options[QuarkMass]]]](* // (Sign[Re[#]]Abs[#])&*);)
];
\[Tau]hinitial = OptionValue[InitialGuess];
mqat\[Tau]hinitial = Quiet[expression /. \[Tau]hv-> \[Tau]hinitial];
PrintV[StringForm["mq at \[Tau]hinitial = `1`", mqat\[Tau]hinitial], "Debug"];

(*Look for a point where mq(\[Tau]h) < mq recursively, starting from the interval {0, \[Tau]hinitial}*)
If[!TrueQ[mqat\[Tau]hinitial <mq],
(\[Tau]hmin = NestWhile[Module[{mqnew, \[Tau]hnew},
If[TrueQ[#[[3]] > mq],
(
\[Tau]hnew = (#[[1]] + #[[2]])/2.;
mqnew = Quiet[(expression /. \[Tau]hv -> (#[[1]] + #[[2]])/2)];
PrintV[StringForm["Evaluated mq at \[Tau]h = `1`, mq = `2`", \[Tau]hnew, mqnew], "Debug"];
If[mqnew === Undefined,(*If it's undefined at the midpoint, make it the left end of the interval*)
{\[Tau]hnew, #[[2]], #[[3]]},
If[mqnew > mq, (*If it's larger than mq, make it the right end of the interval*)
{#[[1]], \[Tau]hnew, mqnew},
{\[Tau]hnew, \[Tau]hnew , mqnew} (*if it's less than mq, signal completion by shrinking the interval to a point*)
]
]),
PrintV[StringForm["Not true that #[[3]] > mq, #[[3]] = `1`, mq = `2`, doubling to `3`.", #[[3]], mq, 2*#[[2]]], "Debug"];{#[[1]], 2*#[[2]], Quiet[expression /. \[Tau]hv -> #[[2]]*2]}
]
]&, {0, \[Tau]hinitial, mqat\[Tau]hinitial},
(#[[1]] =!= #[[2]]) || (#[[1]] === Undefined)& (*End recursion when the interval has shrunk to a point*),
1, OptionValue[MaxRecursion]
]),
\[Tau]hmin = {\[Tau]hinitial, \[Tau]hinitial, mqat\[Tau]hinitial}];

(*Return if the area where mq(\[Tau]h) < mq was not found*)
If[(\[Tau]hmin[[1]] =!= \[Tau]hmin[[2]]) || (\[Tau]hmin[[1]] === Undefined) || \[Tau]hmin[[3]] > mq,
PrintV[StringForm["Could not find \[Tau]h where mq(\[Tau]h) < mq, best interval =  [`1`, `2`], mq[`1`] = `3`", \[Tau]hmin[[1]], \[Tau]hmin[[2]], \[Tau]hmin[[3]]],"Checks"];
Return[OptionValue[NotFound]];,
PrintV[StringForm["Found where mq(\[Tau]h) < mq, \[Tau]hmin = `1`", \[Tau]hmin], "Debug"];
];

PrintV["Looking for a point where mq(\[Tau]h) > mq", "Debug"];

(*Then find a point where mq(\[Tau]h) > mq*)
\[Tau]hmax = NestWhile[Module[{mqnew},
 PrintV[StringForm["Evaluating at \[Tau]h = `1`, prev mq  = `2`", #[[1]]*1.5, #[[2]]], "Debug"];
 mqnew = Quiet[expression /.\[Tau]hv -> #[[2]]]; 
 If[mqnew === Undefined,
 {#[[1]], (#[[1]] + #[[2]])/2, mqnew},(*It's undefined, so we need to retry between the two branch points*)
  If[mqnew > mq,
   {#[[2]], #[[2]], mqnew}, (*We're done*)
   {#[[2]], #[[2]] * 1.5, mqnew}(*Need to try higher*)
 ]
 ]
 ]&, 
 {\[Tau]hmin[[1]], \[Tau]hmin[[1]] * 1.5, mqat\[Tau]hinitial}, 
 !(#[[1]] == #[[2]])&, 1, OptionValue[MaxRecursion]];

If[!TrueQ[\[Tau]hmax[[3]] >mq],
(
PrintV["Could not find \[Tau]h where mq(\[Tau]h) > mq","Checks"];
Return[OptionValue[NotFound]];
)
];
PrintV[StringForm["\[Tau]hmin = `1`, \[Tau]hmax = `2`, mq at \[Tau]h initial = `3`", \[Tau]hmin, \[Tau]hmax, mqat\[Tau]hinitial], "Debug"];

(*Then count the number of nodes*)
If[OptionValue[NodeCountSubdivision] > 0,
(*Localize this, since a number of variables are needed*)
Module[{centers, left, right, l, n, sign, point, node, pbmax, mqleft, lastnodel, lastnoder, last2node, steplen, steps23, leftzero, rightzero, subdiv, oldsl, lastextend},
left = First[\[Tau]hmin];
right = First[\[Tau]hmax];
PrintV[StringForm["Subdividing interval [`1`, `2`]", left, right], "Debug"];
n = OptionValue[NodeCountSubdivision];
l = (right - left)/n;
(*Make a table of the centers of the subdivision intervals*)
(*An exponential progression seems to be more useful than linear*)
(*centers = Table[left*(right/left)^((1/2 + k)/n), {k, 0, n-1}];*)
node = 0;

mqleft = Quiet[expression /. \[Tau]hv -> left];

lastnodel = left;
lastnoder = left;
steplen = N[left*(right/left)^(1/2/n)-left];
PrintV[StringForm["Initial steplen `1`", steplen], "Debug"];

oldsl = steplen;

(*Now compute two more steps, if they're not of the same sign subdivide, rinse, and repeat*)
subdiv[] := Module[{masses}, steplen = steplen/4; PrintV[StringForm["New steplen = `1`", steplen],"Debug"];
	masses = {Re[Quiet[expression /. \[Tau]hv -> (left + steplen)]], Re[Quiet[expression /. \[Tau]hv -> (left + 2 steplen)]]};
	PrintV[StringForm["New attempt: `1`", masses], "Debug"]; 
	masses];
steps23 = NestWhile[subdiv[]&
,

{Quiet[expression /. \[Tau]hv -> (left + steplen)], Quiet[expression /. \[Tau]hv -> (left + 2 steplen)]}, (*Initial condition: count the first two points*)

((Sign[#[[1]]] == 1 ||Sign[#[[2]]] == 1) && (#[[1]] =!= #[[2]]))& (*Continue as long as at least one of the signs is positive. If the two numbers are precisely equal, assume that we've come so close that they are same anyway*)
];

PrintV[StringForm["Initial estimate points are `1`, `2`, `3`", mqleft, steps23[[1]], steps23[[2]]], "Debug"];

If[(mqleft != steps23[[1]]) && (steps23[[1]] != steps23[[2]]),
(
(*Now the array {-maxmq, steps23} contains three negative masses spaced at steplen intervals*)
(*The intersections of the quadratic passing through the three mq points with the constant c, at a coordinate system where left = 0, steplen = 1*)
leftzero = (3*y1 - 4*y2 + y3 - Sqrt[(-3*y1 + 4*y2 - y3)^2 - 4*(-2*c + 2*y1)*(y1 - 2*y2 + y3)])/
 (2*(y1 - 2*y2 + y3))/.{y1 -> mqleft, y2 -> steps23[[1]], y3 -> steps23[[2]], c -> mq};
rightzero = (3*y1 - 4*y2 + y3 + Sqrt[(-3*y1 + 4*y2 - y3)^2 - 4*(-2*c + 2*y1)*(y1 - 2*y2 + y3)])/
 (2*(y1 - 2*y2 + y3))/.{y1 -> mqleft, y2 -> steps23[[1]], y3 -> steps23[[2]], c -> mq};
pbmax = y1 - ((-3*y1)/2 + 2*y2 - y3/2)^2/(2*(y1 - 2*y2 + y3)) /. {y1 -> mqleft, y2 -> steps23[[1]], y3 -> steps23[[2]], c -> mq};
PrintV[StringForm["Maximum estimate from parabolic fit: `1`", pbmax], "Debug"];
If[Im[leftzero] != 0|| Im[rightzero] != 0, 
PrintV[StringForm["Imaginary values encountered at lah = `1`, leftzero = `2`, rightzero = `3`, the points are {`4`, `5`, `6`}, stepsize is `7`",
\[Lambda]h, leftzero, rightzero, mqleft, steps23[[1]], steps23[[2]], steplen], "Debug"]
(*Return[Undefined]*),
PrintV[StringForm["leftzero = `1`, rightzero = `2`, stepsize is `3`",
leftzero, rightzero, steplen], "Debug"];
];

(*Guess that our parabolic estimate is correct within a factor of defined by NodeCountSteps. If the zeros are in the wrong order,
that means they're in the wrong interval, the coef of the n^2 term is negative. In that case, use the absolute value and hope that
it still gives an order-of-magnitude estimate on the curvature of mq.*)
steplen = Abs[(rightzero - leftzero)]*steplen/OptionValue[NodeCountSteps];
PrintV[StringForm["Parabolic estimate for steplength is `1`, `2` steps", steplen, rightzero - leftzero], "Debug"];

)
,
(
steplen = oldsl;
pbmax = 0;
PrintV[StringForm["All masses are equal, reverting to steplen `1`", steplen],"Debug"];
)
];



(* Params in the function: {Sign, tauh, maxmq, mq, oldmq, steps}*)

lastextend = 0;
(*Look for a change of sign*)
point = NestWhile[
	Module[{next\[Tau]h = #[[2]] + steplen, cmq, maxmq = #[[3]]},
   (
	cmq = Re[Quiet[expression /. \[Tau]hv -> next\[Tau]h]];
	If[-cmq > maxmq, maxmq = -cmq];
	PrintV[StringForm["Step: # = `1`, \[Tau]h = `2`, lastextend = `3`", #, #[[2]], lastextend], "Debug"];
	If[cmq === Undefined, (*Sometimes, very near to the lower limit, the solver fails. Try to save things by assuming that this was the same as the previous point*)
		(PrintV["There was an undefined mq","Debug"];
		{#[[1]], next\[Tau]h, maxmq, #[[4]], #[[5]], #[[6]] - 1})
	,
	If[(cmq < mq) != #[[1]], 
		((*We found a node*)
			node++;
			last2node = lastnodel;
			lastnodel = #[[2]];
			lastnoder = next\[Tau]h;
			steplen = Max[(lastnoder - last2node)*OptionValue[NodeStretch], steplen];
			lastextend = 0;
			PrintV[StringForm["Found node at `1`, new steplen `2`", next\[Tau]h, steplen], "Debug"];

			(*If there were very few steps in the node, we might have missed the maximum by a large factor. Compute a few points to get a better estimate*)
			If[cmq > mq && OptionValue[NodeCountMinSteps] - #[[6]] < 3, 
				Module[{start = #[[2]], delta = (next\[Tau]h - #[[2]])/10, newmqs},
				newmqs =  Quiet[((expression /. \[Tau]hv -> (start + # delta))& /@ {1, 2, 3, 4, 5, 6, 7, 8, 9})];
				PrintV[StringForm["New mqs = `1`", newmqs], "Debug"];
				maxmq = Quiet[Max[Append[newmqs, maxmq]]];
				];
			];

			{cmq < mq, next\[Tau]h, maxmq, cmq, #[[4]], OptionValue[NodeCountMinSteps]}),
		(
			If[(lastextend - #[[6]]) > 4, (steplen = 2 steplen; lastextend = #[[6]]; PrintV["Changing slowly, doubling stepsize","Debug"])];
			{#[[1]], next\[Tau]h, maxmq, cmq, #[[4]], #[[6]] - 1}
		)
		] (*Return a list {whether we're above or below mq, the current subinterval}*)
	]
	)]&, 
   {True, left, Max[-mqleft, -pbmax, -steps23[[1]], -steps23[[2]]], mqleft, mqleft, OptionValue[NodeCountMinSteps]}, (*Start from below*)
	AnalyzeConditionalV[((#[[3]]*OptionValue[NodeCountHeuristic] > #[[4]] || #[[1]] == True || (#[[5]] >= #[[4]]) || #[[6]] >= 0 ) && #[[2]] <=  right), "Debug"]& 
	];

(*(*Look for a change of sign*)
point = NestWhile[Module[{cmq = Quiet[expression /. \[Tau]hv -> centers[[#[[2]]]]]},
 (PrintV[StringForm["Step: # = `1`, mq = `2`, maxmq = `3`", #, cmq, maxmq], "Debug"];
	If[-cmq > maxmq, maxmq = -cmq];
	If[(cmq < mq) != #[[1]], (node++; lastnode = #[[2]];PrintV[StringForm["Found node at `1`, point `2`", centers[[lastnode]], lastnode], "Debug"];)];
	 {cmq < mq, #[[2]] + 1, cmq, If[lastnode != #[[2]], #[[4]] + 1, 0]})]&, (*Return a list {whether we're above or below mq, the current subinterval}*)
   {True, 1, maxmq, 0}, (*Start from below*)
	#[[4]] < OptionValue[NodeCountMinPoints] || (maxmq*OptionValue[NodeCountHeuristic] > #[[3]] && #[[2]] <=  n)& (*End recursion when we've either seen enough sign changes or gone through the list*)
	];
*)
PrintV[StringForm["point `1`, node = `2`", point, node], "Debug"];
(*If[point[[2]] <= n, (\[Tau]hmin = {If[lastnode > 1, centers[[lastnode-1]], left]}; \[Tau]hmax = {centers[[lastnode]]}),
  (*Error control for the multinode solutions. Does NOT work correctly yet.*)
  If[node <= OptionValue[NodeNumber] && First[point] == True, PrintV["Could not find the correct node"]; Return[OptionValue[NotFound]]];
 ];*)
 If[lastnodel != lastnoder,
 \[Tau]hmin = {lastnodel}; \[Tau]hmax = {lastnoder};
 ];
];

];

PrintV[StringForm["Done subdividing, searching for solution in interval [`1`, `2`]", \[Tau]hmin, \[Tau]hmax], "Debug"];

(*Now that we've guaranteed the existence of a solution, find it by FindRoot*)
PrintV["Solving for a better root", "Debug"];
(*The check is needed to quiet any error messages (we'll check the result itself later),
 due to a bug in mathematica: http://mathematica.stackexchange.com/questions/20099/quiet-doesnt-work-with-findroot-when-using-brent-method*)
(*On the other hand, if the result is still a number, we'll want to use it.*)
\[Tau]hsol = Module[{\[Tau]hsol1}, Quiet@Check[\[Tau]hsol1 = \[Tau]hv /. FindRoot[Quiet[expression] == mq, {\[Tau]hv,First[\[Tau]hmin], First[\[Tau]hmax]}, Method-> "Brent"(*, AccuracyGoal -> 30, PrecisionGoal -> 29*)], If[NumericQ[\[Tau]hsol1],\[Tau]hsol1, Undefined, Undefined]]];

(*Error control*)
(*If this produces errors, it should be displayed, since they'll affect the the reliability of the error estimate*)
(*We need to scale here even when mq = 0, since the unscaled error is not relevant compared to dimensionful variables.*)
mqact = QuarkMass[SolveAndScaleVQCDBH[\[Lambda]h, \[Tau]hsol, nt, pots, Evaluate[FilterRules[{opts}, Options[SolveAndScaleVQCDBH]]]], \[Gamma], Evaluate[FilterRules[{opts}, Options[QuarkMass]]]]; 
PrintV[StringForm["Actual mq(`2`) = `1`", mqact, \[Tau]hsol], "Checks"];

(*This is a bit of a hack, but we'll have to see how it works in practice:*)
(*If the actual value of mq is not a number, look for the lower limit of existence. It is
often, in limiting cases where the solution barely exists, a very good approximation in itself
Whether this is actually the case, gets of course checked.*)
If[!(NumericQ[mqact] && NumericQ[\[Tau]hsol]),
 PrintV[StringForm["Solution not numeric. \[Tau]hmin = `1`, \[Tau]hmax = `2`", \[Tau]hmin, \[Tau]hmax], "Debug"]; 
 If[NumericQ[\[Tau]hsol],
		(*The \[Tau]hsol obtained is numeric, but mq is not: we'll just go towards larger th *)
		\[Tau]hsol = FindNumLimit[(Quiet[expression] /. \[Tau]hv -> #)&, {\[Tau]hsol, First[\[Tau]hmax]}];,

		(*\[Tau]hsol is also not numeric, i.e., mq was not defined at either end of the interval.*)
		(*Try looking for a numeric transition between \[Tau]hmax and say 1.5 \[Tau]hmax. It might*)
		(*be good enough.*)
		\[Tau]hsol = FindNumLimit[(Quiet[expression] /. \[Tau]hv -> #)&, {0, First[\[Tau]hmax]}];
 ];
 (*Recompute mqact*)	
 mqact = QuarkMass[SolveAndScaleVQCDBH[\[Lambda]h, \[Tau]hsol, nt, pots, Evaluate[FilterRules[{opts}, Options[SolveAndScaleVQCDBH]]]], \[Gamma], Evaluate[FilterRules[{opts}, Options[QuarkMass]]]];
 PrintV[StringForm["mq was not numeric, retried: actual mq(`2`) = `1`", mqact, \[Tau]hsol], "Checks"]; 
];

(*The next line should be modified, it now always produces an error for mq = 0*)
If[mq != 0, 
If[Abs[mqact - mq] > (OptionValue[MqErrorLimit]/.mqv -> mq)mq,Message[\[Tau]hFromQuarkMass::largeerror, mqact, Abs[mqact-mq]/mq, OptionValue[MqErrorLimit]]];,
(*mq zero error message*)
If[Abs[mqact - mq] > (OptionValue[MqErrorLimit]/.mqv -> mq), Message[\[Tau]hFromQuarkMass::largeabserror, mqact, OptionValue[MqErrorLimit]/.mqv-> mq]];
];

(*Return the best guess anyway*)
\[Tau]hsol
]
]


qFromSols[{q_, f_, \[Lambda]_, \[Tau]d_, \[Tau]_, Amin_, Amax_, \[CapitalLambda]scale_, fscale_, success_, bcoefs_, ell_}] :=
If[success == True,
q
,
Undefined];
\[Lambda]FromSols[{q_, f_, \[Lambda]_, \[Tau]d_, \[Tau]_, Amin_, Amax_, \[CapitalLambda]scale_, fscale_, success_, bcoefs_, ell_}] :=
If[success == True,
\[Lambda]
,
Undefined];
fFromSols[{q_, f_, \[Lambda]_, \[Tau]d_, \[Tau]_, Amin_, Amax_, \[CapitalLambda]scale_, fscale_, success_, bcoefs_, ell_}] :=
If[success == True,
f
,
Undefined];
bFromSols[{q_, f_, \[Lambda]_, \[Tau]d_, \[Tau]_, Amin_, Amax_, \[CapitalLambda]scale_, fscale_, success_, bcoefs_, ell_}] :=
If[success == True,
Function[A, Exp[A]/\[CapitalLambda]scale]
,
Undefined];
\[Tau]FromSols[{q_, f_, \[Lambda]_, \[Tau]d_, \[Tau]_, Amin_, Amax_, \[CapitalLambda]scale_, fscale_, success_, bcoefs_, ell_}] :=
If[success == True,
\[Tau]
,
Undefined];
Options[zFromSols] = Options[NIntegralFunction];
zFromSols[{q_, f_, \[Lambda]_, \[Tau]d_, \[Tau]_, Amin_, Amax_, \[CapitalLambda]scale_, fscale_, success_, bcoefs_, ell_}, opts : OptionsPattern[]] := NIntegralFunction[\[CapitalLambda]scale Exp[-A]q[A], {A, Amax, Amin}, opts, AccuracyGoal -> 80]


ARangeFromSols[{q_, f_, \[Lambda]_, \[Tau]d_, \[Tau]_, Amin_, Amax_, \[CapitalLambda]scale_, fscale_, success_, bcoefs_, ell_}] := {Amin, Amax};
AMinFromSols[{q_, f_, \[Lambda]_, \[Tau]d_, \[Tau]_, Amin_, Amax_, \[CapitalLambda]scale_, fscale_, success_, bcoefs_, ell_}] := Amin;
AMaxFromSols[{q_, f_, \[Lambda]_, \[Tau]d_, \[Tau]_, Amin_, Amax_, \[CapitalLambda]scale_, fscale_, success_, bcoefs_, ell_}] := Amax;


\[CapitalLambda]ScaleFromSols[{q_, f_, \[Lambda]_, \[Tau]d_, \[Tau]_, Amin_, Amax_, \[CapitalLambda]scale_, fscale_, success_, bcoefs_, ell_}]:=
If[success == True,
\[CapitalLambda]scale
,
Undefined];
fScaleFromSols[{q_, f_, \[Lambda]_, \[Tau]d_, \[Tau]_, Amin_, Amax_, \[CapitalLambda]scale_, fscale_, success_, bcoefs_, ell_}]:=
If[success == True,
fscale
,
Undefined];
qhFromSols[{q_, f_, \[Lambda]_, \[Tau]d_, \[Tau]_, Amin_, Amax_, \[CapitalLambda]scale_, fscale_, success_, bcoefs_, ell_}]:=
If[success == True,
q[Amin]
,
Undefined];

qhFromThermoData[fscale_, \[Lambda]h_, \[Tau]h_, nt_, {Vg_, Vf_, \[Kappa]_, \[Omega]_}] := -Sqrt[3] fscale/(Sqrt[Vg[\[Lambda]h] - Vf[\[Lambda]h, \[Tau]h]Sqrt[1+ nt^2/(\[Omega][\[Lambda]h]^2 Vf[\[Lambda]h, \[Tau]h]^2)]]);

TemperatureFromThermoData[fscale_, \[CapitalLambda]_, \[Lambda]h_, \[Tau]h_, nt_, pots_] := -1/(4 Pi) fscale^2 / \[CapitalLambda] / qhFromThermoData[fscale, \[Lambda]h, \[Tau]h, nt, pots];

TemperatureFromSols[{q_, f_, \[Lambda]_, \[Tau]d_, \[Tau]_, Amin_, Amax_, \[CapitalLambda]scale_, fscale_, success_, bcoefs_, ell_}]:=
If[success == True,
-1/(4 Pi )fscale^2/q[Amin]/\[CapitalLambda]scale
,
Undefined];

s4G5FromSols[{q_, f_, \[Lambda]_, \[Tau]d_, \[Tau]_, Amin_, Amax_, \[CapitalLambda]scale_, fscale_, success_, bcoefs_, ell_}]:=
If[success == True,
1/\[CapitalLambda]scale^3
,
Undefined];

ntildeFromSols[{q_, f_, \[Lambda]_, \[Tau]d_, \[Tau]_, Amin_, Amax_, \[CapitalLambda]scale_, fscale_, success_, bcoefs_, ell_}, nt_]:=
If[success == True,
nt/\[CapitalLambda]scale^3
,
Undefined];

APrimeFromSols[{q_, f_, \[Lambda]_, \[Tau]d_, \[Tau]_, Amin_, Amax_, \[CapitalLambda]scale_, fscale_, success_, bcoefs_, ell_}, nt_, {_, Vf_, \[Kappa]_, _}] :=
If[success == True,
Function[A, -q[A] Exp[A] Sqrt[(1+ \[Tau]'[A]^2 f[A] \[Kappa][\[Lambda][A]]/q[A]^2)nt^2 \[CapitalLambda]scale^(-6)/(nt^2 \[CapitalLambda]scale^(-6) + Exp[6 A] \[CapitalLambda]scale^(-6) \[Kappa][\[Lambda][A]]^2 Vf[\[Lambda][A], \[Tau][A]]^2)]/(\[CapitalLambda]scale \[Kappa][\[Lambda][A]])]

];

Options[AFromSols] = Join[{Method -> "NIntegralFunction", Densification -> 4}, Options[NIntegralFunction]];
AFromSols[{q_, f_, \[Lambda]_, \[Tau]d_, \[Tau]_, Amin_, Amax_, \[CapitalLambda]scale_, fscale_, success_, bcoefs_, ell_}, nt_, pots_List, opts : OptionsPattern[]] := Module[{Ap},Block[{$vcontext = "AFromSols"},
 PrintV[StringForm["Integrating A_0, options are `1`", opts], "Debug"];
 Ap = APrimeFromSols[{q, f, \[Lambda], \[Tau]d, \[Tau], Amin, Amax, \[CapitalLambda]scale, fscale, success, bcoefs, ell}, nt, pots];
 If[OptionValue[Method] === "NIntegralFunction",
 NIntegralFunction[Ap[A], {A, Amin, Amax}, opts],
 ((*The other method is FromInterpolatingFunction*)
  Module[{grid = InterpolatingFunctionGrid[Head[\[Lambda][A]]], data, intp, Afun},
  data = {#, Ap[#]}& /@ Densify[grid, OptionValue[Densification]];
  intp = Interpolation[data];
  Afun[Av_]  = Integrate[intp[A], {A, Amin, Av}];
  Afun
 ]
 )
 ]
 
 ]
];

Options[AAndMuFromSols] = Options[AFromSols];
AAndMuFromSols[{q_, f_, \[Lambda]_, \[Tau]d_, \[Tau]_, Amin_, Amax_, \[CapitalLambda]scale_, fscale_, success_, bcoefs_, ell_}, nt_, {Vg_, Vf_, \[Kappa]_, \[Omega]_}, opts : OptionsPattern[]] := 
Module[{Afunc}, Afunc = AFromSols[{q, f, \[Lambda], \[Tau]d, \[Tau], Amin, Amax, \[CapitalLambda]scale, fscale, success, bcoefs, ell}, nt, {Vg, Vf, \[Kappa], \[Omega]}, opts];
{Function[A, Afunc[A] - Afunc[Amin]], Afunc[Amax] - Afunc[Amin]}
];



End[ ]
EndPackage[ ]
