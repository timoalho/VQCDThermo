(* ::Package:: *)

(***********************************************
ThermoFactory.m

Scaffolding for easy computation of thermodynamics in the VQCD -model.
Uses VQCDCore.m as the core solver.

Written by Timo Alho

Package created on 2013/9/17
************************************************)
BeginPackage["VQCDThermo`", 
{"VQCDThermo`VQCDCore`", "VQCDThermo`NumericsExtra`", "VQCDThermo`Verbosity`"}]

FindThermoComputationBoxes::usage = "Splits the physical region defined by the potentials in to the \[Lambda]h intervals that we need to compute."
FindThermoTachyonicComputationBoxes::usage = ""
MakentconstLists::usage = "Creates the list with lah ranges for various nt's"
Make\[Lambda]hRangesFromnt::usage = "Gives the lah ranges for handling a given nt"
Make\[Lambda]hconstLists::usage = "Creates a list of \[Lambda]h's with search ranges for nt"

ParamListToComputationList::usage = "Converts a list of nt values and \[Lambda]h ranges to the format suitable for ComputeTunedThermo."

ComputeParametricRawThermo::usage = "Computes the raw variables, which can then be used to form the thermo observables, along a parametrized path in the (nt1, \[Lambda]h) -plane"

fScalePointsFromData::usage = "Returns the data points for fscale from the saved data"
\[CapitalLambda]PointsFromData::usage = "Returns the \[CapitalLambda] datapoints"
\[Mu]PointsFromData::usage = "Returns the \[Mu] datapoints"
PlotPointsAndInterpolation::usage = "Given data, prints diagnostic plots with both the data points and interpolations"
PlotPointsAndInterpolationDirectory::usage = "Prints diagnostic plots for for all data in a given directory"

ComputeThermoCurve::usage = "Computes the thermo along a curve, finding the limits automatically etc."
ComputeTunedThermo::usage = "Given a thermo tuning, computes thermo"
ComputeAndSaveThermo::usage = "Computes and saves a thermo, given tuning"

ComputeStandardThermo::usage = "The most straightforward way to compute the thermo. Attempts to automatically analyze your potentials etc."

LoadThermo::usage = ""
LoadThermoDirectory::usage = ""

(*Export symbol names that are to be saved to file. Without these, the saved variables would be prepended with ThermoFactory`Private`, making
loading dependent on internal details of this package.
Ugly solution, suggestions for a better method are welcome.*)
CurveType;
Param;
Potentials;
Results;
\[Lambda]hFunction;
ntFunction;
\[Tau]hFunction;

Begin["`Private`"]


Options[ParamListToComputationList] = Options[ComputeTunedThermo]
ParamListToComputationList[curveType_String, list_List, pots_List, opts : OptionsPattern[]] := {curveType, pots, #[[1]], #[[2]], #[[3]], #[[4]], opts}&/@ list;


Options[FindThermoComputationBoxes] = Join[{\[Lambda]hrange -> {0.01, 100}}, Options[FindFunctionRoots]]
FindThermoComputationBoxes[pots_List, opts : OptionsPattern[]] := Module[{roots, Veff, ntfun, ntfun1, ntroots, \[Lambda]hnummax},
Veff[\[Lambda]h_, nt_, \[Tau]h_]= VEffective[{Vg, Vf, \[Kappa], \[Omega]}][\[Lambda]h, nt, \[Tau]h];
ntfun[\[Lambda]h_] = (nt /. Last[Quiet[Solve[Veff[\[Lambda]h, nt, 0] == 0, nt]]])^2 /. {Vf -> pots[[2]], Vg -> pots[[1]], \[Kappa]-> pots[[3]], \[Omega] -> pots[[4]]};
(*ntfun[\[Lambda]h_] = ntfun1[\[Lambda]h] ;*)
\[Lambda]hnummax = First[FindFunctionRoots[ntfun[\[Lambda]h], {\[Lambda]h, First[OptionValue[\[Lambda]hrange]], Last[OptionValue[\[Lambda]hrange]]}]];
(*Find the potential (no pun intended) extrema*)
roots = {Sqrt[ntfun[#]], #, ntfun''[#] < 0}&/@  FindFunctionRoots[ntfun'[\[Lambda]h], {\[Lambda]h, First[OptionValue[\[Lambda]hrange]], \[Lambda]hnummax}, Evaluate[FilterRules[{opts},
	Options[FindFunctionRoots]]]];
AppendTo[roots, {Sqrt[ntfun[First @ OptionValue[\[Lambda]hrange]]], First[OptionValue[\[Lambda]hrange]], ntfun'[First @ OptionValue[\[Lambda]hrange]] < 0}];

FindConsecutiveRegions[0, First[OptionValue[\[Lambda]hrange]], \[Lambda]hnummax, roots]
]


Options[FindThermoTachyonicComputationBoxes] = Options[FindThermoComputationBoxes]
FindThermoTachyonicComputationBoxes[pots_List, opts : OptionsPattern[]] := Module[{ntfun, roots},
	(*Approximate the zero mass limit by tauh = infty*)
	ntfun[\[Lambda]h_] = Limit[ntCritical[pots, \[Tau]h][\[Lambda]h], \[Tau]h -> Infinity];
	roots = {ntfun[#], #, ntfun''[#] < 0}& /@ FindFunctionRoots[ntfun'[\[Lambda]h], {\[Lambda]h, Sequence @@ OptionValue[\[Lambda]hrange]}, NumPoints -> 10000, Evaluate[FilterRules[{opts},
	Options[FindFunctionRoots]]]];
	(*Add the extrema at the limits of the interval*)
	roots = Join[roots, {ntfun[#], #, ntfun'[#] < 0}& /@ OptionValue[\[Lambda]hrange]];
	FindConsecutiveRegions[0, Sequence @@ OptionValue[\[Lambda]hrange], roots]

]


Options[Make\[Lambda]hRangesFromnt] = Options[FindThermoComputationBoxes]
Options[MakentconstLists] = Join[Options[Make\[Lambda]hRangesFromnt], {SpreadFunction -> Identity, MinCompBoxHeight -> 0.05}]

Make\[Lambda]hRangesFromnt[n_?NumericQ, compboxes_List?(ListQ[#[[1]]]&)] := {n, #[[3]], #[[4]], #[[5]]}& /@ Select[compboxes, #[[1]] <= n < #[[2]]&];

Make\[Lambda]hRangesFromnt[n_?NumericQ, pots_List] := Module[{compboxes}, compboxes = FindThermoComputationBoxes[pots, opts];
	Make\[Lambda]hRangesFromnt[n, compboxes]];

MakentconstLists[compboxes_List?(ListQ[#[[1]]]&), nCurves_Integer, opts : OptionsPattern[]] := Module[{nmax, nvals},
nmax = Max[#[[2]]&/@compboxes];
nvals = Table[OptionValue[SpreadFunction] @ n, {n, InverseFunction[OptionValue[SpreadFunction]] @ 0, 
 InverseFunction[OptionValue[SpreadFunction]] @ nmax, ((#[[2]] - #[[1]])/(nCurves + 1))& @  InverseFunction[OptionValue[SpreadFunction]] [{0, nmax}]}]; (*+1 since the limit case never gets included*)
Flatten[Make\[Lambda]hRangesFromnt[#, compboxes, opts]&/@nvals,1]
]

MakentconstLists[tachyon : "Tachyonic" | "NonTachyonic", pots_List, nCurves_Integer, opts : OptionsPattern[]] := Module[{compboxes, nmax, nvals},
compboxes = If[tachyon === "NonTachyonic", FindThermoComputationBoxes[pots, opts], FindThermoTachyonicComputationBoxes[pots, opts]];
MakentconstLists[Select[compboxes, (#[[5]] - #[[3]] >=  OptionValue[MinCompBoxHeight])&], nCurves]
]



Options[Make\[Lambda]hconstLists] = {\[Lambda]hrange -> {0.01, 100},
	\[Lambda]hmargin -> 10^-3, (*Points _very _ close to the limit usually behave badly numerically.*)
	SpreadFunction -> (Exp) (*The function with which to spread the points. Identity gives a linear spread*)}

Make\[Lambda]hconstLists["Tachyonic", pots_List, nCurves_Integer, opts : OptionsPattern[]] := Module[{Veff, ntfun, \[Lambda]hvals, \[Lambda]hnummax, \[Lambda]hr},
ntfun[\[Lambda]h_] = Limit[ntCritical[pots, \[Tau]h][\[Lambda]h], \[Tau]h -> Infinity];

\[Lambda]hr = OptionValue[\[Lambda]hrange];
\[Lambda]hvals = Table[OptionValue[SpreadFunction] @ i, Evaluate[{i,Sequence @@ InverseFunction[OptionValue[SpreadFunction]] @ \[Lambda]hr, ((#[[2]] - #[[1]])/nCurves)& @  InverseFunction[OptionValue[SpreadFunction]] [\[Lambda]hr]}]];

{#, 0, 0, ntfun[#]}& /@ \[Lambda]hvals
]

Make\[Lambda]hconstLists["NonTachyonic", pots_List, nCurves_Integer, opts : OptionsPattern[]] := Module[{Veff, ntfun, \[Lambda]hvals, \[Lambda]hnummax, \[Lambda]hr},
Veff[\[Lambda]h_, nt_, \[Tau]h_]= VEffective[{Vg, Vf, \[Kappa], \[Omega]}][\[Lambda]h, nt, \[Tau]h];
ntfun[\[Lambda]h_] = (nt /. Last[Quiet[Solve[Veff[\[Lambda]h, nt, 0] == 0, nt]]])^2 /. {Vf -> pots[[2]], Vg -> pots[[1]], \[Kappa]-> pots[[3]], \[Omega] -> pots[[4]]};

(*Compute the maximum \[Lambda]h by using the fact that at nt = 0 it is precisely when Veff' = 0*)
\[Lambda]hnummax = (\[Lambda]h(1 - OptionValue[\[Lambda]hmargin]))/. FindRoot[Derivative[1, 0, 0][Veff][\[Lambda]h, 0, 0] /. {Vf -> pots[[2]], Vg -> pots[[1]], \[Kappa]-> pots[[3]], \[Omega] -> pots[[4]]}, {\[Lambda]h, Sequence@@OptionValue[\[Lambda]hrange]}];

\[Lambda]hr = {First[OptionValue[\[Lambda]hrange]], \[Lambda]hnummax};
\[Lambda]hvals = Table[OptionValue[SpreadFunction] @ i, Evaluate[{i,Sequence @@ InverseFunction[OptionValue[SpreadFunction]] @ \[Lambda]hr, ((#[[2]] - #[[1]])/nCurves)& @  InverseFunction[OptionValue[SpreadFunction]] [\[Lambda]hr]}]];

{#, 0, 0, Sqrt[ntfun[#]]}& /@ \[Lambda]hvals
]


Options[ComputeParametricRawThermo] = Join[Options[FunctionInterpolation], Options[SolveAndScaleVQCDBH],
	{AIntegralAccuracy -> 60, (*Accuracy for computing the A, \[Mu] -integral*)
	InterpolationOptions -> {}, (*Options to be passed to FunctionInterpolation*)
	IndicatorFunction -> Function[{fscale, \[CapitalLambda], qh, \[Mu], \[Lambda]h, nt1, \[Tau]h}, -fscale^2 /\[CapitalLambda]/qh], (*This is the function that FunctionInterpolation is actually applied to.*)
	Interpolate\[Lambda]h -> False,
	Interpolatent -> False,
	Interpolate\[Tau]h -> False
}];
ComputeParametricRawThermo[umin_?NumericQ, umax_?NumericQ, nt_, \[Lambda]h_, \[Tau]h_, pots_List, opts : OptionsPattern[]] := 
 Module[{\[CapitalLambda]fun, fscalefun,tachyonopts, intpopts, unummax, unummin, pointlist, mufun, indicatorfun, outValues, fscaleidx, \[CapitalLambda]idx, \[Mu]idx, \[Lambda]hidx = Undefined, ntidx = Undefined, \[Tau]hidx = Undefined},

Block[{$vcontext = "ComputeParametricRawThermo"},

 PrintV[StringForm["umin = `1`, umax = `2`", umin, umax], "Debug"];
 PrintV[StringForm["IndicatorFunction is `1`", OptionValue[IndicatorFunction]], "Debug"];
 
(*Set up options for SolveAndScaleFiniteTTachyons*)
tachyonopts = FilterRules[{opts}, Options[SolveAndScaleVQCDBH]];

 (*Set up the filter for producing the table*)
 outValues[u_, fscale_, \[CapitalLambda]_, \[Mu]_, \[Lambda]hv_, ntv_, \[Tau]hv_] = Module[{cidx = 5, list = {u, fscale, \[CapitalLambda], \[Mu]}},
 fscaleidx = 2;
 \[CapitalLambda]idx = 3;
 \[Mu]idx = 4;
 If[OptionValue[Interpolate\[Lambda]h], AppendTo[list, \[Lambda]hv]; \[Lambda]hidx = cidx; cidx = cidx+1;];
 If[OptionValue[Interpolatent], AppendTo[list, ntv]; ntidx = cidx; cidx = cidx+1;];
 If[OptionValue[Interpolate\[Tau]h], AppendTo[list, \[Tau]hv]; \[Tau]hidx = cidx; cidx = cidx+1;];
 list];

 PrintV[StringForm["Listfun is `1`", DownValues[outValues]], "Debug"];

(*Compute the InterpolatingFunctions*)
pointlist = Reap[
PrintV["Constructing fscale[u]...", "Progress"];
PrintV[StringForm["Interpolation options are `1`", OptionValue[InterpolationOptions]], "Debug"];
unummax = umax;
unummin = umin;
indicatorfun = FunctionInterpolation[Module[{\[Tau]hval, fscale, sol, Afun, mu, \[CapitalLambda], qh, \[Lambda]hval, ntval},
\[Lambda]hval = \[Lambda]h[u];
ntval = nt[u];
\[Tau]hval = \[Tau]h[\[Lambda]hval, ntval];
PrintV[StringForm["Computing thermopoint at \[Lambda]h = `1`, nt = `2`, \[Tau]h = `3`", \[Lambda]hval, ntval, \[Tau]hval], "All"];
sol = SolveAndScaleVQCDBH[\[Lambda]hval, \[Tau]hval, ntval, pots, tachyonopts];
fscale = fScaleFromSols[sol];
\[CapitalLambda] = \[CapitalLambda]ScaleFromSols[sol];
qh = qhFromBoundaryData[fscale, \[Lambda]hval, \[Tau]hval, ntval, pots];
mu = AAndMuFromSols[sol, ntval, pots, AccuracyGoal -> OptionValue[AIntegralAccuracy]][[2]];
PrintV[StringForm["Thermo done, indicatorfun = `1`", OptionValue[IndicatorFunction][fscale, \[CapitalLambda], qh, \[Lambda]hval, nt1val, \[Tau]hval]], "All"];
If[NumericQ[fscale], (Sow[outValues[u, fscale, \[CapitalLambda], mu, \[Lambda]hval, ntval, \[Tau]hval]];OptionValue[IndicatorFunction][fscale, \[CapitalLambda], qh, mu, \[Lambda]hval, ntval, \[Tau]hval]),
   (If[unummax>u, unummax = u]; If[unummin < u, unummin = u];0)]
], {u, umin, umax}, Evaluate[Sequence @@ OptionValue[InterpolationOptions]]];

][[2]][[1]]; (*These index the list of sowed values*)
 
PrintV[StringForm["List length before removing duplicate points: `1`", Length[pointlist]], "Debug"]; 
pointlist = Union[pointlist, SameTest -> (#1[[1]] == #2[[1]]&)];
PrintV[StringForm["List length after removing duplicate points: `1`", Length[pointlist]], "Debug"]; 

PrintV[StringForm["Done, actual limits for u = [`1`, `2`], (nt1, \[Lambda]h) = [(`3`, `4`), (`5`, `6`\.1d)]", unummax, unummin, First[pointlist][[ntidx]], First[pointlist][[\[Lambda]hidx]], Last[pointlist][[ntidx]], Last[pointlist][[\[Lambda]hidx]]], "Debug"];

(*Return: the full list of computed points*)
{pointlist, {fscaleidx, \[CapitalLambda]idx, \[Mu]idx, \[Lambda]hidx, ntidx, \[Tau]hidx}}

]

]


ComputeThermoCurve::nosolution = "The solution is not numerical at the given midpoint, u = `1`."
Options[ComputeThermoCurve] := 
 Join[{IndicatorFunction -> Function[{fscale, \[CapitalLambda], qh, \[Mu], \[Lambda]h, nt1, \[Tau]h}, -fscale^2 /\[CapitalLambda]/qh], (*By default, use temperature as the interpolation criterion.*)
  InterpolationOptions -> {}, (*Additional/alternative options for the interpolation*)
 SolveAndScaleOptions -> {}, (*Additional/alternative options to the core solver*)
 \[Tau]hOptions -> {}, (*Additional/alternative options to pass to computing \[Tau]h*)
 \[Tau]hExistenceOptions -> {}, (*Additional/alternative options to pass to computing \[Tau]h when only testing it's existence*)
 NumLimitOptions -> {} (*Additional/alternative option to pass to FindNumLimit*)
 }, Options[ComputeParametricRawThermo]]
 ComputeThermoCurve[ntfun_, \[Lambda]hfun_, tachyonic_?((# == True) || (# == False)&), umin_?NumericQ, umid_?NumericQ, umax_?NumericQ, pots_List, opts : OptionsPattern[]] := Module[
  {numres, uminlim, umaxlim, \[Tau]hFun, time, existfun},
 If[!tachyonic,
 (\[Tau]hFun[\[Lambda]h_, nt_, existenceOnly_ : False] = 0;)
 ,
 (
  \[Tau]hFun[\[Lambda]h_?NumericQ, nt_?NumericQ, existenceOnly_ : False] = \[Tau]hFromQuarkMass[0, \[Lambda]h, nt, pots, Evaluate[If[existenceOnly, Options[\[Tau]hExistenceOptions], Options[\[Tau]hOptions]]], OptionValue[SolveAndScaleOptions], AccuracyGoal -> 60, NodeCountSubdivision-> 0, ARange -> 60, MassAUV -> 60]
 )
 ];

 existfun[u_?NumericQ] = \[CapitalLambda]ScaleFromSols[SolveAndScaleVQCDBH[\[Lambda]hfun[u], \[Tau]hFun[\[Lambda]hfun[u], ntfun[u], True], ntfun[u], pots, OptionValue[SolveAndScaleOptions]]];

 time = First[AbsoluteTiming[

 If[!(existfun[umid]), (Message[ComputentConstThermo::nosolution, umid]; Return[Undefined])];

 If[Quiet[!(NumericQ[existfun[umin]] === True)], uminlim = FindNumLimit[Quiet[existfun[u]], {u, umin, umid}, OptionValue[NumLimitOptions]], uminlim = umin];
 PrintV[StringForm["Found u lower limit `1`", uminlim], "Progress"];

 If[Quiet[!(NumericQ[existfun[umax]] === True)], umaxlim = FindNumLimit[Quiet[existfun[u]], {u, umid, umax}, OptionValue[NumLimitOptions]], umaxlim = umax];
 PrintV[StringForm["Found u upper limit `1`", umaxlim], "Progress"];

 ]];

 PrintV[StringForm["Found limits for solution in `1` seconds...", time], "Progress"];

 PrintV[StringForm["Starting to compute thermo interpolation. Started at `1`:`2`... ", DateList[][[4]], DateList[][[5]]], "Progress"];
 time = First[AbsoluteTiming[ 
 numres =  ComputeParametricRawThermo[uminlim, umaxlim, ntfun, \[Lambda]hfun, \[Tau]hFun[#1, #2]&, pots, 
               OptionValue[SolveAndScaleOptions], 
               Interpolate\[Tau]h -> tachyonic,
			   Evaluate[FilterRules[{opts},
	Options[ComputeParametricRawThermo]]]];
 ]];
 PrintV[StringForm["Done in `1` seconds", time], "Progress"];
 
 numres
]


ThermoTune["\[Lambda]hconstant", param_?NumericQ, pots_List] := {
	\[Lambda]hfun -> (param&), 
	ntfun -> (#&),
	tachyonic -> False,
	umap -> (#&),
	Interpolate\[Lambda]hDefault -> False,
	InterpolatentDefault -> False,
	NumLimitDefaults -> {PrecisionGoal -> 7, MaxRecursion -> 30},
	SolveAndScaleDefaults -> {AccuracyGoal -> 60, ARange -> 120, MaxSteps -> 100000},
	(*The proper tuning actually depends on which part of the la* curve we are. The condition below is _very _ heuristic, 
	but it works for the potentials tested this far.*)
    Sequence @@ If[Module[{ntfun = ntCritical[pots]}, param > Last[(FindFunctionRoots[D[ntfun[\[Lambda]h]^2, \[Lambda]h], {\[Lambda]h, 0, 20}])]],
	{
		InterpolationDefaults ->  {PrecisionGoal -> 6, AccuracyGoal -> 6, MaxRecursion -> 35},
		IndicatorDefaults -> Function[{fscale, \[CapitalLambda], qh, \[Mu], \[Lambda]h, nt1, \[Tau]h}, -fscale^2 /\[CapitalLambda]/qh]
	},
	{
		InterpolationDefaults -> {PrecisionGoal-> 5.3, AccuracyGoal -> 20, MaxRecursion-> 35},
		IndicatorDefaults -> Function[{fscale, \[CapitalLambda], qh, \[Mu], \[Lambda]h, nt1, \[Tau]h}, fscale^2  \[CapitalLambda] ]
	}] 
}


ThermoTune["ntconstant", param_?NumericQ, pots_List] = {
	\[Lambda]hfun -> (Exp[#]&), 
	ntfun -> (param&),
	tachyonic -> False,
	umap -> (Log[#]&),
	Interpolate\[Lambda]hDefault -> False,
	InterpolatentDefault -> False,
	NumLimitDefaults -> {},
	SolveAndScaleDefaults -> {AccuracyGoal -> 60, ARange -> 120, MaxSteps -> 100000},
	InterpolationDefaults -> {PrecisionGoal-> 5, AccuracyGoal -> 60, MaxRecursion-> 25},
	IndicatorDefaults -> Function[{fscale, \[CapitalLambda], qh, \[Mu], \[Lambda]h, nt1, \[Tau]h}, -fscale^2 /\[CapitalLambda]/qh]
}


ThermoTune["\[Lambda]hconstant\[Tau]h", param_?NumericQ, pots_List] = {
	\[Lambda]hfun -> (param&), 
	ntfun -> (#&),
	tachyonic -> True,
	umap -> (#&),
	Interpolate\[Lambda]hDefault -> False,
	InterpolatentDefault -> False,
	NumLimitDefaults -> {PrecisionGoal-> 5, MaxRecursion -> 25},
	SolveAndScaleDefaults -> {\[CapitalLambda]relerrLimit -> 1/500, AccuracyGoal -> 60, ARange -> 120},
	InterpolationDefaults -> {PrecisionGoal-> 6, AccuracyGoal -> 6, MaxRecursion-> 6},
	IndicatorDefaults -> Function[{fscale, \[CapitalLambda], qh, \[Mu], \[Lambda]h, nt1, \[Tau]h}, fscale^2 \[CapitalLambda] \[Mu]/\[Tau]h]
}


ThermoTune["ntconstant\[Tau]h", param_?NumericQ, pots_List] = {
	\[Lambda]hfun -> (Exp[#]&), 
	ntfun -> (param&),
	tachyonic -> True,
	umap -> (Log[#]&),
	Interpolate\[Lambda]hDefault -> False,
	InterpolatentDefault -> False,
	NumLimitDefaults -> {MaxRecursion -> 25},
	SolveAndScaleDefaults -> {\[CapitalLambda]relerrLimit -> 1/500, AccuracyGoal -> 60, ARange -> 120},
	InterpolationDefaults -> {PrecisionGoal-> 6, AccuracyGoal -> 10, MaxRecursion-> 6},
	IndicatorDefaults -> Function[{fscale, \[CapitalLambda], qh, \[Mu], \[Lambda]h, nt1, \[Tau]h}, -fscale^2 /\[CapitalLambda]/qh]
}


ThermoTune["\[Lambda]end", param_?NumericQ, pots_List] := Module[{\[Lambda]hfunloc},
	\[Lambda]hfunloc[u_?NumericQ] := FindNumLimit[\[Tau]hFromQuarkMass[0, \[Lambda]h, u, pots, AccuracyGoal -> 60, NodeCountSubdivision-> 0, ARange -> 120, MassAUV -> 110, MaxRecursion -> 80], {\[Lambda]h, 0.01, 100.}, MaxRecursion -> 120];
	{
	\[Lambda]hfun -> \[Lambda]hfunloc,
	ntfun -> (#&),
	tachyonic -> True,
	umap -> (#&),
	Interpolate\[Lambda]hDefault -> True,
	InterpolatentDefault -> False,
	NumLimitDefaults -> {},
	SolveAndScaleDefaults -> {\[CapitalLambda]relerrLimit -> 1/500},
	InterpolationDefaults -> {PrecisionGoal -> 2, AccuracyGoal -> 2, MaxRecursion -> 6},
	IndicatorDefaults -> Function[{fscale, \[CapitalLambda], qh, \[Mu], \[Lambda]h, nt1, \[Tau]h}, \[CapitalLambda]]
}
]


ThermoTune["ntend", param_?NumericQ, pots_List] := Module[{ntfunloc, ntcrit},
	ntcrit[\[Lambda]h_] = Limit[ntCritical[pots, \[Tau]h][\[Lambda]h], \[Tau]h -> Infinity];
	ntfunloc[u_?NumericQ] := (Print[StringForm["ntfun[`1`] = `2`", Exp[u], ntcrit[Exp[u]]]];FindNumLimit[\[Tau]hFromQuarkMass[0, Exp[u], nt, pots, AccuracyGoal -> 60, NodeCountSubdivision-> 0, ARange -> 120, MassAUV -> 110, MaxRecursion -> 80], {nt, 0, ntcrit[Exp[u]]}, MaxRecursion -> 120]);
	{
	\[Lambda]hfun -> (Exp[#]&),
	ntfun -> ntfunloc,
	tachyonic -> True,
	umap -> (Log[#]&),
	Interpolate\[Lambda]hDefault -> False,
	InterpolatentDefault -> True,
	NumLimitDefaults -> {},
	SolveAndScaleDefaults -> {\[CapitalLambda]relerrLimit -> 1/500},
	InterpolationDefaults -> {PrecisionGoal -> 2, AccuracyGoal -> 2, MaxRecursion -> 6},
	IndicatorDefaults -> Function[{fscale, \[CapitalLambda], qh, \[Mu], \[Lambda]h, nt1, \[Tau]h}, \[CapitalLambda]]
}
]


ComputeTunedThermo::nosolution = "The solution is not numerical at the given midpoint, \[Lambda]h = `1`."
Options[ComputeTunedThermo] = 
 {IndicatorFunction -> Function[{fscale, \[CapitalLambda], qh, \[Mu], \[Lambda]h, nt1, \[Tau]h}, fscale^2  \[CapitalLambda] ], (*By default, use temperature as the interpolation criterion.*)
  InterpolationOptions -> {}, (*Additional/alternative options for the interpolation*)
 SolveAndScaleOptions -> {}, (*Additional/alternative options to the core solver*)
 \[Tau]hOptions -> {}, (*Additional/alternative options to pass to computing \[Tau]h*)
 \[Tau]hExistenceOptions -> {}, (*Additional/alternative options to pass to computing \[Tau]h when only testing it's existence*)
 NumLimitOptions -> {} (*Additional/alternative options to finding the numerical limits*)
 }
 ComputeTunedThermo[type_String, pots_List, param_?NumericQ, umin_?NumericQ, umid_?NumericQ, umax_?NumericQ, opts : OptionsPattern[]] := 
  (
	PrintV[StringForm["Computing `1` thermo, param = `2`", type, param], "Progress"];
	ComputeThermoCurve[ntfun, \[Lambda]hfun, tachyonic, umap[umin], umap[umid], umap[umax], pots,
   NumLimitOptions -> Join[OptionValue[NumLimitOptions], NumLimitDefaults],
   SolveAndScaleOptions -> Join[OptionValue[SolveAndScaleOptions], SolveAndScaleDefaults],
   InterpolationOptions -> Join[OptionValue[InterpolationOptions], InterpolationDefaults],
   opts,
   IndicatorFunction -> IndicatorDefaults,
   Interpolate\[Lambda]h -> Interpolate\[Lambda]hDefault,
   Interpolatent -> InterpolatentDefault
 ] /. ThermoTune[type, param, pots])


ComputeAndSaveThermo::nonascii = "Warning: The filename `1` contains non-ASCII characters, and may be invalid on some filesystems!";
Options[ComputeAndSaveThermo] = Options[ComputeTunedThermo]
ComputeAndSaveThermo[type_String, pots_List, param_?NumericQ, umin_?NumericQ, umid_?NumericQ, umax_?NumericQ, opts : OptionsPattern[]] :=
Module[{nums, hash, filename},
 Block[{$vcontext = "ComputeAndSaveThermo"},
 nums = ComputeTunedThermo[type, pots, param, umin, umid, umax, opts];
 Block[{CurveType = type, Param = param, Potentials = pots, Results, \[Lambda]hFunction, ntFunction, \[Tau]hFunction},
	(*Compress the numerical results, as this typically save O(80%) space *)
	Results = Compress[nums];
	(*Now that everything is named equally, we need to construct unique filenames. Use a SHA160 hash.*)
	hash = Hash[ToString[{Results, pots, param}]]; (*There's a peculiar bug that Hash[...] doesn't seem to evaluate occasionally, possibly related to running in a parallel kernel.  Hopefully the ToString heals it...*)
	filename = ToString[StringForm["`1``2`.m", StringReplace[type, {"\[Lambda]" -> "la", "\[Tau]" -> "tau"}], hash]];
	If[!StringFreeQ[filename, RegularExpression["[^[:ascii:]]"]], Message[ComputeAndSaveThermo::nonascii, filename]];
	PrintV[StringForm["Saving results to `1`\.1d", filename], "Progress"];
	Quiet[DeleteFile[filename]]; (*Delete the file if it already exists. No space for safety here, as we typically want to run in batches.*)

	Save[filename, CurveType];
	Save[filename, Param];
	Save[filename, Potentials];
	PrintV[StringForm["Result size `1` bytes, compressed to `2` bytes", ByteCount[nums], ByteCount[Results]], "Debug"];
	PrintV[StringForm["Results: `1`", nums], "All"];
	Save[filename, Results];
	\[Lambda]hFunction = \[Lambda]hfun /. ThermoTune[type, param, pots];
	ntFunction = ntfun /. ThermoTune[type, param, pots];
	\[Tau]hFunction = If[!tachyonic, (0&), \[Tau]hfun] /. ThermoTune[type, param, pots]; (*Currently we do not pass the actual \[Tau]h function,
																		so just use 0 for non-tachyonic and a dummy for tachyonic.*)
	Save[filename, \[Lambda]hFunction];
	Save[filename, ntFunction];
	Save[filename, \[Tau]hFunction];
 ];
 PrintV[StringForm["Done `1` at `2`.", type, param], "Progress"];
 ];
 filename
]


Options[ComputeStandardThermo] = {ntcurves -> 30,
							       \[Lambda]hcurves -> 30,
								NonTachyonicOnly -> False,
							       nt\[Tau]hcurves -> 20,
							      \[Lambda]h\[Tau]hcurves -> 20,
								\[Lambda]h\[Tau]hrange -> {0.01, 10^5},
								ExtraSetup :>  Null, 
								Parallel -> True,
								\[Lambda]hOptions -> {},
								ntOptions -> {},
								\[Lambda]h\[Tau]hOptions -> {},
								nt\[Tau]hOptions -> {},
								\[Lambda]h\[Tau]hmargin -> 10^-3
								}

ComputeStandardThermo[pots_List, destinationDirectory_String, opts: OptionsPattern[]] := Module[{\[Lambda]hClist, ntClist, \[Lambda]h\[Tau]hClist, nt\[Tau]hClist, func, comps, \[Lambda]endnt0, list},
Block[{$vcontext = "ComputeStandardThermo"}, 

(*Symmetric phase*)
ntClist = ParamListToComputationList["ntconstant", MakentconstLists["NonTachyonic", pots, OptionValue[ntcurves]], pots, OptionValue[\[Lambda]hOptions]];
\[Lambda]hClist  = ParamListToComputationList["\[Lambda]hconstant", Make\[Lambda]hconstLists["NonTachyonic", pots, OptionValue[\[Lambda]hcurves]], pots, OptionValue[ntOptions]];

list = Join[ntClist, \[Lambda]hClist];

If[!OptionValue[NonTachyonicOnly],
(*Tachyonic phase*)
nt\[Tau]hClist = ParamListToComputationList["ntconstant\[Tau]h", MakentconstLists["Tachyonic", pots, OptionValue[nt\[Tau]hcurves], \[Lambda]hrange-> OptionValue[\[Lambda]h\[Tau]hrange]], pots, OptionValue[\[Lambda]h\[Tau]hOptions]];
(*There's no clear way to get the lower limit where to start the \[Lambda]h const curves from. Heuristically, we just find the num limit for nt = 0, and use that*)
PrintV["Computing \[Lambda]end(nt = 0).", "Progress"];
\[Lambda]endnt0 = FindNumLimit[\[Tau]hFromQuarkMass[0, \[Lambda]h, 0, pots, AccuracyGoal -> 60, NodeCountSubdivision-> 0, ARange -> 120, MassAUV -> 110, MaxRecursion -> 80], {\[Lambda]h, 0.01, 100.}, MaxRecursion -> 120];
PrintV[StringForm["\[Lambda]end(nt = 0) = `1`", \[Lambda]endnt0], "Progress"];
\[Lambda]h\[Tau]hClist = ParamListToComputationList["\[Lambda]hconstant\[Tau]h", Make\[Lambda]hconstLists["Tachyonic", pots, OptionValue[\[Lambda]h\[Tau]hcurves], \[Lambda]hrange-> ({Min[#], Max[#]}&@IntervalIntersection[Interval[OptionValue[\[Lambda]h\[Tau]hrange]], Interval[{(1+ OptionValue[\[Lambda]h\[Tau]hmargin])\[Lambda]endnt0, Infinity}]])], pots, OptionValue[nt\[Tau]hOptions]];
list = Join[list, nt\[Tau]hClist, \[Lambda]h\[Tau]hClist];
];

PrintV[StringForm["Generated lists `1`, `2`, `3` and `4`", ntClist, \[Lambda]hClist, nt\[Tau]hClist, \[Lambda]h\[Tau]hClist], "Checks"];

(*Set up parallel or non-parallel computation*)
If[OptionValue[Parallel],
(
func = Composition[ParallelSubmit, (SetDirectory[destinationDirectory];OptionValue[ExtraSetup]; ComputeAndSaveThermo[Sequence @@ #]; ResetDirectory[])&];
ParallelNeeds["VQCDThermo`"];
DistributeDefinitions[pots];
DistributeDefinitions[Options[ComputeStandardThermo]];(*This needs to be done if any of the options to this function are to be used in the parallel kernels.*)
),
(
func = ComputeTunedThermo;
)
];

(*TODO: set up \[Tau]hcurve computations.*)

(*Start the computation*)
comps = Map[func, list];
(*Wait for them to finish*)
WaitAll[comps];

PrintV["Done!", "Progress"];
]
]


LoadThermo[filename_String] := 
Module[{nums},
 Block[{CurveType, Param, Potentials, Results, \[Lambda]hFunction, ntFunction, \[Tau]hFunction},
 Get[filename];
 nums = Uncompress[Results];
 {CurveType, Param, Potentials, nums, \[Lambda]hFunction, ntFunction, \[Tau]hFunction}
]
]


LoadThermoDirectory[directory_String] := (Print[StringForm["Loading file `1`", #]];LoadThermo[#])& /@ FileNames["*.m", {directory}];


PlotPointsAndInterpolation[in : {data_, {fscaleidx_?NumberQ, \[CapitalLambda]idx_?NumberQ, \[Mu]idx_?NumberQ, \[Lambda]hidx_, ntidx_, \[Tau]hidx_}}] := 
Module[{\[CapitalLambda]points, fscalepoints, \[Mu]points, \[CapitalLambda]intp, fsintp, \[Mu]intp, \[CapitalLambda]plot, fplot, \[Mu]plot, \[Tau]hpoints, \[Tau]hintp, \[Tau]hplot},
	\[CapitalLambda]points = \[CapitalLambda]PointsFromData[in];
	fscalepoints = fScalePointsFromData[in];
	\[Mu]points = \[Mu]PointsFromData[in];
	\[CapitalLambda]intp = Interpolation[\[CapitalLambda]points];
	\[CapitalLambda]plot = Show[ListPlot[\[CapitalLambda]points, AxesLabel -> {"u", "\[CapitalLambda]"}], Plot[\[CapitalLambda]intp[x], {x, \[CapitalLambda]points[[1, 1]], \[CapitalLambda]points[[-1, 1]]}, PlotRange -> All]];
	fsintp = Interpolation[fscalepoints];
	fplot = Show[ListPlot[fscalepoints, AxesLabel -> {"u", "fscale"}, PlotRange -> All], Plot[fsintp[x], {x, fscalepoints[[1, 1]], fscalepoints[[-1, 1]]}, PlotRange -> All]];
	\[Mu]intp = Interpolation[\[Mu]points];
	\[Mu]plot = Show[ListPlot[\[Mu]points, AxesLabel -> {"u", "\[Mu]"}, PlotRange -> All], Plot[\[Mu]intp[x], {x, \[Mu]points[[1, 1]], \[Mu]points[[-1, 1]]}, PlotRange-> All]];
	Print[\[CapitalLambda]plot];
	Print[fplot];
	Print[\[Mu]plot];
	If[!(\[Tau]hidx === Undefined), 
		\[Tau]hpoints = \[Tau]hPointsFromData[in];
		\[Tau]hintp = Interpolation[\[Tau]hpoints];
		\[Tau]hplot = Show[ListPlot[\[Tau]hpoints, AxesLabel -> {"u", "\[Tau]h"}, PlotRange -> All], Plot[\[Tau]hintp[u], {u, \[Tau]hpoints[[1, 1]], \[Tau]hpoints[[-1, 1]]}, PlotRange-> All]];
		Print[\[Tau]hplot];
	];
	{\[CapitalLambda]plot, fplot, \[Mu]plot, \[Tau]hplot}
]


PlotPointsAndInterpolation[x : {ct_String, param_?NumericQ, _List, data_, ___}] := (Print[StringForm["Printing curvetype `1`, param `2`", ct, param]];PlotPointsAndInterpolation[data]);


PlotPointsAndInterpolation[file_String] := PlotPointsAndInterpolation[LoadThermo[file]];


PlotPointsAndInterpolationDirectory[dir_String] := PlotPointsAndInterpolation /@ LoadThermoDirectory[dir];


fScalePointsFromData[{data_, {fscaleidx_?NumberQ, ___}}] := {#[[1]], #[[fscaleidx]]}& /@ data;


\[CapitalLambda]PointsFromData[{data_, {_, \[CapitalLambda]idx_?NumberQ, ___}}] := {#[[1]], #[[\[CapitalLambda]idx]]}& /@ data;


\[Mu]PointsFromData[{data_, {_, _, \[Mu]idx_?NumberQ, ___}}] := {#[[1]], #[[\[Mu]idx]]}& /@ data;


\[Tau]hPointsFromData[{data_, {_, _, _, _, _, \[Tau]hidx_?NumberQ}}] := {#[[1]], #[[\[Tau]hidx]]}& /@ data;


End[ ]
EndPackage[ ]
