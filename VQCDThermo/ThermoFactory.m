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
\[Tau]hPointsFromData::usage = "Returns the \[Tau]h datapoints"
\[Lambda]hPointsFromData::usage = "Returns the \[Lambda]h datapoints"
ntPointsFromData::usage = "Returns the nt datapoints"
paramGridFromData::usage = "Returns the parameter values along which the data is defined"
\[CapitalLambda]funFromData::usage = "Returns \[CapitalLambda] as a function"
fScalefunFromData::usage = "Returns fScale as a function"
\[Mu]funFromData::usage = "Returns \[Mu] as a function"
\[Tau]hfunFromData::usage = "Returns \[Tau]h as a function"
\[Lambda]hfunFromData::usage = "Returns \[Lambda]h as a function"
ntfunFromData::usage = "Returns nt as a function"
PlotPointsAndInterpolation::usage = "Given data, prints diagnostic plots with both the data points and interpolations"
PlotPointsAndInterpolationDirectory::usage = "Prints diagnostic plots for for all data in a given directory"

ComputeThermoCurve::usage = "Computes the thermo along a curve, finding the limits automatically etc."
ComputeTunedThermo::usage = "Given a thermo tuning, computes thermo"
ComputeAndSaveThermo::usage = "Computes and saves a thermo, given tuning"

ComputeStandardThermo::usage = "The most straightforward way to compute the thermo. Attempts to automatically analyze your potentials etc."

LoadThermo::usage = ""
LoadThermoDirectory::usage = ""
DataFromThermo::usage = "Returns the raw data, given thermo loaded from file"

MakeThermoFunctions::usage = "Given the basic thermo data, builds the parametric thermo functions and returns them in an associative array"
MakeThermoGrid::usage = "Constructs the thermogrid used for the analysis, given a list containing the standard set of thermocurves, as computed for example by ComputeStandardThermo"
ListThermo::usage = "A very simple helper for quickly checking which curves are present in a list of thermo results"

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

\[Lambda]h;
nt;
\[Tau]h;
T;
p0;
\[Mu];
fScale;
\[CapitalLambda];
n;
min;
max;
lattice;
proots;

(*Protect the above variables, since accidentally overwriting them could prevent access
to the associative array returned by the thermo analysis functions.*)
(Protect[#];)&/@{\[Lambda]h,
nt,
\[Tau]h,
T,
p0,
p,
\[Mu],
fScale,
\[CapitalLambda],
n,
min,
max,
lattice,
proots};

Begin["`Private`"]


Options[ParamListToComputationList] = Options[ComputeTunedThermo]
ParamListToComputationList[curveType_String, list_List, pots_List, opts : OptionsPattern[]] := {curveType, pots, #[[1]], #[[2]], #[[3]], #[[4]], opts}&/@ Sort[list, #1[[1]] < #2[[1]]&];


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
	(*Note that the condition at the minimum is ntfun' < 0, and at the maximum ntfun' > 0*)
	roots = Join[roots, {{ntfun[#], #, ntfun'[#] < 0}& @ First @ OptionValue[\[Lambda]hrange], 
						 {ntfun[#], #, ntfun'[#] > 0}& @ Last @ OptionValue[\[Lambda]hrange]}];
	FindConsecutiveRegions[0, Sequence @@ OptionValue[\[Lambda]hrange], roots]

]


Options[Make\[Lambda]hRangesFromnt] = Options[FindThermoComputationBoxes]
Options[MakentconstLists] = Join[Options[Make\[Lambda]hRangesFromnt], {SpreadFunction -> (ConditionalExpression[#^2, # >= 0]&), MinCompBoxHeight -> 0.05, ntrange -> {0, Infinity}, ntlist -> Automatic}]

Make\[Lambda]hRangesFromnt[n_?NumericQ, compboxes_List?(ListQ[#[[1]]]&), opts : OptionsPattern[]] := {n, #[[3]], #[[4]], #[[5]]}& /@ Select[compboxes, #[[1]] <= n < #[[2]]&];

Make\[Lambda]hRangesFromnt[n_?NumericQ, pots_List, opts : OptionsPattern[]] := Module[{compboxes}, compboxes = FindThermoComputationBoxes[pots, Evaluate[FilterRules[{opts},
	Options[FindThermoComputationBoxes]]]];
	Make\[Lambda]hRangesFromnt[n, compboxes]];

MakentconstLists[compboxes_List?(ListQ[#[[1]]]&), nCurves_Integer, opts : OptionsPattern[]] := Module[{nmax, nvals, ifun},
Block[{$vcontext = "MakentconstLists"},
nmax = Min[Max[#[[2]]&/@compboxes], Last[OptionValue[ntrange]]];

PrintV[StringForm["compboxes `1`, maximum `2`", compboxes, nmax], "Debug"];

(*Build the inverse function, simplifying with the knowledge nt > 0*)
ifun[nt_] = Simplify[InverseFunction[OptionValue[SpreadFunction]][nt], nt > 0];

nvals = If[OptionValue[ntlist] === Automatic,
	(*No explicit list of values given, so generate a table with the requested number of curves*)
	Table[OptionValue[SpreadFunction] @ n, {n, ifun @ First[OptionValue[ntrange]], 
ifun @ nmax, ((#[[2]] - #[[1]])/(nCurves + 1))& @ Thread @ ifun[{0, nmax}]}] (*+1 since the limit case never gets included*)
	,
	(*An explicit list was given, so use that. Note that both nCurves and ntrange lower limit are ignored!*)
	Select[OptionValue[ntlist], # <=  nmax&]
	];
Flatten[Make\[Lambda]hRangesFromnt[#, compboxes,  Evaluate[FilterRules[{opts},
	Options[Make\[Lambda]hRangesFromnt]]]]&/@nvals,1]
]
]

MakentconstLists[tachyon : "Tachyonic" | "NonTachyonic", pots_List, nCurves_Integer, opts : OptionsPattern[]] := Module[{compboxes, nmax, nvals},
compboxes = If[tachyon === "NonTachyonic", FindThermoComputationBoxes[pots, Evaluate[FilterRules[{opts},
	Options[FindThermoComputationBoxes]]]], FindThermoTachyonicComputationBoxes[pots, Evaluate[FilterRules[{opts},
	Options[FindThermoComputationBoxes]]]]];
MakentconstLists[Select[compboxes, (#[[5]] - #[[3]] >=  OptionValue[MinCompBoxHeight])&], nCurves, opts]
]



Options[Make\[Lambda]hconstLists] = {\[Lambda]hrange -> {0.01, 10^5},
	\[Lambda]hlist -> Automatic,
	ntrange -> {0, 100},
	\[Lambda]hmargin -> 10^-3, (*Points _very _ close to the limit usually behave badly numerically.*)
	SpreadFunction -> (Exp) (*The function with which to spread the points. Identity gives a linear spread*)}

Make\[Lambda]hconstLists["Tachyonic", pots_List, nCurves_Integer, opts : OptionsPattern[]] := Module[{Veff, ntfun, \[Lambda]hvals, \[Lambda]hnummax, \[Lambda]hr},
ntfun[\[Lambda]h_] = Limit[ntCritical[pots, \[Tau]h][\[Lambda]h], \[Tau]h -> Infinity];

\[Lambda]hr = OptionValue[\[Lambda]hrange];
\[Lambda]hvals = If[OptionValue[\[Lambda]hlist] === Automatic,
	(*No explicit list of points given, generate it*)
	Table[OptionValue[SpreadFunction] @ i, Evaluate[{i,Sequence @@ InverseFunction[OptionValue[SpreadFunction]] @ \[Lambda]hr, ((#[[2]] - #[[1]])/nCurves)& @  InverseFunction[OptionValue[SpreadFunction]] [\[Lambda]hr]}]]
	,
	(*Use the explicit list, but limit it to the specified range*)
	Select[OptionValue[\[Lambda]hlist], First[\[Lambda]hr] <=  # <= Last[\[Lambda]hr]&]
];

{#, First[OptionValue[ntrange]], First[OptionValue[ntrange]], Min[ntfun[#], Last[OptionValue[ntrange]]]}& /@ \[Lambda]hvals
]

Make\[Lambda]hconstLists["NonTachyonic", pots_List, nCurves_Integer, opts : OptionsPattern[]] := Module[{Veff, ntfun, \[Lambda]hvals, \[Lambda]hnummax, \[Lambda]hr},
Veff[\[Lambda]h_, nt_, \[Tau]h_]= VEffective[{Vg, Vf, \[Kappa], \[Omega]}][\[Lambda]h, nt, \[Tau]h];
ntfun[\[Lambda]h_] = (nt /. Last[Quiet[Solve[Veff[\[Lambda]h, nt, 0] == 0, nt]]])^2 /. {Vf -> pots[[2]], Vg -> pots[[1]], \[Kappa]-> pots[[3]], \[Omega] -> pots[[4]]};

(*Compute the maximum \[Lambda]h by using the fact that at nt = 0 it is precisely when Veff' = 0*)
\[Lambda]hnummax = (\[Lambda]h(1 - OptionValue[\[Lambda]hmargin]))/. FindRoot[Derivative[1, 0, 0][Veff][\[Lambda]h, 0, 0] /. {Vf -> pots[[2]], Vg -> pots[[1]], \[Kappa]-> pots[[3]], \[Omega] -> pots[[4]]}, {\[Lambda]h, Sequence@@OptionValue[\[Lambda]hrange]}];

\[Lambda]hr = {First[OptionValue[\[Lambda]hrange]], \[Lambda]hnummax};

\[Lambda]hvals = If[OptionValue[\[Lambda]hlist] === Automatic,
	(*No explicit list given, generate the table*)
	Table[OptionValue[SpreadFunction] @ i, Evaluate[{i,Sequence @@ InverseFunction[OptionValue[SpreadFunction]] @ \[Lambda]hr, ((#[[2]] - #[[1]])/nCurves)& @  InverseFunction[OptionValue[SpreadFunction]] [\[Lambda]hr]}]]
	,
	Select[OptionValue[\[Lambda]hlist], First[\[Lambda]hr] <= # <= Last[\[Lambda]hr]&]
	];

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
qh = qhFromThermoData[fscale, \[Lambda]hval, \[Tau]hval, ntval, pots];
mu = AAndMuFromSols[sol, ntval, pots, AccuracyGoal -> OptionValue[AIntegralAccuracy], IntegrationGrid -> Head[\[Lambda]FromSols[sol][A]]][[2]];
PrintV[StringForm["Thermo done, indicatorfun = `1`", OptionValue[IndicatorFunction][fscale, \[CapitalLambda], qh, \[Lambda]hval, mu, \[Lambda]hval, ntval, \[Tau]hval]], "All"];
If[NumericQ[fscale], 
	(Sow[outValues[u, fscale, \[CapitalLambda], mu, \[Lambda]hval, ntval, \[Tau]hval]];OptionValue[IndicatorFunction][fscale, \[CapitalLambda], qh, mu, \[Lambda]hval, ntval, \[Tau]hval]),
	(*Not numeric, but check if u is numeric either. If not, then FunctionInterpolation is just trying to simplify the 
	expression before the actual evaluation.*)
   If[NumericQ[u], 
	(If[unummax>u, unummax = u]; If[unummin < u, unummin = u];PrintV[StringForm["Non-numeric value `1` at u = `5`, \[Lambda]h = `2`, nt = `3`, \[Tau]h = `4`", fscale, \[Lambda]hval, ntval, \[Tau]hval, u], "Progress"]; 0),
	fscale
	]]
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
 ComputeThermoCurve[ntfun_, \[Lambda]hfun_, tachyonic_?((# == True) || (# == False)&), umin_?NumericQ, umid0_?NumericQ, umax_?NumericQ, pots_List, opts : OptionsPattern[]] := Module[
  {numres, uminlim, umaxlim, \[Tau]hFun, time, existfun, umid},
 If[!tachyonic,
 (\[Tau]hFun[\[Lambda]h_, nt_, existenceOnly_ : False] = 0;)
 ,
 (
  \[Tau]hFun[\[Lambda]h_?NumericQ, nt_?NumericQ, existenceOnly_ : False] = \[Tau]hFromQuarkMass[0, \[Lambda]h, nt, pots, Evaluate[If[existenceOnly, Options[\[Tau]hExistenceOptions], Options[\[Tau]hOptions]]], OptionValue[SolveAndScaleOptions], AccuracyGoal -> 60, NodeCountSubdivision-> 0, ARange -> 60, MassAUV -> 60]
 )
 ];

 existfun[u_?NumericQ] := Module[{\[Tau]h = \[Tau]hFun[\[Lambda]hfun[u], ntfun[u], True]},
	If[NumericQ[\[Tau]h],
		\[CapitalLambda]ScaleFromSols[SolveAndScaleVQCDBH[\[Lambda]hfun[u], \[Tau]h, ntfun[u], pots, Evaluate[OptionValue[SolveAndScaleOptions]]]],
		Undefined (*Return undefined if even \[Tau]h is not defined*)
	]
	];

 time = First[AbsoluteTiming[

 umid = If[!(NumericQ[existfun[umid0]] === True),
 (*(Message[ComputentConstThermo::nosolution, umid]; Return[Undefined])*)
 (*The solution does not exist at the midpoint. However, try to see if it exists at either the upper or lower limit, and use that as the new midpoint:*)
  Module[{existmax, existmin},
	existmax = NumericQ[existfun[umax]] === True;
	
	If[existmax,
		umax,
		(existmin = NumericQ[existfun[umin]] === True;
		 If[existmin,
			umin,
			(Message[ComputentConstThermo::nosolution, umid]; Return[Undefined])
		])	
	  ]
	],
	umid0
 ];

 If[Quiet[!(NumericQ[existfun[umin]] === True)], uminlim = FindNumLimit[Quiet[existfun[u]], {u, umin, umid}, Evaluate[OptionValue[NumLimitOptions]]
], uminlim = umin];
 PrintV[StringForm["Found u lower limit `1`", uminlim], "Progress"];

 If[Quiet[!(NumericQ[existfun[umax]] === True)], umaxlim = FindNumLimit[Quiet[existfun[u]], {u, umid, umax}, Evaluate[OptionValue[NumLimitOptions]]], umaxlim = umax];
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

 PrintV[StringForm["Done `6`curve from {nt, \[Lambda]h} = {`2`, `3`} to {`4`, `5`} in `1` seconds",
		time,
		If[OptionValue[Interpolatent],
			ntPointsFromData[numres][[1, 2]],
			ntfun[umin]
		],
		If[OptionValue[Interpolate\[Lambda]h],
			\[Lambda]hPointsFromData[numres][[1, 2]],
			\[Lambda]hfun[umin]
		],
		If[OptionValue[Interpolatent],
			ntPointsFromData[numres][[-1, 2]],
			ntfun[umax]
		],
		If[OptionValue[Interpolate\[Lambda]h],
			\[Lambda]hPointsFromData[numres][[-1, 2]],
			\[Lambda]hfun[umax]
		],
		If[tachyonic, "tachyonic ", ""]], "Progress"];
 
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
								\[Lambda]hlist -> Automatic,
								ntlist -> Automatic,
								NonTachyonicOnly -> False,
							       nt\[Tau]hcurves -> 20,
							      \[Lambda]h\[Tau]hcurves -> 20,
								\[Lambda]h\[Tau]hrange -> {0.01, 10^5},
								nt\[Tau]hrange -> {0, Infinity},
								\[Lambda]h\[Tau]hlist -> Automatic,
								nt\[Tau]hlist -> Automatic,
								ExtraSetup :>  Null, 
								Parallel -> True,
								\[Lambda]hOptions -> {},
								ntOptions -> {},
								\[Lambda]h\[Tau]hOptions -> {},
								nt\[Tau]hOptions -> {},
								\[Lambda]h\[Tau]hmargin -> 10^-3,
								\[Lambda]h\[Tau]hlowerlimit -> Automatic
								}

ComputeStandardThermo[pots_List, destinationDirectory_String, opts: OptionsPattern[]] := Module[{\[Lambda]hClist, ntClist, \[Lambda]h\[Tau]hClist, nt\[Tau]hClist, func, comps, \[Lambda]endnt0, list},
Block[{$vcontext = "ComputeStandardThermo"}, 

(*Symmetric phase*)
ntClist = ParamListToComputationList["ntconstant", MakentconstLists["NonTachyonic", pots, OptionValue[ntcurves], ntlist -> OptionValue[ntlist]], pots, OptionValue[\[Lambda]hOptions]];
PrintV[StringForm["Symmetric nt -values to compute: `1`", ntClist[[All, 3]]], "Progress"];
\[Lambda]hClist  = ParamListToComputationList["\[Lambda]hconstant", Make\[Lambda]hconstLists["NonTachyonic", pots, OptionValue[\[Lambda]hcurves], \[Lambda]hlist -> OptionValue[\[Lambda]hlist]], pots, OptionValue[ntOptions]];
PrintV[StringForm["Symmetric \[Lambda]h -values to compute: `1`", \[Lambda]hClist[[All, 3]]], "Progress"];

list = Join[ntClist, \[Lambda]hClist];

If[!OptionValue[NonTachyonicOnly],
(*Tachyonic phase*)
nt\[Tau]hClist = ParamListToComputationList["ntconstant\[Tau]h", MakentconstLists["Tachyonic", pots, OptionValue[nt\[Tau]hcurves], \[Lambda]hrange-> OptionValue[\[Lambda]h\[Tau]hrange], ntrange -> OptionValue[nt\[Tau]hrange], ntlist -> OptionValue[nt\[Tau]hlist]], pots, OptionValue[\[Lambda]h\[Tau]hOptions]];
PrintV[StringForm["Broken nt -values to compute: `1`", nt\[Tau]hClist[[All, 3]]], "Progress"];
(*There's no clear way to get the lower limit where to start the \[Lambda]h const curves from. Heuristically, we just find the num limit for nt = 0, and use that*)
If[(OptionValue[\[Lambda]h\[Tau]hlowerlimit] === Automatic) && !(OptionValue[\[Lambda]h\[Tau]hlist] === {}),
(
	PrintV["Computing \[Lambda]end(nt = 0).", "Progress"];
	\[Lambda]endnt0 = FindNumLimit[\[Tau]hFromQuarkMass[0, \[Lambda]h, 0, pots, AccuracyGoal -> 60, NodeCountSubdivision-> 0, ARange -> 120, MassAUV -> 110, MaxRecursion -> 80], {\[Lambda]h, 0.01, 100.}, MaxRecursion -> 120];
),
	\[Lambda]endnt0 = OptionValue[\[Lambda]h\[Tau]hlowerlimit];
]
	
PrintV[StringForm["\[Lambda]end(nt = 0) = `1`", \[Lambda]endnt0], "Progress"];
\[Lambda]h\[Tau]hClist = ParamListToComputationList["\[Lambda]hconstant\[Tau]h", Make\[Lambda]hconstLists["Tachyonic", pots, OptionValue[\[Lambda]h\[Tau]hcurves], \[Lambda]hrange-> ({Min[#], Max[#]}&@IntervalIntersection[Interval[OptionValue[\[Lambda]h\[Tau]hrange]], Interval[{(1+ OptionValue[\[Lambda]h\[Tau]hmargin])\[Lambda]endnt0, Infinity}]]), ntrange -> OptionValue[nt\[Tau]hrange], \[Lambda]hlist -> OptionValue[\[Lambda]h\[Tau]hlist]], pots, OptionValue[nt\[Tau]hOptions]];
PrintV[StringForm["Broken \[Lambda]h -values to compute: `1`", \[Lambda]h\[Tau]hClist[[All, 3]]], "Progress"];
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
func = (SetDirectory[destinationDirectory];OptionValue[ExtraSetup]; ComputeAndSaveThermo[Sequence @@ #]; ResetDirectory[])&;
)
];

(*TODO: set up \[Tau]hcurve computations.*)

(*Start the computation*)
comps = Map[func, list];
(*Wait for them to finish*)
If[OptionValue[Parallel],
	WaitAll[comps];
]

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


LoadThermoDirectory[directory_String] := (PrintV[StringForm["Loading file `1`", #], "Progress"];LoadThermo[#])& /@ FileNames["*.m", {directory}];


DataFromThermo[{ct_, param_, pots_, data_, \[Lambda]hfun_, ntfun_, \[Tau]hfun_}] := data;


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


\[Lambda]hPointsFromData[{data_, {_, _, _, \[Lambda]hidx_?NumberQ, ___}}] := {#[[1]], #[[\[Lambda]hidx]]}& /@ data;


ntPointsFromData[{data_, {_, _, _, _, ntidx_?NumberQ, _}}] := {#[[1]], #[[ntidx]]}& /@ data;


(*Convenience wrappers for the above*)
\[CapitalLambda]PointsFromData[{ct_String, param_?NumericQ, pots_List, points : {data_, {fscaleidx_?NumberQ, \[CapitalLambda]idx_?NumberQ, \[Mu]idx_?NumberQ, \[Lambda]hidx_, ntidx_, \[Tau]hidx_}}, \[Lambda]hFunction_, ntFunction_, \[Tau]hFunction_}] :=
	\[CapitalLambda]PointsFromData[points];
fScalePointsFromData[{ct_String, param_?NumericQ, pots_List, points : {data_, {fscaleidx_?NumberQ, \[CapitalLambda]idx_?NumberQ, \[Mu]idx_?NumberQ, \[Lambda]hidx_, ntidx_, \[Tau]hidx_}}, \[Lambda]hFunction_, ntFunction_, \[Tau]hFunction_}] :=
	fScalePointsFromData[points];
\[Mu]PointsFromData[{ct_String, param_?NumericQ, pots_List, points : {data_, {fscaleidx_?NumberQ, \[CapitalLambda]idx_?NumberQ, \[Mu]idx_?NumberQ, \[Lambda]hidx_, ntidx_, \[Tau]hidx_}}, \[Lambda]hFunction_, ntFunction_, \[Tau]hFunction_}] :=
	\[Mu]PointsFromData[points];
\[Tau]hPointsFromData[{ct_String, param_?NumericQ, pots_List, points : {data_, {fscaleidx_?NumberQ, \[CapitalLambda]idx_?NumberQ, \[Mu]idx_?NumberQ, \[Lambda]hidx_, ntidx_, \[Tau]hidx_?NumberQ}}, \[Lambda]hFunction_, ntFunction_, \[Tau]hFunction_}] :=
		\[Tau]hPointsFromData[points]
\[Lambda]hPointsFromData[{ct_String, param_?NumericQ, pots_List, points : {data_, {fscaleidx_?NumberQ, \[CapitalLambda]idx_?NumberQ, \[Mu]idx_?NumberQ, \[Lambda]hidx_?NumberQ, ntidx_, \[Tau]hidx_}}, \[Lambda]hFunction_, ntFunction_, \[Tau]hFunction_}] :=
		\[Lambda]hPointsFromData[points]
ntPointsFromData[{ct_String, param_?NumericQ, pots_List, points : {data_, {fscaleidx_?NumberQ, \[CapitalLambda]idx_?NumberQ, \[Mu]idx_?NumberQ, \[Lambda]hidx_, ntidx_?NumberQ, \[Tau]hidx_}}, \[Lambda]hFunction_, ntFunction_, \[Tau]hFunction_}] :=
		ntPointsFromData[points]


(*Convenience wrapper for extracting the basic functions*)
\[CapitalLambda]funFromData[{ct_String, param_?NumericQ, pots_List, points : {data_, {fscaleidx_?NumberQ, \[CapitalLambda]idx_?NumberQ, \[Mu]idx_?NumberQ, \[Lambda]hidx_, ntidx_, \[Tau]hidx_}}, \[Lambda]hFunction_, ntFunction_, \[Tau]hFunction_}] := 
 Interpolation[\[CapitalLambda]PointsFromData[points]];
fScalefunFromData[{ct_String, param_?NumericQ, pots_List, points : {data_, {fscaleidx_?NumberQ, \[CapitalLambda]idx_?NumberQ, \[Mu]idx_?NumberQ, \[Lambda]hidx_, ntidx_, \[Tau]hidx_}}, \[Lambda]hFunction_, ntFunction_, \[Tau]hFunction_}] := 
 Interpolation[fScalePointsFromData[points]];
\[Mu]funFromData[{ct_String, param_?NumericQ, pots_List, points : {data_, {fscaleidx_?NumberQ, \[CapitalLambda]idx_?NumberQ, \[Mu]idx_?NumberQ, \[Lambda]hidx_, ntidx_, \[Tau]hidx_}}, \[Lambda]hFunction_, ntFunction_, \[Tau]hFunction_}] := 
 Interpolation[\[Mu]PointsFromData[points]];
\[Tau]hfunFromData[{ct_String, param_?NumericQ, pots_List, points : {data_, {fscaleidx_?NumberQ, \[CapitalLambda]idx_?NumberQ, \[Mu]idx_?NumberQ, \[Lambda]hidx_, ntidx_, \[Tau]hidx_}}, \[Lambda]hFunction_, ntFunction_, \[Tau]hFunction_}] := 
	If[\[Tau]hidx === Undefined,
		If[NumericQ @ \[Tau]hFunction[data[[1, 1]]], \[Tau]hFunction, 0&],
		Interpolation[\[Tau]hPointsFromData[points]]
	];
\[Lambda]hfunFromData[{ct_String, param_?NumericQ, pots_List, points : {data_, {fscaleidx_?NumberQ, \[CapitalLambda]idx_?NumberQ, \[Mu]idx_?NumberQ, \[Lambda]hidx_, ntidx_, \[Tau]hidx_}}, \[Lambda]hFunction_, ntFunction_, \[Tau]hFunction_}] := 
	If[\[Lambda]hidx === Undefined,
		\[Lambda]hFunction,
		Interpolation[\[Lambda]hPointsFromData[points]]
	];
ntfunFromData[{ct_String, param_?NumericQ, pots_List, points : {data_, {fscaleidx_?NumberQ, \[CapitalLambda]idx_?NumberQ, \[Mu]idx_?NumberQ, \[Lambda]hidx_, ntidx_, \[Tau]hidx_}}, \[Lambda]hFunction_, ntFunction_, \[Tau]hFunction_}] := 
	If[ntidx === Undefined,
		ntFunction,
		Interpolation[ntPointsFromData[points]]
	];


paramGridFromData[{data_, __}] := First[#]& /@ data;
paramGridFromData[{ct_String, param_?NumericQ, pots_List, points : {data_, {fscaleidx_?NumberQ, \[CapitalLambda]idx_?NumberQ, \[Mu]idx_?NumberQ, \[Lambda]hidx_, ntidx_, \[Tau]hidx_}}, \[Lambda]hFunction_, ntFunction_, \[Tau]hFunction_}] :=
	paramGridFromData[points];


(*The following are the data analysis functions. They might be moved to a separate package
in the future, but as long as it's just two functions, that doesn't seem warranted.*)


Options[MakeThermoFunctions] = {InterpolationDensityFactor -> 8 (*The factor by which to increase the density of the interpolation grid for forming dp, which helps with the fact that dp tends to be less regular
than any of the individual functions.*)}
MakeThermoFunctions[in : {ct_String, param_?NumericQ, pots_List, {data_, {fscaleidx_?NumberQ, \[CapitalLambda]idx_?NumberQ, \[Mu]idx_?NumberQ, \[Lambda]hidx_, ntidx_, \[Tau]hidx_}}, \[Lambda]hFunction_, ntFunction_, \[Tau]hFunction_}, opts: OptionsPattern[]] := Module[{\[Lambda]fun, sfun, fsfun, nfun,ntfun,  \[Mu]fun, \[CapitalLambda]fun, dp, Tfun, \[Tau]fun, pfun, dplist, dpint, grid, proots, res},
Print[StringForm["Making thermo `1` = `2`", ct, param]];
fsfun = fScalefunFromData[in];
\[Mu]fun = \[Mu]funFromData[in];
\[CapitalLambda]fun = \[CapitalLambda]funFromData[in];
\[Lambda]fun = \[Lambda]hfunFromData[in];
\[Tau]fun = \[Tau]hfunFromData[in];
ntfun = ntfunFromData[in];

nfun = (ntfun[#]/\[CapitalLambda]fun[#]^3 / (4 Pi))&;

Tfun = TemperatureFromThermoData[fsfun[#], \[CapitalLambda]fun[#], \[Lambda]fun[#], \[Tau]fun[#], ntfun[#], pots]&;


(*Print[Plot[Tfun[u], {u, data[[1, 1]], data[[-1, 1]]}]];*)

(*The best way this far to compute the integral: form the derivatives at the original grid points and interpolate, which can be integrated analytically.*)
dp = 1/\[CapitalLambda]fun[#]^3 Tfun'[#] + nfun[#] \[Mu]fun'[#]&;
grid = paramGridFromData[DataFromThermo[in]];
pfun = NIntegralFunction[dp[x], {x, Last[grid], First[grid]}, IntegrationGrid -> grid, Densification -> OptionValue[InterpolationDensityFactor]];


(*Print[Show[Plot[{pfun[u]}, {u, data[[1, 1]], data[[-1, 1]]}]]];*)

(*proots = Select[FindFunctionRoots[pfun[x], {x, First[grid], Last[grid]}, NumPoints -> grid], # != Last[grid]&];
Print[StringForm["Roots found: `1`", proots]];*)

res[\[Lambda]h] = \[Lambda]fun;
res[nt] = ntfun;
res[\[Tau]h] = \[Tau]fun;
res[T] = Tfun;
res[p0] = pfun;
res[\[Mu]] = \[Mu]fun;
res[fScale] = fsfun;
res[\[CapitalLambda]] = \[CapitalLambda]fun;
res[n] = nfun;
res[min] = First[grid];
res[max] = Last[grid];
res[lattice] = grid; (*Also return the set of points where the functions have actually been calculated.*)

res

];
MakeThermoGrid[thermo_List, plotThermo : (True | False) : True] := Module[{symmlahfuns, brokelahfuns, symmntfuns, brokentfuns, \[Lambda]hinvert, ntinvert, pressuresFromCrossings, plot, plotnt},

(*Define a helper for inserting the inverses of the parameter functions, which
we can now do since we know we can invert lambdah*)
\[Lambda]hinvert[in_] := Module[{out, inv},
	inv = InverseFunction[in[\[Lambda]h]];
	out[\[Lambda]h][x_] := x;
	out[nt] = in[nt][Null]; (*nt is constant*)

	(*insert the inverse function to the physical variables*)
	(out[#][\[Lambda]h_] := in[#][inv[\[Lambda]h]])&/@ {\[Tau]h, T, p0, \[Mu], fScale, \[CapitalLambda], n};

	out[min] = in[\[Lambda]h][in[min]];
	out[max] = in[\[Lambda]h][in[max]];

	out[lattice] = in[\[Lambda]h]/@ in[lattice];

	out
];

(*Helper function for inserting the necessary inverses to the \[Lambda]h constant curves*)
ntinvert[in_] := Module[{out, inv},
	inv = InverseFunction[in[nt]]; (*this is probably just the identity function,
								but it allows a bit more generality to use this
								here explicitly*)
	out[\[Lambda]h] = in[\[Lambda]h][Null]; (*\[Lambda]h is constant*)
	out[nt][x_] := x; 

	(*insert the inverse function to the physical variables*)
	(out[#][nt_] := in[#][inv[nt]])&/@ {\[Tau]h, T, p0, \[Mu], fScale, \[CapitalLambda], n};

	out[min] = in[nt][in[min]];
	out[max] = in[nt][in[max]];

	out[lattice] = in[nt]/@ in[lattice];

	out
];

(*Helper for plotting the data*)
If[!plotThermo, (plot[in_] := Null),
plot[in_] := Module[{style = If[in[\[Tau]h][in[max]] =!= 0, Blue, Red],
					 Tlimits = {Min[#], Max[#]}&@{in[T][#]& /@ in[lattice]},
					 \[Mu]limits = {Min[#], Max[#]}&@{in[\[Mu]][#]& /@ in[lattice]},
					 plimits = {Min[#], Max[#]}&@{in[p][#]& /@ in[lattice]}},
	Print @ Show[
	(*Set up the canvas*)
	LogPlot[Undefined, Evaluate[{t, Evaluate @ Sequence @@ \[Mu]limits}], 
		PlotRange -> Tlimits,
		AxesLabel -> {"\[Mu]", "T"}],
	(*Plot the interpolation*)
	ParametricPlot[{in[\[Mu]][\[Lambda]h], Log @ in[T][\[Lambda]h]}, {\[Lambda]h, in[min], in[max]},
		AspectRatio -> 1(*, PerformanceGoal -> "Quality"*), MaxRecursion -> 15(*, Method -> {"Refinement" -> {ControlValue -> .01\[Degree]}}*), PlotPoints -> {Automatic, 500}, PlotStyle -> style, PlotRange -> All],
	(*Plot the points*)
	ListPlot[{in[\[Mu]][#], Log @ in[T][#]}& /@ in[lattice], PlotStyle -> style, PlotRange -> All],
	(*Plot the zeroes of pressure*)
	If[in[proots] =!= {},
		ListPlot[{in[\[Mu]][#], Log @ in[T][#]}&/@in[proots], PlotStyle -> Black, PlotMarkers -> {Automatic, Small}],
		{}
		]
	];
	 Print @ Show[
	(*Set up the canvas*)
	LogLinearPlot[Undefined, Evaluate[{t, Evaluate @ Sequence @@ Tlimits}], 
		PlotRange -> plimits,
		AxesLabel -> {"T", "p"}],
		(*Plot the interpolation*)
		ParametricPlot[{Log @ in[T][\[Lambda]h], in[p][\[Lambda]h]}, {\[Lambda]h, in[min], in[max]},
		AspectRatio -> 1, PerformanceGoal -> "Quality", MaxRecursion -> 15(*, Method -> {"Refinement" -> {ControlValue -> .01\[Degree]}}*), PlotPoints -> {500, Automatic}, PlotStyle -> style, PlotRange -> All],
		(*Plot the points*)
		ListPlot[{Log@in[T][#], in[p][#]}& /@ in[lattice], PlotStyle -> style, PlotRange -> All],
		(*Plot the pressure zeroes*)
		If[in[proots] =!= {}, 
			ListPlot[{Log @ in[T][#], in[p][#]}&/@in[proots], PlotStyle -> Black, PlotMarkers -> {Automatic, Small}, PlotRange-> All],
			{}
		]
	];
 ];
];
(*First, compute the basic thermo for the broken phase, nt = 0, since this fixes the zero of the pressure*)
brokelahfuns = Module[{temp},{temp = \[Lambda]hinvert[MakeThermoFunctions[First @ Cases[thermo, {"ntconstant\[Tau]h", 0., _List, _List, _, _, _}]]];
temp[p] = temp[p0];
temp[proots] = Select[FindFunctionRoots[temp[p][x], {x, temp[min], temp[max]}, NumPoints -> temp[lattice]], # != Last[temp[lattice]]&];
PrintV[StringForm["Roots found: `1`", temp[proots]], "Progress"];

temp}];

(*Then, compute the broken phase \[Lambda]h constant functions, since now we can fix the pressure constants*)
brokentfuns = Module[{\[Lambda]hval = #[[2]], pconst, temp},temp = ntinvert[MakeThermoFunctions[#]];
Print[StringForm["Constructing at `1`", #[[2]]]];
pconst = -temp[p0][0.] + brokelahfuns[[1]][p][\[Lambda]hval];
temp[p] = Function[x, temp[p0][x] +  pconst];

temp[proots] = Select[FindFunctionRoots[temp[p][x], {x, temp[min], temp[max]}, NumPoints -> temp[lattice]], # != Last[temp[lattice]]&];
Print[StringForm["Roots found: `1`", temp[proots]]];
 plot[temp];

temp
]&/@ Sort[Cases[thermo, {"\[Lambda]hconstant\[Tau]h", _?NumericQ, _List, _List, _, _, _}], #1[[2]] < #2[[2]]&];

(*Symmetric phase, nt = 0, can be done when we know the broken phase:*)
symmlahfuns = Module[{temp},{(temp = \[Lambda]hinvert[MakeThermoFunctions[First@ Cases[thermo, {"ntconstant", 0., _List, _List, _, _, _}]]];
temp[p] = Function[x, temp[p0][x] - temp[p0][brokelahfuns[[1]][min]] + brokelahfuns[[1]][p][brokelahfuns[[1]][min]]
				];

temp[proots] = Select[FindFunctionRoots[temp[p][x], {x, temp[min], temp[max]}, NumPoints -> temp[lattice]], # != Last[temp[lattice]]&];
PrintV[StringForm["Roots found: `1`", temp[proots]], "Progress"];

temp
)}];

(*symmetric phase \[Lambda]h constant functions, since now we can fix the pressure constants*)
symmntfuns = Module[{\[Lambda]hval = #[[2]], pconst, temp},temp = ntinvert[MakeThermoFunctions[#]];
PrintV[StringForm["Constructing at `1`", #[[2]]], "Progress"];
pconst = -temp[p0][0.] + symmlahfuns[[1]][p][\[Lambda]hval];
temp[p] = Function[x, temp[p0][x] +  pconst];

temp[proots] = Select[FindFunctionRoots[temp[p][x], {x, temp[min], temp[max]}, NumPoints -> temp[lattice]], # != Last[temp[lattice]]&];
PrintV[StringForm["Roots found: `1`", temp[proots]], "Progress"];
 plot[temp];

temp
]&/@ Sort[Cases[thermo, {"\[Lambda]hconstant", _?NumericQ, _List, _List, _, _, _}],#1[[2]] < #2[[2]]&];

(*Rest of the broken phase nt constant functions, since now we can find the crossings to set the constants*)
pressuresFromCrossings[therm_, ntf_] :=
Module[{ntval = therm[[2]], crossings, temp, pconst},
temp = \[Lambda]hinvert[MakeThermoFunctions[therm]];
crossings = {{ntval, #[\[Lambda]h]}, (temp[T][#[\[Lambda]h]] - #[T][ntval])/(temp[T][#[\[Lambda]h]] + #[T][ntval]), #[p][ntval], temp[p0][#[\[Lambda]h]], #[p][ntval] - temp[p0][#[\[Lambda]h]], (#[p][ntval] - temp[p0][#[\[Lambda]h]])/temp[p0][#[\[Lambda]h]]}&/@ Select[ntf,#[min] <= temp[nt] <=  #[max]&];
PrintV[StringForm["Crossingslist, {\[Lambda]h, p, p0, \[CapitalDelta]p, \[CapitalDelta]p/p0}, `1`, mean = `2`, stddeviation = `3`", crossings, Mean[crossings[[All, 4]]], StandardDeviation[crossings[[All, 4]]]], "Checks"];
pconst = Last[crossings][[5]];
temp[p] = Function[x, temp[p0][x] + pconst];

(*Now we can also find zeroes of pressure, which are potential points
for the phase transitions*)
temp[proots] = Select[FindFunctionRoots[temp[p][x], {x, temp[min], temp[max]}, NumPoints -> temp[lattice]], # != Last[temp[lattice]]&];
PrintV[StringForm["Roots found: `1`", temp[proots]], "Progress"];
 plot[temp];

temp
];

brokelahfuns = Join[brokelahfuns,  pressuresFromCrossings[#, brokentfuns]&/@ Sort[Cases[thermo, {"ntconstant\[Tau]h", _?(# != 0. && NumericQ[#]&), _List, _List, _, _, _}], #1[[2]] < #2[[2]]&]];

symmlahfuns = Join[symmlahfuns, pressuresFromCrossings[#, symmntfuns]&/@ Sort[Cases[thermo, {"ntconstant", _?(# != 0. && NumericQ[#]&), _List, _List, _, _, _}],#1[[2]] < #2[[2]]&]];


{symmlahfuns,
symmntfuns,
brokelahfuns,
brokentfuns
}
]
ListThermo[list_] := Print[StringForm["`1` at `2`", #[[1]], #[[2]]]]&/@ list


End[ ]
EndPackage[ ]
