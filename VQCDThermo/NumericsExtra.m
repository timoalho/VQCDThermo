(* ::Package:: *)

BeginPackage["VQCDThermo`NumericsExtra`", {"VQCDThermo`Verbosity`", "DifferentialEquations`InterpolatingFunctionAnatomy`"}]

RemoveImaginaries::usage = "RemoveImaginaries[obj] replaces all complex numbers in obj with their real parts. A warning is issued (but the replacement still takes place), if Im/Re exceeds the bound set by the option Tolerance (default 10^-9)"

ExtractInterpolatingFunctionScalings::usage = "ExtractInterpolatingFunctionScalings[f[z]], where f[z] is a scaled and shifted InterpolatingFunction, i.e. f[z] = scale1 * f0[offset + scale2 * z] for real numbers scale1, scale2 and offset, returns a list {scale1, scale2, offset}"

InterpolatingFunctionToList::usage = "InterpolatingFunctionToList[f[z]] converts an InterpolatingFunction, or a shifted and scaled interpolating function, to a corresponding list of points {z, f[z]}.

If f[z] is of the form f[z] = scale1 * f0[offset + scale2 * z], where f0 is an InterpolatingFunction and scale1, scale2 and offset are numbers, the scalings and offsets are applied to the list.
"

NIntegralFunction::usage = "NIntegralFunction[fun, min, max] returns an InterpolatingFunction which is is the definite integral of fun from min to x, x in the range min to max."

IntegrateWithInterpolatingFunctions::usage = "IntegrateWithInterpolatingFunctions[fun, {x, xmin, xmax}] attempts to integrate the expression fun, replacing any InterpolatingFunction
objects therein by their corresponding interpolating polynomials, integrating the resulting expression analytically, and then doing the piecewise sums. Returns a CompiledFunction object,
which evaluates the resulting new interpolation function.  Note that the piecewise interpolations are no longer necessarily polynomials.

There are at least the following limitations to this process: 
-The expression resulting from replacing each InterpolatingFunction in fun with the corresponding polynomials must be analytically integrable with Integrate.
-any arguments of InterpolatingFunctions in fun must have inverses solvable by Solve[...]. This excludes, among others, nested
InterpolatingFunctions.
"

RemoveImaginaries::LargeIm = "There are large imaginary parts relative to the real part present, Re = `1`, Im = `2`";
RemoveImaginaries::LargestIms = "Largest relative imaginary part `2`, largest absolute imaginary part `1` (which is relatively `3`)"

ExpFunctionInterpolation::usage = "ExpFunctionInterpolation[fun[x], {x, xmin, xmax}] returns an interpolation of fun[x]. The interpolation is a polynomial in Exp[x], making this better suited to functions which have a logarithmic behaviour. For log(x), the interpolation becomes exact.

Options: the same as FunctionInterpolation
"

FindNumLimit::usage = "FindNumLimit[f, {x, xmin, xmax}]: when f or f, but not both, are non-numeric (technically NumericQ is false), FindNumLimit finds the xc where the transition from numeric to non-numeric happens. The returned point is always such that f is numeric there.

Options and their defaults:
PrecisionGoal -> 5, the number of digits of precision wanted
MaxRecursion -> 20, the maximum number of recursive subdivisions
"

FindZeroCrossings::usage = "FindZeroCrossings[data_List, valueFun :  #[[2]]&] looks
for consecutive elements of the list data s.t. valueFun[element] changes sign between them."

FindFunctionRoots::usage = "Finds all the zeroes of a function, up to a finite search accuracy set by NumPoints"

FindConsecutiveRegions::usage = "Given a list of extrema for a function, splits an area in the I quadrant of R^2 into regions that are consecutive in y. This is what you'd need for example in doing integrals in the y direction for a function which is not an injection."

Densify::usage = "Densify[list, n] takes a list of numbers and returns a list which has n times more points, linearly interpolated between the original points."

FmTune::usage = "Plays a selected tune. Useful for signaling completion of long calculations etc"

IndyTune::usage = "A tune for FmTune, suitable for archeological discoveries and other pompous events."

BassLine::usage = "A bassline, suitable for mounting suspension"

Begin["`Private`"]

RemoveImaginaries[obj_, opts : OptionsPattern[{Tolerance->10^(-9)}] (*Options and their defaults*)] := Module[{maxRelativeIm = 0, maxAbsoluteIm = 0, maxAbsoluteImRelative = 0, ret},
Quiet[
ret = obj /. Complex[x_, y_] :> (If[Abs[y] > OptionValue[Tolerance] Abs[x], Message[RemoveImaginaries::LargeIm, x, y]; 
							If[Abs[y] > maxAbsoluteIm, (maxAbsoluteIm = Abs[y]; maxAbsoluteImRelative = Abs[y/x])]; If [Abs[y/x] > maxRelativeIm, maxRelativeIm = Abs[y/x]]];
x (*In the end, just return x*)),{Power::infy}];
If[maxAbsoluteIm > 0, Message[RemoveImaginaries::LargestIms, maxAbsoluteIm, maxRelativeIm, maxAbsoluteImRelative]];
ret
]

ExtractInterpolatingFunctionScalings[_InterpolatingFunction[Plus[ofs_,Times[mult_,_Symbol]]]] := {1, mult, ofs};

ExtractInterpolatingFunctionScalings[Times[InterpolatingFunction_[Plus[ofs_,Times[mult_,_Symbol]]], scale1_]]:= {scale1, mult, ofs}


(*InterpolatingFunctionToList: the first form operates on an actual pure InterpolatingFunction, and the two other forms extract the scaling parameters and apply them a step at a time, when necessary.*)
InterpolatingFunctionToList[fun_InterpolatingFunction] := Transpose[{Flatten[InterpolatingFunctionGrid[fun]], InterpolatingFunctionValuesOnGrid[fun]}];
InterpolatingFunctionToList[fun_InterpolatingFunction[Plus[ofs_,Times[mult_,_Symbol]]]]:=  (({(#[[1]] - ofs)/mult, #[[2]] }&) /@ InterpolatingFunctionToList[fun]);
InterpolatingFunctionToList[Times[fun_InterpolatingFunction[Plus[ofs_,Times[mult_,x_Symbol]]], scale1_]]:=  ({#[[1]], #[[2]]* scale1}&)/@InterpolatingFunctionToList[fun[ofs + mult x]];
(*The form below is needed to make this also work when the argument in fun[x] is also in the function call. In addition, the code authors limited skills are not sufficient
 to figure out a way to extract the InterpolatingFunction itself and call the first definition above, so unfortunately code repetition could not be avoided here.
 
 Any suggestions to make this more elegant are welcome.
*)
InterpolatingFunctionToList[fun_InterpolatingFunction[_Symbol]] :=  Transpose[{Flatten[InterpolatingFunctionGrid[fun]], InterpolatingFunctionValuesOnGrid[fun]}];


SetAttributes[NIntegralFunction, HoldAll]
SyntaxInformation[NIntegralFunction]  = {"LocalVariables" -> {"Plot", {2, 2}}}
Options[NIntegralFunction] = Join[Options[NDSolve], Options[Integrate],
	{IntegrationGrid -> None,(*This option sets Switches to a different method of integration, where an 
							InterpolatingFunction is formed from the expression, which is then integrated by built-in methods.
							The grid is given by this option either directly as a list of points, or indirectly as an InterpolatingFunction
							which's grid is extracted.
							*)
	EnforceEndpointOrder -> True, (*If this is true, an affine map is used to make sure that Integrate does the numerical sum in
									the interpolating function in the order determined by the endpoints. This may be significant if there are
									several orders of magnitude changes in the function along the interval*)
	Densification -> 8 (*Controls whether the grid in the previous option is made denser by adding points*)}];
NIntegralFunction[expr_, intv_, opts : OptionsPattern[]] := Module[{expr2, x0, x1, var = ReleaseHold[Hold[intv] /. {x_, y__} :> HoldPattern[x]], locvar, y},
expr2 = ReleaseHold[Hold[expr] /. var :> locvar];
{x0, x1} = Rest[intv];
If[OptionValue[IntegrationGrid] === None,
   (*The default choice, use NDSolve*)
	Module[{sols},
		sols = NDSolve[{y'[x] == expr2 /. locvar -> x, y[x0] == 0}, {y},{x, x0, x1}, Evaluate[FilterRules[{opts}, Options[NDSolve]]]];
		First[y/.sols]
	],
	(*In this case, we're given a list of points to evaluate at*)
	If[Head[OptionValue[IntegrationGrid]] === List,
		Module[{gridpoints = Densify[OptionValue[IntegrationGrid],OptionValue[Densification]], ifun, affine, invaffine},
			If[OptionValue[EnforceEndpointOrder] && x0 > x1,
				affine[x_] := x0 - x,
				affine[x_] := x,
				affine[x_] := x (*Just in case*)
			];
			PrintV[StringForm["Integrating `2` on grid `1`", gridpoints, expr2], "Debug"];
			ifun[t_] = N[Integrate[Interpolation[{affine[#], affine'[#](expr2 /. locvar -> #)}& /@ gridpoints][t], t, Sequence @@ Evaluate[FilterRules[{opts}, Options[Integrate]]]],
						OptionValue[WorkingPrecision]];
			Function[x, ifun[affine[x]] - ifun[affine[x0]]] (*affine is self-inverse*)
		],
		(*In this case, the grid is an interpolatingfunction. The complication is though, that there's no
		guarantee that the argument is our integration variable. Do some pattern matching to attempt to divine that*)
		Module[{integrate},
			(*Option1: the grid is just an InterpolatingFunction. In this case, let's just hope that the argument is indeed our integration variable.*)
			integrate[gridfun_InterpolatingFunction] := Module[{grid = Flatten[InterpolatingFunctionGrid[gridfun]]},
														If[Head[grid] === List,
															NIntegralFunction[expr2 /. locvar -> x, {x, x0, x1}, IntegrationGrid -> grid, opts],
															Undefined (*Return undefined if we couldn't extract the grid*)
														]];
			(*Option 2: the grid is an InterpolatingFunction which has a more complicated argument than just the integration variable.*)
			integrate[gridfun : InterpolatingFunction[__][inexpr_]] := Module[{fun, invfun, xpoints, ifun},
				fun = ReleaseHold[Hold[inexpr] /. var :> locvar];
				PrintV[StringForm["Function is `1`\.1d", fun], "Debug"];
				xpoints = Flatten[InterpolatingFunctionGrid[Head[gridfun]]];
				(*In order to support non-monotonic functions, we take the Union of the real intervals produced by all branches of the solution.*)
				invfun[x1_] = locvar /. Solve[fun == x1, locvar];

				PrintV[StringForm["Inverse function is `1`", invfun[y]], "Debug"];

				ifun = NIntegralFunction[D[First[invfun[x]], x] (expr2 /. locvar -> First[invfun[x]]), {x, fun /. locvar -> x0, fun /. locvar -> x1}, IntegrationGrid -> xpoints, opts];
				(*Function[x, ifun[x]]*)
				Function[x, Evaluate[ifun[fun /. locvar -> x]]]
			];
			integrate[OptionValue[IntegrationGrid]]
		]
	]
]

]


(*The following are helpers needed by IntegrateWithInterpolatingFunctions*)
Options[PolyListFromInterpolatingFunction] = {PolyOrder -> Automatic};
PolyListFromInterpolatingFunction[intp_, var_, opts : OptionsPattern[]] := Module[{points, order, intvs, x, polys},
points = InterpolatingFunctionToList[intp];
order = If[OptionValue[PolyOrder] === Automatic, First[InterpolatingFunctionInterpolationOrder[intp]], OptionValue[PolyOrder]];
Print[StringForm["Creating a polynomial of order `1`, PolyOrder = `2`", order, OptionValue[PolyOrder]]];
intvs = Partition[points, order + 1, 1];
polys = {{First[#[[Ceiling[order/2]]]], First[#[[Ceiling[order/2] + 1]]]}, InterpolatingPolynomial[#, var]}&/@intvs;
(*We're almost finished, but the poly intervals don't contain the ends now. Extend them and return the result.*)
MapAt[{{points[[1, 1]], #[[1, 2]]}, #[[2]]}&, MapAt[{{#[[1, 1]], points[[-1, 1]]}, #[[2]]}&, polys, -1], 1]
]
PolyListToCoefList[polys_, var_] := Module[{},
{#[[1]], CoefficientList[#[[2]], var]}&/@ polys
]

Options[CoefListFromInterpolatingFunction] = Options[PolyListFromInterpolatingFunction];
CoefListFromInterpolatingFunction[intp_, var_, opts : OptionsPattern[]] := PolyListToCoefList[PolyListFromInterpolatingFunction[intp, var, opts], var];

Options[InterpolatingFunctionToGenericPolynomial] = {PolyOrder -> Automatic}
InterpolatingFunctionToGenericPolynomial[intp_, coefname_, opts : OptionsPattern[]] := Module[{order, poly, coef},
order = If[OptionValue[PolyOrder] === Automatic, First[InterpolatingFunctionInterpolationOrder[Head[intp]]], OptionValue[PolyOrder]];
poly[x_] :=  Sum[coefname[n] x^n, {n, 0, order}];
poly @@ intp
]

SetAttributes[IntegrateWithInterpolatingFunctions, HoldAll];
Options[IntegrateWithInterpolatingFunctions] = Join[{CompileResult -> True, WorkingPrecision -> $MachinePrecision}, Options[PolyListFromInterpolatingFunction]];
IntegrateWithInterpolatingFunctions::outofrange = "Warning: Argument `1` is out of range, using extrapolation";
SyntaxInformation[IntegrateWithInterpolatingFunctions] = {"LocalVariables" -> {"Plot", {2, 2}}, "ArgumentsPattern" -> {_, _, OptionsPattern[]}};
IntegrateWithInterpolatingFunctions[expr_, intv_, opts: OptionsPattern[]] := Module[{polyexpr, intgexpr, intps = {}, coefnum = 1, coefs, xpoints,xpoints1,  xintvs, localexprs, combinedrules, accu, fun1, fun, var = ReleaseHold[Hold[intv] /. {x_, y__} :> HoldPattern[x]], expr2, xmin,  xmax, localvar, x0, x1, wp = OptionValue[WorkingPrecision]}, 

Block[{$vcontext = "IntegrateWithInterpolatingFunctions"},

expr2 = ReleaseHold[Hold[expr] /. var :> localvar];

{x0, x1} = Rest[intv];
{xmin, xmax} = Sort[{x0, x1}]; (*Sort the endpoints.*)

(*First we go through the expression. Anything which's Head's Head is an InterpolatingFunction (the form is InterpolatingFunction[params][var], which is why we take the double head) is processed by appending it to a list of interpolating functions, and extracting the piecewise coefficients of the polynomials.*)
polyexpr =  (If[Head[Head[#]] === InterpolatingFunction, 
(intps = Append[intps, {Head[#], CoefListFromInterpolatingFunction[Head[#], localvar,  Evaluate[Sequence @@ FilterRules[{opts}, Options[CoefListFromInterpolatingFunction]]]],#[[1]](*This is the variable involved*)}];
InterpolatingFunctionToGenericPolynomial[#, coefs[coefnum++], Evaluate[Sequence @@ FilterRules[{opts}, Options[InterpolatingFunctionToGenericPolynomial]]]]) , #])&//@ expr2;

(*Take the list of points, which will determine the final intervals*)
(*If the argument of an interpolating function is more complicated, we need to map
the interpolation points with the inverse function of that, to find the points
where we need to change the interpolating functions.*)
xpoints1 = Module[{invfun, xpoints2, infun},
infun = Last[#];
xpoints2 = Flatten[InterpolatingFunctionGrid[First@#]];
(*In order to support non-monotonic functions, we take the Union of the real intervals produced by all branches of the solution.*)
invfun[x1_] = localvar /. Solve[infun == x1, localvar];

PrintV[StringForm["Inverse function is `1`", invfun[y]], "Debug"];

PrintV[StringForm["Local points list `1`", Union[Select[Flatten[invfun/@ xpoints2], Function[point, Element[point, Reals]]]]], "Debug"];

Union[Select[Flatten[N[invfun/@ xpoints2, wp]], Function[point, Element[point, Reals]]]]

]&/@ intps;
(*Then take the union of all these points, add the initial and final points, and only take the points in the interval.*)
xpoints = Select[Union[Flatten[Join[xpoints1, {xmin, xmax}]]], xmin <= # <= xmax&];
xintvs = Partition[xpoints, 2, 1];

(*Turn the coefficients into a list of replacement rules*)
intps = Function[intpn, {#[[1]], Function[coefnum, coefs[intpn][coefnum - 1] -> #[[2, coefnum]]]/@ Range[1, Length[#[[2]]]], intps[[intpn, 3]]}&/@ intps[[intpn, 2]]] /@ Range[1, Length[intps]];

(*Combine the intervals and replacement rules*)
combinedrules = Function[interval, {interval, Flatten[Function[rule, Select[rule, (#[[1, 1]] <=  (#[[3]] /. localvar ->((interval[[1]] + interval[[2]])/ 2))<= #[[1, 2]])&, 1][[1, 2]]] /@ intps]}]/@ xintvs;

intgexpr = HornerForm[Integrate[polyexpr, localvar]]; (*The horner form is faster and more accurate to evaluate. 
													The problem is that the result is not necessarily a polynomial. What to do?*)

PrintV[StringForm["Polynomial expression is `1`, which integrates to `2`", polyexpr, intgexpr], "Debug"];

accu = N[0, wp];
 localexprs = Function[nintv, 
Module[{localexpr, localexprshifted, first, last},
 PrintV[StringForm["evaluating interval `1`", combinedrules[[nintv, 1]]], "Debug"];
localexpr = N[intgexpr /. combinedrules[[nintv, 2]], wp];
(*The order of evaluating the endpoints depends on the order the integration endpoints. We deal with that here.*)
 {first, last} = If[x0 < x1, {1, 2}, {2, 1}];
localexprshifted = N[localexpr - (localexpr /. localvar -> combinedrules[[nintv, 1, first]]) + accu, wp];
accu = N[localexprshifted /. localvar ->combinedrules[[nintv, 1, last]], wp]; 
{combinedrules[[nintv, 1]],  localexprshifted}
]
]/@ (If[x0 < x1, #, Reverse[#]]& @ Range[1, Length[xintvs]]);

(*Construct the result to a function. For compatibility with InterpolatingFunction, issue a warning but extrapolate the result if out of range*)
fun1[x_?NumericQ] := Module[{express = Select[localexprs, #[[1, 1]] <=  x <=  #[[1, 2]]&, 1]},
If[express === {},
Message[IntegrateWithInterpolatingFunctions::outofrange, x];
express = Sort[{First[localexprs], Last[localexprs]}, Function[{expre1, expre2}, Abs[expre1[[1, 1]] - x] < Abs[expre2[[1, 2]] -x]]];
];
N[express[[1, 2]] /. localvar-> x, wp]]
;

(*Attempt to compile the function, unless forbidden by the user.*)
If[OptionValue[CompileResult] === True,
fun = Compile[{x}, fun1[x]],
fun = fun1
];

Derivative[n_][fun][x_] := D[Select[localexprs, #[[1, 1]] <=  x <=  #[[1, 2]]&, 1][[1, 2]], {localvar, n}] /. localvar-> x;
Coefficients[fun][x_] := Select[localexprs, #[[1, 1]] <=  x <=  #[[1, 2]]&, 1][[1, 2]] /. localvar -> "x";

fun
]
]


Options[ExpFunctionInterpolation] := Options[FunctionInterpolation];
SetAttributes[ExpFunctionInterpolation, HoldAll]; 
ExpFunctionInterpolation[fun_, {var_Symbol, min_?NumericQ, max_?NumericQ}, opts : OptionsPattern[]] := Module[{intlog, int},
intlog = FunctionInterpolation[Evaluate[fun /. var -> Exp[logvar]], {logvar, Log[min], Log[max]}, Evaluate[Sequence @@ FilterRules[{opts}, Options[FunctionInterpolation]]]];
intlog[Log[#]]&
]


FindNumLimit::PrecisionNotReached = "The requested precision could not be reached within MaxRecursion subdivisions, estimated error is `1`";
FindNumLimit::notbracketed = "The given endpoints do not bracket a transition from numeric to non-numeric, left value = `1`, right value = `2`."
Options[FindNumLimit] = {PrecisionGoal -> 5, MaxRecursion -> 20};

SetAttributes[FindNumLimit, HoldAll]
SyntaxInformation[FindNumLimit]  = {"LocalVariables" -> {"Plot", {2, 2}}}

FindNumLimit[expr_, intv_, opts: OptionsPattern[]] := Module[{fun, x0, x1, subDivide, interval, xminval, xmaxval, var = ReleaseHold[Hold[intv] /. {x_, y__} :> HoldPattern[x]],  locvar, expr2},
Block[{$vcontext = "FindNumLimit"},

 expr2 = ReleaseHold[Hold[expr] /. var :> locvar];
{x0, x1} = Rest[intv];

 fun[x_] := expr2 /. locvar -> x;

xminval = Quiet[fun[x0]];
xmaxval = Quiet[fun[x1]];

If[NumericQ[xminval] == NumericQ[xmaxval],
 Message[FindNumLimit::notbracketed, xminval, xmaxval];
 Return[Undefined];
];

subDivide[{left_, right_, leftval_, rightval_}]:= Module[{midpoint, midval},
midpoint= (left+right)/2;  
midval = Quiet[fun[midpoint]]; (*No need to print out the error messages, since failure is a valid result here*)
If[NumericQ[midval] == NumericQ[leftval],
 (
  PrintV[StringForm["Going right, fun[`1`] = `2`, fun[`3`] = `4`, fun[`5`] = `6`", left, leftval, midpoint, midval, right, rightval],"Debug"];
 {midpoint, right, midval, rightval}
 ),
 (
 PrintV[StringForm["Going left, fun[`1`] = `2`, fun[`3`] = `4`, fun[`5`] = `6`", left, leftval, midpoint, midval, right, rightval],"Debug"];
 {left, midpoint, leftval, midval}
 )
]];

interval = NestWhile[subDivide[#]&, {x0, x1, xminval, xmaxval}, Abs[#[[1]]-#[[2]]] > 10^(-OptionValue[ PrecisionGoal])&, 1, OptionValue[MaxRecursion]];
	(*Best guess at the midpoint*)

(*Error control*)
If[Abs[interval[[1]] - interval[[2]]]  > 10^(-OptionValue[PrecisionGoal]), Message[FindNumLimit::PrecisionNotReached, Abs[interval[[1]] - interval[[2]]]]];

(*Return the numeric side of the interval*)
If[NumericQ[interval[[3]]], interval[[1]], interval[[2]]]
]
]


FindZeroCrossings[data_List, valueFun_ : (#[[2]]&)] := (First /@ #)&/@ Select[Partition[data, 2, 1], (Sign[valueFun[#[[1]]]] != Sign[valueFun[#[[2]]]]) && (valueFun[#[[2]]] \[Element] Reals) && (valueFun[#[[1]]] \[Element] Reals)&]


SetAttributes[FindFunctionRoots, HoldAll]
Options[FindFunctionRoots] = Join[{NumPoints -> 300}, Options[FindRoot]]
SyntaxInformation[FindFunctionRoots]  = {"LocalVariables" -> {"Plot", {2, 2}}}
FindFunctionRoots[expr_, intv_, opts: OptionsPattern[]] := Module[{x0, x1, table, crossings, x,  var = ReleaseHold[Hold[intv] /. {x_, y__} :> HoldPattern[x]],  locvar, expr2},

expr2 = ReleaseHold[Hold[expr] /. var :> locvar];
{x0, x1} = Rest[intv];

table = Table[{i, expr2 /. locvar -> i}, {i, x0, x1, (x1 - x0)/OptionValue[NumPoints]}];

crossings = FindZeroCrossings[table];

(If[(#[[1]] != 0) && (#[[2]] != 0), x/. FindRoot[N[expr2 /. locvar -> x], {x, #[[1]], #[[2]]}, Method ->  "Brent", Evaluate[Sequence @@ FilterRules[{opts}, Options[FindRoot]]]], First[Sort[#]]]&) /@ crossings
]


FindConsecutiveRegions[xmin_?NumericQ, miny_?NumericQ, maxy_?NumericQ, extrema_List] := Module[{maxs, mins, zeros, pos, firstmin, midy},
(*First always determine the mid lambda. A valid choice is any local maximum.*)
maxs = Select[extrema, (#[[3]] == True)&];
midy = If[maxs === {}, Return[{xmin, Infinity, miny, miny, maxy}], First[maxs][[2]]];

(*Get the minimum with the smallest x*)
pos=Ordering[Select[extrema, (#[[3]] == False) &][[All,1]]];
If[pos === {}, Return[{xmin, maxs[[1, 1]], miny,midy, maxy}]];
     firstmin = extrema[[First[pos]]];

(*Divide the further boxes at the local minimum.*)
{FindConsecutiveRegions[firstmin[[1]], firstmin[[2]], maxy, Select[extrema, (#[[2]] > firstmin[[2]])&]], FindConsecutiveRegions[firstmin[[1]], miny, firstmin[[2]], Select[extrema, (#[[2]] < firstmin[[2]])&]], {xmin, firstmin[[1]], miny, midy, maxy}}
]


Densify[list_List, n_?NumericQ] := Flatten[{Table[((n - i)#[[1]] + i #[[2]])/n, {i, 0, n - 1}]& /@ Partition[list, 2, 1], Last[list]}]


NoteFreq[semisfromA440_]:= 440. 2^(semisfromA440/12);
FmNote[carrier_, modf_, len_, index_]:=Play[ Exp[-3t ](1-HeavisideTheta[t/len-0.9](t/len - 0.9)10) Sin[2 Pi carrier  t +  index  Exp[-2t] Sin[2 Pi modf carrier  t]], {t, 0, len}, PlayRange-> {-1, 1}, SampleRate-> 44100, SampleDepth-> 24]
FmTune[Notes_List, mod_, transpose_ : 0, bpm_ : 50]:= Sound[FmNote[NoteFreq[#[[1]] + transpose], mod, #[[2]]60/bpm, #[[3]]]]& /@ Notes;


BassLine = {{-24, 0.25, 1}, {-24, 0.25, 2}, {-12, 0.25, 2}, {-24, 0.25, 1},{-24, 0.25, 1}, {-17, 0.25, 1.4}, {-24, 0.25, 1.3}, {-24, 0.25, 1}, {-16, 0.25, 1.2}, {-24, 0.25, 1}, {-24, 0.25, 1},  {-19, 0.25, 1.4}, {-24, 0.25, 1}, {-24, 0.25, 1.9}, {-21, 0.25, 2.6}, {-19, 0.25, 3.6}, {-21, 0.25, 1.8}, {-26, 0.25, 1.7}, {-24, 2, 4.6}};


IndyTune = {{2, 0.75, 2.5}, {3, 0.25, 2.5}, {5, 0.5, 2.5}, {10, 2, 4.5},{0, 0.75, 2.5}, {2, 0.25, 2.5}, {3, 2, 3.5}, {5, 0.75, 2.5}, {7, 0.25, 2.5}, {9, 0.5, 3.5}, {15, 2, 1.5}, {7, 0.75, 2}, {9, 0.25, 3.5}, {10, 1, 4}, {12, 1, 4.5}, {14, 1.5, 5}};


End[ ]

EndPackage[ ]
