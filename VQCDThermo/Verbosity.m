(* ::Package:: *)

BeginPackage["VQCDThermo`Verbosity`"]

(*Set up defaults for global parameters that determine how much debug/intermediate output is generated.
 Currently the allowed values are:
	"All": the maximal amount of intermediate results is printed and plotted.
  "Debug": many things, including checks that must hold always if the code is correct, are printed.
"Checks": values and functions where the calculation could silently fail are printed/plotted such as fit results against the fitted points etc.
"None": Nothing is printed or plotted.
*)

PrintV::usage = "PrintV[output, level] conditionally prints output, depending on level and the current setting of SetVerbosity. If level <= GetVerbosity[], PrintV[output, level] is equivalent to Print[output], otherwise it does nothing.

The first parameter is only evaluated after checking the level, so lengthy Plot[...]'s etc will not take time if the current verbosity level is below the level determined in PrintV."
EvaluateV::usage = "EvaluateV[level, expr1, expr2] conditionally evaluates expr1 or expr2, depending on level and the current setting of SetVerbosity. If level <= GetVerbosity[], expr1 is evaluated, otherwise expr2."

(*The next two functions may be at a slightly wrong place, or perhaps this package could be renamed DebugTools etc.*)
AnalyzeConditional::usage = "Subdivides conditional expressions at operators And, Or and Xor, and prints out truth values of the subexpressions. Returns the value of the full expression";
AnalyzeConditionalV::usage = "Evaluates the given conditional, or applies AnalyzeConditional, depending on verbosity"

SetVerbContext::usage = "Sets the context for the verbosity, to allow varying the verbosity level locally."
ResetVerbContext::usage = "Resets the verb context to the previous value."

SetVerbosity::usage = "SetVerbosity[level] sets the maximum level of verbosity that is printed out.

The allowed values of level are:
\"All\": the maximal amount of intermediate results is printed and plotted.
\"Debug\": many things, including checks that must hold always if the code is correct, are printed.
\"Checks\": values and functions where the calculation could silently fail are printed/plotted such as fit results against the fitted points etc.
\"None\": Nothing is printed or plotted.
"

$vcontext::usage = "Sets the name of the current verbosity context";

GetVerbosity::usage  = "GetVerbosity[] returns the current verbosity level."

Begin["`Private`"]

AllowedVerbs = "All"|"Debug"| "Checks"| "Progress" |  "None";
Verbtest = (StringMatchQ[#, AllowedVerbs]&);
VerbLevels["All"] = 40;
VerbLevels["Debug"] = 30;
VerbLevels["Checks"] = 20;
VerbLevels["Progress"] = 10;
VerbLevels["None"] = 0;
SetVerbosity[lev_Integer]:= (Verbosity = lev)
SetVerbosity[lev_?Verbtest] := (Verbosity = VerbLevels[lev]);
SetVerbosity[lev_Integer, context_] := (vcontextlevels[ToString[context]] = lev);
SetVerbosity[lev_?Verbtest, context_] := (vcontextlevels[ToString[context]] = VerbLevels[lev]; Print[StringForm["Set vlevel for context `1`", ToString[context]]]);
Verbosity = 10;

$vcontext = "Global";
vcontextlevels["Global"] = -1;

GetVContextLevel[context_] := vcontextlevels[ToString[context]];

GetVerbosity[] := Verbosity;

SetAttributes[PrintV, HoldFirst](*This prevents evaluation of the first argument of PrintV. This way calling for example PrintV[Plot[Slowfunction ... ], "Debug"] will not spend time calculating the plot if $Verbosity is below "Debug".*)
SetAttributes[EvaluateV, HoldRest]

(*Evaluates expr1 or expr2 based on verbosity.*)
EvaluateV[level_?Verbtest, expr1_, expr2_ : Null] := If[(VerbLevels[level] >  Verbosity) && !TrueQ[(VerbLevels[level] <= vcontextlevels[$vcontext])],
	expr2,
	expr1,
	expr1
];

PrintV[out_, level_?Verbtest] := EvaluateV[level, Print[out], Null];

SetAttributes[AnalyzeConditional, HoldFirst]
AnalyzeConditional[expr_, topLevel_ : True]:= Module[{
comps = (Hold/@Unevaluated[expr]),
result = expr},
If[topLevel,
	Print[StringForm["Expression `1`\n\t is `2` ", HoldForm[expr], result]];
];
Function[x, AnalyzeConditional[Unevaluated[x], False], HoldFirst] @@@ comps;
result
]/;MatchQ[Unevaluated@expr,_And|_Or|_Xor|_Not];
AnalyzeConditional[expr_, False] := Print[StringForm["Subexpression `1`\n\t is `2`", HoldForm[expr], expr]];

SetAttributes[AnalyzeConditionalV, HoldFirst]
AnalyzeConditionalV[expr_, level_?Verbtest] := EvaluateV[level, AnalyzeConditional[expr], expr];

(*PrintV[out_, level_?Verbtest] := ( 
(*Print[StringForm["Current vcontextlevel = `1`", vcontextlevels[Last[vcontext]]]];*)
(*Print[StringForm["vcontextlevels = `1`, current vcontext = `2`", DownValues[vcontextlevels], $vcontext]];*)
If[(VerbLevels[level] >  Verbosity) && !TrueQ[(VerbLevels[level] <= vcontextlevels[$vcontext])],
Return[]];
(*Print[StringForm["Current evaluation stack `1`",stack]];*)
Print[out];
);*)

End[ ]

EndPackage[ ]












