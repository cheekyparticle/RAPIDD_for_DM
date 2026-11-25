(* Export the polynomial coefficients of the nuclear responses of DMFormFactor v6 (arXiv:1308.6288).

   wolframscript -file export_responses.wls <dmformfactor-V6.m> <Elastic folder> <model> <output> [nuclei]

   model = Anand2013    : default one-body density matrices (OBDMs) of DMFormFactor, SetIsotope[Z, A, "default", "default"]
   model = McElvain2022 : OBDMs of the Elastic repository (github.com/Berkeley-Electroweak-Physics/Elastic,
                          files <Elastic folder>/<name>.txt), injected into DMFormFactor after SetIsotope
                          (which sets b, J and T = |N - Z|/2).
   nuclei (optional)    : comma-separated A values (e.g. 19,23) to export only some nuclei of the model.

   Each output line is
     A Z response t1 t2 J a0 ... a10
   with response = M, SigmaPP, SigmaP, PhiPP, PhiTildeP, Delta, M_PhiPP, SigmaP_Delta (ResponseNuclear j = 1..8),
   t1, t2 = 0, 1 (isospin +, -) and a_k the coefficients of y^k of ResponseNuclear[y, j, t1, t2] E^(2 y) (exact
   polynomials of DMFormFactor, not fits). Machine numbers are written in InputForm (x*^n for x 10^n). *)

args = $ScriptCommandLine;
If[Length[args] < 5, Print["usage: export_responses.wls <dmformfactor-V6.m> <Elastic folder> <model> <output> [A,A,...]"]; Exit[1]];
{pkg, elastic, model, out} = args[[2 ;; 5]];
only = If[Length[args] >= 6, ToExpression /@ StringSplit[args[[6]], ","], All];

Get[pkg];

(* {Z, A, Elastic file} *)
nuclei = Switch[model,
  "Anand2013", {{9, 19, ""}, {11, 23, ""}, {14, 28, ""}, {14, 29, ""}, {14, 30, ""},
    {32, 70, ""}, {32, 72, ""}, {32, 73, ""}, {32, 74, ""}, {32, 76, ""}, {53, 127, ""},
    {54, 128, ""}, {54, 129, ""}, {54, 130, ""}, {54, 131, ""}, {54, 132, ""}, {54, 134, ""}, {54, 136, ""}},
  "McElvain2022", {{6, 12, "C12_ck"}, {6, 13, "C13_ck"}, {8, 16, "O16_4hw"}, {8, 18, "O18_usdb"},
    {9, 19, "F19_usdb"}, {11, 23, "Na23_usdb"}, {13, 27, "Al27_usdb"},
    {14, 28, "Si28_usdb"}, {14, 29, "Si29_usdb"}, {14, 30, "Si30_usdb"}, {18, 40, "Ar40_Sussex"},
    {20, 40, "Ca40_kbp"}, {32, 70, "Ge70_GCN2850"}, {32, 72, "Ge72_GCN2850"}, {32, 73, "Ge73_GCN2850"},
    {32, 74, "Ge74_GCN2850"}, {32, 76, "Ge76_GCN2850"}, {53, 127, "I127_GCN5082"}, {55, 133, "Cs133_GCN5082"},
    {54, 128, "Xe128_GCN5082"}, {54, 129, "Xe129_GCN5082"}, {54, 130, "Xe130_GCN5082"},
    {54, 131, "Xe131_GCN5082"}, {54, 132, "Xe132_GCN5082"}, {54, 134, "Xe134_GCN5082"}, {54, 136, "Xe136_GCN5082"}},
  _, Print["unknown model ", model]; Exit[1]];
If[only =!= All, nuclei = Select[nuclei, MemberQ[only, #[[2]]] &]];

(* Elastic files: blocks with the header 2Ji 2Ti 2Jf 2Tf 2Jop 2Top, then rows N 2j N' 2j' value, closed by -1 / -2.
   Same convention as the internal tables of DMFormFactor {value, Jop, Top, N, 2j, N', 2j'}. *)
loadElastic[z_, a_, file_] := Module[{lines, J, T, rows = {}, Ji = None, Ti = None},
  lines = Select[StringSplit /@ ReadList[file, String], Length[#] > 0 &];
  Do[Which[Length[l] == 6, J = ToExpression[l[[5]]]/2; T = ToExpression[l[[6]]]/2;
             If[Ji === None, Ji = ToExpression[l[[1]]]/2; Ti = ToExpression[l[[2]]]/2],
           Length[l] == 5, AppendTo[rows, Join[{ToExpression[l[[5]]], J, T}, ToExpression /@ l[[1 ;; 4]]]]], {l, lines}];
  SetIsotope[z, a, "default", "default"];
  If[Ji =!= dmformfactor`Private`JIso,
    Print["ERROR: J of ", file, " = ", Ji, ", J of DMFormFactor = ", dmformfactor`Private`JIso]; Exit[2]];
  If[Ti =!= Abs[a - 2 z]/2,
    Print["ERROR: T of ", file, " = ", Ti, ", |N - Z|/2 = ", Abs[a - 2 z]/2]; Exit[2]];
  Clear[dmformfactor`Private`myDM];
  Do[dmformfactor`Private`myDM[i] = rows[[i]], {i, Length[rows]}];
  dmformfactor`Private`DensityMatrix = dmformfactor`Private`myDM;
  SetDMmatrixLength[Length[rows]];];

names = {"M", "SigmaPP", "SigmaP", "PhiPP", "PhiTildeP", "Delta", "M_PhiPP", "SigmaP_Delta"};
str[x_] := If[x === 0, "0", ToString[N[x], InputForm]];
stream = OpenWrite[out];
Do[{z, a, f} = n;
  t = First@AbsoluteTiming[
    If[model == "Anand2013", SetIsotope[z, a, "default", "default"],
       loadElastic[z, a, FileNameJoin[{elastic, f <> ".txt"}]]];
    jn = dmformfactor`Private`JIso;
    Do[c = CoefficientList[Expand[ResponseNuclear[yy, j, t1, t2] E^(2 yy)], yy];
       If[Length[c] > 11, Print["ERROR: degree > 10 for ", a, " ", names[[j]]]; Exit[3]];
       If[! FreeQ[c, yy], Print["ERROR: not a polynomial for ", a, " ", names[[j]]]; Exit[3]];
       WriteString[stream, StringRiffle[Join[{ToString[a], ToString[z], names[[j]], ToString[t1], ToString[t2],
           ToString[N[jn], InputForm]}, str /@ PadRight[c, 11]], " "], "\n"],
       {j, 1, 8}, {t1, 0, 1}, {t2, 0, 1}]];
  m0 = ResponseNuclear[0, 1, 0, 0] 16 Pi/(2 jn + 1);
  Print[a, " (Z = ", z, ", ", If[f == "", "default", f], "): J = ", jn, ", T = ", dmformfactor`Private`TIso,
    ", time ", Round[t, 0.1], " s, M++(0) 16pi/(2J+1) = ", ToString[N[m0], InputForm], " (A^2 = ", a^2, ")"],
  {n, nuclei}];
Close[stream];
