#!/usr/bin/env wolframscript
(* ::Package:: *)

Print[$ScriptCommandLine]
filebase=$ScriptCommandLine[[2]]
Print[AbsoluteTiming[edges=#[[1]]->#[[2]]&/@Import[filebase<>"edges.csv"][[2;;-1]];]];
Print[AbsoluteTiming[g=Graph[Range[Length[components]],edges,GraphLayout->{"VertexLayout"->"LayeredDigraphEmbedding"},ImageSize->{Automatic,500},AspectRatio->1/2,VertexSize->0.002*Length[components],PerformanceGoal->"Quality",VertexStyle->{Directive[Black]},EdgeStyle->Directive[Opacity[0.1],Arrowheads[0]]];]]
Print[AbsoluteTiming[p=LayeredGraphPlot[g];]]
Print[AbsoluteTiming[vcs=p[[1,1,2,3;;-1]]/.Disk[u_,r_]:>u/.DynamicName[u_, v_]:>u/.Point[u_]:>u;]]
Print[AbsoluteTiming[ecs=p[[1,1,1,4;;-1]]/.{Arrow[{DynamicLocation[u_,Automatic,Center],DynamicLocation[v_,Automatic,Center]}]:>{vcs[[ToExpression[StringSplit[u,"$"][[2]]]]],vcs[[ToExpression[StringSplit[v,"$"][[2]]]]]},Arrow[BezierCurve[{DynamicLocation[u_,Automatic,Center],l__,DynamicLocation[v_,Automatic,Center]}]]:>Join[{vcs[[ToExpression[StringSplit[u,"$"][[2]]]]]},{l},{vcs[[ToExpression[StringSplit[v,"$"][[2]]]]]}]};]]
BinaryWrite[filebase<>"nodes.dat",N[Flatten[vcs]],"Real64"]
BinaryWrite[filebase<>"edges.dat",N[Flatten[ecs]],"Real64"]
BinaryWrite[filebase<>"edgelengths.dat",Length/@ecs,"Integer64"]
