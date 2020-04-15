#!/usr/bin/env wolframscript
(* ::Package:: *)

Print[$ScriptCommandLine]

filebase=$ScriptCommandLine[[2]]
dat=Import[filebase<>"out.dat"];
{dim,ns,runtime} = dat[[1,1;;3]]
multiindices=Partition[BinaryReadList[filebase<>"multiindices.npy","Integer64"][[17;;-1]],Length[dat[[2]]]];
elements=dat[[3]];
spatoms=Partition[BinaryReadList[filebase<>"spatoms.npy","Integer64"][[17;;-1]],Length[elements]];
temp=BinaryReadList[filebase<>"temperatures.npy","Real64"][[17;;-1]];
pres=BinaryReadList[filebase<>"pressures.npy","Real64"][[17;;-1]];
rows=BinaryReadList[filebase<>"rows.npy","Integer64"][[17;;-1]];
columns=BinaryReadList[filebase<>"columns.npy","Integer64"][[17;;-1]];
data=BinaryReadList[filebase<>"data.npy","Real64"][[17;;-1]];
states=Partition[Import[filebase<>"states.npy","Real64"][[17;;-1]],dim];
A=SparseArray[Join[Table[{rows[[i]]+1,columns[[i]]+1}->data[[i]],{i,1,Length[data]}]],{dim,dim}];
thrs=10^-3;
vals=A["NonzeroValues"];
pos=A["NonzeroPositions"];
Do[If[A[[pos[[i,1]],pos[[i,2]]]]<thrs*A[[pos[[i,2]],pos[[i,1]]]],vals[[i]]=0],{i,1,Length[pos]}]
A2=SparseArray[#[[1]]->#[[2]]&/@Transpose[{pos,vals}]];
spforms=Table[DisplayForm[RowBox[Table[If[spatoms[[i,j]]>0,If[spatoms[[i,j]]>1,Subscript[elements[[j]],spatoms[[i,j]]],elements[[j]]],""],{j,1,Length[elements]}]]],{i,1,Length[spatoms]}];

Print[AbsoluteTiming[g=Graph[#[[1]]->#[[2]]&/@A2["NonzeroPositions"]];]]
Print[AbsoluteTiming[components=ConnectedComponents[g];]]

Print[AbsoluteTiming[edges=Flatten[Table[ParallelTable[If[Max[Abs[A2[[components[[i]],components[[j]]]]]]>0&&i!=j,i->j,{}],{j,1,i}],{i,1,Length[components]}]];]]
Print[AbsoluteTiming[g=Graph[Range[Length[components]],edges,GraphLayout->{"VertexLayout"->"LayeredDigraphEmbedding"},ImageSize->{Automatic,500},AspectRatio->1/2,VertexSize->0.002*Length[components],PerformanceGoal->"Quality",VertexStyle->{Directive[Black]},EdgeStyle->Directive[Opacity[0.1],Arrowheads[0]]];]]
Print[AbsoluteTiming[p=Show[g];]]
(*Print[AbsoluteTiming[Export[filebase<>"connected.pdf",p]]]*)

Print[AbsoluteTiming[vcs=p[[1,1,2,3;;-1]]/.Disk[u_,r_]:>u/.DynamicName[u_, v_]:>u/.Point[u_]:>u;]]
Print[AbsoluteTiming[ecs=p[[1,1,1,4;;-1]]/.{Arrow[{DynamicLocation[u_,Automatic,Center],DynamicLocation[v_,Automatic,Center]}]:>{vcs[[ToExpression[StringSplit[u,"$"][[2]]]]],vcs[[ToExpression[StringSplit[v,"$"][[2]]]]]},Arrow[BezierCurve[{DynamicLocation[u_,Automatic,Center],l__,DynamicLocation[v_,Automatic,Center]}]]:>Join[{vcs[[ToExpression[StringSplit[u,"$"][[2]]]]]},{l},{vcs[[ToExpression[StringSplit[v,"$"][[2]]]]]}]};]]
BinaryWrite[filebase<>"nodes.dat",N[Flatten[vcs]],"Real64"]
BinaryWrite[filebase<>"edges.dat",N[Flatten[ecs]],"Real64"]
BinaryWrite[filebase<>"edgelengths.dat",Length/@ecs,"Integer64"]
