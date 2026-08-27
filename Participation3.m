%six bar linkage
%static equilibrium

%define the joints
A = [7 4 0];
B = [5 16 0];
E = [18 35 0];
F = [43 32 0];
C = [25 25 0];
D = [23 10 0];
G = [45 17 0];

%center of mass of each link
S1 = (A+B)/2;
S2 = (B+C+E)/3; %not accurate!!!
S3 = (C+D)/2;
S4 = (E+F)/2;
S5 = (F+G)/2;


%define the lengths of bars
lAB= norm(B-A);
lBC= norm(C-B);
lCD= norm(D-C);
lBE = norm(E-B);
lEF = norm(F-E);
lFG = norm(G-F);

WAB = [0 -1 0];

WBEC = [0 -1 0];

WCD = [0 -1 0];

WEF = [0 -1 0];

WFG = [0 -1 0];

syms FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin

ForceA =  [FAx FAy 0];
ForceB = [FBx FBy 0];
ForceC = [FCx FCy 0];
ForceE =  [FEx FEy 0];
ForceF = [FFx FFy 0];
ForceG =  [FGx FGy 0];
ForceD = [FDx FDy 0];
InputTorque = [0 0 Tin];

AppliedForce = [50 0 0];

%Static Equilibrium Conditions for Link AB
%ForceAB
eqn1 = ForceA + ForceB + WAB == 0;

%momentAB
eqn2 = cross(A-S1,ForceA) + cross(B-S1,ForceB) + InputTorque == 0;

%SEC's BEC
%ForceBEC
eqn3 = -ForceB + ForceC + ForceE + WBEC == 0;

%momentBEC
eqn4 = cross(B-S2,-ForceB) + cross(C-S2,ForceC)+ cross(E-S2,ForceE) == 0;

%Static Equilibrium Conditions for Link CD
%ForceCD
eqn5 = -ForceC + ForceD + WCD == 0;

%momentCD
eqn6 = cross(C-S3,-ForceC) + cross(D-S3,ForceD) == 0;

%Static Equilibrium Conditions for Link EF
%ForceEF
eqn7 = -ForceE + ForceF + WEF == 0;

%momentEF
eqn8 = cross(E-S4,-ForceE) + cross(F-S4,ForceF) == 0;

%Static Equilibrium Conditions for Link FG
%ForceFG
eqn9 = -ForceF + ForceG + WFG + AppliedForce == 0;

%momentFG
eqn10 = cross(F-S5,-ForceF) + cross(G-S5,ForceG) == 0;



eqnMatrix= [eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10];

StaticSolution=solve(eqnMatrix,[FAx, FAy, FBx, FBy, FCx, FCy, FDx, FDy, FEx, FEy, FFx, FFy, FGx, FGy, Tin]);
ForceAx = double(StaticSolution.FAx) %[output:09062f1d]
ForceAy = double(StaticSolution.FAy) %[output:1b72e175]
ForceBx = double(StaticSolution.FBx) %[output:52b7ed10]
ForceBy = double(StaticSolution.FBy) %[output:3dc57759]
ForceCx = double(StaticSolution.FCx) %[output:02368dd7]
ForceCy = double(StaticSolution.FCy) %[output:1b178119]
ForceDx = double(StaticSolution.FDx) %[output:61570b74]
ForceDy = double(StaticSolution.FDy) %[output:343aa1e5]
ForceEx = double(StaticSolution.FEx) %[output:4747de2e]
ForceEy = double(StaticSolution.FEy) %[output:1357f864]
ForceFx = double(StaticSolution.FFx) %[output:0cfc0298]
ForceFy = double(StaticSolution.FFy) %[output:661aad6f]
ForceGx = double(StaticSolution.FGx) %[output:0105cfe7]
ForceGy = double(StaticSolution.FGy) %[output:7f8e97c6]
Torque = double(StaticSolution.Tin) %[output:114f8990]



%Angular Velocity Calculations
%Loop ABCDA

syms wBEC wCD
omega_AB=[0 0 1];
omega_BEC=[0 0 wBEC];
omega_CD=[0 0 wCD];

eqn11 = cross(omega_AB,B-A) + cross(omega_BEC,C-B) + cross(omega_CD,D-C) == 0;

loop1Solution = solve(eqn11,[wBEC wCD]);
angularVelocity_CD = double(loop1Solution.wCD) %[output:54cca0c3]
angularVelocity_BEC = double(loop1Solution.wBEC) %[output:0b11af7f]

%Loop DCEFGD

syms wGF wEF

vectoromegaCD=[0 0 angularVelocity_CD];
vectoromegaBEC=[0 0 angularVelocity_BEC];

omega_FE=[0 0 wEF];
omega_GF=[0 0 wGF];

eqn12 = cross(vectoromegaCD, C-D) + cross(vectoromegaBEC,E-C) + cross(omega_FE,F-E) + cross(omega_GF, G-F) == 0;

loop2Solution = solve(eqn12,[wEF wGF]);
angularVelocity_GF = double(loop2Solution.wGF) %[output:24fcd856]
angularVelocity_EF = double(loop2Solution.wEF) %[output:358ba2a3]

vectoromegaGF=[0 0 angularVelocity_GF] %[output:975b4d35]
vectoromegaEF= [0 0 angularVelocity_EF] %[output:4b3cc6a0]

%Angular Acceleration
%Loop 1 ABCDA
%Angular Acceleration for Loop 1 ABCDA
syms aBEC aCD
alpha_AB = [0 0 0];
alpha_BEC= [0 0 aBEC];
alpha_CD = [0 0 aCD];
a_B_A = cross(alpha_AB,B-A)+cross(omega_AB,cross(omega_AB,B-A));
a_C_B = cross(alpha_BEC,C-B) + cross(vectoromegaBEC,cross(vectoromegaBEC,C-B));
a_D_C = cross(alpha_CD,D-C) + cross(vectoromegaCD,cross(vectoromegaCD,D-C));

eqn13 = a_B_A + a_D_C + a_C_B == 0;

loop1AccSolution = solve(eqn13, [aBEC aCD]);

alphaBEC = double(loop1AccSolution.aBEC) %[output:7b337843]
alphaCD = double(loop1AccSolution.aCD) %[output:84ab136d]
vectoralphaCD = [0 0 alphaCD];
vectoralphaBEC = [0 0 alphaBEC];

%Loop 2 DCEFGD
%Angular Acceleration for Loop 2 DCEFGD

syms aEF aFG
alpha_EF = [0 0 aEF];
alpha_FG = [0 0 aFG];


a_E_F = cross(alpha_EF,F-E) + cross(vectoromegaEF,cross(vectoromegaEF,F-E));
a_C_E = cross(vectoralphaBEC,E-C) + cross(vectoromegaBEC,cross(vectoromegaBEC,E-C));
a_F_G = cross(alpha_FG,G-F) + cross(vectoromegaGF,cross(vectoromegaGF,G-F));
a_C_D = cross(vectoralphaCD,C-D) + cross(vectoromegaCD,cross(vectoromegaCD,C-D));

eqn14 = a_E_F + a_C_D + a_C_E + a_F_G == 0;

loop2AccSolution = solve(eqn14, [aEF, aFG]);

alphaEF = double(loop2AccSolution.aEF) %[output:4d5a52dd]
alphaFG = double(loop2AccSolution.aFG) %[output:32b55e49]

vectoralphaFG = [0 0 alphaFG] %[output:6d1ea489]
vectoralphaEF = [0 0 alphaEF] %[output:41e6487f]

%Velocity at Joints with respect to ground

%Velocity at B
vB_A = cross(omega_AB,B-A) %[output:0789c55c]

%Velocity at E
vE_B = cross(vectoromegaBEC,E-B);
vE_A = vE_B+vB_A %[output:6c38fab0]

%Velocity at C
vC_D = cross(vectoromegaCD,C-D)  %[output:2c3b3875]

%Velocity at F
vF_G = cross(vectoromegaGF, F-G) %[output:0e8a01b5]


%Accelerations at Joints with respect to ground

%Acceleration at E_A
a_E_A = a_B_A + cross(vectoralphaBEC,E-B)+cross(vectoromegaBEC,cross(vectoromegaBEC,E-B)) %[output:12408de1]

%Acceleration at B_A
vectora_BA = cross(alpha_AB,B-A)+cross(omega_AB,cross(omega_AB,B-A)) %[output:7126fb9a]

%Acceleration at C_D
vectora_CD = cross(vectoralphaCD,C-D)+cross(vectoromegaCD,cross(vectoromegaCD,C-D)) %[output:3fee64e9]

%Acceleration at F_G
vectora_FG = cross(vectoralphaFG,F-G)+cross(vectoromegaGF,cross(vectoromegaGF,F-G)) %[output:0ac5d130]

%Angular velocities and accelerations for each cm

%V_S4/G = V_S4_F + V_F_G

v_S4_F = cross(vectoromegaEF,S4-F);
v_S4_G = v_S4_F + vF_G %[output:51c4b2b9]

%V_S1/A = cross(omega_AB,S1-A)

v_S1_A = cross(omega_AB, S1-A) %[output:4e8940a6]

%v_S2_A = V_S2_B + V B_A

v_S2_B = cross(vectoromegaBEC,S2-B);
v_S2_A = v_S2_B+vB_A %[output:4a27f52c]

%v_S3_D = cross(vectoromegaCD,S3-D)
v_S3_D = cross(vectoromegaCD,S3-D) %[output:2733b56e]

%v_S5_G= cross(vectoromegaGF,S5-G)
v_S5_G = cross(vectoromegaGF,S5-G) %[output:924ef89b]

%a_S1_A
a_S1_A = cross(omega_AB,cross(omega_AB,S1-A)) %[output:3aa7c1fd]

%a_S2_A= vectora_BA+
a_S2_A = vectora_BA + cross(vectoralphaBEC,S2-B) + cross(vectoromegaBEC,cross(vectoromegaBEC,S2-B)) %[output:860a53f9]

%a_S3_D
a_S3_D = cross(vectoralphaCD,S3-D) + cross(vectoromegaCD,cross(vectoromegaCD,S3-D)) %[output:3984e0a5]

%a_S4_G
a_S4_G = vectora_FG + cross(vectoralphaEF,S4-F) + cross(vectoromegaEF,cross(vectoromegaEF,S4-F)) %[output:11a9d9bd]

%a_S5_G
a_S5_G = cross(vectoralphaFG,S5-G)+ cross(vectoromegaGF,cross(vectoromegaGF,S5-G)) %[output:3c53e447]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":24.8}
%---
%[output:09062f1d]
%   data: {"dataType":"textualVariable","outputData":{"name":"ForceAx","value":"-29.1509"}}
%---
%[output:1b72e175]
%   data: {"dataType":"textualVariable","outputData":{"name":"ForceAy","value":"-23.0670"}}
%---
%[output:52b7ed10]
%   data: {"dataType":"textualVariable","outputData":{"name":"ForceBx","value":"29.1509"}}
%---
%[output:3dc57759]
%   data: {"dataType":"textualVariable","outputData":{"name":"ForceBy","value":"24.0670"}}
%---
%[output:02368dd7]
%   data: {"dataType":"textualVariable","outputData":{"name":"ForceCx","value":"3.8799"}}
%---
%[output:1b178119]
%   data: {"dataType":"textualVariable","outputData":{"name":"ForceCy","value":"28.5996"}}
%---
%[output:61570b74]
%   data: {"dataType":"textualVariable","outputData":{"name":"ForceDx","value":"3.8799"}}
%---
%[output:343aa1e5]
%   data: {"dataType":"textualVariable","outputData":{"name":"ForceDy","value":"29.5996"}}
%---
%[output:4747de2e]
%   data: {"dataType":"textualVariable","outputData":{"name":"ForceEx","value":"25.2710"}}
%---
%[output:1357f864]
%   data: {"dataType":"textualVariable","outputData":{"name":"ForceEy","value":"-3.5325"}}
%---
%[output:0cfc0298]
%   data: {"dataType":"textualVariable","outputData":{"name":"ForceFx","value":"25.2710"}}
%---
%[output:661aad6f]
%   data: {"dataType":"textualVariable","outputData":{"name":"ForceFy","value":"-2.5325"}}
%---
%[output:0105cfe7]
%   data: {"dataType":"textualVariable","outputData":{"name":"ForceGx","value":"-24.7290"}}
%---
%[output:7f8e97c6]
%   data: {"dataType":"textualVariable","outputData":{"name":"ForceGy","value":"-1.5325"}}
%---
%[output:114f8990]
%   data: {"dataType":"textualVariable","outputData":{"name":"Torque","value":"396.9454"}}
%---
%[output:54cca0c3]
%   data: {"dataType":"textualVariable","outputData":{"name":"angularVelocity_CD","value":"0.9149"}}
%---
%[output:0b11af7f]
%   data: {"dataType":"textualVariable","outputData":{"name":"angularVelocity_BEC","value":"0.1915"}}
%---
%[output:24fcd856]
%   data: {"dataType":"textualVariable","outputData":{"name":"angularVelocity_GF","value":"1.0635"}}
%---
%[output:358ba2a3]
%   data: {"dataType":"textualVariable","outputData":{"name":"angularVelocity_EF","value":"-0.1047"}}
%---
%[output:975b4d35]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"vectoromegaGF","rows":1,"type":"double","value":[["0","0","1.0635"]]}}
%---
%[output:4b3cc6a0]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"vectoromegaEF","rows":1,"type":"double","value":[["0","0","-0.1047"]]}}
%---
%[output:7b337843]
%   data: {"dataType":"textualVariable","outputData":{"name":"alphaBEC","value":"-0.0328"}}
%---
%[output:84ab136d]
%   data: {"dataType":"textualVariable","outputData":{"name":"alphaCD","value":"-0.2158"}}
%---
%[output:4d5a52dd]
%   data: {"dataType":"textualVariable","outputData":{"name":"alphaEF","value":"-0.1596"}}
%---
%[output:32b55e49]
%   data: {"dataType":"textualVariable","outputData":{"name":"alphaFG","value":"0.0578"}}
%---
%[output:6d1ea489]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"vectoralphaFG","rows":1,"type":"double","value":[["0","0","0.0578"]]}}
%---
%[output:41e6487f]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"vectoralphaEF","rows":1,"type":"double","value":[["0","0","-0.1596"]]}}
%---
%[output:0789c55c]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"vB_A","rows":1,"type":"double","value":[["-12","-2","0"]]}}
%---
%[output:6c38fab0]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"vE_A","rows":1,"type":"double","value":[["-15.6383","0.4894","0"]]}}
%---
%[output:2c3b3875]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"vC_D","rows":1,"type":"double","value":[["-13.7234","1.8298","0"]]}}
%---
%[output:0e8a01b5]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"vF_G","rows":1,"type":"double","value":[["-15.9523","-2.1270","0"]]}}
%---
%[output:12408de1]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"a_E_A","rows":1,"type":"double","value":[["2.1474","-13.1237","0"]]}}
%---
%[output:7126fb9a]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"vectora_BA","rows":1,"type":"double","value":[["2","-12","0"]]}}
%---
%[output:3fee64e9]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"vectora_CD","rows":1,"type":"double","value":[["1.5623","-12.9870","0"]]}}
%---
%[output:0ac5d130]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"vectora_FG","rows":1,"type":"double","value":[["1.3948","-17.0806","0"]]}}
%---
%[output:51c4b2b9]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"v_S4_G","rows":1,"type":"double","value":[["-15.7953","-0.8188","0"]]}}
%---
%[output:4e8940a6]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"v_S1_A","rows":1,"type":"double","value":[["-6","-1","0"]]}}
%---
%[output:4a27f52c]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"v_S2_A","rows":1,"type":"double","value":[["-13.7872","0.1064","0"]]}}
%---
%[output:2733b56e]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"v_S3_D","rows":1,"type":"double","value":[["-6.8617","0.9149","0"]]}}
%---
%[output:924ef89b]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"v_S5_G","rows":1,"type":"double","value":[["-7.9761","-1.0635","0"]]}}
%---
%[output:3aa7c1fd]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"a_S1_A","rows":1,"type":"double","value":[["1","-6","0"]]}}
%---
%[output:860a53f9]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"a_S2_A","rows":1,"type":"double","value":[["1.9032","-12.7036","0"]]}}
%---
%[output:3984e0a5]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"a_S3_D","rows":1,"type":"double","value":[["0.7811","-6.4935","0"]]}}
%---
%[output:11a9d9bd]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"a_S4_G","rows":1,"type":"double","value":[["1.7711","-15.1022","0"]]}}
%---
%[output:3c53e447]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"a_S5_G","rows":1,"type":"double","value":[["0.6974","-8.5403","0"]]}}
%---
