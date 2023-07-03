clear;
close all;
clc;
yalmip('clear')

rng('default');

%% Setup
% Sampling time
T = 24;
Ts = 1;
DT = T/Ts;
TimeLine = 0:Ts:T-Ts;

% Number of batteries
Np1 = 4000; % 4*10^(6)
Np2 = 6000; % 6*10^(6)

% Rated charging power
P = 0.005; % 5*10^(-6)

% Maximum power flow transmission
PFmax = 20;

% Initialization
[Obj1,BatObj1] = Area1_Initialization(Ts,DT,Np1,P);
[Obj2,BatObj2] = Area2_Initialization(Ts,DT,Np2,P);

% load('Initialization.mat');

%% Individual agent model based optimization

[d1_M1,d2_M1,g1_M1,g2_M1,pf_M1,fval_M1,SelectionTime_M1,OptimizationTime_M1] = TwoArea_AgentBasedModel(DT,Np1,Np2,PFmax,Obj1.A,Obj1.P,Obj1.E,Obj1.InflexDemand,Obj1.Gmax,Obj1.Gmin,Obj1.a,Obj1.b,Obj2.A,Obj2.P,Obj2.E,Obj2.InflexDemand,Obj2.Gmax,Obj2.Gmin,Obj2.a,Obj2.b,Obj1.A_Max,Obj2.A_Max);
% load('M1_Cap0.mat');
% load('M1_Cap5.mat');
% load('M1_Cap10.mat');
% load('M1_Cap15.mat');
% load('M1_Cap20.mat');
MP1_M1 = 2*Obj1.a*(g1_M1) + Obj1.b;
MP2_M1 = 2*Obj2.a*(g2_M1) + Obj2.b;
fprintf('\nIndividual agent model optimization: finished! \n');

%% Aggregate demand model based optimization - Full time window

d1_M2 = d1_M1; d2_M2 = d2_M1; g1_M2 = g1_M1; g2_M2 = g2_M1; pf_M2 = pf_M1; MP1_M2 = MP1_M1; MP2_M2 = MP2_M1;
fval_M2 = fval_M1;
% [d1_M2,d2_M2,g1_M2,g2_M2,pf_M2,fval_M2,SelectionTime_M2,OptimizationTime_M2] = TwoArea_AggregateDemand_FT(DT,Np1,Np2,PFmax,Obj1.A,Obj1.P,Obj1.E,Obj1.InflexDemand,Obj1.Gmax,Obj1.Gmin,Obj1.a,Obj1.b,Obj2.A,Obj2.P,Obj2.E,Obj2.InflexDemand,Obj2.Gmax,Obj2.Gmin,Obj2.a,Obj2.b);
fprintf('\nAggregate demand model based optimization - Full time window: finished! \n');

%% Aggregate demand model based optimization - Individual selected time window

[d1_M3,d2_M3,g1_M3,g2_M3,pf_M3,fval_M3,SelectionTime_M3,OptimizationTime_M3] = TwoArea_AggregateDemand_IST(DT,Np1,Np2,PFmax,Obj1.A,Obj1.P,Obj1.E,Obj1.InflexDemand,Obj1.Gmax,Obj1.Gmin,Obj1.a,Obj1.b,Obj2.A,Obj2.P,Obj2.E,Obj2.InflexDemand,Obj2.Gmax,Obj2.Gmin,Obj2.a,Obj2.b,Obj1.A_Max,Obj2.A_Max);
MP1_M3 = 2*Obj1.a*(g1_M3) + Obj1.b;
MP2_M3 = 2*Obj2.a*(g2_M3) + Obj2.b;
fprintf('\nAggregate demand model based optimization - Individual selected time window: finished! \n');

%% Aggregate demand model based optimization - Merged selected time window

[d1_M4,d2_M4,g1_M4,g2_M4,pf_M4,fval_M4,SelectionTime_M4,OptimizationTime_M4] = TwoArea_AggregateDemand_MST(DT,Np1,Np2,PFmax,Obj1.A,Obj1.P,Obj1.E,Obj1.InflexDemand,Obj1.Gmax,Obj1.Gmin,Obj1.a,Obj1.b,Obj2.A,Obj2.P,Obj2.E,Obj2.InflexDemand,Obj2.Gmax,Obj2.Gmin,Obj2.a,Obj2.b,Obj1.A_Max,Obj2.A_Max);
MP1_M4 = 2*Obj1.a*(g1_M4) + Obj1.b;
MP2_M4 = 2*Obj2.a*(g2_M4) + Obj2.b;
fprintf('\nAggregate demand model based optimization - Merged selected time window: finished! \n');

%% Export to excel

% Output data

ExcelDataInput = table(d1_M1,d2_M1,g1_M1,g2_M1,pf_M1,MP1_M1,MP2_M1);
writetable(ExcelDataInput,'Results.xlsx','Sheet','AgentModel')

ExcelDataInput = table(d1_M2,d2_M2,g1_M2,g2_M2,pf_M2,MP1_M2,MP2_M2);
writetable(ExcelDataInput,'Results.xlsx','Sheet','AggregateModel_FT')

ExcelDataInput = table(d1_M3,d2_M3,g1_M3,g2_M3,pf_M3,MP1_M3,MP2_M3);
writetable(ExcelDataInput,'Results.xlsx','Sheet','AggregateModel_IST')

ExcelDataInput = table(d1_M4,d2_M4,g1_M4,g2_M4,pf_M4,MP1_M4,MP2_M4);
writetable(ExcelDataInput,'Results.xlsx','Sheet','AggregateModel_MST')


% Input data
Di1 = Obj1.InflexDemand';
Di2 = Obj2.InflexDemand';

d_Max1 = Obj1.A_Max';
d_Max2 = Obj2.A_Max';

PF_Max = PFmax*ones(DT,1);
PF_Min = -PFmax*ones(DT,1);

ExcelDataInput = table(Di1,Di2,d_Max1,d_Max2,PF_Max,PF_Min);
writetable(ExcelDataInput,'Results.xlsx','Sheet','Inflexible')