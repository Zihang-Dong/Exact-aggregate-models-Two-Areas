function [Obj,BatObj] = Area2_Initialization(Ts,DT,Np,P)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


% Object initialization
Obj = struct;
BatObj(1:Np) = struct;

% Area 1
Obj.InflexDemand = 10^(-3)*[26078	25697	25433	25134	25046	24781	24902	25244	25835	26576	27613	28541	29218	29321	29612	29738	29955	29318	28338	27111	25627	24199	22766	21758	20925	20527	20513	20956	20887	20520	20111	19756	19539	19439	19757	20105	21184	22002	23095	23954	24625	25125	25149	25246	24758	24307	23725	23338];

Obj.InflexDemand = Obj.InflexDemand(1:Ts/(0.5):end);

% Obj.InflexDemand = smooth(Obj.InflexDemand,'sgolay')';

Obj.Gmax = 50;
Obj.Gmin = 0;

Obj.P = P*ones(Np,1);
Obj.A = zeros(Np,DT);
Obj.E = zeros(Np,1);

Obj.A_Max = zeros(1,DT);

Obj.a = 0.02*(10^(3))^(2);
Obj.b = 14*10^(3);

% Batteries
for j = 1:Np
    
    % Power
    BatObj(j).Pr = P;

%     % Availability window
%     BatObj(j).A = zeros(1,DT);
%     C = nchoosek(1:DT,randi([1 DT]));
%     BatObj(j).A(C(randi([1 size(C,1)]),:)) = 1;
%     Obj.A(j,:) = BatObj(j).A;

    % Availability window
    BatObj(j).A = zeros(1,DT);
    BatObj(j).ts = max(0,round(5 + 1*randn(1),0));
    BatObj(j).td = min(DT-BatObj(j).ts,round(10 + 2*randn(1),0));
    BatObj(j).A(BatObj(j).ts+1:BatObj(j).ts+BatObj(j).td) = 1;
    Obj.A(j,:) = BatObj(j).A;
    
    % Maximum capacity
    Obj.A_Max = Obj.A_Max + Obj.A(j,:)*BatObj(j).Pr;
    
    % Energy
    BatObj(j).E = randi([1 sum(BatObj(j).A)])*BatObj(j).Pr;
    Obj.E(j) = BatObj(j).E;
end

end

