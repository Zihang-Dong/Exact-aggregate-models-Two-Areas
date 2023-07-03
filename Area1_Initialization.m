function [Obj,BatObj] = Area1_Initialization(Ts,DT,Np,P)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


% Object initialization
Obj = struct;
BatObj(1:Np) = struct;

% Area 1
Obj.InflexDemand = 10^(-3)*[30545	31311	31754	32008	31969	32095	32134	32315	32914	34107	35786	36198	36037	35341	34197	33098	31960	30948	30085	29026	28031	27000	25866	24265	23358	23242	23572	23254	22820	22427	21980	21651	21264	20983	20829	21161	21422	22421	23147	24067	24600	25736	27025	28249	29049	29570	29962	30224];

Obj.InflexDemand = Obj.InflexDemand(1:Ts/(0.5):end);

% Obj.InflexDemand = smooth(Obj.InflexDemand,'sgolay')';

Obj.Gmax = 50;
Obj.Gmin = 0;

Obj.P = P*ones(Np,1);
Obj.A = zeros(Np,DT);
Obj.E = zeros(Np,1);

Obj.A_Max = zeros(1,DT);

Obj.a = 0.01*(10^(3))^(2);
Obj.b = 15*10^(3);

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

