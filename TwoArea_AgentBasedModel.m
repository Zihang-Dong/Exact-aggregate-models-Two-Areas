function [d1,d2,g1,g2,pf,fval,SelectionTime,OptimizationTime] = TwoArea_AgentBasedModel(T,Np1,Np2,PFmax,A1,P1,E1,Di1,Gmax1,Gmin1,a1,b1,A2,P2,E2,Di2,Gmax2,Gmin2,a2,b2,d1_max,d2_max)
%%

% optimization
% options = sdpsettings('solver','fmincon','verbose',0);
% options.fmincon.Algorithm = 'sqp';
% options.fmincon.MaxIter = 10^(5);
% options.fmincon.MaxFunEvals = 10^(6);
% options.fmincon.TolCon = 1.0000e-12;
% options.fmincon.TolFun = 1.0000e-12;
% options.fmincon.TolFunValue = 1.0000e-12;

% options = sdpsettings('solver','gurobi','verbose',0);
% options.gurobi.IntFeasTol = 1.0000e-9;
% options.gurobi.FeasibilityTol = 1.0000e-9;
% options.gurobi.OptimalityTol = 1.0000e-9;
% options.gurobi.MIPGap = 1.0000e-9;
% options.gurobi.PSDTol = 1.0000e-9;
% options.gurobi.PerturbValue = 1.0000e-9;

options = sdpsettings('solver','quadprog','verbose',0);
% options.quadprog.TolCon = 1.0000e-12;
% options.quadprog.TolFun = 1.0000e-12;
% options.quadprog.TolFunValue = 1.0000e-12;

SelectionTime = 0;

%% Optimization problem

% declare variables
d1 = sdpvar(1,T);
d2 = sdpvar(1,T);
g1 = sdpvar(1,T);
g2 = sdpvar(1,T);
pf = sdpvar(1,T);
u1 = sdpvar(Np1,T);
u2 = sdpvar(Np2,T);

% objective function
Objective = a1*(g1*g1')+b1*sum(g1) + a2*(g2*g2')+b2*sum(g2);

% constraints
Constraints = [d1 >= 0, d2 >= 0];

% real time constraints
for t = 1:T
    % Power transmission
    Constraints = [Constraints, -PFmax <= pf(t) <= PFmax];
    
    % power balance
    Constraints = [Constraints, g1(t) + pf(t) == Di1(t)+d1(t), g2(t) - pf(t) == Di2(t)+d2(t)];
    
    % generator constraints
    Constraints = [Constraints, g1(t) >= Gmin1, g1(t) <= Gmax1, g2(t) >= Gmin2, g2(t) <= Gmax2];
    
    % aggregate demand
    Constraints = [Constraints, d1(t) == sum(u1(:,t)), d2(t) == sum(u2(:,t))];
    
    % power
    for j = 1:Np1
        Constraints = [Constraints, 0 <= u1(j,t) <= P1(j).*A1(j,t)];
    end
    for j = 1:Np2
        Constraints = [Constraints, 0 <= u2(j,t) <= P2(j).*A2(j,t)];
    end
    
    % maximum power
    Constraints = [Constraints, d1(t) <= d1_max(t), d2(t) <= d2_max(t)];
    
end

% energy constraint
for j = 1:Np1
    Constraints = [Constraints, u1(j,:)*ones(T,1) == E1(j)];
end
for j = 1:Np2
    Constraints = [Constraints, u2(j,:)*ones(T,1) == E2(j)];
end
        
% optimization problem
OptimizationStart= tic;
diagnostics = optimize(Constraints,Objective,options);
OptimizationTime = toc(OptimizationStart);

if diagnostics.problem ~= 0
    error('Something else happened')
end


%% optimal solutions
d1 = value(d1);
d1 = d1';
g1 = value(g1);
g1 = g1';
d2 = value(d2);
d2 = d2';
g2 = value(g2);
g2 = g2';
pf = value(pf);
pf = pf';
fval = value(Objective);

end

