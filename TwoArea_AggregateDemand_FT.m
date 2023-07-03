function [d1,d2,g1,g2,pf,fval,SelectionTime,OptimizationTime] = TwoArea_AggregateDemand_FT(T,Np1,Np2,PFmax,A1,P1,E1,Di1,Gmax1,Gmin1,a1,b1,A2,P2,E2,Di2,Gmax2,Gmin2,a2,b2)
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

% 2^T - 1 sub time intervals
W = [];
for t = 1:T
    Ind = nchoosek(1:T,t);
    WW = zeros(size(Ind,1),T);
    for tt = 1:size(Ind,1)
        WW(tt,Ind(tt,:)) = 1;
    end
    W = [W; WW];
end

% declare variables
d1 = sdpvar(T,1);
d2 = sdpvar(T,1);
g1 = sdpvar(T,1);
g2 = sdpvar(T,1);
pf = sdpvar(T,1);

% objective function
Objective = a1*(g1'*g1)+b1*sum(g1) + a2*(g2'*g2)+b2*sum(g2);

% constraints
Constraints = [d1 >= 0, d2 >= 0];

% total energy
Constraints = [Constraints, ones(1,T)*d1 == ones(1,Np1)*E1, ones(1,T)*d2 == ones(1,Np2)*E2]; 

% real time constraints
for t = 1:T
    % Power transmission
    Constraints = [Constraints, -PFmax <= pf(t) <= PFmax];
    
    % power balance
    Constraints = [Constraints, g1(t) + pf(t) == Di1(t)+d1(t), g2(t) - pf(t) == Di2(t)+d2(t)];
    
    % generator constraints
    Constraints = [Constraints, g1(t) >= Gmin1, g1(t) <= Gmax1, g2(t) >= Gmin2, g2(t) <= Gmax2];
end

% necessary and sufficient conditions
for i = 1:size(W,1)
    RHS1 = 0;
    for j = 1:Np1
        RHS1 = RHS1 + min(sum(A1(j,:).*W(i,:)),E1(j)/P1(j))*P1(j);
    end
    Constraints = [Constraints, W(i,:)*d1 <= RHS1];
    
    RHS2 = 0;
    for j = 1:Np2
        RHS2 = RHS2 + min(sum(A2(j,:).*W(i,:)),E2(j)/P2(j))*P2(j);
    end
    Constraints = [Constraints, W(i,:)*d2 <= RHS2];
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
g1 = value(g1);
d2 = value(d2);
g2 = value(g2);
pf = value(pf);
fval = value(Objective);

end

