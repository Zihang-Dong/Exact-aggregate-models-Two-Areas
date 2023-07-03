function [WK,SelectionTime] = SelectWindow(T,N,A,P,E,Di)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    Wk = ones(1,T);
    WK = Wk;
    WK_cost = Di*WK';
    for j = 1:N
        WK_cost = WK_cost + min(sum(A(j,:).*WK),E(j)/P(j))*P(j);
    end
    WK_cost = WK_cost/T;

    SelectionStart = tic;

    for i = T-1:-1:1

        Wk_Temp = [];
        Wk_AveCost = Inf;

        for n = 1:size(Wk,1)
            Wkk = Wk(n,:);

            TSet = find(Wkk == 1);

            TSubSet = nchoosek(TSet,i);
            N_row = size(TSubSet,1);
            Wk0 = zeros(N_row,T);

            if N_row > T
                Row = 1:1:N_row;
                for Column = 1:i
                    Wk0(sub2ind([N_row,T],Row,TSubSet(:,Column)')) = 1;
                end
            else
                for tt = 1:N_row
                    Wk0(tt,TSubSet(tt,:)) = 1;
                end
            end

            D_Sum = Wk0*Di';

            N0 = N;
            for k = 1:N/N0
                D_Sum_Temp = min(Wk0*A((k-1)*N0+1:(k-1)*N0+N0,:)'.*kron(ones(N_row,1),P((k-1)*N0+1:(k-1)*N0+N0)'), kron(ones(N_row,1),E((k-1)*N0+1:(k-1)*N0+N0)'));
                D_Sum = D_Sum + sum(D_Sum_Temp,2);
            end

            Ave_Cost = D_Sum./sum(Wk0,2);

            if min(Ave_Cost) < Wk_AveCost
                Wk_AveCost = min(Ave_Cost);

                Wk_Temp = Wk0(Ave_Cost == min(Ave_Cost),:);
            elseif min(Ave_Cost) == Wk_AveCost
                Wk_Temp = [Wk_Temp; Wk0(Ave_Cost == min(Ave_Cost),:)];
            else
                continue;
            end

        end

        Wk = unique(Wk_Temp,'row');

        WK = [WK; Wk];
        WK_cost = [WK_cost; Wk_AveCost];

    end

    SelectionTime = toc(SelectionStart);


end

