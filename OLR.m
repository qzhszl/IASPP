function [W,A,S_P,R] = OLR(D,k)% OLR algorithm
% Modified in 04/09/2023
% Author: IVAN JOKIC
% D demand matrix,the input parameter k is defined as an allowed deviation of the norm (D - S) from 0. For example, b = 0.7 means allowing the norm of (D - S) >= 1 - 0.7 = 0.3. This interpretation fits the IASPP problem and provides a better understanding. As b gets smaller, the deviation of the norm from 0 increases while we obtain fewer links in the resulting graph 
% W weighted adj, A is the adj, S_P distance matrix of W, 

N = size(D,1); val = 1;                                 % Define number of nodes
%% If we want to start from the graph obtained by DOR
W = ISPP_LR(D);
A = (W > 0);
%% If we start from the complete graph
% A = ones(N) - eye(N);                                   % Start from the complete graph
W = k*A.*D;
W_tilde = (W + (ones(N) - A)).^-1 - (ones(N) - A);      % Compute W tilde
S_P = distances(graph(W, 'lower'));                     % Initialise the shortest path matrix
while(nnz(S_P > D) == 0 && val > 0)                     % Remove links one by one until we exceed the constraints
    Omega = Omega_G(W_tilde);                           % Compute the effective resistance Omega
    R = A.*((Omega+eye(N)).^-1 - W_tilde).*(D-S_P);     % Compute R
    [val,~] = max(max(R));                              % Identify the maximum element
    [row,col] = find(R == val);                         % Identify the link
    A(row(1),col(1)) = 0; A(col(1),row(1)) = 0;         % Remove the link
    W = k*A.*D;
    W_tilde = (W+(ones(N)-A)).^-1-(ones(N)-A);          % Compute W tilde
    S_P = distances(graph(W, 'lower'));                 % Update the shortest path weight matrix
end
if(row(1) ~= col(1))
    A(row(1),col(1)) = 1; A(col(1),row(1)) = 1;         % Return lastly removed link
    W = k*A.*D;
    S_P = distances(graph(W, 'lower'));                 % Update the shortest path matrix
end
end


function Omega = Omega_G(W)
N = size(W,1);
Q = diag(W*ones(N,1)) - W;
[X,Mu] = eig(Q);
c=find(abs(diag((Mu)))<1e-6);
X(:,c)=[];
Mu(:,c) =[];
Mu(c,:) =[];
Mu_inv = diag(1./diag(Mu));
Qp= X*Mu_inv*X';
zeta = diag(Qp);
Omega = zeta*ones(1,N)+ones(N,1)*zeta'-2*Qp;
end

