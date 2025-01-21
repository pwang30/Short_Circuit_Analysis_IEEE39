clear;
clc;                                                                       % Clear data and command window
mpc = IEEE39;                                                              % Reading data files 
tic;                                                                       % Start the timer       
f = input('Enter the node f where the short circuit fault will occur:');   
branch = sparse(mpc.branch);                                               % Get branch data and create a sparse matrix
generator = sparse(mpc.gen);                                               % Get generator data and create a sparse matrix
a = -0.5+sqrt(3)*1i/2;                                                     % Twiddle Factor，a=e^(j120)
[b,c] = size(branch);                                                      % Get the number of rows and columns of the branch data
nodes = 39;                                                                % Number of nodes in IEEE39 system39

%% Forming the positive sequence admittance matrix
bBij1 = -1./branch(:,5);                                                   % Positive sequence admittance of branch
gBij1 = -1./generator(:,2);                                                % Positive sequence admittance of generators
B1 = sparse(branch(:,3),branch(:,4),-bBij1,nodes,nodes);                   % Upper triangular non-object element mutual admittance
B1 = B1 + sparse(branch(:,4),branch(:,3),-bBij1,nodes,nodes);              % Upper triangular non-object element mutual admittance
B1 = B1 + sparse(branch(:,3),branch(:,3),bBij1,nodes,nodes);               % Main diagonal element self-admittance
B1 = B1 + sparse(branch(:,4),branch(:,4),bBij1,nodes,nodes);
B1 = B1 + sparse(generator(:,1),generator(:,1),gBij1,nodes,nodes);
Y1 = 1i*B1;

%% Calculate the fth column element of the node impedance matrix related to the short-circuit point f using the positive-sequence node admittance matrix
If1 = zeros(nodes,1);                                                      % 39 rows and 1 column of all zero matrix
If1(f) =1 ;                                                                % The short-circuit node is set to 1
Z1 = Y1\If1;                                                               % Apply left division to the matrix and obtain the U vector, which gives the Z vector +

%% Three-phase short circuit
    disp('---------------Three-phase short circuit---------------');
    % Calculate the three-phase short-circuit current at the fault point
    If_1 = 1/Z1(f);
    disp('Three-phase short-circuit current at fault:');
    disp(If_1);
    % Non-fault node short circuit voltage
    U_1=1-Z1*If_1;
    disp('The voltage at the non-fault node is:');
    [d1,f1] = size(U_1);
    for i = 1:d1
        fprintf('  %d    %f\n',i,U_1(i));
    end
    %disp(U_1);
    % Calculate the short-circuit current of each branch
    U_ij1 = sparse (diag(U_1)*ones(nodes,nodes)-ones(nodes,nodes)*diag(U_1));
    I_ij1 = U_ij1.*Y1;
    disp('The current of each branch at the non-fault position of three-phase short circuit is:')
    disp(sparse(I_ij1));    

toc;