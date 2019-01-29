% example_4_2_02_270119.m
% Consider the struction composed of two quadratic bars as shown in Fig.
% 4.3. Given E = 210GPa and A = 0.003m^2, P = 10kN, and node 5 is displaced
% to the right by 0.002m, determine:
% 1. The global stiffness matrix for the structure.
% 2. The displacements at nodes 2,3, and 4.
% 3. The reaction at node 1 and 5.
% 4. The element stresses.
%
%%
clear
clc

%----------STEP 1: Declaration, assignment of given parameters:------

% Modulus of Elasticity:
E = input('Enter Modulus of Elasticity of Bar, E = ');
% If user entry is empty, it assigns a default value:
if isempty(E)
    E = 210e6;
end

% Declare no. of elements, n:
nE = input('Enter no. of quadratic elements, n = ');
% If user entry is empty, it assigns a default value:
if isempty(nE)
    nE = 1;
end

% Declare vector for varying length of bar, L, as vector:
L = ones(1,nE);
for i = 1:nE
    L(i) = input(['Enter length of element ', num2str(i), ' = ']);
end

% Declare cross-sectional area, A, as vector:
A = input('Enter Cross-sectional area, A (meter-sqr) = ');

% Display parameters:
disp('*******************');                    % Demarcation of line
disp(['Modulus of Elasticity, E = ',num2str(E),' Pa']);
disp(['Number of Quadratic elements, n = ', num2str(nE)]);
xL = 1:length(L);
msgSpec = 'Length of element, L%d = %d meters\n';
D = [xL; L];         % For organised display
fprintf(msgSpec,D);
disp(['Cross-sectional Area, A = ', num2str(A),' m^2']);

%----------STEP 2: Generating the Element Stiffness Matrix, k:-------------

k = zeros(3,3*nE);   % Preallocate each elemental stiffness matrix. However,
                    % note that this process is optional. Also note that
                    % the parameters for pre-allocation are different from
                    % those used in linear elemental stiffness matrix.
j = 1;              % Initialize iter (Important!)
for i = 1:3:3*nE
    k(1:3,i:i+2) = quadraticBarElementStiffness(E,A,L((i+2)/3));
    % Also simultaneously display each elemental stiffness matrix as though
    % they were seperate array blocks:
    disp(['k',num2str(j),' = ']); disp(k(1:3,i:i+2));
    j = j + 1;
end

%-----------STEP 3: Assembling the Global Stiffness Matrix---------
% The size of the global stiffness matrix is (n+1) x (n+1). Thus, intialize
% a zero matrix and make calls to the appropriate assemble function

% Determine dimension of the Global Stiffness Matrix using Arithmetic
% Progression expression:
% (HINT: in future program, at the begining of the selection of the no.of
% element, deploy a switch statement to pre-select the following variables
% for Linear, Quadratic, etc elements
a1 = 3;        % For quadratic element the first term, a1, corresponds to
               % the no. of node for one quadratic element
d = 2;

nN = a1 + (nE-1)*d;   % Using A.P formula to compute total no. of nodes,
                      % nN, for global element
disp(['Total no. of nodes for the quadratic element = ',num2str(nN)]);

K = zeros(nN,nN);     % Pre-allocation

c = 1;           % Initialize column range index
j = 1;           % Initialize iter

for i = 1:nE
    K = quadraticBarAssemble(K,k(1:3,c:c+2),j,j+2,j+1);
    c = c + 3; j = j + 2;
end
fprintf('\n');              % NOTE: Escape sequence (aka
                            % Control Character in MATLAB term) for newline
                            % i.e. '/n' does NOT work with disp()
disp('Global Stiffness matrix, K = '); disp(K);

%-------------STEP 4: Applying the Boundary Conditions:-----------------

% Define the Global Displacement matrix, U:
U = zeros(nN,1);        % Pre-allocation
UUnknown = U;          % Store the no. of unknown U variables

% Define the Global Nodal Force matrix, F:
F = zeros(nN,1);        % Pre-allocation
FUnknown = F;           % Store the no. of unknown F variable

% Iter
count = 0; countZeros = 0;

% Prompt user for U and F and simultaneously deploy sorting mechanism:

for i = 1:length(U)
    val = input(['Enter value of Global Displcmt., U(',num2str(i),')'...
        ' if known, or simply press "Enter" if unknown: ']);
    if isempty(val)             % Track Unknown U variables
        UUnknown(i) = i;        % Store position of unknown
    else
        U(i) = val;
        if U(i) ~= 0
            countZeros = countZeros + 1;    % Flags non-zero values present
        end
    end
end
for i = 1:length(F)
    val = input(['Enter value of Global Force, F(',num2str(i),')'...
        ' if known, or simply press "Enter" if unknown: ']);
    if isempty(val)             % Track unknown F variables
        FUnknown(i) = i;        % Store position of unknown
        count = count + 1;      % Store no. of unknown F
    else
        F(i) = val;
    end
end

fprintf('\tU\t\tF\n');
disp([U,F]);
disp('U unknowns    F unknowns')
disp([UUnknown,FUnknown]);

%-------------STEP 5: Solving the Equations:-------------------
% Using partitioning
% Test to see which parameter should be used between UUnknown and FUnknown:

% Initialize and populate the partition matrix:
rem = length(FUnknown) - count;     % Calc the remainder
disp(rem);

if (count > 0) && (countZeros > 0)
    disp('Special partitioning with superimposition required:');
    [kP,f0] = specialPartMatrix(rem, FUnknown, U, F, K);
    % Compute the corresponding partition displacement vector, uP:
    uP = kP\f0;
elseif (count > 0)
    disp('Normal parititioning required:');
    [kP,fP] = normalPartMatrix(rem, FUnknown, F, K);
    % Compute the corresponding partition displacement vector, uP:
    uP = kP\fP;
else
    disp('Error in specified boundary conditions!');
end

% Display parititioned nodal displacement:
disp('uP = '); disp(uP);

% ------------ STEP 6: Post-processing ----------------:

% ------Re-construct the original displacement vector dimension and
% integrate with the newly computed results accordingly:
uR = zeros(nN,1);
disp(uR);
rP = 1; % no reset in this context
rR = 1;

% TO DO: Investigate to see how UUnknown can be included in the partition
% matrix
for i = 1:length(uR)
    if UUnknown(i) == 0
        uR(i) = U(i);
    elseif UUnknown(i) ~= 0
        uR(i) = uP(rP);
        rP = rP + 1;
    end
end

% ------ Replace U with the final results of the Global nodal
% displacement vector:

U = uR;                 % Re-assignment

disp('Global Displacement, U = '); disp(U);

% Now, generate the Global Force vector by simply computing matrix
% multiplication:

F = K*U;
disp('Global Force, F = '); disp(F);

% ------ Setup the element nodal displacement vector by making calls to
% the function linearBarElementForces().m
% NOTE: This node positions follow the i, j, m, where m represents the
% node at the middle.

c = 1;           % Reset column range index
j = 1;           % Reset iter

u = zeros(3,nE);   % Initilize matrix to store elemental node variables
                  % NOTE: Quadratic elements always have 3 nodes hence the
                  % dimension of the row of the elemental matrix u is set
                  % as 3. Please take note.

% Compute the element stress vectors (sigma)
sigma = zeros(3,nE);    % Same as above.

% Generate the element nodal vectors (u, f) and simultaneously compute the
% element stress vectors (sigma):
disp(['Element Nodal Displacement (u) & Stress (sigma) vectors are '...
    'given as follows: ']);

for i = 1:nE
    u(:,i) = [U(j); U(j+2); U(j+1)]; % order of nodal points: i, j, m
    % I used the all member selection operator, ':', b/c MATLAB won't let
    % me use a whole number without throwing in the following error:
    % 'Assignment has more non-singleton rhs dimensions than non-singleton
    % subscripts'. IMO: This is b/c the LHS is meant to iterate over i wrt
    % to the RHS using some kind of internal range-based for loop system
    % similar to that in Python and perhaps, C++. Howover, when a constant
    % a constant value is specified, e.g. 'u(3:i) = '... MATLAB has no way
    % of performing this loop, hence the error message.
    sigma(:,i) = quadraticBarElementStresses(k(1:3,c:c+2),u(:,i),A);
    % Simultaneously display results:
    disp(['u',num2str(i),' = ']); disp(u(:,i));
    disp(['sigma',num2str(i),' = ']); disp(sigma(:,i));

    c = c + 3; j = j + 2;
end