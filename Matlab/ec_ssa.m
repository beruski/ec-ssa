function res = ec_ssa(t,X,R,P,C,BV,U,I,K)
% ec_ssa(t,X,R,P,C,BV,U,I,K)
%
% Electrochemical SSA v5. Evaluates the stochastic time evolution of a
% system with M reaction channels involving N species. Based on the work
% described in:
%
% Beruski, O. Stochastic Electrochemical Kinetics. arXiv:1608.07507
%
% The following input parameters are expected:
%
% t is the time vector: it may be a scalar, 2x1 or 3x1 vector.
% Options in the code.
%
% X is the initial species population: it is a Nx1 vector.
%
% R and P are the reactant and product stoichiometry matrix, respectively:
% they are MxN matrices, with the corresponding stoichiometric coefficient
% for species j and channel i at R(i,j). R(i,j) and P(i,j) >= 0.
%
% C is the transition constants vectors: it is a Mx1 vector.
%
% BV is the Butler-Volmer coefficients matrix: its is an 3xM matrix.
% Line 1 contains the electron stoichiometry coeffients, with BV(1,i) > 0
% for oxidation, and BV(1,i) < 0 for reduction channels. BV(1,i) = 0
% signals a non-electrochemical channel.
% Line 2 contains the transfer coefficients, usually 0 < BV(2,i) <= 1.0.
% Line 3 contains the formal potentials.
%
% U is the applied electrode potential: it is a scalar. Relevant for K = 1.
%
% I is the applied electric current: it is a scalar. Relevant for K = 2.
%
% K is the flag determining the type of run: (1) for potentiostatic, (2)
% for galvanostatic, and (0) for non-electrochemical.
%
% Implemented by Beruski, O.

if length(t) == 1
    % If t is a scalar, a default linear space is created for the sampling
    % times.
    st = linspace(0,t);
    clear t;
elseif length(t) == 2
    % If t has 2 components, t(1) gives the stopping time, while
    % t(2) gives the length of the sampling time vector.
    st = linspace(0,t(1),t(2));
    clear t;
elseif length(t) >= 3
    % If length(t) >= 3, the function treats it as the sampling time
    % vector.
    st = t;
    clear t;
else
    % In the case t is invalid.
    disp('Problem with input t: invalid value.');
    disp('Aborting...');
    res = -1;
    return;
end

[M N] = size(R);
S = length(st);
t_max = st(S);

% Parameters

q_e = 96485/6.022E23;
T = 298.15;
k_B = 8.3144621/6.022E23;
f = q_e/(k_B*T);

THRESH = 1E-6;

% Checking for input problems

if size(P) ~= [M N]
    disp('Problem with inputs R and/or P: different sizes.');
    disp('Aborting...');
    res = -1;
    return;
end

if size(X) ~= [1 N]
    disp('Problem with input X: wrong size.');
    disp('Aborting...');
    res = -1;
    return;
end

if size(C) ~= [1 M]
    disp('Problem with input C: wrong size.');
    disp('Aborting...');
    res = -1;
    return;
end

if K == 1 || K == 2
    if size(BV) ~= [3 M]
        disp('Problem with input BV: wrong size.');
        disp('Aborting...');
        res = -1;
        return;
    end
end

% Declaring and initializing variables

s = 1;   % Sample number
t = 0;   % Total time
tau = 0; % Sampled time for next reaction
mu = 0;  % Sampled channel fired
A = 0;   % Total reaction propensity

a = zeros(1,M); % Individual reaction propensities
T = zeros(S,1); % Sampled times
Y = zeros(S,N); % Sampled population numbers
Z = zeros(S,M); % Sampled reaction propensities

if K == 1
    J = zeros(S,1); % Sample electric current
elseif K == 2
    E = zeros(S,1); % Sampled electrode potential
end

% Dummy variables

d = 0;
x = 0;
dx = 0;
DU = 0;

if K == 0
    while t <= t_max
        % Resetting variables
        a(:) = C(:);
        A = 0;
        d = 0;
        % Calculating reaction propensities
        for i = 1:M
            for j = 1:N
                if R(i,j) == 1
                    a(i) = a(i)*X(j);
                elseif R(i,j) == 2
                    a(i) = a(i)*X(j)*(X(j)-1)/2;
                elseif R(i,j) == 3
                    a(i) = a(i)*X(j)*(X(j)-1)*(X(j)-2)/6;
                elseif R(i,j) > 3
                    disp('Reaction order not implemented for indeces:');
                    disp([i j]);
                    disp('Consider doing so after line 113');
                end
            end
            A = A + a(i);
        end
        % Writing initial conditions
        if s == 1
            Y(s,:) = X(:);
            Z(s,:) = a(:);
        end
        % Checking for equillibrium
        if A == 0
            T(s) = t;
            Y(s,:) = X(:);
            disp('Zero total propensity. Terminating SSA.');
            break;
        end
        % Sampling (tau,mu)
        tau = (1/A)*log(1/rand);
        mu = rand*A;
        i = 0;
        while d < mu
            i = i + 1;
            d = d + a(i);
        end
        mu = i;
        % Applying modifications
        t = t + tau;
        X = X - R(mu,:) + P(mu,:);
        % Checking for sampling time
        if t  >= st(s)
            s = s + 1;
            T(s) = t;
            Y(s,:) = X(:);
            Z(s,:) = a(:);
        end
    end
elseif K == 1
    while t <= t_max
        % Resetting variables
        a(:) = C(:);
        A = 0;
        I = 0;
        d = 0;
        % Calculating reaction propensities
        for i = 1:M
            for j = 1:N
                if R(i,j) == 1
                    a(i) = a(i)*X(j);
                elseif R(i,j) == 2
                    a(i) = a(i)*X(j)*(X(j)-1)/2;
                elseif R(i,j) == 3
                    a(i) = a(i)*X(j)*(X(j)-1)*(X(j)-2)/6;
                elseif R(i,j) > 3
                    disp('Reaction order not implemented for indeces:');
                    disp([i j]);
                    disp('Consider doing so after line 175');
                end
            end
            a(i) = a(i)*exp(BV(1,i)*BV(2,i)*f*(U-BV(3,i)));
            A = A + a(i);
            I = I + BV(1,i)*q_e*a(i);
        end
        % Writing initial conditions
        if s == 1
            Y(s,:) = X(:);
            Z(s,:) = a(:);
        end
        % Checking for equillibrium
        if A == 0
            T(s) = t;
            Y(s,:) = X(:);
            J(s) = I;
            disp('Zero total propensity. Terminating SSA.');
            break;
        end
        % Sampling (tau,mu)
        tau = (1/A)*log(1/rand);
        mu = rand*A;
        i = 0;
        while d < mu
            i = i + 1;
            d = d + a(i);
        end
        mu = i;
        % Applying modifications
        t = t + tau;
        X = X - R(mu,:) + P(mu,:);
        % Checking for sampling time
        if t  >= st(s)
            s = s + 1;
            T(s) = t;
            Y(s,:) = X(:);
            Z(s,:) = a(:);
            J(s) = I;
        end
    end
elseif K == 2
    while t <= t_max
        % Resetting variables
        a(:) = C(:);
        A = 0;
        U = sign(I)*1.0;
        DU = 1.0;
        d = 0;
        % Calculating reaction quocient
        for i = 1:M
            for j = 1:N
                if R(i,j) == 1
                    a(i) = a(i)*X(j);
                elseif R(i,j) == 2
                    a(i) = a(i)*X(j)*(X(j)-1)/2;
                elseif R(i,j) == 3
                    a(i) = a(i)*X(j)*(X(j)-1)*(X(j)-2)/6;
                elseif R(i,j) > 3
                    disp('Reaction order not implemented for indeces:');
                    disp([i j]);
                    disp('Consider doing so after lines 236');
                end
            end
        end
        % Calculating electrode potential
        while abs(DU) >= THRESH
            x = -I/q_e;
            dx = 0.0;
            for i = 1:M
                x = x + BV(1,i)*a(i)*exp(BV(1,i)*BV(2,i)*f*(U-BV(3,i)));
                dx = dx + f*BV(1,i)^2*BV(2,i)*a(i)*exp(BV(1,i)*BV(2,i)*f*(U-BV(3,i)));
            end
            DU = x/dx;
            U = U - DU;
        end
        % Calculating reaction propensities
        for i = 1:M
            a(i) = a(i)*exp(BV(1,i)*BV(2,i)*f*(U-BV(3,i)));
            A = A + a(i);
        end
        % Writing initial conditions
        if s == 1
            Y(s,:) = X(:);
            Z(s,:) = a(:);
        end
        % Checking for equillibrium
        if A == 0 || isnan(A) == 1
            T(s) = t;
            Y(s,:) = X(:);
            E(s) = U;
            disp('Zero total propensity. Terminating SSA.');
            break;
        end
        % Sampling (tau,mu)
        tau = (1/A)*log(1/rand);
        mu = rand*A;
        i = 0;
        while d < mu
            i = i + 1;
            d = d + a(i);
        end
        mu = i;
        % Applying modifications
        t = t + tau;
        X = X - R(mu,:) + P(mu,:);
        % Checking for sampling time
        if t  >= st(s)
            s = s + 1;
            T(s) = t;
            Y(s,:) = X(:);
            Z(s,:) = a(:);
            E(s) = U;
        end
    end
end

if K == 1 || K == 2
    subplot(2,1,1)
end
plot(T(1:s-1),Y(1:s-1,:),'o-')
xlabel('\bf{Time (s)}')
ylabel('\bf{Population Number}')
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
if K == 0
    res = [T Y Z];
elseif K == 1
    subplot(2,1,2)
    plot(T(1:s-1),J(1:s-1),'ko-')
    xlabel('\bf{Time (s)}')
    ylabel('\bf{Electric Current (A)}')
    set(gca,'XMinorTick','on')
    set(gca,'YMinorTick','on')
    res = [T J Y Z];
elseif K == 2
    subplot(2,1,2)
    plot(T(1:s-1),E(1:s-1),'ko-')
    xlabel('\bf{Time (s)}')
    ylabel('\bf{Electrode Potential (V)}')
    set(gca,'XMinorTick','on')
    set(gca,'YMinorTick','on')
    res = [T E Y Z];
end
end