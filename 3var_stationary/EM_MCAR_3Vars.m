function [K_em, gamma_em, tau2_em, alpha_em, T, epsilon, time] = EM_MCAR_3Vars(S, ...
			A, num_vars, r, n, z, H, sig2eps, Delta, maxit, avgtol, opts)       
%% input:
% S: is an nxr basis function matrix.
% num_vars: number of variables used in MCAR (number of diagonal blocks in S)
% r: gives detailed information of blocks in S (r = [ r1 r2 r3 r4 ... ])
% n: gives detailed information of blocks in S (n = [ n1 n2 n3 n4 ...])
% z: is the detailed residuals.
% H: is the proximity matrix with dimension N x N
% sig2eps: is the measurem error. 
% Delta: is the NxN diagonal matrix. So far, Delta is speye(N). Users can
% generalize it to other forms.
% maxit: is the maximum number of iterations in EM algorithm.
% avgtol: the tolerance error when estimating parameters.
% opts: is an optimzation object specifying the optimzation algorithm.

N = size(H,1); % number of area units in discretized domain

% default values for the last two parameters
if nargin<11, avgtol=1e-4; end
if nargin<10, maxit=2000; end
if nargin<9 % homogeneous CAR 
	Delta = speye(N);
	DeltaInv = speye(N);
	reord = symrcm(H);
	large = eigs(H(reord,reord), 1, 'LA');
	small = eigs(H(reord,reord), 1, 'SA');
	lbg = 1/small + 1e-3*abs(1/small);
	ubg = 1/large - 1e-3*abs(1/large);
	gamma = (lbg + ubg)/2 + zeros(num_vars,1);
else % weighted CAR
	DeltaInv = spdiags(diag(Delta).^(-1), 0, N, N);
	lbg = -1+1e-3;
	ubg = 1-1e-3;
	gamma = 0.16 + zeros(num_vars,1);
end	

% initial values
for i=1:num_vars
    if i==1
        varest(i)=var( z(1:n(1)), 1 );
    else
        varest(i) =   var( z( ( sum( n(1:(i-1))) + 1):sum( n(1:i) ) ), 1);
    end         
end

for j=1:num_vars
    if j==1
       K_old = 0.95*varest(1)*eye(r(1));
    else
       K_old = blkdiag( K_old, 0.95*varest(j)*eye(r(j)) );
    end 
end

tau2 = 0.05*sum(varest)/num_vars + zeros(num_vars,1);
alpha = [ -0.1677; 0.5335; -0.0021; -0.0114; -0.0011; -0.0001];  

t = 1;
done = 0;

update = 0;
lbt = 0.01;
ubt = 10;

%[gamma; tau2; alpha]
lb = [lbg; lbg; lbg; lbt; lbt; lbt; -5; -5; -5; -5; -5; -5];
ub = [ubg; ubg; ubg; ubt; ubt; ubt; 5; 5; 5; 5; 5; 5];

if nargin < 12
	 opts = optimoptions(@fmincon,'TolX', 1e-3, 'Algorithm', ...
	   	'active-set','Display','off');
    %opts = optimset('TolX', 1e-3,'Display','off'); 
end

fun = @neg2Q;

% help term
Veps = sig2eps(1) .* speye(N);
for k=1:num_vars
    if k==1
        continue
    end
    Veps_block = sig2eps(k) .* speye(N);
    Veps = blkdiag( Veps, Veps_block);
end

AVA = A*Veps*A';
AVAInv = spdiags((diag(AVA)).^(-1), 0, sum(n), sum(n));   
A2A = A'*AVAInv*A;
A2z = A'*AVAInv*z;
A2AS = A'*AVAInv*S;

while done == 0
ts(t)=tic;

Q1 = (DeltaInv-gamma(1,t)*H)/tau2(1,t);
Q2 = (DeltaInv-gamma(2,t)*H)/tau2(2,t);
Q3 = (DeltaInv-gamma(3,t)*H)/tau2(3,t);
P1 = alpha(1,t) * speye(N) + alpha(2,t) * H;
P2 = alpha(3,t) * speye(N) + alpha(4,t) * H;
P3 = alpha(5,t) * speye(N) + alpha(6,t) * H;

M = [ Q1+P1'*Q2*P1+P2'*Q3*P2 , P2'*Q3*P3-P1'*Q2, -P2'*Q3;
      (P2'*Q3*P3-P1'*Q2)' ,      Q2+P3'*Q3*P3,   -P3'*Q3;
      -(P2'*Q3)' ,              -(P3'*Q3)',          Q3 ];

mid = M + A2A; 

SDS = S' * AVAInv*S - S'*AVAInv*A* cholmod2(mid,A2AS);    
temp = K_old\speye( sum(r) ) + SDS;

mid3z = cholmod2(mid, A2z);    

Dz = AVAInv*z - AVAInv*(A*mid3z);

% update K
muEta = K_old*(S'*Dz) - K_old*(SDS*(temp\(S'*Dz)));
SigEta = K_old - K_old*SDS*K_old' + K_old*SDS*(temp\SDS)*K_old';
K_new = SigEta + muEta*muEta';

% find numerical solution for tau2, gamma
if update == 0
	theta0 = [gamma(:,t); tau2(:,t); alpha(:,t)];
	theta = fmincon(fun, theta0, [], [], [], [], lb, ub, [], opts);
    gamma(:,t+1) = theta(1:3);
    tau2(:,t+1) = theta(4:6);
    alpha(:,t+1) = theta(7:12);
    
elseif update == 1
    gamma(:,t+1) = gamma(:,t);
    tau2(:,t+1) = tau2(:,t);
    alpha(:,t+1) = alpha(:,t);
end

% check convergence 
if sum( abs(gamma(:,t+1)-gamma(:,t)) ) + sum( abs(tau2(:,t+1)-tau2(:,t)) ) ...
        + sum( abs(alpha(:,t+1)-alpha(:,t)) ) < 1e-3 
	update = 1;
end

diff = sum(sum((K_new-K_old).^2, 1), 2) + sum( (gamma(:,t+1)-gamma(:,t)).^2 ) ...
        + sum( (tau2(:,t+1)-tau2(:,t)).^2 ) + sum( (alpha(:,t+1)-alpha(:,t)).^2 ) ;
if diff < min(avgtol*sum(r)^2,1)
	done = 1;
end

if t > maxit
	done = 1;
	disp(strcat('Algorithm did not converge after ', num2str(maxit),' iterations')); 
end
te(t) = toc(ts(t));

K_old = K_new;
t = t+1;
disp(strcat(' t= ', num2str(t)))
end  % end while


% check positive definiteness of K
E = eig(K_new);
if nnz(E<=0) > 0
	disp('Error: K is NOT Positive definite!')
end

K_em = K_new;
gamma_em = gamma(:,t);
tau2_em = tau2(:,t);
alpha_em = alpha(:,t);

if nargout > 4
	T = t;
end

if nargout > 5
	epsilon = diff;
end	

if nargout > 6
	time = te;
end

function f = neg2Q(theta)   
  
    f_Q1 = (DeltaInv-theta(1)*H)/theta(4);
    f_Q2 = (DeltaInv-theta(2)*H)/theta(5);
    f_Q3 = (DeltaInv-theta(3)*H)/theta(6);
    f_P1 = theta(7) * speye(N) + theta(8) * H;
    f_P2 = theta(9) * speye(N) + theta(10) * H;
    f_P3 = theta(11) * speye(N) + theta(12) * H;
    
    f_M = [ f_Q1+f_P1'*f_Q2*f_P1+f_P2'*f_Q3*f_P2 , f_P2'*f_Q3*f_P3-f_P1'*f_Q2, -f_P2'*f_Q3;
           (f_P2'*f_Q3*f_P3-f_P1'*f_Q2)' ,         f_Q2+f_P3'*f_Q3*f_P3,       -f_P3'*f_Q3;
            -(f_P2'*f_Q3)' ,                       -(f_P3'*f_Q3)',              f_Q3 ];
  
    Tmid = f_M + A2A; 
 
    %trick
    [Lq1, ~, sq1] = chol( f_Q1 , 'lower','vector');
    logQ1 = 2*sum(log(diag(Lq1))); clear Lq1 sq1;
    
    [Lq2, ~, sq2] = chol( f_Q2 , 'lower','vector');
    logQ2 = 2*sum(log(diag(Lq2))); clear Lq2 sq2;
    
    [Lq3, ~, sq3] = chol( f_Q3 , 'lower','vector');
    logQ3 = 2*sum(log(diag(Lq3))); clear Lq3 sq3;
    
    logM = logQ1 + logQ2 + logQ3;
    
    [Lt, ~, st] = chol(Tmid, 'lower','vector');  
    
    logT = 2*sum(log(diag(Lt)));
	part1 = -logM+logT;
	
    Tmid4S_and_Tmid3Z = cholmod2(Tmid,[A2AS A2z]); 
    
    SDnewS = S' * AVAInv*S - S'*AVAInv*A* Tmid4S_and_Tmid3Z(:, 1:size(A2AS,2) ); 
    zDnewS = z' * AVAInv * S - z' * AVAInv*A*Tmid4S_and_Tmid3Z(:, 1:size(A2AS,2) );
        
	DnewZ = AVAInv*z - AVAInv*(A*Tmid4S_and_Tmid3Z(:,end));
    clear Tmid4S_and_Tmid3Z;

	part2 = z'*DnewZ - 2*zDnewS*muEta;

	part3 = sum(sum(SDnewS.*SigEta', 2)) + muEta'*(SDnewS*muEta);
	f = part1 + part2 + part3;
  
end  % end nested function



end
