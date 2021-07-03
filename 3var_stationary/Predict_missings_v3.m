function [Y_pred] = Predict_missings_v3(S, Sp, ...
			A, Ap, z, sig2eps, n, r, n_pred, X_trend, Sigma_CAR_inv, K, aprop) 
    %% input:
    % X_trend is an (N-n)*1 vector, where N-n is the number of predicted
    % values.
    % A: (n1+n2+n3)x(N+N+N) is constraint, and is the same as what is used in EM algorithms.
    % Ap: corresponding predicted positions for X.
    % S: is a basis function matrix.
    % Sp: is the corresponding basis function matrix for prediction.
    % H: is the proximity matrix.
    % sig2eps: is the measurement error.
    % aprop: is the adjustment on sig2eps.
    % Sigma_CAR_inv, K: can be obtained directly/indirectly from EM algorithms.
   
    N = size(A,2)/3;
    
    Veps = sig2eps(1) .* speye(N);
    for k=1:3
        if k==1
            continue
        end
        Veps_block = sig2eps(k) .* speye(N);
        Veps = blkdiag( Veps, Veps_block);
    end

    AVA = A*Veps*A';
    AVAInv = spdiags((diag(AVA)).^(-1), 0, sum(n), sum(n));   %AVA is diagonal
    A2A = A'*AVAInv*A;
    A2z = A'*AVAInv*z;
    A2AS = A'*AVAInv*S;

    mid = Sigma_CAR_inv + A2A; 

    [L1, ~, s1] = chol(mid, 'lower','vector');  
    mid4S(s1, :) = L1'\(L1\A2AS(s1, :));
    DS = AVAInv*S - AVAInv*(A*mid4S);
    SDS = S'*DS; clear mid4S;
    temp = K\speye( sum(r) ) + SDS;

    mid3z(s1, 1) = L1'\( L1\A2z(s1,1) ); clear R1;            
    Dz = AVAInv*z - AVAInv*(A*mid3z);

    muEta = K*(S'*Dz) - K*(SDS*(temp\(S'*Dz)));

   
    % calculate muXi  
    tranz1 = A'*( DS*(temp\(S'*Dz)) );
    tranz2 = A' * Dz - tranz1;
    muXi = Sigma_CAR_inv\tranz2;
    
    
    
    Y1_pred = X_trend(1:n_pred(1)) + (1/aprop(1)) * Sp(1:n_pred(1),1:r(1)) * muEta(1:r(1)) + (1/aprop(1)) * Ap(1:n_pred(1),1:N) * muXi(1:N);
    Y2_pred = X_trend( (1+n_pred(1)):(n_pred(1)+n_pred(2)) ) + (1/aprop(2)) * Sp( (1+n_pred(1)):(n_pred(1)+n_pred(2)), (1+r(1)):(r(1)+r(2))) * muEta((1+r(1)):(r(1)+r(2))) + (1/aprop(2)) * Ap( (1+n_pred(1)):(n_pred(1)+n_pred(2)),(1+N):(2*N)) * muXi((1+N):(2*N));
    Y3_pred = X_trend( (1+n_pred(1)+n_pred(2)):end ) + (1/aprop(3)) * Sp( (1+n_pred(1)+n_pred(2)):end, (1+r(1)+r(2)):end) * muEta((1+r(1)+r(2)):end) + (1/aprop(3)) * Ap( (1+n_pred(1)+n_pred(2)):end,(1+2*N):(3*N)) * muXi((1+2*N):(3*N));
  
    Y_pred = [Y1_pred; Y2_pred; Y3_pred];
end