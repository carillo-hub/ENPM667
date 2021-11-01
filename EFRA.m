function [y_k_new, K_k, Phi_k_new, P_k_new] = EFRA(X_k, P_k_curr, Phi_k_curr,y_k_curr)

    alpha = 0.5;
    lambda = 0.98;
    beta = 0.005;
    gamma = beta;

    y_k_new = X_k * Phi_k_curr'
    K_k = (alpha.*P_k_curr.*X_k) / (alpha + (X_k'.*P_k_curr.*X_k));
    Phi_k_new = Phi_k_curr + K_k.*(y_k_new-y_k_curr);
    P_k_new = ( (1/lambda)*( eye - (K_k*X_k') ) ) + (beta*eye) - ( gamma*(P_k_curr*P_k_curr) );

end 