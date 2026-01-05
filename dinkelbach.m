function Xopt = dinkelbach(C, lambda0, tol)
%Dinkelbach  solves min Tr(CX)/(det(X)^{1/n}) s.t. Tr(X) = 1, X > 0
%   
% INPUTS: C positive definite
%         lambda0 > 0
%         tol

    [n, ~] = size(C);
    fprintf('NEW Dinkelbach: n = %i, lambda0 = %6f, tol = %e\n',n,lambda0,tol)

    c = eig(C);
    assert(all(c > 0), "C is not positive definite");

    iter = 1;

    lambda_k = lambda0;
    p = @(mu, lambda) prod(c + mu) - (lambda/n)^n;
    p_prime = @(mu) prod(c + mu) * sum(1 ./ (c + mu));

    mu_k = 0;
    % find a root mu_k of p using Newton's method
    while abs(p(mu_k, lambda_k)) > tol
        mu_k = mu_k - p(mu_k, lambda_k) / p_prime(mu_k);

        if mu_k <= -min(c) % to ensure that C + mu I is PD still
            mu_k = -min(c) + 0.05;
            break
        end
    end

    M_X = 1/(sum(1 ./ (c + mu_k)));

    X_k = M_X * inv(C + mu_k * eye(n));
    
    k_lambda = trace(C * X_k) - lambda_k * (det(X_k)^(1/n));

    fprintf('%12s','iter','mu_k','lambda_k', 'k_lambda', 'M_X');
    fprintf('\n')
    fprintf('%12i',iter)
	fprintf('%12f', mu_k, lambda_k, k_lambda, M_X)
    fprintf('\n')

    while abs(k_lambda) > tol
        iter = iter + 1;

        lambda_k = trace(C * X_k)/(det(X_k)^(1/n));

        if lambda_k < 0
            lambda_k = abs(lambda_k);
            warning("lambda_k is negative; corrected to be positive");
        end

        mu_k = 0;
        % find a root mu_k of p using Newton's method
        while abs(p(mu_k, lambda_k)) > tol
            mu_k = mu_k - p(mu_k, lambda_k) / p_prime(mu_k);
    
            if mu_k < -min(c) % to ensure that C + mu I is PD still
                mu_k = -min(c) + 0.05;
            end
        end

        M_X = 1/(sum(1 ./ (c + mu_k)));

        X_k = M_X * inv(C + mu_k * eye(n));
    
        k_lambda = trace(C * X_k) - lambda_k * (det(X_k)^(1/n));

        fprintf('%12i',iter)
	    fprintf('%12f', mu_k, lambda_k, k_lambda, M_X)
        fprintf('\n')

    end
    Xopt = X_k;
end
