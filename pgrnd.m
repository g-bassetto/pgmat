function X = pgrnd(b, z)
    if b < 1
        % generate samples using the truncated sum approximation
        X = pgrndsum(b, z);
    else
        intb = floor(b);
        if intb == b
            % if b is integer use the fast algoritm
            if b > 1 % sum b independet PG(1,z)
                X = 0;
                for i = 1:b
                    X = X + pgrnd1z(z);
                end
            else % sample from the fast algorithm
                X = pgrnd1z(z);
            end
        else
            frab = b - intb;
            X = pgrnd(intb, z) + pgrnd(frab, z);
        end
    end
end

function [X, iter] = pgrnd1z(z)
    t = 0.64;
    z = z / 2;
    K = 0.125 * pi^2 + 0.5 * z^2;
    p = 0.5 * pi * exp(-K*t) / K;
    q = 2 * exp(-z) * pigauss(t, 1/z);
    iter = 0;
    while true
        iter = iter + 1;
        % Generate uniform random U and V
        U = rand(1);
        V = rand(1);
        if U < p / (p + q)
            % Truncated Exponential
            X = t + rndexp() / K;
        else
            % Truncated Inverse-Gaussian
            X = rnd_trunc_ig(1/z, t);
        end
        % compare with coefficients
        S = coeff(X, 0);
        Y = V * S;
        n = 0;
        while true
            n = n + 1;
            if mod(n, 2) == 1 % n in odd
                S = S - coeff(X, n);
                if Y < S
                    X = 0.25 * X;
                    return;
                end
            else % n is even
                S = S + coeff(X, n);
                if Y > S
                    break;
                end
            end
        end % INNER WHILE
    end % OUTER WHILE
    
    function a = coeff(x, n)
        c = pi * (n + 0.5);
        if x <= t
            a_ = exp(-2*(n + 0.5)^2 / x) * sqrt(8 / (pi * x)^3);
        else
            a_ = exp(-0.5 * x * c^2);
        end
        a = c * a_;
    end
end

function X = pgrndsum(b, z)
    K = 200; % number of terms in the sum
    k = (1:K)';
    % the weights in the sum
    w = 1./((k - 0.5).^2 + 0.25 * (z / pi)^2);
    % the iid gamma variables
    g = gamrnd(b, 1, K, 1);
    X = 0.5 * (w'* g) / pi^2;
end

function X = rnd_trunc_ig(mu, t)
    %RND_TRUNC_IG samples a truncate inverse-Gaussian distribution IG(m, 1)
    %in the interval (0, t)
    if mu > t
        X = rnd_ig_exp(mu, t);
    else
        X = rnd_ig_chi(mu, t);
    end
end

function X = rnd_ig_exp(mu, t)
    %RND_IG_EXP samples a truncate inverse-Gaussian distribution IG(m, 1)
    %in the interval (0, t) when m > t
    
    z = 1 / mu;
    exitflag = false;
    
    while ~exitflag
        exitflag2 = false;
        while ~exitflag2
            E = rndexp();
            E1 = rndexp();
            exitflag2 = E^2 <= 2 * E1 / t;
        end
        X = t / (1 + t*E)^2;
        a = exp(-0.5 * X * z^2);
        U = rand(1);
        exitflag = U <= a;
    end
end

function X = rnd_ig_chi(mu, t)
    %RND_IG_CHI samples a truncate inverse-Gaussian distribution IG(m, 1)
    %in the interval (0, t) when m <= t
    exitflag = false;
    while ~exitflag
        Y = randn(1)^2;
        mY = mu * Y;
        X = mu + 0.5 * mu * (mY - sqrt(4*mY + mY^2));
        U = rand(1);
        if U > mu / (mu + X)
            X = mu^2 / X;
        end
        exitflag = X <= t;
    end
end

function y = pigauss(x, mu)
    x05 = sqrt(x);
    y = normcdf((x/mu - 1) / x05) + exp(2/mu) * normcdf(-(x/mu + 1) / x05);
end