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
