N = 1e4;

% tic
% [X1, M1] = arrayfun(@(~) rndpg1z(0.5), 1:N);
% toc

[B, C] = ndgrid(logspace(-1,1,3), logspace(-1.3,1.3,3));
T = zeros(9,1);
for i = 1:9
    b = B(i);
    c = C(i);
    X = zeros(N, 1);
    t = tic();
    for n = 1:N
        X(n) = pgrnd(b, c);
    end
    t(i) = toc(t);
    subplot(3,3,i);
    hist(X, linspace(0, max(X(:)), 50));
    xlim([0, max(X(:))]);
    title(sprintf('b = %0.2f, c = %0.2f', b, c))
end