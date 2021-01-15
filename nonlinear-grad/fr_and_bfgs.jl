using Optim
using LinearAlgebra

function f(x)
    return x[1]^2 + 4 * x[2]^2 - 4 * x[1] - 8 * x[2]
end

function g(x)
    return [2 * x[1] - 4, 8 * x[2] - 8]
end

function fr(func, grad, x0, eps)
    k = 0
    x = x0
    d = -grad(x)
    x1 = x .+ 100*eps

    while norm(x - x1) > eps
        @show k
        alpha = Optim.minimizer(optimize(a -> func(x + a .* d), [0.], BFGS()))
        @show alpha
        @show x
        @show grad(x)
        x1 = copy(x)
        x += alpha .* d
        beta = grad(x)' * grad(x) / (grad(x1)' * grad(x1))
        d = -grad(x) + beta .* d
        @show beta
        @show d
        k += 1
    end

    return x
end

fr(f, g, [0., 0.], 0.001)
