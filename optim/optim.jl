module Opt

using Optim
using LinearAlgebra

function sd_optimize(func, grad, x0, eps)
    k = 0
    x = copy(x0)   
    x1 = x .+ 100 * eps

    while norm(x - x1) > eps
        x1 = copy(x)
        d = -grad(x)
        alpha = Optim.minimizer(optimize(a -> func(x + a .* d), [0.], BFGS()))
        x += alpha .* d
        k += 1
    end

    return x, func(x)
end


function newton_optimize(func, grad, jacob, x0, eps)
    k = 0
    x = copy(x0)   
    x1 = x .+ 100 * eps

    while norm(x - x1) > eps
        x1 = copy(x)
        d = -inv(jacob(x)) * grad(x) 
        x += d
        k += 1
    end
    
    return x, func(x)
end


function quasi_newton_optimize(func, grad, x0, eps, kind)
    @assert kind in ["SR1", "BGFS", "DFP"]

    k = 0
    x = copy(x0)
    H = I(size(x)[1]) # here we take H0 to be identity
    x1 = x .+ 100 * eps
    g1 = zeros(size(x)[1])

    while norm(x - x1) > eps
        x1 = copy(x); g = copy(g1)
        d = -H * g
        alpha = Optim.minimizer(optimize(a -> func(x + a .* d), [0.], BFGS()))

        s = alpha .* d
        x += s

        g1 = grad(x)
        y = g1 - g

        # update H
        if kind == "SR1"
            H += (s - H * y) * (s - H * y)' / ((s - H * y)' * y)
        elseif kind == "BGFS"
            H += (1 + y' * H * y / (y' * s)) * (s * s') / (y' * s) - (s * y' * H + H * y * s') / (y' * s)
        else # DFP
            H += s * s' / (s' * y) - H * y * y' * H / (y' * H * y)
        end

        k += 1
    end
    
    return x, func(x)
end


function bb_optimize(func, grad, x0, eps, kind)
    @assert kind in ["BB1" "BB2"]
    
    k = 0
    x = copy(x0)   
    x1 = x .+ 100 * eps
    g1 = zeros(size(x)[1])

    while norm(x - x1) > eps
        x1 = copy(x)
        d = -grad(x)
        
        if kind == "BB1"
            alpha = g1' * g / (g' * G * g)
        else
            alpha = g1' * G * g / (g' * G * G * g)
        end

        x += alpha .* d
        g1 = -copy(d)

        k += 1
    end
    
    return x, func(x)
end


function fr_optimize(func, grad, x0, eps)
    k = 0
    x = x0
    d = -grad(x)
    x1 = x .+ 100 * eps

    while norm(x - x1) > eps
        alpha = Optim.minimizer(optimize(a -> func(x + a .* d), [0.], BFGS()))
        x1 = copy(x)
        x += alpha .* d
        beta = grad(x)' * grad(x) / (grad(x1)' * grad(x1))
        d = -grad(x) + beta .* d
        k += 1
    end

    return x
end

end
