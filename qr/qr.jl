module QR

using LinearAlgebra
using BenchmarkTools

# Householder transformation
function householder(xx)
    x = copy(xx)
    n = length(x)
    eta = norm(x, Inf)
    x /= eta
    sigma = x[2:n]' * x[2:n]
    v = zeros(n)
    v[2:n] = x[2:n]

    if sigma==0
        beta = 0
    else
        alpha = sqrt(x[1]^2 + sigma)
        if x[1]<=0
            v[1] = x[1] - alpha
        else
            v[1] = -sigma / (x[1]+alpha)
        end
        beta = 2 * v[1]^2 / (sigma+v[1]^2)
        v = v/v[1]
    end

    return v, beta
end

# apply Householder transformation
function apply_house(A, v, beta)
    w = beta * A' * v
    B = A - v * w'
    return B
end

# QR decomposition using Householder transformation, m>=n
function qr_decomp!(A)
    m, n = size(A)
    @assert m>=n
    d = zeros(min(n, m-1))
    for j=1:n
        if j<m
            v, beta = householder(A[j:m, j])
            A[j:m, j:n] = apply_house(A[j:m, j:n], v, beta)
            d[j] = beta
            A[j+1:m, j] = v[2:m-j+1]
        end
    end
    return A, d
end

function qr_solver(A, b)
    m, n = size(A)

    # do qr_decomp to A
    vr = copy(A)
    vr, d = qr_decomp!(vr)
    y = copy(b)

    # c1= Q1T b
    y = copy(b)
    for j=1:length(d)
        v = vr[j:end,j]
        v[1] = 1.0
        beta = d[j]
        y[j:end] = apply_house(y[j:end], v, beta)
    end

    
    # solve Rx=c1
    for j=n:-1:2
        y[j] = y[j] / vr[j,j]
        y[1:j-1] = y[1:j-1] - y[j]*vr[1:j-1, j]
    end
    y[1] = y[1] / vr[1, 1]

    return y[1:n]
end

end
