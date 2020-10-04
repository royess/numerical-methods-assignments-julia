module LinearEq

using LinearAlgebra

function lu_decomp!(A, col_select::Bool=true)
    n = size(A)[1]
    u = zeros(Int, n-1)

    for k in 1:n-1
        if col_select
            p = findmax(abs.(A[k:n, k]))[2] + k-1
            A[[k p], :] = A[[p k], :]
            u[k] = p

            if A[k, k] == 0
                error("A is sigular.")
            end
        end

        A[k+1:n, k] = A[k+1:n, k] / A[k, k]
        A[k+1:n, k+1:n] = A[k+1:n, k+1:n] - A[k+1:n, k] * A[k, k+1:n]'
    end
    
    return A, u
end

function lu_solve_linear_eqs(A, b, col_select::Bool=true)
    # decompose A to get P, L, U
    n = size(A)[1]
    LU = copy(A)
    y = copy(b)
    LU, u = lu_decomp!(LU, col_select)

    # solve Ly=Pb

    if col_select
        for k=1:n-1
            y[[k u[k]]] = y[[u[k] k]]
        end
    end
    
    for j=1:n-1
        y[j+1:n] = y[j+1:n] - y[j]*LU[j+1:n, j]
    end

    # solve Ux=y
    for j=n:-1:2
        y[j] = y[j] / LU[j,j]
        y[1:j-1] = y[1:j-1] - y[j]*LU[1:j-1, j]
    end
    
    y[1] = y[1] / LU[1, 1]

    return y

end

function cholesky_decomp!(A)
    n = size(A)[1]

    for k=1:n
        A[k, k] = sqrt(A[k, k])
        
        A[k+1:n, k] = A[k+1:n, k] / A[k, k]

        for j=k+1:n
            A[j:n, j] = A[j:n, j] - A[j:n, k]*A[j, k] 
        end
    end

    return A
end

function ldlt_decomp!(A)
    n = size(A)[1]
    v = ones(n)

    for j=1:n
        for i=1:j-1 
            v[i] = A[j, i] * A[i, i]
        end

        A[j, j] = A[j, j] - A[j, 1:j-1]'*v[1:j-1]
        A[j+1:n, j] = (A[j+1:n, j] - A[j+1:n, 1:j-1]*v[1:j-1]) / A[j, j]
    end
    return A
end

function sqroot_solve_linear_eqs(A, b, improved=true)
    n = size(A)[1]
    y = copy(b)
    L = copy(A)
    U = zeros(n ,n)

    if improved
        L = ldlt_decomp!(L)
        D = Diagonal(L)
        L = L - D + I(n)
        U = D * L'
    else
        L = cholesky_decomp!(L)
        U = L'
    end

    for j=1:n-1
        y[j] = y[j] / L[j, j]
        y[j+1:n] = y[j+1:n] - y[j]*L[j+1:n, j]
    end
    y[n] = y[n] / L[n, n]

    for j=n:-1:2
        y[j] = y[j] / U[j,j]
        y[1:j-1] = y[1:j-1] - y[j]*U[1:j-1, j]
    end
    y[1] = y[1] / U[1, 1]

    return y
end

end