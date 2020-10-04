using LinearAlgebra

include("linear_equations.jl")

A1 = Matrix{Float64}(Tridiagonal(8.0*ones(84-1), 6.0*ones(84), 1.0*ones(84-1)))
b1 = vcat([7.0], 15.0*ones(84-2), [14.0])

print("习题 1.1\n")

print("准确结果:\n")
@show x1_accurate = ones(84)
print("\n")

print("选列主元的结果:\n")
@show x1_with_selection =LinearEq.lu_solve_linear_eqs(A1, b1, true)
print("\n")

print("不选列主元的结果:\n")
@show x1_without_selection = LinearEq.lu_solve_linear_eqs(A1, b1, false)
print("\n")

print("----------\n")

A2 = Matrix{Float64}(Tridiagonal(1.0*ones(100-1), 10.0*ones(100), 1.0*ones(100-1)))
b2 = 10*rand(Float64, 100)

print("习题 1.2 (1)\n")

print("改进的平方根法的结果:\n")
@show x2_ldlt = LinearEq.sqroot_solve_linear_eqs(A2, b2, true)
@show b2-A2*x2_ldlt
print("\n")

print("平方根法的结果:\n")
@show x2_chol = LinearEq.sqroot_solve_linear_eqs(A2, b2, false)
@show b2-A2*x2_chol
print("\n")

print("习题 1.2 (2)\n")
print("由于 Hilbert 矩阵的高度病态, 这里使用 128 位的浮点数.")

A3 =  [BigFloat(1) / BigFloat(i+j-1) for i in 1:40, j in 1:40]
b3 = sum(A3, dims=2)

print("改进的平方根法的结果:\n")
@show x3_ldlt = LinearEq.sqroot_solve_linear_eqs(A3, b3, true)
@show b3-A3*x3_ldlt
print("\n")

print("平方根法的结果:\n")
@show x3_chol = LinearEq.sqroot_solve_linear_eqs(A3, b3, false)
@show b3-A3*x3_chol
print("\n")