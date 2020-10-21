using LinearAlgebra
using BenchmarkTools

include("../lu/lu.jl")
include("qr.jl")

A1 = Matrix{BigFloat}(Tridiagonal(8.0*ones(84-1), 6.0*ones(84), 1.0*ones(84-1)))
b1 = vcat([BigFloat(7.0)], BigFloat(15.0)*ones(84-2), [BigFloat(14.0)])

print("习题 1.1\n")

x1_accurate = ones(84)

println("LU分解:")
x1_with_selection = LU.lu_solver(A1, b1, true)

println("性能:")
show(stdout, MIME"text/plain"(), @benchmark LU.lu_solver(A1, b1, true))
println()

println("误差:")
show(stdout, MIME"text/plain"(), norm(x1_with_selection-x1_accurate, 2))
println()

println("QR分解:")
x1_qr = QR.qr_solver(A1, b1)

println("性能:")
show(stdout, MIME"text/plain"(), @benchmark QR.qr_solver(A1, b1))
println()

println("误差:")
show(stdout, MIME"text/plain"(), norm(x1_qr-x1_accurate, 2))
println()

println("----------")


A2 = Matrix{Float64}(Tridiagonal(1.0*ones(100-1), 10.0*ones(100), 1.0*ones(100-1)))
b2 = 10*rand(Float64, 100)

print("习题 1.2 (1)\n")

println("改进的平方根法:")
x2_ldlt = LU.sqroot_solver(A2, b2, true)

println("性能:")
show(stdout, MIME"text/plain"(), @benchmark LU.sqroot_solver(A2, b2, true))
println()

println("误差:")
show(stdout, MIME"text/plain"(), norm(b2-A2*x2_ldlt, 2))
println()

println("QR分解:")
x2_qr = QR.qr_solver(A2, b2)

println("性能:")
show(stdout, MIME"text/plain"(), @benchmark QR.qr_solver(A2, b2))
println()

println("误差:")
show(stdout, MIME"text/plain"(), norm(b2-A2*x2_qr, 2))
println()

println("----------")

print("习题 1.2 (2)\n")
print("由于 Hilbert 矩阵的高度病态, 这里使用 128 位的浮点数.")

A3 =  [BigFloat(1) / BigFloat(i+j-1) for i in 1:40, j in 1:40]
b3 = sum(A3, dims=2)

println("改进的平方根法:")
x3_ldlt = LU.sqroot_solver(A3, b3, true)

println("性能:")
show(stdout, MIME"text/plain"(), @benchmark LU.sqroot_solver(A3, b3, true))
println()

println("误差:")
show(stdout, MIME"text/plain"(), norm(b3-A3*x3_ldlt, 2))
println()

println("QR分解:")
x3_qr = QR.qr_solver(A3, b3)

println("性能:")
show(stdout, MIME"text/plain"(), @benchmark QR.qr_solver(A3, b3))
println()

println("误差:")
show(stdout, MIME"text/plain"(), norm(b3-A3*x3_qr, 2))
println()



