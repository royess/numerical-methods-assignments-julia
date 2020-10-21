using LinearAlgebra

include("qr.jl")

t_arr = [-1.0, -0.75, -0.5, 0.0, 0.25, 0.5, 0.75]
y_arr = [1.00, 0.8125, 0.75, 1.00, 1.3125, 1.75, 2.3125]

A = hcat(ones(7), t_arr, t_arr.^2)

show(stdout, MIME"text/plain"(), QR.qr_solver(A, y_arr))
println()