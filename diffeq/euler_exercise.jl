ENV["GKSwstype"] = "nul"

using Plots
include("euler.jl")

f(x, y) = -1. / x^2 - y / x - y ^ 2
y_truth(x) = -tan.(pi / 4. .+ log.(x)) ./ x

x1_arr, y1_arr = Euler.euler_solver(f, 0.01, 100, 1, -1)
x2_arr, y2_arr = Euler.improved_euler_solver(f, 0.01, 100, 1, -1)

plot!(x1_arr, y_truth(x1_arr), label="truth")
plot!(x1_arr, y1_arr, label="euler")
plot!(x2_arr, y2_arr, label="improved euler")
png("euler")
