ENV["GKSwstype"] = "nul"
using Plots

include("runge_kutta.jl")

function generator(σ, ρ, β)
    f(t, y) = [
        σ * (y[2] - y[1]),
        ρ * y[1] - y[2] - y[1] * y[3],
        y[1] * y[2] - β * y[3]]
    return f
end


t_arr, y_arr = RungeKutta.runge_kutta_solver(generator(10,28,8/3), 5000, 0.01, 0.0, [1.0, 1.0, 1.0])
p1 = scatter(y_arr[:,1], y_arr[:,2], y_arr[:,3], marker=1.6, label="[1.0, 1.0, 1.0]")

t_arr, y_arr = RungeKutta.runge_kutta_solver(generator(10,28,8/3), 5000, 0.01, 0.0, [3.0, -1.0, -1.0])
p2 = scatter(y_arr[:,1], y_arr[:,2], y_arr[:,3], marker=1.6, label="[3.0, -1.0, -1.0]")

t_arr, y_arr = RungeKutta.runge_kutta_solver(generator(10,28,8/3), 5000, 0.01, 0.0, [-1.0, 3.0, -1.0])
p3 = scatter(y_arr[:,1], y_arr[:,2], y_arr[:,3], marker=1.6, label="[-1.0, 3.0, -1.0]")

t_arr, y_arr = RungeKutta.runge_kutta_solver(generator(10,28,8/3), 5000, 0.01, 0.0, [-1.0, -1.0, 3.0])
p4 = scatter(y_arr[:,1], y_arr[:,2], y_arr[:,3], marker=1.6, label="[-1.0, -1.0, 3.0]")

f1 = plot(p1, p2, p3, p4, layout=4, size=(800,800))

png(f1, "lorenz-location")

t_arr, y_arr = RungeKutta.runge_kutta_solver(generator(10,30,8/3), 5000, 0.01, 0.0, [1.0, 1.0, 1.0])
p5 = scatter(y_arr[:,1], y_arr[:,2], y_arr[:,3], marker=1.6, label="ρ->30")

t_arr, y_arr = RungeKutta.runge_kutta_solver(generator(10,26,8/3), 5000, 0.01, 0.0, [1.0, 1.0, 1.0])
p6 = scatter(y_arr[:,1], y_arr[:,2], y_arr[:,3], marker=1.6, label="ρ->26")

t_arr, y_arr = RungeKutta.runge_kutta_solver(generator(10,28,9/3), 5000, 0.01, 0.0, [1.0, 1.0, 1.0])
p7 = scatter(y_arr[:,1], y_arr[:,2], y_arr[:,3], marker=1.6, label="β->9/3")

t_arr, y_arr = RungeKutta.runge_kutta_solver(generator(10,28,7/3), 5000, 0.01, 0.0, [1.0, 1.0, 1.0])
p8 = scatter(y_arr[:,1], y_arr[:,2], y_arr[:,3], marker=1.6, label="β->7/3")

t_arr, y_arr = RungeKutta.runge_kutta_solver(generator(11,28,8/3), 5000, 0.01, 0.0, [1.0, 1.0, 1.0])
p9 = scatter(y_arr[:,1], y_arr[:,2], y_arr[:,3], marker=1.6, label="σ->11")

t_arr, y_arr = RungeKutta.runge_kutta_solver(generator(9,28,8/3), 5000, 0.01, 0.0, [1.0, 1.0, 1.0])
p10 = scatter(y_arr[:,1], y_arr[:,2], y_arr[:,3], marker=1.6, label="σ->9")

f2 = plot(p5, p6, p7, p8, p9, p10, layout=6, size=(1200,800))

png(f2, "lorenz-parameter")
