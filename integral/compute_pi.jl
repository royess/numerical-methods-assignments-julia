ENV["GKSwstype"] = "nul"
using Plots

f(x) = 4 / (1 + x ^ 2)

function composite_midpoint(func, a, b, n)
    h = (b - a) / n
    return h * sum([func(a + (i + 1/2) * h) for i in 0:n-1])
end

function composite_trapezoidal(func, a, b, n)
    h = (b - a) / n
    return h/2. * sum([func(a + i * h) + func(a + (i+1) * h) for i in 0:n-1])
end

function composite_simpson(func, a, b, n)
    h = (b - a) / n
    return h/6. * sum([func(a + i * h) + 4 * func(a + (i + 1/2) * h) + func(a + (i+1) * h) for i in 0:n-1])
end

@show h = 1/100
@show composite_midpoint(f, 0, 1, 100)
@show composite_trapezoidal(f, 0, 1, 100)
@show composite_simpson(f, 0, 1, 100)

n_arr = [200, 500, 1000, 2000, 5000, 10000]
h_arr = 1. ./ n_arr

scatter!(h_arr, [pi - composite_midpoint(f, 0, 1, n) for n in n_arr], label="midpoint")
scatter!(h_arr, [pi - composite_trapezoidal(f, 0, 1, n) for n in n_arr], label="trapezoidal")
scatter!(h_arr, [pi - composite_simpson(f, 0, 1, n) for n in n_arr], label="simpson")
xlabel!("h")
ylabel!("error")

png("hscaling")
