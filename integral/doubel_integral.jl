function double_simpson(func, a, b, c, d, n)
    h = (b - a) / n
    k = (d - c) / n

    return h * k / 36 * sum([(
        16 * func(a + (i + 0.5) * h, c + (j + 0.5) * k)
        + 4 * (
            func(a + i * h, c + (j + 0.5) * k) + func(a + (i + 1) * h, c + (j + 0.5) * k)
            + func(a + (i + 0.5) * h, c + j * k) + func(a + (i + 0.5) * h, c + (j + 1) * k))
        + func(a + i * h, c + j * k) + func(a + (i + 1) * h, c + j * k)
        + func(a + i * h, c + (j + 1) * k) + func(a + (i + 1) * h, c + (j + 1) * k))
        for i in 0:n-1 for j in 0:n-1])
end

f(x, y) = exp(-x * y)

# for the second problem, change the integral into polar coordinates
g(r, θ) = exp(-r^2 * cos(θ) * sin(θ)) * r

@show double_simpson(f, 0., 1., 0., 1., 1000)
@show double_simpson(g, 0., 1., 0., pi/4., 1000)
