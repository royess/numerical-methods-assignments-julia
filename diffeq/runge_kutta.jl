module RungeKutta
    function runge_kutta_solver(func, n, h, x0, y0)
        x = x0
        y = copy(y0)

        x_arr = zeros(n)
        y_arr = zeros(n, size(y)[1])

        for i in 0:n-1
            k1 = func(x, y)
            k2 = func(x + 0.5 * h, y + (0.5 * h) .* k1)
            k3 = func(x + 0.5 * h, y + (0.5 * h) .* k2)
            k4 = func(x + h, y + h .* k3)
            y += (h / 6.0) .* (k1 + 2 * k2 + 2 * k3 + k4)
            x += h
            x_arr[i + 1] = x
            y_arr[i + 1, :] = y
        end

        return x_arr, y_arr
    end
end