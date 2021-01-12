module Euler

function euler_solver(func, h, n, x0, y0)
    x = x0
    y = y0

    x_arr = zeros(n)
    y_arr = zeros(n)

    for i in 0:n-1 
        y += h * func(x, y)
        x += h
        x_arr[i + 1] = x
        y_arr[i + 1] = y
    end

    return x_arr, y_arr
end

function improved_euler_solver(func, h, n, x0, y0)
    x = x0
    y = y0

    x_arr = zeros(n)
    y_arr = zeros(n)

    for i in 0:n-1
        k1 = func(x, y)
        k2 = func(x + h, y + h * k1)
        y += 0.5 * h * (k1 + k2)
        x += h
        x_arr[i + 1] = x
        y_arr[i + 1] = y
    end

    return x_arr, y_arr
end

end