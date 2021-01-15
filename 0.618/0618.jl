function search0618(func, a0, b0, eps, tau=0.618)
    a = a0
    b = b0
    i = 0
    while b - a > eps
        al = a + (1. - tau) * (b - a)
        ar = a + tau * (b - a)
        if func(al) < func(ar)
            b = ar
        else
            a = al
        end
        i += 1
        @show i
        @show a
        @show b
    end

    return (a + b) / 2.
end

function func1(alpha)
    return 1 - alpha * exp(-alpha^2)
end

print(search0618(func1, 0, 1, 0.01))
