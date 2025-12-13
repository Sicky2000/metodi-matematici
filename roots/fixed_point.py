def fixed_point(func, x0, tol, max_iter):
    xr = x0
    iter_count = 0
    ea = 100.0

    while iter_count < max_iter:
        xrold = xr
        xr = func(xrold)
        iter_count+=1

        if xr != 0:
            ea = abs((xr - xrold) / xr) * 100

        if ea < tol or iter_count >= max_iter:
            break

    return xr, iter_count
