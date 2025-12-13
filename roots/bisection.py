def bisect(func, lower_bound, upper_bound, tolerance, max_iterations):

    if func(lower_bound) * func(upper_bound) >= 0:
        print("Error: The root is not bracketed. f(lower) and f(upper) must have opposite signs.")
        return None

    iteration_count = 0
    root_estimate = 0.0
    current_error = 100.0

    print(f"{'Iter':<10} {'Root Estimate':<20} {'Current Error (%)':<20}")
    print("-" * 50)

    while iteration_count < max_iterations:
        previous_estimate = root_estimate
        root_estimate = (lower_bound + upper_bound) / 2
        iteration_count += 1

        if iteration_count > 1:
            if root_estimate != 0:
                current_error = abs((root_estimate - previous_estimate) / root_estimate) * 100
            else:
                current_error = abs(upper_bound - lower_bound)

        test_val = func(lower_bound) * func(root_estimate)

        if test_val < 0:
            upper_bound = root_estimate
        elif test_val > 0:
            lower_bound = root_estimate
        else:
            current_error = 0.0

        print(f"{iteration_count:<10} {root_estimate:<20.5f} {current_error:<20.5f}")

        if current_error < tolerance or iteration_count >= max_iterations:
            break

    return root_estimate