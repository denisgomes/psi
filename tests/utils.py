from math import fabs


def compare(a, sol, error_percent=3.0):
    """Compare the absolute value of two values and return True if both are
    within the error_percent bounds.
    """
    upper_bound = sol*(1 + error_percent/100)
    lower_bound = sol*(1 - error_percent/100)

    if sol < 0.:
        upper_bound, lower_bound = lower_bound, upper_bound

    return (a <= upper_bound and a >= lower_bound)
