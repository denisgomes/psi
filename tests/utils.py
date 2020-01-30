from math import fabs


def compare(a, sol, error_percent=1.0):
    """Compare the absolute value of two values and return True if both are
    within the error_percent bounds.
    """
    upper_bound = fabs(sol)*(1 + error_percent/100)
    lower_bound = fabs(sol)*(1 - error_percent/100)

    return (fabs(a) <= upper_bound and fabs(a) >= lower_bound)
