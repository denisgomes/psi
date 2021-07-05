from math import fabs


def compare(a, sol, error_percent=3.0, sigfig=4):
    """Compare the absolute value of two values and return True if both are
    within the error_percent bounds.

    The values are compared upto a predefined number of significant figures.
    """
    upper_bound = sol*(1 + error_percent/100)
    lower_bound = sol*(1 - error_percent/100)

    if sol < 0.:
        upper_bound, lower_bound = lower_bound, upper_bound

    return (round(a, sigfig) <= round(upper_bound, sigfig) and
            round(a, sigfig) >= round(lower_bound, sigfig))
