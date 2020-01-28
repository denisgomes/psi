
def compare(a, sol, percent=3):
    upper_bound = sol*(1 + percent/100)
    lower_bound = sol*(1 - percent/100)

    return a <= upper_bound and a >= lower_bound
