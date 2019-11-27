"""Various implementations of vectorized numerical solvers."""

import numpy as np


def gauss(a, b):
    """Gaussian elimination with partial row pivoting to reduce roundoff error.

    The a vector is the global stiffness matrix.

    The b vector is the global force matrix. It is mutated into the solution
    vector to save storage space.
    """
    n = len(b)
    tempa = np.zeros(n, dtype=np.float64)
    for k in range(n-1):
        # Determine where the largest pivot exists
        indexbig = np.argmax(np.abs(a[k:n, k]))
        # Move up the row with the largest pivot (entire row)
        tempa[:] = a[k, :]
        a[k, :] = a[k+indexbig, :]
        a[k+indexbig, :] = tempa
        tempb = b[k]
        b[k] = b[k+indexbig]
        b[k+indexbig] = tempb
        # Forward gaussian elimination
        for i in range(k+1, n):
            fac = a[i, k]/a[k, k]
            # a[i,k] = 0.0  # do not need to set explicitly
            a[i, k+1:n] = a[i, k+1:n] - fac*a[k, k+1:n]
            b[i] = b[i] - fac*b[k]

    for k in range(n-1, -1, -1):
        b[k] = (b[k] - np.dot(a[k, k+1:n], b[k+1:n])) / a[k, k]

    return b


if __name__ == "__main__":
    # Solution is [5, 1, -2]
    a = np.array([[1, 4, 1],
                  [1, 6, -1],
                  [2, -1, 2]])

    b = np.array([7, 13, 5])

    x = gauss(a, b)

    print(x)
