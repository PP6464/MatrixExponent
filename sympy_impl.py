from sympy import *

def mat_exp(M):
    # Throw error if non-square
    m, n = shape(M)
    if m != n:
        raise RuntimeError("M must be square")
    # First check if M is diagonal
    is_diagonal = True
    for i in range(n):
        if not is_diagonal:
            break
        for j in range(n):
            if i != j and M[i, j] != 0:
                is_diagonal = False
                break
    if is_diagonal:
        # Element-wise exponentiate along diagonals
        res_mat = Matrix([[0 for _ in range(n)] for _ in range(n)])
        for i in range(n):
            res_mat[i, i] = exp(M[i, i])
        return res_mat
    # Check if matrix is diagonalizable
    try:
        change_basis, L = M.diagonalize()
        # Matrix is diagonalizable
        return change_basis * mat_exp(L) * change_basis ** -1
    except:
        # Matrix is not diagonalizable, so check for nilpotency
        current_powers = [eye(n)]
        k = 0
        is_nilpotent = False
        zero_mat = eye(n) - eye(n)
        while k < n:
            k += 1
            current_powers.append(M * current_powers[-1])
            if zero_mat == current_powers[-1]:
                # M is nilpotent with degree k
                is_nilpotent = True
                break
        if is_nilpotent:
            res_mat = zero_mat
            for i in range(k):
                res_mat += current_powers[i] / factorial(i)
            return res_mat
    return exp(M) # Just use default exponentiation otherwise