"""Post-comparison exact checks, using only the reviewer's local primitives."""
from fractions import Fraction as F
from math import comb
import json
import resource
import time


def matvec(a, v):
    return [sum(x*y for x, y in zip(row, v)) for row in a]


def rowmat(v, a):
    return [sum(v[i]*a[i][j] for i in range(len(v))) for j in range(len(v))]


def solve(a, b):
    aug = [list(row)+[v] for row, v in zip(a, b)]
    n = len(b)
    for k in range(n):
        p = next(i for i in range(k, n) if aug[i][k])
        aug[k], aug[p] = aug[p], aug[k]
        pivot = aug[k][k]
        aug[k] = [v/pivot for v in aug[k]]
        for i in range(n):
            if i != k:
                v = aug[i][k]
                aug[i] = [x-v*y for x, y in zip(aug[i], aug[k])]
    return [row[-1] for row in aug]


def main():
    begin = time.perf_counter()
    cases = []
    for n in range(1, 7):
        rate = F(3, 2)
        for weighted in (False, True):
            w = [F(1)]+[F(j+1, j) if weighted else F(1) for j in range(2, n+1)]
            a = [[rate*(int(i == j+1)-int(i == j)) for j in range(n)] for i in range(n)]
            initial = [F(j+1, 2*n*(n+1)) for j in range(n)]
            inputs = [F(k+2, k+3) for k in range(n)]
            # Forward jets from the state equation, independent of observation formula.
            jets = [initial]
            for k in range(n-1):
                step = matvec(a, jets[-1])
                step[0] += inputs[k]
                jets.append(step)
            obs = [w]
            for k in range(1, n):
                obs.append(rowmat(obs[-1], a))
            rhs = []
            for k in range(n):
                measured = sum(wj*ij for wj, ij in zip(w, jets[k]))
                known_input = sum(obs[k-1-h][0]*inputs[h] for h in range(k))
                rhs.append(measured-known_input)
            recovered = solve(obs, rhs)
            assert recovered == initial
            r = F(1)-F(1, 4)-sum(recovered)
            assert r > 0
            cases.append({'n': n, 'weighted': weighted, 'initial_state_recovered': True})
        # All-unit weights: compare coefficients in X, including constant term.
        for k in range(n):
            lhs = sum(rate**(j-1)*comb(n-j, k)*rate**(n-j-k)
                      for j in range(1, n-k+1))
            rhs = F(comb(n, k+1))*rate**(n-1-k)
            assert lhs == rhs
    print(json.dumps({'state_recovery_cases':cases, 'unit_weight_Q_coefficients_checked':21,
                      'elapsed_seconds':time.perf_counter()-begin,
                      'ru_maxrss_platform_units':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                      'limit':'synthetic input jets test the linear recovery identity, not epidemic simulation'},indent=2))


if __name__ == '__main__':
    main()
