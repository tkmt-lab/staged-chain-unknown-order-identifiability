"""Standalone exact polynomial checks; no producer code or data imported."""
from fractions import Fraction as F
from math import factorial
import json
import resource
import time


def add(*ps):
    d = {}
    for p in ps:
        for k, v in p.items():
            d[k] = d.get(k, F(0)) + v
    return {k: v for k, v in d.items() if v}


def scale(p, a):
    return {k: v*a for k, v in p.items() if v*a}


def mul(p, q):
    d = {}
    for k, a in p.items():
        for l, b in q.items():
            e = tuple(x+y for x, y in zip(k, l))
            d[e] = d.get(e, F(0)) + a*b
    return {k: v for k, v in d.items() if v}


def diff(p, i):
    d = {}
    for k, v in p.items():
        if k[i]:
            l = list(k)
            l[i] -= 1
            d[tuple(l)] = v*k[i]
    return d


def apply(p, vector):
    return add(*(mul(diff(p, i), v) for i, v in enumerate(vector)))


def determinant(a):
    a = [[F(v) for v in row] for row in a]
    d = F(1)
    for j in range(len(a)):
        i = next((i for i in range(j, len(a)) if a[i][j]), None)
        if i is None:
            return F(0)
        if i != j:
            a[i], a[j] = a[j], a[i]
            d = -d
        pivot = a[j][j]
        d *= pivot
        for i in range(j+1, len(a)):
            ratio = a[i][j]/pivot
            for k in range(j+1, len(a)):
                a[i][k] -= ratio*a[j][k]
    return d



def main():
    begin = time.perf_counter()
    results = []
    for n in range(1, 7):
        weights = [F(1)] + [F(j+1, j) for j in range(2, n+1)]
        total = sum(weights)
        zeta = F(3, 2)
        a = F(1)/(2*total)
        # Physical observability with the non-unit output row.
        amat = [[zeta*(int(i == j+1)-int(i == j)) for j in range(n)] for i in range(n)]
        row = weights[:]
        obs = []
        for r in range(n):
            obs.append(row)
            row = [sum(row[k]*amat[k][j] for k in range(n)) for j in range(n)]
        det = determinant(obs)
        expected = (-1)**(n*(n-1)//2)*zeta**(n*(n-1)//2)*weights[-1]**n
        assert det == expected and det
        # Independently evaluate determinant and explicit characteristic polynomial.
        char_checks = 0
        for lam in range(n+1):
            mat = [[F(lam+1)*int(i == j)-int(i == j+1)
                    -a*int(i == 0)*weights[j] for j in range(n)] for i in range(n)]
            rhs = F(lam+1)**n-a*sum(weights[j]*F(lam+1)**(n-1-j) for j in range(n))
            assert determinant(mat) == rhs
            char_checks += 1
        # Q evaluated in the shifted polynomial coordinate L=X+zeta.
        q_shift = [weights[n-1-r]*zeta**(n-1-r) for r in range(n)]
        recovered = [q_shift[n-j]/zeta**(j-1) for j in range(1, n+1)]
        assert recovered == weights
        assert q_shift[0] == weights[-1]*zeta**(n-1) and q_shift[-1] == 1
        tail = [sum(weights[i:]) for i in range(n)]
        assert tail[0] == total
        assert all(-tail[i]+(tail[i+1] if i+1<n else 0) == -weights[i] for i in range(n))
        residuals = 0
        if n >= 2:
            one = {(0,)*(n+1): F(1)}
            vs = []
            for i in range(n+1):
                e = [0]*(n+1); e[i] = 1
                vs.append({tuple(e): F(1)})
            x,s,t = vs[:3]
            h = [one, {}]+vs[3:]
            tp = [one]
            for i in range(1,n): tp.append(mul(tp[-1],t))
            q = [add(*(scale(mul(h[r-j],tp[j]),F(1,factorial(j))) for j in range(r+1))) for r in range(n)]
            p = add(*(scale(q[i],weights[i]) for i in range(n)))
            xp = mul(x,p)
            vector = [scale(mul(xp,s),-1),add(scale(s,-1),mul(xp,s)),add(one,scale(mul(xp,t),-1))]
            for r in range(2,n):
                ar = add(scale(tp[r],F((-1)**r,factorial(r))),scale(h[r],-1),mul(t,h[r-1]))
                vector.append(mul(xp,ar))
            u = [mul(s,qi) for qi in q]
            jw = add(*(scale(u[i],weights[i]) for i in range(n)))
            assert not add(apply(x,vector),mul(x,jw))
            for i in range(n):
                target=add(mul(x,jw) if i==0 else u[i-1],scale(u[i],-1))
                assert not add(apply(u[i],vector),scale(target,-1))
            residuals = n+1
        results.append({'n':n,'weights':list(map(str,weights)),
                        'characteristic_evaluation_checks':char_checks,
                        'moving_coordinate_identities':residuals,
                        'observability_determinant':str(det),
                        'recovered_weights_match':True,'EL_linear_identity':True})
    # Pole numerator recurrence q^k+1 D^k(y/q), Dq=-y.
    one = {(0,0):F(1)}
    q={(1,0):F(1)}; y={(0,1):F(1)}
    p=y
    residues=[]
    # Dy is set to q+y^2; remainder recursion must be independent of this choice.
    vector=[scale(y,-1),add(q,mul(y,y))]
    for k in range(7):
        rem={e:v for e,v in p.items() if e[0]==0}
        assert rem=={(0,k+1):F(factorial(k))}
        residues.append(k)
        p=add(mul(q,apply(p,vector)),scale(mul(y,p),k+1))
    print(json.dumps({'cases':results,'pole_remainder_checks_k':residues,
                      'elapsed_seconds':time.perf_counter()-begin,
                      'ru_maxrss_platform_units':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                      'scope':'finite exact checks; arbitrary-n claims use FIRST-STAGE proof'},indent=2))


if __name__ == '__main__':
    main()
