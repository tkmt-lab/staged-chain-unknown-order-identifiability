> Public edition: historical references and earlier execution limits are mapped in [the local reference guide](../../HISTORICAL-REFERENCES.md); use the current bounded replay entry there. The mathematical text below is retained.

# Independent blind lower bound: B_min(1,3) >= 11

Input: `studies/R044/blind/INPUT.json`, SHA256
`86dca53623397e754b047ffa8dc8a139316986f08d1092cc826433f8c9510a86`,
read at 2026-11-08 20:25:49 UTC. No prior research, notes, mappings,
reports, code, or proposed witness was consulted. After deriving the
contact equations and finding the candidate independently, this agent
used the root agent's newly written `certify_interval.py` and its newly
written `lower_root_certificate.json` to complete rigorous certification.

The result is a rigorous existence proof with an exactly specified
algebraic witness, not a claim that rounded decimal parameters give exact
jet equality. The 384-bit dyadic arithmetic uses integer endpoint
operations with outward rounding. Floating point values printed below
are only illustrations.

## Exact algebraic specification

Put T=1/25 and d=-59/200. Define eight real numbers
z=(i,g,p2,p1,p0,nb,nq1,nq0) as the unique zero in the rational box
`root_box_exact` in `lower_root_certificate.json` of the eight polynomials
F_0,...,F_7 below. This box is centered within one dyadic unit of the
112-digit decimal `root_center` in `lower_candidate.json`, with each
radius equal to 10^-35 rounded upward to the 384-bit dyadic grid.

The polynomial construction uses formal Taylor coefficients in the
auxiliary common variable tau=t/T:

    S[0]=1, I[0]=i,
    Y[k] = sum_{j=0}^k S[j] I[k-j],
    S[k+1] = -T Y[k]/(k+1),
    I[k+1] = T (Y[k]-g I[k])/(k+1),
    V[k] = (Y[k] - sum_{j=1}^k S[j] V[k-j])/(141/200).

Compute through index 11. For k=0,...,8 define

    F_k = (k+3)(k+2)(k+1) V[k+3]
        + p2 (k+2)(k+1) V[k+2] + p1 (k+1) V[k+1] + p0 V[k]
        + nb (k+2)(k+1) Y[k+2] + nq1 (k+1) Y[k+1] + nq0 Y[k].

These are polynomials with rational coefficients. The exact integer
interval contraction certificate proves existence and uniqueness of z:
the fixed-point image lies strictly inside the box, and the infinity norm
of its derivative is at most 1.104e-24. It also proves F_8<0, approximately
-4.617899619429018e-6. Contraction is applied to z-Y F(z), with a fixed
rational matrix Y. The bound ||Id-Y J||<1 implies Y is nonsingular, so
its fixed point is a zero of F.

For orientation, i is approximately 12.60350355547678 and g approximately
10.90581427513914. All subsequent model parameters are exact functions of
this uniquely defined z and the roots specified next.

Set P(x)=x^3+p2 x^2+p1 x+p0, q1=nq1/nb, q0=nq0/nb, Q(x)=x^2+q1 x+q0,
and beta_R=-nb/T. Define a<b<c as the three positive roots of
r^3-p2 r^2+p1 r-p0. Their exact rational isolating brackets are
`rates_fixed_T_exact` in `lower_realization_certificate.json`. Exact
interval signs at the endpoints prove that there are three distinct
positive roots, so P is Hurwitz.

Put u=a/2, v=b/2 and define the right physical-time parameters by

    K_R = (1/T) [[-a, 0, 0], [u, -b, 0], [0, v, -c]],
    C_R = [1, (q1-b-c)/u, (c^2-q1*c+q0)/(u*v)],
    S_R(0)=141/200.

The transfer function is exactly Q(x)/P(x) in the scaled time coordinate:

    C_R (x Id - T K_R)^(-1) e1
      = 1/(x+a) + (q1-b-c)/((x+a)(x+b))
        + (c^2-q1*c+q0)/((x+a)(x+b)(x+c))
      = Q(x)/P(x).

The left parameters are beta_L=1, K_L=[-g], C_L=[1], S_L(0)=1,
X_L(0)=i. Let y be its analytic output. Define h=S_L+d and
w=y/(beta_R h). Let w_j=w^(j)(0). The right initial state is the unique
solution of

    C_R X0 = w_0,
    C_R K_R X0 = w_1 - y(0),
    C_R K_R^2 X0 = w_2 - y'(0) - (C_R K_R e1) y(0).

This specification is exact. The certificate calculates these quantities
without numerical differentiation: y0=i, y1=i(1-g-i),
y2=y1(1-g-i)+i[-i-(1-g)i], v0=y0/h0,
v1=(y1+i*v0)/h0, v2=(y2+2*i*v1+y1*v0)/h0, and w_j=v_j/beta_R.

## Admissibility and minimality

`lower_realization_certificate.py` imports only the newly written exact
interval arithmetic, the candidate, and the root certificate. Its output
`lower_realization_certificate.json` contains exact rational enclosures
for every claimed strict inequality and each nonzero elimination pivot.

For orientation the physical model is approximately

    beta_R = 1.800765530240113,
    K_R = [[-0.094075196609234, 0, 0],
           [ 0.047037598304617, -7.298373017203031, 0],
           [0, 3.649186508601516, -37.179675025823259]],
    C_R = [1, 803.377762420834, 283.659088082825],
    X_R(0) = [0.233986375164218, 0.0117006702169091, 0.00103494691734808].

The exact enclosures prove beta_R>0, all entries of C_R and X_R(0)>0,
and i,g>0. Both C_1 values are identically 1. The chain matrix is Metzler
and has strictly negative eigenvalues. Its column sums are
-a/(2T), -b/(2T), -c/T, all strictly negative. Its controllability
determinant is u^2 v/T^3>0; the three interval Gaussian elimination
pivots of its observability matrix exclude zero. Thus the right model
is controllable and observable. The one-dimensional left model is also
controllable and observable.

For a baseline common population N0=100 choose recovery populations by
conservation. The certificates give R_L(0)>0 (approximately 86.39649644)
and R_R(0)>0 (approximately 99.04827801). All initial populations are
strictly positive. Polynomial ODEs give analytic solutions; strict
positivity persists in a sufficiently small neighborhood of zero.
Defining R'= -1^T KX gives the required conservation equation.

## Why precisely the first eleven jet entries agree

Let D_tau=T D_t. Along the left solution, V=y/h and Y=y. The polynomial
system says the residual

    r = P(D_tau) V - beta_R*T Q(D_tau) y

has its derivatives in tau of orders 0,...,7 equal to zero, while
r^(8)(0)=8! F_8 is strictly negative. In physical time let
P_t(x)=T^-3 P(Tx), Q_t(x)=T^-2 Q(Tx). Then

    r_phys = P_t(D_t) V - beta_R Q_t(D_t)y = T^-3 r,
    r_phys^(8)(0)=T^-11 8! F_8.

Drive the right linear system by the target left output: let
Xbar'=K_R Xbar+e1*y, Xbar(0)=X0. Put e=C_R Xbar-w. The transfer identity
and the chosen initial state give

    P_t(D_t)e = -r_phys/beta_R,
    e(0)=e'(0)=e''(0)=0.

It follows successively that e^(j)(0)=0 for j=0,...,10, and

    e^(11)(0) = -T^-11 8! F_8/beta_R != 0.

The forced path (h,Xbar) has output beta_R h C_R Xbar=y+beta_R h e.
Its state-equation defect relative to the actual right nonlinear system
is O(t^11). Analytic ODE comparison, or coefficient recursion, gives a
state difference O(t^12). Consequently the actual right output satisfies

    y_R-y_L = beta_R h e + O(t^12),
    y_R^(k)(0)=y_L^(k)(0) for k=0,...,10,
    y_R^(11)(0)-y_L^(11)(0) = -(141/200) T^-11 8! F_8 > 0.

The exact interval certificate encloses this last difference near
312963884883025.8. There is no one-sided time shift or change of time:
the constructed K matrices and outputs use the same physical t.
The auxiliary tau merely supplies an exact common algebraic scale.

## Every prescribed common population

For arbitrary N>0 take lambda=N/100, multiply every initial compartment
of both baseline models by lambda, leave K and C unchanged, and divide
both beta values by lambda. The corresponding state trajectories and
outputs are multiplied by lambda. All conditions, minimality, exact
jet equalities through order 10, and the strict order-11 difference
are preserved. Therefore the claimed lower bound holds for every N>0.

## Reproduction

From the research repository root, in sequence:

    python studies/R044/verify-b13-blind-2026-11-09/lower_candidate.py
    python studies/R044/verify-b13-blind-2026-11-09/certify_interval.py \
      studies/R044/verify-b13-blind-2026-11-09/lower_candidate.json \
      studies/R044/verify-b13-blind-2026-11-09/lower_root_certificate.json
    python studies/R044/verify-b13-blind-2026-11-09/lower_realization_certificate.py

The candidate generator uses mpmath only to locate the box. The two
certification steps use exact rational/integer arithmetic for all proof
assertions. Each computation was sequential, single process/thread, and
completed well below 30 minutes and 3 GiB. The exploratory search scripts
are not required for the proof.
