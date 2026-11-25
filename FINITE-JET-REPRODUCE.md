# Finite-jet evidence and reproduction

For every real common population N>0, the positive compartment feedback class in section 7 has sharp output-germ orders 4, 8, 11, 15 for (1,m), m=1,2,3,4. Upper bounds include nonminimal realizations; lower witnesses are minimal. Exact derivatives use the same physical time. These are output-germ orders, not general recovery orders for arbitrary compartment realizations.

```sh
sha256sum -c MANIFEST.sha256
python3 -B finite-jet/reproduce.py --output build/finite-jet-checks
python3 -B upper/check.py --verify
```

The portable command uses Python 3.12's standard library and checks upper four and lower four/eight/eleven/fifteen in five sequential jobs. Its arithmetic, root contraction and positive-realization checks must be combined with the exact feedback induction and population scaling in section 7. Interval zero containment alone does not prove equality of derivatives. The next command checks the integrity of the complete adopted upper-bound supplement. Full mathematical replay commands, dependencies, proofs and resource limits are in [upper/README.md](upper/README.md).

The weighted-chain consequence uses germ identification theorem W after the finite-jet implication. It concerns a true one-stage chain and competitors of fixed order 1,2,3,4, or their union with order fifteen. Sharpness for the general compartment class does not imply weighted-chain sharpness at m=2,3,4. There is no uniform fifteen bound for unbounded competitor order or for arbitrary pairs of orders up to four.

The seven-item nonblind integration review was accepted; [INTEGRATION-REVIEW.md](INTEGRATION-REVIEW.md) records its precise scope and fixed source hash. Lower eleven's independent first stage and upper eleven's later review after opening the proof remain distinct. Reviews for four/eight carry exposure qualifications; fifteen was checked nonblind. The original article's five checks do not certify the addition. Packaging and smoke replay do not constitute a new full global certification, human peer review, formal proof or priority assessment. The degree-1954 global contact classification is not required for the local lower-fifteen witness.

General finite Taylor reduction, elimination and multiplicity methods are established work. The literature comparisons retain their original scope; this supplement adds no new priority search or claim of a new general method.
