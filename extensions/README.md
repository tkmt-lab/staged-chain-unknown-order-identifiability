# Additional finite-jet bounds

This extension collects fixed-order results for the positive compartment feedback class. The output is y=beta S C X, with S'=-y, X'=KX+e1 y and R'=-1^T KX. Both models have beta>0, C>0, Ce1=1, Hurwitz Metzler K with nonpositive column sums, strictly positive initial S,X,R, a common known population N>0 and the same physical time. B_min is the least highest derivative order B for which equality of output derivatives 0 through B implies equality of analytic output germs. Upper bounds include nonminimal realizations. Lower witnesses are minimal on both sides.

|Pair / assertion|Bound|Verification stage and independence|Proof and executable evidence|
|---|---|---|---|
|(1,5), lower|18 <= B_min|Two stages: separately implemented first-stage witness certification from two fixed inputs, followed by a nonblind comparison. Prior result and cyclic-realization summaries were visible; not completely unseen blind review.|[input](l18/INPUT.json), [witness](l18/CANDIDATE.json), [proof](l18/STAGE1.md), `l18/independent.py`|
|(1,5), upper|B_min <= 42|Stage one left 42 undecided; acceptance was in stage two after reading the proof, with algebraic rechecks and a separate resultant implementation. Nonblind. The review recorded a possible brief third computation process; mathematical acceptance did not certify strict resource compliance.|[proof](u42/PROOF.md), `u42/certificate42.json`, `projection42.json`, `alternate42.json`, `check42.py`, `check_projection42.py`|
|(2,2), lower|10 <= B_min|Single nonblind verification; separate rational/interval implementation and perturbed-coefficient/state controls.|[proof](l22/PROOF.md), `l22/check.py`, `b22data/lower_PAIR10.json`, `lower_system_u.json`|
|(2,2), upper|B_min <= 75|Nonblind separate algebraic recalculation. The final order transfer and return to germs were checked by reading the argument, not recalculated.|[proof](u22/PROOF.md), [analytic bridge](u22/ANALYTIC.md), `u22/check.py`, `algebra.py`, `eliminate.py`, three `b22data/upper_degree8*.json` certificates|
|(2,3), upper|B_min <= 278|Nonblind verification of universal identities and support bounds. Irreducibility, invariant-factor branches and final analytic uniqueness depend on prose arguments read in review; the full universal resultant was not expanded.|[proof](u23/PROOF.md), `u23/checker.py` (generates its universal certificates without external data)|
|(2,3), common S0 branch|30 suffices|Same nonblind review and prose qualifications as the previous row. The branch requires equal initial susceptible states in the original common physical coordinates.|Same proof and checker|

These give the intervals 18 <= B_min(1,5) <= 42 and 10 <= B_min(2,2) <= 75, not optimal values. For (2,3), the older proposed lower bound 10 has not been rechecked here and is not adopted as an accepted lower bound. This extension establishes no parameter/state identification theorem for the broader compartment class and no uniform finite order for unbounded competitors. No new novelty claim is made.

## Reproduction

Install the versions in `requirements.txt` in the environment used for reproduction, then from the distribution root run:

```sh
python3 -B reproduce.py --extensions --output build/extended-checks
```

This preserves the five existing standard-library baseline checks and then runs nine sequential extension jobs. To run only the extension, from the distribution root:

```sh
python3 -B extensions/reproduce.py --output build/extensions
```

Expected: nine zero exit codes and final `status: PASS`. Each numerical job is one process and one thread, with 7 GiB address-space ceiling, 7000 CPU seconds and 7100 wall seconds; the small runner is the second process. Checkers may impose tighter limits. Outputs and writable certificate copies are confined to the chosen fresh output directory. Choose an output outside `extensions/`. The upper-75 universal expansion runs before interpolation and elimination; its output is consumed by those jobs. It can be the slowest step. These limits are ceilings, not claimed runtimes.

The existing `finite-jet/reproduce.py` and `upper/check.py` remain available for the previous results. Their large CRT regeneration jobs and inherited verification limits remain as documented in `../FINITE-JET-REPRODUCE.md` and `../upper/README.md`; the extension does not relabel inherited runs as new reproductions.

The proofs retain mathematical excerpts of the original derivations and verification. Historical program names in the upper-75 structural notes identify previous calculations; the current adopted complete route is `check.py bridge`, `check.py universal`, `algebra.py`, `eliminate.py`. It independently reconstructs the identities used by those notes. The corrected h=0 argument and germ bridge are in `u22/ANALYTIC.md`. The upper-278 checker directly recomputes all the universal identities needed by its bundled proof; separate historical generators are not dependencies.

`SOURCE-HASHES.json` connects fixed original hashes to distributed hashes. Portability changes alter paths, metadata or resource ceilings, not arithmetic. This packaging run is a single-author nonblind integration check, not a new independent verification-B acceptance, formal proof or human peer review.

## Public derived versions

The following five inputs/checkers are accepted as public derived versions, not byte-identical copies of their originals. The complete original-to-public SHA256 mapping and the permitted changes are recorded below and in `SOURCE-HASHES.json`. The five distributed files match the public SHA256 values in the table. Their numerical values and arithmetic are unchanged. The other thirteen inputs/checkers retain their original bytes.

|Public file|Original SHA256|Public SHA256|Change type|
|---|---|---|---|
|`l18/INPUT.json`|`9783b9fb21f5f3b67651b408c7809f0fa1b8307001ee23110a694351a3b292db`|`3e8c9e4d5247650fbece3681834925069a745653699b919d01ddb8fa8b7216c6`|Remove only the study metadata field; all remaining JSON values are unchanged.|
|`l18/independent.py`|`ff85a8ef542edf824a1ad517a69c3f99a152307e68e0ceb0b9997eb70ea002b4`|`3be86a82792b9af49755fccf5949dba0bfc13db899836a2e5f24791261794786`|Use the bundled input directory and its metadata-adjusted INPUT hash; reduce address-space limit from 8 to 6 GiB and CPU limit from 7200 to 7000 seconds. Arithmetic is unchanged.|
|`l22/check.py`|`92e56839d47a5576a9cc37245f104485c93dc4a3b4fab9ea1e5a06954f26451a`|`58a2aab4eb23b032bb192f14d44b1e7bbeff79f9da29c9587fa260a10dfea01a`|Replace the historical proof name with PROOF.md and use the bundled b22data input directory. Arithmetic is unchanged.|
|`u22/algebra.py`|`7d81590131e7a93c5596da21b361a541e5c8313cf71747382b8eb19be84d142b`|`51a11a3c587d0832215e868b10e39bb23a1ee1e156f02abb49260bb8661d98af`|Use the bundled b22data input directory. Arithmetic is unchanged.|
|`u22/eliminate.py`|`f749fd11adb90a4cc560622d59e0ed53388e71a51e73c88d1d296d94f20fdcd9`|`d39582376abc701e65e335de9f0fc662b42524d6e2193a17061f106cfb67d158`|Use the bundled b22data input directory. Arithmetic is unchanged.|
