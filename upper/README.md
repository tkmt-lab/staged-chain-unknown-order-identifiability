# Complete upper-bound evidence for orders 8, 11 and 15

This supplement contains the adopted proofs, exact certificate inputs and checkers for the three universal upper bounds. All numeric inputs are local. A private repository, credentials and remote data downloads are unnecessary after extraction. `SOURCE-MAP.json` fixes every selected input and checker by SHA256 and lists path adaptations and the explicitly recorded public-metadata edits. Files with historical claims retain their original date and scope; the accepted integration record determines the current combined claim.

From the package root:

```sh
sha256sum -c MANIFEST.sha256
python3 -B upper/check.py --verify
python3 -B upper/check.py --group smoke --output build/upper-smoke
```

The integrity command checks every selected file. The smoke command runs seven bounded checks; it does not recertify all global exclusions. Each replay uses a fresh copy, refuses numeric inputs outside that copy, runs sequentially with one numerical thread, CPU 590 s / wall 595 s and address space at most 6 GiB per child (several original checkers lower this to 3 or 4 GiB). The parent plus child are two Python processes. A timeout is INCOMPLETE, never PASS. Do not use Python `-O`.

Python 3.12, SymPy 1.14.0, NumPy and SciPy are required. The accepted upper-eleven review used python-flint 0.8.0; upper-fifteen used 0.9.0. This packaging test uses 0.9.0. Install dependencies in an isolated environment; dependency binaries are not included. Example: `python3 -m pip install -r upper/requirements.txt`. Exact installed versions used here are recorded in the delivery test record. Numerical assertions use exact arithmetic, not floating-point tolerances.

`JOBS.json` is the ordered executable index. Use `--group u8`, `--group u11`, `--group u15` or `--job <id>` with a new `--output` directory. `--group u15-without-crt` reruns the eight other U15 jobs while inheriting the three saved raw CRT certificates. Single jobs can use the bundled intermediate certificates; a whole group regenerates them in dependency order. The final accepted certificates are preserved in the source copy. Timestamps and output hashes can change on replay.

|Bound|Complete adopted reasoning|Executable route|
|---|---|---|
|8|`sources/next16-upper8/PROOF.md`, sections 1–6: residual support, general branch, g=1, degree-drop cubic, four real roots, h=0 and ODE uniqueness|u8-general, u8-exception, u8-classify_exact; algebra/ contains the original parallel coefficient certificates|
|11|`sources/verify-b13-blind-2026-11-09/stage2/reduction/AUDIT.md`, `coverage/AUDIT.md`, `roots/AUDIT.md`; `sources/b13-close-2026-11-07/upper_contact10_PROOF.md`, `structure_G22_PROOF.md`, `upper_nine_RUR_PROOF.md`|u11 group: original-ODE reduction and determinants; 17 final-route elimination checks; bridge/new-prime checks; complete 28-root cover and positivity exclusion|
|15|`sources/verify-b15-2026-11-09/THEORY.md`, `VERDICTS.json` (P01–P08); `sources/b15-close-2026-11-09/PROOF.md`, `LIFTING.md`, local-cluster prose|u15 group: determinant, three full raw CRT checks, all fibers, both local valuations, both complete finite-field projection evaluations, boundary and negative controls|

The U11 route includes `sources/b13-lower3-2026-11-07/elim_general.json` and the U15 route includes all four finite-field projection polynomials in `sources/b14-pilot-2026-11-09/upper/104729/`. Those are indispensable numerical inputs despite their historical directory names. Obsolete exploration branches, prime caches and the degree-1954 global lower-bound classification are not dependencies of these adopted proofs.

Prose and code divide the obligations: all-rank contact reduction, ODE uniqueness, finite projection/closure, lifting from characteristic p, and precision-aware local factors must be read with the exact certificates. A file hash or stored PASS field is not a mathematical proof. The full U15 CRT/projection route may time out at 595 seconds on a slower machine; this packaging task does not assert completion of those jobs. Keep the saved certificates and log the unfinished job; splitting by primes/evaluation points is a possible future engineering step, not a completed check here.

The lower-witness proofs are included at `sources/next16-lower8/verify_ALGEBRAIC_PROOF.md`, `sources/verify-b13-blind-2026-11-09/lower_proof.md`, and `sources/verify-contact11-2026-11-10/PROOF.md`. Their portable checkers are in `../finite-jet/`. Historical relative links to exploration reports are provenance, not additional inputs; the complete adopted routes are the table above and `JOBS.json`.

Historical scripts/prose may mention larger time limits. The bounded entry above is authoritative for this delivery and never extends the 595-second per-job wall limit. The original CRT03 replay took about 599.64 seconds in its accepted run; a slower bounded replay can therefore be incomplete without invalidating that historical certificate.

See [historical references](HISTORICAL-REFERENCES.md) for old filenames, output labels and execution examples. [REPAIR-RECORD.json](../REPAIR-RECORD.json) records the source-to-public editing boundary.
