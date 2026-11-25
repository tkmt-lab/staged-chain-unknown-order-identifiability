# Accepted finite-jet integration review

The seven-item integration review accepted the synthesis for every real N>0: B_N(1,1)=4, B_N(1,2)=8, B_N(1,3)=11, B_N(1,4)=15, and the restricted implication for weighted chains with true order one and competitor order one through four. It was a single-reviewer nonblind review of the research proof and code, inheriting the fixed component certificates. It was not an unseen independent derivation or a new complete global recomputation.

|Item|Accepted result and evidence boundary|
|---|---|
|1 Sharp list|Upper/lower components match 4/8/11/15. All lower witnesses and upper four were replayed; the global U11/U15 arithmetic was inherited. See upper/README.md for the complete adopted evidence.|
|2 Quantifiers|Both models vary all allowed coefficients/states; same physical time and known common real N>0. Upper bounds do not assume minimality; lower witnesses are minimal. Scaling multiplies all states and output by a=N/Nbase and beta by 1/a, leaving K,C,time fixed.|
|3 Contact lemma|Output equality through B implies residual order at least B+1-m; initial q and z through m-1 plus residual equations give the reverse feedback induction. The next difference occurs at m+ell+1. Exact equality is not inferred from interval zero containment.|
|4 Upper compression|U8 keeps h=0, g=1, all denominator and cubic degree-drop branches. U11 keeps all-rank reduction, G22 closure/intersection, full projection and unique fibers, all 28 roots. U15 retains P01–P08, raw CRT bound M>B+C, H80, precision-aware valuations, full projection and all fibers.|
|5 W connection|K=zeta(L-I), C=w^T satisfies the larger model assumptions. Finite jet implies germ; W then gives identification only for true order one, competitors 1–4. General lower witnesses at 2/3/4 do not prove weighted-chain sharpness.|
|6 Preservation|Old sections 1–6, old claim entries, original five computations/baselines and runners were retained. The integration patch tree and PDFs were checked in the accepted review. The evidence package examined at that stage retained its mathematical TeX and PDF byte-for-byte.|
|7 Literature|Existing fixed-version comparisons retained; general finite Taylor reduction, elimination and multiplicity bounds are established methods. No new priority or general-method novelty claim.|

Lower eleven was established from the dedicated first-stage statement; upper eleven was accepted after opening the research proof in the second stage. Exposure qualifications for four/eight, and nonblind verification for fifteen, remain. The old article's eleven-point review does not automatically cover the new section.

The accepted review identified one separate delivery deficiency: the complete upper-bound archives were absent from the integration patch. Version 1.1.0 supplies the adopted proof routes, all their numeric inputs and their checkers in upper/. Packaging completeness and replay results are documented separately; this review's historical execution scope is unchanged. The review covered the integration patch.

## Fixed identities

- REPORT.md: SHA256 `6c87f2d425a4bcf2864b326ec8ebc39e5774600e676956676760b48e3e776c24`.
- INPUT-HASHES.json: SHA256 `78a8be9a886dd685c74dabd1224ae17c9783222d125b50f1e69c36c2b55da455`.
- finite-jet-section.tex: SHA256 `58106eb7f1c48c8aea276f3f63d0035df353fde6c9bbaa538893852f9d8b4b2a`.
- finite-jet-integration.patch: SHA256 `c81c0c5bea606f5df59bf4a744a4bea4778fcaa38cc2d879d3bcd37459ab0c3a`.
- Acceptance receipt: NOTE at commit `aefd3e513ddefaee6c513e9cf65717b017ff61ae`, SHA256 `dc929bceaf8df15b7041a2dc92b64e2c4e61acceb42805589ff694aa0b103a7f`; see [receipt excerpt](INTEGRATION-ACCEPTANCE.md).

Packaging source snapshot: `84f6a857a90ad4e7cd103dc4006357da8b9f191a`. This is the evidence source snapshot, before the acceptance receipt; its NOTE SHA256 is `4f7750eec542f28e41601bc516e2186f36a1c040f83714180e032288b2c1dd39`. Acceptance is recorded separately at the receipt commit above; these hashes identify historical review/receipt documents and are not a request for readers to obtain a private repository. The mathematical evidence and executable inputs needed by the adopted routes are bundled locally. AI assistance was used for the synthesis, verification and this packaging; this record does not assert human peer review.

The accompanying `INTEGRATION-REVIEW-EXCERPT.md` preserves the original review's seven verdicts and mathematical discussion verbatim (Japanese). Its opening packaging-deficiency statement describes the earlier integration patch, before this supplement was prepared. Historical internal paths in that excerpt identify provenance; the local adopted evidence map is upper/README.md. The excerpt's hash is covered by MANIFEST.sha256.
