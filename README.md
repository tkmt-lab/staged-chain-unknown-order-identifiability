# Identifying infection stages, infectivity weights, and initial states from incidence

Version 1.0.0 — release date: 2026-11-16. Read [the paper](paper.pdf), [source](paper.tex), [claims](claims.json), and [verification record](VERIFICATION-RECORD.md).

## First run

With Python 3.10+ (standard library only), run from this directory:

```sh
python3 -B reproduce.py --output build/checks
```

Expected: JSON with `"status": "PASS"`, `"manifest_entries": 29`, and five `"baseline_match": true` results (`darboux`, `equilibrium`, `boundary`, `weighted`, `recovery`). Timing and peak memory vary. The same result is saved to `build/checks/summary.json`. No network or installation is needed for this check.

For a single staged SIR chain with positive rates and initial states, common progression rate, positive infectivity weights normalized by w1=1, and known total population N0, the complete noiseless incidence germ in the same physical time identifies the finite stage number, both rates, all weights, and all initial states against every finite positive competing order. No extra generic condition or upper order bound is imposed. Equal weights give the SInR corollary.

This is a prose theorem with bounded exact checks and staged AI review. It is not formally certified or human peer reviewed. It gives no guarantee for finite samples, finite derivative data, noisy estimation, unknown population or time scale, nonpositive competitors, unequal progression rates, branched chains, or epidemiological adequacy.

## Reproduce

Python 3.10+ on Linux/macOS, standard library only; no network, installation, sibling repository, or Git history is needed. Keep assertions enabled.

```sh
python3 -B reproduce.py --output build/checks
python3 -B check_packet.py
python3 -B publication_privacy.py
```

Five checks run serially with 60-second CPU and wall limits per child. Memory is measured, not capped. They compare scientific fields against frozen exact baselines, excluding timing and memory from equality. These finite checks do not certify the all-order analytic and valuation proofs. Outputs remain under ignored build/.

The PDF is included. Optional offline rebuild with pdfLaTeX (TeX Live with the `lmodern`, AMS, `booktabs`, `geometry` and `hyperref` packages installed). Run twice to resolve cross-references:

```sh
mkdir -p build/pdf
pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error -output-directory=build/pdf paper.tex
pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error -output-directory=build/pdf paper.tex
```

## Attribution

Pant, Saucedo, and Pogudin formulated the SInR model and fixed-order identifiability question in equation (3.3) and Conjecture 1 of their version-2 preprint, *Uncovering identifiability of epidemiological models: basic reproduction number and complementary data streams*, doi:10.64898/2026.01.16.26344284. Their question and model are credited in the paper.

The preceding same-author work [Global identifiability for the staged-infection SIR model](https://github.com/tkmt-lab/sinr-incidence-identifiability) (R042) concerns each fixed known order, generic true points, and arbitrary real competitors. Here the order and normalized infectivity weights are unknown and every positive true point is covered, but competitors are restricted to the positive family. Neither competitor domain subsumes the other claim. The companion is not external independent evidence.

Hurtado–Kirosingh provide linear-chain background; Ritschel–Wyller provide a comparison for mixed-Erlang kernel estimation. Full references are in the article. Priority is not asserted.

## 日本語

版 1.0.0 — 公開日：2026-11-16。

最初の実行例（Python 3.10 以上、標準ライブラリのみ。一式のディレクトリで実行）：

```sh
python3 -B reproduce.py --output build/checks
```

期待する出力は JSON の `"status": "PASS"`、`"manifest_entries": 29` と、5 件（`darboux`、`equilibrium`、`boundary`、`weighted`、`recovery`）すべての `"baseline_match": true` です。時間・メモリ値は実行ごとに変わり、同じ結果が `build/checks/summary.json` に保存されます。この計算に通信や追加インストールは不要です。続けて上の Reproduce 節の残りの二つの検査を実行してください。

正の一列の段階感染モデルについて、総人口が既知で、遷移率が全段で共通、感染力の重みが正で先頭の重みが1の場合を扱います。同じ物理時間の完全な無雑音の罹患率関数から、段数・二つの率・全重み・全初期状態が、任意の有限段数の正の競合の中で一意に決まります。全重み1ではSInRの結果に戻ります。

有限個の測定や有限階微分、雑音下の推定精度、未知の総人口・時間尺度、非正の競合、段別の遷移率や分岐鎖は対象外です。一般の証明は論文の文章にあり、付属の小規模な有理数計算はその機械認証ではありません。人間の査読・形式認証・優先権の認定はありません。

問いとモデルはPant–Saucedo–Pogudinに帰属します。同著者の先行R042は既知の固定段数・genericな真点と任意の実競合を扱い、今回とは競合の範囲が違います。先行論文は独立の外部検証ではありません。再現手順は上記と同じです。

## License

Code and CI: Apache-2.0. Original article, text and generated data: CC BY 4.0. See LICENSE and the full license texts. Third-party papers are not redistributed.
