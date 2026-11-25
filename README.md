# Identifying infection stages, infectivity weights, and initial states from incidence

Version 1.1.0 — 2026-11-26 — DOI: to be assigned on publication.

For a single staged SIR chain with a common progression rate, a complete noiseless incidence germ and known population determine the unknown finite stage number, both rates, normalized infectivity weights, and all initial states.
This gives researchers studying staged infection models a uniqueness result even when the stage structure and infectivities are unknown.
The true model and every finite competing model have strictly positive rates, weights and initial states, w1=1, and use the same physical time.

## First run

With Python 3.10+ (standard library only), run from this directory:

```sh
python3 -B reproduce.py --output build/checks
```

Expected: JSON with `"status": "PASS"`, `"manifest_entries": 282`, and five `"baseline_match": true` results (`darboux`, `equilibrium`, `boundary`, `weighted`, `recovery`). Timing and peak memory vary. The same result is saved to `build/checks/summary.json`.

Read [the paper](paper.pdf), [source](paper.tex), [claims](claims.json), and [verification record](VERIFICATION-RECORD.md). The result covers every positive true point and all finite positive competing orders. Equal weights give the SInR corollary. The complete analytic germ means the incidence function on a neighborhood of one time, in the original physical time.

## Reproduce

Python 3.10+ on Linux/macOS, standard library only. The checks run offline from this directory with assertions enabled.

```sh
python3 -B reproduce.py --output build/checks
python3 -B check_packet.py
python3 -B publication_privacy.py
```

The package is self-contained; the checks use only Python and the files in this directory. Five checks run serially with 60-second CPU and wall limits per child. Memory use is reported; the resource limits cover CPU and wall time. They compare scientific fields against frozen exact baselines, excluding timing and memory from equality. Outputs remain under ignored `build/`.

The PDF is included. Optional offline rebuild with pdfLaTeX (TeX Live with the `lmodern`, AMS, `booktabs`, `geometry` and `hyperref` packages installed). Run twice to resolve cross-references:

```sh
mkdir -p build/pdf
pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error -output-directory=build/pdf paper.tex
pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error -output-directory=build/pdf paper.tex
```

## Attribution

Pant, Saucedo, and Pogudin formulated the SInR model and fixed-order identifiability question in equation (3.3) and Conjecture 1 of their version-2 preprint, *Uncovering identifiability of epidemiological models: basic reproduction number and complementary data streams*, doi:10.64898/2026.01.16.26344284. Hurtado–Kirosingh provide linear-chain background; Ritschel–Wyller provide a comparison for mixed-Erlang kernel estimation. Full references are in the article. I am grateful to these authors for the model, question, and mathematical background on which this article builds.

The preceding same-author work [Global identifiability for the staged-infection SIR model](https://github.com/tkmt-lab/sinr-incidence-identifiability) concerns each fixed known order, generic true points, and arbitrary real competitors. Here the order and normalized infectivity weights are unknown and every positive true point is covered, with competitors restricted to the positive family. The two results address different competitor domains and true-point conditions.

Author: Tatsuhiko Takimoto. Research contact: research@tatsuhikotakimoto.com

## 日本語：発生数から感染の段数・感染力の重み・初期状態を識別する

版 1.1.0 — 2026-11-26 — DOI: 公開時に付く。

一列の段階感染モデルで遷移率が全段に共通の場合に、完全な無雑音の発生数の芽と既知の総人口から、未知の有限段数・二つの率・正規化した感染力の重み・全初期状態が一意に決まることを示します。
段階感染モデルを調べる研究者にとって、段数や感染力も未知の場合に、観測からモデルが一意に定まる条件を与える結果です。
真のモデルと任意の有限段数の競合モデルの双方で、率・重み・全初期状態が厳密に正、先頭の重みが1で、同じ物理時間を用いることが条件です。

最初の実行例（Python 3.10 以上、標準ライブラリのみ。一式のディレクトリで実行）：

```sh
python3 -B reproduce.py --output build/checks
```

期待する出力は JSON の `"status": "PASS"`、`"manifest_entries": 282` と、5 件（`darboux`、`equilibrium`、`boundary`、`weighted`、`recovery`）すべての `"baseline_match": true` です。時間・メモリ値は実行ごとに変わり、同じ結果が `build/checks/summary.json` に保存されます。続けて上の Reproduce 節の残りの二つの検査を実行してください。

発生数の芽とは、ある時刻の近くの解析関数全体です。結果は全ての正の真点と、任意の有限段数の正の競合を対象とし、全重み1の場合にはSInRの系が得られます。[論文](paper.pdf)に証明、[検証記録](VERIFICATION-RECORD.md)に有限検算とAI点検の範囲を記しています。PDFの再生成も上の Reproduce 節の手順で行えます。検査は Python と一式のファイルだけで完結し、CPU時間・実時間に上限を設け、メモリ使用量を記録します。

### 先行研究への感謝

Pant、Saucedo、Pogudin は上記の原著第2版の式 (3.3) と予想1でSInRのモデルと固定段数の識別の問いを示しました。Hurtado–Kirosingh は鎖表現の構成と一般化を、Ritschel–Wyller は比較対象となる混合Erlangカーネルの推定を扱っています。詳しい書誌は論文の参考文献にあります。モデルと問い、数学的な背景を与えてくださったこれらの研究の著者の方々に感謝します。

同著者の先行論文 [Global identifiability for the staged-infection SIR model](https://github.com/tkmt-lab/sinr-incidence-identifiability) は、既知の固定段数について、genericな真点と任意の実競合を扱います。今回は段数と正規化した重みも未知とし、全ての正の真点を対象に、競合を正のモデルに限ります。二つの結果では、真点と競合の条件が異なります。

著者：滝本竜彦。研究の連絡先：research@tatsuhikotakimoto.com

## Scope / 範囲

The theorem rests on a prose proof, supplemented by bounded exact checks and staged AI review. The finite checks do not certify the all-order analytic and valuation proofs. Human peer review, formal proof certification, and priority are not claimed; the same-author preceding article is not external independent evidence.

The unrestricted-order weighted-chain theorem uses a complete incidence germ and supplies no uniform finite derivative order for unbounded competitors. Section 7 supplies a restricted finite-jet implication: with true order one and competitor order fixed at 1, 2, 3 or 4, exact derivatives through order 4, 8, 11 or 15 suffice; order 15 covers their union. This gives no guarantee for finite samples, noisy estimation, unknown population or time scale, nonpositive competitors, unequal progression rates, branched chains, or epidemiological adequacy. Neither result subsumes the other's combination of true-point and competitor conditions.

定理の根拠は散文の証明で、有限の正確計算と段階的なAI点検がこれを補います。有限検算は全次数の解析的議論・付値論の機械認証ではなく、人間の査読・形式認証・優先権も認定していません。同著者の先行論文は独立の外部検証ではありません。

競合段数を制限しない重み付き鎖の定理は完全な発生数の芽を用い、無制限の競合に対する一律の有限微分階数を与えません。第7節の限定された系では、真の段数1・競合固定1、2、3、4に対して、それぞれ0〜4、8、11、15階の厳密微分値で十分です。競合1〜4の和集合には15階で十分です。有限個の測定、雑音下の推定精度、未知の総人口・時間尺度、非正の競合、段別の遷移率、分岐鎖、実際の疾病への適合性は対象外です。真点と競合の条件が異なるため、今回の結果と先行論文のいずれも他方を包含しません。

## License

Code and CI: Apache-2.0. Original article, text and generated data: CC BY 4.0. See LICENSE and the full license texts. The distribution consists of the original article, supporting files and license texts; cited third-party papers are linked in the article.

## Finite jets against one infectious compartment

Section 7 gives sharp output-germ orders 4, 8, 11, 15 for a one-compartment model compared with orders one through four in the broader positive compartment class. Its upper bounds restrict to common-rate weighted chains; together with the weighted identification theorem they identify a one-stage true model against competitors of at most four stages. Sharpness in the broader class does not establish sharpness in the chain subclass for orders two through four.

```sh
python3 -B finite-jet/reproduce.py --output build/finite-jet-checks
```

Expected output: `PASS`, five bounded certificate jobs. See [finite-jet reproduction](FINITE-JET-REPRODUCE.md) for the precise checked subset and the separate global upper-bound evidence. This command does not certify all four universal upper bounds.

第7節は一般の正の区画フィードバックで、一感染状態と1〜4状態を比較する出力germの最適階数4・8・11・15を示す。Wへの系は真の段数1・競合段数4以下に限る。上の5検査は小さな証明書の再現で、8・11・15の大域上界の再証明ではない。

The complete adopted upper-bound evidence is bundled in [upper/README.md](upper/README.md). Verify every delivered file with `sha256sum -c MANIFEST.sha256`. The accepted integration review is in [INTEGRATION-REVIEW.md](INTEGRATION-REVIEW.md).

## Additional fixed-order intervals

Version 1.1.0 adds 18 <= B_min(1,5) <= 42, 10 <= B_min(2,2) <= 75, and B_min(2,3) <= 278 (30 suffices on the common-initial-susceptibility branch). These are bounds, not optimal values. The older lower bound 10 for (2,3) has not been rechecked and is not adopted here as verified. Verification stages and limitations for every row, proof texts, certificates and the bounded replay entry are in [extensions/README.md](extensions/README.md).

追加は (1,5) の18〜42、(2,2) の10〜75、(2,3) の上界278と共通初期感受性の枝30。最適値ではない。(2,3) の旧下界10は再点検しておらず、成立した下界として採用しない。各行の検証段・独立性の留保と証拠は上記入口にある。

```sh
python3 -B reproduce.py --extensions --output build/extended-checks
```

The `studies/…` paths in the evidence are historical provenance names from the research stage, not paths inside this distribution; the historical note in `extensions/u42/check42.py` likewise records provenance, not an execution dependency.

証拠の中の `studies/…` は研究時の来歴の名前で、公開物の中の道ではない。`extensions/u42/check42.py` の歴史的注記も来歴であり、実行時の依存先ではない。
