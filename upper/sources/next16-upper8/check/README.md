# 同作業内の別実装検算

入力の next15-upper/PROOF.md と next13/PROOF.md、および前回の independent_check.py を読んだ。同じ今回の代数担当コードを読む前に general.py を作り、物理座標 (S,X) の時間微分から F を導出した。W は掲載された式を転記せず、11×11 初期微分行列の逆から生成した。今回の algebra/ のコード・結果JSONは import していない。担当間では最高係数式と終結式・実根分類の速報を交換した。したがって同作業内の別実装検算であり、盲検や形式証明ではない。

- exception.py/json: g=1 の新例外枝を全係数で排除。
- general.py/json: QQ(a,b) で全係数消去。rank 7 の初期jet条件を持つ4次元高接触空間。
- classify_exact.py/json: 終結式を再計算。C3 の零枝は元二式の剰余から有限解なし。H12 の全4実根を Sturm で数え、各候補の必要正性を厳密有理区間演算で棄却。
- realroots.py/json: 探索記録のみ。C3 の根を一般 b 式に代入した正に見える数値は有限解ではない。元二式の残差を満たさず、厳密な C3 枝処理で排除される。これを正の実現候補と採用していない。

再現は repo root で、各スクリプトに `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout 1800 python3 -B` を前置する。exception.py → general.py → classify_exact.py の順。Python標準FractionとSymPyの有理数・多項式演算を用い、確定判定に浮動小数を用いない。成果JSONは決定的に再生成する。

## 計算中の修正

classify_exact の最初の版は幅10^-12の有理区間のHorner評価で、第3根の値の符号を確定できず assert に停止した。根区間を有理数のまま幅10^-35以下に細分して確定した。これは区間演算の過大評価の解消で、判定条件や根を変更していない。次の版は検査を通過したが SymPy Integer を JSON に書けず停止した。Sturm の符号値を Python int に変換して保存を完了した。上記の C3 の探索上の見かけ候補も含め、これらを成立証明に使っていない。

最終点検では全11単項式、h=0の4単項式、SL/Lの行列式、v/d0の係数−1/1、消去後の残る9係数の恒等零をassertした。係数−1/1の初回assertはSymPy式の構文一致を比較して停止したため、`cancel(差)==0`という有理式恒等一致へ直し、全検査が通過した。
