# GTEx Step10–12 ERV Figure Guide（日本語メモ）

このメモは、GTEx側のERV解析で作成したFigureが何を示しているかを、後で見返して理解できるように整理したものです。

## 前提

GTEx側のWGCNAは **gene-only** です。

```text
GTEx gene expression
→ WGCNA
→ AGE-associated candidate modules
→ そのmodule内のgene座標にERV annotationを後付け
```

したがって、GTExでは「ERVがWGCNA moduleに直接入っている」のではなく、

- module geneのpromoter付近にERVがあるか
- gene TSSからnearest ERVまでどれくらい近いか
- 特定ERV familyが特定moduleで多いか／少ないか

を調べています。

---

# Step10 promoter ERV overlap

## Figure: Step10_promoter_ERV_percent

### 何を見ているか

各AGE-associated candidate moduleについて、**promoter領域にERV fragmentが重なるgeneの割合**を示しています。

横軸：promoter ERV overlapを持つgeneの割合（%）  
縦軸：WGCNA module

破線は、解析対象となった全WGCNA geneでの基準値です。

```text
All eligible WGCNA genes
promoter ERV overlap = 23.86%
```

### 読み方

- 破線より右：全体平均よりpromoter ERV overlapが多い
- 破線より左：全体平均よりpromoter ERV overlapが少ない

おおよその傾向：

```text
greenyellow  約24–25%
blue         約24%
turquoise    約23%
red          約20–21%
midnightblue 約19–20%
darkgreen    約19%
pink         約17%
```

これは **module間でpromoter ERVを持つgeneの割合が違うか** を可視化したFigureです。

---

# Step10 promoter ERV enrichment / depletion

## Figure: Step10_promoter_ERV_forest

### 何を見ているか

各moduleでpromoter ERV overlapが、全体背景に比べて多いか少ないかを **odds ratio (OR)** で評価しています。

横軸：Odds ratio（log scale）  
縦軸：WGCNA module

```text
OR = 1
→ 全体背景と同程度

OR > 1
→ promoter ERV overlapが多い方向

OR < 1
→ promoter ERV overlapが少ない方向
```

点：odds ratio  
横線：95% confidence interval

### 読み方

95% CIがOR = 1をまたぐ場合、明確なenrichment / depletionとは言いにくいです。

このFigureは、単純な「割合」だけでなく、**背景に対する相対的な偏り**を確認するためのものです。

---

# Step10 ERV family-specific promoter overlap

## Figure: Step10_ERV_family_heatmap

### 何を見ているか

promoter ERV overlapをERV familyごとに分けて、各moduleで多いか少ないかを比較しています。

横軸：ERV family

```text
ERV1
ERVK
ERVL
ERVL-MaLR
```

縦軸：WGCNA module

色：log2 odds ratio

```text
log2 OR > 0
→ そのERV familyがそのmoduleで多い方向

log2 OR < 0
→ 少ない方向

log2 OR = 0
→ 背景と同程度
```

`*` は FDR < 0.05 を示します。

### 読み方

このFigureは、

> promoter ERV全体が多いか

だけではなく、

> ERV1 / ERVK / ERVL / ERVL-MaLR のどのfamilyが特定moduleに偏っているか

を見るためのものです。

GSE164471側ではERV subfamily expressionを直接定量しているため、後で統合するときに重要なFigureです。

---

# Step11 / Step12 nearest ERV distance

## Figure: Step10_nearest_ERV_distance_boxplot

### 何を見ているか

各moduleのgeneについて、**annotated TSSから最も近いERV fragmentまでの距離**を比較しています。

横軸：

```text
log10(nearest ERV distance + 1 bp)
```

縦軸：WGCNA module

破線：全eligible genesのmedian

```text
median nearest ERV distance = 3,803 bp
```

### 読み方

- 左側：ERVがTSSの近くにあるgeneが多い
- 右側：ERVまで遠いgeneが多い

log10表示なので、軸は通常のbp距離ではありません。

例：

```text
log10(distance + 1) = 2
→ 約100 bp

= 3
→ 約1,000 bp

= 4
→ 約10,000 bp
```

このFigureは、**module geneがERVの近くに位置する傾向があるか**を見るものです。

---

# 4つのFigureの違い

```text
promoter_ERV_percent
→ promoter ERVを持つgeneの割合

promoter_ERV_forest
→ その割合が背景より多い／少ないか（OR + 95% CI）

ERV_family_heatmap
→ どのERV familyがどのmoduleで偏っているか

nearest_ERV_distance_boxplot
→ gene TSSからnearest ERVまでの距離
```

---

# GTEx側でERVについて実際に言えること

GTEx WGCNAにはERV expressionそのものは入っていません。

解析したのは、AGE-associated candidate moduleに属するgeneについての **genomic ERV context** です。

```text
AGE-associated candidate modules
→ promoter ERV overlap
→ ERV-family-specific overlap
→ nearest ERV distance
```

この情報はGSE164471側の

```text
TEcount ERV subfamily expression
Telescope ERV locus expression
ERV–gene co-expression
```

と後で統合するために使います。

特に、

```text
GTExでERV近傍にあるAGE-related gene
＋
GSE164471で同じERV family/subfamilyと共発現するgene
＋
Telescopeで近傍locusが確認できる
```

という形で一致すれば、統合解析で優先候補として扱いやすくなります。
