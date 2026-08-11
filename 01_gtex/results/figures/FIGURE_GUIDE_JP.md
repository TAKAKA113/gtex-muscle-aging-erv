# GTEx Figure guide（日本語）

このファイルは、`01_gtex/results/figures/` にある主要Figureを、**「何を見ている図か → 軸の意味 → どう読むか → 何が言えるか」**の順で簡単に説明するためのメモです。

GTEx側では、**ERV発現そのものをWGCNAに入れていません。** まずgene-only WGCNAでAGE-associated candidate modulesとhub genesを作り、その後にRepeatMasker/RoadmapなどのERV・epigenomic annotationを重ねています。

---

## 1. `06_GenAge_enrichment_light.pdf`

### 何を見ている図か

WGCNAで得たmoduleの中に、既知の老化関連遺伝子である**GenAge genesが多く含まれているか**を見ています。

```text
WGCNA module genes
      ↓
GenAge gene listと照合
      ↓
moduleごとのoverlap
      ↓
Fisher's exact test / FDR
```

### 見方

- 各行または各点：WGCNA module
- enrichmentの指標：GenAge genesがそのmoduleにどれくらい偏っているか
- Odds ratioが使われている場合：
  - `OR > 1` = GenAge genesが多い方向
  - `OR = 1` = 背景と同程度
  - `OR < 1` = GenAge genesが少ない方向
- FDR / `*` がある場合：多重検定補正後に有意かを見る

### この図から言えること

この図は、**「AGEと相関したmoduleが、既知の老化遺伝子も多く含むか」**を確認する図です。

ただし、GenAge enrichmentがないmoduleでも、新しい老化関連networkである可能性は残ります。GenAgeは既知遺伝子との一致を見る外部validationです。

### 関連する数字

- WGCNA input: **20,000 genes**
- WGCNA内のGenAge genes: **191 genes**
- 7 target modulesから抽出したhub candidates: **210 genes**
- その中のGenAge genes: **3 genes**
  - `BSCL2`（red）
  - `S100B`（darkgreen）
  - `PML`（pink）

---

# Step10: module全体とERVのゲノム位置関係

Step10では、7つのAGE-associated candidate modulesについて、**moduleに属するgeneの近くにERVが多いか**を見ています。

重要なのは、ここで見ているERVは**expressionではなくgenomic annotation**だという点です。

---

## 2. `Step10_ERV_family_heatmap.pdf`

### 何を見ている図か

各WGCNA moduleについて、gene promoterに重なるERVをfamily別に分け、**どのERV familyが多い／少ないか**を見ています。

対象family：

```text
ERV1
ERVK
ERVL
ERVL-MaLR
```

### 軸

- 縦軸：WGCNA module
- 横軸：ERV family
- 色：`log2 Odds Ratio`
- `*`：FDR < 0.05

### 色の読み方

```text
log2 OR > 0
→ そのmoduleのgene promoterに、そのERV familyが多い

log2 OR = 0
→ backgroundと同程度

log2 OR < 0
→ そのERV familyが少ない
```

例えば `log2 OR = 1` なら、Odds Ratioでは約2倍です。

### この図から言えること

この図で見ているのは、

> **AGE-associated candidate moduleごとに、promoter付近に存在するERV familyの構成が違うか**

です。

`*` が付いているcellは、そのfamilyの偏りがFDR補正後も有意だったことを示します。

### 言えないこと

- ERVがそのgeneを発現制御している、とはまだ言えない
- そのERVがRNAとして発現している、とは言えない
- ERVが老化を引き起こしている、とは言えない

これは**genomic proximity / enrichment**の結果です。

---

## 3. `Step10_nearest_ERV_distance_boxplot.pdf`

### 何を見ている図か

各moduleのgeneについて、**TSS（transcription start site）から最も近いERV fragmentまでの距離**を比較しています。

### 軸

- 縦軸：WGCNA module
- 横軸：`log10(nearest ERV distance + 1 bp)`
- 破線：全eligible genesの中央値

全体のmedian：

```text
3,803 bp
```

### 横軸の意味

log10変換されているので、だいたい次のように読みます。

```text
0  → 0 bp付近
1  → 10 bp
2  → 100 bp
3  → 1,000 bp
4  → 10,000 bp
5  → 100,000 bp
```

したがって、**左にあるほどERVがgene TSSに近い**です。

### Boxplotの見方

- 太線：中央値
- 箱：中央50%の範囲
- whisker：さらに広い分布
- 点：外れ値

### この図から言えること

この図では、

> **AGE-associated candidate moduleのgeneが、全体と比べてERVの近くに配置されている傾向があるか**

を見ています。

moduleの中央値が破線より左なら「ERVまでの距離が全体より短い傾向」、右なら「遠い傾向」です。

ただし、この図だけでは統計的有意差は判断できません。対応する統計検定結果が必要です。

---

## 4. `Step10_promoter_ERV_forest.pdf`

### 何を見ている図か

各moduleのgene promoterにERVが重なる割合が、**全WGCNA backgroundと比べてenrichmentしているか／depletionしているか**をOdds Ratioで示した図です。

### 軸

- 縦軸：WGCNA module
- 横軸：Odds Ratio（log scale）
- 点：Odds Ratio
- 横棒：95% confidence interval
- `OR = 1`：backgroundと差がない基準

### 読み方

```text
OR > 1
→ promoter ERV overlapが多い

OR < 1
→ promoter ERV overlapが少ない

95% CIが1をまたぐ
→ 明確なenrichment/depletionとは言いにくい

95% CIが1より完全に右
→ enrichment方向

95% CIが1より完全に左
→ depletion方向
```

### この図から言えること

`Step10_promoter_ERV_percent.pdf`が単純な割合なのに対し、このforest plotは、

> **その割合の差がbackgroundに対してどのくらい大きいか、また不確実性がどれくらいあるか**

を示しています。

---

## 5. `Step10_promoter_ERV_percent.pdf`

### 何を見ている図か

各moduleについて、**gene promoterに1つ以上ERV fragmentが重なるgeneの割合**をそのまま%で表示しています。

### 軸

- 縦軸：WGCNA module
- 横軸：Genes with promoter ERV overlap (%)
- 破線：全eligible WGCNA genesの割合

全体のreference：

```text
23.86%
```

### 読み方

```text
bar > 23.86%
→ そのmoduleではpromoter ERVを持つgeneが全体より多い

bar < 23.86%
→ 全体より少ない
```

### Step10 forestとの違い

```text
promoter_ERV_percent
= 実際に何%なのかを見る

promoter_ERV_forest
= backgroundに対してどの程度enrichment/depletionかを見る
```

したがって、この2枚は同じ情報を別の角度から見ています。

---

# Step11: hub genesに絞ったERV解析

Step10は**module内の全gene**を見ました。

Step11では、そこからさらに絞った**210 hub-gene candidates**を対象にERVとの関係を見ています。

```text
7 AGE-associated candidate modules
        ↓
各module top30 hub genes
        ↓
210 hub candidates
        ↓
ERV genomic contextを確認
```

---

## 6. `Step11_hub_MM_GS_ERV_scatter.pdf`

### 何を見ている図か

hub candidate genesについて、

```text
MM = moduleの中心性
GS_AGE = ageとの関連
```

を同時に見ながら、**そのgeneがERV genomic annotationを持つか**を重ねたscatter plotです。

### MMとは

```text
MM = Module Membership
```

gene expressionとmodule eigengeneの相関です。

```text
|MM|が高い
→ module全体の発現patternによく一致
→ module内で中心的なgene候補
```

### GS_AGEとは

```text
GS_AGE = Gene Significance for AGE
```

gene expressionとAGEの相関です。

```text
GS_AGE > 0
→ 年齢とともに発現が上がる方向

GS_AGE < 0
→ 年齢とともに発現が下がる方向

|GS_AGE|が大きい
→ 年齢との関連が強い
```

### Scatterの読み方

基本的には、

```text
|MM|が大きい
＋
|GS_AGE|が大きい
```

geneほど、**module hubとして中心的で、かつAGEとの関連も強い候補**です。

そこにERV statusが色・形などで重ねられています。

### この図で重要なgene

この図の目的は、

> **network中心性が高くAGEとも関連するhub geneの中に、ERV genomic contextを持つgeneが存在するか**

を見ることです。

GSE164471との統合時には、ここでERV contextを持つhub geneが特に重要な候補になります。

---

## 7. `Step11_hub_promoter_ERV_percent.pdf`

### 何を見ている図か

Step10のpromoter ERV overlapを、**module全体ではなくhub genesだけに限定して**見ています。

### 読み方

各moduleのtop hub genesについて、

```text
promoterにERVあり
        /
hub genes
        ×100
```

の割合です。

### Step10との違い

```text
Step10_promoter_ERV_percent
= module内の全genes

Step11_hub_promoter_ERV_percent
= moduleのhub genesだけ
```

ここでhub genesのERV overlapが高ければ、

> **ERV近傍geneが単にmoduleに存在するだけでなく、module networkの中心部分にも存在する**

という情報になります。

ただし、これも因果関係を意味しません。

---

# Step12A: hub geneとERVの位置関係をさらに分解

Step12Aでは、hub genesに関連するERVが、**geneのどの位置にあるか**を詳しく分類しています。

---

## 8. `Step12A_hub_ERV_context_by_module.pdf`

### 何を見ている図か

hub genesのERV genomic contextを、**WGCNA moduleごとに比較**した図です。

つまり、

```text
red hub genes
pink hub genes
blue hub genes
...
```

について、ERVがどのgenomic contextに存在するかを比較しています。

### 見方

図のlegendに表示されているgenomic-context categoryごとに、各moduleのhub genesがどの程度含まれるかを見ます。

主に確認したいのは、

- promoter近傍が多いmodule
- gene body内のERVが多いmodule
- より遠位のERVが多いmodule

のようなmodule間の違いです。

### この図から言えること

この図は、

> **AGE-associated networkのhub genesとERVの位置関係が、moduleによって異なるか**

を見る図です。

後のGSE164471/Telescope統合では、特定moduleのhub gene近傍にあるERV locusを優先して照合できます。

---

## 9. `Step12A_hub_ERV_genic_context.pdf`

### 何を見ている図か

Step12Aの情報を、module別ではなく、**hub genes全体としてERVがどのgenic contextに存在するか**にまとめた図です。

### `context_by_module`との違い

```text
hub_ERV_context_by_module
= moduleごとの違いを見る

hub_ERV_genic_context
= hub全体で、ERVがどこに位置するかを見る
```

### 見方

図中のcategoryの比率を比較します。

たとえばpromoter / intragenic / otherなどのcategoryが表示されている場合、どのcontextがhub-associated ERVで最も多いかを確認します。

### この図から言えること

この図は、

> **AGE-related hub candidatesに関連しているERV annotationが、geneのどの領域に多いか**

を俯瞰するための図です。

これはERVがgene regulationに関わる可能性を考えるための**位置情報**であり、regulatory activityそのものを証明する図ではありません。

---

# Step12B: Roadmap epigenomeによる機能的annotation

Step10–12Aまでは主に「ERVがgeneのどこにあるか」という**genomic position**を見ていました。

Step12BではRoadmap Epigenomics annotationを使い、ERVが存在する領域が**どのようなchromatin stateを持つか**を確認しています。

このstepによって、単なる「近いERV」から、

```text
近いERV
        ↓
その領域がenhancerなどのregulatory chromatin stateか？
```

という情報を追加しています。

---

## 10. `Step12B_Roadmap_ERV_family_by_epigenome.pdf`

### 何を見ている図か

ERV familyごとに、Roadmap Epigenomicsでannotateされたepigenomic/chromatin stateとの関係を比較しています。

対象ERV familyはStep10と同様に主に、

```text
ERV1
ERVK
ERVL
ERVL-MaLR
```

です。

### 見方

縦軸・横軸はFigure上のlabelに従いますが、基本的には、

```text
ERV family
×
epigenomic state / epigenome
```

の組み合わせを比較する図です。

色や値がenrichmentを表している場合、値が大きい組み合わせほど、そのERV familyがそのchromatin contextに多く見られることを意味します。

### この図から言えること

この図の目的は、

> **ERV familyによって、存在するepigenomic contextが異なるか**

を見ることです。

GTEx側ではERV expressionを直接測定していないため、このRoadmap情報は、ERV locusに**regulatory potentialがありそうか**を見る補助annotationになります。

---

## 11. `Step12B_Roadmap_hub_enhanc...pdf`

### 何を見ている図か

ファイル名から分かる通り、Roadmap annotationを用いて、**hub genesに関連するERVとenhancer annotationの関係**を見ています。

このFigureは、Step12Bの中でも特に、

> **hub gene近傍のERVが、enhancerとしてannotateされたchromatin領域に重なるか**

を確認するためのものです。

### 見方

Figure上でmodule、hub gene、ERV family、epigenomeなどが比較されている場合は、enhancer overlapが多いgroupを確認します。

重要なのは、

```text
ERVがhub geneに近い
```

だけでなく、

```text
ERVがhub geneに近い
＋
Roadmap上でenhancer contextに重なる
```

という候補を分けられる点です。

### この図から言えること

このFigureでenhancer overlapを持つERV–hub gene pairは、**GSE164471/Telescopeで優先的に確認する候補**になります。

ただしRoadmap annotationは別sample/epigenome由来のreference annotationなので、GTEx sampleでそのenhancerが実際にactiveだったことを直接証明するものではありません。

---

# Figure全体のつながり

このFigure群は、次の順番で理解すると整理しやすいです。

```text
WGCNA
↓
AGE-associated candidate modules
↓
06 GenAge
既知のaging genesとの一致を確認
↓
Step10
module全体とERV genomic context
↓
Step11
hub genesだけに絞ってERV contextを確認
↓
Step12A
hub-associated ERVのgenic positionを細分化
↓
Step12B
Roadmap epigenomeでregulatory contextを追加
```

つまり、後半に行くほど、

```text
AGE-associated gene
        ↓
network hub gene
        ↓
ERVが近いhub gene
        ↓
特定genic contextにERVがあるhub gene
        ↓
Roadmap上でもregulatory contextを持つ候補
```

と候補を絞っています。

---

# GSE164471との統合で重要になる部分

GTExとGSE164471ではERV情報の種類が違います。

```text
GTEx
Gene expression WGCNA
+ ERV genomic position
+ Roadmap epigenomic annotation

GSE164471
Gene + ERV subfamily expression
+ ERV–gene co-expression
+ Telescope locus expression
```

したがって最終統合では、例えば次のようなcandidateが強くなります。

```text
GTEx
AGE-associated module / hub gene
＋
promoter or nearby ERV
＋
Roadmap regulatory annotation

AND

GSE164471
同じgeneがERV subfamilyと強くco-express
＋
Telescopeで対応ERV locusがgene近傍に存在
```

このようなcandidateは、**gene network、ERV expression、ERV locus、genomic proximity、epigenomic context**の複数の独立した情報で支持されることになります。

---

# 最低限の覚え方

| Figure | 一言でいうと |
|---|---|
| `06_GenAge_enrichment_light.pdf` | moduleに既知aging genesが多いか |
| `Step10_ERV_family_heatmap.pdf` | moduleごとにどのERV familyがpromoterに多いか |
| `Step10_nearest_ERV_distance_boxplot.pdf` | module genesとERVの距離 |
| `Step10_promoter_ERV_forest.pdf` | promoter ERV overlapのenrichment/depletion |
| `Step10_promoter_ERV_percent.pdf` | promoter ERVを持つgeneの割合 |
| `Step11_hub_MM_GS_ERV_scatter.pdf` | hub中心性・AGE関連・ERV contextを同時表示 |
| `Step11_hub_promoter_ERV_percent.pdf` | hub genesだけでpromoter ERV率を見る |
| `Step12A_hub_ERV_context_by_module.pdf` | hub-associated ERV contextをmodule別に比較 |
| `Step12A_hub_ERV_genic_context.pdf` | hub-associated ERVがgeneのどこにあるか |
| `Step12B_Roadmap_ERV_family_by_epigenome.pdf` | ERV familyとchromatin stateの関係 |
| `Step12B_Roadmap_hub_enhanc...pdf` | hub-associated ERVのenhancer context |

## 一番重要な注意

GTEx側のFigureは、基本的に

> **ERVがどこにあるか、どのgene/module/hubと位置的・annotation上関係するか**

を示しています。

**ERV expressionそのものはGTEx WGCNAに入っていません。**

ERV expressionとgene expressionの直接的な共発現を見るのはGSE164471側で、そのlocusを細かく分解するのがTelescopeです。
