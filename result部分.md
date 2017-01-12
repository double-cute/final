SpliceNet的主要优点：
1. 可以处理高维数据.
  - 特别当[${exons per gene}:${sample size}]比例值特别大时表现良好.


## 模拟实验：

1. varying number of exons(dimensions) and samples.
  - 以exon和sample的数量作为变量进行实验.
    - exon数量决定了矩阵的维数.
2. 论文给出的R软件包中的RNASeqNet模块的性能也是基于以上数据测试的.

1. Table1显示出了RNASeqNet从exon表达数据中抽象出依赖关系的突出能力.
2. RNASeqNet和SpliceNet的测试和运行都是基于特定癌症的基因：ERBB2、MAPK
  - 数据来自KEGG的singaling pathways.
  - 采用不同数量的sample.
3. Figure3显示出了在低sample量的数据上SpliceNet比RNASeqNet的表现更加出众.
4. Figure4/5是非常细节的SpliceNet在Bcl-x和EGFR中心网络上的结果.
  - Differential edges inferred by SpliceNet converged to cancer-specific splice variants reported in literature
  - 以上用SpliceNet预测出的normal-cancer diff的边在已有的医学文献中已经被证实.
    - 说明预测准确.

## 真实数据上的实验：

1. SpliceNet和RNASeqNet都在在真实的RNA-Seq数据上跑过.
  - 器官有lung、kidney、liver.
2. 结果在Table4上，显示出SpliceNet表现要比RNASeqNet好很多.

## signaling pathways验证：

- 为了可比性，SpliceNet和RNASeqNet一样，是基于non-small cell lung cancer-specific pathway验证的.
- 验证用的ERBB2/MAPK的signaling pathway从KEGG上下载下来的.

1. 45个LUSC sample用来预测边，结果在Figure3a中给出.


## isoform-diff-net在non-small-cell LUADsample上的结果总结

- 为了深入理解diff-net的优势以及应用价值，一个更细致的结果：SpliceNet在Bcl-x和EGFR中心网络上的结果
  - Figure4/5

1. 简单介绍Bcl-x基因：
  - 和绝大多数non-small-cell lung cancer有关.
  - 两个isform：Bcl-xL（抑制细胞凋亡）、Bcl-xS（促进细胞凋亡）

2. 和Bcl-x有关的两个基因SIVA1和CFLAR，三者之间的关系在normal和cancer中都显著存在.
  - 特别是在肿瘤发生（凋亡循环）过程中扮演这非常重要的作用.
  - Figure4a that there is no difference between the networks derived from cancer and normal samples, and is difficult to explain tumorigenesis.
    - 由于从基因交互网络中看不出正常和癌症之间的差异（两张图完全一样），因此必须从更深入的isoform层面去解释.

3. 同样，EGFR中心网络中，在cancer和normal状况下都有和CD44和CEACAM1有交互，从基因交互网络中无法描述cancer/normal之间的差异.
  - 但是从diff-net中可以明显的看出来.
  - Figure5/a/

- Figure4和5的a都是KEGG上的signaling pathways，
