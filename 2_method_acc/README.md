## 寻找最佳参数时同时计算准确率使用的数据

### 仅展示bulk数据

1. TCGA-LUAD
2. GSE42568
3. TCGA-BRCA
4. GSE9891 (恶性 vs.低恶性潜能) (不视为肿瘤vs.非肿瘤)

查下图以查看数据匹配：

![data_match](paper_plot/plot/benchmarkdata_sankey.png)

### 真实标签

1. GSE161529 HER2 有肿瘤标签
2. GSE161529 TNBC 有肿瘤标签
3. GSE123902 有肿瘤标签
4. GSE165897 有肿瘤标签

### 说明

因为真实判断标签标识的均为肿瘤，所以二元标签也用肿瘤。

与生存与预后差相关的细胞不一定是肿瘤细胞，所以不纳入。