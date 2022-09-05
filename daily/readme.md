## 待解决问题
### DNBelab C4单细胞平台，双细胞率问题：
#### 1 双细胞检测软件测评
1 常用的检测软件<br>
#### 2 软件测评指标
1 去除双细胞前后的umap/tsne图：双细胞在umap图上的分布<br>
2 各个软件检测结果的一致性<br>
3 去除双细胞前后的轨迹图，也是看双细胞的在轨迹图上的位置<br>
4 去除前后DE基因分析<br>
#### 3 自动化

###部分结果
1. DoubletFinder scrublet doubletdetection scDblFinder鉴定出的双细胞交集很少，不到20%，solo使用CPU耗时，没测评
2. cluster内部的doublet分布没有规律,推断为同种类型细胞组成的双细胞, 
