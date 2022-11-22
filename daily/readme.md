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


### 2022/09/13
#### 学东西就是为了用，所以就是要有明确的目的性
咋逐渐意识到我自个很菜呢？线性代数很菜，矩阵的秩，分解，乘法，学单细胞才意识到，机器学习，看了很多遍入门，半途放弃，随机森立呀，决策树呀，SVM呀，很感兴趣就是工作中没机会用，也没法精进；RNA速率呀回忆起微积分，重拾吧？python基础语法看了n遍了？<br>
我今后的生信之路就是学好数学和编程，然后知识的归纳提炼和输出，一方面是工作，另一方面和自己的对知识的理念相契合，做一个有不断探索原理的人，持续输出能力的人，因为热爱这个吧？<br>
https://github.com/Avik-Jain/100-Days-Of-ML-Code

http://lihuaxi.xjx100.cn/news/102615.html
### 读源码，搞懂单细胞数据整合的原理，重点SCTransform，harmony和liger

### 多样本单细胞批次效应，怎么定量评估，umap图欠缺定量的描述

###
sed -i 's/[[:space:]]\{1,\}/,/g' gene_counts.txt
#横坐标轴标签倾斜
axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5)

## 了解蛋白组
https://bioincloud.tech/cloudir/reports/TMT_report/TMT_Report.html<br>
