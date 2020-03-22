# python-处理hmp文件用于QTL分析

## 原始文件：

## rawdata_hmp.xlsx	基因型文件

## all.filter.RQTL.geno.txt	含有遗传距离的文件

## MR_phenotype.xlsx	表型文件



### 用R转换为txt格式

```R
#!/usr/bin/env R

install.packages("openxlsx")
library(openxlsx)
rawdata <- read.xlsx("rawdata_hmp.xlsx",sheet = 1)
write.table(rawdata,"output.txt",quote = FALSE,row.names = FALSE) #不输出行名以及字符串的双引号

######usage#####
write.table(x, file = "", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
```



### 用python进行数据分析时使用的是R的结果文件：output.txt，也可直接用python读取原始xlsx文件

```python
#!/usr/bin/env python
# -*- coding: utf-8 -*-
#fdy 2020-3-22

import pandas as pd  #导入模块并命名为pd
import collections as col

################################将遗传距离添加到基因型文件#######################
data = pd.read_table(r'D:\yanglab_data\QTL\dachucao\output.txt',sep=' ') #读取文件
#print(data)
data2 = pd.read_table(r'D:\yanglab_data\QTL\dachucao\all.filter.RQTL.geno.txt',nrows=4,header=None) #读取文件前四行
data2_T = data2.T #转置数据
data2_T.iloc[0,3] = "fdy" #将cM单元格重命名为fdy，cM为遗传距离
#print(data2_T)
data2_T.to_csv(r'D:\yanglab_data\QTL\dachucao\rawdata_cM.csv',index = None,header=None) #将结果输出到csv文件，并删除行列索引
data3 = pd.read_csv(r'D:\yanglab_data\QTL\dachucao\rawdata_cM.csv',sep=',')
#print(data3)
data4 = data3.loc[:,['rs','fdy']] #读取rs和fdy两列
#print(data4)
data5 = col.OrderedDict(zip(data4.iloc[:,0],data4.iloc[:,1])) #将rs，fdy构建为字典，形式为{rs：fdy}
data6 = dict(data5)
col_name = data.columns.tolist() #提取data文件的列索引
#print(col_name)
col_name.insert(3,'fdy') #在第四列插入一列，名为fdy
#print(col_name)
data = data.reindex(columns=col_name) #重建data的索引
data['fdy'] = data['rs'].map(data6) #为fdy这一列添加值，以data中rs列的名为key，在字典data6中寻找其对应的值，类似excel的vlookup功能
#print(data)
data = data.dropna(subset=["fdy"]) #删除fdy这一列中为NaN的行
#print(data)
data_T = data.T
#print(data_T)
data_T.to_csv(r'D:\yanglab_data\QTL\dachucao\final_hmp.csv',header=None)
hmp = pd.read_csv(r'D:\yanglab_data\QTL\dachucao\final_hmp.csv')
#print(hmp)

#####################################插入表型数据#############################
pt = pd.read_excel(r'D:\yanglab_data\QTL\dachucao\MR_phenotype.xlsx') #读取表型文件
#print(pt)
pt_dic = col.OrderedDict(zip(pt.iloc[:,0],pt.iloc[:,1])) 
data7 = dict(pt_dic)
#print(data7)
col_name = hmp.columns.tolist()
#print(col_name)
col_name.insert(1,'phenotype_water')
#print(col_name)
hmp = hmp.reindex(columns=col_name)
hmp['phenotype_water'] = hmp['rs'].map(data7)
#print(hmp)
hmp.drop([0, 3],inplace = True) #删除第1，4行，设置inplcae则可直接修改原文件
hmp.iloc[0,1] = " " #设置单元格为空
hmp.iloc[1,1] = " "
hmp.drop(['rs'],axis=1,inplace=True) #删除rs列
hmp = hmp.dropna(subset=["phenotype_water"])
hmp_rep = hmp.replace({'M':'AC','C':'CC','A':'AA','N':'-'}) #将M替换为AC
print(hmp_rep)
hmp_rep.to_csv(r'D:\yanglab_data\QTL\dachucao\results.csv',index = None)


```

```python
#################################其它可能用到的一些方法##############################

import pandas as pd
import collections as col

data = pd.read_table(r'D:\yanglab_data\QTL\dachucao\all.filter.RQTL.geno.txt',nrows=4)
data2 = pd.read_csv(r'D:\yanglab_data\QTL\dachucao\output.csv',sep=',',header=None,low_memory=False)

####获取M00001-M10000的列####
file = open(r'D:\yanglab_data\QTL\dachucao\fdy.txt','w')
file.write('num'+'\n')
for i in range(len(data2)):
    file.write('m'+("{:0>5d}".format(i))+'\n')
file.close()


data3 = pd.read_table(r'D:\yanglab_data\QTL\dachucao\fdy.txt')
print(data3)
index = data3.pop('num')
data2['num'] = index
#cols = data.shape[1] #计算列数
#rows = data.shape[0]
#print(rows)
#print(cols)
data_T = data.T
print(data_T)
#data2['num'] = range(len(data2)) #添加行
#print(data)
print(data2)


#############################构建字典的其它方法#######################
result_dic = data4.set_index('rs').T.to_dict()
result_dic = data4.groupby('rs')['fdy'].apply(list).to_dict()
```

