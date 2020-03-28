# python脚本汇总

### 一些python脚本汇总，由JYX编写

### 总结：数据分析常用for循环和def函数进行批量操作

#### add_title.py

```python
#!usr/bin/python3
#此脚本用于给txt文件添加列名
output = open('p_morethan_5.annotation','w') 
output.write('chr' + '\t' + 'pos'+'\t'+ 'log10(P)' + '\t' + 'compound_num' + '\t' + 'gene_start' + '\t' + 'gene_end'+'\t'+'gene'+'\t'+'description'+'\n' )

for i in open('p_morethan_5.intersect'):
    ii = i.strip().split('\t') #将文件按空格分列
    output.write(ii[0]+'\t'+ii[1]+'\t'+ii[3]+'\t'+ii[4]+'\t'+ii[6]+'\t'+ii[7]+'\t'+ii[8]+'\t'+ii[9]+'\n' ) #写入列名所需的内容

print('DONE')

```

#### count_gene_number.py

```python
#!usr/bin/python3
dic_num = {} #写一个空字典，用于存放数据
for i in open('p_morethan_5.annotation'):
    ii = i.strip().split('\t')
    if ii[6] in dic_num:
        dic_num[ii[6]] = dic_num[ii[6]] + 1
    else:
        dic_num[ii[6]] = 1
output = open('gene_replicate.txt','w')
for i in dic_num:
    output.write(i+'\t'+str(dic_num[i])+'\n')
print('DONE')

```

#### draw_scatter.py

```python
#!usr/bin/python3
#coding=utf-8
#这个脚本用来将“get_target_range.py”找出的一系列点画成散点图
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math

def Draw_figure(fl):
    fl_name = fl + '.txt.output'
    fig_name = fl + '.png'
    pos = []
    p = []
    for i in open(fl_name):
        ii = i.strip().split('\t')
        pos.append(int(ii[3]) - 154974052)
        p.append(math.log10(float(ii[6]))*-1)
    plt.figure(figsize=(20,8),dpi = 80)
    plt.scatter(pos,p, c = 'orange')
    plt.ylabel('-log(p)',fontsize=18)
    plt.title(fl,fontsize=20)  #图片标题
    plt.savefig(fig_name)  #图片名称
    plt.close(fig_name)
    print('DONE',fig_name,'!')

fl_ls = ['125','167','28','320','35','36']
for j in fl_ls:
    Draw_figure(j)

print('ALL DONE')

```

#### find_P_for_specific.py

```python
#!usr/bin/python3
#coding=utf-8
import math

def get_p(n):
    fl_name ='/data/FDY_analysis/mGWAS/analysis/GWAS/test/' + n + '.txt'
    fl = open(fl_name)
    fl.readline()
    output_fl = n + '.output'
    output = open(output_fl,'w')
    for i in fl:
        ii = i.strip().split('\t')
        p = math.log10(float(ii[6])) * -1
        if p > 5:
            output.write(ii[2] + '\t' + ii[3] + '\t' + str(int(ii[3]) + 1) +'\t' +str(p) + '\n')
    output.close()
    print('finish ',output_fl)

for j in open('list.txt'):
    j = j.strip()
    get_p(j)

```

#### format_tassel_results_batch.py

```python
#!usr/bin/python3
#coding=utf-8
#这个脚本用来将tassel的结果作图，这样就不用再把数据传回win电脑上作图了。做出的图不如那个R脚本好看，但可以做参考。
#运行过程中会产生一个tmp.txt文件用于临时存放一些数据。
#180529修改，可以用于批量作图。
import datetime
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

start = datetime.datetime.now()
def format_data(n):
    txt = str(n) + '.txt'
    path = '/data/GWAS/170515_daixie/warm/results/' + txt  #批量作图时，注意修改数据文件所在路径
    res = open(path)
    print('reading data, please waite')
    output = open('tmp.txt','w')
    a = res.readline()
    a = res.readline()
    end1 = 301284077
    end2 = 237032434 + end1
    end3 = 232132502 + end2
    end4 = 241403054 + end3
    end5 = 217807808 + end4
    end6 = 169151979 + end5
    end7 = 176396028 + end6
    end8 = 175788001 + end7
    end9 = 156637819 + end8
    end10 = 150177341 + end9
    for i in res:
        ii = i.strip().split('\t')
        p = str(math.log10(float(ii[6])) * -1)
        if ii[2] == '1':
            pos = ii[3]
            output.write(ii[2]+'\t'+pos+'\t'+p+'\n')
        elif ii[2] == '2':
            pos = str(int(ii[3]) + end1)
            output.write(ii[2] + '\t' + pos + '\t' + p + '\n')
        elif ii[2] == '3':
            pos = str(int(ii[3]) + end2)
            output.write(ii[2] + '\t' + pos + '\t' + p + '\n')
        elif ii[2] == '4':
            pos = str(int(ii[3]) + end3)
            output.write(ii[2] + '\t' + pos + '\t' + p + '\n')
        elif ii[2] == '5':
            pos = str(int(ii[3]) + end4)
            output.write(ii[2] + '\t' + pos + '\t' + p + '\n')
        elif ii[2] == '6':
            pos = str(int(ii[3]) + end5)
            output.write(ii[2] + '\t' + pos + '\t' + p + '\n')
        elif ii[2] == '7':
            pos = str(int(ii[3]) + end6)
            output.write(ii[2] + '\t' + pos + '\t' + p + '\n')
        elif ii[2] == '8':
            pos = str(int(ii[3]) + end7)
            output.write(ii[2] + '\t' + pos + '\t' + p + '\n')
        elif ii[2] == '9':
            pos = str(int(ii[3]) + end8)
            output.write(ii[2] + '\t' + pos + '\t' + p + '\n')
        elif ii[2] == '10':
            pos = str(int(ii[3]) + end9)
            output.write(ii[2] + '\t' + pos + '\t' + p + '\n')
    output.close()
    print('finish format the data of ',txt)
#    print('generating figure')
data = ''
def get_data(n):
    global data
    fl = data[data.chr == n]
    x = fl['pos']
    y = fl['p']
    return(x,y)
def generate_figure(n):
    txt = str(n) + '.txt'
    path = '/data/GWAS/170515_daixie/warm/results/' + txt  #批量作图时，注意修改数据文件所在路径
    global data
    data = pd.read_table('tmp.txt',names=['chr','pos','p'])
    x1,y1 = get_data(1)
    x2,y2 = get_data(2)
    x3,y3 = get_data(3)
    x4,y4 = get_data(4)
    x5,y5 = get_data(5)
    x6,y6 = get_data(6)
    x7,y7 = get_data(7)
    x8,y8 = get_data(8)
    x9,y9 = get_data(9)
    x10,y10 = get_data(10)
    plt.figure(figsize=(20,8),dpi = 80)
    plt.scatter(x1,y1, c = 'orange',alpha=0.5)
    plt.scatter(x2,y2, c = 'blue',alpha=0.5)
    plt.scatter(x3,y3, c = 'orange',alpha=0.5)
    plt.scatter(x4,y4, c = 'blue',alpha=0.5)
    plt.scatter(x5,y5, c = 'orange',alpha=0.5)
    plt.scatter(x6,y6, c = 'blue',alpha=0.5)
    plt.scatter(x7,y7, c = 'orange',alpha=0.5)
    plt.scatter(x8,y8, c = 'blue',alpha=0.5)
    plt.scatter(x9,y9, c = 'orange',alpha=0.5)
    plt.scatter(x10,y10, c = 'blue',alpha=0.5)
    plt.title(txt,fontsize=20)  #图片标题
    plt.ylabel('-log(p)',fontsize=18)
    figname = txt + '.png'
    plt.savefig(figname)  #图片名称
    print('finish generating figure ',txt)

i = 1
while(i <= 1407):
    format_data(i)
    generate_figure(i)
    i += 1

end = datetime.datetime.now()
used = (end-start).seconds
print(used,'seconds were used')
print('done')

```

#### get_a_region_GWAS_results.py

```python
#!usr/bin/python3
import math
fl = ['40','50','195']
start_site = 172159342
end_site = 172164107
chro = 7
start = start_site - 5e5
end = end_site + 5e5
output = open('GRMZM2G081682_cold.txt','w')
for i in fl:
    fl_name = i + '_output2.txt'
    res = open(fl_name)
    res.readline()
    res.readline()
    for j in res:
        jj = j.strip().split('\t')
        if chro == int(jj[2]):
            if int(jj[3]) > start and int(jj[3]) < end:
                p = math.log10(float(jj[6])) * -1
                output.write(jj[0] + '\t' + jj[2] + '\t' + jj[3] + '\t' + str(p) + '\n')
    print('Finish',fl_name)
print('DONE')

```

#### get_annotation.py

```python
#!usr/bin/python3
#coding=utf-8
dic_anno = {}
for i in open('/data/RNA-Seq/maize_annotation/v3_gene_annotation.txt'):
    ii = i.strip().split('\t')
    dic_anno[ii[0]] = ii[1]
print('Dic_anno done.')

output = open('/data/JYX/get_annotation_180517/gene_list.annotation.txt','w')
for i in open('/data/JYX/get_annotation_180517/gene_list.txt'):
    i = i.strip()
    if '？' in i:
        ii = i.split('？')
        ii[0] = ii[0].strip()
        output.write(i + '\t' + ii[0] + '\t' + dic_anno.get(ii[0], 'NA') + '\n')
    elif 'GRMZM' in i and '_' in i:
        ii = i.split('_')
        ii[0] = ii[0].strip()
        output.write(i + '\t' + ii[0] + '\t' + dic_anno.get(ii[0], 'NA') + '\n')
    elif 'AC' in i and 'T' in i:
        ii = i.split('T')
        j = ''.join(ii)
        output.write(i + '\t' + j + '\t' + dic_anno.get(j, 'NA') + '\n')
    elif ' ' in i:
        ii = i.split(' ')
        ii[0] = ii[0].strip()
        output.write(i + '\t' + ii[0] + '\t' + dic_anno.get(ii[0], 'NA') + '\n')
    else:
        output.write(i + '\t' + i + '\t' + dic_anno.get(i, 'NA') + '\n')
print('DONE')

```

#### get_cold_gene_GWAS_results.py

```python
#!usr/bin/python3
import math

def get_results(n):
    nn = n.strip().split('\t')
    fl = nn[0] + '_cold.txt'
    output = open(fl,'w')
    chro = int(nn[1])
    start_site = int(nn[2])
    end_site = int(nn[3])
    start = start_site - 5e5
    end = end_site + 5e5
    meta_ls = nn[4:]
    for i in meta_ls:
        fl_name = i + '_output2.txt'
        res = open(fl_name)
        res.readline()
        res.readline()
        for j in res:
            jj = j.strip().split('\t')
            if chro == int(jj[2]):
                if int(jj[3]) > start and int(jj[3]) < end:
                    p = math.log10(float(jj[6])) * -1
                    output.write(jj[0] + '\t' + jj[2] + '\t' + jj[3] + '\t' + str(p) + '\n')
        print('Finish',fl_name)
    print('Finish',fl,'===========')
    output.close()

for m in open('cold_gene_region.txt'):
    get_results(m)

print('DONE')

```

#### get_cold_genes_region.py

```python
#!usr/bin/python3
dic_pos = {}
for i in open('/data/maize_gff_and_bed/v2/ZmB73_5a_gene_annotation.bed'):
    ii = i.strip().split('\t')
    gene_info = ii[0] + '\t' + ii[1] + '\t' + ii[2]
    dic_pos[ii[3]] = gene_info

output = open('cold_gene_region.txt','w')
for i in open('cold.txt'):
    ii = i.strip().split()
    try:
        gene_info = dic_pos[ii[0]]
        meta = '\t'.join(ii[1:])
        output.write(ii[0] + '\t' + gene_info + '\t' + meta + '\n')
    except:
        pass
output.close()
print('DONE')

```

#### get_gene_description.py

```python
#!usr/bin/python3
dic_des = {}
for i in open('TAIR10_functional_descriptions_20140331.txt'):
    ii = i.strip().split('\t')
    ii[0] = ii[0].split('.')[0]
    des = '\t'.join(ii[1:])
    dic_des[ii[0]] = des
output = open('diff.anno','w')
for i in open('diff'):
    i = i.strip()
    output.write(i+'\t'+dic_des[i]+'\n')
print('done')

```

#### get_intersectBed_batch.py

```python
#!usr/bin/python3
output = open('intersect.sh','w')
for i in open('list.txt'):
    i = i.strip()
    fl = i + '.output'
    fl_inter = i + '.intersect'
    output.write('intersectBed -wo -a '+ fl + ' -b /data/maize_gff_and_bed/v2/ZmB73_5a_gene_annotation.bed > ' + fl_inter + '\n')

```

#### get_p_morethan_5.py

```python
#!usr/bin/python3
#coding=utf-8
#这个脚本用来从这一堆GWAS结果里，把P大于5的点都找出来，并以BED格式输出，方便后续与基因对应并加注释
import math

output = open('p_morethan_5.txt','w')
output.write('chr' + '\t' + 'pos1' + '\t' + 'pos2' + '\t' + 'log10(P)' + '\t' + 'compound_num' + '\n')
def get_p(n):
    fl_name =str(n) + '_output2.txt'
    print('Processing',fl_name)
    fl = open(fl_name)
    fl.readline()
    fl.readline()
    for i in fl:
        ii = i.strip().split('\t')
        p = math.log10(float(ii[6])) * -1
        if p > 5:
            output.write(ii[2] + '\t' + ii[3] +'\t' +  str(int(ii[3]) + 1) + '\t' + str(p) +'\t' + str(n) +  '\n')

i = 1
while i <= 350:
    get_p(i)
    i += 1

print('DONE')

```

#### get_target_metabolism.py

```python
#!usr/bin/python3
output = open('target_metabolism','w')
dic_meta = {}
meta_total = open('metabolism-JYX.txt')
output.write(meta_total.readline())
for i in meta_total:
    ii = i.strip().split('\t')
    data = '\t'.join(ii[1:])
    dic_meta[ii[0]] = data
print('Finish dic_meta')

for i in open('metabolism_list'):
    i = i.strip()
    meta = dic_meta.get(i,'NA')
    output.write(i + '\t' + meta + '\n')
print('DONE')

```

#### get_the_metabolism_name.py

```python
#!usr/bin/python3
#coding=utf-8
#从代谢组GWAS的结果中找出ICE1所在区间P大于5的点，并将这个点和所对应的代谢产物输出。
output = open('cold_907.txt','w')
def get_name(name):
    target = []
    fl = str(name)+'.txt'
    print('Now processing',fl)
    for i in open(fl):
        if 'chr3' in i:
            ii = i.strip().split('\t')
            if int(ii[3]) >= 156362362 and int(ii[3]) <= 156365307:
                if float(ii[6]) <= 0.00001:
                    target.append(fl+'\t'+i)
        elif 'chr4' in i:
            return target
            break

res = []
for i in range(1,908):
    res.extend(get_name(i))
#print(res)
dic_name = {}
for i in open('metabolism_name_907.txt'):
    ii = i.strip().split('\t')
    dic_name[ii[0]] = ii[1]
#print(dic_name)
for i in res:
    ii = i.strip().split('\t')
    name = ii[0].split('.')
    meta = dic_name[name[0]]
    output.write(meta+'\t'+i)
print('Done')

```

#### change_to_log.py

```python
#!usr/bin/python3
import math
output = open('wc2.log10_p','w')
for i in open('wc2.ps'):
    ii = i.split('\t')
    pos = ii[0].split('_')
    logp = math.log(float(ii[2]),10) * -1
    output.write('\t'.join(pos)+'\t'+str(logp)+'\n')
print('done')

```

#### generate_figure_for_resequencing_GWAS.py

```python
#!usr/bin/python3
#conding=uft-8
#这个脚本用来将重测序GWAS做的结果一一画图
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
def gener_figure(prefix,name_ls):
    for i in name_ls:
        name_fl = prefix + '_' + i + '_output2.txt'
        fig_name = prefix+'_'+i+'.png'
        fig_titl = prefix+'_'+i
        fl = open(name_fl)
        a = fl.readline()
        a = fl.readline()
        pos = []
        p = []
        for a in fl:
            aa = a.strip().split('\t')
            pos.append(int(aa[3]))
            p.append(math.log10(float(aa[6]))*-1)
        plt.figure(figsize=(20,8),dpi = 80)
        plt.scatter(pos,p, c = 'orange')
        plt.ylabel('-log(p)',fontsize=18)
        x_ticks = np.arange(100,5100,100)
        y_ticks = np.arange(0,9,0.5)
        plt.xticks(x_ticks,rotation=45)
        plt.yticks(y_ticks)
#        plt.ylim((0,10))
        plt.title(fig_titl,fontsize=20)  #图片标题
        plt.savefig(fig_name)  #图片名称
        plt.close(fig_name)
        print('DONE',fig_name,'!')

name_ls1 = ['19','28','109','167','223','320','344','384','576','596','681','768','788','792']
name_ls2 = ['137','196','218','294','318','374','380','507','540','607','633','721','854','855','907','918','921','964','986','990','1069','1077','1132','1175','1280','1306','1339']
name_ls3 = ['541','742','793','925','1035','1362','1407']

gener_figure('907_cold',name_ls1)
gener_figure('1407_cold',name_ls2)
gener_figure('1407_warm',name_ls3)

print('ALL DONE')

```

#### wig_normalised.py

```python
import pandas as pd
names = ['position','reads']
df =  pd.read_table('/data/FDY/ChIPSeq/181125_A00403_0097_AHFGTCDSXX/macs_result2/maize_ChIP_MACS_wiggle/treat/maize_ChIP_treat_afterfiting_all.wig',sep='\t',names=names,low_memory=False)
sum=df['reads'].sum(skipna=True)
df['reads']=df['reads']*1000000/sum
df.to_csv('ZmCT.wig',sep='\t',index=False,header=False)
df.describe()

names = ['position','reads']
df = pd.read_table('/data/FDY/ChIPSeq/181125_A00403_0097_AHFGTCDSXX/macs_result2/maize_ChIP_MACS_wiggle/control/maize_ChIP_control_afterfiting_all.wig',sep='\t',names=names,low_memory=False)
sum=df['reads'].sum(skipna=True)
df['reads']=df['reads']*1000000/sum
df.to_csv('WT.wig',sep='\t',index=False,header=False)
df.describe()
```

