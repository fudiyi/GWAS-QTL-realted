# python-统计数据是否符合正态分布

#### 适用于：

1. #####  检测表型数据

2. #####  共表达网络计算 PCC 后进行统计

```python
#!/usr/bin/env ipython

import pandas as pd
import numpy as np
import matplotlib as plt

data = pd.read_csv('file.csv')
data1 = np.array(data['phenotype or PCC']) # 选取相对应的数据读入numpy
bins = np.linspace(-1,1,20) # 设置直方图区间
plt.hist(data1,bins,histype='step')
plt.show()
```

i