dotgroupplot
======

### What is it?


### Installation

1. Clone the [GitHub repo](https://github.com/mourisl/dotgroupplot), e.g. with `git clone https://github.com/mourisl/dotgroupplot`.
2. Copy "dotgroupplot" folder to your project folder or run "python3 setup.py install" to install.

I will try to add dotgroupplot to PyPi in future.

### Usage
Here is a minimal example:

```python
import math
import pandas as pd
from dotgroupplot import dotgroupplot

n = 100 # Draw 100 groups
idlist = []
conditionlist = []
for i in range(n):
    idlist += [str(i)] * (i + 1)
    conditionlist += ["Healthy"] * math.floor((i+1)/2) + ["Disease"] * math.ceil((i+1)/2)
df = pd.DataFrame({"CDR3":idlist, "Condition":conditionlist})    

dotgroupplot.dotgroupplot(df, "CDR3", hue="Condition")
```

See more examples in [example.ipynb](https://github.com/mourisl/dotgroupplot/blob/main/example.ipynb).

### Requirements
+ Python >= 3.6
+ seaborn >= 0.9
+ matplotlib >= 2.2.2
