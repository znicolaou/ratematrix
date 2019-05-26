#!/usr/bin/env python
import numpy as np
import tlist
import timeit

start=timeit.default_timer()
test=np.zeros(4);
test=test+1
tlist.list(test.astype(int))
print(timeit.default_timer()-start)
