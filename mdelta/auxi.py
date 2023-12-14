#!/usr/bin/env python
# -*- coding: utf-8 -*-
def loopindex(row, col):
    ll = []
    min_ = min(row, col)
    for i in range(min_):
        for j in range(i):
            ll.append([j, i])
            ll.append([i, j])
        ll.append([i, i])
    if row > col:
        for j in range(col+1, row+1):
            for i in range(col):
                ll.append([j-1, i])
    elif col > row:
        for j in range(row+1, col+1):
            for i in range(row):
                ll.append([i, j-1])
    return ll
