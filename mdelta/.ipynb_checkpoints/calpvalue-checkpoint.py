#!/usr/bin/env python
# -*- coding: utf-8 -*-

from copy import *
import random

def RandomlyShuffle(Seq: str, times: int):
    result = []
    if Seq.find('(') == -1:
        for i in range(times):
            result.append(Seq)
        return result
    Seq_list = []
    index_list = []
    index_dict = {}
    brackets_num = 0  # 括号个数
    node_num = -1  # 节点序号
    node_tmp = ''
    
    for i in list(Seq):
        if i == '(':
            brackets_num += 1
            Seq_list.append(i)
            node_num += 1
        elif i == ',':
            if brackets_num != 0 and node_tmp:
                Seq_list.append(node_tmp)
                node_num += 1
                index_list.append(node_num)
                index_dict[node_num] = node_tmp
                node_tmp = ''
            Seq_list.append(i)
            node_num += 1
        elif i == ')':
            if node_tmp:
                Seq_list.append(node_tmp)
                node_num += 1
                index_list.append(node_num)
                index_dict[node_num] = node_tmp
                node_tmp = ''
            brackets_num -= 1
            Seq_list.append(i)
            node_num += 1
        else:
            node_tmp += i
            if i == ';':
                Seq_list.append(i)
                node_num += 1
    # print(Seq_list)
    # print(index_dict)
    # print(index_list)
    Seq_list_tmp = deepcopy(Seq_list)
    index_list_tmp = deepcopy(index_list)
    #print('原来的序列: \n',"".join(Seq_list))
    #print('改变的序列: \n')
    for i in range(times):
        random.shuffle(index_list_tmp)
        # print(index_list_tmp)
        for i, j in zip(index_list_tmp, index_list):
            Seq_list_tmp[j] = index_dict[i]
        result.append("".join(Seq_list_tmp))
    # print(result)
    return(result)

# 计算随机打乱位置的，定性计算的，得分矩阵
def get_rmatrix(tree1: str, tree2: str, rqueue):
    print(tree1)
    print(tree2)
    rqueue.put(tree2)
