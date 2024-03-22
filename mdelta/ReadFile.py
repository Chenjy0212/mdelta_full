# -*- coding: utf-8 -*- 
import itertools
import csv
from itertools import islice
import math
import re
import pandas as pd
from pandas.api.types import is_string_dtype
from pandas.api.types import is_numeric_dtype
import numpy as np
import os
from Bio import Phylo

# 去除中间节点和枝长信息
def simplify_tree(clade):
    if clade.is_terminal():
        return clade.name
    else:
        # 递归处理子节点
        children_names = [simplify_tree(child) for child in clade]
        return "(" + ",".join(children_names) + ")"
# 读取 .nwk类文件并把谱系树中的节点名称转换为类型
# 例如   (a1,b7); ---> (a,b);
def ReadTreeSeq_Name2Type(TreeSeqFilePath, Name2TypeFilePath):
    tree = Phylo.read(TreeSeqFilePath, "newick")
    newick_tree = simplify_tree(tree.clade)
    
    file2= open(Name2TypeFilePath,encoding='utf-8')
    for line in file2:
        line = line.replace('"',"")
        line = line.replace("","")
        spl = line.strip().split(',')
        newick_tree = newick_tree.replace(spl[0], spl[1])
    # print(newick_tree)
    return newick_tree.replace(';', '')

def ReadTreeSeq(TreeSeqFilePath):
    tree = Phylo.read(TreeSeqFilePath, "newick")
    newick_tree = simplify_tree(tree.clade)
    return newick_tree

def Scoredict(lllleaf, llllleaf, mav:float, miv:float):
    #如果是自动生成的话
    #可以用到笛卡尔积
    score_dict = {}
    for i in itertools.product(set(lllleaf), set(llllleaf)):
        score_dict[i[0].nodeobj+'_'+i[1].nodeobj] = float(miv)

    #或者用到相同节点才匹配
    #score_dict = {}
    for i in lllleaf:
        score_dict[i.nodeobj+'_'+i.nodeobj] = float(mav)
    for i in llllleaf:
        score_dict[i.nodeobj+'_'+i.nodeobj] = float(mav)
    return score_dict


#定性的
def QualitativeScoreFile(lllleaf, llllleaf, mav:float, miv:float, ScoreFile):
    print('Building Quantitative Score_Dictionary. Loading...正在构建定性矩阵得分字典')
    score_dict={}
    for i in itertools.product(set(lllleaf), set(llllleaf)):
        score_dict[i[0].nodeobj+'_'+i[1].nodeobj] = float(miv)
    for i in lllleaf:
        score_dict[i.nodeobj+'_'+i.nodeobj] = float(miv)
    for i in llllleaf:
        score_dict[i.nodeobj+'_'+i.nodeobj] = float(miv)
    
    df = pd.read_csv(ScoreFile, na_values='NAN', low_memory = False, header = 0)#, index_col = 0)
    
    #定性处理
    for i in range(0, len(df)):
        score_dict[str(df.iat[i,0])+ '_' + str(df.iat[i,1])] = float(df.iat[i,2])
        score_dict[str(df.iat[i,1])+ '_' + str(df.iat[i,0])] = float(df.iat[i,2])
        
    print('Succeed!!!')
    return score_dict

#定量的
def QuantitativeScoreFile(lllleaf, llllleaf, mav:float, miv:float, LScoreFile):
    print('Building Quantitative Score_Dictionary. Loading...正在构建定量矩阵得分字典')
    score_dict={}
    for i in itertools.product(set(lllleaf), set(llllleaf)):
        score_dict[i[0].nodeobj+'_'+i[1].nodeobj] = float(miv)
    for i in lllleaf:
        score_dict[i.nodeobj+'_'+i.nodeobj] = float(miv)
    for i in llllleaf:
        score_dict[i.nodeobj+'_'+i.nodeobj] = float(miv)
    
    df = pd.read_csv(LScoreFile, na_values='NAN', header = 0, index_col = 0)
    row = df.index.tolist()
    col = df.columns.tolist()
    for i in range(len(row)):
        for j in range(len(col)):
            score_dict[str(row[i])+ '_' + str(col[j])] = df.values[i][j]
    print('Succeed!!!')
    return score_dict

## ================== add by lzz ========================
def leafLable_to_celltype_info(node_list):
    # get root1 leaves' new label to celltype infos
    labels = []
    celltypes = []
    for each_node in node_list:
        if "(" not in each_node.nodeobj:
            labels.append(each_node.label)
            celltypes.append(each_node.nodeobj)
    ## conver to strings
    
    return labels, celltypes

#创建文件夹
def mymkdir(path):
    # os.path.exists 函数判断文件夹是否存在
    folder = os.path.exists(path)
    if len(path)!=0:
    # 判断是否存在文件夹如果不存在则创建为文件夹
        if not folder:
            if path[0] == '/':#为绝对路径
                os.makedirs(path)  # makedirs 创建文件时如果路径不存在会创建这个路径
            # os.makedirs 传入一个path路径，生成一个递归的文件夹；如果文件夹存在，就会报错,因此创建文件夹之前，需要使用os.path.exists(path)函数判断文件夹是否存在；
            else:
                b = os.getcwd()
                os.makedirs(b +'/' + path)