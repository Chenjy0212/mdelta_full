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


# 读取 .nwk类文件并把谱系树中的节点名称转换为类型
# 例如   (a1,b7); ---> (a,b);
def ReadTreeSeq_Name2Type(TreeSeqFilePath, Name2TypeFilePath):
    file1= open(TreeSeqFilePath,encoding='utf-8') #读取到文本的所有内容
    content=file1.read()
    content=re.sub("\:\d+", "", content)
    content=re.sub("\)\d+", ")", content)
    content=re.sub("\.\d+", "", content)
    content=content.strip()
    #print(content)
    file2= open(Name2TypeFilePath,encoding='utf-8')
    #file2.readline() #跳过第一行
    for line in file2:
        spl = line.strip().split(',')
        content = content.replace(spl[0], spl[1])
    return content.replace(';', '')

def ReadTreeSeq(TreeSeqFilePath):
    file1= open(TreeSeqFilePath,encoding='utf-8') #读取到文本的所有内容
    content=file1.read()
    content=re.sub("\:\d+", "", content)
    content=re.sub("\)\d+", ")", content)
    content=re.sub("\.\d+", "", content)
    content=content.strip()
    #print(content)
    return content.replace(';', '')

def Scoredict(lllleaf, llllleaf, mav:float, miv:float):
    #如果是自动生成的话
    #可以用到笛卡尔积
    score_dict = {}
    for i in itertools.product(set(lllleaf), set(llllleaf)):
        score_dict[i[0].nodeobj+'_'+i[1].nodeobj] = float(miv)
    #score_dict[i[0].nodeobj+'_'+i[1].nodeobj] = random.random()/10
    #print(score_dict)

    #或者用到相同节点才匹配
    #score_dict = {}
    for i in lllleaf:
        score_dict[i.nodeobj+'_'+i.nodeobj] = float(mav)
    for i in llllleaf:
        score_dict[i.nodeobj+'_'+i.nodeobj] = float(mav)
    return score_dict

#def reverseScore(Score, matchScore:float):
#    #Score =matchScore-log(Score+1,math.e)
#    Score =matchScore-(Score**0.5)
#    return Score

#定性的
def QuantitativeScoreFile(lllleaf, llllleaf, mav:float, miv:float, ScoreFile):
    print('Building Score_Dictionary. Loading...')
    score_dict={}
    for i in itertools.product(set(lllleaf), set(llllleaf)):
        score_dict[i[0].nodeobj+'_'+i[1].nodeobj] = float(miv)
    for i in lllleaf:
        score_dict[i.nodeobj+'_'+i.nodeobj] = float(mav)
    for i in llllleaf:
        score_dict[i.nodeobj+'_'+i.nodeobj] = float(mav)
    #typeXn_dict = {}
    
    df = pd.read_csv(ScoreFile, na_values='NAN', low_memory = False, header = 0)#, index_col = 0)
    #print(df)
    #col_name = df.columns.tolist()
    #col_name = list(df)
    #print(col_name)
    
    #定性处理
    for i in range(0, len(df)):
        #if is_numeric_dtype(df.iat[i,2]):
        #print(str(df.iat[i,0]), str(df.iat[i,1]), float(df.iat[i,2]))
        score_dict[str(df.iat[i,0])+ '_' + str(df.iat[i,1])] = float(df.iat[i,2])
        score_dict[str(df.iat[i,1])+ '_' + str(df.iat[i,0])] = float(df.iat[i,2])
        #elif is_numeric_dtype(df.iat[i,1]):
        #    score_dict[str(df.iat[i,0])+ '_' + str(df.iat[i,2])] = float(df.iat[i,1])
        #    score_dict[str(df.iat[i,2])+ '_' + str(df.iat[i,0])] = float(df.iat[i,1])
        #elif is_numeric_dtype(df.iat[i,0]):
        #    score_dict[str(df.iat[i,1])+ '_' + str(df.iat[i,2])] = float(df.iat[i,0])
        #    score_dict[str(df.iat[i,2])+ '_' + str(df.iat[i,1])] = float(df.iat[i,0])
    '''
    csv_reader=csv.reader(open(ScoreFile,'r', encoding="utf-8"))
    #for row in islice(csv_reader, 1, None): #跳过第一行名称信息
    for row in csv_reader:
        #print(row)
        if len(row) == 3:
            if row[2].isdigit():
                score_dict[row[0]+ '_' + row[1]] = float(row[2])
                score_dict[row[1]+ '_' + row[0]] = float(row[2])
        else:
            typeXn_dict[row[0]]=row[1:]
        #print(list(typeXn_dict.keys()))
        for i in itertools.product(list(typeXn_dict.keys()), list(typeXn_dict.keys())):
            #print(i[0]+ '_' + i[1])
            # cmath.sqrt() 返回的是complex复数形式，不利于计算
            score_dict[i[0]+ '_' + i[1]] = reverseScore(sum((abs(float(a)**2-float(b)**2)**0.5) for a,b in zip(typeXn_dict[i[0]],typeXn_dict[i[1]])),math.ceil(len(row)**0.5) if matchScore==-999. else matchScore)
            '''
    return score_dict

#定量的
def QuantitativeScoreFile2(lllleaf, llllleaf, mav:float, miv:float, LScoreFile):
    print('Building Score_Dictionary. Loading...')
    score_dict={}
    for i in itertools.product(set(lllleaf), set(llllleaf)):
        score_dict[i[0].nodeobj+'_'+i[1].nodeobj] = float(miv)
    for i in lllleaf:
        score_dict[i.nodeobj+'_'+i.nodeobj] = 0  #float(mav)
    for i in llllleaf:
        score_dict[i.nodeobj+'_'+i.nodeobj] = 0  #float(mav)
    
    df = pd.read_csv(LScoreFile, na_values='NAN', low_memory = False)
    row = df.index.tolist()
    col = df.columns.tolist()
    for i in range(len(row)):
        for j in range(len(col)):
            score_dict[str(row[i])+ '_' + str(col[j])] = df.values[i][j]
    
    return score_dict
'''
def QuantitativeScoreFile2(lllleaf, llllleaf, mav:float, miv:float, LScoreFile):
    print('Building Score_Dictionary. Loading...')
    score_dict={}
    myp = 0.
    for i in itertools.product(set(lllleaf), set(llllleaf)):
        score_dict[i[0].nodeobj+'_'+i[1].nodeobj] = float(miv)
    for i in lllleaf:
        score_dict[i.nodeobj+'_'+i.nodeobj] = 0  #float(mav)
    for i in llllleaf:
        score_dict[i.nodeobj+'_'+i.nodeobj] = 0  #float(mav)
    
    #把不对等的行数处理
    df1 = pd.read_csv(ScoreFile1, na_values='NAN', low_memory = False, header = 0, index_col = 0)
    df2 = pd.read_csv(ScoreFile2, na_values='NAN', low_memory = False, header = 0, index_col = 0)
    #df2.rename(columns=lambda x:x+'@#$%', inplace=True)
    #print(df2.columns.values)
    
    df = pd.concat([df1,df2],axis=1).fillna(0) #不存在的行设为0
    df1 = df.iloc[:,:df1.shape[1]]
    df2 = df.iloc[:,df1.shape[1]:]
    #print(df1.shape[1], df2.shape[1], df.shape[1])
    
    dis_matrix = pd.DataFrame([[0.0 for i in range(len(df2.columns.values))] for j in range(len(df1.columns.values))],
                        index=[i for i in df1.columns.values],
                        columns=[j for j in df2.columns.values])
    
    for loci, i in enumerate(df1.columns.values):
        for locj, j in enumerate(df2.columns.values):
            #print(sum(df1[i] - df2[j]))
            dissst = np.linalg.norm(df1[i]-df2[j])
            score_dict[str(i)+ '_' + str(j)] = dissst
            dis_matrix.values[loci][locj] = dissst
            #score_dict[str(i)+ '_' + str(j)] = sum((df1[i]-df2[j])^2)#np.sqrt(sum(np.power((df1[i]-df2[j]), 2)))
            #score_t = np.corrcoef(df1[i], df1[j])
            #score_dict[str(i)+ '_' + str(j)] = np.corrcoef(df1[i], df1[j])
            #print(score_t, end = " ")
    
    # 保存距离矩阵
    # dis_matrix.to_csv(TreeSeqFileName + '_' + TreeSeqFileName2 + '_dis_matrix.csv')
    
    score_list = list(score_dict.values())
    #print(score_list)
    score_list.sort()
    #print(len(score_list))
    #print(round(len(score_list) * 0.05))
    myp = score_list[round(len(score_list) * (1-pp/100))]
    
    for key,val in score_dict.items():
        score_dict[key] = math.log((myp - val) + 1) if myp >= val else -math.log((val - myp) + 1) #myp/(val + 1) if myp > val else -(val/(myp + 1)) + pv
        #score_dict[key] = myp - val
    #score_list = list(score_dict.values())
    #print(score_list)
    
    return score_dict
'''
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