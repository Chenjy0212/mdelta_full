#!/usr/bin/env python
# coding: utf-8
# def install(package):
#     subprocess.check_call(["pip", "install", package])

# # 要安装的包列表
# packages = ['matplotlib']

# # 循环遍历包列表，如果没有安装，则使用pip安装该包
# for package in packages:
#     try:
#         __import__(package)
#     except ImportError:
#         print(f"{package} not found. Installing...")
#         install(package)

import sys
import json
import networkx as nx
from networkx.algorithms import bipartite
import matplotlib.pyplot as plt
import os


# 保存的文件夹
folder = os.path.abspath(sys.argv[2]) + '/network/'
if not os.path.exists(folder):
    os.makedirs(folder)

with open(sys.argv[1],'r',encoding='utf8') as fp1:
    json_dataA1G11 = json.load(fp1) 
if len(sys.argv) == 4 and sys.argv[2] != 'non' and sys.argv[3] != 'non':
    with open(sys.argv[2],'r',encoding='utf8') as fp2:
        json_dataA1A1 = json.load(fp2)
    with open(sys.argv[3],'r',encoding='utf8') as fp3:
        json_dataG11G11 = json.load(fp3)

    for i in json_dataA1A1:
        i['Root1_label_tmp'] = i.pop('Root2_label')
    for i in json_dataG11G11:
        i['Root2_label_tmp'] = i.pop('Root1_label')

    json_data = json_dataA1G11 + json_dataA1A1[1:] + json_dataG11G11[1:]
else:
    json_data = json_dataA1G11

X = []
Y = []
leaves = {}
B = nx.Graph()
for i in json_data:
    if i.get('Root1_label', False) and i.get('Root2_label', False):
        T1_label = 'T1_T2_' + str(i['Root1_label'])
        T2_label = 'T2_T1_' + str(i['Root2_label'])
        X.append(T1_label)
        Y.append(T2_label)
        leaves[T1_label] = i['Root1_node'].count(',') + 1
        leaves[T2_label] = i['Root2_node'].count(',') + 1
        B.add_edge(T1_label, T2_label, score = i['Score'])
        
    elif i.get('Root1_label_tmp', False):
        T1_label = 'T1_T1_' + str(i['Root1_label'])
        T2_label = 'T1_T1_' + str(i['Root1_label_tmp'])
        #print(T1_label, T2_label)
        X.append(T1_label)
        X.append(T2_label)
        leaves[T1_label] = i['Root1_node'].count(',') + 1
        leaves[T2_label] = i['Root2_node'].count(',') + 1
        B.add_edge(T1_label, T2_label, score = i['Score'])
    elif i.get('Root2_label_tmp', False):
        T1_label = 'T2_T2_' + str(i['Root2_label_tmp'])
        T2_label = 'T2_T2_' + str(i['Root2_label'])
        Y.append(T1_label)
        Y.append(T2_label)
        leaves[T1_label] = i['Root1_node'].count(',') + 1
        leaves[T2_label] = i['Root2_node'].count(',') + 1
    else:
        print(i)
    B.add_edge(T1_label, T2_label, score = i['Score'])
    
X = list(set(X))
Y = list(set(Y))  
node_colorX = [B.degree(v) for v in X]
node_colorY = [B.degree(v) for v in Y]
node_sizeX = [(20 * leaves[i]) for i in X]
node_sizeY = [(20 * leaves[i]) for i in Y]
node_colorX = [-(20 * leaves[i]) for i in X] # more dark color,less size of subtree 
node_colorY = [-(20 * leaves[i]) for i in Y]
# node_sizeX = [B.degree(v) for v in X]
# node_sizeY = [B.degree(v) for v in Y]
width = 1
edge_color = 'black'

#最多连接的前十位
#sorted(((node,degree) for node,degree in B.degree), key= lambda d:d[1],reverse = True)[:10]
#print(X)
#B.nodes()
# print(node_sizeX, node_sizeY )

Edges = []
pos = dict()
pos.update( (n, (1, i)) for i, n in enumerate(X) )
pos.update( (n, (2, i+0.5)) for i, n in enumerate(Y) )
mz = len(X)+len(Y)
plt.figure(figsize=(min(mz, 50), min(2*mz, 100)))

nx.draw_networkx_nodes(B, pos, nodelist=X, node_color=node_colorX,alpha=0.95, node_size = node_sizeX, cmap = plt.cm.Reds)
nx.draw_networkx_nodes(B, pos, nodelist=Y, node_color=node_colorY,alpha=0.95, node_size = node_sizeY, cmap = plt.cm.Reds,)
#nx.draw_networkx_labels(B,pos, font_size=1, font_color='yellow')
#colors = [ B.edges[u,i]['score'] for u,i in B.edges]
edges = nx.draw_networkx_edges(B, pos = pos, edge_color = edge_color,
        width=width, edge_cmap=plt.cm.Blues, edge_vmin = 0, alpha=0.9)
plt.savefig(folder + "Line.pdf")
plt.close()


# In[ ]:


plt.figure(figsize=(min(mz, 50), min(mz, 50)))
nx.draw_networkx(B,
                 pos=nx.circular_layout(B),
                 nodelist = X, 
                 node_color = node_colorX, 
                 node_size = node_sizeX,   
                 node_shape= 'o', 
                 width = width, 
                 # linewidths = 2,
                 with_labels = False, 
                 cmap = plt.cm.Reds,
                 edge_color = edge_color
                )
nx.draw_networkx(B,
                 pos=nx.circular_layout(B), #nx.shell_layout(B), 反过来 
                 nodelist = Y, 
                 node_color = node_colorY, 
                 node_size = node_sizeY,  
                 node_shape= 'v', 
                 width = width, 
                 with_labels = False, 
                 cmap = plt.cm.Reds,
                 edge_color = edge_color
                )
plt.savefig(folder + "Circle.pdf")
plt.close()

plt.figure(figsize=(min(mz, 50), min(mz, 50)))
nx.draw_networkx(B,
                 pos=nx.spring_layout(B),
                 nodelist = X, 
                 node_color = node_colorX, 
                 node_size = node_sizeX, 
                 node_shape= 'o', 
                 width = width, 
                 # linewidths = 2,
                 with_labels = False, 
                 cmap = plt.cm.Reds,
                 edge_color = edge_color
                )
nx.draw_networkx(B,
                 pos=nx.spring_layout(B),
                 nodelist = Y, 
                 node_color = node_colorY, 
                 node_size = node_sizeY,  
                 node_shape= 'v', 
                 width = width, 
                 with_labels = False, 
                 cmap = plt.cm.Reds,
                 edge_color = edge_color
                )
plt.savefig(folder + "Cluster.pdf")
plt.close()

print('network ok!!!')

