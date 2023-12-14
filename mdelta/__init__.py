# -*- coding: utf-8 -*- 

from .ReadFile import *
from .mytest import *
from .mymath import *
from .auxi import *
from .Pvalue import *
from .calpvalue import *
import pandas as pd
from copy import *


class MultiTree:
    def __init__(self, nodeobj:str, level:'int'=0, label:str="root"):
        self.nodeobj = nodeobj if nodeobj[-1] != ';' else nodeobj[:-1]
        self.left = None
        self.right = None
        self.level = level
        self.label = label
    # 前序遍历,获取按照层级排序的全节点序列对 [node:level]
    def nodes(self,ll={}):
        lll= ll 
        if self.nodeobj is not None:
            lll[self]=self.level
            #print(self.nodeobj,' Level:',self.level,' Label:',self.label)
        if self.left is not None:
            self.left.nodes(lll)
        if self.right is not None:
            self.right.nodes(lll)
        
        lll = sorted(lll.items(), key=lambda item:item[1])
        llist = [[] for i in range(lll[-1][1]+1)]
        for i in lll: 
            llist[i[1]].append(i[0])
        return llist
    # 后序遍历_得到nodeobj:level对
    def postorder(self):
        if self.left is not None:
            self.left.postorder()
        if self.right is not None:
            self.right.postorder()
        if self.nodeobj is not None:
            print(self.nodeobj,' Level:',self.level,' Label:',self.label)
    # 前序遍历_得到nodeobj:level对
    def preorder(self):
        if self.nodeobj is not None:
            print(self.nodeobj,' Level:',self.level,' Label:',self.label)
        if self.left is not None:
            self.left.preorder()
        if self.right is not None:
            self.right.preorder()
    # 层序遍历
    def levelorder(self):
        # 返回某个节点的左孩子
        def LChild_Of_Node(node):
            return node.left if node.left is not None else None
        # 返回某个节点的右孩子
        def RChild_Of_Node(node):
            return node.right if node.right is not None else None
        # 层序遍历列表
        level_order = []
        # 是否添加根节点中的数据
        if self.nodeobj is not None:
            level_order.append([self])
        # 二叉树的高度
        height = self.height()
        if height >= 1:
            # 对第二层及其以后的层数进行操作, 在level_order中添加节点而不是数据
            for _ in range(2, height + 1):
                level = []  # 璇ュ眰鐨勮妭鐐?
                for node in level_order[-1]:
                    # 如果左孩子非空，则添加左孩子
                    if LChild_Of_Node(node):
                        level.append(LChild_Of_Node(node))
                    # 如果右孩子非空，则添加右孩子
                    if RChild_Of_Node(node):
                        level.append(RChild_Of_Node(node))
                # 如果该层非空，则添加该层
                if level:
                    level_order.append(level)

             # 取出每层中的数据
            for i in range(0, height):  # 层数
                for index in range(len(level_order[i])):
                    level_order[i][index] = level_order[i][index].nodeobj
        return level_order
    # 二叉树的高度
    def height(self):
        # 空的树高度为0, 只有root节点的树高度为1
        if self.nodeobj is None:
            return 0
        elif self.left is None and self.right is None:
            return 1
        elif self.left is None and self.right is not None:
            return 1 + self.right.height()
        elif self.left is not None and self.right is None:
            return 1 + self.left.height()
        else:
            return 1 + max(self.left.height(), self.right.height())
    # 原来的多叉树的叶子节点
    def leaves(self,ll=[],flag = 1): #意味 0 这个是叶子节点
        if self is None:
            return None
        if flag == 1:
            if self.left is None:
                ll.append(self)
                return ll
            else:
                self = self.left
        
        if self.left is not None:
            self.left.leaves(ll,flag = 0)
        if self.nodeobj is not None:
            if self.left is None:
                ll.append(self)
        if self.right is not None:
            self.right.leaves(ll,flag = 0)
        return ll
        
    def leaves_nodeobj(self,ll=[],flag = 1): #意味 0 这个是叶子节点
        if self is None:
            return None
        if flag == 1:
            if self.left is None:
                ll.append(self.nodeobj)
                return ll
            else:
                self = self.left
        
        if self.left is not None:
            self.left.leaves_nodeobj(ll,flag = 0)
        if self.nodeobj is not None:
            if self.left is None:
                ll.append(self.nodeobj)
        if self.right is not None:
            self.right.leaves_nodeobj(ll,flag = 0)
        return ll
    
    def leaves_label(self,ll=[],flag = 1): #意味 0 这个是叶子节点
        if self is None:
            return None
        if flag == 1:
            if self.left is None:
                ll.append(self.label)
                return ll
            else:
                self = self.left
        
        if self.left is not None:
            self.left.leaves_label(ll,flag = 0)
        if self.nodeobj is not None:
            if self.left is None:
                ll.append(self.label)
        if self.right is not None:
            self.right.leaves_label(ll,flag = 0)
        return ll
    #生成树
    def CreatTree(self):
        if(self.nodeobj[0] == '('): #存在括号意味着还没达到叶子结点
            node_list = []
            brackets_num = 0 #括号个数
            node_num = 0 #节点序号
            node_tmp = ''
            for i in self.nodeobj[1:-1]:
                if i == '(':
                    brackets_num +=1
                    node_tmp += i
                elif i == ')':
                    brackets_num -=1
                    node_tmp += i
                elif i == ',' and brackets_num == 0:
                    node_list.append(node_tmp)
                    node_num +=1
                    node_tmp = ''
                else:
                    node_tmp += i
            node_list.append(node_tmp)   
            #print(node_list)    #查看拆分后的节点序列
            self_tmp = self
            for index,item in enumerate(node_list):
                if index == 0:
                #第一个子节点都是父节点的左节点
                #后续子节点就是上一个子节点的右节点 
                #赋予label
                    label = str(index) if self.label == 'root' else self_tmp.label+'_'+str(index)    
                    self.left = MultiTree(item,label=label)
                    self = self.left
                    self.CreatTree()
                else:
                    label = str(index) if not '_' in self.label else self_tmp.label+'_'+str(index)    
                    self.right = MultiTree(item,label=label)
                    self = self.right
                    self.CreatTree()
    #赋层左子树后序遍历找出最大层级d
    def Posorder_Max_Level(self):
        level = self.level;
        while self.right:
            #有右子树就是有兄弟结点
            level = max(level, self.right.level)
            self = self.right
        return level
    #赋予层级level
    def Level(self):
        if self.left == None: #叶子节点
            self.level = 0
        else:
            self.level = self.left.Posorder_Max_Level()+1
    def Postorder_Level(self):
        if self.left is not None:
            self.left.Postorder_Level()
        if self.right is not None:
            self.right.Postorder_Level()
        if self.nodeobj is not None:
            self.Level()
    #叶子节点个数,需要的是节点下的左节点才正确
    def leaf_count(self, flag = 1):
        if self is None:
            return 0
        if flag == 1:
            if self.left is None:
                return 1
            else:
                self = self.left
        if self.left is None and self.right is None:
            return 1
        elif self.right is not None and self.left is not None:
            return self.left.leaf_count(flag = 0) + self.right.leaf_count(flag = 0)
        elif self.left is None and self.right is not None:
            return 1 + self.right.leaf_count(flag = 0)
        elif self.right is None and self.left is not None:
            return self.left.leaf_count(flag = 0)
        
    #总节点个数
    def node_count(self):
        if self is None:
            return 0
        elif self.left is None and self.right is None:
            return 1
        elif self.left is not None and self.right is not None:
            return 1 + self.left.node_count() + self.right.node_count()
        elif self.left is None and self.right is not None:
            return 1 + self.right.node_count()
        elif self.right is None and self.left is not None:
            return 1 + self.left.node_count()
        #else:
            #return self.left.node_count()+ self.right.node_count()
    #获取当前节点一级子节点个数
    def son_count(self):
        self_tmp = self
        son_num = 0
        if self_tmp is None:
            return son_num
        elif self_tmp.left is None:
            return son_num
        elif self_tmp.left is not None:
            son_num += 1
            self_tmp = self_tmp.left
            while self_tmp.right:
                self_tmp = self_tmp.right
                son_num += 1
        return son_num    
    #获取当前节点一级子节点self
    def son(self):
        self_tmp = self
        son_tmp = []
        if self_tmp is None:
            return None
        elif self_tmp.left is None:
            son_tmp.append(self_tmp)
            return son_tmp
        elif self_tmp.left is not None:
            self_tmp = self_tmp.left
            son_tmp.append(self_tmp)
            while self_tmp.right:
                self_tmp = self_tmp.right
                son_tmp.append(self_tmp)
        return son_tmp
    # 中序遍历，得到生成树数据
    def inorder(self,data=[]):
        if self.left is not None:
            self.left.inorder(data)
        if self.nodeobj is not None:
            data.append(self)
        if self.right is not None:
            self.right.inorder(data)
        return data

def label_leaves_list_to_tree(label_list, tree_str):
    str_tmp = ''
    flag = 0
    for i in tree_str:
        if i == '(':
            str_tmp += i
        elif i == ';':
            str_tmp += i
        elif i == ',':
            if str_tmp[-1] == ')':
                str_tmp += ','
            else:
                str_tmp += label_list[flag]
                str_tmp += ','
                flag +=1
        elif i == ')':
            if str_tmp[-1] == ')':
                str_tmp += i
            else:
                str_tmp += label_list[flag]
                str_tmp += i
                flag +=1
    return str_tmp

def scoremat(TreeSeqFile:str, 
            TreeSeqFile2:str, 
            Name2TypeFile:str = '',
            Name2TypeFile2:str = '', 
            ScoreDictFile:str = '', 
            ScoreDictFile1:str = '',
            ScoreDictFile2:str = '',
            top:int = 0,
            pv:float = -1.,
            pp:float = 5,
            mav:float = 2.,
            miv:float = -1.,
            Alg:str = 'KM',
            Tqdm:int = 1,
            merge:float = 10,
            notebook:int = 0,
            diff:int = 0,
            ):
    if Name2TypeFile != '':
        TreeSeqType = ReadTreeSeq_Name2Type(TreeSeqFile,Name2TypeFile)
        TreeSeqOri = ReadTreeSeq(TreeSeqFile)
    else:
        TreeSeqType = ReadTreeSeq(TreeSeqFile)
        TreeSeqOri = ReadTreeSeq(TreeSeqFile)

    if Name2TypeFile2 != '':
        TreeSeqType2 = ReadTreeSeq_Name2Type(TreeSeqFile2,Name2TypeFile2)
        TreeSeqOri2 = ReadTreeSeq(TreeSeqFile2)
    else:
        TreeSeqType2 = ReadTreeSeq(TreeSeqFile2)
        TreeSeqOri2 = ReadTreeSeq(TreeSeqFile2)

    root1 = MultiTree(TreeSeqType)
    root2 = MultiTree(TreeSeqType2)
    oroot1 = MultiTree(TreeSeqOri)
    oroot2 = MultiTree(TreeSeqOri2)

    root1.CreatTree()
    root1.Postorder_Level()
    lll = root1.nodes({}) #二维表示
    lllnode = [j for i in lll for j in i]
    lllnode_obj = [j.nodeobj for i in lll for j in i]
    lllnode_label = [j.label for i in lll for j in i]
    # get root1 leaves' new label to celltype infos
    root1_label2celltype = leafLable_to_celltype_info(lllnode)

    oroot1.CreatTree()
    oroot1.Postorder_Level()
    olll = oroot1.nodes({}) #二维表示
    olllnode = [j for i in olll for j in i]
    oroot1_label2celltype = leafLable_to_celltype_info(olllnode)

    lllloop = []
    for i in lll:
        lllloop.append(len(i))
    llldict = {}
    for index,iter in enumerate(lllnode):
        llldict[index] = [lllnode.index(i) for i in iter.son()]

    root2.CreatTree()
    root2.Postorder_Level()
    llll = root2.nodes({}) #二维表示
    llllnode = [j for i in llll for j in i]
    llllnode_obj = [j.nodeobj for i in llll for j in i]
    llllnode_label = [j.label for i in llll for j in i]
    # get root2 leaves' new label to celltype infos
    root2_label2celltype = leafLable_to_celltype_info(llllnode)

    oroot2.CreatTree()
    oroot2.Postorder_Level()
    ollll = oroot2.nodes({}) #二维表示
    ollllnode = [j for i in ollll for j in i]
    oroot2_label2celltype = leafLable_to_celltype_info(ollllnode)

    llllloop = []
    for i in llll:
        llllloop.append(len(i))
    lllldict = {}
    for index,iter in enumerate(llllnode):
        lllldict[index] = [llllnode.index(i) for i in iter.son()]

    if ScoreDictFile != '':
        score_dict = QuantitativeScoreFile(root1.leaves([]),root2.leaves([]),mav, miv, ScoreDictFile)
    elif ScoreDictFile1 != '' and ScoreDictFile2 != '':
        result_tmp = QuantitativeScoreFile2(root1.leaves([]),root2.leaves([]),mav, miv, ScoreDictFile1, ScoreDictFile2, pv, pp)
        score_dict = result_tmp
        #print(pv)
    else:
        score_dict = Scoredict(root1.leaves([]),root2.leaves([]), mav, miv)

    mmatrix = pd.DataFrame([[0.0 for i in range(len(llllnode))] for j in range(len(lllnode))],
                        index=[i.label for i in lllnode],
                        columns=[j.label for j in llllnode])
    mmatrix.index.name = 'Root1' 
    mmatrix.columns.name = 'Root2'
    matrix_values = mmatrix.values

    ttrace = pd.DataFrame([[[] for i in range(len(llllnode))] for j in range(len(lllnode))],
                        index=[i.label for i in lllnode],
                        columns=[j.label for j in llllnode])
    ttrace.index.name = 'Root1' 
    ttrace.columns.name = 'Root2'
    trace_value = ttrace.values
    if Tqdm == 1:
        if notebook == 1:
            from tqdm.notebook import tqdm
        else:
            from tqdm import tqdm
        with tqdm(total=(root1.node_count()*root2.node_count()),desc='Matrix Node') as pbar:
            #执行循有规律就是(0,0)(0,1)(1,0)(1,1)(0,2)(2,0)(1,2)(2,1)(2,2)(0,3)(3,0)...(n,m)
            for loop_index in loopindex(root1.level+1,root2.level+1):
                #print(loop_index)
                for i in range(lllloop[loop_index[0]]):
                    for j in range(llllloop[loop_index[1]]): 
                        i_index = 0
                        j_index = 0
                        for i_tmp in range(loop_index[0]):
                            i_index += lllloop[i_tmp]
                        i_index+=i
                        for j_tmp in range(loop_index[1]):
                            j_index += llllloop[j_tmp]
                        j_index+=j
                        matrix_tmp = matrix_values[llldict[i_index],:]
                        matrix_tmp = matrix_tmp[:,lllldict[j_index]]

                        matrix_values[i_index][j_index] = GetMaxScore(trace=trace_value,
                                                                        root1=lll[loop_index[0]][i],
                                                                        root2=llll[loop_index[1]][j],
                                                                        allmatrix = matrix_values,
                                                                        root1_index=i_index, 
                                                                        root2_index=j_index, 
                                                                        local_matrix=matrix_tmp,
                                                                        local_matrix_root1_index = llldict[i_index],
                                                                        local_matrix_root2_index = lllldict[j_index], 
                                                                        dict_score=score_dict,
                                                                        Algorithm=Alg,
                                                                        prune=pv,
                                                                        lll_label = [i.label for i in lll[0]],
                                                                        llll_label = [i.label for i in llll[0]],
                                                                        merge = merge,)
                        pbar.update(1)
    else:
        for loop_index in loopindex(root1.level+1,root2.level+1):
                #print(loop_index)
                for i in range(lllloop[loop_index[0]]):
                    for j in range(llllloop[loop_index[1]]): 
                        i_index = 0
                        j_index = 0
                        for i_tmp in range(loop_index[0]):
                            i_index += lllloop[i_tmp]
                        i_index+=i
                        for j_tmp in range(loop_index[1]):
                            j_index += llllloop[j_tmp]
                        j_index+=j
                        matrix_tmp = matrix_values[llldict[i_index],:]
                        matrix_tmp = matrix_tmp[:,lllldict[j_index]]

                        matrix_values[i_index][j_index] = GetMaxScore(trace=trace_value,
                                                                        root1=lll[loop_index[0]][i],
                                                                        root2=llll[loop_index[1]][j],
                                                                        allmatrix = matrix_values,
                                                                        root1_index=i_index, 
                                                                        root2_index=j_index, 
                                                                        local_matrix=matrix_tmp,
                                                                        local_matrix_root1_index = llldict[i_index],
                                                                        local_matrix_root2_index = lllldict[j_index], 
                                                                        dict_score=score_dict,
                                                                        Algorithm=Alg,
                                                                        prune=pv,
                                                                        lll_label = [i.label for i in lll[0]],
                                                                        llll_label = [i.label for i in llll[0]],
                                                                        merge = merge,)
    
    def changemat(rac, tracemat, tracemat_value, mat, list_tmp1, list_tmp2, ttused):
        for i in tracemat_value[tuple(rac)]:
            if isinstance(i[0],int) and isinstance(i[1],int):
                list_tmp1.append(tracemat.index[i[0]])
                #print(tracemat.index[i[0]])
                list_tmp2.append(tracemat.columns[i[1]])
                #print(tracemat.index[i[1]])
            elif not isinstance(i[0],int) and not isinstance(i[1],int):    
                list_tmp1.append(tracemat.index[rac[0]])
                #print(tracemat.index[rac[0]])
                list_tmp2.append(tracemat.columns[rac[1]])
                return
                
            #print(i[0],i[0])
        for i in tracemat_value[tuple(rac)]:
            ttused.append(str(i))
            if isinstance(i[0],int) and isinstance(i[1],int):
                #print(i[0],root1.leaf_count(), i[1],root2.leaf_count())
                if i == [] or (i[0]<root1.leaf_count() and i[1] <root2.leaf_count()):
                    mat[tuple(i)] = -99999.
                else:
                    mat[tuple(i)] = -99999.
                    #print(i)
                    changemat(i, tracemat,tracemat_value, mat, list_tmp1, list_tmp2, ttused)
            #return mat
    def getmatchtree(rac, lllnode_obj, llllnode_obj, tracemat_value, mat, tree_tmp1, tree_tmp2):
        for i in tracemat_value[tuple(rac)]:
            if isinstance(i[0],int) and isinstance(i[1],int):
                if i[0]<root1.leaf_count() and i[1] <root2.leaf_count():
                    mat[tuple(i)] = -99999.
                    if tree_tmp1 == []:
                        tree_tmp1.append('('+str(lllnode_obj[i[0]]) + ')')
                        tree_tmp2.append('('+str(llllnode_obj[i[1]]) + ')')
                    elif tree_tmp1[-1] == '(': #if tree_tmp2[-1] == '(':
                        tree_tmp1.append(str(lllnode_obj[i[0]]))
                        tree_tmp2.append(str(llllnode_obj[i[1]]))
                    else:
                        tree_tmp1.append(','+str(lllnode_obj[i[0]]))
                        tree_tmp2.append(','+str(llllnode_obj[i[1]]))
                else:
                    mat[tuple(i)] = -99999.
                    if tree_tmp1 == []:
                        pass
                    elif tree_tmp1[-1] != '(':
                        tree_tmp1.append(',')
                        tree_tmp2.append(',')
                    tree_tmp1.append('(')
                    tree_tmp2.append('(')
                    getmatchtree(i, lllnode_obj, llllnode_obj, tracemat_value, mat, tree_tmp1, tree_tmp2)
                    tree_tmp1.append(')')
                    tree_tmp2.append(')')
            else:
                tree_tmp1.append(lllnode_obj[rac[0]])
                tree_tmp2.append(llllnode_obj[rac[1]])

    def where_prune(match_list:list, leaves_list:list):
        leaves_list_tmp = deepcopy(leaves_list)
        for i in leaves_list:
            if i in match_list:
                leaves_list_tmp.remove(i)
        return leaves_list_tmp
    
    if top == 0: #默认情况下
        T1root_T2root = []
        
        mat_tmp = deepcopy(matrix_values)
        list_tmp1 = []
        list_tmp2 = []
        mat_tmp[-1,-1] = -99999.
        ttused = [str([root1.node_count()-1, root2.node_count()-1])]
        changemat([-1,-1],ttrace, trace_value,mat_tmp, list_tmp1, list_tmp2, ttused)
        
        mat_tmp2 = deepcopy(matrix_values)
        tree_tmp1 = []
        tree_tmp2 = []
        mat_tmp2[-1,-1] = -99999.
        getmatchtree([-1,-1],lllnode_obj, llllnode_obj, trace_value,mat_tmp2, tree_tmp1, tree_tmp2)
        #tree_tmp1.append(';')
        #tree_tmp2.append(';')

        mat_tmp3 = deepcopy(matrix_values)
        tree_tmp3 = []
        tree_tmp4 = []
        mat_tmp3[-1,-1] = -99999.
        getmatchtree([-1,-1],lllnode_label, llllnode_label, trace_value,mat_tmp3, tree_tmp3, tree_tmp4)
        #tree_tmp3.append(';')
        #tree_tmp4.append(';')
        
        #list_tmp1.insert(0,root1.label)
        #list_tmp2.insert(0,root2.label)
        T1root_T2root.append({'Score':matrix_values[-1][-1],
                            'Root1_label':root1.label + ';', 
                            'Root1_node':root1.nodeobj + ';',
                            'Root1_seq':oroot1.nodeobj + ';',
                            'Root1_label_node': label_leaves_list_to_tree(root1.leaves_label([]), root1.nodeobj) + ';',
                            'Root2_label':root2.label + ';', 
                            'Root2_node':root2.nodeobj + ';', 
                            'Root2_seq':oroot2.nodeobj + ';',
                            'Root2_label_node': label_leaves_list_to_tree(root2.leaves_label([]), root2.nodeobj) + ';',
                            'Root1_match': list_tmp1,
                            'Root2_match': list_tmp2,
                            'Root1_match_tree': '(' + ''.join(tree_tmp1) + ');',
                            'Root2_match_tree': '(' + ''.join(tree_tmp2) + ');',
                            'Root1_match_label_tree': '(' + ''.join(tree_tmp3) + ');',
                            'Root2_match_label_tree': '(' + ''.join(tree_tmp4) + ');',
                            'Root1_prune':where_prune(list_tmp1, list(map(lambda x:x.label,root1.leaves([])))),
                            'Root2_prune':where_prune(list_tmp2, list(map(lambda x:x.label,root2.leaves([])))),
                            'row':root1.node_count()-1, 
                            'col':root2.node_count()-1,
                            'used':ttused
                             })

        return({'matrix':mmatrix, 
                'tree1_leaves_nodename': oroot1_label2celltype[1],
                'tree1_leaves_label': root1_label2celltype[0],
                'tree1_leaves_celltype': root1_label2celltype[1],
                'tree2_leaves_nodename': oroot2_label2celltype[1],
                'tree2_leaves_label': root2_label2celltype[0],
                'tree2_leaves_celltype': root2_label2celltype[1],
                'score_dict':score_dict,
                'T1root_T2root':T1root_T2root})
        
    elif top > 0 and top < min(root1.node_count(), root2.node_count()):

        mat_tmp = deepcopy(matrix_values)
        mat_tmp2 = deepcopy(matrix_values)
        mat_tmp3 = deepcopy(matrix_values)
        scorelist=[]
        for jjj in range(top):
            percent = 100.
            #print(mat_tmp)
            maxscore = np.max(mat_tmp)
            del_i_index = np.where(mat_tmp==np.max(mat_tmp))[0][0]
            del_j_index = np.where(mat_tmp==np.max(mat_tmp))[1][0]
            
            list_tmp1 = []
            list_tmp2 = []
            mat_tmp[del_i_index,del_j_index] = -99999.
            ttused = [str([del_i_index, del_j_index])]
            changemat([del_i_index,del_j_index],ttrace, trace_value,mat_tmp, list_tmp1, list_tmp2, ttused)
            
            tree_tmp1 = []
            tree_tmp2 = []
            mat_tmp2[del_i_index,del_j_index] = -99999.
            getmatchtree([del_i_index,del_j_index],lllnode_obj, llllnode_obj, trace_value,mat_tmp2, tree_tmp1, tree_tmp2)
           
            tree_tmp3 = []
            tree_tmp4 = []
            mat_tmp3[del_i_index,del_j_index] = -99999.
            getmatchtree([del_i_index,del_j_index],lllnode_label, llllnode_label, trace_value,mat_tmp3, tree_tmp3, tree_tmp4)
            
            #if list_tmp1[0] != lllnode[del_i_index].label:
            #    list_tmp1.insert(0,lllnode[del_i_index].label)
            #if list_tmp2[0] != llllnode[del_j_index].label:    
            #    list_tmp2.insert(0,llllnode[del_j_index].label)
            if jjj > 0 and diff > 0:
                for i in range(len(scorelist)):
                    percent = min(percent, len(set(scorelist[i]['used']) - set(ttused)) / len(scorelist[i]['used'])*100.)
                
                while percent < float(diff):
                    maxscore = np.max(mat_tmp)
                    if (maxscore + 99999.) < 0.01:
                        break
                    del_i_index = np.where(mat_tmp==np.max(mat_tmp))[0][0]
                    del_j_index = np.where(mat_tmp==np.max(mat_tmp))[1][0]
                    
                    list_tmp1 = []
                    list_tmp2 = []
                    mat_tmp[del_i_index,del_j_index] = -99999.
                    ttused = [str([del_i_index, del_j_index])]
                    changemat([del_i_index,del_j_index],ttrace, trace_value,mat_tmp, list_tmp1, list_tmp2, ttused)
                    
                    tree_tmp1 = []
                    tree_tmp2 = []
                    mat_tmp2[del_i_index,del_j_index] = -99999.
                    getmatchtree([del_i_index,del_j_index],lllnode_obj, llllnode_obj, trace_value,mat_tmp2, tree_tmp1, tree_tmp2)
                    
                    tree_tmp3 = []
                    tree_tmp4 = []
                    mat_tmp3[del_i_index,del_j_index] = -99999.
                    getmatchtree([del_i_index,del_j_index],lllnode_label, llllnode_label, trace_value,mat_tmp3, tree_tmp3, tree_tmp4)
                    percent = 100
                    for i in range(len(scorelist)):
                        percent = min(percent, len(set(scorelist[i]['used']) - set(ttused)) / len(scorelist[i]['used'])*100.)
                    
            #tree_tmp1.append(';')
            #tree_tmp2.append(';')
            #tree_tmp3.append(';')
            #tree_tmp4.append(';')     
                
            scorelist.append({'Score':maxscore,
                            'Root1_label':lllnode[del_i_index].label + ';', 
                            'Root1_node':lllnode[del_i_index].nodeobj + ';',
                            'Root1_seq':olllnode[del_i_index].nodeobj + ';',
                            'Root1_label_node': label_leaves_list_to_tree(lllnode[del_i_index].leaves_label([]), lllnode[del_i_index].nodeobj) + ';',
                            'Root2_label':llllnode[del_j_index].label + ';', 
                            'Root2_node':llllnode[del_j_index].nodeobj + ';', 
                            'Root2_seq':ollllnode[del_j_index].nodeobj + ';', 
                            'Root2_label_node': label_leaves_list_to_tree(llllnode[del_j_index].leaves_label([]), llllnode[del_j_index].nodeobj) + ';',
                            'Root1_match': list_tmp1,
                            'Root2_match': list_tmp2,
                            'Root1_match_tree': '(' + ''.join(tree_tmp1) + ');',
                            'Root2_match_tree': '(' + ''.join(tree_tmp2) + ');',
                            'Root1_match_label_tree': '(' + ''.join(tree_tmp3) + ');',
                            'Root2_match_label_tree': '(' + ''.join(tree_tmp4) + ');',
                            'Root1_prune':where_prune(list_tmp1, list(map(lambda x:x.label,lllnode[del_i_index].leaves([])))),
                            'Root2_prune':where_prune(list_tmp2, list(map(lambda x:x.label,llllnode[del_j_index].leaves([])))),
                            'row':del_i_index, 
                            'col':del_j_index,
                            'used':ttused
                             })

        return({'matrix':mmatrix, 
                'tree1_leaves_nodename': oroot1_label2celltype[1],
                'tree1_leaves_label': root1_label2celltype[0],
                'tree1_leaves_celltype': root1_label2celltype[1],
                'tree2_leaves_nodename': oroot2_label2celltype[1],
                'tree2_leaves_label': root2_label2celltype[0],
                'tree2_leaves_celltype': root2_label2celltype[1],
                'score_dict':score_dict,
                'TopScoreList':scorelist})
        
    else:
        print("Parameter top cannot be negative, or it might be out of range.")
        print("root1_nodecount:",root1.node_count(), "root2_nodecount:",root2.node_count())