import numpy as np
from munkres import Munkres
from math import log
import copy


class Greedy_Algorithm:
    def __init__(self, Local_matrix,
                 local_matrix_root1_index=[],
                 local_matrix_root2_index=[],
                 prune=-1):
        self.Local_matrix = Local_matrix
        self.local_matrix_root1_index = local_matrix_root1_index
        self.local_matrix_root2_index = local_matrix_root2_index
        self.prune = prune

    def calculate(self):
        mat_tmp = self.Local_matrix
        sum = 0
        for i in range(min(self.Local_matrix.shape)):
            # print(mat_tmp)
            sum += np.max(mat_tmp)
            del_i_index = np.where(mat_tmp == np.max(mat_tmp))[0][0]
            del_j_index = np.where(mat_tmp == np.max(mat_tmp))[1][0]
            mat_tmp = np.delete(mat_tmp, del_i_index, 0)
            mat_tmp = np.delete(mat_tmp, del_j_index, 1)
        return sum

    def trace_back(self):
        mat_tmp = self.Local_matrix
        trace = []
        for i in range(min(self.Local_matrix.shape)):
            if np.max(mat_tmp) <= self.prune:
                break
            del_i_index = np.where(mat_tmp == np.max(mat_tmp))[0][0]
            del_j_index = np.where(mat_tmp == np.max(mat_tmp))[1][0]
            trace.append([self.local_matrix_root1_index[del_i_index],
                         self.local_matrix_root2_index[del_j_index]])
            mat_tmp[del_i_index, :] = -99999
            mat_tmp[:, del_j_index] = -99999
        return trace


'''
def topscore(mat, root1addr, root2addr, tracemat, num:int=-1):
    mat_tmp = copy.deepcopy(mat.values)
    scorelist=[]
    if num == -1:
        loopt = min(5,min(mat_tmp.shape))
    else:
        loopt = min(min(mat_tmp.shape), num)
    for _ in range(loopt):
        #print(mat_tmp)      
        del_i_index = np.where(mat_tmp==np.max(mat_tmp))[0][0]
        del_j_index = np.where(mat_tmp==np.max(mat_tmp))[1][0]
        #scorelist.append([np.max(mat_tmp),mat.index[del_i_index],mat.columns[del_j_index]])
        #scorelist.append([np.max(mat_tmp), del_i_index, del_j_index])
        scorelist.append({'Score':np.max(mat_tmp), 
                          'Root1_label':root1addr[del_i_index].label, 
                          'Root1_node':root1addr[del_i_index].nodeobj,
                          'Root2_label':root2addr[del_j_index].label, 
                          'Root2_node':root2addr[del_j_index].nodeobj, 
                          'row':del_i_index, 
                          'col':del_j_index})
        mat_tmp[del_i_index] = -999.
        mat_tmp[:,del_j_index] = -999.
    return scorelist
'''


def GetMaxScore(trace,
                root1,
                root2,
                allmatrix,
                root1_index: int,
                root2_index: int,
                local_matrix,
                local_matrix_root1_index,
                local_matrix_root2_index,
                lll_label,
                llll_label,
                prune=-1.,
                dict_score={},
                Algorithm: str = 'KM',
                merge: float = 100.
                ):
    # 如果两棵树都是，即是根节点也是叶节点（单一元素），直接查字典，查不到就划归为罚分值prune
    root1node = root1
    root2node = root2
    score = prune
    key = root1node.nodeobj+'_'+root2node.nodeobj
    key2 = root2node.nodeobj+'_'+root1node.nodeobj
    if root1node.left == None and root2node.left == None:
        # 路径只能是节点直接匹配
        #if key in dict_score:
            if float(dict_score[key]) > prune:
                trace[root1_index][root2_index].append([root1.nodeobj, root2.nodeobj])
            return float(dict_score[key])
        #else:
            #dict_score[key] = score
            #dict_score[key2] = score
            #return score
    # 如果只有一个树表示为根节点
    elif root1node.left == None and root2node.left is not None:
        # if key in dict_score:
        #    return dict_score[key]
        if (dict_score.get(root1node.nodeobj+'_'+root1node.nodeobj, prune)+prune*(root2node.leaf_count()-1) < prune):
            #dict_score[key] = score
            #dict_score[key2] = score
            return score
        else:
            dearr = local_matrix.flatten()
            maxnum = np.max(dearr)
            score = max(score, maxnum+(root2.leaf_count()-1)*prune)
            if score > prune:
                trace[root1_index][root2_index].append(
                    [root1_index, local_matrix_root2_index[dearr.tolist().index(maxnum)]])
            # root2node = root2node.left #来到子树下
            #key_tmp = root1node.nodeobj+'_'+root2node.nodeobj
            #score = max(score, dict_score[key_tmp])
            # while root2node.right: #子树还有节点
            #    root2node = root2node.right
            #    key_tmp = root1node.nodeobj+'_'+root2node.nodeobj
            #    score = max(score, dict_score[key_tmp])
            #score = max(score+((root2.leaf_count()-1)*prune),prune)
            #dict_score[key] = score
            #dict_score[key2] = score
            # return root2.leaf_count()
            return score

    elif root1node.left is not None and root2node.left == None:
        # if key in dict_score:
        #    return dict_score[key]
        # print(max(local_matrix))
        if (dict_score.get(root2node.nodeobj+'_'+root2node.nodeobj, -1.)+prune*(root1node.leaf_count()-1) < prune):
            #dict_score[key] = score
            #dict_score[key2] = score
            return score
        else:
            dearr = local_matrix.flatten()
            maxnum = np.max(dearr)
            score = max(score, maxnum+(root1.leaf_count()-1)*prune)
            if score > prune:
                trace[root1_index][root2_index].append(
                    [local_matrix_root1_index[dearr.tolist().index(maxnum)], root2_index])
            # root1node = root1node.left #来到子树下
            #key_tmp = root1node.nodeobj+'_'+root2node.nodeobj
            #score = max(score, dict_score[key_tmp])
            # while root1node.right: #子树还有节点
            #    root1node = root1node.right
            #    key_tmp = root1node.nodeobj+'_'+root2node.nodeobj
            #    score = max(score, dict_score[key_tmp])
            #score = max(score+((root1.leaf_count()-1)*prune),prune)
            #dict_score[key] = score
            #dict_score[key2] = score
            # return root2.leaf_count()
            return score
    # 两个树都是非叶子结点
    else:
        summ = 0.
        r1leaf = root1node.leaf_count()
        r2leaf = root2node.leaf_count()
        if r1leaf >= r2leaf:
            for i in root2node.leaves([]):
                summ += dict_score.get(i.nodeobj+'_'+i.nodeobj, prune)
            if (summ + prune*(r1leaf-r2leaf) < r2leaf*prune):
                return r2leaf*prune
        else:
            for i in root1node.leaves([]):
                summ += dict_score.get(i.nodeobj+'_'+i.nodeobj, prune)
            if (summ + prune*(r2leaf-r1leaf) < r1leaf*prune):
                return r1leaf*prune

        if(Algorithm == 'GA'):
            ga = Greedy_Algorithm(
                local_matrix, local_matrix_root1_index, local_matrix_root2_index, prune)
            score = ga.calculate()+abs(r1leaf-r2leaf)*prune
            if score > min(r1leaf, r2leaf)*prune:
                trace[root1_index][root2_index] = ga.trace_back()
            else:
                score = min(r1leaf, r2leaf)*prune
        elif(Algorithm == 'GAR'):
            pass

        else:  # 默认为KM计算
            # github包版本
            trace_tmp = []
            cost_matrix = []
            for row in local_matrix:
                # print(row)
                cost_row = []
                for col in row:
                    cost_row += [99999 - col]
                cost_matrix += [cost_row]
            m = Munkres()
            indexes = m.compute(cost_matrix)
            # print(indexes)
            #print_matrix(local_matrix, msg='Highest profit through this local_matrix:')
            total = 0
            for row, column in indexes:
                #print(row, column)
                value = local_matrix[row][column]
                total += value
                if value > prune:
                    trace_tmp.append(
                        [local_matrix_root1_index[row], local_matrix_root2_index[column]])
                #print(f'({row}, {column}) -> {value}')
            #print('total profit= ',total)
            score = total + \
                abs(local_matrix.shape[0] - local_matrix.shape[1])*prune
            if score > min(r1leaf, r2leaf)*prune:
                trace[root1_index][root2_index] = trace_tmp
            else:
                score = min(r1leaf, r2leaf)*prune

        # 遍历
        for i in zip(root1node.son(), local_matrix_root1_index):
            if i[0].left is not None:
                score_tmp = allmatrix[i[1], root2_index] + \
                    prune*(r1leaf - i[0].leaf_count())
                if score < score_tmp:
                    score = score_tmp
                    trace[root1_index][root2_index].clear()
                    trace[root1_index][root2_index].append([i[1], root2_index])
                #score = max(allmatrix[i[1],root2_index] + prune*(r1leaf - i[0].leaf_count()), score)
        for j in zip(root2node.son(), local_matrix_root2_index):
            if j[0].left is not None:
                score_tmp = allmatrix[root1_index, j[1]] + \
                    prune*(r2leaf - j[0].leaf_count())
                if score < score_tmp:
                    score = score_tmp
                    trace[root1_index][root2_index].clear()
                    trace[root1_index][root2_index].append([root1_index, j[1]])
                #score = max(allmatrix[root1_index,j[1]] + prune*(r2leaf - j[0].leaf_count()), score)

        # print()
        if (abs(merge - (100)) < 1e-10):
            #print('不merge')
            pass
        else:
            # print("merge:",merge)
            #print(root1.label, root1.node_count(), root2.label, root2.node_count())
            #print(root1.node_count() - root1.leaf_count() - 1, root2.node_count() - root2.leaf_count() - 1)
            #print(root1.nodeobj, root2.nodeobj)
            root1_leaves_nodeobj = []
            root2_leaves_nodeobj = []
            root1_leaves_label = []
            root2_leaves_label = []
            for i in root1.leaves([]):
                root1_leaves_nodeobj.append(i.nodeobj)
                root1_leaves_label.append(i.label)
                # print('Leav：',i.label)
            for j in root2.leaves([]):
                root2_leaves_nodeobj.append(j.nodeobj)
                root2_leaves_label.append(j.label)
                # print('Leav：',j.label)
            # print(root1_leaves_label_nodeobj)
            # print(root1_leaves_label)
            # print(root2_leaves_label_nodeobj)

            a = [x for x in root1_leaves_nodeobj if x in root2_leaves_nodeobj]
            b = [root1_leaves_label[x]
                 for x, xi in enumerate(root1_leaves_nodeobj) if xi in a]
            c = [root2_leaves_label[x]
                 for x, xi in enumerate(root2_leaves_nodeobj) if xi in a]
            # print(a)
            # print(b)
            # print(c)

            score_tmpp = ((root1.node_count() - root1.leaf_count() - 1) + (root2.node_count() -
                          root2.leaf_count() - 1) + abs(root1.leaf_count() - root2.leaf_count())) * merge
            #print(abs(root1.leaf_count() - root2.leaf_count()))
            for x in a:
                score_tmpp += float(dict_score[x+'_'+x])

            if score_tmpp > score:
                score = score_tmpp
                trace[root1_index][root2_index].clear()
                for i in zip(b, c):
                    trace[root1_index][root2_index].append(
                        [lll_label.index(i[0]), llll_label.index(i[1])])
                    #print(lll_label.index(i[0]), llll_label.index(i[1]))
            #print(lll_label, llll_label)

        return score
