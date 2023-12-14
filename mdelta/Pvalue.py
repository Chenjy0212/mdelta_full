from copy import *
import random
from multiprocessing import *
import multiprocessing as mp
from .mymath import *
import psutil
from .ReadFile import *
from .auxi import *
import pandas as pd


class MultiTree:
    def __init__(self, nodeobj: str, level: 'int' = 0, label: str = "root"):
        self.nodeobj = nodeobj if nodeobj[-1] != ';' else nodeobj[:-1]
        self.left = None
        self.right = None
        self.level = level
        self.label = label
    # 鍓嶅簭閬嶅巻,鑾峰彇鎸夌収灞傜骇鎺掑簭鐨勫叏鑺傜偣搴忓垪瀵? [node:level]

    def nodes(self, ll={}):
        lll = ll
        if self.nodeobj is not None:
            lll[self] = self.level
            #print(self.nodeobj,' Level:',self.level,' Label:',self.label)
        if self.left is not None:
            self.left.nodes(lll)
        if self.right is not None:
            self.right.nodes(lll)

        lll = sorted(lll.items(), key=lambda item: item[1])
        llist = [[] for i in range(lll[-1][1]+1)]
        for i in lll:
            llist[i[1]].append(i[0])
        return llist
    # 鍚庡簭閬嶅巻_寰楀埌nodeobj:level瀵?

    def postorder(self):
        if self.left is not None:
            self.left.postorder()
        if self.right is not None:
            self.right.postorder()
        if self.nodeobj is not None:
            print(self.nodeobj, ' Level:', self.level, ' Label:', self.label)
    # 鍓嶅簭閬嶅巻_寰楀埌nodeobj:level瀵?

    def preorder(self):
        if self.nodeobj is not None:
            print(self.nodeobj, ' Level:', self.level, ' Label:', self.label)
        if self.left is not None:
            self.left.preorder()
        if self.right is not None:
            self.right.preorder()
    # 灞傚簭閬嶅巻

    def levelorder(self):
        # 杩斿洖鏌愪釜鑺傜偣鐨勫乏瀛╁瓙
        def LChild_Of_Node(node):
            return node.left if node.left is not None else None
        # 杩斿洖鏌愪釜鑺傜偣鐨勫彸瀛╁瓙

        def RChild_Of_Node(node):
            return node.right if node.right is not None else None
        # 灞傚簭閬嶅巻鍒楄〃
        level_order = []
        # 鏄惁娣诲姞鏍硅妭鐐逛腑鐨勬暟鎹?
        if self.nodeobj is not None:
            level_order.append([self])
        # 浜屽弶鏍戠殑楂樺害
        height = self.height()
        if height >= 1:
            # 瀵圭浜屽眰鍙婂叾浠ュ悗鐨勫眰鏁拌繘琛屾搷浣?, 鍦╨evel_order涓坊鍔犺妭鐐硅€屼笉鏄暟鎹?
            for _ in range(2, height + 1):
                level = []  # 璇ュ眰鐨勮妭鐐?
                for node in level_order[-1]:
                    # 濡傛灉宸﹀瀛愰潪绌猴紝鍒欐坊鍔犲乏瀛╁瓙
                    if LChild_Of_Node(node):
                        level.append(LChild_Of_Node(node))
                    # 濡傛灉鍙冲瀛愰潪绌猴紝鍒欐坊鍔犲彸瀛╁瓙
                    if RChild_Of_Node(node):
                        level.append(RChild_Of_Node(node))
                # 濡傛灉璇ュ眰闈炵┖锛屽垯娣诲姞璇ュ眰
                if level:
                    level_order.append(level)

            # 鍙栧嚭姣忓眰涓殑鏁版嵁
            for i in range(0, height):  # 灞傛暟
                for index in range(len(level_order[i])):
                    level_order[i][index] = level_order[i][index].nodeobj
        return level_order
    # 浜屽弶鏍戠殑楂樺害

    def height(self):
        # 绌虹殑鏍戦珮搴︿负0, 鍙湁root鑺傜偣鐨勬爲楂樺害涓?1
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
    # 鍘熸潵鐨勫鍙夋爲鐨勫彾瀛愯妭鐐?

    def leaves(self, ll=[], flag=1):  # 意味 0 这个是叶子节点
        if self is None:
            return None
        if flag == 1:
            if self.left is None:
                ll.append(self)
                return ll
            else:
                self = self.left

        if self.left is not None:
            self.left.leaves(ll, flag=0)
        if self.nodeobj is not None:
            if self.left is None:
                ll.append(self)
        if self.right is not None:
            self.right.leaves(ll, flag=0)
        return ll

    def leaves_nodeobj(self, ll=[], flag=1):  # 意味 0 这个是叶子节点
        if self is None:
            return None
        if flag == 1:
            if self.left is None:
                ll.append(self.nodeobj)
                return ll
            else:
                self = self.left

        if self.left is not None:
            self.left.leaves_nodeobj(ll, flag=0)
        if self.nodeobj is not None:
            if self.left is None:
                ll.append(self.nodeobj)
        if self.right is not None:
            self.right.leaves_nodeobj(ll, flag=0)
        return ll

    def leaves_label(self, ll=[], flag=1):  # 意味 0 这个是叶子节点
        if self is None:
            return None
        if flag == 1:
            if self.left is None:
                ll.append(self.label)
                return ll
            else:
                self = self.left

        if self.left is not None:
            self.left.leaves_label(ll, flag=0)
        if self.nodeobj is not None:
            if self.left is None:
                ll.append(self.label)
        if self.right is not None:
            self.right.leaves_label(ll, flag=0)
        return ll

    def CreatTree(self):
        if(self.nodeobj[0] == '('):  # 瀛樺湪鎷彿鎰忓懗鐫€杩樻病杈惧埌鍙跺瓙缁撶偣
            node_list = []
            brackets_num = 0  # 鎷彿涓暟
            node_num = 0  # 鑺傜偣搴忓彿
            node_tmp = ''
            for i in self.nodeobj[1:-1]:
                if i == '(':
                    brackets_num += 1
                    node_tmp += i
                elif i == ')':
                    brackets_num -= 1
                    node_tmp += i
                elif i == ',' and brackets_num == 0:
                    node_list.append(node_tmp)
                    node_num += 1
                    node_tmp = ''
                else:
                    node_tmp += i
            node_list.append(node_tmp)
            # print(node_list)    #鏌ョ湅鎷嗗垎鍚庣殑鑺傜偣搴忓垪
            self_tmp = self
            for index, item in enumerate(node_list):
                if index == 0:
                    # 绗竴涓瓙鑺傜偣閮芥槸鐖惰妭鐐圭殑宸﹁妭鐐?
                    # 鍚庣画瀛愯妭鐐瑰氨鏄笂涓€涓瓙鑺傜偣鐨勫彸鑺傜偣
                    # 璧嬩簣label
                    label = str(
                        index) if self.label == 'root' else self_tmp.label+','+str(index)
                    self.left = MultiTree(item, label=label)
                    self = self.left
                    self.CreatTree()
                else:
                    label = str(
                        index) if not ',' in self.label else self_tmp.label+','+str(index)
                    self.right = MultiTree(item, label=label)
                    self = self.right
                    self.CreatTree()
    # 璧嬪眰宸﹀瓙鏍戝悗搴忛亶鍘嗘壘鍑烘渶澶у眰绾

    def Posorder_Max_Level(self):
        level = self.level
        while self.right:
            # 鏈夊彸瀛愭爲灏辨槸鏈夊厔寮熺粨鐐?
            level = max(level, self.right.level)
            self = self.right
        return level
    # 璧嬩簣灞傜骇level

    def Level(self):
        if self.left == None:  # 鍙跺瓙鑺傜偣
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
    # 鍙跺瓙鑺傜偣涓暟,闇€瑕佺殑鏄妭鐐逛笅鐨勫乏鑺傜偣鎵嶆纭?

    def leaf_count(self, flag=1):
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
            return self.left.leaf_count(flag=0) + self.right.leaf_count(flag=0)
        elif self.left is None and self.right is not None:
            return 1 + self.right.leaf_count(flag=0)
        elif self.right is None and self.left is not None:
            return self.left.leaf_count(flag=0)

    # 总节点个数
    def node_count(self, flag=0):
        if self is None:
            return 0
        elif self.left is None and self.right is None:
            return 1
        elif self.left is not None and self.right is not None:
            if flag == 1:
                return 1 + self.left.node_count(1) + self.right.node_count()
            else:
                return 1 + self.left.node_count()
        elif self.left is None and self.right is not None:
            return 1 + self.right.node_count()
        elif self.right is None and self.left is not None:
            return 1 + self.left.node_count(1)
        # else:
            # return self.left.node_count()+ self.right.node_count()

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
    # 鑾峰彇褰撳墠鑺傜偣涓€绾у瓙鑺傜偣self

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
    # 涓簭閬嶅巻锛屽緱鍒扮敓鎴愭爲鏁版嵁

    def inorder(self, data=[]):
        if self.left is not None:
            self.left.inorder(data)
        if self.nodeobj is not None:
            data.append(self)
        if self.right is not None:
            self.right.inorder(data)
        return data


def FindNode(Seq: str, times: int) -> list:
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
    return result


class OP:
    def __init__(self,
                 seq1_list_result,
                 seq2_list_result,
                 ScoreDictFile,
                 poolnum=1,
                 mav: float = 2.,
                 miv: float = -1.,
                 pv: float = -1.,
                 notebook: int = 0,
                 Tqdm: int = 1,
                 merge: float = 10,
                 ):
        # 直接调用 Manager 提供的 list() 和 dict()
        self.manager = mp.Manager
        self.mp_lst = self.manager().list()
        self.Seq1_list = seq1_list_result
        self.Seq2_list = seq2_list_result
        self.scoredictfile = ScoreDictFile
        self.mav = mav
        self.miv = miv
        self.pv = pv
        self.Tqdm = Tqdm
        self.notebook = notebook
        self.poolnum = min(poolnum, max(psutil.cpu_count(False), 1))
        self.length = len(seq1_list_result)
        self.merge = merge

    def Foo(self, i, j, scoredictfile):
        root1_tmp = MultiTree(i)
        root1_tmp.CreatTree()
        root1_tmp.Postorder_Level()
        lll_tmp = root1_tmp.nodes({})  # 二维表示
        lllnode_tmp = [j for i in lll_tmp for j in i]
        lllloop_tmp = []
        for i in lll_tmp:
            lllloop_tmp.append(len(i))
        llldict_tmp = {}
        for index, iter in enumerate(lllnode_tmp):
            llldict_tmp[index] = [lllnode_tmp.index(i) for i in iter.son()]

        root2_tmp = MultiTree(j)
        root2_tmp.CreatTree()
        root2_tmp.Postorder_Level()
        llll_tmp = root2_tmp.nodes({})  # 二维表示
        llllnode_tmp = [j for i in llll_tmp for j in i]
        llllloop_tmp = []
        for i in llll_tmp:
            llllloop_tmp.append(len(i))
        lllldict_tmp = {}
        for index, iter in enumerate(llllnode_tmp):
            lllldict_tmp[index] = [llllnode_tmp.index(i) for i in iter.son()]

        score_dict = {}
        if scoredictfile == '':
            score_dict = Scoredict(root1_tmp.leaves(
                []), root2_tmp.leaves([]), self.mav, self.miv)
        else:
            score_dict = QuantitativeScoreFile(root1_tmp.leaves(
                []), root2_tmp.leaves([]), self.mav, self.miv, scoredictfile)

        mmatrix = pd.DataFrame([[0.0 for i in range(len(llllnode_tmp))] for j in range(len(lllnode_tmp))],
                               index=[i.nodeobj for i in lllnode_tmp],
                               columns=[j.nodeobj for j in llllnode_tmp])
        mmatrix.index.name = 'Root1'
        mmatrix.columns.name = 'Root2'
        matrix_values = mmatrix.values

        ttrace = pd.DataFrame([[[] for i in range(len(llllnode_tmp))] for j in range(len(lllnode_tmp))],
                              index=[i.label for i in lllnode_tmp],
                              columns=[j.label for j in llllnode_tmp])
        ttrace.index.name = 'Root1'
        ttrace.columns.name = 'Root2'
        trace_value = ttrace.values
        # 执行循有规律就是(0,0)(0,1)(1,0)(1,1)(0,2)(2,0)(1,2)(2,1)(2,2)(0,3)(3,0)...(n,m)
        for loop_index in loopindex(root1_tmp.level+1, root2_tmp.level+1):
            # print(loop_index)
            for i in range(lllloop_tmp[loop_index[0]]):
                for j in range(llllloop_tmp[loop_index[1]]):
                    i_index = 0
                    j_index = 0
                    for i_tmp in range(loop_index[0]):
                        i_index += lllloop_tmp[i_tmp]
                    i_index += i
                    for j_tmp in range(loop_index[1]):
                        j_index += llllloop_tmp[j_tmp]
                    j_index += j
                    matrix_tmp = matrix_values[llldict_tmp[i_index], :]
                    matrix_tmp = matrix_tmp[:, lllldict_tmp[j_index]]

                    matrix_values[i_index][j_index] = GetMaxScore(trace=trace_value,
                                                                  root1=lll_tmp[loop_index[0]][i],
                                                                  root2=llll_tmp[loop_index[1]][j],
                                                                  allmatrix=matrix_values,
                                                                  root1_index=i_index,
                                                                  root2_index=j_index,
                                                                  local_matrix=matrix_tmp,
                                                                  local_matrix_root1_index=llldict_tmp[i_index],
                                                                  local_matrix_root2_index=lllldict_tmp[j_index],
                                                                  dict_score=score_dict,
                                                                  prune=self.pv,
                                                                  Algorithm='',
                                                                  lll_label=[
                                                                      i.label for i in lll_tmp[0]],
                                                                  llll_label=[
                                                                      i.label for i in llll_tmp[0]],
                                                                  merge=self.merge)
        # print(mmatrix)
        # print(ttrace)
        # print(matrix_values[len(lllnode)-1][len(llllnode)-1])
        self.mp_lst.append(matrix_values[-1][-1])

    def flow(self):
        pool = mp.Pool(self.poolnum)
        if self.Tqdm == 1:
            if self.notebook == 1:
                from tqdm.notebook import tqdm
            else:
                from tqdm import tqdm
            pbar = tqdm(total=self.length)
            pbar.set_description(' Pvalue ')
            update = lambda *args: pbar.update()
        else:
            update = None

        for i, j in zip(self.Seq1_list, self.Seq2_list):
            pool.apply_async(func=self.Foo, args=(
                i, j, self.scoredictfile), callback=update)

        pool.close()
        pool.join()


def pvalue(times: int,
           topscorelist,
           ScoreDictFile: str = '',
           CPUs: int = 50,
           mav: float = 2.,
           miv: float = -1.,
           pv: float = -1.,
           notebook: int = 0,
           Tqdm: int = 1,
           merge: float = 10,
           ):
    Seq1_list_result_max = []
    Seq2_list_result_max = []
    for i in topscorelist:
        Seq1_list_result_max.append(FindNode(i['Root1_node'], times))
        Seq2_list_result_max.append(FindNode(i['Root2_node'], times))

    score_list_max = []
    # print(len(topscorelist))
    for i in range(len(topscorelist)):
        op_max = OP(
            Seq1_list_result_max[i], Seq2_list_result_max[i], ScoreDictFile, CPUs, mav, miv, pv, notebook, Tqdm, merge)
        op_max.flow()
        score_list_max.append(op_max.mp_lst + [topscorelist[i]['Score']])

    # score_list_max.append(float(MaxScore))
    return(score_list_max)
    '''
    top_rate = []
    for i in score_list_max:
        per_less_max = 0
        for j in i:
            if j <= i[-1]:
                per_less_max += 1
        top_rate.append(per_less_max*100/len(i))
    return(top_rate)
    '''
