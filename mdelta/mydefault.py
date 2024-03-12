# -- coding: utf-8 --

import ipywidgets as widgets
from ipywidgets import Button, Layout, jslink, IntText, IntSlider
from ipywidgets import TwoByTwoLayout
from ipywidgets import interact, interactive, fixed, interact_manual
from IPython.display import display
from multiprocessing import cpu_count
import os


top_left_text = IntText(description='Top left', layout=Layout(width='auto', height='auto'))
top_right_text = IntText(description='Top right', layout=Layout(width='auto', height='auto'))
bottom_left_slider = IntSlider(description='Bottom left', layout=Layout(width='auto', height='auto'))
bottom_right_slider = IntSlider(description='Bottom right', layout=Layout(width='auto', height='auto'))

top_left_file = widgets.Textarea(layout=Layout(width='auto', height='auto'))
top_right_file = widgets.Textarea(layout=Layout(width='auto', height='auto'))
bottom_left_upload = widgets.Combobox(layout=Layout(width='auto', height='auto'))
bottom_right_upload = widgets.Combobox(layout=Layout(width='auto', height='auto')) 

def f(x):
    return x

def get_default():
    #展示logo
    img1 = open("./image/mlogo.png", "rb")
    logoimage = img1.read()
    logo = widgets.Image(
        value= logoimage,
        format='png',
        width='auto',
        height='auto',
    )
    display(logo)
    
    #匹配流程图
    img2 = open("./image/mDELTA.png", "rb")
    mDELTAimage = img2.read()
    img_mDELTA = widgets.Image(
        value=mDELTAimage,
        format='png',
        width='auto',
        height='auto',
    )
    #中山大学
    img3 = open("./image/mdelta_algorithm.png", "rb")
    algorithmimage = img3.read()
    mdelta_algorithm = widgets.Image(
        value = algorithmimage,
        format='png',
        width='auto',
        height='auto',
    )
    how = widgets.HTML(
        value = """
<h1>Positional arguments</h1>
<table border="1">
<tr>
<th>Parameter</th> 
<th>Type</th> 
<th>Description</th>
</tr>
<td> TreeSeqFile </td>
<td> [path/filename] </td> 
<td> A text file storing cell lineage tree #1 in newick format. Tips can be labeled by name or cell type. Branch lengths should be removed. </td>
</tr>
<tr>
<td> TreeSeqFile2 </td>
<td> [path/filename] </td> 
<td> A text file storing cell lineage tree #2 in newick format. Tips can be labeled by name or cell type. Branch lengths should be removed. </td>
</tr>
</table>


<h1>Optional arguments</h1>
<table border="1">
<tr>
<th>Parameter</th> 
<th>Type</th> 
<th>Description</th>
</tr>
<td> -h, --help </td>
<td></td> 
<td> show this help message and exit </td>
</tr>

<tr>
<td> -nt NAME2TYPEFILE, <br>--Name2TypeFile NAME2TYPEFILE </td>
<td> [path/filename] </td> 
<td> List of correspondance between tip name and cell type for cell lineage tree #1. </td>
</tr>

<tr>
<td> -nt2 NAME2TYPEFILE, <br>--Name2TypeFile2 NAME2TYPEFILE </td>
<td> [path/filename] </td> 
<td> List of correspondance between tip name and cell type for cell lineage tree #2. </td>
</tr>

<tr>
<td> -xsd XSCOREDICTFILE, <br>--ScoreDictFile SCOREDICTFILE </td>
<td> [path/filename] </td> 
<td> A comma-delimited text file used to
                        determine similarity scores between cells. If there
                        are exactly three columns, they will be interpreted as
                        (1) the cell (name or type) in Tree #1, (2) the cell
                        in Tree #2, and (3) the similarity score. If
                        otherwise, the first column will be interpreted as the
                        cell (name or type) and the remaining columns as
                        features of the cell (e.g. expression of a gene). The
                        similarity scores will be estimated between all pairs
                        of cells based on the Euclidean distance calculated
                        using all the features. Overrides `-ma` and `-mi`. </td>
</tr>

<tr>
<td> -lsd LSCOREDICTFILE1, <br>--ScoreDictFile1 SCOREDICTFILE1 </td>
<td> [path/filename] </td> 
<td> Calculate the Euclidean distance based on the characteristics of the cell (such as gene expression) to obtain the result x, and obtain the corresponding score through a series of operations such as - ln (x+1). </td>
</tr>


<tr>
<td> -t TOP, --top TOP </td>
<td> [int > 0] </td> 
<td> Performs local (instead of global)
                        alignment, and output the top NUM local alignments
                        with the highest score (e.g. `-t 10`). In the case of
                        global alignment, this parameter should be omitted. </td>
</tr>

<tr>
<td> -ma MAV, --mav MAV </td>
<td> [float] </td> 
<td> Default=2. Score of two paired end nodes with the same name or type.  </td>
</tr>

<tr>
<td> -mi MIV, --miv MIV </td>
<td> [float] </td> 
<td> Default=-1. Shorthand for a simple matching
                        score scheme, where the matching score between a pair
                        of the same cell types is MAV and all other pairs are
                        MIV. (e.g. `-ma 2 -mi -2`). Overridden by `-sd`. </td>
</tr>

<tr>
<td> -p PV, --pv PV  </td>
<td> [float] </td> 
<td> Default=-1. The score for pruning a tip of the tree (e.g.
                        `-p -2`). Default to -1. </td>
</tr>

<tr>
<td> -T TQDM, <tqdm>--Tqdm TQDM  </td>
<td> [0(off) or 1(on)] </td> 
<td> Whether to display the running Progress bar </td>
</tr>

<tr>
<td> -n NOTEBOOK, <br>--notebook NOTEBOOK  </td>
<td> [0(off) or 1(on)] </td> 
<td> Toggle for the jupyter notebook
                        environment. </td>
</tr>

<tr>
<td> -P PERM, --PERM PERM  </td>
<td> [int > 0] </td> 
<td> Toggle for the statistical significance. For
                        each observed alignment, the aligned trees will be
                        permuted PERM times to generate a null distribution of
                        alignment scores, with which a P value can be
                        calculated for the observed alignment score. </td>
</tr>

<tr>
<td> -c CPUS, --CPUs CPUS  </td>
<td> [int > 0] </td> 
<td> Number of threads for multi-processing.
                        Default to 50., it can reach the maximum number of
                        local CPU cores - 1. </td>
</tr>

<tr>
<td> -o OUTPUT, <br>--output OUTPUT  </td>
<td> [path] </td> 
<td> Output path, eg:'/home/username' </td>
</tr>

<tr>
<td> -mg MERGE, <br>--merge MERGE  </td>
<td> [float] </td> 
<td> This is the scaling factor for calculating the
                        score of merging an internal node (e.g. -mg -1), which
                        is multiplied by the number of tips of the internal
                        node to be merged. Default to 0. </td>
</tr>

<tr>
<td> -x DIFF, <br>--diff DIFF  </td>
<td> [int > 0] </td> 
<td> Alignment must consist of a minimal of DIFF
                        percent aligned cell pairs that are different from
                        previous(better) local alignments in order to be
                        considered as another new alignment (e.g. `-x 20`
                        means 20 percent). </td>
</tr>

</table>

<br>
<div> 
    More details on :
    <a href='https://github.com/Chenjy0212/mdelta_full' style='display: flex' target="_blank"> 
        <h1 style='color: red'>M </h1>
        <h1 style='color: orange'>D </h1>
        <h1 style='color: green'>E </h1>
        <h1 style='color: lightblue'>L </h1>
        <h1 style='color: blue'>T </h1>
        <h1 style='color: purple'>A </h1>
    </a>
</div>
        """,
        # placeholder='Some HTML',
        # description='Some HTML',
    )
    '''
    positional = TwoByTwoLayout(top_left=top_left_file,
               top_right=top_right_file,
               bottom_left=bottom_left_upload,
               bottom_right=bottom_right_upload)
    TreeSeqFileList = ['/home/ee_while/TreeFile_nwk/Example_CLT_1.nwk',
                      '/home/ee_while/TreeFile_nwk/Example_CLT_2.nwk',
                      '/home/ee_while/TreeFile_nwk/a1_cbrad5_onecell.nwk',
                      '/home/ee_while/TreeFile_nwk/e5_hesc_onecell.nwk',
                      '/home/ee_while/TreeFile_nwk/f11_hesc_onecell.nwk',
                      '/home/ee_while/TreeFile_nwk/g11_cbrad5_onecell.nwk',
                      '/home/ee_while/TreeFile_nwk/a2_cbrad5_onecell.nwk',
                      '/home/ee_while/TreeFile_nwk/gs_hesc_onecell.nwk',]
    #控件的标题
    positional.top_left.description = 'TreeSeqFile'
    positional.bottom_left.description = '注:删除后选择'
    positional.top_right.description='TreeSeqFile2'
    positional.bottom_right.description= '注:删除后选择'
    #空白时候的提示
    positional.top_right.placeholder=\
    positional.bottom_right.placeholder=\
    positional.bottom_left.placeholder=\
    positional.top_left.placeholder='请在此输入文件名或文件路径'
    #默认值
    positional.top_left.value = '/home/ee_while/TreeFile_nwk/Example_CLT_1.nwk'
    positional.top_right.value = '/home/ee_while/TreeFile_nwk/Example_CLT_2.nwk'
    #默认值的选项目列表
    positional.bottom_right.options = TreeSeqFileList
    positional.bottom_left.options = TreeSeqFileList

    jslink((positional.top_left, 'value'), (positional.bottom_left, 'value'))
    jslink((positional.top_right, 'value'), (positional.bottom_right, 'value'))
    '''
    #必选参数
    #TreeSeqFile,TreeSeqFile2
    TreeSeqFile = widgets.Combobox(layout=Layout(width='auto', height='auto'))
    TreeSeqFile2 = widgets.Combobox(layout=Layout(width='auto', height='auto'))
    #控件的标题
    TreeSeqFile.description = 'TreeSeqFile'
    TreeSeqFile2.description = 'TreeSeqFile2'
    #空白时候的提示
    TreeSeqFile.placeholder = '请在此输入文件名或文件路径，或者在最右方点击选择'
    TreeSeqFile2.placeholder = '请在此输入文件名或文件路径，或者在最右方点击选择'
    #默认值
    # TreeSeqFile.value = '/home/ee_while/TreeFile_nwk/Example_CLT_1.nwk'
    # TreeSeqFile2.value = '/home/ee_while/TreeFile_nwk/Example_CLT_2.nwk'
    #默认值的选项目列表
    # TreeSeqFileList = ['/home/ee_while/TreeFile_nwk/Example_CLT_1.nwk',
    #                   '/home/ee_while/TreeFile_nwk/Example_CLT_2.nwk',
    #                   '/home/ee_while/TreeFile_nwk/a1_cbrad5_onecell.nwk',
    #                   '/home/ee_while/TreeFile_nwk/e5_hesc_onecell.nwk',
    #                   '/home/ee_while/TreeFile_nwk/f11_hesc_onecell.nwk',
    #                   '/home/ee_while/TreeFile_nwk/g11_cbrad5_onecell.nwk',
    #                   '/home/ee_while/TreeFile_nwk/g2_cbrad5_onecell.nwk',
    #                   '/home/ee_while/TreeFile_nwk/gs_hesc_onecell.nwk',]
    TreeSeqFileList = ['ExampleFile/tree1.nwk',
                      'ExampleFile/tree2.nwk',
                      ]
    TreeSeqFile.options = TreeSeqFileList
    TreeSeqFile2.options = TreeSeqFileList
    #合在一起
    positional = widgets.VBox([TreeSeqFile, TreeSeqFile2])
    
    #option 可选参数\
    #名称转类型
    Name2TypeFile = widgets.Combobox(layout=Layout(width='auto', height='auto'))
    Name2TypeFile2 = widgets.Combobox(layout=Layout(width='auto', height='auto'))
    #控件的标题
    Name2TypeFile.description = 'Name2TypeFile'
    Name2TypeFile2.description = 'Name2TypeFile2'
    #空白时候的提示
    Name2TypeFile.placeholder = '请在此输入文件名或文件路径，或者在最右方点击选择'
    Name2TypeFile2.placeholder = '请在此输入文件名或文件路径，或者在最右方点击选择'
    #默认值
    #Name2TypeFile.value = "/mnt/data5/disk/ee_while/mdelta/ExampleFile/Name2Type.csv"
    #Name2TypeFile2.value = "/mnt/data5/disk/ee_while/mdelta/ExampleFile/Name2Type.csv"
    #默认值的选项目列表
    # Name2TypeFileist = ['/mnt/data5/disk/ee_while/mdelta/ExampleFile/Name2Type.csv',
    #                   '/home/ee_while/JOB220429/a1_cbrad5/A1_CBRAD5_newNAME2CLASS.csv',
    #                   '/home/ee_while/JOB220429/g11_cbrad5/G11_CBRAD5_newNAME2CLASS.csv',
    #                   '/home/ee_while/JOB220429/g2_cbrad5/G2_CBRAD5_newNAME2CLASS.csv',
    #                   '/home/ee_while/JOB220429/gs_hesc/GS_HESC_newNAME2CLASS.csv',]
    Name2TypeFileist = ['ExampleFile/Name2Type.csv',
                      ]
    Name2TypeFile.options = Name2TypeFileist
    Name2TypeFile2.options = Name2TypeFileist
    #合在一起
    Name2Type = widgets.VBox([Name2TypeFile, Name2TypeFile2])
    
    #定性计算
    XScoreDictFile = widgets.Combobox(layout=Layout(width='auto', height='auto'))
    #控件的标题
    XScoreDictFile.description = 'XScoreDictFile'
    #空白时候的提示
    XScoreDictFile.placeholder = '请在此输入文件名或文件路径，或者在最右方点击选择'
    #默认值的选项目列表
    ScoreDictFileist = ['ExampleFile/Xscorefile.csv',
                       ]
    XScoreDictFile.options = ScoreDictFileist
    #合在一起
    dingxing = widgets.VBox([XScoreDictFile,])
    
    #定量计算
    #默认值的选项目列表
    LScoreDictFileist = ['ExampleFile/Lscorefile.csv']
                        # '/home/ee_while/JOB220901/A1.csv',
                        # '/home/ee_while/JOB220901/G2.csv',
                        # '/home/ee_while/JOB220901/G11.csv',
                        # '/home/ee_while/JOB220901/E5.csv',
                        # '/home/ee_while/JOB220901/F11.csv',
                        # '/home/ee_while/JOB220901/GS.csv',]
    LScoreDictFile = widgets.Combobox(layout=Layout(width='auto', height='auto'), description = 'LScoreDictFile', placeholder = '请在此输入文件名或文件路径，或者在最右方点击选择', options = LScoreDictFileist)
#     QScoreDictFile2 = widgets.Combobox(layout=Layout(width='auto', height='auto'), description = 'ScoreDictFile2', placeholder = '请在此输入文件名或文件路径，或者在最右方点击选择', options = QScoreDictFileist)
    
    #定量计算截取罚分的percent
    # PP1 = widgets.FloatRangeSlider(value = [5.0, 5.0], description='pp', min = 0.0, max= 100.0, layout=Layout(width='auto', height='auto'))
    # PPstep = widgets.FloatSlider(value = 5.0, description='pp_step', min = 0.01, max= 100.0, layout=Layout(width='auto', height='auto'))
    #合在一起
    # dingliang = widgets.VBox([QScoreDictFile, QScoreDictFile2, widgets.VBox([PP1, PPstep])])
    dingliang = widgets.VBox([LScoreDictFile])
    
    #top
    top1 = widgets.IntSlider(description='top', min = 0, max= 2000, layout=Layout(width='auto', height='auto'))
    top2 = widgets.IntText(layout=Layout(width='auto', height='auto'))
    top1.value = 0
    jslink((top1, 'value'), (top2, 'value'))
    #合在一起
    top = widgets.VBox([top1, top2])
    
    #mav
    mav1 = widgets.FloatRangeSlider(value=(2.0, 2.0), description='mav', min = 0.0, max= 10.0, layout=Layout(width='auto', height='auto'))
    mavstep = widgets.FloatSlider(value=1.0, min = 0.1, max= 10.0, description='mav_step', layout=Layout(width='auto', height='auto'))
    #合在一起
    mav = widgets.VBox([mav1, mavstep])
    
    #miv
    miv1 = widgets.FloatRangeSlider(value = (-1.0, -1.0), description='miv', min= -10.0, max = 0.0, layout=Layout(width='auto', height='auto'))
    mivstep = widgets.FloatSlider(value=1.0, min = 0.1, max= 10.0, description='miv_step', layout=Layout(width='auto', height='auto'))
    #合在一起
    miv = widgets.VBox([miv1, mivstep])
    
    #merge
    merge1 = widgets.FloatRangeSlider(value = (100., 100.), description='merge', min= 0.0, max = 100.0, layout=Layout(width='auto', height='auto'))
    mergestep = widgets.FloatSlider(value=1.0, min = 0.1, max= 10.0, description='merge_step', layout=Layout(width='auto', height='auto'))
    #合在一起
    merge = widgets.VBox([merge1, mergestep])
    
    #prune
    p1 = widgets.FloatRangeSlider(value = (-1.0, -1.0), description='prune | p', min= -10.0, max = 0.0, layout=Layout(width='auto', height='auto'))
    pstep = widgets.FloatSlider(value=1.0, min = 0.1, max= 10.0, description='p_step', layout=Layout(width='auto', height='auto'))
    #合在一起
    p = widgets.VBox([p1, pstep])
    
    #tqdm 进度条
    mytqdm = widgets.Checkbox(
        # value=True,
        value = False,
        description='Tqdm 循环进度条',
        disabled=False,
        indent=False)
    #jupyter 编译环境
    jupyter = widgets.Checkbox(
        value=True,
        description='Jupyter lab/nptebook编译环境',
        disabled=False,
        indent=False)
    
    #输出的路径
    output = widgets.Text(description = '保存路径', value='result', placeholder = '默认为当前目录下的 result文件夹下，可输入文件夹绝对或相对路径名称，eg： /home/user/result/ or /home/user/result or result or result/', layout=Layout(width='auto', height='auto'))
    
    #不同于前N项最优匹配的序列百分比
    diff1 = widgets.IntRangeSlider(value =(0, 0), description='diff | x', min = 0.0, max= 100.0, layout=Layout(width='auto', height='auto'))
    diffstep = widgets.IntSlider(value=10, min = 0.1, max= 100.0, description='diff_step', layout=Layout(width='auto', height='auto'))
    #合在一起
    diff = widgets.VBox([diff1, diffstep])
    
    #随机打乱的次数
    P1 = widgets.IntSlider(description='PERM', min = 0, max= 2000, layout=Layout(width='auto', height='auto'))
    P2 = widgets.IntText(layout=Layout(width='auto', height='auto'))
    P1.value = 0
    jslink((P1, 'value'), (P2, 'value'))
    #合在一起
    P = widgets.VBox([P1, P2])
    
    #随机打乱的次数
    P1 = widgets.IntSlider(description='PERM', min = 0, max= 2000, layout=Layout(width='auto', height='auto'))
    P2 = widgets.IntText(layout=Layout(width='auto', height='auto'))
    P1.value = 0
    jslink((P1, 'value'), (P2, 'value'))
    cpus = widgets.IntSlider(description='Threads', min = 1, max= cpu_count(), layout = Layout(width='auto', height='auto'))
    cpus2 = widgets.IntText(layout=Layout(width='auto', height='auto'))
    cpus.value = cpu_count()
    jslink((cpus, 'value'), (cpus2, 'value'))
    #合在一起
    P = widgets.VBox([P1, P2, cpus, cpus2])
    
    #可选参数
    options = widgets.Accordion(children=[
        Name2Type,
        dingxing,
        dingliang,
        top,
        mav,
        miv,
        p,
        merge,
        diff,
        P,
    ])
    
    title_list = ['谱系树终末节点名称转为节点类型文件',
                  '定性计算: 终末节点类型匹配得分文件',
                  '定量计算：终末节点特征文件（如细胞基因表达量）',
                  'LOCAL alignment：最优匹配的前N项结果（top=0 则为GLOBAL alignment）',
                  '终末节点一致匹配得分',
                  '终末节点错配减分',
                  '剪枝罚分',
                  '中间节点融合罚分',
                  '不同于前N项最优匹配的序列百分比',
                  '随机打乱终末节点所在位置的次数（不改变拓扑结构）',
                   ]
    
    for i in range(len(title_list)):
        options.set_title(i, title_list[i])
    
    #菜单—类型1
    accordion = widgets.Accordion(children=[
        positional,
        options,
    ])
    accordion.set_title(0, '必选参数 CLT')
    accordion.set_title(1, '可选参数')
    
    menu1 = widgets.VBox([accordion, widgets.HBox([mytqdm, jupyter]), output])
    menu2 = widgets.VBox([widgets.HTML(value="<h2><b>必选参数</b></h3>"),
                          positional, 
                          widgets.HTML(value="<h2><b>可选参数</b></h3>"),
                          widgets.Label(value="定性计算"),
                          Name2Type,
                          dingxing,
                          widgets.Label(value="定量计算"),
                          dingliang,
                          widgets.Label(), #相当于分隔符
                          top,
                          mav,
                          miv,
                          p,
                          merge,
                          diff,
                          P,
                          widgets.HBox([mytqdm, jupyter]), 
                          output])
    
    #最外层
    tab_nest = widgets.Tab()
    tab_nest.children = [menu1,
                         menu2, 
                         img_mDELTA,
                         mdelta_algorithm,
                         how
                        ]
    tab_nest.set_title(0, '菜单—类型1')
    tab_nest.set_title(1, '菜单—类型2')
    tab_nest.set_title(2, 'mdleta 匹配示例图')
    tab_nest.set_title(3, 'mDELTA 得分矩阵图')
    tab_nest.set_title(4, '参数解析')
    display(tab_nest)
    
    #if positional.top_left.value in TreeSeqFileList:
    #    positional.top_left.value = '/home/ee_while/TreeFile_nwk/' + positional.top_left.value
    #if positional.top_right.value in TreeSeqFileList:
    #    positional.top_right.value = '/home/ee_while/TreeFile_nwk/' + positional.top_left.value
    return {
        'TreeSeqFile': TreeSeqFile,
        'TreeSeqFile2': TreeSeqFile2,
        'Name2TypeFile': Name2TypeFile,
        'Name2TypeFile2': Name2TypeFile2,
        'XScoreDictFile1': XScoreDictFile,
        'LScoreDictFile2': LScoreDictFile,
        # 'ScoreDictFile1': QScoreDictFile,
        # 'ScoreDictFile2': QScoreDictFile2,
        # 'ScoreDictFile': ScoreDictFile1,
        # 'pp': PP1,
        # 'ppstep': PPstep,
        'top': top1,
        'mav':mav1,
        'mavstep':mavstep,
        'miv':miv1,
        'mivstep':mivstep,
        'p':p1,
        'pstep':pstep,
        'tqdm':mytqdm,
        'n':jupyter,
        'mg':merge1,
        'mgstep':mergestep,
        'x':diff1,
        'xstep':diffstep,
        'o':output,
        'PERM':P1,
        'CPU':cpus,
        'mdelta':'./mdelta/mDELTA.py',
        # 'MND': './feature/draw.py',
        'match_tree': './feature/match_tree.r',
        'network': './feature/network.py',
        'densitree':'./feature/densitree.r',
        'da':'./feature/densitreeALL.r'
           }

def forlist(start, end, step):
    l=[start]
    while start < end:
        start += step
        l.append(min(end,start))
    return l

def get_listvalue(l):
    l_new = []
    for i in l:
        if hasattr(i, 'value'):
            if i.value != '':
                l_new.append(i.value)
            else:
                l_new.append('non')
        else: l_new.append(i)
    return l_new

def get_output(o):
    if o.strip() == '' or o.strip() == '/' or o.strip() == 'non': # 如果为空或者为 / 的集合，默认为当前目录
        output = os.getcwd() + '/'
    elif o[0] == '/': # 如果是绝对路径
        if o[-1] != '/': 
            output = o + '/'
        else: output = o
    else: # 相对路径
        if o[-1] != '/': 
            output = os.getcwd() + '/' + o + '/'
        else: output = os.getcwd() + '/' + o
    # print(output)
    return output

def TF_to_10(n, t):
    if n == True: n = 1
    else: n = 0
    if t == True: t = 1
    else: t = 0
    # print(n, t)
    return n,t
    