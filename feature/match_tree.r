# if(!require("rjson", quietly=TRUE)){
# 	install.packages("rjson", repos = "http://cran.us.r-project.org")
# }
# if (!require("BiocManager", quietly = TRUE)){
#     install.packages("BiocManager", repos = "http://cran.us.r-project.org")
# }
# if(!require("ggtree", quietly=TRUE)){
# 	BiocManager::install("ggtree")
# }
# if(!require("ggtreeExtra", quietly=TRUE)){
# 	BiocManager::install("ggtreeExtra")
# }
options (warn = -1)

args = commandArgs(T)
# print(args)

# 导入包
suppressMessages(library(rjson))
suppressMessages(library(jsonlite))
suppressMessages(library(ggtree))
suppressMessages(library(ggplot2))
suppressMessages(library(treeio))
suppressMessages(library(dplyr))

if(!is.na(args[1])){
    # print(args[1])
	data <- jsonlite::stream_in( file(args[1]) )
}
mav = as.numeric(args[4])
#类型分数
if(!is.na(args[2]) && args[2]!= 'non'){
	data2 <- read.csv(args[2])#header = FALSE)
	data2$pair = paste(data2[,1], data2[,2], sep="--")
    
    # 完善所有的情况
    rowss = nrow(data2)
    for (i in 1:rowss) {
        rpair = paste(data2[i,2], data2[i,1], sep="--")
        # print(rpair)
        # 判断某列的值是否满足条件
        if(!(rpair %in% data2$pair)) {
            # 如果满足条件，则在新的dataframe中添加一行数据
                new_row <- data.frame(CellType1 = data2[i,2], CellType2 = data2[i,1], Scores = data2[i,3], pair = rpair)
            # 将新行添加到新的dataframe中
                data2 <- rbind(data2, new_row)
        }
    }
}else{
    data2 <- data.frame(
        CellType1 = c('xxx'),
        CellType2 = c('xxx'),
        Scores = c(mav),
        pair = c('xxx'))
}
# print(data2)

if(!is.na(args[3]) && args[3] != '' && args[3] != 'non'){
	output = args[3]
}else{
    output = ''
}
folder_path <- paste(output,'Match_tree', sep="")
if (!dir.exists(folder_path)) {
    # 如果文件夹不存在，则创建文件夹
    dir.create(folder_path)
    # print("文件夹已创建")
}

# 定义函数
#查找字符最后出现位置的“前一位”函数
LastIndexOf = function(str,str2){ #前为字符串，后为查找字符
    cd = nchar(str);
    cd2=nchar(str2);
    if(is.na(cd) ||is.na(cd2)){
        return(0);
    }
    for(i in cd:1){
        t = substr(str,i,i);
        if(t == str2){
            return(i-1)
        }
    }
    return(0)
}

TableLabel = function(atable){
    for (j in c(1:nrow(atable))){
        i = j
        while(is.na(atable[atable[i,]$parent,]$label)){
            newlabel = atable[i,]$label
            index = LastIndexOf(newlabel,"_")
            if(index != 0){
                tmpstr = substr(newlabel, 1 ,index)
                atable[atable[i,]$parent,]$label = tmpstr
            }
            
            if(atable[i,]$parent == atable[i,]$node){
                #atable[i,]$label = 'root'
                #如果为根节点则跳出
                break
            }
            i = atable[i,]$parent
        }
        
        #父节点为当前节点的，为根节点
    }
    rootindex = which(is.na(atable$label) == "TRUE")
    #print(rootindex)
    atable[rootindex,]$label = 'root'
    atable
}

TablePrune = function(btable){
    for (j in c(1:nrow(filter(btable, isTip == TRUE)))){
        i = j
        while(btable$parent[i] != btable$node[i]){
            btable$prune[btable$parent[i]] = btable$prune[btable$parent[i]] + btable$prune[j]
            
            i = btable$parent[i]
        }
    }
    btable
}

Mysize = 1
#错配的分数
mav = as.numeric(args[4])     

layoutList =list('rectangular', 
                 'dendrogram',
                 'slanted', 
                 'ellipse', 
                 'roundrect', 
                 'fan', 
                 'circular', 
                 'inward_circular', 
                 'radial', 
                 'equal_angle', 
                 'daylight' ,
                 'ape')
mylayout = layoutList[[1]]
mylayout2 = layoutList[[1]]

## 第一个位置：新建一个其实进度条
pb <- txtProgressBar(style=3, char = '>')

for(Position in 1:nrow(data)){
    #选取第几个
    #Position=5
    #print(Position)
    #data$Root1_match_label_tree[Position]
    #data$Root1_label_node[Position]
    #data$Root1_node[Position]
    
    aaa <- tryCatch({
        x = read.tree(text = data$Root1_node[Position])
        y = read.tree(text = data$Root2_node[Position])
    },error=function(e){
        # cat("ERROR :",conditionMessage(e),"\n")
        # cat(Position, " 只有单个叶子节点，跳过~\n")
    })
    if (is.null(aaa)){
        next
    }
    
    x = read.tree(text = data$Root1_node[Position])
    y = read.tree(text = data$Root2_node[Position])
    x_match = read.tree(text = data$Root1_match_tree[Position])
    #length(x_match$tip.label)
    y_match = read.tree(text = data$Root2_match_tree[Position])
    x_label = read.tree(text = data$Root1_label_node[Position])
    y_label = read.tree(text = data$Root2_label_node[Position])
    x_label_match = read.tree(text = data$Root1_match_label_tree[Position])
    y_label_match = read.tree(text = data$Root2_match_label_tree[Position])
    x_label_tree <- ggtree(x_label)
    x_label_table = TableLabel(x_label_tree$data)
    x_label_table <- as.data.frame(x_label_table)
    #x_label_table
    x_label_table$prune = 0
    y_label_tree <- ggtree(y_label)
    y_label_table = TableLabel(y_label_tree$data)
    y_label_table <- as.data.frame(y_label_table)
    #y_label_table
    y_label_table$prune = 0
    x_label_match_tree <- ggtree(x_label_match)
    y_label_match_tree <- ggtree(y_label_match)

    if(length(data$Root1_prune[Position][[1]]) == 0){
        data$Root1_prune[Position][[1]] = list('no')
    }
    #x_p <- data$Root1_prune[Position][[1]]
    x_prune <- as.data.frame(data$Root1_prune[Position])
    colnames(x_prune) <- "x_p"
    #x_prune
    #x_label_table[c(1,2,3)]
    x_prune <- merge(x_prune, x_label_table[c(1,2,3)], by.x= "x_p", by.y = "label", all.x=TRUE)
    for (i in data$Root1_prune[Position][[1]]){
        x_label_table$prune[which(x_label_table$label == i)] = 1
    }
    x_label_table <- TablePrune(x_label_table)
    x_label_table <- x_label_table %>% filter(prune > 0)
    #x_label_table
    #y_p <- data$Root2_prune[Position][[1]]
    #y_prune <- data.frame(y_p)
    if(length(data$Root2_prune[Position][[1]]) == 0){
        data$Root2_prune[Position][[1]] = list('no')
    }
    #print(length(data$Root2_prune[Position][[1]]))
    y_prune <- as.data.frame(data$Root2_prune[Position])
    colnames(y_prune) <- "y_p"
    #y_prune
    y_prune <- merge(y_prune, y_label_tree$data[c(1,2,3)], by.x= "y_p", by.y = "label", all.x=TRUE)#[[2]]
    for (i in data$Root2_prune[Position][[1]]){
        y_label_table$prune[which(y_label_table$label == i)] = 1
    }
    y_label_table <- TablePrune(y_label_table)
    y_label_table <- y_label_table %>% filter(prune > 0)
    #y_label_table

    pg1<-ggtree(x, layout = mylayout2)+
        geom_tiplab(hjust = -.3)+
        geom_tippoint(mapping=aes(subset = node %in% x_prune[[3]],
                         node = node),
                      color="#FF0000",size = 5,shape = 4,stroke = 5,)+
        geom_tiplab(hjust = -.3, mapping=aes(subset = node %in% x_prune[[3]],
                         node = node),geom = 'text', colour = alpha('firebrick', .9))
                    #geom = 'shadowtext', bg.colour = alpha('firebrick', .5))

    pg2<-ggtree(y)
    d1 <- pg1$data
    d2 <- pg2$data
    d2$x <- max(d2$x) - d2$x + max(d1$x) + 2
    pp1 <- pg1+ 
        geom_tree(data = d2, layout = mylayout2)+
        geom_tiplab(data = d2, hjust=1.5)+
        geom_tippoint(data = d2,
                      mapping=aes(subset = node %in% y_prune[[3]],
                         node = node),
                      color="#FF0000",size = 5,shape = 4,stroke = 5,)+
        geom_tiplab(data = d2,hjust = 1.5, mapping=aes(subset = node %in% y_prune[[3]],
                    node = node),geom = 'text', colour = alpha('firebrick', .9))
                    #geom = 'shadowtext', bg.colour = alpha('firebrick', .5))
    #pp1
    
    pg<-ggtree(x_match, layout = mylayout)
    #geom_tippoint(
    #    mapping = aes(
    #        color = label,
    #    ),
    #    size = Mysize,
    #)
    pg2 <- ggtree(y_match)

    d1 <- pg$data
    d1 <- as.data.frame(d1)
    d2 <- pg2$data
    d2 <- as.data.frame(d2)
    d1_label <- TableLabel(x_label_match_tree$data)
    d1_label <- as.data.frame(d1_label)
    d2_label <- TableLabel(y_label_match_tree$data)
    d2_label <- as.data.frame(d2_label)
    #d1_label
    #d2_label
    
    tree1_match_prune <- merge(d1_label, x_label_table[c(3,9)], by.x= "label", by.y = "label", all.x=TRUE)
    #tree1_match_prune
    tree2_match_prune <- merge(d2_label, y_label_table[c(3,9)], by.x= "label", by.y = "label", all.x=TRUE)
    #tree2_match_prune
    ## reverse x-axis and 
    ## set offset to make the tree on the right-hand side of the first tree

    d2$x <- max(d2$x) - d2$x + max(d1$x) + 4
    tree2_match_prune$x <- max(tree2_match_prune$x) - tree2_match_prune$x + max(tree1_match_prune$x) + 4
    
    dddd1 = d1 %>% filter(!is.na(label))
    dddd2 = d2 %>% filter(!is.na(label))
    #dddd1
    #dddd2

    sc_ <- c()
    #print(data2$pair)
    for(i in 1:nrow(dddd1)) {
        dict_index = which(data2$pair == paste(dddd1$label[i], dddd2$label[i], sep="--"))
        if (length(dict_index) == 0){
            sc_ = append(sc_, mav)
        }
        else{
            sc_ = append(sc_, data2[dict_index, 3])
        }
    }
    data_frame_sc_ <- data.frame(Score = sc_)

    pp <- pg + 
    geom_nodelab(data = tree1_match_prune,
                mapping = aes(label = prune,
                             fill = prune,
                             ),
                geom = "label")+
    geom_tiplab(hjust = -.1)+
    
    geom_tree(data=d2 , layout = mylayout)+
    geom_nodelab(data = tree2_match_prune,
                mapping = aes(label = prune,
                             fill = prune,
                             ),
                geom = "label")+
    geom_tiplab(data = d2, hjust = 1.3)+
    geom_text(data = data_frame_sc_, x = (max(d1$x)+min(d2$x))/2, y = dddd1$y - .3, label = data_frame_sc_$Score, angle = 0, )+
    geom_point(data = data_frame_sc_, x = (max(d1$x)+min(d2$x))/2, y = dddd1$y,  mapping = aes(color = Score, size = Score),)+
    scale_fill_gradientn(colors = RColorBrewer::brewer.pal(3, "YlGnBu")) +
    theme(legend.position = "right")
    

    #dddd1 = d1 %>% filter(!is.na(label))
    #dddd2 = d2 %>% filter(!is.na(label))
    #dddd1
    #dddd2
    #pp1 /  pp

    #保存的画幅大小
    mylength = max(max(length(x$tip.label)+length(x_match$tip.label),length(y$tip.label)+length(y_match$tip.label)), 10)
    #mylength

    filename = paste(output,"Match_tree/",Position, ".pdf", sep="")
    # print(filename)
    ggsave(filename, plot = pp1 /  pp,
           width = min(mylength*10,1000) ,
           height = min(mylength*12,1000),
           units = 'mm',
           limitsize = FALSE)

    ## 第二个位置：实时反映进度
    setTxtProgressBar(pb, Position/nrow(data))
    # setTxtProgressBar(pb, Position)
}

## 第三个位置关闭进度条
close(pb) 

print('matchtree ok!!!')
