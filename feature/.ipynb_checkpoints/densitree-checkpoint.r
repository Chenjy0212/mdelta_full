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

suppressMessages(library(rjson))
suppressMessages(library(jsonlite))
suppressMessages(library(ggtree))
suppressMessages(library(ggplot2))
suppressMessages(library(treeio))
suppressMessages(library(dplyr))

if(!is.na(args[1]) && args[1] != ''){
	data<-jsonlite::stream_in(file(args[1]))
    # print(args[1])
}
if(!is.na(args[2]) && args[2] != '' && args[2] != 'non'){
	output = args[2]
}else{
    output = ''
}

folder_path <- paste(output,'DensitreeBEST', sep="")
if (!dir.exists(folder_path)) {
    # print(paste('创建新文件夹: ',output,'DensitreeBEST', sep=""))
  # 如果文件夹不存在，则创建文件夹
  dir.create(folder_path)
  # print("文件夹已创建")
}

sorttree1 <- sort(table(data$Root1_label), decreasing = TRUE)
sorttree2 <- sort(table(data$Root2_label), decreasing = TRUE)
# print(sorttree1)

whichtree = 'Tree1'
if(max(sorttree2) >= max(sorttree1)){
    sorttree1 = sorttree2
    whichtree = 'Tree2'
    data[c('Root1_label', 'Root1_node', 'Root1_seq', 'Root1_label_node', 'Root1_match', 'Root1_match_tree', 'Root1_match_label_tree', 'Root1_prune', 'Root2_label', 'Root2_node', 'Root2_seq', 'Root2_label_node', 'Root2_match', 'Root2_match_tree', 'Root2_match_label_tree', 'Root2_prune')]<-data[c('Root2_label', 'Root2_node', 'Root2_seq', 'Root2_label_node', 'Root2_match', 'Root2_match_tree', 'Root2_match_label_tree', 'Root2_prune', 'Root1_label', 'Root1_node', 'Root1_seq', 'Root1_label_node', 'Root1_match', 'Root1_match_tree', 'Root1_match_label_tree', 'Root1_prune')]
}

Chosen_Root1 = names(sorttree1[1])

l_tmp = which(data$Root1_label==Chosen_Root1, arr.ind = TRUE)
Root1_Root2 <- data[l_tmp,]

mylayout = "slanted"
mysize = 3

tree1<- Root1_Root2$Root1_node[1] #取第一个
tree1_seq <- Root1_Root2$Root1_seq[1] #取第一个
# 设置文件路径和文件名
file_path <- paste(output,'DensitreeBEST/',"BEST_",nrow(Root1_Root2),'_',whichtree,'_',Chosen_Root1,"_sub_tree_celltype.nwk", sep="")
file_path2 <- paste(output,'DensitreeBEST/',"BEST_",nrow(Root1_Root2),'_',whichtree,'_',Chosen_Root1,"_sub_tree.nwk", sep="")

# 将字符串保存到文件
writeLines(tree1, file_path)
writeLines(tree1_seq, file_path2)

tree1_label <-Root1_Root2$Root1_label_node[1]

tree2 <- Root1_Root2$Root2_node
x = read.tree(text = tree1)
x_label = read.tree(text = tree1_label)

if(exists('mytype1')){
	for(i in 1:length(x$tip)) {
		index = which(mytype1[,2] == x$tip.label[i], arr.in=TRUE)
		if(length(index) == 0){
            x$tip.label[i] <- 'UK'
			}
		
	}
}

pg<-ggtree(x,layout=mylayout,
            alpha=0.01, 
            size = mysize,
           #colour='steelblue',
          )
pg <- pg+ 
geom_tippoint(
    mapping = aes(
        color = label,
        x = x
    ),
    size = 10,
    shape = 15,
)+theme(legend.text = element_text(size = 30), 
         legend.title = element_text(size = 30),
        legend.position="bottom", # legend至于底部
         legend.box="horizontal", # legend水平放置
        )

pg_label<-ggtree(x_label, layout=mylayout,)
#pg_label$data
#pg

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

TabelLabel = function(atabel){
    for (j in c(1:nrow(atabel))){
        i = j
        while(is.na(atabel[atabel[i,]$parent,]$label)){
            newlabel = atabel[i,]$label
            index = LastIndexOf(newlabel,"_")
            if(index != 0){
                tmpstr = substr(newlabel, 1 ,index)
                atabel[atabel[i,]$parent,]$label = tmpstr
            }
            
            if(atabel[i,]$parent == atabel[i,]$node){
                #atabel[i,]$label = 'root'
                #如果为根节点则跳出
                break
            }
            i = atabel[i,]$parent
        }
        
        #父节点为当前节点的，为根节点
    }
    rootindex = which(is.na(atabel$label) == "TRUE")
    #print(rootindex)
    atabel[rootindex,]$label = 'root'
    atabel
}

mytable <-TabelLabel(pg_label$data)
# 消除! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.
mytable <- as.data.frame(mytable)

# print(mytable)
leaves <- pg$data[pg$data$isTip == TRUE, ]
match_dataframe <- data.frame(y = leaves$y, base = leaves$label)

for(i in 1:nrow(Root1_Root2)) {
    Root2_match_tree = Root1_Root2$Root2_match_tree[i]
    # print(Root2_match_tree)
    Root1_match_label_tree = Root1_Root2$Root1_match_label_tree[i]
    #Root1_match_label_tree = '(0,(1_0,1_1));'
    x_match = read.tree(text = Root2_match_tree)
  
    if(exists('mytype2')){
        # print(mytype2)
        for(j in 1:length(x_match$tip)) {
            index = which(mytype2[,2] == x_match$tip.label[j], arr.in=TRUE)
            if(length(index) == 0){
                x_match$tip.label[j] <- 'UK'
            }
        }
    }
    x_match_label = read.tree(text = Root1_match_label_tree)
    pg_match <- ggtree(x_match, layout=mylayout,)
    pg_match_label <- ggtree(x_match_label, layout=mylayout,)
    #pg_match_label$data
    
    mymatchtable <- TabelLabel(pg_match_label$data)
    # print(mymatchtable)
    match_table_label <- merge(mytable[c(-1,-2,-4)], mymatchtable[c(1,2,3,4)], by = "label", all.y=TRUE)
    # print(match_table_label)
    match_table  <- merge(match_table_label[,-1], pg_match$data[c(2,3)], by = 'node', all.x=TRUE)
    # print(match_table)
    leaves_tmp <- match_table[match_table$isTip == TRUE, ][c('y','label')]
    # print(leaves_tmp)
    colnames(leaves_tmp)[2] <- i 
    
    match_dataframe <- merge(match_dataframe, leaves_tmp, by = 'y', all.x = TRUE)
    
    match_table$x_index = match_table$x + .5*i

    pg <- pg + geom_tree(data = match_table,
                         layout= mylayout,
                         alpha=.03,
                         size = mysize,
                         colour='steelblue',
                        )+
                geom_tippoint(data = match_table,
                  mapping = aes(
                  color = label,
                  x = x_index
                  ),
                  size = 10,
                shape = 15,)
    #theme(legend.position = 'none')
}
pg + scale_color_manual(values = c('C3' = '#A6CEE3',
                              'C1' = '#1F78B4',
                              'C2' = '#B2DF8A',
                              'C4' = '#33A02C',
                              'R2' = '#FB9A99',
                              'C5' = '#E31A1C',
                              'C7' = '#FDBF6F',
                              'C6' = '#FF7F00',
                              'C10' = '#CAB2D6',
                              'C9' = '#6A3D9A',
                              'R1' = '#FFFF99',
                              'C8' = '#B15928',
                              # 'UK' = 'grey',
                              'A' = 'red', 'B' = 'blue', 'C' = 'green', 'D' = 'yellow', 'E' = 'pink', 'F' = 'orange'))

width = max((2.2 * max(sorttree1)), 10)
# print(width)
height = max(nrow(match_dataframe), 15)/1.5  
# print(height)
width = max(width, height)
height = max(width, height)
# print(width)
# print(height)
ggsave(filename = paste(output,'DensitreeBEST/',"BEST_",nrow(Root1_Root2),"_dt.pdf", sep=""), width=width, height=height , limitsize = FALSE)

fileout = paste(output,'DensitreeBEST/',"BEST_",nrow(Root1_Root2),"_densitree_leaves_match.csv", sep="")
write.csv(match_dataframe[order(match_dataframe$y, decreasing = TRUE), ][c(-1)], fileout, row.names = FALSE)

print('densitreeBEST ok!!!')
