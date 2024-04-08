options(warn = -1)

args <- commandArgs(T)

suppressMessages(library(rjson))
suppressMessages(library(jsonlite))
suppressMessages(library(ggtree))
suppressMessages(library(ggplot2))
suppressMessages(library(treeio))
suppressMessages(library(dplyr))

if (!is.na(args[1]) && args[1] != "") {
    suppressMessages(data <- jsonlite::stream_in(file(args[1]), verbose = FALSE))
    # cat(args[1])
}
if (!is.na(args[2]) && args[2] != "" && args[2] != "non") {
    output <- args[2]
} else {
    output <- ""
}

# celltype color
if (!is.na(args[3]) && args[3] != "non") {
    ccfile <- read.csv(file = args[3])
    cc <- c()
    # 遍历数据框的列
    for (i in 1:nrow(ccfile)) {
        cc[ccfile[i, 1]] <- ccfile[i, 2]
    }
}

folder_path <- paste(output, "DensitreeALL", sep = "")
if (!dir.exists(folder_path)) {
    # 如果文件夹不存在，则创建文件夹
    dir.create(folder_path)
    # cat("文件夹已创建")
}

# 定义函数
# 查找字符最后出现位置的“前一位”函数
LastIndexOf <- function(str, str2) { # 前为字符串，后为查找字符
    cd <- nchar(str)
    cd2 <- nchar(str2)
    if (is.na(cd) || is.na(cd2)) {
        return(0)
    }
    for (i in cd:1) {
        t <- substr(str, i, i)
        if (t == str2) {
            return(i - 1)
        }
    }
    return(0)
}

TabelLabel <- function(atabel) {
    for (j in c(1:nrow(atabel))) {
        i <- j
        while (is.na(atabel[atabel[i, ]$parent, ]$label)) {
            newlabel <- atabel[i, ]$label
            index <- LastIndexOf(newlabel, "_")
            if (index != 0) {
                tmpstr <- substr(newlabel, 1, index)
                atabel[atabel[i, ]$parent, ]$label <- tmpstr
            }

            if (atabel[i, ]$parent == atabel[i, ]$node) {
                # atabel[i,]$label = 'root'
                # 如果为根节点则跳出
                break
            }
            i <- atabel[i, ]$parent
        }

        # 父节点为当前节点的，为根节点
    }
    rootindex <- which(is.na(atabel$label) == "TRUE")
    # cat(rootindex)
    atabel[rootindex, ]$label <- "root"
    atabel
}

# 保存 densitree info
info <- data.frame(var1 = "", var2 = "", var3 = "")[-1, ]
# cat(info)

if (nrow(table(data$Root1_label)) == 1) {
    cat("只有单一匹配，结果请查看DensitreeBEST~")
    # 结束R会话
    q()
}

sorttree1 <- sort(table(data$Root1_label), decreasing = TRUE)
sorttree2 <- sort(table(data$Root2_label), decreasing = TRUE)
whichtree <- "Tree1"

match_times1 <- as.numeric(sorttree1)
match_times2 <- as.numeric(sorttree2)

## 第一个位置：新建一个其实进度条
pb1 <- txtProgressBar(style = 3, char = ">")
for (iii in 1:nrow(sorttree1)) {
    Chosen_Root1 <- names(sorttree1[iii])

    l_tmp <- which(data$Root1_label == Chosen_Root1, arr.ind = TRUE)
    Root1_Root2 <- data[l_tmp, ]
    # 如果 densitree 过于少，例如少于 2，则不画了
    if (nrow(Root1_Root2) < 2) {
        setTxtProgressBar(pb1, iii / nrow(sorttree1))
        next
    }

    mylayout <- "slanted"
    mysize <- 3

    tree1 <- Root1_Root2$Root1_node[1] # 取第一个
    tree1_seq <- Root1_Root2$Root1_seq[1] # 取第一个
    # cat(tree1)
    # 设置文件路径和文件名
    file_path <- paste(output, "DensitreeALL/", nrow(Root1_Root2), "_", whichtree, "_", Chosen_Root1, "_sub_tree_celltype.nwk", sep = "")
    file_path2 <- paste(output, "DensitreeALL/", nrow(Root1_Root2), "_", whichtree, "_", Chosen_Root1, "_sub_tree.nwk", sep = "")

    # 将字符串保存到文件
    # writeLines(tree1, file_path)
    writeLines(tree1_seq, file_path2)

    tree1_label <- Root1_Root2$Root1_label_node[1]

    tree2 <- Root1_Root2$Root2_node

    x <- read.tree(text = tree1)
    if (length(x$tip.label) <= 1) {
        setTxtProgressBar(pb1, iii / nrow(sorttree1))
        next()
    }
    x_label <- read.tree(text = tree1_label)

    if (exists("mytype1")) {
        for (i in 1:length(x$tip)) {
            index <- which(mytype1[, 2] == x$tip.label[i], arr.in = TRUE)
            if (length(index) == 0) {
                x$tip.label[i] <- "UK"
            }
        }
    }

    pg <- ggtree(x,
        layout = mylayout,
        alpha = 0.01,
        size = mysize,
        # colour='steelblue',
    )
    pg <- pg +
        geom_tippoint(
            mapping = aes(
                color = label,
                x = x
            ),
            size = 10,
            shape = 15,
        ) + theme(
            legend.text = element_text(size = 30),
            legend.title = element_text(size = 30),
            legend.position = "bottom", # legend至于底部
            legend.box = "horizontal", # legend水平放置
        )

    pg_label <- ggtree(x_label, layout = mylayout, )

    mytable <- TabelLabel(pg_label$data)
    # 消除! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.
    mytable <- as.data.frame(mytable)

    # cat(mytable)
    leaves <- pg$data[pg$data$isTip == TRUE, ]
    match_dataframe <- data.frame(y = leaves$y, base = leaves$label)

    base_prune_count <- c(0, "base_prune") # 统计base剪枝个数
    prune_count <- c(0, "prune") # 统计对应的剪枝个数
    score_count <- c(0, "Score") # 统计对应的得分
    tree_seq <- c(0, tree1_seq) # 子树的原结构
    Serial_Number <- c(0, "Serial Number")

    # 对齐标志位
    alig_index <- 1
    for (i in 1:nrow(Root1_Root2)) {
        Root2_match_tree <- Root1_Root2$Root2_match_tree[i]
        # cat(Root2_match_tree)
        Root1_match_label_tree <- Root1_Root2$Root1_match_label_tree[i]
        # cat(Root1_match_label_tree)

        x_match <- read.tree(text = Root2_match_tree)
        if (length(x_match$tip.label) <= 1) {
            next()
        }
        if (exists("mytype2")) {
            # cat(mytype2)
            for (j in 1:length(x_match$tip)) {
                index <- which(mytype2[, 2] == x_match$tip.label[j], arr.in = TRUE)
                if (length(index) == 0) {
                    x_match$tip.label[j] <- "UK"
                }
            }
        }
        x_match_label <- read.tree(text = Root1_match_label_tree)
        pg_match <- ggtree(x_match, layout = mylayout, )
        # cat(pg_match$data)
        pg_match_label <- ggtree(x_match_label, layout = mylayout, )
        # cat(pg_match_label$data)

        mymatchtable <- TabelLabel(pg_match_label$data)
        # cat(mymatchtable)
        match_table_label <- merge(mytable[c(-1, -2, -4)], mymatchtable[c(1, 2, 3, 4)], by = "label", all.y = TRUE)
        # cat(match_table_label)
        match_table <- merge(match_table_label[, -1], pg_match$data[c(2, 3)], by = "node", all.x = TRUE)
        # cat(match_table)
        leaves_tmp <- match_table[match_table$isTip == TRUE, ][c("y", "label")]
        # cat(leaves_tmp)
        colnames(leaves_tmp)[2] <- i

        match_dataframe <- merge(match_dataframe, leaves_tmp, by = "y", all.x = TRUE)

        match_table$x_index <- match_table$x + .5 * i

        pg <- pg + geom_tree(
            data = match_table,
            layout = mylayout,
            alpha = .03,
            size = mysize,
            colour = "steelblue",
        ) +
            geom_tippoint(
                data = match_table,
                mapping = aes(
                    color = label,
                    x = x_index
                ),
                size = 10,
                shape = 15,
            )
        # theme(legend.position = 'none')
        base_prune_count[alig_index + 2] <- lengths(Root1_Root2$Root1_prune[i])
        prune_count[alig_index + 2] <- lengths(Root1_Root2$Root2_prune[i])
        score_count[alig_index + 2] <- Root1_Root2$Score[i]
        tree_seq[alig_index + 2] <- Root1_Root2$Root2_seq[i]
        Serial_Number[alig_index + 2] <- l_tmp[i]
        alig_index <- alig_index + 1
    }
    match_dataframe <- rbind(match_dataframe, base_prune_count, prune_count, score_count, tree_seq, Serial_Number)
    # pg + scale_color_manual(values = c(
    #     "C3" = "#A6CEE3",
    #     "C1" = "#1F78B4",
    #     "C2" = "#B2DF8A",
    #     "C4" = "#33A02C",
    #     "R2" = "#FB9A99",
    #     "C5" = "#E31A1C",
    #     "C7" = "#FDBF6F",
    #     "C6" = "#FF7F00",
    #     "C10" = "#CAB2D6",
    #     "C9" = "#6A3D9A",
    #     "R1" = "#FFFF99",
    #     "C8" = "#B15928",
    #     # 'UK' = 'grey',
    #     "A" = "#a11f1f", "B" = "#0e0ea6", "C" = "#246f24", "D" = "#afaf1e", "E" = "#aa707a", "F" = "#c58204", "X" = "#5c5450"
    # ))
    if (!is.na(args[3]) && args[3] != "non") {
        pg <- pg + scale_color_manual(values = cc)
    }

    width <- max((2.3 * match_times1[iii]), 15)
    # cat(width)
    height <- max(nrow(match_dataframe), 15)
    # cat(height)
    width <- max(width, height)
    height <- max(width, height)
    filename <- paste(nrow(Root1_Root2), "_", whichtree, "-", Chosen_Root1, sep = "")
    filenamepic <- paste(output, "DensitreeALL/", filename, "_dt.pdf", sep = "")
    ggsave(filename = filenamepic, width = width, height = width, limitsize = FALSE)

    fileout <- paste(output, "DensitreeALL/", filename, "_densitree_leaves_match.csv", sep = "")
    write.csv(match_dataframe[order(match_dataframe$y, decreasing = TRUE), ][c(-1)], fileout, row.names = FALSE)


    info_row <- c(match_times1[iii], tree1, filename)
    info <- rbind(info, info_row)

    setTxtProgressBar(pb1, iii / nrow(sorttree1))
    # setTxtProgressBar(pb, Position)
}
close(pb1)

# 到tree2
whichtree <- "Tree2"
sorttree1 <- sorttree2
data[c("Root1_label", "Root1_node", "Root1_seq", "Root1_label_node", "Root1_match", "Root1_match_tree", "Root1_match_label_tree", "Root1_prune", "Root2_label", "Root2_node", "Root2_seq", "Root2_label_node", "Root2_match", "Root2_match_tree", "Root2_match_label_tree", "Root2_prune")] <- data[c("Root2_label", "Root2_node", "Root2_seq", "Root2_label_node", "Root2_match", "Root2_match_tree", "Root2_match_label_tree", "Root2_prune", "Root1_label", "Root1_node", "Root1_seq", "Root1_label_node", "Root1_match", "Root1_match_tree", "Root1_match_label_tree", "Root1_prune")]
## 第一个位置：新建一个进度条
pb2 <- txtProgressBar(style = 3, char = ">")
for (iii in 1:nrow(sorttree1)) {
    Chosen_Root1 <- names(sorttree1[iii])

    l_tmp <- which(data$Root1_label == Chosen_Root1, arr.ind = TRUE)
    Root1_Root2 <- data[l_tmp, ]
    if (nrow(Root1_Root2) < 2) {
        setTxtProgressBar(pb2, iii / nrow(sorttree1))
        next
    }

    mylayout <- "slanted"
    mysize <- 3

    tree1 <- Root1_Root2$Root1_node[1] # 取第一个
    tree1_seq <- Root1_Root2$Root1_seq[1] # 取第一个

    # 设置文件路径和文件名
    file_path <- paste(output, "DensitreeALL/", nrow(Root1_Root2), "_", whichtree, "_", Chosen_Root1, "_sub_tree.nwk", sep = "")

    # 将字符串保存到文件
    writeLines(tree1, file_path)

    tree1_label <- Root1_Root2$Root1_label_node[1]

    tree2 <- Root1_Root2$Root2_node

    x <- read.tree(text = tree1)
    if (length(x$tip.label) <= 1) {
        setTxtProgressBar(pb2, iii / nrow(sorttree1))
        next()
    }
    x_label <- read.tree(text = tree1_label)

    if (exists("mytype1")) {
        for (i in 1:length(x$tip)) {
            index <- which(mytype1[, 2] == x$tip.label[i], arr.in = TRUE)
            if (length(index) == 0) {
                x$tip.label[i] <- "UK"
            }
        }
    }

    pg <- ggtree(x,
        layout = mylayout,
        alpha = 0.01,
        size = mysize,
        # colour='steelblue',
    )
    pg <- pg +
        geom_tippoint(
            mapping = aes(
                color = label,
                x = x
            ),
            size = 10,
            shape = 15,
        ) + theme(
            legend.text = element_text(size = 30),
            legend.title = element_text(size = 30),
            legend.position = "bottom", # legend至于底部
            legend.box = "horizontal", # legend水平放置
        )

    pg_label <- ggtree(x_label, layout = mylayout, )

    mytable <- TabelLabel(pg_label$data)
    # 消除! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.
    mytable <- as.data.frame(mytable)

    # cat(mytable)
    leaves <- pg$data[pg$data$isTip == TRUE, ]
    match_dataframe <- data.frame(y = leaves$y, base = leaves$label)

    base_prune_count <- c(0, "base_prune") # 统计base剪枝个数
    prune_count <- c(0, "prune") # 统计对应的剪枝个数
    score_count <- c(0, "Score") # 统计对应的得分
    tree_seq <- c(0, tree1_seq) # 子树的原结构
    Serial_Number <- c(0, "Serial Number")

    # 对齐标志位
    alig_index <- 1
    for (i in 1:nrow(Root1_Root2)) {
        Root2_match_tree <- Root1_Root2$Root2_match_tree[i]
        # cat(Root2_match_tree)
        Root1_match_label_tree <- Root1_Root2$Root1_match_label_tree[i]

        x_match <- read.tree(text = Root2_match_tree)
        if (length(x_match$tip.label) <= 1) {
            next()
        }

        if (exists("mytype2")) {
            # cat(mytype2)
            for (j in 1:length(x_match$tip)) {
                index <- which(mytype2[, 2] == x_match$tip.label[j], arr.in = TRUE)
                if (length(index) == 0) {
                    x_match$tip.label[j] <- "UK"
                }
            }
        }
        x_match_label <- read.tree(text = Root1_match_label_tree)
        pg_match <- ggtree(x_match, layout = mylayout, )
        pg_match_label <- ggtree(x_match_label, layout = mylayout, )
        # pg_match_label$data

        mymatchtable <- TabelLabel(pg_match_label$data)
        # cat(mymatchtable)
        match_table_label <- merge(mytable[c(-1, -2, -4)], mymatchtable[c(1, 2, 3, 4)], by = "label", all.y = TRUE)
        # cat(match_table_label)
        match_table <- merge(match_table_label[, -1], pg_match$data[c(2, 3)], by = "node", all.x = TRUE)
        # cat(match_table)
        leaves_tmp <- match_table[match_table$isTip == TRUE, ][c("y", "label")]
        # cat(leaves_tmp)
        colnames(leaves_tmp)[2] <- i

        match_dataframe <- merge(match_dataframe, leaves_tmp, by = "y", all.x = TRUE)

        match_table$x_index <- match_table$x + .5 * i

        pg <- pg + geom_tree(
            data = match_table,
            layout = mylayout,
            alpha = .03,
            size = mysize,
            colour = "steelblue",
        ) +
            geom_tippoint(
                data = match_table,
                mapping = aes(
                    color = label,
                    x = x_index
                ),
                size = 10,
                shape = 15,
            )
        # theme(legend.position = 'none')
        base_prune_count[alig_index + 2] <- lengths(Root1_Root2$Root1_prune[i])
        prune_count[alig_index + 2] <- lengths(Root1_Root2$Root2_prune[i])
        score_count[alig_index + 2] <- Root1_Root2$Score[i]
        tree_seq[alig_index + 2] <- Root1_Root2$Root2_seq[i]
        Serial_Number[alig_index + 2] <- l_tmp[i]
        alig_index <- alig_index + 1
    }
    match_dataframe <- rbind(match_dataframe, base_prune_count, prune_count, score_count, tree_seq, Serial_Number)
    # pg + scale_color_manual(values = c(
    #     "C3" = "#A6CEE3",
    #     "C1" = "#1F78B4",
    #     "C2" = "#B2DF8A",
    #     "C4" = "#33A02C",
    #     "R2" = "#FB9A99",
    #     "C5" = "#E31A1C",
    #     "C7" = "#FDBF6F",
    #     "C6" = "#FF7F00",
    #     "C10" = "#CAB2D6",
    #     "C9" = "#6A3D9A",
    #     "R1" = "#FFFF99",
    #     "C8" = "#B15928",
    #     # 'UK' = 'grey',
    #     "A" = "#a11f1f", "B" = "#0e0ea6", "C" = "#246f24", "D" = "#afaf1e", "E" = "#aa707a", "F" = "#c58204", "X" = "#5c5450"
    # ))
    if (!is.na(args[3]) && args[3] != "non") {
        pg <- pg + scale_color_manual(values = cc)
    }

    width <- max((2.3 * match_times2[iii]), 15)
    # cat(width)
    height <- max(nrow(match_dataframe), 15)
    # cat(height)
    width <- max(width, height)
    height <- max(width, height)
    filename <- paste(nrow(Root1_Root2), "_", whichtree, "-", Chosen_Root1, sep = "")
    filenamepic <- paste(output, "DensitreeALL/", filename, "_dt.pdf", sep = "")
    ggsave(filename = filenamepic, width = width, height = width, limitsize = FALSE)

    fileout <- paste(output, "DensitreeALL/", filename, "_densitree_leaves_match.csv", sep = "")
    write.csv(match_dataframe[order(match_dataframe$y, decreasing = TRUE), ][c(-1)], fileout, row.names = FALSE)

    info_row <- c(match_times2[iii], tree1, filename)
    info <- rbind(info, info_row)

    setTxtProgressBar(pb2, iii / nrow(sorttree1))
    # setTxtProgressBar(pb, Position)
}
close(pb2)

colnames(info) <- c("Times", "SubTree", "Filename")
# cat(info)
write.csv(info, paste(output, "DensitreeALL/densitrees_info.csv", sep = ""), row.names = FALSE, col.names = TRUE)

cat("densitreeALL ok!!!")
