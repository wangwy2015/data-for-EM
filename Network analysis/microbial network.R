
##读入OTU丰度和分类单元数据
otutab<-read.csv("otu_bacteria.csv", header=T, row.names= 1)
otutax<-read.csv("bacteria_tax.csv", header=T, row.names= 1)

library(amplicon)
otu<-data.frame(otutab)
Group<-data.frame(name=rownames(otutab))
format2lefse(otutab = otu,
             taxonomy = otutax,
             metadata = Group,
             groupID = "name",
             output = "LEFSE_bacteria.txt")

#筛选相对丰度大于0.01%的OTU###########################################################################
otutab<-read.csv("otu_bacteria.csv", header=T, row.names= 1)
#OTU表行名为物种，列名为处理
# 2. 读取物种注释
tax<-read.csv("bacteria_tax.csv", header=T, row.names= 1)
# OTU丰度筛选阈值，默认0.01%，0为来筛选
thre<- 0.01
# 输出文件名前缀
prefix = "tax"

# 生成各分类级汇总特征表
library(usethis)
library(devtools)#提前安装好再加载 install.packages("devtools")
devtools::install_github("microbiota/amplicon")
suppressWarnings(suppressMessages(library(amplicon)))
format2stamp(otutab, tax, thre, prefix)

otutab<-read.csv("otu_fungi.csv", header=T, row.names= 1)
#OTU表行名为物种，列名为处理
# 2. 读取物种注释
tax<-read.csv("fungi_tax.csv", header=T, row.names= 1)
# OTU丰度筛选阈值，默认0.01%，0为来筛选
thre<- 0.01
# 输出文件名前缀
prefix = "tax"

# 生成各分类级汇总特征表
library(usethis)
library(devtools)#提前安装好再加载 install.packages("devtools")
devtools::install_github("microbiota/amplicon")
suppressWarnings(suppressMessages(library(amplicon)))
format2stamp(otutab, tax, thre, prefix)
############################################################################################
#网络构建############################################################################
### 计算相关性系数，建议使用spearman系数。
library(tidyverse)
library(igraph)
library(psych)
data<-read.csv("bacteria_0.01.csv", header=T, row.names= 1)
data<-t(data)
cor <- corr.test(data, use = "pairwise",method="spearman",adjust="holm", alpha=.05,ci=FALSE) # ci=FALSE,不进行置信区间计算，数据量较大时，可以加快计算速度。
cor.r <- data.frame(cor$r) # 提取R值
cor.p <- data.frame(cor$p) # 提取p值
colnames(cor.r) = rownames(cor.r)
colnames(cor.p) = rownames(cor.p) # 变量名称中存在特殊字符，为了防止矩阵行名与列名不一致，必须运行此代码。
write.csv(cor.r,"cor.r_bacteria.csv",quote = FALSE,row.names = TRUE) # 保存结果到本地
write.csv(cor.p,"cor.p_bacteria.csv",quote = FALSE,row.names = TRUE)

knitr::kable(
  head(cor.r),
  caption = "cor.r"
)
head(cor.p)

### 确定相关性关系 保留p<=0.05且abs(r)>=0.6的变量间相关关系
cor.r[abs(cor.r) < 0.6 | cor.p > 0.05] = 0
cor.r[abs(cor.r) < 0.5] = 0
cor.r[abs(cor.r) = 1] = 0
cor.r[abs(cor.r) > 0.5] = 1
cor.r = as.matrix(cor.r)
write.csv(cor.r,"cor.r_bacteria_1.csv",quote = FALSE,row.names = TRUE)
g = graph_from_adjacency_matrix(cor.r,mode = "undirected",weighted = TRUE,diag = FALSE)

#### 将数据转换为long format进行过滤，后面绘制网络图需要节点和链接数据，在这一步可以完成格式整理
cor.r$node1 = rownames(cor.r) 
cor.p$node1 = rownames(cor.p)

r = cor.r %>% 
  gather(key = "node2", value = "r", -node1) %>%
  data.frame()

p = cor.p %>% 
  gather(key = "node2", value = "p", -node1) %>%
  data.frame()
head(r)
head(p)

#### 将r和p值合并为一个数据表
cor.data <- merge(r,p,by=c("node1","node2"))
cor.data

#### 保留p<=0.05且abs(r)>=0.6的变量间相关关系，并添加网络属性
cor.data <- cor.data %>%
  filter(abs(r) >= 0.6, p <= 0.05, node1 != node2) %>%
  mutate(
    linetype = ifelse(r > 0,"positive","negative"), # 设置链接线属性，可用于设置线型和颜色。
    linesize = abs(r) # 设置链接线宽度。
  ) # 此输出仍有重复链接，后面需进一步去除。
head(cor.data)
##############网络图节点属性整理
#### 计算每个节点具有的链接数
c(as.character(cor.data$node1),as.character(cor.data$node2)) %>%
  as_tibble() %>%
  group_by(value) %>%
  summarize(n=n()) -> vertices
colnames(vertices) <- c("node", "n")
head(vertices)
write.csv(vertices,"vertices_bacteria.csv",quote = FALSE,row.names = TRUE)
#### 添加变量分类属性
type<-read.csv("bacteria_0.01_tax.csv", header=T, row.names= 1)
vertices %>%
  left_join(type,by="node") -> vertices 

#### 网络图中节点会按照节点属性文件的顺序依次绘制，为了使同类型变量位置靠近，按照节点属性对节点进行排序。
vertices$type = factor(vertices$type,levels = unique(vertices$type))
vertices = arrange(vertices,type)
write.csv(vertices,"vertices_bacteria_node.csv",quote = FALSE,row.names = FALSE)
head(vertices)
#####################################################################################################
############构建真菌网络########################################################################

##读入OTU丰度和分类单元数据
#筛选相对丰度大于0.01%的OTU#################
otutab<-read.csv("otu_fungi.csv", header=T, row.names= 1)
#OTU表行名为物种，列名为处理
# 2. 读取物种注释
tax<-read.csv("fungi_tax.csv", header=T, row.names= 1)
# OTU丰度筛选阈值，默认0.01%，0为来筛选
thre<- 0.01
# 输出文件名前缀
prefix = "tax"

# 生成各分类级汇总特征表
library(usethis)
library(devtools)#提前安装好再加载 install.packages("devtools")
devtools::install_github("microbiota/amplicon")
suppressWarnings(suppressMessages(library(amplicon)))
format2stamp(otutab, tax, thre, prefix)

otutab<-read.csv("otu_fungi.csv", header=T, row.names= 1)
#OTU表行名为物种，列名为处理
# 2. 读取物种注释
tax<-read.csv("fungi_tax.csv", header=T, row.names= 1)
# OTU丰度筛选阈值，默认0.01%，0为来筛选
thre<- 0.01
# 输出文件名前缀
prefix = "tax"

# 生成各分类级汇总特征表
library(usethis)
library(devtools)#提前安装好再加载 install.packages("devtools")
devtools::install_github("microbiota/amplicon")
suppressWarnings(suppressMessages(library(amplicon)))
format2stamp(otutab, tax, thre, prefix)
############################################################################################
#网络构建############################################################################
### 计算相关性系数，建议使用spearman系数。
library(tidyverse)
library(igraph)
library(psych)
data<-read.csv("fungi_0.01.csv", header=T, row.names= 1)
data<-t(data)
cor <- corr.test(data, use = "pairwise",method="spearman",adjust="holm", alpha=.05,ci=FALSE) # ci=FALSE,不进行置信区间计算，数据量较大时，可以加快计算速度。
cor.r <- data.frame(cor$r) # 提取R值
cor.p <- data.frame(cor$p) # 提取p值
colnames(cor.r) = rownames(cor.r)
colnames(cor.p) = rownames(cor.p) # 变量名称中存在特殊字符，为了防止矩阵行名与列名不一致，必须运行此代码。
write.csv(cor.r,"cor.r_fungi.csv",quote = FALSE,row.names = TRUE) # 保存结果到本地
write.csv(cor.p,"cor.p_fungi.csv",quote = FALSE,row.names = TRUE)

knitr::kable(
  head(cor.r),
  caption = "cor.r"
)
head(cor.p)

### 确定相关性关系 保留p<=0.05的变量间相关关系
cor.r[abs(cor.r) < 0.3] = 0
cor.r = as.matrix(cor.r)
write.csv(cor.r,"cor.r_fungi_1.csv",quote = FALSE,row.names = TRUE)
g = graph_from_adjacency_matrix(cor.r,mode = "undirected",weighted = TRUE,diag = FALSE)

#### 将数据转换为long format进行过滤，后面绘制网络图需要节点和链接数据，在这一步可以完成格式整理
cor.r$node1 = rownames(cor.r) 
cor.p$node1 = rownames(cor.p)

r = cor.r %>% 
  gather(key = "node2", value = "r", -node1) %>%
  data.frame()

p = cor.p %>% 
  gather(key = "node2", value = "p", -node1) %>%
  data.frame()
head(r)
head(p)

#### 将r和p值合并为一个数据表
cor.data <- merge(r,p,by=c("node1","node2"))
cor.data

#### 保留p<=0.05且abs(r)>=0.6的变量间相关关系，并添加网络属性
cor.data <- cor.data %>%
  filter( p <= 0.05, node1 != node2) %>%
  mutate(
    linetype = ifelse(r > 0,"positive","negative"), # 设置链接线属性，可用于设置线型和颜色。
    linesize = abs(r) # 设置链接线宽度。
  ) # 此输出仍有重复链接，后面需进一步去除。
head(cor.data)
##############网络图节点属性整理
#### 计算每个节点具有的链接数
c(as.character(cor.data$node1),as.character(cor.data$node2)) %>%
  as_tibble() %>%
  group_by(value) %>%
  summarize(n=n()) -> vertices
colnames(vertices) <- c("node", "n")
head(vertices)
write.csv(vertices,"vertices_fungi.csv",quote = FALSE,row.names = TRUE)
#### 添加变量分类属性
type<-read.csv("fungi_0.01_tax.csv", header=T, row.names= 1)
vertices %>%
  left_join(type,by="node") -> vertices 

#### 网络图中节点会按照节点属性文件的顺序依次绘制，为了使同类型变量位置靠近，按照节点属性对节点进行排序。
vertices$type = factor(vertices$type,levels = unique(vertices$type))
vertices = arrange(vertices,type)
write.csv(vertices,"vertices_fungi_node.csv",quote = FALSE,row.names = FALSE)
head(vertices)

############################### 构建graph结构数据################################################
g <- graph_from_data_frame(cor.data, vertices = vertices, directed = FALSE )
g
write.csv(g,"g_fungi.csv",quote = FALSE,row.names = FALSE)
vcount(g) # 节点数目：82
ecount(g) # 链接数:297
#get.vertex.attribute(g) # 查看graph中包含的节点属性
#get.edge.attribute(g) # 查看graph中包含的链接属性

### 3.3 简单图
is.simple(g) # 非简单图，链接数会偏高，所以需要转换为简单图。
E(g)$weight <- 1
g <- igraph::simplify(g,
                      remove.multiple = TRUE,
                      remove.loops = TRUE,
                      edge.attr.comb = "first")
#g <- delete.vertices(g,which(degree(g) == 0)) # 删除孤立点
is.simple(g)
E(g)$weight <- 1
is.weighted(g)
vcount(g) # 节点数目：21
ecount(g) # 链接数:33
### 3.4 计算节点链接数
V(g)$degree <- degree(g)
#vertex.attributes(g)
#edge.attributes(g) 
g
### graph保存到本地
write.graph(g,file = "all.gml_bacteria",format="gml") # 直接保存graph结构，gml能保存的graph信息最多。
net.data  <- igraph::as_data_frame(g, what = "both")$edges # 提取链接属性
write.csv(net.data,"net.data_fungi.csv",quote = FALSE,row.names = FALSE) # 保存结果到本地。
head(net.data)

vertices  <- igraph::as_data_frame(g, what = "both")$vertices # 提取节点属性
write.csv(vertices,"vertices_bacteria.csv",quote = FALSE,row.names = FALSE)
head(vertices) # 直接读入之前


############################### 构建graph结构数据################################################
g <- graph_from_data_frame(cor.data, vertices = vertices, directed = FALSE )
g
write.csv(g,"g_bacteria.csv",quote = FALSE,row.names = FALSE)
vcount(g) # 节点数目：82
ecount(g) # 链接数:297
#get.vertex.attribute(g) # 查看graph中包含的节点属性
#get.edge.attribute(g) # 查看graph中包含的链接属性

### 3.3 简单图
is.simple(g) # 非简单图，链接数会偏高，所以需要转换为简单图。
E(g)$weight <- 1
g <- igraph::simplify(g,
                      remove.multiple = TRUE,
                      remove.loops = TRUE,
                      edge.attr.comb = "first")
#g <- delete.vertices(g,which(degree(g) == 0)) # 删除孤立点
is.simple(g)
E(g)$weight <- 1
is.weighted(g)
vcount(g) # 节点数目：21
ecount(g) # 链接数:33
### 3.4 计算节点链接数
V(g)$degree <- degree(g)
#vertex.attributes(g)
#edge.attributes(g) 
g
### graph保存到本地
write.graph(g,file = "all.gml_bacteria",format="gml") # 直接保存graph结构，gml能保存的graph信息最多。
net.data  <- igraph::as_data_frame(g, what = "both")$edges # 提取链接属性
write.csv(net.data,"net.data_bacteria.csv",quote = FALSE,row.names = FALSE) # 保存结果到本地。
head(net.data)

vertices  <- igraph::as_data_frame(g, what = "both")$vertices # 提取节点属性
write.csv(vertices,"vertices_bacteria.csv",quote = FALSE,row.names = FALSE)
head(vertices) # 直接读入之前保存的链接和节点属性文件，之后可直接生成graph或用于其他绘图软件绘图。
#################################################################################################################


### ####准备网络图布局数据
#?layout_in_circle # 帮助信息中，有其它布局函数。
layout1 <- layout_in_circle(g) # 径向布局适合节点较少的数据。
layout2 <- layout_with_fr(g) # fr布局。
layout3 <- layout_on_grid(g) # grid布局。
head(layout1)

### #############设置绘图颜色
#?rgb() # 可以将RGB值转换为16进制颜色名称。
#### 设置节点与分组背景色
color <- c(rgb(65,179,194,maxColorValue = 255),
           rgb(255,255,0,maxColorValue = 255),
           rgb(201,216,197,maxColorValue = 255))

names(color) <- unique(V(g)$type) # 将颜色以节点分类属性命名
V(g)$point.col <- color[match(V(g)$type,names(color))] # 设置节点颜色。
#names(color2) <- unique(V(g)$type) # 如果想要节点颜色与背景颜色不一致，则可以为节点单独设置一个颜色集。
#V(g)$point.col <- color2[match(V(g)$type,names(color2))]

#### 边颜色按照相关性正负设置
#E(g)$color <- ifelse(E(g)$linetype == "positive",rgb(255,215,0,maxColorValue = 255),"gray50")
E(g)$color <- ifelse(E(g)$linetype == "positive","red",rgb(0,147,0,maxColorValue = 255))

###  绘制fr布局网络图-不添加背景色
pdf("network_bacteria.pdf",family = "Times",width = 10,height = 12)
par(mar=c(5,2,1,2))
plot.igraph(g, layout=layout2,#更多参数设置信息?plot.igraph查看。
            
            ##节点颜色设置参数##
            vertex.color=V(g)$point.col,
            vertex.frame.color ="black",
            vertex.border=V(g)$point.col,
            ##节点大小设置参数##
            vertex.size=V(g)$degree*2,
            ##节点标签设置参数##
            vertex.label=g$name,
            vertex.label.cex=0.8,
            #vertex.label.dist=0, # 标签距离节点中心的距离，0表示标签在节点中心。
            #vertex.label.degree = 0, # 标签对于节点的位置，0-右，pi-左，-pi/2-上，pi/2-下。
            vertex.label.col="black",
            
            ##链接属性参数以edge*开头##
            edge.arrow.size=0.5,
            edge.width=abs(E(g)$r)*2,
            edge.curved = TRUE
)
##设置图例，与plot.igraph()函数一起运行##
legend(
  title = "Type",
  list(x = min(layout1[,1])-0.2,
       y = min(layout1[,2])-0.17), # 图例的位置需要根据自己的数据进行调整，后面需要使用AI手动调整。
  legend = c(unique(V(g)$type)),
  fill = color,
  #pch=1
)

legend(
  title = "|r-value|",
  list(x = min(layout1[,1])+0.4,
       y = min(layout1[,2])-0.17),
  legend = c(0.6,0.8,1.0),
  col = "black",
  lty=1,
  lwd=c(0.6,0.8,1.0)*2,
)

legend(
  title = "Correlation (±)",
  list(x = min(layout1[,1])+0.8,
       y = min(layout1[,2])-0.17),
  legend = c("positive","negative"),
  col = c("red",rgb(0,147,0,maxColorValue = 255)),
  lty=1,
  lwd=1
)

legend(
  title = "Degree",
  list(x = min(layout1[,1])+1.2,
       y = min(layout1[,2])-0.17),
  legend = c(1,seq(0,8,2)[-1]),# max(V(g)$degree)
  col = "black",
  pch=1,
  pt.lwd=1,
  pt.cex=c(1,seq(0,8,2)[-1])
)
dev.off()
