# ============================================================
# File: 02_transmission_network_analysis.R
# Original file: 2016.1.26传播网络模型.R
# Purpose: Transmission network construction and node-level metrics for onward transmission
# Study: Imported malaria analysis
#
# Inputs:
# - table1_primary.csv
# - table2_edges.csv
#
# Outputs:
# - node_metrics1.27更新后.csv
#
# Notes:
# - Personal working-directory paths have been removed.
# - Run from the repository root using source("R/02_transmission_network_analysis.R").
# - Original unmodified script is retained in original_scripts/.
# ============================================================

source(file.path("R", "00_setup.R"))

# setwd() removed; project paths are defined in R/00_setup.R
table1 <- read.csv(input_file("table1_primary.csv"), stringsAsFactors = FALSE)
table2 <- read.csv(input_file("table2_edges.csv"),   stringsAsFactors = FALSE)
head(table1)
head(table2)

##数据整理
# 检查字段名
names(table1)
names(table2)
missing_primary <- setdiff(
  unique(table2$from_id),
  unique(table1$primary_id)
)

missing_primary

dup_edges <- table2[
  duplicated(table2[, c("from_id", "to_id")]),
]

dup_edges

self_loop <- table2[table2$from_id == table2$to_id, ]
self_loop

library(dplyr)

n_secondary_df <- table2 %>%
  count(from_id, name = "n_secondary")

table1 <- table1 %>%
  left_join(
    n_secondary_df,
    by = c("primary_id" = "from_id")
  )
table1$n_secondary[is.na(table1$n_secondary)] <- 0
setdiff(table2$from_id, table1$primary_id)
table(table1$primary_id)[table(table1$primary_id) > 1]

summary(table1$n_secondary)
table(table1$n_secondary)

###指标获取
library(igraph)
library(dplyr)

# 构建有向传播网络
g <- graph_from_data_frame(
  d = table2,          # from_id -> to_id
  vertices = table1,   # 节点属性表
  directed = TRUE
)

g <- graph_from_data_frame(
  d = table2 %>% select(from_id, to_id),
  vertices = NULL,
  directed = TRUE
)


vcount(g)   # 节点数
ecount(g)   # 边数

##网络层面（Global-level）拓扑指标
network_metrics <- data.frame(
  total_nodes = vcount(g),
  total_edges = ecount(g),
  network_density = edge_density(g, loops = FALSE),
  avg_path_length = average.path.length(
    g,
    directed = TRUE,
    unconnected = TRUE
  ),
  num_weak_components = components(g, mode = "weak")$no
)

network_metrics

V(g)$outdegree <- degree(g, mode = "out")
table(V(g)$outdegree == V(g)$n_secondary)

V(g)$indegree <- degree(g, mode = "in")

V(g)$betweenness <- betweenness(
  g,
  directed = TRUE,
  normalized = TRUE
)

super_spreader_threshold <- quantile(
  V(g)$outdegree,
  probs = 0.95
)
super_spreader_threshold

V(g)$is_superspreader <- V(g)$outdegree > super_spreader_threshold

ss_contribution <- sum(V(g)$outdegree[V(g)$is_superspreader]) /
  sum(V(g)$outdegree)

ss_contribution

node_metrics <- data.frame(
  primary_id       = V(g)$name,
  outdegree        = V(g)$outdegree,
  indegree         = V(g)$indegree,
  betweenness      = V(g)$betweenness,
  is_superspreader = V(g)$is_superspreader
)

head(node_metrics)

table1 <- table1 %>%
  left_join(node_metrics, by = "primary_id")

write.csv(
  node_metrics,
  file = output_table_file("node_metrics1.27更新后.csv"),
  row.names = FALSE
)

##绘图代码
library(igraph)
library(ggraph)
library(ggplot2)

# 只保留至少参与过一次传播的节点
g_sub <- induced_subgraph(
  g,
  vids = V(g)[degree(g, mode = "all") > 0]
)
V(g_sub)$outdegree <- degree(g_sub, mode = "out")
V(g_sub)$is_spreader <- V(g_sub)$outdegree > 0

p_network <- ggraph(g_sub, layout = "fr") +
  
  # 边：弱化处理，突出结构而非路径
  geom_edge_link(
    edge_width = 0.4,
    alpha = 0.3,
    colour = "grey70"
  ) +
  
  # 节点：大小 = outdegree，颜色 = 是否传播
  geom_node_point(
    aes(
      size = outdegree,
      color = is_spreader
    ),
    alpha = 0.8
  ) +
  
  scale_size_continuous(
    range = c(2, 10),
    name = "Number of secondary cases"
  ) +
  
  scale_color_manual(
    values = c("FALSE" = "grey50", "TRUE" = "firebrick"),
    labels = c("No secondary cases", "Caused secondary cases"),
    name = "Transmission role"
  ) +
  
  theme_void() +
  theme(
    legend.position = "bottom"
  )

p_network


###考虑疟疾类型的代码
# primary -> species 映射
primary_species <- table1[, c("primary_id", "Plasmodium_species")]

# 初始化
V(g)$species <- primary_species$Plasmodium_species[
  match(V(g)$name, primary_species$Primary_id)
]

# 二代继承一代疟型
edge_df <- as_data_frame(g, what = "edges")

for (i in seq_len(nrow(edge_df))) {
  from <- edge_df$from[i]
  to   <- edge_df$to[i]
  V(g)$species[to] <- V(g)$species[from]
}

# 转 factor
V(g)$species <- as.factor(V(g)$species)

p_network_species <- ggraph(g_sub, layout = "fr") +
  geom_edge_link(alpha = 0.3, colour = "grey70") +
  geom_node_point(
    aes(
      size = outdegree,
      color = species
    ),
    alpha = 0.85
  ) +
  scale_color_brewer(
    palette = "Set1",
    name = "Plasmodium species"
  ) +
  scale_size_continuous(
    range = c(2, 10),
    name = "Number of secondary cases"
  ) +
  theme_void() +
  theme(legend.position = "bottom")

p_network_species

###二代病例柱状图
outdeg <- degree(g, mode = "out")
super_spreader_threshold <- quantile(outdeg, 0.95)
library(ggplot2)

deg_df <- data.frame(outdegree = outdeg)

p_deg <- ggplot(deg_df, aes(x = outdegree)) +
  geom_histogram(
    binwidth = 1,
    fill = "steelblue",
    color = "white",
    boundary = 0
  ) +
  geom_vline(
    xintercept = super_spreader_threshold,
    linetype = "dashed",
    color = "red",
    linewidth = 1
  ) +
  labs(
    x = "Out-degree (Number of secondary cases)",
    y = "Number of cases",
    title = "Out-degree distribution of the transmission network"
  ) +
  theme_minimal(base_size = 13)

p_deg

p_deg +
  coord_cartesian(
    xlim = c(1, 45),
    ylim = c(0, 50)
  )
p_deg +
  coord_cartesian(
    xlim = c(1, 45),
    ylim = c(0, 50)
  ) +
  theme_classic(base_size = 13)


###颜色渐变饼图
ss_df <- panelC_sum[panelC_sum$group != "Other cases", ]
other_df <- panelC_sum[panelC_sum$group == "Other cases", ]
ss_df <- ss_df |>
  dplyr::arrange(desc(transmissions))
panelC_plot <- rbind(ss_df, other_df)

panelC_plot$group <- factor(
  panelC_plot$group,
  levels = panelC_plot$group   # 保留当前顺序
)

# 超级传播者数量
n_ss <- nrow(ss_df)

# 橙色 → 深蓝色 渐变
ss_colors <- colorRampPalette(
  c("#FDAE61", "#F46D43", "#D73027", "#4575B4", "#313695")
)(n_ss)

fill_colors <- c(
  ss_colors,
  "grey80"
)

names(fill_colors) <- panelC_plot$group

p_panelC <- ggplot(
  panelC_plot,
  aes(x = "", y = transmissions, fill = group)
) +
  geom_bar(
    stat = "identity",
    width = 0.6,
    color = "black",
    linewidth = 0.3
  ) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = fill_colors) +
  labs(
    title = "Contribution of individual super-spreaders to total transmission",
    fill = "Case ID"
  ) +
  theme_void(base_size = 13) +
  theme(
    legend.position = "right"
  )

p_panelC

