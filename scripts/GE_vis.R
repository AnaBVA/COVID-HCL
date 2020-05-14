library(ggplot2)

ace2 <- read.csv("ace2_tmprss2/ace2_hcl.txt")
tmprss2 <- read.csv("ace2_tmprss2/tmprss2_hcl.txt")

gene <- ace2

g <- ggplot(gene,aes(x=ncells_without_0/ncells,
                     y = mean_without_0,
                     label = ntissue
                     )
            ) +
    geom_point() +
    labs(title = "ACE2", x = "% cells with expression", y = " ln(CPM/100 + 1)")

g + geom_label()


ggplotly(g)
