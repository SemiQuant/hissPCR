# Need to make it so it collapses reads that are close

args = commandArgs(trailingOnly=TRUE)
f_in <- args[1]

name <- gsub(".counts.tsv", "", basename(f_in))

filter_thresh <- 200
require(tidyverse)
require(ggplot2)
require(cowplot)
require(plotly)
require(htmlwidgets)
dat <- read_tsv(f_in, col_names = F)

# colnames(dat) <- c("Gene", "Start", "End", "Count", "Strand")


dat$count <- 1
dat <- dat %>%
  # select(-X1) %>%
  group_by_all() %>%
  summarise(count = sum(count)) %>%
  filter(count > filter_thresh) %>%
  ungroup()


colnames(dat) <- c("Gene", "Start", "End", "Strand", "Count")

write.table(file = paste0(name, "_summarized_counts.tsv"), x = dat, sep = "\t", col.names = T, row.names = F)


round_any <- function(x, accuracy, f=round){f(x/ accuracy) * accuracy}
dat$Start <- round_any(dat$Start, 10, f = floor)
dat$End <- round_any(dat$End, 10, f = floor)
dat <- dat %>%
  group_by(Gene, Start, End, Strand) %>%
  summarise(Count = sum(Count)) %>%
  filter(Count >= filter_thresh)



plots <- list()
genes <- unique(dat$Gene)
for (i in 1:length(genes)){
  tmp <- dat %>% 
    filter(Gene == genes[i])
  plots[[i]] <- ggplot(tmp, aes(x = Start, xend = End, y = Count, yend = Count, color = Strand)) + 
    xlab("Position") + ylab("Depth") +
    geom_segment() + expand_limits(y=0)
}

# p <- ggplot(dat, aes(x = Start, xend = End, y = Count, yend = Count, color = Strand)) + 
#   geom_segment() + expand_limits(y=0)

pdf(paste0(name, ".pdf"))
ggpubr::ggarrange(plotlist = plots, labels = genes, common.legend = T, legend = "bottom")
dev.off()

dpth_plt <- read_tsv(args[2], col_names = F) %>% 
  plot_ly(x = ~X2, y=~X3, color = ~X1) %>% 
  layout(xaxis = list(title = "Postion"),
         yaxis = list(title = "Depth"),
         title = name)
try({
  saveWidget(as_widget(dpth_plt), paste0(name, "_cov.html"), selfcontained = T)
})
