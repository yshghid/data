# Pathway enrichment bubble plot of GO terms

## Load package

```R
library(ggplot2)
```

## Set path
```R
setwd("/data-blog/bi3")
getwd()
```
```plain text
'/data-blog/bi3'
```

## Draw bubble plot
```R
condition <- '150_con'
gpsource <- 'REAC'

df <- read.csv(paste0("gprofiler_",condition,"_termsize.csv"))
df <- df[df$source == gpsource, ]
df$reg_type <- 'up'
df$nlog <- abs(df$negative_log10_of_adjusted_p_value)
df <- df[order(df$negative_log10_of_adjusted_p_value), ]

ggplot(df, aes(x = reorder(term_name, nlog), y = negative_log10_of_adjusted_p_value, size = intersection_size, color = nlog)) +
  geom_point(alpha = 0.6) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1)) +
  labs(title = paste0("Bubble Plot - ",gpsource," / ",condition),
       x = "Term",
       y = "-log10(p-adj)",
       size = "Intersection Size",
       color = "-log10(p-adj)") +
  scale_size(range = c(1,10)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  coord_flip()
ggsave(filename = paste0("bubble_plot_",condition,"_",gpsource,".png"), width =15, height = 9)
```
