# library(dplyr)
# library(tidyr)
library(stringr)
library(ggplot2)
library(tidyverse)
library(ggh4x)

df <- read.csv("data/all.csv")
nrow(df)
theme_set(theme_bw())

datadir <- function(ntaxa, ngt, ils, nbp, m, r) {
    return(paste0(paste("data", ntaxa, ngt, ils, nbp, m, r, sep = "/"), "/"))
}

ci_by_group <- df %>%
    filter(gtee < 0.9) %>%
    mutate(across(c(ils, ngt, nbp), as.factor)) %>%      # ensure factors
    group_by(ils, ngt, nbp) %>%
    summarise(
        n    = sum(!is.na(gtee)),
        mean = mean(gtee, na.rm = TRUE),
        ci_low  = quantile(gtee, 0.1),
        ci_high = quantile(gtee, 0.9),
        .groups = "drop"
    ) %>%
    select(ils, ngt, nbp, n, mean, ci_low, ci_high) %>%
    complete(ils, ngt, nbp) %>%
    arrange(ils, ngt, nbp)
ci_by_group




p <- df %>%
    mutate(
        ngt = paste0(ngt, " gene trees"),
        nbp = paste0(nbp, " base pairs"),
        ils = factor(ils, levels=c("low", "high"), labels=c("Low ILS", "High ILS")),
        ntaxa = factor(ntaxa, levels=c(25, 50, 100, 200), labels=c("25 taxa", "50 taxa", "100 taxa", "200 taxa"))
    ) %>%
    group_by(nbp, ils) %>%
    summarise(
        y = median(gtee),
        ymin = quantile(gtee, 0.0),
        ymax = quantile(gtee, 0.75),
        # ymin = max(0, y - sd(gtee)),
        # ymax = y + sd(gtee)
    ) %>%
    ggplot(aes(x = ils, y = y, ymin = ymin, ymax = ymax, color = nbp, group = interaction(nbp, ils))) +
    geom_point(position = position_dodge(width = 1)) +
    geom_errorbar(position = position_dodge(width = 1), width = 0.5) +
    theme_bw() +
    theme(
      panel.grid.minor  = element_blank(),
      legend.position   = "bottom",
      legend.box        = "vertical",
      legend.spacing.y  = unit(-0.5, "lines")
    ) +
    guides(
      color = guide_legend(order = 2),
      shape = guide_legend(order = 1)
    ) +
    labs(
        x = "",
        y = "Gene Tree Estimation Error",
        color = "Number of Gene Trees"
    ) +
    scale_y_continuous(limits = c(0, 1))
p

pdf("figs/gtee.pdf", width=10, height=5)
p
dev.off()



# \textcolor{cyan}{percentages}
# squirrel has smaller values overall, other methods are very comparable.
cyandf <- df %>% mutate(leq = hwcd <= input_error) #%>% filter(imethod != "squirrel")
paste0(
    round(100*mean(filter(cyandf, ntaxa == 25)$leq), digits=1), "%, ",
    round(100*mean(filter(cyandf, ntaxa == 50)$leq), digits=1), "%, ",
    round(100*mean(filter(cyandf, ntaxa == 100)$leq), digits=1), "%, and ",
    round(100*mean(filter(cyandf, ntaxa == 200)$leq), digits=1), "%"
)
paste0(
    round(100*mean(filter(cyandf, m == 10)$leq), digits=1), "%, ",
    round(100*mean(filter(cyandf, m == 20)$leq), digits=1), "%"
)
paste0(
    round(100*mean(filter(cyandf, m == 10)$leq), digits=1), "%, ",
    round(100*mean(filter(cyandf, m == 20)$leq), digits=1), "%"
)
paste0(
    round(100*mean(filter(cyandf, ils == "low")$leq), digits=1), "%, ",
    round(100*mean(filter(cyandf, ils == "high")$leq), digits=1), "%"
)


df$gtee_bin <- if_else(df$gtee < 0.25, "low", if_else(df$gtee < 0.5, "mid", if_else(df$gtee < 0.75, "high", "ultra-high")))
df$gtee_bin <- if_else(df$gtee < 0.35, "low", if_else(df$gtee < 0.7, "mid", "high"))
ggplot(df, aes(x = input_error, y = hwcd, color = gtee_bin)) + geom_point() + facet_wrap(~gtee_bin)

summary(aov(hwcd ~ input_error + gtee_bin, df))




for(idx in which(df$gtee > 0.75 | (df$gtee > 0.6 & df$ils == "Low ILS"))) {
    row <- df[idx,]
    cat(paste("rm -r data", row$ntaxa, row$ngt, row$ils, row$nbp, row$m, row$r, "*", sep="/"), ";  ")
    cat(paste("./scripts/perform_simulation.sh", row$ntaxa, row$ngt, row$ils, row$nbp, row$m, row$r, "fake", sep=" "))
}
