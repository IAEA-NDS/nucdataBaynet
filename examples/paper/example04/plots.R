library(ggplot2)
library(data.table)

savepath <- "<save-dir>"

plotdt <- data.table::melt(node_dt[NODE=="truexs"], variable.name="variable", measure.vars = patterns("^ZPOST", cols=names(node_dt)))
plotdt[, variable := as.factor(smoothval[variable])]

ggp <- ggplot() + theme_bw()
ggp <- ggp + geom_errorbar(aes(x=ENERGY, ymin=OBS-UNC, ymax=OBS+UNC), width=1000, size=0.5, data=node_dt[NODE=="exp"])
ggp <- ggp + geom_point(aes(x=ENERGY, y=OBS), data=node_dt[NODE=="exp"], size=0.7)
ggp <- ggp + geom_line(aes(x=ENERGY, y=value, col=variable), data=plotdt)
ggp <- ggp + theme(legend.title=element_blank(), legend.position=c(0.9,0.9))
ggp <- ggp + xlab("energy") + ylab("variation")
#ggp <- ggp + scale_x_continuous(expand=c(0,0))
ggp <- ggp + theme(legend.position="none")
ggp

#ggsave(paste0(savepath, "sparsegp-regularization-example.svg"), width=8.7, height=5, units="cm")
