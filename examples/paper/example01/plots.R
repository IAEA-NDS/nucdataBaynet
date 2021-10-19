################################################################################
#                 Produce plots showing results
################################################################################

library(igraph)
library(ggplot2)
library(plotly)

savepath <- "<path to savedir>"

# plot the Bayesian network
grph <- get_network_structure(compmap$getMaps(), node_dt$NODE, node_dt$OBS)
layout <- structure(c(-0.0184285714285714, 0.845529758060234, -0.881099020449749,
                      -1, 0.521711364697923, -0.656405973332188, -0.708792782471717,
                      0.994589018496853, 0.0552857142857142, 1, -0.0368571428571428,
                      0.589714285714286, 0.850551332833744, 1, 0.450735438405833, -0.436548549405845,
                      -0.349894207552815, -0.902673212457804, -1, -0.760178571428571,
                      -0.147428571428571, -0.133607142857143), .Dim = c(11L, 2L))
# svg(file.path(savepath, "U5nf-baynet.svg"), width=2.75, height=2.75, onefile=TRUE, antialias="subpixel")
par(mar=c(0,0,0,0), oma=rep(0,4))
plot.igraph(grph, layout = layout)
# dev.off()

energysel <- node_dt[, ENERGY > 7000 & ENERGY < 9000]

# plot posterior of normalization uncertainties
expsel <- node_dt[, grepl("Weston", EXPREF) | grepl("Blons", EXPREF)]
sel <- node_dt[NODE=="relsyserr" & energysel & expsel, IDX]
ggp <- ggplot(data=node_dt[sel]) + theme_bw()
ggp <- ggp + geom_ribbon(aes(x=ENERGY, ymin=PRED-PREDUNC, ymax=PRED+PREDUNC, fill=EXPREF), alpha=0.2)
ggp <- ggp + geom_line(aes(x=ENERGY, y=PRED, col=EXPREF))
ggp <- ggp + theme(legend.position="none")
ggp <- ggp + coord_cartesian(ylim=c(-0.07, 0.05), expand=0)
ggp <- ggp + xlab("energy [eV]") + ylab("relative systematic error")
ggp <- ggp + scale_y_continuous(labels = scales::percent_format(suffix="\\%"))
#ggplotly(ggp)
ggp
# ggsave(paste0(savepath, "U5nf-relnormerr.svg"), width=8.75, height=7, units="cm")

ggp <- ggplot() + theme_bw()
ggp <- ggp + geom_point(aes(x=ENERGY, y=OBS, col=EXPREF), data=node_dt[NODE=="expdata" & energysel])

ggp <- ggp + geom_ribbon(aes(x=ENERGY, ymin=PRED-PREDUNC, ymax=PRED+PREDUNC), data=node_dt[NODE=="truexs_avg"],
                         fill="green", alpha=0.8)
ggp <- ggp + geom_line(aes(x=ENERGY, y=PRED), data=node_dt[NODE=="truexs_avg"], col="black")

ggp <- ggp + geom_line(aes(x=ENERGY, y=PRED), data=node_dt[NODE=="truexs"], col="brown", alpha=0.4)
ggp <- ggp + theme(legend.position = "none")
ggp <- ggp + xlab("energy [eV]") + ylab("cross section [barn]")
ggp <- ggp + coord_cartesian(xlim=c(7000,9000), ylim=c(0,6), expand=0)
ggp
#ggplotly(ggp)
# ggsave(paste0(savepath, "U5nf-evaluation.svg"), width=8.75, height=7, units="cm")

# prior specification as table output for latex

library(xtable)
print_dt <- unique(node_dt[, list(
  PRIOR,
  UNC,
  EMIN = min(ENERGY),
  EMAX = max(ENERGY),
  NUM = .N
), by=NODE])
print_dt

# for latex output
print.xtable(xtable(print_dt,
                    display = c("s","s","g","g","d","d","d"),
                    align = c("l","l","|","c","c","r","r","r")),
             include.rownames=FALSE,
             math.style.exponents=TRUE)


# table with experimental datasets for latex
expinfo <- list(
  "12877" = list(
    citekey = "westonSubthresholdFissionCross1984",
    author = "Weston et al",
    year = 1984
  ),
  "20483" = list(
    citekey = "blonsMeasurementAnalysisFission1971",
    author = "Blons et al",
    year = 1971
  ),
  "20783" = list(
    citekey = "mignecoIntermediateStructureKeV1975",
    author = "Migneco et al",
    year = 1975
  ),
  "20826" = list(
    citekey = "wagemansNeutronInducedFission1976",
    author = "Wagemans et al",
    year = 1976
  ),
  "23294" = list(
    citekey = "paradelaHighAccuracy2352016",
    author = "Paradela et al",
    year = 2016
  )
)

expinfodt <- node_dt[grepl("^expdata",NODE), list(NUM=.N), by=c("EXFOR")]
expinfodt[,  AUTHOR := unlist(lapply(expinfo[as.character(EXFOR)], function(x) x$author))]
expinfodt[, YEAR := as.integer(unlist(lapply(expinfo[as.character(EXFOR)], function(x) x$year)))]
expinfodt[, REF := unlist(lapply(expinfo[as.character(EXFOR)], function(x) {
  if (!is.na(x$citekey)) paste0("\\cite{",x$citekey,"}") else ""
}))]
setcolorder(expinfodt, neworder=c("EXFOR", "NUM", "AUTHOR", "YEAR", "REF"))
expinfodt <- expinfodt[order(YEAR)]
expinfodt[, NUM := as.integer(NUM)]
expinfodt[, EXFOR := as.character(EXFOR)]

library(xtable)
expinfo_xtable <- xtable(expinfodt,
  caption = paste0("Experimental datasets used in the example evaluation of $^{235}$U(n,f) reaction between seven and 9 keV. ",
                   "The column EXFOR contains the EXFOR accession number, NUM the number of datapoints of the datasets in the energy range considered. "),
  label = "tbl:u5nf-datasets"
)

print.xtable(expinfo_xtable,
             sanitize.text.function = identity,
             include.rownames = FALSE)
