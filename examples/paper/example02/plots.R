################################################################################
#                  Plots of results
################################################################################

savepath <- "<path to savedir>"

# plot the Bayesian network
layout <- structure(c(0.454824742268041, 0.875334306075898, 0.0877731958762888,
                      0.435646987881181, 0.817886597938145, 0.058666007864398,
                      0.519061704242949,  0.251823319039584, -0.479173386424051,
                      -0.474869034245737, 1,  -1, -0.585642035911904, 0.987497304880946,
                      -0.949292576293966,  0.936601620325883, 0.663615675894183,
                      0.0997422680412374, 0.454824742268041,  -0.498711340206185,
                      -0.654309278350515, -0.848570236806739, -0.546587628865979,
                      -0.160525791948715, -0.482752577319587, -0.213853389296632, -1,
                      0.227412371134021, 0.291247422680413, -0.0398969072164947,
                      0.299226804123711,  0.355082474226804, -0.622389033626708,
                      -0.0558556701030927, -0.187515463917526,  0.686226804123711, 1,
                      0.801927835051547, 0.490731958762887, 0.750061855670103 ), .Dim = c(20L, 2L))

#svg(paste0(savepath,"baynet3.svg"), width=17/2.54, height=14/2.54, onefile=TRUE, antialias="subpixel")
grph <- get_network_structure(compmap$getMaps(), node_dt$NODE, obs)
par(oma=rep(0,4), mar=rep(0,4))
plot.igraph(grph, layout=layout)
#dev.off()


# compute binned averages for compariosn
binsize <- 0.2
breaks <- seq(1,2,by=binsize)
mid <- (breaks[-1] + head(breaks, n=-1)) / 2
node_dt[, EBINS:=cut(ENERGY, breaks), by=REAC]
avg_dt <- node_dt[grepl("^expdata_", NODE), list(OBSMEAN = mean(OBS)), by=c("REAC","EBINS")]
avg_dt[, ECENTER := mid[match(EBINS, levels(EBINS))]]
avg_dt[, BINSIZE := binsize]

create_plotlayers <- function(ggp, node_dt, compgrep, Emin, Emax, expalpha=0.3, expsize=0.5) {
  ensel <- node_dt[, ENERGY >= Emin & ENERGY <= Emax]
  sel <- node_dt[, grepl(compgrep, NODE) & ensel]
  ggp <- ggp + geom_errorbar(aes(x=ENERGY, ymin=OBS-UNC, ymax=OBS+UNC, col=EXPREF),
                             data=node_dt[grepl("^expdata_", NODE) & ensel], alpha=expalpha)
  ggp <- ggp + geom_point(aes(x=ENERGY, y=OBS, col=EXPREF), data=node_dt[grepl("^expdata_", NODE) & ensel], size=expsize)
  ggp <- ggp + facet_wrap(~REAC, scales="free")
  # avg component prediction
  ggp <- ggp + geom_ribbon(aes(x=ENERGY, ymin=PRED-PREDUNC, ymax=PRED+PREDUNC), data=node_dt[sel], alpha=0.2)
  ggp <- ggp + geom_line(aes(x=ENERGY, y=PRED), data=node_dt[sel])
  # plot a sample from average posterior distribution
  ggp <- ggp + geom_line(aes(x=ENERGY, y=V1), data=node_dt[grepl(compgrep, NODE) & ensel], alpha=0.4, linetype="dashed")
  ggp <- ggp + geom_line(aes(x=ENERGY, y=V2), data=node_dt[grepl(compgrep, NODE) & ensel], alpha=0.4, linetype="dashed")
  # plot a sample from average posterior distribution
  ggp <- ggp + xlim(c(Emin,Emax)) #coord_cartesian(xlim=c(1,2)) + xlim(c(1,2))
  ggp <- ggp + xlab("energy [MeV]") + ylab("cross section [mbarn]")
  return(ggp)
}

ggf <- ggplot() + theme_bw() + theme(legend.position="none")

# plot the average components
compgrep <- "^truexs_avg_[A-Z]+"
ggp <- create_plotlayers(ggf, node_dt, compgrep, 1, 2)
# binned average
ggp <- ggp + geom_segment(aes(x=ECENTER-BINSIZE/2, xend=ECENTER+BINSIZE/2, y=OBSMEAN, yend=OBSMEAN),
                          col = "blue", data=avg_dt)
ggp
#ggsave(paste0(savepath, "truexs_avg2.svg"), width=17, height=6, units="cm")

# plot the compound truexs
compgrep <- "^truexs_[A-Z]+"
ggp <- create_plotlayers(ggf, node_dt, compgrep, 1.47, 1.52, expalpha=0.8, expsize=1)
ggp
#ggsave(paste0(savepath, "truexs2.svg"), width=17, height=6, units="cm")

# create an example posterior covariance block
sel1 <- node_dt[grepl("^truexs_avg_TOT", NODE) & ENERGY >= 1 & ENERGY <= 2, IDX]
sel2 <- node_dt[grepl("^truexs_hires_EL", NODE) & ENERGY >= 1 & ENERGY <= 2, IDX]
example_covmat <- get_posterior_cov(compmap, zpost, U, obs, c(sel1,sel2) , c(sel1,sel2), ret.dep=TRUE)
example_cormat <- cov2cor(as.matrix(example_covmat))
example_cormat <- example_cormat[seq_along(sel1), length(sel1) + seq_along(sel2)]
rownames(example_cormat) <- node_dt[sel1, ENERGY]
colnames(example_cormat) <- node_dt[sel2, ENERGY]

library(reshape2)
excordt <- as.data.table(melt(example_cormat))
ggp <- ggplot() + theme_bw()
ggp <- ggp + geom_tile(aes(x=Var1, y=Var2, fill=value), data=excordt)
ggp <- ggp + xlab("truexs_avg_TOT energy [MeV]") + ylab("truexs_hires_EL energy [MeV]")
ggp <- ggp + scale_fill_continuous(type="viridis")
ggp <- ggp + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
ggp <- ggp + labs(fill="corr")
ggp <- ggp + theme(axis.title=element_text(size=9), legend.title=element_text(size=9))
ggp
#ggsave(paste0(savepath, "corr-totxs-avg-TOT-totxs-hires-EL.png"), width=8.7, height=6, units="cm", dpi=300)

# print the table with the Bayesian network information
library(xtable)
print_dt <- unique(node_dt[, list(
  PRIOR,
  UNC,
  EMIN = min(ENERGY),
  EMAX = max(ENERGY),
  NUM = .N
), by=NODE])
print_dt

# table with node prior for latex output
print.xtable(xtable(print_dt,
                    display = c("s","s","g","g","f","f","d"),
                    align = c("l","l","|","c","c","r","r","r")),
             include.rownames=FALSE,
             math.style.exponents=TRUE)

# table with experimental datasets
expinfo <- list(
  "10529004" = list(
    citekey = "pereyHighResolutionInelastic1971",
    author = "Perey et al",
    year = 1971
  ),
  "11700002" = list(
    citekey = "barrowsStudyNgReactions1965",
    year = 1965,
    author = "Barrows"
  ),
  "13764002" = list(
    citekey = NA,
    author = "Harvey",
    year = 1987
  ),
  "22316003" = list(
    citekey = NA,
    author = "Rohr et al",
    year = 1995
  ),
  "23134005" = list(
    citekey = "beyerInelasticScatteringFast2014",
    author = "Beyer et al",
    year = 2014
  ),
  "32201002" = list(
    citekey = "korzhStudyCrossSections1994",
    author = "Korzh et al",
    year = 1994
  ),
  "40532014" = list(
    citekey = "korzhStudyCrossSections1977",
    author = "Korzh et al",
    year = 1977
  )
)

expinfodt <- node_dt[grepl("^expdata_",NODE), list(NUM=.N), by=c("REAC","EXPREF")]
expinfodt[, AUTHOR := unlist(lapply(expinfo[EXPREF], function(x) x$author))]
expinfodt[, YEAR := as.integer(unlist(lapply(expinfo[EXPREF], function(x) x$year)))]
expinfodt[, REF := unlist(lapply(expinfo[EXPREF], function(x) {
  if (!is.na(x$citekey)) paste0("\\cite{",x$citekey,"}") else ""
}))]
setcolorder(expinfodt, neworder=c("REAC", "EXPREF", "NUM", "AUTHOR", "YEAR", "REF"))
expinfodt <- expinfodt[order(REAC, YEAR)]
setnames(expinfodt, c("EXPREF"), c("EXFOR"))

library(xtable)
expinfo_xtable <- xtable(expinfodt,
  caption = paste0("Experimental datasets used in the example evaluation of neutron-induced reactionf of $^{56}$Fe between one and two MeV. ",
                   "The column EXFOR contains the EXFOR accession number, NUM the number of datapoints of the datasets in the energy range considered. ",
                   "A missing reference means that no accessible publication is known to the authors."),
  label = "tbl:fe56-slow-datasets"
)

print.xtable(expinfo_xtable,
             sanitize.text.function = identity,
             include.rownames = FALSE)
