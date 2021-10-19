library(data.table)
library(nucdataBaynet)
library(igraph)

expDt <- readRDS("inpdata/expDt.rda")
refParamDt <- readRDS("inpdata/refParamDt.rda")
extNeedsDt <- readRDS("inpdata/extNeedsDt.rda")
fullSensDt <- readRDS("inpdata/fullSensDt.rda")

energy_thresholds <- list(
  "(N,D)" = 8.10259,
  "(N,2N)" = 1.13990E+01,
  "(N,P)" = 2.96572E+00,
  "(N,INL)" = 8.62048E-01,
  "(N,N+P)" = 1.03673E+01,
  "(N,EL)" = 0,
  "(N,T)" = 1.21430E+01,
  "(N,A)" = 0,
  "(N,TOT)" = 0
)

# rename the experimental reactions
exp_reacdic <- list(
  "(26-FE-56(N,P)25-MN-56,,SIG)" = "(N,P)",
  "(26-FE-56(N,EL)26-FE-56,,SIG)" = "(N,EL)",
  "(26-FE-56(N,TOT),,SIG)" = "(N,TOT)",
  "(26-FE-56(N,INL)26-FE-56,,SIG)" = "(N,INL)",
  "(26-FE-56(N,D)25-MN-55,,SIG)" = "(N,D)",
  "(26-FE-56(N,A)24-CR-53,,SIG)" = "(N,A)",
  "(26-FE-56(N,2N)26-FE-55,,SIG)" = "(N,2N)",
  "(26-FE-56(N,T)25-MN-54,,SIG)" = "(N,T)",
  "(26-FE-56(N,N+P)25-MN-55,,SIG)" = "(N,N+P)"
)

talys_reacdic <- list(
  "CS/EL" = "(N,EL)",
  "CS/REAC/000001/TOT" = "(N,A)",
  "CS/REAC/000100/TOT" = "(N,T)",
  "CS/REAC/001000/TOT" = "(N,D)",
  "CS/REAC/010000/TOT" = "(N,P)",
  "CS/REAC/100000/TOT" = "(N,INL)",
  "CS/REAC/110000/TOT" = "(N,N+P)",
  "CS/REAC/200000/TOT" = "(N,2N)",
  "CS/TOT" = "(N,TOT)"
)

expDt[, REAC := unlist(exp_reacdic[REAC])]
extNeedsDt[, REAC := unlist(talys_reacdic[REAC])]

# bin the total cross section for simplicity
totxs_dt <- expDt[REAC == "(N,TOT)"]
totxs_dt <- totxs_dt[, {
  breaks <- seq(min(L1), max(L1)+0.2, by=0.2)
  mid <- (breaks[-1] + head(breaks, n=-1)) / 2
  Efact <- cut(L1, breaks , include.lowest=TRUE)
  mymean <- tapply(DATA, Efact, mean)
  list(L1=mid, DATA=as.vector(mymean))
}, by=c("REAC","EXPID")]
totxs_dt <- totxs_dt[!is.na(DATA)]
expDt <- expDt[REAC != "(N,TOT)"]
expDt <- rbindlist(list(expDt,totxs_dt), fill=TRUE)


talysinfo <- list(
  yref = extNeedsDt$V1,
  pref = unlist(refParamDt[ADJUSTABLE==TRUE, PARVAL]),
  S = sparseMatrix(i = fullSensDt$IDX1, j = fullSensDt$IDX2, x = fullSensDt$X,
                   dims = c(nrow(extNeedsDt), nrow(refParamDt)))[, refParamDt$ADJUSTABLE],
  refParamDt = refParamDt,
  extNeedsDt = extNeedsDt,
  fullSensDt = fullSensDt,
  energy_thresholds = energy_thresholds,
  all_reacs = unique(extNeedsDt$REAC),
  indep_reacs = setdiff(unique(extNeedsDt$REAC), "(N,TOT)"),
  dep_reacs = "(N,TOT)"
)
# remove the total channel because we will explicitly define it as a sum
should_keep_row <- extNeedsDt[, REAC %in% talysinfo$indep_reacs]
talysinfo$yref <- talysinfo$yref[should_keep_row]
talysinfo$S <- talysinfo$S[should_keep_row,,drop=FALSE]
