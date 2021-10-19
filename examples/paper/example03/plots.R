################################################################################
#                  Produce relevant plots
################################################################################

library(ggplot2)
#library(plotly)

savepath <- "<path to savedir>"

# plot the Bayesian network

graphlayout <- structure(c(
  0.948258100558659, 0.872601226517603, 0.340317651389185,
  -0.917058322346441, -0.558661370212771, 0.896369832402235, 0.115802316756472,
  -0.831914018494573, -0.356731843575419, 0.914530726256983, 0.6322523360211,
  0.565967576000719, -0.887191909668635, -0.268957167238054, 0.927502793296089,
  -0.125865348732544, -0.796374655922659, 0.354137430167598, 0.683627932960894,
  0.530954277280521, 0.0438200554549677, -1, -0.788632747831876,
  0.740705027932961, -0.154296838790076, -0.542860536243004, 0.0402134078212292,
  -0.118045810055866, -0.258144134078212, 0.489046927374301, -0.651247249513969,
  0.455319553072626, 0.565967576000719, 0.517585474860335, 0.241374731177236,
  -0.362794432545306, -0.104691897510122, 0.924908379888268, -0.719939998672922,
  0.686222346368715, 0.120640223463687, 0.722544134078212, 0.302976292968554,
  -0.54759911791926, -0.343623198932569, 1, -0.933415963304882,
  0.841887150837989, 0.133612290502793, 0.943069273743017, 0.366947145597999,
  -0.796374655922659, -0.540741522606088, 0.758054054072224, -0.574662569832402,
  0.673250279329609, 0.364577854759872, 0.535746368715084, 0.0162921015551126,
  -0.49784401031858, -0.337649916397008, 0.167339664804469, -0.229605586592179,
  -0.97160782122905, 1, 0.279784685062276, 0.965607606536671, 0.499424581005587,
  -0.999853641996033, -0.648671131848133, 0.247766480446927, -0.439753072625698,
  -0.997173128351175, 0.904524756411861, 0.522215113118341, 0.982348618432271,
  0.185500558659218, -0.955898199628966, -0.410910696374057, 0.118045810055866,
  -0.476074860335195, -0.7822156424581, 0.996003901308808, 0.0118352645792558,
  0.729556926587344, 0.102479329608939, -0.769114103939879, -0.302367888875023,
  0.348948603351955, -0.502018994413408, -0.12582905027933, -0.585040223463687,
  0.474367002317802, 0.387864804469274, 0.653842558140281, -0.0479966480446926,
  -0.595417877094972, -0.873509804524704, 0.789355807855007, -0.595417877094972,
  -0.0934305791819305, 0.522774301675978, 0.553907262569833, -0.162150837988827,
  -0.769243575418994, -0.966546496666734, 0.463239285373228, -0.84448156424581,
  -0.262493904010503, 0.701788826815643, 0.787404469273743, -0.0454022346368714,
  -1, -0.870925451965203, 0.618225355463578, -0.748488268156424,
  0.149178770949721, 0.29446592178771, 0.770138423317818, -0.3178156424581,
  -0.670655865921788, -0.627996311372126, 0.757067043252851, -0.216633519553072
), .Dim = c(62L, 2L))

# svg(paste0(savepath,"fe56-fast-baynet.svg"), width=17/2.54, height=17/2.54, onefile=TRUE, antialias="subpixel")
par(oma=rep(0,4), mar=rep(0,4), cex=0.7)
igraph_options(vertex.size = 12)
grph <- get_network_structure(compmap$getMaps(), node_dt$NODE, node_dt$OBS)
plot(grph, layout=graphlayout)
# dev.off()


# plot the reaction channels
xsel <- node_dt[, !grepl("^truexs_",NODE) | PRED < 6000]
ggp <- ggplot() + theme_bw()

# prior uncertainty (model pars + model defect)
ggp <- ggp + geom_ribbon(aes(x=ENERGY, ymin=PRIORPRED-PRIORUNC, ymax=PRIORPRED+PRIORUNC),
                         data=node_dt[xsel & grepl("^truexs_",NODE)], alpha=0.01, fill="green", col="#eeeeee")

# experimental data
ggp <- ggp + geom_errorbar(aes(x=ENERGY, ymin=OBS-PRIORUNC, ymax=OBS+PRIORUNC, col=EXPREF),
                           data=node_dt[grepl("^exp_",NODE)], alpha=0.4)
ggp <- ggp + geom_point(aes(x=ENERGY, y=OBS, col=EXPREF), data=node_dt[grepl("^exp_",NODE)])

# prior estimate (model pars + model defect)
ggp <- ggp + geom_line(aes(x=ENERGY, y=PRIORPRED), data=node_dt[xsel & grepl("^truexs_",NODE)], col="green")

ggp <- ggp + geom_ribbon(aes(x=ENERGY, ymin=PRED-PREDUNC, ymax=PRED+PREDUNC), data=node_dt[xsel & grepl("^truexs_",NODE)], alpha=0.2)
ggp <- ggp + geom_line(aes(x=ENERGY, y=PRED), data=node_dt[xsel & grepl("^truexs_",NODE)])

ggp <- ggp + theme(legend.position="none")
ggp <- ggp + facet_wrap(~ REAC, scales="free")
ggp <- ggp + xlab("energy [MeV]") + ylab("cross section [mbarn]")
# ggplotly(ggp)
ggp
# ggsave(paste0(savepath, "fe56-fast-truexs-eval.svg"), width=17, height=12, units="cm")

# plot the model defects
ggp <- ggplot() + theme_bw()
ggp <- ggp + geom_ribbon(aes(x=ENERGY, ymin=PRED-PREDUNC, ymax=PRED+PREDUNC), data=node_dt[xsel & grepl("^def_\\(N,P\\)$",NODE)], alpha=0.2)
ggp <- ggp + geom_line(aes(x=ENERGY, y=PRED), data=node_dt[xsel & grepl("^def_\\(N,P\\)$",NODE)])
ggp <- ggp + xlab("energy [MeV]") + ylab("relative model defect")
ggp <- ggp + scale_y_continuous(label = scales::percent_format(suffix="\\%"))
ggp
# ggsave(paste0(savepath, "fe56-fast-example-relative-defect-n-p.svg"), width=8.7, height=5, units="cm")

# plot the model parameters and their uncertainties
plotDt <- node_dt[grepl("^auxmodpar$", NODE)]
plotDt[, DATAMIN := pmax(0.9, PRED-PREDUNC)]
plotDt[, DATAMAX := pmin(1.1, PRED+PREDUNC)]
plotDt[, DIFF := abs(DATAMAX-DATAMIN)]
plotDt[, PRED := pmin(pmax(0.9, PRED), 1.1)]
plotDt <- plotDt[order(abs(DATAMAX-DATAMIN))]
# otherwise Latex complains:
plotDt[, PARNAME := gsub("_"," ",PARNAME)]
plotDt[, PARNAME := factor(PARNAME, levels=rev(PARNAME), ordered=TRUE)]


ggp <- ggplot(data=plotDt[1:10])
ggp <- ggp + theme_bw() #+ theme(axis.text.x = element_text(angle=90))
ggp <- ggp + theme(axis.text=element_text(size=9),
                   axis.title=element_text(size=10),
                   axis.title.x=element_text(hjust=1),
                   plot.title=element_text(size=12))
ggp <- ggp + xlab('parameter value relative to default') + ylab('')
ggp <- ggp + geom_errorbarh(aes(y=PARNAME, xmin=DATAMIN, xmax=DATAMAX), size=0.5, height=0.3)
ggp <- ggp + geom_point(aes(y=PARNAME, x=PRED), col='red', size=0.75)
ggp
# ggsave(paste0(savepath, "fe56-fast-modelpars-posterior.svg"), width=8.7, height=6, units="cm")
# ggsave(paste0(savepath, "fe56-fast-modelpars-posterior-more-converged.svg"), width=8.7, height=6, units="cm")


# print the table with the Bayesian network information
library(xtable)
print_dt <- unique(node_dt[, list(
  PRIOR,
  UNC,
  EMIN = min(ENERGY),
  EMAX = max(ENERGY),
  NUM = .N
), by=NODE])
print_dt <- print_dt[!grepl("^truexs_",NODE)]
print_dt <- print_dt[!grepl("^pred_",NODE)]
print_dt

# for latex output
print.xtable(xtable(print_dt,
                    display = c("s","s","g","g","f","f","d"),
                    align = c("l","l","|","c","c","r","r","r")),
             include.rownames=FALSE,
             math.style.exponents=TRUE)

# print table with experimental datasets for latex

expinfo <- list(
  "13132002" = list(
    citekey = "bowersAnalysisLonglivedIsotopes1989",
    author = "Bowers and Greenwood",
    year = 1989L
  ),
  "20091004" = list(
    citekey = "wenusch2nCrosssectionMeasurements1962",
    author = "Wenusch and Vonach",
    year= 1962L
  ),
  "20416003" = list(
    citekey = "frehautStatus2nCross1980",
    author = "Frehaut et al",
    year = 1980L
  ),
  "20721021" = list(
    citekey = "mollaSystematicStudyReactions1977",
    author = "Molla and Qaim",
    year = 1977L
  ),
  "20854015" = list(
    citekey = "corcalciucStudyNeutronInduced1978",
    author = "Corcalciuc",
    year = 1978L
  ),
  "23171003" = list(
    citekey = "wallnerProductionLonglivedRadionuclides2011",
    author = "Wallner et al",
    year = 2011L
  ),
  "12812012" = list(
    citekey = "sarafCrossSectionsSpectra1991",
    author = "Saraf et al",
    year = 1991L
  ),
  "32737002" = list(
    citekey = "wangCrossSectionsFe2015",
    author = "Wang et al",
    year = 2015L
  ),
  "10827031" = list(
    citekey = "grimesChargedparticleEmissionReactions1979",
    author = "Grimes et al",
    year = 1979L
  ),
  "10037004" = list(
    citekey = "boschungScatteringFastNeutrons1971",
    author = "Boschung et al",
    year = 1971L
  ),
  "10958012" = list(
    citekey = "el-kadiElasticInelasticScattering1982",
    author = "El-Khadi et al",
    year = 1982L
  ),
  "11708003" = list(
    citekey = "kinneyNeutronElasticInelastic1968",
    author = "Kinney",
    year = 1968L
  ),
  "14462004" = list(
    citekey = "ramirezNeutronScatteringCross2017",
    author = "Ramirez et al",
    year = 2017
  ),
  "23134005" = list(
    citekey = "beyerInelasticScatteringFast2014a",
    author = "Beyer et al",
    year = 2014L
  ),
  "30656021" = list(
    citekey = "xiaminMeasurementsInducedGamma1982",
    author = "Xiamin et al",
    year = 1982L
  ),
  "41156006" = list(
    citekey = "simakov14MeVFacilityResearch1992",
    author = "Simakov et al",
    year = 1992L
  ),
  "10022010" = list(
    citekey = "barrallHighEnergyNeutron1969",
    author = "Barrall et al",
    year = 1969L
  ),
  "10031005" = list(
    citekey = "barrallCrossSectionsReactions1969",
    author = "Barrall et al",
    year = 1969L
  ),
  "10289002" = list(
    citekey = "dyer56Fe58FeCross1972",
    author = "Dyer and Hamilton",
    year = 1972L
  ),
  "10309004" = list(
    citekey = "singhNeutronReactionCross1972",
    author = "Singh",
    year = 1972L
  ),
  "10417007" = list(
    citekey = "grundlStudyFissionNeutronSpectra1967",
    author = "Grundl",
    year = 1967L
  ),
  "10835002" = list(
    citekey = "sothrasStudySystematics2n1977",
    author = "Sothras",
    year = 1977L
  ),
  "11274031" = list(
    citekey = "paulCrossSectionMeasurements1953",
    author = "Paul and Clarke",
    year = 1953L
  ),
  "11464006" = list(
    citekey = "kernCrossSectionsReactions1959",
    author = "Thompson and Ferguson",
    year = 1959L
  ),
  "11494009" = list(
    citekey = "gabbardCrossSectionsCharged1962",
    author = "Gabbard and Kern",
    year = 1962L
  ),
  "11696007" = list(
    citekey = "crossActivationCrossSections1963",
    author = "Cross et al",
    year = 1963L
  ),
  "11701002" = list(
    citekey = "santryExcitationCurvesReactions1964",
    author = "Santry and Butler",
    year = 1964L
  ),
  "11701003" = list(
    citekey = "santryExcitationCurvesReactions1964",
    author = "Santry and Butler",
    year = 1964L
  ),
  "11701004" = list(
    citekey = "santryExcitationCurvesReactions1964",
    author = "Santry and Butler",
    year = 1964L
  ),
  "11703002" = list(
    citekey = "mcclureInelasticScattering141955",
    author = "Mcclure and Kent",
    year = 1955L
  ),
  "11715003" = list(
    citekey = "terrellExcitationFunctionFe1958",
    author = "Terrell and Holm",
    year = 1958L
  ),
  "11715004" = list(
    citekey = "terrellExcitationFunctionFe1958",
    author = "Terrell and Holm",
    year = 1958L
  ),
  "11718005" = list(
    citekey = "chittendenNewIsotopeManganese1961",
    author = "Chittenden et al",
    year = 1961L
  ),
  "12812010" = list(
    citekey = "sarafCrossSectionsSpectra1991",
    author = "Saraf et al",
    year = 1991L
  ),
  "12956012" = list(
    citekey = "spangler14MeVCrossSection1975",
    author = "Spangler et al",
    year = 1975L
  ),
  "12969013" = list(
    citekey = "meadowsMeasurement14MeV1987",
    author = "Meadows et al",
    year = 1987L
  ),
  "20280004" = list(
    citekey = "yasumiNuclearReactionsInduced1957",
    author = "Yasumi",
    year = 1957L
  ),
  "20377002" = list(
    citekey = "liskienCrosssectionMeasurementThreshold1965",
    author = "Liskien and Paulsen",
    year = 1965
  ),
  "20387004" = list(
    citekey = "liskienCrossSectionsCu631966",
    author = "Liskien and Paulsen",
    year = 1966L
  ),
  "20721094" = list(
    citekey = "mollaSystematicStudyReactions1977",
    author = "Molla and Qaim",
    year = 1977L
  ),
  "20772003" = list(
    citekey = "ryvesCrossSectionMeasurements1978",
    author = "Ryves et al",
    year = 1978L
  ),
  "20798003" = list(
    citekey = "robertson56Fe56Mn27Al1973",
    author = "Robertson et al",
    year = 1973L
  ),
  "20815014" = list(
    citekey = "vonachPrecisionMeasurementsExcitation1968",
    author = "Vonach et al",
    year = 1968L
  ),
  "20887015" = list(
    citekey = "bormannExcitationFunctionsNeutron1965",
    author = "Bormann et al",
    year = 1965L
  ),
  "20888004" = list(
    citekey = "bonazzolaMeasurementActivationCross1964",
    author = "Bonazzola",
    year = 1964L
  ),
  "20890004" = list(
    citekey = "cuzzocreaExcitationFunctionsNeutroninduced1968",
    author = "Cuzzocrea et al",
    year = 1968L
  ),
  "20993002" = list(
    citekey = "kudoAbsoluteMeasurement56Fe1977",
    author = "Kudo",
    year = 1977L
  ),
  "21049005" = list(
    citekey = "mostafaMeasurementsRelativeNeutron1976",
    author = "Mostafa",
    year = 1976L
  ),
  "21339005"= list(
    citekey = "bormannUeberWirkungsquerschnitteEiniger1962",
    author = "Bormann et al",
    year = 1962L
  ),
  "21352002" = list(
    citekey = "pollehnBestimmungWirkungsquerschnittenEiniger1961",
    author = "Pollehn and Neuert",
    year = 1961L
  ),
  "21352007" = list(
    citekey = "pollehnBestimmungWirkungsquerschnittenEiniger1961",
    author = "Pollehn and Neuert",
    year = 1961L
  ),
  "21372003" = list(
    citekey = "hemingwayDeterminationCrossSections1966",
    author = "Hemingway et al",
    year = 1966L
  ),
  "21419006" = list(
    citekey = "deprazMesureSectionsEfficaces1960",
    author = "Depraz et al",
    year = 1960L
  ),
  "21487008" = list(
    citekey = "allanProtonsInteraction141957",
    author = "Allan",
    year = 1957L
  ),
  "21868002" = list(
    citekey = "kudo56Fe56MnCross1982",
    author = "Kudo",
    year = 1982L
  ),
  "21923003" = list(
    citekey = NA,
    author = "Kudo",
    year = 1984L
  ),
  "22089042" = list(
    citekey = "ikedaActivationCrossSection1988",
    author = "Ikeda et al",
    year = 1988L
  ),
  "22093011" = list(
    citekey = "kimuraCalibratedFissionFusion1990",
    author = "Kimura and Kobayashi",
    year = 1990L
  ),
  "22312004" = list(
    citekey = "ikedaAbsoluteMeasurementsActivation1993",
    author = "Ikeda et al",
    year = 1993L
  ),
  "22338048" = list(
    citekey = "ercan14MeVNeutron1991",
    author = "Ercan et al",
    year = 1991L
  ),
  "22414017" = list(
    citekey = "fesslerNeutronActivationCrossSection2000",
    author = "Fessler et al",
    year = 2000L
  ),
  "22497004" = list(
    citekey = "coszachNeutroninducedReactionsContributing2000",
    author = "Coszach et al",
    year = 2000L
  ),
  "22976017" = list(
    citekey = "mannhartMeasurementNeutronActivation2007",
    author = "Mannhart and Schmidt",
    year = 2007L
  ),
  "30483002" = list(
    citekey = "chi-chouCrossSectionMeasurement1978",
    author = "Chi-Chou et al",
    year = 1978L
  ),
  "30483003" = list(
    citekey = "chi-chouCrossSectionMeasurement1978",
    author = "Chi-Chou et al",
    year = 1978L
  ),
  "30562019" = list(
    citekey = "ngocInvestigations2nReactions1980",
    author = "Ngoc et al",
    year = 1980L
  ),
  "30644006" = list(
    citekey = "viennotExcitationFunctionsReactions1982",
    author = "Viennot et al",
    year = 1982L
  ),
  "30676002" = list(
    citekey = "sharmaAbsoluteMeasurementsFe561978",
    author = "Sharma et al",
    year = 1978L
  ),
  "30707013" = list(
    citekey = "guptaPreequilibriumEmissionEffect1985",
    author = "Gupta et al",
    year = 1985L
  ),
  "30755003" = list(
    citekey = "muyaoShellEffectCross1987",
    author = "Muyao et al",
    year = 1987L
  ),
  "30802002" = list(
    citekey = "ngocNeutronActivationCross1983",
    author = "Ngoc et al",
    year = 1983L
  ),
  "30807008" = list(
    citekey = "garleaNeutronCrossSections1985",
    author = "Garlea et al",
    year = 1985L
  ),
  "30978022" = list(
    citekey = "viennotCrossSectionMeasurementsNp1991",
    author = "Viennot et al",
    year = 1991L
  ),
  "30993003" = list(
    citekey = "zongyuAbsoluteMeasurementCross1993",
    author = "Zongyu et al",
    year = 1993L
  ),
  "31479002" = list(
    citekey = "fugaStudyExcitationFunction1991",
    author = "Fuga",
    year = 1991L
  ),
  "31479003" = list(
    citekey = "fugaStudyExcitationFunction1991",
    author = "Fuga",
    year = 1991L
  ),
  "31524008" = list(
    citekey = "belgaidMeasurement14MeV1992",
    author = "Belgaid et al",
    year = 1992L
  ),
  "33045003" = list(
    citekey = "mulikMeasurement56Fe56Mn2013",
    author = "Mulik et al",
    year = 2013L
  ),
  "40485002" = list(
    citekey = "nemilovCrosssectionsReactionsNi581978",
    author = "Nemilov and Tofimov",
    year = 1978L
  ),
  "41118012" = list(
    citekey = "klochkovaInvestigationReactionsAl271992",
    author = "Klochkova et al",
    year = 1992L
  ),
  "41240012" = list(
    citekey = "filatenkovSystematicMeasurementActivation1999",
    author = "Filatenkov et al",
    year = 1999L
  ),
  "41313002" = list(
    citekey = "ramendikDetermination56Fe56Mn1977",
    author = "Ramendik et al",
    year = 1977L
  ),
  "41614019" = list(
    citekey = "filatenkovNeutronActivationCross2016",
    author = "Filatenkov",
    year = 2016L
  ),
  "20669004" = list(
    citekey = "qaimInvestigationReactions141976",
    author = "Qaim and Stocklin",
    year = 1976L
  ),
  "10037005" = list(
    citekey = "boschungScatteringFastNeutrons1971",
    author = "Boschung et al",
    year = 1971L
  ),
  "13764002" = list(
    citekey = NA,
    author = "Harvey",
    year = 1987
  ),
  "22316003" = list(
    citekey = NA,
    author = "Rohr et al",
    year = 1995L
  ),
  "41325003" = list(
    citekey = "tutubalinTotalNeutronCross1973",
    author = "Tutubalin et al",
    year = 1973L
  ),
  "41316002" = list(
    citekey = "kozyrRadiativeTransitionsUnbound1978",
    author = "Kozyr and Prokopets",
    year = 1978L
  ),
  "41118013" = list(
    citekey = "klochkovaInvestigationReactionsAl271992",
    author = "Klochkova",
    year = 1992L
  ),
  "32737003" = list(
    citekey = "wangCrossSections562015",
    author = "Wang et al",
    year = 2015L
  ),
  "20377002" = list(
    citekey = NA,
    author = "Liskien et al",
    year = 1965L
  ),
  "41240012" = list(
    citekey = "filatenkovSystematicMeasurementActivation1999",
    author = "Filatenkov et al",
    year = 1999
  )
)

print_dt <- node_dt[grepl("^exp_", NODE),
                    list(EMIN=min(ENERGY), EMAX=max(ENERGY)),
                    by=c("REAC","EXPREF")]
print_dt <- print_dt[order(REAC,EXPREF)]
print_dt[, AUTHOR := unlist(lapply(EXPREF, function(x) if (is.null(expinfo[[x]])) NA else expinfo[[x]]$author))]
print_dt[, YEAR := unlist(lapply(EXPREF, function(x) if (is.null(expinfo[[x]])) NA else expinfo[[x]]$year))]
print_dt[, YEAR := as.integer(YEAR)]
print_dt[, REF := unlist(lapply(expinfo[EXPREF], function(x) {
  if (!is.na(x$citekey)) paste0("\\cite{",x$citekey,"}") else ""
}))]
print_dt[, EMIN := NULL]
print_dt[, EMAX := NULL]
print_dt <- print_dt[order(REAC, YEAR)]
setnames(print_dt, "EXPREF", "EXFOR")

print_dt <- cbind(print_dt[1:45], print_dt[46:90])

library(xtable)
expinfo_xtable <- xtable(print_dt[1:45,],
  caption = paste0("Experimental datasets used in the example evaluation of neutron-induced reactionf of $^{56}$Fe between five and thirty MeV. ",
                   "The column EXFOR contains the EXFOR accession number. ",
                   "A missing reference means that no accessible publication is known to the authors."),
  label = "tbl:fe56-fast-datasets"
)

print.xtable(expinfo_xtable,
             sanitize.text.function = identity,
             include.rownames = FALSE)
