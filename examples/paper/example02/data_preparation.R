################################################################################
#                  Data preparation
################################################################################

rawdata <- copy(n_fe56_xsdata)
rawdata[grepl("TOT", REAC), REAC := "TOT"]
rawdata[grepl("EL", REAC), REAC := "EL"]
rawdata[grepl("INL", REAC), REAC := "INL"]
