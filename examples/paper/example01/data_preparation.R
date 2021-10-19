################################################################################
#                  Data preparation
################################################################################

rawFilePath <- system.file("extdata", "U5nf_data.c5m", package = "nucdataBaynet")

tableFormat <- rbind.data.frame(
  c("Proj", "I5"), c("Target", "I6"), c(NA, "X"), c("MF", "I3"), c("MT", "I4"),
  c(NA, "X3"), c("Energy", "F9"), c("dEnergy", "F9"), c("Data", "F9"),
  c("dData", "F9"), c(NA, "X39"), c("Reference", "A25"), c("Entry", "F5"),
  c("SubP", "F3"), c(NA, "X"), c("dSys", "F9"), c("dStat", "F9"), c("dOther", "F9"),
  c("dTot", "F9"), c("dSysP", "F9"), c("dStatP", "F9"), c("dOtherP", "F9"),
  c("dTotP", "F9"), c("dDataP", "F9")
)
names(tableFormat) <- c("varname", "format")

rawdt <- as.data.table(read.fortran(rawFilePath, tableFormat$format,
                                    col.names=with(tableFormat, varname[!is.na(varname)])))

rawdt[Energy > 7000 & Energy < 12000, Reference]
rawdt[["Reference"]] <- gsub("[, ]+"," ",rawdt[["Reference"]])
