## ------------------------------------------------------------------------
library(ggpubr)
library(ggplot2)
library(readr)
library(tidyr)
library(hexbin)

library(xlsx)



## ------------------------------------------------------------------------
defaultInput = "C:/Users/sebas/Desktop/testAnalysis/data/"
txtList <- list("cAng_avg.txt",
                "cAng_ampl.txt",
                "cAng_avg_abs.txt",
                "cAng_f_avg_f1.txt",
                "cAng_f_avg_f2.txt",
                "y_avg.txt",
                "y_ampl.txt",
                "y_avg_abs.txt",
                "z_avg.txt",
                "z_ampl.txt",
                "z_avg_abs.txt",
                "z_f_avg_f1.txt",
                "z_f_avg_f2.txt")
message("Insert raw data Path: (press space for default path)")
inputPath = readline("")
if(inputPath == "") inputPath = defaultInput
setwd(inputPath)
message(paste0("processing files at:   ",getwd()))

dir <- "data"
nonWTLabel <- "LAPD"
outputFolder <- c(paste0(file.path(Sys.getenv("USERPROFILE"),"Desktop"), "/plots_SpermQ_", gsub(":","_",date()), "/"))

xLim = c(0,90) #shown part of the flagellum in ?m
zoom = 2.2

dir.create(outputFolder, showWarnings = T)

input <- data.frame(ID = integer(),
                    path = character(),
                    stringsAsFactors = F,
                    date = character(),
                    mouse = character(),
                    exp = character(),
                    type = character(),
                    moreThan1Analysis = logical()
)
ID = 0
dateList <- list.files()
i = 0
for(i in 1:length(dateList)){
  dDir = paste0(dateList[i])
  miceList = list.files(dDir)
  j = 0
  for(j in 1:length(miceList)){
    mDir = paste0(dDir,"/", miceList[j])
    experimentList = list.files(mDir)
    h = 0
    for(h in 1:length(experimentList)){
      expDir = paste0(mDir, "/", experimentList[h])
      path = paste0(expDir, "/", list.files(expDir)[length(list.files(expDir))], "/")
      ID = ID +1
      input[ID,] = list(ID, path, dateList[i], miceList[j], experimentList[h], "WT",F)
      if(!grepl("WT",input$mouse[ID])){
        input$type[ID] = nonWTLabel
      }
      if(length(list.files(expDir))>1){input$moreThan1Analysis[ID] = T}
    }
  }
}
input$type = factor(input$type, levels = c("WT", nonWTLabel))
rownames(input) = input$ID
rm(h,i,j)
rm(dDir,mDir,expDir,experimentList,ID,path, miceList)
write.xlsx(input, paste0(outputFolder, "expList(input).xlsx"))

path_activity <- paste0(outputFolder,"activityProfiles","/")
dir.create(path_activity, showWarnings = F)

## ------------------------------------------------------------------------
saveActivityProfiles <- function(rawList, txtName = "txtName", path = "Desktop"){
  if(path == "Desktop")path = paste0(file.path(Sys.getenv("USERPROFILE"),"Desktop"),"/")
  data = data.frame(ID = numeric(),
                    time = factor(levels = c("before", "after", "1min", "3min", "5min")),
                    average = numeric(),
                    type = factor(levels = c("WT", "LAPD")))
  i = 0
  for(i in 1:length(rawList)){
    data[nrow(data) + 1,] = list(rawList[[i]]$ID[1], "before", mean(rawList[[i]]$before, na.rm = T), rawList[[i]]$type[1])
    data[nrow(data) + 1,] = list(rawList[[i]]$ID[1], "after", mean(rawList[[i]]$after, na.rm = T), rawList[[i]]$type[1])
    data[nrow(data) + 1,] = list(rawList[[i]]$ID[1], "1min", mean(rawList[[i]]$'1min', na.rm = T), rawList[[i]]$type[1])
    data[nrow(data) + 1,] = list(rawList[[i]]$ID[1], "3min", mean(rawList[[i]]$'3min', na.rm = T), rawList[[i]]$type[1])
    data[nrow(data) + 1,] = list(rawList[[i]]$ID[1], "5min", mean(rawList[[i]]$'5min', na.rm = T), rawList[[i]]$type[1])
  }
  svg(paste0(path, gsub(".txt","",txtName),"_singleLines", ".svg"))
  lines <- ggplot(data, aes(x = time, y = average, colour = type)) + geom_line(aes(group = ID)) + theme_minimal() + facet_wrap(~ID) 
  #+ labs(title = paste0("activity profile (based on ",txtName,") for each sperm"), x = "time", y = "activity")
  print(lines)
  dev.off()
  
  svg(paste0(path_activity, gsub(".txt","",txtName),"_boxplots", ".svg"))
  boxplot = ggboxplot(raw, x = "time", y = "value", color = "type") + theme_minimal() +
    ggtitle(paste0("activity profile (based on ",txtName,") for each sperm"))
  print(boxplot)
  dev.off()
  save(boxplot, lines, file = paste0(path, txtName, ".RData")) 
  return()
}

getLinePlot <- function(data, type = "_raw"){
  name = paste0(gsub(".txt","",txtName),type, "_line")
  svg(paste0(path, name, ".svg"))
  plot <- ggplot(raw, aes(x=arcLength, y=value, colour = type)) + geom_line(aes(group = ID))+ facet_grid(time ~type) +
    theme_minimal() + xlab("arc Length (?m)") + xlim(xLim) + ylim(yLim) +ggtitle(name)
  plots[["line"]] <<- plot
  print(plot)
  dev.off()
}

getBoxplot <- function(data, type = "_raw"){
  name = paste0(gsub(".txt","",txtName), type,"_boxplot")
  svg(paste0(path, name, ".svg"))
  plot <- ggboxplot(raw, x="bin", y="value", color = "type", add = "boxplot") + theme_minimal() + xlab("arc Length (?m)") + ggtitle(name)
  plots[["boxplot"]] <<- plot
  print(plot)
  dev.off()
  
  name = paste0(gsub(".txt","",txtName),type, "_boxplot_facetWrap")
  svg(paste0(path, name, ".svg"))
  print(plot + facet_wrap(~time, ncol = 1) )
  dev.off()
}

getDensityPlot <- function(data, type = "_raw"){
  name = paste0(gsub(".txt","",txtName),type, "_density")
  svg(paste0(path, name, ".svg"))
  plot <- ggplot(raw, aes(x=arcLength, y=value)) + facet_grid(time~type) + xlim(xLim) + theme_minimal() + 
    geom_hex(bins = 20) +  scale_fill_distiller(palette= "Spectral", direction=1) + ggtitle(name)
  plots[["density"]] <<- plot
  print(plot)
  dev.off()
}

getFrame <- function(data){
  wtMean = mean(subset(data, data$type == "WT")$value, na.rm = T)
  nWtMean = mean(subset(data, data$type != "WT")$value, na.rm = T)
  wtSd= sd(subset(data, data$type == "WT")$value, na.rm = T)
  nWtSd = sd(subset(data, data$type != "WT")$value, na.rm = T)
  if(wtSd > nWtSd){
    sd = wtSd
  }else{
    sd = nWtSd
  }
  if(wtMean > nWtMean){
    yLim <<- c(nWtMean-zoom*sd, wtMean+zoom*sd)
  }else{
    yLim <<- c(wtMean-zoom*sd, nWtMean+zoom*sd)
  }
  if(all(data$value>=0,na.rm = T)) yLim[1] <<- 0
}

notToRemove = list()
notToRemove = ls()

## ------------------------------------------------------------------------

message("---------- parameter in process:--------------------")


for(task in 1:length(txtList)){
  txtName <- txtList[[task]]
  message(paste(task, "/", length(txtList), "-",txtList[[task]]))
  path <- paste0(outputFolder,gsub(".txt","",txtName), "/")
  dir.create(path, showWarnings = F)
  
  
  ## ------------------------------------------------------------------------
  rawList <- vector("list", 1)
  i = 0
  for(i in 1:length(input$path)){
    dataFrame <- read_delim(paste0(input$path[i],txtName), 
                            "\t", escape_double = FALSE, locale = locale(), 
                            na = "empty", trim_ws = TRUE, skip = 5, col_types = cols())
    colnames(dataFrame) = c("arcLength","before", "after", "1min", "3min", "5min")
    dataFrame$arcLength = factor(dataFrame$arcLength, levels = unique(dataFrame$arcLength))
    dataFrame$ID = rep(input$ID[i], nrow(dataFrame))
    dataFrame$type = rep(input$type[i], nrow(dataFrame))
    rawList[[i]] <- dataFrame
  }
  names(rawList) <- input$ID
  rm(dataFrame)
  
  ## ------------------------------------------------------------------------
  i = 0
  raw = data.frame()
  for(i in 1:length(rawList)){
    raw <- rbind(raw,gather(rawList[[i]], key = "time", value = "value", before:'5min', factor_key=F, na.rm = F, convert = T))
  }
  raw$time = factor(raw$time, levels = c("before", "after", "1min", "3min", "5min"))
  raw$ID = factor(raw$ID, levels = unique(raw$ID))
  raw$arcLength = as.numeric(as.character(raw$arcLength))
  raw$bin <- cut(raw$arcLength, c(-1,10,30,50,70,90,110,120,200))
  raw$bin <- factor(raw$bin, levels = as.character(unique(raw$bin)))
  
  ## ------------------------------------------------------------------------
  getFrame(raw)
  
  plots <- list()
  
  saveActivityProfiles(rawList, txtName = txtName, path = path_activity)
  
  getLinePlot(data = raw, type = "_raw")
  
  getBoxplot(data = raw, type  = "_raw")
  
  getDensityPlot(data = raw, type = "_raw")
  
  save(file = paste0(path, gsub(".txt", "", txtName),"_raw.RData"),plots)
  
  ## ------------------------------------------------------------------------
  i = 0
  subtracted = data.frame()
  subList <- rawList
  for(i in 1:length(subList)){
    subList[[i]]$after = subList[[i]]$after - subList[[i]]$before
    subList[[i]]$'1min' = subList[[i]]$'1min' - subList[[i]]$before
    subList[[i]]$'3min' = subList[[i]]$'3min' - subList[[i]]$before
    subList[[i]]$'5min' = subList[[i]]$'5min' - subList[[i]]$before
    subList[[i]]$before <- NULL
    subtracted <- rbind(subtracted,gather(subList[[i]], key = "time", value = "value", after:'5min', factor_key=F, na.rm = F, convert = T))
  }
  subtracted$ID = factor(subtracted$ID, levels = unique(subtracted$ID))
  subtracted$arcLength = as.numeric(as.character(subtracted$arcLength))
  subtracted$bin <- cut(subtracted$arcLength, c(-1,10,30,50,70,90,110,120,200))
  subtracted$bin <- factor(subtracted$bin, levels = as.character(unique(subtracted$bin)))
  subtracted$time = factor(as.character(subtracted$time), levels = c("after", "1min", "3min", "5min"))
  
  ## ------------------------------------------------------------------------
  getFrame(subtracted)
  
  plots <- list()
  
  getLinePlot(data = subtracted, type = "_subtracted")
  
  getBoxplot(data = subtracted, type  = "_subtracted")
  
  getDensityPlot(data = subtracted, type = "_subtracted")
  
  save(file = paste0(path, gsub(".txt", "", txtName),"_subtracted.RData"),plots)
  
  ## ------------------------------------------------------------------------
  
  i = 0
  divided = data.frame()
  devList <- rawList
  for(i in 1:length(devList)){
    devList[[i]]$after = devList[[i]]$after / devList[[i]]$before
    devList[[i]]$'1min' = devList[[i]]$'1min' / devList[[i]]$before
    devList[[i]]$'3min' = devList[[i]]$'3min' / devList[[i]]$before
    devList[[i]]$'5min' = devList[[i]]$'5min' / devList[[i]]$before
    devList[[i]]$before = NULL
    divided <- rbind(divided,gather(devList[[i]], key = "time", value = "value", after:'5min', factor_key=F, na.rm = F, convert = T))
  }
  divided$ID = factor(divided$ID, levels = unique(divided$ID))
  divided$arcLength = as.numeric(as.character(divided$arcLength))
  divided$bin <- cut(divided$arcLength, c(-1,10,30,50,70,90,110,120))
  divided$bin <- factor(divided$bin, levels = as.character(unique(divided$bin)))
  divided$time = factor(as.character(divided$time), levels = c("after", "1min", "3min", "5min"))
  
  ## ------------------------------------------------------------------------
  plots <- list()
  
  getLinePlot(data = divided, type = "_divided")
  
  getBoxplot(data = divided, type  = "_divided")
  
  getDensityPlot(data = divided, type = "_divided")
  
  save(file = paste0(path, gsub(".txt", "", txtName),"_divided.RData"),plots)
  
  
  ## ------------------------------------------------------------------------
  
  wideTable <- spread(raw, time, value)
  wideTable <- wideTable[order(as.numeric(wideTable$ID)),]
  write.xlsx(wideTable, paste0(path, gsub(".txt", "", txtName),"_rawData.xlsx"))
  
}

