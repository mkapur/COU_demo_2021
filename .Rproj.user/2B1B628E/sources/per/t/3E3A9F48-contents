require(dplyr)
require(readr)
require(purrr)
# devtools::install_github("mkapur/kaputils", dependencies = F)
library(kaputils)
# devtools::install_github("r4ss/r4ss@2663227")
# devtools::install_github("r4ss/r4ss") ## post writeforecast merge
library(r4ss)

bg.statevals <- data.frame(matrix(NA, ncol = 3, nrow = 2))
colnames(bg.statevals) <- c('low','base','high')
row.names(bg.statevals) <- c('Fem','Mal')
bg.statevals$low <- c(0.046,0.048)
bg.statevals$base <- c(0.063,0.065)
bg.statevals$high <- c(0.086,0.089)
compname <- c('mkapur',"Maia Kapur")[1]

for(state in c('base','low','high')){
  
  rootdir.temp <- paste0("C:/Users/",compname,"/Dropbox/UW/assessments/blackgill-2019-update/ABC_base")
  if(state == 'base'){ ## only read in once per region
    catch_projections <- read.csv(paste0(rootdir.temp,"/blackgill_proj.csv"))
  }
  SS_autoForecast(rootdir = rootdir.temp,
                  basedir = "base_2015",
                  catch_proportions = catch_projections[catch_projections$YEAR==2021,5:ncol(catch_projections)],
                  state = state,
                  statesex = 2,
                  statevals = bg.statevals,
                  forecast_start = 2021,
                  forecast_end = 2031,
                  fixed_catches = catch_projections[catch_projections$YEAR<2021,5:ncol(catch_projections)],
                  Flimitfraction = catch_projections$PSTAR_0.45[catch_projections$YEAR >2020])
} ## end ABC runs

## Run Const/Upper catch x 3 states ----
## Did these separately because the constant nature of the catch 
rootdir <- paste0("C:/Users/",compname,"/Dropbox/UW/assessments/blackgill-2019-update/")
forecast_start <- 2021; forecast_end <- 2031; t = 10
for(catch in c('constant','upper')){ ## loop catch scen
  for(state in c('low','base','high')){
    
    df<-data.frame()
    catch_projections <- read.csv(paste0(rootdir,"/ABC_base/blackgill_proj.csv"))## from base dir
    Flimitfraction <- catch_projections$PSTAR_0.45[catch_projections$YEAR == 2030] ## doesn't really matter, already have catch vals
    catch_proportions <- catch_projections[catch_projections$YEAR == 2021,5:ncol(catch_projections)]
    const.catch <- mean(rowSums(catch_projections[catch_projections$TYPE == 'PROJECTION',5:ncol(catch_projections)])) ## avg 2019/2020
    fixed_catches <- catch_projections[catch_projections$TYPE == 'ACTUAL',5:ncol(catch_projections)]
    
    if(state != 'base'){
      lastrun <- paste0(rootdir,"ABC_",state)
    } else if(state == 'base'){
      lastrun <- paste0(rootdir,"ABC_",state,"/forecasts/forecast2030")
      
    }
    
    mod1 <- SS_output(lastrun, covar = FALSE, hidewarn = T, verbose = F) ## just load once for structure this hasn't executed yet
    
    newdir.temp <- paste0(rootdir,catch,"_",state)
    dir.create(newdir.temp) ## make special folder and copy files
    file.copy(list.files(lastrun,
                         full.names = TRUE,
                         recursive = TRUE),
              to = newdir.temp, overwrite = TRUE)
    setwd(newdir.temp) ## now forecast2030 appropriate to state is replicated here
    
    ## only need to change catches in forecast file
    fore <- SS_readforecast(file = './forecast.ss',
                            Nareas = mod1$nareas,
                            Nfleets = mod1$nfishfleets,
                            nseas = 1,
                            version = paste(mod1$SS_versionNumeric),
                            readAll = TRUE)
    
    
    fore$vals_fleet_relative_f <- paste(paste0(catch_proportions, collapse = " "))
    
    if(catch == 'constant'){
      ## apply 2019/2020 average to all yrs
      tempForeCatch <- SS_ForeCatch(mod1,
                                    yrs = 2021:(2021+(t-2)),
                                    average = FALSE,
                                    total = const.catch)
      fore$ForeCatch[min(which((fore$ForeCatch$Year==2021))):nrow(fore$ForeCatch),]<- tempForeCatch[,1:4]
      
      writecatch <- fore$ForeCatch %>% filter(Year > 2020) %>% group_by(Year) %>% dplyr::summarise(Catch_Used = sum(Catch_or_F))
      idx = nrow(writecatch)
      writecatch[idx+1,'Year'] <- 2030
      writecatch[idx+1,'Catch_Used'] <- const.catch
      write.csv(writecatch,
                file = "./tempForeCatch.csv",row.names = FALSE) ## save constant catched used       
      
    } else if (catch == 'upper'){
      ## apply 50% over 2021 to all years
      upperStream <- 1.5*mod1$derived_quants[grep("ForeCatch_2021", mod1$derived_quants$Label),"Value"]
      tempForeCatch <- SS_ForeCatch(mod1,
                                    yrs = 2021:(2021+(t-2)),
                                    average = FALSE,
                                    total = upperStream)
      
      fore$ForeCatch[min(which((fore$ForeCatch$Year==2021))):nrow(fore$ForeCatch),]<- tempForeCatch[,1:4]
      
      writecatch <- fore$ForeCatch %>% filter(Year > 2020) %>% group_by(Year) %>% dplyr::summarise(Catch_Used = sum(Catch_or_F))
      idx = nrow(writecatch)
      writecatch[idx+1,'Year'] <- 2030
      writecatch[idx+1,'Catch_Used'] <- upperStream
      write.csv(writecatch,
                file = "./tempForeCatch.csv",row.names = FALSE) ## save upperstream catched used
    }
    ## save file
    SS_writeforecast(fore, file = './forecast.ss', overwrite = TRUE)
    ## execute this model
    setwd(newdir.temp); system('ss3 -nohess') ## works
  } ## end states of nature
} ## end catch scenarios

## Build decision table (not in SS_executive summary) ----
YOI <- 2021:2030
dec_table <- matrix(NA, nrow = length(YOI)*3, ncol = 9)
dec_table <- data.frame(dec_table)
names(dec_table) <- c('Scenario','Year','catch',paste(c("spawnbio","depl"),rep(c('low','base','high'),each = 2)))
dec_table$Year <- rep(YOI,3)
idxr <- idxc <- 1
for(catch in c('constant','ABC','upper')){ ## loop catch scen
  idxc <- 1 ## reset to initial column for new catch scenario
  for(state in c('low','base','high')){
    if(catch != 'ABC'  | state != 'base'){
      tempdir <- paste0("C:/Users/",compname,"/Dropbox/UW/assessments/blackgill-2019-update/",catch,"_",state)
    } else if(catch == 'ABC'){
      tempdir <- paste0("C:/Users/",compname,"/Dropbox/UW/assessments/blackgill-2019-update/",catch,"_",state,"/forecasts/forecast2030")
    }
    mod <- SS_output(tempdir, covar = F)
    if(catch == 'constant' & idxc ==2){
      catch_projections <- read.csv(paste0("C:/Users/",compname,"/Dropbox/UW/assessments/blackgill-2019-update/ABC_base/blackgill_proj.csv"))
      const.catch <- mean(rowSums(catch_projections[catch_projections$TYPE == 'PROJECTION',5:ncol(catch_projections)])) ## avg 2019/2020
      dec_table$catch[idxr:(idxr+length(YOI)-1)] <- round(const.catch,2)
      
    } else if(catch == 'upper' & idxc ==2){
      upperStream <- read.csv(paste0(tempdir,"/tempforecatch.csv"))
      #1.5*mod$derived_quants[grep("ForeCatch_2021", mod$derived_quants$Label),"Value"]
      
      dec_table$catch[idxr:(idxr+length(YOI)-1)] <- round(upperStream[1,2],2)
      
    } else if (catch == 'ABC' &  idxc ==2){
      catchvals <- read.csv(paste0(rootdir,"/ABC_base/forecasts/forecast2030/tempforecatch.csv"))
      
      # mod$timeseries[, grepl('Yr|dead[(]B', names(mod$timeseries))] %>% 
      # filter(Yr %in% YOI) %>%
      # select(-Yr) %>% rowSums(.) %>% round(.,2)
      
      dec_table$catch[idxr:(idxr+length(YOI)-1)] <- round(catchvals$Catch_Used,0)
      
    }
    # read.csv(paste0(tempdir,"/tempForeCatch.csv"))
    
    
    ## input what was given to forecast file
    dec_table$Scenario[idxr:(idxr+length(YOI)-1)] <- rep(catch, length(idxr:(idxr+length(YOI)-1)))
    
    
    dec_table[idxr:(idxr+length(YOI)-1),idxc*2+2] <-  round_any(mod$derived_quants[grep(paste0("SSB_",YOI,collapse = "|"),
                                                                              mod$derived_quants$Label),"Value"],1000)/1000
    dec_table[idxr:(idxr+length(YOI)-1),idxc*2+3] <-  round(mod$derived_quants[grep(paste0("Bratio_",YOI,collapse = "|"),
                                                                              mod$derived_quants$Label),"Value"],2)
    idxc <- idxc+1 ## move to next set of columns as state updates
    # idxc <- idxc+3; idxr <-
    #     df["Depletion",y] <- paste0(round(basemod10$derived_quants[grep(paste0("Bratio_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"],3)*100,"%")
    
    
  } ## end state
  idxr <- idxr+length(YOI) ## jump down to next set of years when catch scenario updates
} ## end catch
## rename to look nice
dec_table$Scenario[dec_table$Scenario == 'constant'] <- c(rep(" ",5),'Constant (2019-2020 Average)',rep(" ",5))
dec_table$Scenario[dec_table$Scenario == 'ABC'] <- c(rep(" ",5),'40-10 Rule',rep(" ",5))
dec_table$Scenario[dec_table$Scenario == 'upper'] <- c(rep(" ",5),'Upper Stream',rep(" ",5))
## save dec_table
write.csv(dec_table, 
          file = paste0("C:/Users/",compname,"/Dropbox/UW/assessments/blackgill-2019-update/decision_table_",
                        Sys.Date(),".csv"),
          row.names = F)


## just generate the plots
bg.mod <- SS_output("C:/Users/",compname,"/Dropbox/UW/assessments/blackgill-2019-update/ABC_base/forecasts/forecast2021")
# save(bg.mod, file = paste0("C:/Users/",compname,"/Dropbox/UW/assessments/blackgill-2019-update/blackgill_2030_base.Rdata"))
View(bg.mod$derived_quants[grep(paste0("ForeCatch_",2021:2030,collapse = "|"), bg.mod$derived_quants$Label),"Value"]) ## start around 50
bg.mod$derived_quants[grep(paste0("OFLCatch_",2021:2030,collapse = "|"), bg.mod$derived_quants$Label),"Value"] ## start around 50

# View(bg.mod$timeseries[, grepl('Yr|dead[(]B', names(bg.mod$timeseries))] %>% filter(Yr %in% 2021:2030) %>% select(-Yr) %>% rowSums(.))
## check that 2017 spbio matches between models
bg.mod <- SS_output(paste0("C:/Users/",compname,"/Dropbox/UW/assessments/blackgill-2019-update/ABC_base/forecasts/forecast2030"))
bg.base <- SS_output(paste0("C:/Users/",compname,"/Dropbox/UW/assessments/blackgill-2019-update/base_2015 - Copy"))
# bg.mod$derived_quants[grep(paste0("SpBio_",2017,collapse = "|"), bg.mod$derived_quants$Label),"Value"] ## start around 50
bg.mod$timeseries$SpawnBio[bg.mod$timeseries$Yr == 2017] == bg.base$timeseries$SpawnBio[bg.base$timeseries$Yr == 2017]
bg.mod$timeseries$SpawnBio[bg.mod$timeseries$Yr == 2018] == bg.base$timeseries$SpawnBio[bg.base$timeseries$Yr == 2018]

bg.upbase <- SS_output(paste0("C:/Users/",compname,"/Dropbox/UW/assessments/blackgill-2019-update/upper_base"))
bg.upbase2 <- SS_output(paste0("C:/Users/",compname,"/Dropbox/UW/assessments/blackgill-2019-update/prediscard_0715/upper_base"))

bg.mod$timeseries$SpawnBio[bg.mod$timeseries$Yr == 2021]

bg.upbase$timeseries$SpawnBio[bg.upbase$timeseries$Yr == 2021]
bg.upbase2$timeseries$SpawnBio[bg.upbase2$timeseries$Yr == 2021]

SS_plots(bg.mod, png = T, dir = rootdir)


bg.mod.upper <- SS_output(paste0("C:/Users/",compname,"/Dropbox/UW/assessments/blackgill-2019-update/upper_base/"))
bg.mod.upper$timeseries[, grepl('Yr|dead[(]B', names(bg.mod.upper$timeseries))] %>% filter(Yr %in% 2021:2030) %>% select(-Yr) %>% rowSums(.)

iterOFL <- NULL; i = 1
for(y in 2021:2030){
  modNother <- SS_output(paste0("C:/Users/",compname,"/Dropbox/UW/assessments/blackgill-2019-update/ABC_base/forecasts/forecast",y))
  iterOFL[i] <- modNother$derived_quants[grep(paste0("OFLCatch_",y,collapse = "|"), modNother$derived_quants$Label),"Value"]
  i = i+1
}
View(iterOFL)


rd <- paste0("C:/Users/",compname,"/Dropbox/UW/assessments/blackgill-2019-update/ABC_base/forecasts/") ## a directory with all 9 model runs

iterOFL <- data.frame('MOD' = NA,'YEAR' = NA, 'OFL' = NA, 'FORECATCH' = NA, 
                      'DEADBIO' = NA,
                      'REALIZEDBUFFER' = NA,
                      'TRUEBUFFER_045' =  NA,
                      'TRUEBUFFER_025' = NA) ## sigma 45)
i <- 1
for(l in seq_along(list.dirs(rd, recursive = F))){
  modNother <- SS_output(list.dirs(rd, recursive = F)[l], covar = FALSE)
  modNother <-  SS_output(paste0("C:/Users/",compname,"/Dropbox/UW/assessments/blackgill-2019-update/ABC_high/"))
  modNother <-  SS_output(paste0("C:/Users/",compname,"/Dropbox/UW/assessments/blackgill-2019-update/ABC_base/forecasts/forecast2030/"))
  
  for(y in 2019:2030){
    # iterOFL[i,'MOD'] <- paste0(basename(list.dirs(rd, recursive = F)[l]))
    iterOFL[i,'YEAR'] <- y
    iterOFL[i,'OFL'] <- modNother$derived_quants[grep(paste0("OFLCatch_",y,collapse = "|"), modNother$derived_quants$Label),"Value"]
    iterOFL[i,'FORECATCH'] <- modNother$derived_quants[grep(paste0("ForeCatch_",y,collapse = "|"), modNother$derived_quants$Label),"Value"]
    iterOFL[i,'DEADBIO'] <-  modNother$timeseries[, grepl('Yr|dead[(]B', names(modNother$timeseries))] %>% filter(Yr == y) %>% select(-Yr) %>% rowSums(.)
    iterOFL[i,'REALIZEDBUFFER'] <-    round(iterOFL[i,'FORECATCH']/iterOFL[i,'OFL'],3)
    iterOFL[i,'TRUEBUFFER_045'] <-    c(NA,NA,0.857,
                                        0.849,
                                        0.841,
                                        0.833,
                                        0.826,
                                        0.818,
                                        0.810,
                                        0.803,
                                        0.795,
                                        0.788)[y-2018]#round(qlnorm(0.45,0,0.5*(1+c(1:10)*0.075)),3)[y-2020]
    iterOFL[i,'TRUEBUFFER_025'] <-   NA # round(qlnorm(0.25,0,0.5*(1+c(1:10)*0.075)),3)[y-2020]
    i <- i+1
  } ## end yrs
} ## end dirs
iterOFL ## th
