######################################################################################################################
# use for mapping
libpath <- '/export/home/leos/libs'; .libPaths(libpath)
library(dplyr); library(ggplot2); library(tidyr); library(plotly)
cpu_mem_model <- function(chunks, dep, var, attribute, version){
    cpu_mem_model <- NULL
    for (chunk in chunks){
        print(chunk)
        check_file <- paste0("code/batch_scripts/mapping/",attribute,'/v',version,"/error_log",chunk,'-',dep,'-',var,".log")
        if (file.exists(check_file)){
            open_file <- readLines(con = check_file)
            cpupercent <- as.numeric(gsub("[^0-9]", "", unlist(strsplit(open_file[grep("cpupercent", open_file)], ","))[5]))
            mem <- as.numeric(gsub("[^0-9]", "", unlist(strsplit(open_file[grep("cpupercent", open_file)], ","))[2]))/1000000
            walltime <- as.numeric(gsub("[^0-9]", "", unlist(strsplit(open_file[grep("cpupercent", open_file)], ","))[3]))
            ncpus <- as.numeric(gsub("[^0-9]", "", unlist(strsplit(open_file[grep("cpupercent", open_file)], ","))[4]))
            hpc_usage <- data.frame('chunk' = chunk,
                                    'cpupercent' = cpupercent/100,
                                    'mem' = mem,
                                    'walltime' = walltime,
                                    'ncpus' = ncpus)
            cpu_mem_model <- rbind(cpu_mem_model, hpc_usage)
        }
    }
    return(cpu_mem_model)
}

cpu_mem_data <- cpu_mem_model(chunks=1:200, dep=1, var=5, attribute=attribute, version=version)
cpu_mem_data
hpc_model <- pivot_longer(cpu_mem_data, cols=c('cpupercent','mem','walltime'), names_to='metric') %>% as.data.frame
p1 <- ggplot(hpc_model, aes(chunk, value))+
    geom_point()+
    facet_wrap(vars(metric), nrow=2, scales='free')+
    geom_smooth(method = "gam", formula = y ~ s(x, k = 4))
ggplotly(p1)
cpu_est <- mgcv::gam(cpupercent ~ s(chunk,k=4), data=cpu_mem_data)

#####################

process <- "mapping"
dir.create(paste0('/export/home/leos/QSRS/code/batch_scripts/', process, '/', attribute), showWarnings=F)
dir.create(paste0('/export/home/leos/QSRS/code/batch_scripts/', process, '/', attribute, '/v', version), showWarnings=F)
if(!map_now){stop('check set up')}
variance = 5; depth = 1; chunk = 1:50

for (var in variance){
    for (dep in depth){
        for (chunk in chunk){
            print(paste0('depth: ',dep, ' | chunk: ',chunk, ' | variance: ',var))
            
            check_file <- paste0("code/batch_scripts/mapping/",attribute,'/v',version,"/error_log",chunk,'-',dep,'-',var,".log")
            if (file.exists(check_file)){}else{
                
                if (var == 50){
                    if (between(chunk, 1, 80)){
                        cpus = 5; cores = 5; mem = 20; walltime = "0:10:00"
                    } else if(chunk > 80){
                        cpus = 10; cores = 5; mem = 20; walltime = "0:10:00"
                    }
                } else if (var != 50){
                    if (between(chunk, 1, 70)){
                        cpus = 5; cores = 5; mem = 20; walltime = "0:10:00"
                    } else if(chunk > 70){
                        cpus = 10; cores = 5; mem = 50; walltime = "0:30:00" # try reducing cores
                    }
                }
                
                lines <- c("#!/bin/bash",
                           paste0("#PBS -l select=1:ncpus=",cpus,":mem=",mem,"gb"),
                           paste0("#PBS -l walltime=",walltime),
                           paste0("#PBS -e /export/home/leos/QSRS/code/batch_scripts/",process,"/",attribute,'/v',version,"/error_log",chunk,'-',dep,'-',var,".log"),
                           paste0("#PBS -o /export/home/leos/QSRS/code/batch_scripts/",process,"/",attribute,'/v',version,"/output_log",chunk,'-',dep,'-',var,".log"),
                           #"#PBS -a 1145",
                           paste0("#PBS -N QSRS",chunk,'-',dep,'-',var),
                           paste0("singularity exec docker://artemis:5050/qld/des/rsc/rserver_4.2.1:latest Rscript /export/home/leos/QSRS/code/setup.R ",chunk," ",cores," ",dep, " ",var))
                
                writeLines(lines, paste0("code/batch_scripts/",process,"/",attribute,'/v',version,"/QSRS",chunk,'-',dep,'-',var,".sh"))
                system(paste0("qsub QSRS/code/batch_scripts/",process,"/",attribute,'/v',version,"/QSRS",chunk,'-',dep,'-',var,".sh"))
                Sys.sleep(1)
            }
        }
    }
}

######################################################################################################################
# use for extracting and modelling
cpus <- 30
mem <- 20
walltime <- "0:30:00"
process <- "modelling" # extracting, modelling
cores <- 3
dir.create(paste0('/export/home/leos/QSRS/code/batch_scripts/', process, '/', attribute), showWarnings=F)
lines <- c("#!/bin/bash",
           paste0("#PBS -l select=1:ncpus=",cpus,":mem=",mem,"gb"),
           paste0("#PBS -l walltime=",walltime),
           paste0("#PBS -e /export/home/leos/QSRS/code/batch_scripts/",process,"/",attribute,"/error_log.log"),
           paste0("#PBS -o /export/home/leos/QSRS/code/batch_scripts/",process,"/",attribute,"/output_log.log"),
           paste0("#PBS -N QSRS"),
           paste0("singularity exec docker://artemis:5050/qld/des/rsc/rserver_4.2.1:latest Rscript /export/home/leos/QSRS/code/setup.R ",cores))

# write the lines to a file
writeLines(lines, paste0("code/batch_scripts/",process,"/",attribute,"/QSRS.sh"))
system(paste0("qsub QSRS/code/batch_scripts/",process,"/",attribute,"/QSRS.sh"))

######################################################################################################################
# use for mosaicing
cpus <- 1
mem <- 20
walltime <- "2:00:00"
process <- "mosaicing"
dir.create(paste0('/export/home/leos/QSRS/code/batch_scripts/', process, '/', attribute), showWarnings=F)

for (dep in 1:6){
    print(dep)
    lines <- c("#!/bin/bash",
               paste0("#PBS -l select=1:ncpus=",cpus,":mem=",mem,"gb"),
               paste0("#PBS -l walltime=",walltime),
               paste0("#PBS -e /export/home/leos/QSRS/code/batch_scripts/",process,"/",attribute,"/error_log",dep,".log"),
               paste0("#PBS -o /export/home/leos/QSRS/code/batch_scripts/",process,"/",attribute,"/output_log",dep,".log"),
               paste0("#PBS -N QSRS",dep),
               paste0("singularity exec docker://artemis:5050/qld/des/rsc/rserver_4.2.1:latest Rscript /export/home/leos/QSRS/code/setup.R ",dep))
    
    # write the lines to a file
    writeLines(lines, paste0("code/batch_scripts/",process,"/",attribute,"/QSRS",dep,".sh"))
    system(paste0("qsub QSRS/code/batch_scripts/",process,"/",attribute,"/QSRS",dep,".sh"))
    Sys.sleep(1)
}