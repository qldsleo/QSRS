libpath <- '/export/home/leos/libs'; .libPaths(libpath)
library(dplyr); library(ggplot2); library(tidyr); library(plotly)
######################################################################################################################
# use for extracting
batch_script_extracting <- function(){
    dir.create(paste0('/export/home/leos/QSRS/code/batch_scripts/extracting/', attribute), showWarnings=F)
    if(!extract_now){stop('check set up')}
    cpus = 5; mem = 100; walltime = '1:00:00'; cores = 5
    lines <- c('#!/bin/bash',
               paste0('#PBS -l select=1:ncpus=',cpus,':mem=',mem,'gb'),
               paste0('#PBS -l walltime=',walltime),
               paste0('#PBS -e /export/home/leos/QSRS/code/batch_scripts/extracting/',attribute,'/error_log.log'),
               paste0('#PBS -o /export/home/leos/QSRS/code/batch_scripts/extracting/',attribute,'/output_log.log'),
               paste0('#PBS -N ',attribute),
               paste0('singularity exec docker://artemis:5050/qld/des/rsc/rserver_4.2.1:latest Rscript /export/home/leos/QSRS/code/QSRS/setup.R ',cores))
    
    # write the lines to a file
    writeLines(lines, paste0('code/batch_scripts/extracting/',attribute,'/QSRS.sh'))
    system(paste0('qsub QSRS/code/batch_scripts/extracting/',attribute,'/QSRS.sh'))
}
batch_script_extracting()
readLines(paste0('/export/home/leos/QSRS/code/batch_scripts/extracting/',attribute,'/error_log.log'))

######################################################################################################################
# use for modelling
batch_script_modelling <- function(){
    dir.create(paste0('/export/home/leos/QSRS/code/batch_scripts/modelling/', attribute), showWarnings=F)
    if(!xv_now){stop('check set up')}
    cpus = 30; mem = 20; walltime = '0:30:00'; cores = 3
    lines <- c('#!/bin/bash',
               paste0('#PBS -l select=1:ncpus=',cpus,':mem=',mem,'gb'),
               paste0('#PBS -l walltime=',walltime),
               paste0('#PBS -e /export/home/leos/QSRS/code/batch_scripts/modelling/',attribute,'/error_log.log'),
               paste0('#PBS -o /export/home/leos/QSRS/code/batch_scripts/modelling/',attribute,'/output_log.log'),
               paste0('#PBS -N ',attribute),
               paste0('singularity exec docker://artemis:5050/qld/des/rsc/rserver_4.2.1:latest Rscript /export/home/leos/QSRS/code/QSRS/setup.R ',cores))
    
    # write the lines to a file
    writeLines(lines, paste0('code/batch_scripts/modelling/',attribute,'/QSRS.sh'))
    system(paste0('qsub QSRS/code/batch_scripts/modelling/',attribute,'/QSRS.sh'))
}
batch_script_modelling()
readLines(paste0('/export/home/leos/QSRS/code/batch_scripts/modelling/',attribute,'/error_log.log'))

######################################################################################################################
# use for mapping
cpu_mem_model <- function(chunks, dep, var, attribute, version){
    cpu_mem_model <- NULL
    for (chunk in chunks){
        print(chunk)
        check_file <- paste0('code/batch_scripts/mapping/',attribute,'/v',version,'/error_log',chunk,'-',dep,'-',var,'.log')
        if (file.exists(check_file)){
            open_file <- readLines(con = check_file)
            cpupercent <- as.numeric(gsub('[^0-9]', '', unlist(strsplit(open_file[grep('cpupercent', open_file)], ','))[5]))
            mem <- as.numeric(gsub('[^0-9]', '', unlist(strsplit(open_file[grep('cpupercent', open_file)], ','))[2]))/1000000
            walltime <- as.numeric(gsub('[^0-9]', '', unlist(strsplit(open_file[grep('cpupercent', open_file)], ','))[3]))
            ncpus <- as.numeric(gsub('[^0-9]', '', unlist(strsplit(open_file[grep('cpupercent', open_file)], ','))[4]))
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
cpu_mem_data <- cpu_mem_model(chunks=1:200, dep=1, var=50, attribute=attribute, version=version)
hpc_model <- pivot_longer(cpu_mem_data, cols=c('cpupercent','mem','walltime'), names_to='metric') %>% as.data.frame
ggplotly(ggplot(hpc_model, aes(chunk, value))+
             geom_point()+
             facet_wrap(vars(metric), nrow=2, scales='free')+
             geom_smooth(method = 'gam', formula = y ~ s(x, k=4)))
#_____________________________________________________________________________________________________________________
batch_script_mapping <- function(variance, depth, chunk){
    
    dir.create(paste0('/export/home/leos/QSRS/code/batch_scripts/mapping/', attribute), showWarnings=F)
    dir.create(paste0('/export/home/leos/QSRS/code/batch_scripts/mapping/', attribute, '/v', version), showWarnings=F)
    if(!map_now){stop('check set up')}
    for (var in variance){
        for (dep in depth){
            for (chunk in chunk){
                print(paste0('depth: ',dep, ' | chunk: ',chunk, ' | variance: ',var))
                
                check_file_exists <- paste0('/export/home/leos/QSRS/code/batch_scripts/mapping/',attribute,'/v',version,'/error_log',chunk,'-',dep,'-',var,'.log')
                if(file.exists(check_file_exists)){
                    check_file <- readLines(paste0('/export/home/leos/QSRS/code/batch_scripts/mapping/',attribute,'/v',version,'/error_log',chunk,'-',dep,'-',var,'.log'))
                    error_status <- check_file[grepl('ExitStatus', check_file)]
                    error_status <- as.numeric(gsub('ExitStatus:', '', error_status))
                } else {error_status <- 1}
                
                if (file.exists(check_file_exists) & error_status == 0){}else{
                    
                    if (var == 50){
                        if (between(chunk, 1, 80)){
                            cpus = 5; cores = 5; mem = 20; walltime = '0:10:00'
                        } else if(chunk > 80){
                            cpus = 10; cores = 5; mem = 20; walltime = '0:10:00'
                        }
                    } else if (var != 50){
                        if (between(chunk, 1, 75)){
                            cpus = 5; cores = 5; mem = 20; walltime = '0:20:00'
                        } else if(chunk > 70){
                            cpus = 5; cores = 2; mem = 20; walltime = '1:00:00' # try reducing cores
                            # cpus = 10; cores = 5; mem = 50; walltime = '0:30:00' # try reducing cores
                        }
                    }
                    
                    lines <- c('#!/bin/bash',
                               paste0('#PBS -l select=1:ncpus=',cpus,':mem=',mem,'gb'),
                               paste0('#PBS -l walltime=',walltime),
                               paste0('#PBS -e /export/home/leos/QSRS/code/batch_scripts/mapping/',attribute,'/v',version,'/error_log',chunk,'-',dep,'-',var,'.log'),
                               paste0('#PBS -o /export/home/leos/QSRS/code/batch_scripts/mapping/',attribute,'/v',version,'/output_log',chunk,'-',dep,'-',var,'.log'),
                               #'#PBS -a 1145',
                               paste0('#PBS -N ',attribute,'-',chunk,'-',dep,'-',var),
                               paste0('singularity exec docker://artemis:5050/qld/des/rsc/rserver_4.2.1:latest Rscript /export/home/leos/QSRS/code/QSRS/setup.R ',chunk,' ',cores,' ',dep, ' ',var))
                    
                    writeLines(lines, paste0('code/batch_scripts/mapping/',attribute,'/v',version,'/QSRS',chunk,'-',dep,'-',var,'.sh'))
                    system(paste0('qsub QSRS/code/batch_scripts/mapping/',attribute,'/v',version,'/QSRS',chunk,'-',dep,'-',var,'.sh'))
                    Sys.sleep(1)
                }
            }
        }
    }
}
batch_script_mapping(variance = 50, depth = 1, chunk = seq(1,200,by=10))

######################################################################################################################
# use for mosaicing
batch_script_mosaicing <- function(variance, depth){
    
    dir.create(paste0('/export/home/leos/QSRS/code/batch_scripts/mosaicing/', attribute), showWarnings=F)
    if(!mosaic_now){stop('check set up')}
    for (var in variance){
        for (dep in depth){
            print(paste0('depth: ',dep, ' | variance: ',var))
            
            cpus = 1; mem = 20; walltime = '1:00:00'
            
            lines <- c('#!/bin/bash',
                       paste0('#PBS -l select=1:ncpus=',cpus,':mem=',mem,'gb'),
                       paste0('#PBS -l walltime=',walltime),
                       paste0('#PBS -e /export/home/leos/QSRS/code/batch_scripts/mosaicing/',attribute,'/error_log',dep,'-',var,'.log'),
                       paste0('#PBS -o /export/home/leos/QSRS/code/batch_scripts/mosaicing/',attribute,'/output_log',dep,'-',var,'.log'),
                       paste0('#PBS -N ',attribute,'-',dep,'-',var),
                       paste0('singularity exec docker://artemis:5050/qld/des/rsc/rserver_4.2.1:latest Rscript /export/home/leos/QSRS/code/QSRS/setup.R ',dep,' ',var))
            
            # write the lines to a file
            writeLines(lines, paste0('code/batch_scripts/mosaicing/',attribute,'/QSRS',dep,'-',var,'.sh'))
            system(paste0('qsub QSRS/code/batch_scripts/mosaicing/',attribute,'/QSRS',dep,'-',var,'.sh'))
            Sys.sleep(1)
        }
    }
}
batch_script_mosaicing(variance = c(50), depth = 1)
