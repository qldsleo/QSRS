libpath <- '/scratch/rsc3/leos/R/libs'; .libPaths(libpath)
library(dplyr); library(ggplot2); library(tidyr); library(plotly)
email <- 'stephen.leo@des.qld.gov.au'

# clay, silt, sand, fs, cs, cat_ca, cat_mg, cat_na, cat_k, cat_cec, tn, tp, clay_act, esp, sar (wb_oc - not done)

######################################################################################################################
# use for extracting
batch_script_extracting <- function(attribute){
    
    for (attr in attribute){
        dir.create(paste0('/scratch/rsc3/leos/QSRS/batch_scripts/extracting/', attr), showWarnings=F)
        cpus = 5; mem = 200; walltime = '6:00:00'; cores = 5
        lines <- c('#!/bin/bash',
                   paste0('#PBS -l select=1:ncpus=',cpus,':mem=',mem,'gb'),
                   paste0('#PBS -l walltime=',walltime),
                   paste0('#PBS -e /scratch/rsc3/leos/QSRS/batch_scripts/extracting/',attr,'/error_log.log'),
                   paste0('#PBS -o /scratch/rsc3/leos/QSRS/batch_scripts/extracting/',attr,'/output_log.log'),
                   paste0('#PBS -N ',attr),
                   paste0('#PBS -M ', email),
                   paste0('#PBS -m ae'),
                   paste0('singularity exec docker://artemis:5050/qld/des/rsc/rserver_4.2.1:latest Rscript /export/home/leos/QSRS/code/QSRS/setup.R ',
                          attr, ' ', 1, ' ', 'extract', ' ', cores))
        
        # write the lines to a file
        writeLines(lines, paste0('/scratch/rsc3/leos/QSRS/batch_scripts/extracting/',attr,'/QSRS.sh'))
        system(paste0('qsub /scratch/rsc3/leos/QSRS/batch_scripts/extracting/',attr,'/QSRS.sh'))
        Sys.sleep(1)
    }
}
batch_script_extracting(attribute = 'col_p')

######################################################################################################################
# use for modelling
batch_script_modelling <- function(attribute, version, depth, cpu){
    
    for (attr in attribute){
        for (vers in version){
            for (dep in depth){
                dir.create(paste0('/scratch/rsc3/leos/QSRS/batch_scripts/modelling/', attr), showWarnings=F)
                cpus = cpu; mem = 20; walltime = '0:30:00'
                
                lines <- c('#!/bin/bash',
                           paste0('#PBS -l select=1:ncpus=',cpus,':mem=',mem,'gb'),
                           paste0('#PBS -l walltime=',walltime),
                           paste0('#PBS -e /scratch/rsc3/leos/QSRS/batch_scripts/modelling/',attr,'/error_log',dep,'.log'),
                           paste0('#PBS -o /scratch/rsc3/leos/QSRS/batch_scripts/modelling/',attr,'/output_log',dep,'.log'),
                           paste0('#PBS -N ',attr,'-',dep),
                           paste0('#PBS -M ', email),
                           paste0('#PBS -m ',ifelse(dep==1, 'ae', 'a')),
                           paste0('singularity exec docker://artemis:5050/qld/des/rsc/rserver_4.2.1:latest Rscript /export/home/leos/QSRS/code/QSRS/setup.R ',
                                  attr, ' ', vers, ' ', 'model', ' ', dep))
                
                # write the lines to a file
                writeLines(lines, paste0('/scratch/rsc3/leos/QSRS/batch_scripts/modelling/',attr,'/QSRS',dep,'-',vers,'.sh'))
                system(paste0('qsub /scratch/rsc3/leos/QSRS/batch_scripts/modelling/',attr,'/QSRS',dep,'-',vers,'.sh'))
                Sys.sleep(1)
            }
        }
    }
}
batch_script_modelling(attribute = c('phw'), version = 5, depth = 1:6, cpu = 15)

######################################################################################################################
# use for model summarising
batch_script_model_sum <- function(attribute, version){
    
    for (attr in attribute){
        for (vers in version){
            dir.create(paste0('/scratch/rsc3/leos/QSRS/batch_scripts/modelling/', attr), showWarnings=F)
            cpus = 15; mem = 20; walltime = '0:30:00'; cores = 5
            
            lines <- c('#!/bin/bash',
                       paste0('#PBS -l select=1:ncpus=',cpus,':mem=',mem,'gb'),
                       paste0('#PBS -l walltime=',walltime),
                       paste0('#PBS -e /scratch/rsc3/leos/QSRS/batch_scripts/modelling/',attr,'/error_log.log'),
                       paste0('#PBS -o /scratch/rsc3/leos/QSRS/batch_scripts/modelling/',attr,'/output_log.log'),
                       paste0('#PBS -N ',attr),
                       paste0('#PBS -M ', email),
                       paste0('#PBS -m ae'),
                       paste0('singularity exec docker://artemis:5050/qld/des/rsc/rserver_4.2.1:latest Rscript /export/home/leos/QSRS/code/QSRS/setup.R ',
                              attr, ' ', vers, ' ', 'model_sum', ' ', cores))
            
            # write the lines to a file
            writeLines(lines, paste0('/scratch/rsc3/leos/QSRS/batch_scripts/modelling/',attr,'/QSRS',vers,'.sh'))
            system(paste0('qsub /scratch/rsc3/leos/QSRS/batch_scripts/modelling/',attr,'/QSRS',vers,'.sh'))
            Sys.sleep(1)
        }
    }
}
batch_script_model_sum(attribute = c('phw'), version = 5)

######################################################################################################################
# use for mapping
cpu_mem_model <- function(attribute, version, chunks, depth, var){
    cpu_mem_model <- NULL
    for (chunk in chunks){
        print(chunk)
        check_file <- paste0('/scratch/rsc3/leos/QSRS/batch_scripts/mapping/',attribute,'/v',version,'/error_log',chunk,'-',depth,'-',var,'.log')
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
cpu_mem_data <- cpu_mem_model(attribute = 'phw', version = 5, chunks = 1:200, depth = 1, var = 50)
hpc_model <- pivot_longer(cpu_mem_data, cols=c('cpupercent','mem','walltime'), names_to='metric') %>% as.data.frame
ggplotly(ggplot(hpc_model, aes(chunk, value))+
             geom_point()+
             facet_wrap(vars(metric), nrow=2, scales='free')+
             geom_smooth(method = 'gam', formula = y ~ s(x, k=4)))
#_____________________________________________________________________________________________________________________
batch_script_mapping <- function(attribute, version, variance, depth, chunk){
    
    for (attr in attribute){
        for (vers in version){
            dir.create(paste0('/scratch/rsc3/leos/QSRS/batch_scripts/mapping/', attr), showWarnings=F)
            dir.create(paste0('/scratch/rsc3/leos/QSRS/batch_scripts/mapping/', attr, '/v', vers), showWarnings=F)
            for (var in variance){
                for (dep in depth){
                    for (chunk in chunk){
                        print(paste0('depth: ',dep, ' | chunk: ',chunk, ' | variance: ',var))
                        
                        check_file_exists <- paste0('/scratch/rsc3/leos/QSRS/batch_scripts/mapping/',attr,'/v',vers,'/error_log',chunk,'-',dep,'-',var,'.log')
                        if(file.exists(check_file_exists)){
                            check_file <- readLines(paste0('/scratch/rsc3/leos/QSRS/batch_scripts/mapping/',attr,'/v',vers,'/error_log',chunk,'-',dep,'-',var,'.log'))
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
                                       paste0('#PBS -e /scratch/rsc3/leos/QSRS/batch_scripts/mapping/',attr,'/v',vers,'/error_log',chunk,'-',dep,'-',var,'.log'),
                                       paste0('#PBS -o /scratch/rsc3/leos/QSRS/batch_scripts/mapping/',attr,'/v',vers,'/output_log',chunk,'-',dep,'-',var,'.log'),
                                       #'#PBS -a 1145',
                                       paste0('#PBS -N ',attr,'-',chunk,'-',dep,'-',var),
                                       paste0('singularity exec docker://artemis:5050/qld/des/rsc/rserver_4.2.1:latest Rscript /export/home/leos/QSRS/code/QSRS/setup.R ',
                                              attr, ' ', vers, ' ', 'map', ' ', chunk,' ', cores,' ', dep, ' ', var))
                            
                            writeLines(lines, paste0('/scratch/rsc3/leos/QSRS/batch_scripts/mapping/',attr,'/v',vers,'/QSRS',chunk,'-',dep,'-',var,'.sh'))
                            system(paste0('qsub /scratch/rsc3/leos/QSRS/batch_scripts/mapping/',attr,'/v',vers,'/QSRS',chunk,'-',dep,'-',var,'.sh'))
                            Sys.sleep(1)
                        }
                    }
                }
            }
        }
    }
}
batch_script_mapping(attribute = 'phw', version = 5, variance = 50, depth = 1, chunk = 2)

######################################################################################################################
# use for mosaicing
batch_script_mosaicing <- function(attribute, version, variance, depth){
    
    for (attr in attribute){
        for (vers in version){
            dir.create(paste0('/scratch/rsc3/leos/QSRS/batch_scripts/mosaicing/', attr), showWarnings=F)
            for (var in variance){
                for (dep in depth){
                    print(paste0('depth: ',dep, ' | variance: ',var))
                    
                    cpus = 1; mem = 20; walltime = '1:00:00'
                    
                    lines <- c('#!/bin/bash',
                               paste0('#PBS -l select=1:ncpus=',cpus,':mem=',mem,'gb'),
                               paste0('#PBS -l walltime=',walltime),
                               paste0('#PBS -e /scratch/rsc3/leos/QSRS/batch_scripts/mosaicing/',attr,'/error_log',dep,'-',var,'.log'),
                               paste0('#PBS -o /scratch/rsc3/leos/QSRS/batch_scripts/mosaicing/',attr,'/output_log',dep,'-',var,'.log'),
                               paste0('#PBS -N ',attr,'-',dep,'-',var),
                               paste0('singularity exec docker://artemis:5050/qld/des/rsc/rserver_4.2.1:latest Rscript /export/home/leos/QSRS/code/QSRS/setup.R ',
                                      attr, ' ', vers, ' ', 'mosiac', ' ', dep,' ',var))
                    
                    # write the lines to a file
                    writeLines(lines, paste0('/scratch/rsc3/leos/QSRS/batch_scripts/mosaicing/',attr,'/QSRS',dep,'-',var,'.sh'))
                    system(paste0('qsub /scratch/rsc3/leos/QSRS/batch_scripts/mosaicing/',attr,'/QSRS',dep,'-',var,'.sh'))
                    Sys.sleep(1)
                }
            }
        }
    }
}
batch_script_mosaicing(attribute = 'ec', version = 3, variance = c(50), depth = 1)
