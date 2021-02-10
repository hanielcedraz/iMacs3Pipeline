#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("tools"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("-f", "--file"), type = "character", default = "samples.txt",
                help = "The filename of the sample file [default %default]",
                dest = "samplesFile"),
    make_option(c('-k', '--callpeak'), type = 'logical', default = TRUE, help = "Main MACS3 Function: Call peaks from alignment results. [default %default",
                dest = 'callpeak'),
    make_option(c("-c", "--column"), type = "character", default = "SAMPLE_ID",
                help = "Column name from the sample sheet to use as sample names [default %default]",
                dest = "column"),
    make_option(c('-b', '--broad'), action = 'store_true', type = 'logical',  default = FALSE, help = "If set, MACS will try to call broad peaks using the --broad-cutoff setting. [default %default", 
                dest = 'broad'),
    make_option(c('-C', '--broad-cutoff'), type = 'numeric', default = 0.05, help = "Cutoff for broad region. This option is not available unless --broad is set. [default %default", 
                dest = 'broadCutoff'),
    make_option(c("-i", "--inputFolder"), type = "character", default = "02-bamFiles",
                help = "Directory where the sequence data is stored [default %default]",
                dest = "inputFolder"),
    make_option(c('-o', '--output'), type = 'character', default = '03-callpeaks_',
                help = 'Folder where narrow output will be stored [default %default]',
                dest = 'outputFolder'),
    make_option(c('-g', '--gsize'), type = 'character', default = 'hs', help = "Effective genome size. It can be 1.0e+9 or 1000000000, or shortcuts:'hs' for human (2.7e9), 'mm' for mouse (1.87e9), 'ce' for C. elegans (9e7) and 'dm' for fruitfly (1.2e8) [default %default]",
                dest = 'gsize'),
    make_option(c('-B', '--bdg'), type = 'logical', default = TRUE, help = "Whether or not to save extended fragment pileup, and local lambda tracks (two files) at every bp into a bedGraph file [default %default]", 
                dest = 'bedGraph'),
    make_option(c("-p", "--processors"), type = "integer", default = 8,
                help = "Number of processors to use [defaults %default]",
                dest = "procs"),
    make_option(c("-a", "--minaqual"), type = "integer", default = 20,
                help = "Skip all reads with alignment quality lower than the given minimum value [defaults %default]",
                dest = "minaQual"),
    make_option(c('-q', '--qvalue'), type = 'numeric', default = 0.05, 
                help = "Minimum FDR (q-value) cutoff for peak detection. [default %default",
                dest = 'qValue'),
    make_option(c("-s", "--sampleprocs"), type = "integer", default = 2,
                help = "Number of samples to process at time [default %default]",
                dest = "mprocs"),
    make_option(c("-x", "--external"), type  =  "character", default = 'FALSE', help = "A space delimeted file with a single line contain several external parameters from MACS3 [default %default]",
                dest = "externalParameters")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list = option_list, description =  paste('Authors: OLIVEIRA, H.C.', 'Version: 0.1', 'E-mail: hanielcedraz@gmail.com', sep = "\n", collapse = '\n')))




######################################################################
## loadSampleFile
loadSamplesFile <- function(file, reads_folder, column){
    ## debug
    file = opt$samplesFile; reads_folder = opt$inputFolder; column = opt$samplesColumn
    ##
    if (!file.exists(file) ) {
        write(paste("Sample file",file,"does not exist\n"), stderr())
        stop()
    }
    ### column Control should be the control sample name
    ### rows can be commented out with #
    targets <- read.table(file, sep = "", header = TRUE, as.is = TRUE)
    if (!all(c("SAMPLE_ID", "Treated", "Input") %in% colnames(targets))) {
        cat('\n')
        write(paste("Expecting the four columns Control, Treated, FileControl, and FileTreated in samples file (tab-delimited)\n"), stderr())
        stop()
    }
    for (i in seq.int(nrow(targets$Control))) {
        if (targets[i,column]) {
            ext <- unique(file_ext(dir(file.path(reads_folder, targets[i,column]), pattern = "gz")))
            if (length(ext) == 0) {
                write(paste("Cannot locate fastq or sff file in folder",targets[i,column], "\n"), stderr())
                stop()
            }
            # targets$type[i] <- paste(ext,sep="/")
        }
        else {
            ext <- file_ext(grep("gz", dir(file.path(reads_folder,targets[i, column])), value = TRUE))
            if (length(ext) == 0) {
                write(paste(targets[i,column], "is not a gz file\n"), stderr())
                stop()
            }
            
        }
    }
    write(paste("samples sheet contains", nrow(targets), "samples to process", sep = " "),stdout())
    return(targets)
}


prepareCore <- function(opt_procs) {
    # if opt_procs set to 0 then expand to samples by targets
    if (detectCores() < opt$procs) {
        write(paste("number of cores specified (", opt$procs,") is greater than the number of cores available (",detectCores(),")",sep = " "),stdout())
        paste('Using ', detectCores(), 'threads')
    }
}


filesList <- function(samples, reads_folder, column){
    counting_list <- list()
    
    for (i in 1:nrow(samples)) {
        files <- dir(path = file.path(reads_folder), recursive = TRUE, pattern = paste0('.bam$'), full.names = TRUE)
        
        listFiles <- lapply(c("Input_", "treated_"), grep, x = files, value = TRUE)
        names(listFiles) <- c("Input", "treated")
        listFiles$sampleName <- samples[i, column]
        listFiles$Input <- listFiles$Input[i]
        listFiles$treated <- listFiles$treated[i]
        counting_list[[paste(listFiles$sampleName)]] <- listFiles
        counting_list[[paste(listFiles$sampleName, sep = "_")]]
        
    }
    
    write(paste("Setting up", length(counting_list), "jobs"),stdout())
    return(counting_list)
}


samples <- loadSamplesFile(opt$samplesFile, opt$inputFolder, opt$controlColumn)
procs <- prepareCore(opt$procs)
runningFiles <- filesList(samples, opt$inputFolder, opt$column)
#runningFiles

if (opt$broad) {
    outputFolder <- paste0(opt$outputFolder, "Broad")
    if (!file.exists(file.path(outputFolder))) dir.create(file.path(outputFolder), recursive = TRUE, showWarnings = FALSE)
} else {
    outputFolder <- paste0(opt$outputFolder, "Narrow")
    if (!file.exists(file.path(outputFolder))) dir.create(file.path(outputFolder), recursive = TRUE, showWarnings = FALSE)
}


external_parameters <- opt$externalParameters
if (file.exists(external_parameters)) {
    con = file(external_parameters, open = "r")
    line = readLines(con, warn = FALSE, ok = TRUE)
}


####################
### Running macs3
####################


if (opt$callpeak) {
    macs.run <- mclapply(runningFiles, function(index){
        try({
            system(paste('macs3',
                         'callpeak',
                         '-n',
                         ifelse(detectCores() < opt$procs, detectCores(), paste(opt$procs)),
                         '-t',
                         #"scrofaPig5_tested_sorted_aligned_RD_uniq.bam",
                         index$treated,
                         '-c',
                         #"scrofaPig5_control_sorted_aligned_RD_uniq.bam",
                         index$Input,
                         '-f',
                         'BAM',
                         '-g',
                         opt$gsize,
                         if (opt$broad) {
                             paste('--broad', '--broad-cutoff', opt$broadCutoff)
                             #paste()
                         } else {
                             paste('-q', opt$qValue, '-B')
                         },
                         '-n',
                         paste(outputFolder, index$sampleName, sep = "/")
                         
                         
            ))})
    }, mc.cores = opt$mprocs
    )
    
    
    if (!all(sapply(macs.run, "==", 0L))) {
        write(paste("Something went wrong with macs3. Some jobs failed"),stderr())
        stop()
    } else {
        write(paste('All jobs finished successfully'), stderr())
    }
} else {
    if (file.exists(external_parameters)) {
        macs.run <- mclapply(runningFiles, function(index){
            try({
                system(paste('macs3',
                             line)
                       
                )
                
            })
        }, mc.cores = opt$mprocs
        )
        
        
        if (!all(sapply(macs.run, "==", 0L))) {
            write(paste("Something went wrong with macs3. Some jobs failed"),stderr())
            stop()
        }else{
            write(paste('All jobs finished successfully'), stderr())
        }
    } else {
        write(paste("You did not give any parameter for running macs3. Please provide the file containing the parameters for analysis"), stderr())
    }
}
