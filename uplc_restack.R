#-----------------------------------------------#
#' #Function for Lindroth Lab uplc data formatting
##Created by Clay Morrow 2017-March-21
##Version 1.2 Updated: 2017-June-28
#-----------------------------------------------#

### Source this file to add uplc.restack() to your library.
# This function will restack raw outptut .txt file(s) from the uplc into either a horizontally
# or vertically stacked data frame and write the output to a .csv file.

##-----Arguments-----##

### The functon uplc.restack() has 10 arguments, 2 of which are necessary:

#---required arguments---#

## 1. 'path.in' is the path to the raw .txt file or of a folder containing only raw .txt
# files. This path can be absolute or relative to your working directory i.e. "~/folder"

## 2.'path.out' is the path to the folder in which the restacked.csv output will be saved. if
# no such directory exists, it will be created unless supress.dir.create=TRUE

#---optional arguments---#

## 3.'multiple' (default: TRUE) is a logical argument indicating whether 1 or more files are to be
# restacked (single file or a folder). If set to FALSE, the 'path.in' must be a single
# .txt file; if set to TRUE (default) then 'path.in' must be a folder containing only
# raw .txt output files.

## 4. 'wide' (default: TRUE) is a logical argument indicating whether the output is to be horizontally
# stacked (wide). The default is TRUE which will yield one instance of each sample
# in a row and add seperate variables (columns) for the data associated with each compound.
# if FALSE, then the observations will be stacked as rows multiple times with single
# readouts for each column and a single new grouping variable indicating the compound
# being measured in that row.

## 5. 'repeated.cols' is a vector containing the columns to be kept for each compound (i.e. repeated).
# The default is: c("Area", "Response", "Conc").

#---Master file arguments---#

## 6.. 'master' (default: FALSE) is a logical argument indicating whether a master file should
# be created in addition to the restacked files. A column for batch# will be added in the
# master file, to distinguish each batch. These batch numbers must be part of the original
# file names to work properly.
## NOTE: if master=TRUE, arguments 7 and 8 will also be used, so the file naming formats may be
# limiting.

## 7. 'name.delim' (default:"_") is a character object that separates words in the a file names
# for example "W16_BATCH10_2017.txt" is seperated by underscores into "W16" "BATCH10" "2017.txt"
# this step is important for extracting and labeling batch numbers

## 8. 'batch.word' (default: "BATCH") is a character object that indicates the batch identifer.
# for example in the file "W16_BATCH10_2017.txt" the batch word is "BATCH"

##9. 'batch.pos' (default: 2) is an integer that indicates which position of the file name contains
# the batch information. for example, in the file "W16_BATCH10_2017.txt", the batch information
# is in the second position.
## Note: arguments 8 and 9 allow the function to extract "10" as the batch number for this file
# the master file will then have 10 in the batch column for all observations from this file.


#---File Creation---#

## 6. 'supress.dir.create' (default:FALSE) is a logical argument indicating whether to supress the creation of a
# folder in the path.out directory. Even, when this value set to FALSE, a new directory
# is only created if it does not already exist.

##----------------------------------------------------------------------------------------------##

uplc.restack(path.in="~/Google Drive/Competiton Garden/Data/2017/Comp 2017 PGs Raw Data/Text Files/",path.out="~/Google Drive/Competiton Garden/Data/Comp 2017 PGs Raw Data/", multiple=T)

##-------Arguments----------##

# path.in <- "Z:/Lindroth Lab/WisAsp Study/Phytochemistry/2016_WisAsp_PGs_Master/June2016_WIsAsp_Raw_PG_Data/morrowcj/UPLC TXT June 2016 files/"
# path.out <- "C:/Users/morro/Desktop/New folder/"
#
# multiple <- TRUE
# wide <- TRUE
# repeated.cols <-  c("Area","Response","Conc") #c(8,10,12)
#
# master <- TRUE
# name.delim <- "_"
# batch.word <- "BATCH"
# batch.pos <- 2
#
# supress.dir.create <- TRUE

##-----------Function----------##
uplc.restack <- function(path.in,
                         path.out,
                         multiple = TRUE,
                         wide = TRUE,
                         repeated.cols = c("Area", "Response", "Conc"),
                         master = FALSE,
                         name.delim = "_",
                         batch.word = "BATCH",
                         batch.pos = 2,
                         supress.dir.create = FALSE) {
  ##Required Libraries
  require(dplyr, quietly = TRUE)
  require(plyr, quietly = TRUE)

  ##output directory check and create
  if (supress.dir.create != TRUE) {
    if (dir.exists(path.out) == FALSE) {
      dir.create(path.out)
      message("directory \"", path.out, "\" created")
    }
  }
  ##Files to use
  if (multiple == TRUE) {
    Files <- list.files(path.in)
  } else if (multiple == FALSE) {
    Files <- basename(path.in)
  } else {
    stop("multiple must be TRUE/FALSE")
  }

  ##iterate over each file
  for (j in 1:length(Files)) {
    #----------------------#

    ## absolute path to each file
    if (multiple == TRUE) {
      file.location <- paste(path.in, Files[j], sep = "")
    } else if (multiple == FALSE) {
      file.location <- path.in
    }

    #---------Data loading------------#
    ## This format is specific to UPLC output
    Data <- read.table(
      file.location,
      header = FALSE,
      sep = c("\t", "\n"),
      skip = 3,
      blank.lines.skip = TRUE,
      stringsAsFactors = FALSE,
      fill = TRUE,
      na.strings = c("NA", ""),
      # These column names expect that the raw output also contains these columns.
      col.names = c(
        "Obs",
        "num",
        "ID",
        "Name",
        "Type",
        "Std.Conc",
        "RT",
        "Area",
        "IS.Area",
        "Response",
        "Detection.flags",
        "Conc",
        "perc.dev"
      )
    )


    # These are the columns that get added for each compound
    cols.to.add <-
      repeated.cols #c(8,10,12) #c("Area","Response","Conc")

    ##Remove Blank Rows
    Data <- Data[which(!is.na(Data$Obs)),]

    ##--------Name Compounds------------##
    #These are the compounds
    Compound.key <- grep("Compound*", Data$Obs, value = TRUE)
    #This is which row they are in
    Compound.index <- c(grep("Compound*", Data$Obs))
    #This takes the name after the colon
    Compound.names <-
      unlist(strsplit(Compound.key, ": "))[seq(2, length(Compound.key) * 2, by = 2)]
    ## Tells which compound is being analyzed in each row
    Compound.ID <- NULL
    for (i in 1:length(Compound.key)) {
      Compound.ID <-
        c(Compound.ID, rep(
          Compound.names[i],
          length(Compound.index[1]:Compound.index[2]) - 1
        ))
    }
    # Creates a data frame with an added column for the observation's compound
    Data.long <- cbind(Data, Compound.ID)
    # removes the subheaders that divide the raw frames
    Data.long <- Data.long[-Compound.index,]

    ##----------write long form----------##
    # This is what the new file will be called
    new.file.name <- sub(".txt", "_restacked.csv", Files[j])
    # only write these if wide=FALSE
    if (wide == FALSE) {
      write.csv(Data.long, paste(path.out, new.file.name, sep = ""))
    } else if (wide == TRUE) {
      ##---------Wide form----------------##

      # counts the samples run (number of unique samples)
      n <- length(unique(na.omit(Data$num)))
      # builds foundation (of first set of each observation and common columns only )
      Base <- Data.long[1:n,]
      Base <-
        subset(Base, select = c(Obs, num, ID, Name, Type, IS.Area))
      #column names for the first 13 columns
      Base.names <-
        names(Base)#c(names(Base)[1:7],paste(Compound.names[1],names(Base)[8:13],sep = "."))
      # an index of which rows sample #1 occurs in (and also, one row after the end of the file)
      first.index <-
        c(seq(1, nrow(Data.long), by = n), nrow(Data.long) + 1)

      # Creating a new data frame by adding
      New.Data <- Base
      # x=1
      for (z in 1:length(Compound.key)) {
        New.Data <-
          cbind(New.Data, Data.long[first.index[z]:(first.index[z + 1] - 1), cols.to.add])
      }

      ##new names
      comps.to.use <- Compound.names[1:length(Compound.names)]

      names(New.Data) <-
        c(Base.names, paste(rep(comps.to.use, each = length(cols.to.add)), names(New.Data)[cols.to.add], sep = "_"))


      ## File Output
      # new.file.name <- sub(".txt","_restacked.csv",Files[j])

      write.csv(New.Data, paste(path.out, new.file.name, sep = ""))

    } else {
      stop("'wide' must be TRUE/FALSE")
    }
  }
  message("Your file(s) have been restacked. Check your 'path.out' directory to find them.")

  #------Master List-------#
  ## This will concatenate all restacked files into a master copy with a batch column.
  if (master == TRUE) {
    #creates a list of all restacked files from above
    CAT.files <- grep("*.csv", list.files(path.out), value = TRUE)
    #splits each file name by the delimeter
    split.files <- strsplit(CAT.files, name.delim)
    #empty vectors to fill
    Batch <- NULL
    MASTER <- NULL

    # files.to.use <- grep("*.csv",list.files(path.out))
    for (i in 1:length(CAT.files)) {
      #
      Batch <-
        c(Batch, na.exclude(as.numeric(unlist(
          strsplit(split.files[[i]][batch.pos], batch.word, fixed = FALSE)
        )))[1])

      data <-
        read.csv(file = paste(path.out, CAT.files[i], sep = ""),
                 header = TRUE)

      data <- cbind(Batch[i], data)

      #----only if you want each batch column
      # new.name <- paste("withbatchcolumn",CAT.files[i],sep = "_")
      #
      # write.csv(x = data,file = paste(path.out2,new.name))

      # binds rows together even if column numbers are different (assumes you used the same order of compounds each time)
      MASTER <- rbind.fill(MASTER, data)

    }
    ##-----write master list-------##

    if (dir.exists(paste(path.out, "MASTER", sep = "")) == FALSE) {
      dir.create(paste(path.out, "MASTER", sep = ""))
      message("'MASTER/' directory created")
    }

    write.csv(MASTER, paste(path.out, "MASTER/", "MASTER.csv", sep = ""))
    message(
      "Your Master file has been created. Check your 'path.out'/Master directory to find it."
    )
  }
}
