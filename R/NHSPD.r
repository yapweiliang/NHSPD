# NHSPD.r
#' NHS Postcode Directory
#'
#' download/process NHS Postcode Directory into a form I want; reload locally saved RDS data.
#'
#' @section Background:
#' The NHS Postcode Directory is a specially tailored dataset based on ONS Postcode Directory.
#'
#' In particular, each postcode is mapped to a LSOA, so this enables linking to the
#' Index of Multiple Deprivation.
#'
#' Download from:
#' \itemize{
#'  \item \url{http://systems.hscic.gov.uk/data/ods/datadownloads/onsdata}
#'  \item \url{http://systems.hscic.gov.uk/data/ods/datadownloads/onsdata/zip-files/gridall.zip}
#' }
#'
#' This library therefore supercedes my \pkg{codepoint} library.
#'
#' on comparing NHSPD data (November 2015) with OS CodePoint Open 2015.4.0 (October 2015),
#' at least the England postcodes with grid coordinates are identical to OS CodePoint
#' except for 4 postcodes where NHSPD data seems correct if plotted on map
#' (TR210AA, TR210PR, TR210PW, TR210WY labelled as Isles of Scilly on CodePoint vs Cornwall on NHSPD)
#'
#' @section TODO:
#' TODO: error checking for CodePoint stuff
#' @docType package
#' @name NHSPD
NULL

.citations <- function(year = 2015) {
  s1 <- paste0("Contains Ordnance Survey data (C) Crown copyright and database right ", year)
  s2 <- paste0("Contains Royal Mail data (C) Royal Mail copyright and database right ", year)
  s3 <- paste0("Contains National Statistics data (C) Crown copyright and database right ", year)
  return(sprintf("The NHSPD package process data that needs the following citations:\n  %s\n  %s\n  %s\n",s1,s2,s3))
}

# default source data files -----------------------------------------------

.default_gridall_DataFile <- "H:/DATASETS/HSCIC/gridall.csv" # local copy of download
.default_gridall_URL <- "http://systems.hscic.gov.uk/data/ods/datadownloads/onsdata/zip-files/gridall.zip"
.default_gridall_RDS <- "H:/DATASETS/HSCIC/gridall.RDS"

.default_CodePoint_Unzipped_Path <- "H:/DATASETS/OS/CodePoint/tmp"

# constants ---------------------------------------------------------------

.headers.short <- c(
  "PCD2",       #	Unit postcode – 8 character version
  "PCDS",       #	Unit postcode - variable length (e-Gif) version
  "DOINTR",     #	Date of introduction
  "DOTERM",     #	Date of termination
  "OSEAST100M", #	National grid reference - Easting
  "OSNRTH100M", #	National grid reference - Northing
  "OSCTY",      #	County
  "ODSLAUA",    #	Local Authority Organisation
  "OSLAUA",     #	Local authority district (LAD)/unitary authority (UA)/ metropolitan district (MD)/ London borough (LB)/ council area (CA)/district council area (DCA)
  "OSWARD",     #	(Electoral) ward/division
  "USERTYPE",   #	Postcode user type
  "OSGRDIND",   #	Grid reference positional quality indicator
  "CTRY",       #	Country
  "OSHLTHAU",   #	Former Strategic health authority (SHA)/ health board (HB)/ health authority (HA)/ health & social care board (HSCB)
  "GOR",        #	Region
  "OLDHA",      #	Pre-2002 Health Authority
  "NHSCR",      #	Commissioning Region
  "CCG",        #	Clinical Commissioning Group (CCG)/ local health board (LHB)/ community health partnership (CHP)/ local commission-ing group (LCG)/ primary healthcare directorate (PHD)
  "PSED",       #	1991 Census Enumeration District (ED)
  "CENED",      #	1991 Census Enumeration District (ED)
  "EDIND",      #	ED positional quality indicator
  "WARD98",     #	1998 ward
  "OA01",       #	2001 Census Output Area
  "NHSRG",      #	NHS Region (Geography), formerly Area Team
  "HRO",        #	Former Pan SHA
  "LSOA01",     #	2001 Lower Layer Super Output Area (LSOA)
  "UR01IND",    #	2001 Census Urban/rural indicator
  "MSOA01",     #	2001 Census Middle Layer Super Output Area (MSOA)
  "CANNET",     #	Former Cancer Network
  "SCN",        #	Strategic Clinical Network
  "OSHAPREV",   #	‘First wave’ Strategic health authority (SHA)/ health board (HB)/ health authority (HA)/ health and social services board (HSSB)
  "OLDPCT",     #	‘First wave’ Primary Care Trust (PCT)/ Local Health Board (LHB)/ Care Trust (CT)
  "OLDHRO",     #	‘Old’ IT Cluster (ITC)
  "PCON",       #	Westminster parliamentary constituency
  "CANREG",     #	Cancer Registry
  "PCT",        #	‘Second wave’ Primary Care Trust (PCT)/ Care Trust/ Care Trust Plus (CT)/ local health board (LHB)/ community health partnership (CHP)/ local commission-ing group (LCG)/ primary healthcare directorate (PHD)
  "OSEAST1M",   #	National grid reference - Easting
  "OSNRTH1M",   #	National grid reference - Northing
  "OA11",       #	2011 Census output area
  "LSOA11",     #	2011Census  lower layer super output area (LSOA)
  "MSOA11"      #	2011 Census middle layer super output area (MSOA)
)

.headers.wanted <- c(
  "PCD2", #	Unit postcode – 8 character version
  "DOTERM", #	Date of termination
  "OSCTY", #	County
#  "ODSLAUA", #	Local Authority Organisation
  "OSLAUA", #	Local authority district (LAD)/unitary authority (UA)/ metropolitan district (MD)/ London borough (LB)/ council area (CA)/district council area (DCA)
#  "OSWARD", #	(Electoral) ward/division
  "OSGRDIND", #	Grid reference positional quality indicator
  "CTRY", #	Country
#  "OSHLTHAU", #	Former Strategic health authority (SHA)/ health board (HB)/ health authority (HA)/ health & social care board (HSCB)
  "NHSCR", #	Commissioning Region
  "CCG", #	Clinical Commissioning Group (CCG)/ local health board (LHB)/ community health partnership (CHP)/ local commission-ing group (LCG)/ primary healthcare directorate (PHD)
#  "OA01", #	2001 Census Output Area
  "NHSRG", #	NHS Region aka Area Team
#  "HRO", #	Former Pan SHA
  "LSOA01", #	2001 Lower Layer Super Output Area (LSOA)
  "SCN", #	Strategic Clinical Network
#  "CANREG", #	Cancer Registry
#  "PCT", #	‘Second wave’ Primary Care Trust (PCT)/ Care Trust/ Care Trust Plus (CT)/ local health board (LHB)/ community health partnership (CHP)/ local commission-ing group (LCG)/ primary healthcare directorate (PHD)
  "OSEAST1M", #	National grid reference - Easting
  "OSNRTH1M", #	National grid reference - Northing
#  "OA11", #	2011 Census output area
  "LSOA11" #	2011Census  lower layer super output area (LSOA)
#  "MSOA11" #	2011 Census middle layer super output area (MSOA)
)

# CodePoint stuff
.codelist.xl.filename <- "Codelist.xlsx" # 2015.4.0 version uses .xlsx
.codelist.xl.areacodes.sheetname <- "AREA_CODES" # name of sheet within Codelist.xls that has table of the different types of codes
.nhscodelist.xl.filename <- "NHS_Codelist.xls"

# .read_codelist ----------------------------------------------------------

.read_codelist <- function(sheet,path,NHS=FALSE) {
  z <- read_excel(path,sheet,col_names=FALSE)
  if(NHS){
    z[,3] <- "NHS"
  } else {
    z[,3] <- sheet # 3 letter abbreviation for the type of code, e.g. CTW=County
  }
  return(z)
}

# getCodePointAreaCodes ---------------------------------------------------

#' get CodePoint area codes
#'
#' obtain the names of all the area codes from an unzipped OS CodePoint Open zipped file
#' and compile them into format that I want
#'
#' @param codepoint.doc.path full path to where we can find the unzipped xls/xlsx files that contain area code names
#'
#' @return data.frame, compiled list of all the area codes
#' @export
#'
#' @examples \dontshow{getCodePointAreaCodes()}
getCodePointAreaCodes <- function( codepoint.doc.path = paste(.default_CodePoint_Unzipped_Path, "Doc", sep = "/") ) {

  codelist.xl.file    <- paste(codepoint.doc.path, .codelist.xl.filename, sep = "/")
  nhscodelist.xl.file <- paste(codepoint.doc.path, .nhscodelist.xl.filename, sep = "/")

  # TODO error checking of path syntax

  # compile all the area codes
  area_codes <- as.list(read_excel(codelist.xl.file, sheet = .codelist.xl.areacodes.sheetname, col_names = FALSE)[,1])
  area_codes <- lapply( area_codes, .read_codelist,  path = codelist.xl.file )
  area_codes <- do.call(rbind, area_codes)
  names(area_codes) <- c("AreaName", "AreaID", "AreaType")

  nhs_areas <- excel_sheets(nhscodelist.xl.file)
  nhs_areas <- lapply( nhs_areas, .read_codelist, path = nhscodelist.xl.file, NHS = TRUE )
  nhs_areas <- do.call(rbind, nhs_areas)
  names(nhs_areas) <- c("AreaID", "AreaName", "AreaType")

  CodePointAreaCodes <- rbind(area_codes,nhs_areas)
  CodePointAreaCodes  <- CodePointAreaCodes[complete.cases(CodePointAreaCodes),]

  return(CodePointAreaCodes)
}

# getNHSPD ----------------------------------------------------------------

#' get NHS Postcode Directory
#'
#' if RDS file exists, load it
#' else if data file (zip/csv) exists, process it (and also save as RDS for next time)
#' else download, process and save
#'
#' @param gridall.file full path to CSV or ZIP containing gridall.csv
#' @param gridall.URL full URL to download
#' @param gridall.RDS full path to saved (processed) data
#' @param force.download Force download of file from source
#' @param remake.RDS Rebuild and save the data
#' @param codepoint.doc.path if \code{remake.RDS==TRUE} - specify where to find unzipped OS CodePoint documentation - see getCodePointAreaCodes()
#' @param codepoint.compatible if \code{remake.RDS==TRUE} - column names to be compatible with my previous codepoint library
#' @param only.England if \code{remake.RDS==TRUE} - only return rows applicable to England
#'
#' @return data.table
#' @export
#'
#' @examples NHSPD <- getNHSPD()
getNHSPD <- function( gridall.file = .default_gridall_DataFile,
                      gridall.URL = .default_gridall_URL,
                      gridall.RDS = .default_gridall_RDS,
                      force.download = FALSE,
                      remake.RDS = FALSE,
                      codepoint.compatible = TRUE,
                      codepoint.doc.path = paste(.default_CodePoint_Unzipped_Path, "Doc", sep = "/"),
                      only.England = TRUE ) {

  warning(.citations())

  if (force.download == FALSE & remake.RDS == FALSE) {
    if(file.exists(gridall.RDS)) {
      message(sprintf("getNHSPD() loading RDS data file '%s'...", gridall.RDS))
      return(readRDS(gridall.RDS))
    } else {
      stop(sprintf("getNHSPD() unable to load RDS file '%s'.  To try to make it, specificy remake.RDS = TRUE next time.", gridall.RDS))
      # alternatively, don't stop - just set remake.RDS <- TRUE
    }
  }

  # remainder combinations are:
  #   force.download = TRUE, remake.RDS == FALSE (illogical, therefore overide)
  #   force.download = TRUE, remake.RDS == TRUE
  #   force.download = FALSE, remake.RDS == TRUE

  if (force.download) { remake.RDS <- TRUE } # overide the illogical combination

  # by this stage, remake.RDS has to be TRUE

  if (file.exists(gridall.RDS)) file.remove(gridall.RDS) # or rename as .bak / .previous

  gridall.file.zip <- paste0( file_path_sans_ext(gridall.file), ".zip")
  DataFile.missing <- !file.exists(gridall.file)

  if(force.download | DataFile.missing) {
    if (force.download) {
      excuse <- "force.download = TRUE"
    } else {
      excuse <- "missing local copy"
    }
    message(sprintf("getNHSPD() downloading gridall.zip from '%s' to '%s' because: '%s'...", gridall.URL, gridall.file.zip, excuse))
    download.success <- download.file( gridall.URL, destfile = gridall.file.zip, method = "auto") # should we warn if overwrite?
    if (download.success == 0) {
      if (file.exists(gridall.file)) file.remove(gridall.file)
      message(sprintf("getNHSPD() unzipping..."))
      unzip(gridall.file.zip, exdir = dirname(gridall.file.zip) ) # should we warn if overwrite?
    } else {
      if (DataFile.missing == FALSE) {
        # download failure (force.download==TRUE), but option to use existing therefore warn, or stop with error
        message("getNHSPD() Unable to download from source, re-using existing file.")
      } else {
        stop("getNHSPD() Unable to download from source; no existing file to use.")
      }
    }
  }

  # process file (make RDS file)
  message(sprintf("getNHSPD() loading data from '%s'...", gridall.file))
  z <- fread( gridall.file, col.names = .headers.wanted, select = match(.headers.wanted, .headers.short) )

  # clean up the data
  message(sprintf("getNHSPD() cleaning up data..."))

  if (only.England) {
    z <- subset(z, CTRY == "E92000001")
  }

  if (FALSE) {
    # remove non-geographic postcodes
    # TODO ?
  }

  if (TRUE) {
    # remove terminated postcodes
    z <- subset(z, DOTERM == "")
    z[, DOTERM:=NULL] # drop this column
  }

  if (codepoint.compatible) {
    setnames(z, c("OSEAST1M", "OSNRTH1M", "PCD2"), c("Eastings", "Northings", "Postcode"))
    z[, Eastings:=as.integer(Eastings)]
    z[, Northings:=as.integer(Northings)]

    #z$Postcode <- paste0(str_sub(z$PCD8,1,4),str_sub(z$PCD8,-3,-1)) # make it the same format as CodePoint database
    z[, Postcode:=paste0(str_sub(Postcode, 1, 4), str_sub(Postcode, -3, -1))]
    # z$Admin_county_code <- z$OSCTY # ? to factorise // probably can remove
    z[, Admin_district_code:=OSLAUA] # ? to factorise

    # add codelist (names) from OS CodePoint
    message("getNHSPD() obtaining area names from local copy of OS CodePoint...")
    CodePointAreaCodes <- as.data.table(getCodePointAreaCodes(codepoint.doc.path))

    # join
    message("getNHSPD() joining areanames with NHSPD...")
    z <- merge(z, CodePointAreaCodes, by.x = "OSLAUA", by.y = "AreaID", all.x = TRUE)
  }

  # save RDS
  message(sprintf("getNHSPD() saving data to '%s'...", gridall.RDS))
  saveRDS(z, file = gridall.RDS)

  # return
  return(z)
}