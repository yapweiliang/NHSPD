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
#' This library therefore supercedes my \pkg{codepoint} library but it still references
#' Area Code names, which can be obtained from the \bold{ONS} or \bold{OS CodePoint Open}
#'
#' @section Notes:
#' on comparing NHSPD data (November 2015) with OS CodePoint Open 2015.4.0 (October 2015),
#' at least the England postcodes with grid coordinates are identical to OS CodePoint
#' except for 4 postcodes where NHSPD data seems correct if plotted on map
#' (TR210AA, TR210PR, TR210PW, TR210WY labelled as Isles of Scilly on CodePoint vs Cornwall on NHSPD)
#'
#' @section TODO:
#' TODO: error checking for CodePoint stuff
#'
#' @family functions
#'
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

.default_LAD_Datafile <- "H:/DATASETS/OS/Area Codes/LAD_2015_UK_NC.csv" # local copy of extracted download
.default_LAD_URL <- "https://geoportal.statistics.gov.uk/Docs/Names%20and%20Codes/Local_Authority_Districts_(UK)_2015_names_and_codes.zip"

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

# .change_year ------------------------------------------------------------

.change_year <- function(s, yr = as.numeric(format(Sys.time(), "%Y")) - 1) {
  # change the year in the URL or data file name to current year - 1, or any specified year
  return(gsub("_(\\d{4})_",sprintf("_%04d_",as.numeric(yr)),s))
}

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
#' obtain the names of all the area codes from an unzipped \emph{OS CodePoint Open} zipped file
#' and compile them into format that I want.
#'
#' To download:
#' \itemize{
#'  \item request from \url{https://www.ordnancesurvey.co.uk/business-and-government/products/code-point-open.html}
#'  \item download OS CodePoint Open into a temporary folder e.g. \code{"/your/folder/tmp"}
#'  \item set \code{.default_CodePoint_Unzipped_Path <- "/your/folder/tmp"} in this library, or,
#'  \item pass \code{codepoint.doc.path = "/your/folder/tmp/Doc"}
#' }
#'
#' @family functions
#'
#' @param codepoint.doc.path full path to where we can find the unzipped \code{Codelist.xls/xlsx} files that contain area code names
#'
#' @return data.table, compiled list of all the area codes
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

  return(as.data.table(CodePointAreaCodes))
}

# getAreaCodes ------------------------------------------------------------

#' get Area Codes
#'
#' loads the list of Area Codes from local copy or attempts to download and unzip it
#'
#' Download and extract \code{LAD_\emph{yyyy}_UK_NC.csv} from
#' \url{https://geoportal.statistics.gov.uk/Docs/Names\%20and\%20Codes/Local_Authority_Districts_(UK)_2015_names_and_codes.zip}
#' or allow the library to do that for you.
#'
#' @family functions
#'
#' @param LAD.file full path to the unzipped file
#' @param LAD.URL full URL to download
#' @param year \code{NULL} to guess the latest year, \code{YYYY} to specify
#' @param force.download Force download of file from source
#'
#' @return data.table, columns renamed to AreaID and AreaName
#' @export
#'
#' @examples \dontshow{AreaCodes <- getAreaCodes()}
getAreaCodes <- function( LAD.file = .default_LAD_Datafile,
                          LAD.URL = .default_LAD_URL,
                          year = NULL,
                          force.download = FALSE ) {

  if (is.null(year)) {
    # guess the latest year
    LAD.file <- .change_year(LAD.file)
    LAD.URL <- .change_year(LAD.URL)
  } else {
    # specify the year
    LAD.file <- .change_year(LAD.file, year)
    LAD.URL <- .change_year(LAD.URL, year)

    if (year <= 2014) {
      # change ...names_and_codes... to ...Names_and_Codes...
      # change ...L_A_D... to ...L_a_d...
      LAD.URL <- gsub("names_and_codes", "Names_and_Codes", LAD.URL)
      LAD.URL <- gsub("Local_Authority_Districts", "Local_authority_districts", LAD.URL)
    }
  }

  LAD.file.zip <- paste0(dirname(LAD.file), "/", basename(LAD.URL))
  DataFile.missing <- !file.exists(LAD.file)

  if(force.download | DataFile.missing) {
    if (force.download) {
      excuse <- "force.download = TRUE"
    } else {
      excuse <- "missing local copy"
    }
    message(sprintf("getAreaCodes() downloading LAD code names from '%s' to '%s' because: '%s'...", LAD.URL, LAD.file.zip, excuse))
    download.success <- try(download.file( LAD.URL, destfile = LAD.file.zip, method = "auto")) # should we warn if overwrite?
    if (download.success == 0) {
      if (file.exists(LAD.file)) file.remove(LAD.file)
      message(sprintf("getAreaCodes() unzipping..."))
      unzip(LAD.file.zip, files = basename(LAD.file), exdir = dirname(LAD.file.zip) ) # should we warn if overwrite?
    } else {
      if (DataFile.missing == FALSE) {
        # download failure (force.download==TRUE), but option to use existing therefore warn, or stop with error
        message("getAreaCodes() Unable to download from source, re-using existing file.")
      } else {
        stop("getAreaCodes() Unable to download from source; no existing file to use.")
      }
    }
  }

  z <- fread(LAD.file) # column names are LADnnCD, LADnnNM where nn is the year
  setnames(z, c("AreaID", "AreaName")) # these are the names I used in my previous codepoint library

  if(FALSE) {
    # optionally add a 3rd column with NULL values named AreaType to be consistent with my previous codepoint library
  }

  return(z)
}


# getNHSPD ----------------------------------------------------------------

#' get NHS Postcode Directory
#'
#' downloads/compiles/re-loads the NHS Postcode Directory in a form that I want
#'
#' if RDS file exists, load it\cr
#' else if data file (zip/csv) exists, process it, and save as RDS for next time\cr
#' else download, process and save for next time
#'
#' Download from:
#' \itemize{
#'  \item \url{http://systems.hscic.gov.uk/data/ods/datadownloads/onsdata}
#'  \item \url{http://systems.hscic.gov.uk/data/ods/datadownloads/onsdata/zip-files/gridall.zip}
#' }
#' or allow the library to do that for you.
#'
#' The \bold{OS CodePoint Open} area names are fuller, e.g. "Sandwell District (B)", instead of "Sandwell" from \bold{ONS}
#'
#' @family functions
#'
#' @param gridall.file full path to CSV or ZIP containing gridall.csv
#' @param gridall.URL full URL to download
#' @param gridall.RDS full path to saved (processed) data
#' @param force.download Force download of file from source
#' @param remake.RDS Rebuild and save the data
#' @param AreaCodes supply a data.table with area codes.\cr
#' defaults to \code{\link{getAreaCodes}()} for brief names from ONS \cr
#' but could also use \code{\link{getCodePointAreaCodes}()} which has fuller names.
#' @param codepoint.compatible if \code{remake.RDS==TRUE} - column names to be compatible with my previous \pkg{codepoint} library
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
                      AreaCodes = getAreaCodes( force.download = force.download ),
                      codepoint.compatible = TRUE,
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
    download.success <- try(download.file( gridall.URL, destfile = gridall.file.zip, method = "auto")) # should we warn if overwrite?
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

    # join
    message("getNHSPD() joining areanames with NHSPD...")
    z <- merge(z, AreaCodes, by.x = "OSLAUA", by.y = "AreaID", all.x = TRUE)
  }

  # save RDS
  message(sprintf("getNHSPD() saving data to '%s'...", gridall.RDS))
  saveRDS(z, file = gridall.RDS)

  # return
  return(z)
}