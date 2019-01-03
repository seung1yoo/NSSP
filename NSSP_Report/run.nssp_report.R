if(!require(stringr)){
  install.packages("stringr")
}
if(!require(lubridate)){
  install.packages("lubridate")
}
library(stringr)
library(lubridate)

########################################
## Setting Area
# Project informations

project_id <- "TBD180101"
reporting_date <- as.Date("2018-11-01")
data_path <- paste("/Users/Yoo/YoosScripts/NSSP", "NSSP_Data", sep = '/')

# Analysis informations
species_name <- "Fusarium graminearum"
ref_source <- "BioMax"

# report informations
report_path <- "/Users/Yoo/YoosScripts/NSSP/NSSP_Report"


########################################

xmonth <- month(as.Date(reporting_date))
xday <- day(as.Date(reporting_date))
dmonth <- stringr::str_pad(xmonth,2,pad="0")
dday <- stringr::str_pad(xday,2,pad="0")

try(
  rmarkdown::render(paste(report_path, "report_template/nssp_report_typeA.Rmd", sep = '/'),
                    params = list(
                      dmonth = dmonth,
                      dday = dday,
                      project_id = project_id,
                      data_path = data_path,
                      species_name = species_name,
                      ref_source = ref_source,
                      report_path = report_path
                      ),
                    output_format = "pdf_document",
                    output_file = paste(data_path, paste0(project_id, "_nssp_analysis_report", ".pdf"), sep = '/'),
                    encoding = "UTF-8"),
  silent = FALSE
)
