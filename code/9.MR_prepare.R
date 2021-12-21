home_PATH<-"/data/sharedcode/kjkim/data/MDHC_DM"
setwd(home_PATH)
TB <- read_excel("meta.xlsx", # path
#sheet = "cust_profile", # sheet name to read from
#range = "B3:E8", # cell range to read from
col_names = TRUE, # TRUE to use the first row as column names
col_types = "guess", # guess the types of columns
skip=3,
na = "NA")


names(TB)<-c("No","ID","Group","Type","NAME","Nationality","Age","Sex","Date","RegID")
table(TB$Group)