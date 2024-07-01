# (1) Groupings of variables --------------------------------------------------
var_groups <- list()

#demographics
var_groups[["demographics"]] <- list()
  
var_groups$demographics[["population"]] <- c("total_population", "cvap")
var_groups$demographics[["race_composit"]] <- c("white_pct", "black_pct", "hispanic_pct", "otherrace_pct")
var_groups$demographics[["race_other"]] <- c("nonwhite_pct", "foreignborn_pct")
var_groups$demographics[["gender"]] <- c("female_pct", "male_pct")
var_groups$demographics[["age"]] <- c("age29andunder_pct", "age30to64", "age65andolder_pct")
var_groups$demographics[["economic_situation"]] <- c("median_hh_inc", "clf_unemploy_pct", "clf_employed_pct")
var_groups$demographics[["education"]] <- c("lesshs_pct", "onlyhs_pct", "morecollege_pct")
var_groups$demographics[["education_whites"]] <- c("lesshs_whites_pct", "onlyhs_whites_pct", "morecollege_whites_pct")
var_groups$demographics[["rurality"]] <- c("rural_pct", "ruralurban_cc")

#elections results
var_groups[["elections"]] <- list()

var_groups$elections[["2020"]] <- c("gop_pct20", "dem_pct20", "others_pct20")
var_groups$elections[["2016"]] <- c("gop_pct16", "dem_pct16", "others_pct16")
var_groups$elections[["2012"]] <- c("gop_pct12", "dem_pct12", "others_pct12")

# (2) Next ----------------------------------------------------------------



