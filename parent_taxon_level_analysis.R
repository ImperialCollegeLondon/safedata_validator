library(RSQLite)

# connect to the database
mydb <- dbConnect(RSQLite::SQLite(), "backbone-current.sqlite")

# set the taxonomic hierarchy order
tax_hier <- c("kingdom","phylum","class","order","family","genus",
			 "species","subspecies", "variety", "form")

# accepted taxa: are all children less nested than parents in rank
accepted_data <- dbGetQuery(mydb, "
SELECT a.taxonRank AS child, 
       b.taxonRank AS parent, 
       count(*) AS n
    FROM (
        SELECT parentNameUsageID, taxonRank 
            FROM backbone 
            WHERE taxonomicStatus == 'accepted') a 
    JOIN (
        SELECT taxonID, taxonRank 
            FROM backbone 
            WHERE taxonomicStatus == 'accepted') b
    ON a.parentNameUsageID = b.taxonID
    GROUP BY a.taxonRank, b.taxonRank;")

accepted_data$parent <- factor(accepted_data$parent, levels = tax_hier)
accepted_data$child <- factor(accepted_data$child, levels = tax_hier)

accepted_xtab <- xtabs(n ~ child + parent, data=accepted_data)

# answer - yes but not always the next rank up:
diag_index <- row(accepted_xtab) - col(accepted_xtab)
sum(accepted_xtab[diag_index > 0])/sum(accepted_xtab)
sum(accepted_xtab[diag_index > 1])/sum(accepted_xtab)

# unaccepted taxa: how bad is it?
unaccepted_data <- dbGetQuery(mydb, "
SELECT a.taxonRank AS child, 
       b.taxonRank AS parent, 
       count(*) AS n
    FROM (
        SELECT parentNameUsageID, taxonRank 
            FROM backbone 
            WHERE taxonomicStatus != 'accepted') a 
    JOIN (
        SELECT taxonID, taxonRank 
            FROM backbone 
            WHERE taxonomicStatus == 'accepted') b
    ON a.parentNameUsageID = b.taxonID
    GROUP BY a.taxonRank, b.taxonRank;")

unaccepted_data$parent <- factor(unaccepted_data$parent, levels = tax_hier)
unaccepted_data$child <- factor(unaccepted_data$child, levels = tax_hier)

unaccepted_xtab <- xtabs(n ~ child + parent, data=unaccepted_data)

# a lot map to equal or less nested status
sum(unaccepted_xtab[diag_index > 0])/sum(unaccepted_xtab)
sum(unaccepted_xtab[diag_index > 1])/sum(unaccepted_xtab)

# specifically doubtful taxa: how bad is it?
doubtful_data <- dbGetQuery(mydb, "
SELECT a.taxonRank AS child, 
       b.taxonRank AS parent, 
       count(*) AS n
    FROM (
        SELECT parentNameUsageID, taxonRank 
            FROM backbone 
            WHERE taxonomicStatus == 'doubtful') a 
    JOIN (
        SELECT taxonID, taxonRank 
            FROM backbone 
            WHERE taxonomicStatus == 'accepted') b
    ON a.parentNameUsageID = b.taxonID
    GROUP BY a.taxonRank, b.taxonRank;")

doubtful_data$parent <- factor(doubtful_data$parent, levels = tax_hier)
doubtful_data$child <- factor(doubtful_data$child, levels = tax_hier)

doubtful_xtab <- xtabs(n ~ child + parent, data=doubtful_data)

# a lot map to equal or less nested status
diag_index <- row(doubtful_xtab) - col(doubtful_xtab)
sum(doubtful_xtab[diag_index > 0])/sum(doubtful_xtab)
sum(doubtful_xtab[diag_index > 1])/sum(doubtful_xtab)



library(knitr)
kable(accepted_xtab)
kable(unaccepted_xtab)
kable(doubtful_xtab)



