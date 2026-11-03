# =============================================================================
# AOP_selection.R
#
# Purpose : Build the node (key event) and edge (key event relationship) tables
#           for an organ-specific Adverse Outcome Pathway (AOP) network by
#           querying the AOP-Wiki RDF SPARQL endpoint. This copy builds the
#           liver network; the brain, lung and kidney networks are built with
#           the same script by changing the keyword list in the Part 1 query.
#
# Workflow: Part 1 - keyword search of AOP-Wiki for organ-relevant AOPs
#           Part 2 - manual exclusion of AOPs whose adverse outcome is not
#                    organ-relevant
#           Part 3 - retrieval of all KEs of the retained AOPs  (network nodes)
#           Part 4 - retrieval of all KERs of the retained AOPs (network edges)
#
# Inputs  : AOP-Wiki RDF SPARQL endpoint
#             https://aopwiki.rdf.bigcat-bioinformatics.org/sparql
#           selected_liver_aops.tsv  (written by Part 2, read back in Part 3)
#
# Outputs : all_liver_aops.tsv       all AOPs returned by the keyword search
#           selected_liver_aops.tsv  AOPs retained after manual exclusion
#           kes.tsv                  network nodes (key events)
#           ker.tsv                  network edges (key event relationships)
#           AOP_selection.RData      full session image
#
# Note    : Parts 1 and 2 are commented out because they only need to be run
#           once. Their output (selected_liver_aops.tsv) is read back in at the
#           start of Part 3. Uncomment them to regenerate the AOP selection.
#
# =============================================================================

# setting working directory
setwd("")
getwd()
# NOTE: change the path above to the location of this repository on your machine.



# Optional: restore a previously saved session instead of re-running the queries.

# load('AOP_selection.RData')
# call libraries
# httr      sends the POST requests to the SPARQL endpoint
# tidyverse provides read_tsv() and general data handling
# dplyr     data frame manipulation
# jsonlite  called below with :: to parse the JSON responses

library(httr)
library(tidyverse)
library(dplyr)

# All queries in this script are sent to the AOP-Wiki RDF SPARQL endpoint.
# Results depend on the AOP-Wiki release loaded at the endpoint at the time of
# querying, so record the release date when reporting results.

# SPARQL endpoint (API)
# https://aopwiki.rdf.bigcat-bioinformatics.org/sparql


# =============================================================================
# PART 1 - Keyword search: find candidate organ-specific AOPs
# Run once, then left commented out. Output: all_liver_aops.tsv
# =============================================================================

###### get all AOPs associated with Liver based on defined key words ###### 

# # define query
# The query returns, for every AOP whose title OR adverse outcome title matches
# one of the liver keywords: the AOP URI, AOP number, AOP title, AO title, the
# list of its KE identifiers, and its KE count. An AOP with two adverse
# outcomes is returned twice, once per AO.

# query <- '
# SELECT DISTINCT ?aop AS ?link 
#        (CONCAT("AOP ", ?aopNumber) AS ?aopNumber) 
#        ?title 
#        ?aotitle 
#        (GROUP_CONCAT(DISTINCT ?keid; SEPARATOR=", ") AS ?KEids) 
#        (COUNT(DISTINCT ?keid) AS ?totalKEs)
# WHERE {
#     ?aop a aopo:AdverseOutcomePathway ;
#          dc:title ?title ;
#          aopo:has_adverse_outcome ?ao ;
#          aopo:has_key_event ?ke .
#     
#     ?ao dc:title ?aotitle .
#     ?ke rdfs:label ?keid .
# 
#     # Merge title and AO title
#     BIND(CONCAT(LCASE(?title), " ", LCASE(?aotitle)) AS ?searchText)
# 
#     # AOP number
#     BIND(STRAFTER(STR(?aop), "aop/") AS ?aopNumber)
# 
#     # Use regex for multiple keywords at once
# The keyword list below defines organ relevance and is the only part of the
# query that changes between organs (for example "brain|neuro|cognitive|..."
# for the brain network). The "i" flag makes the match case-insensitive.

#     FILTER regex(?searchText, "liver|hepatic|hepato|steatosis|hepatitis|cholestasis|cirrhosis|kupffer cells|hepatic sinusoidal|bile canaliculi|fatty liver|lipidosis", "i")
# }
# GROUP BY ?aop ?aopNumber ?title ?aotitle
# ORDER BY ?aopNumber
# 
# '
# 
# # Get response
# Send the query as a form-encoded POST request and ask for JSON output.

# response <- httr::POST("https://aopwiki.rdf.bigcat-bioinformatics.org/sparql", 
#                        body = list(query = query, format = "json"), 
#                        encode = "form")
# # read the results
# data <- jsonlite::fromJSON(httr::content(response, "text"))
# results <- data$results$bindings
# 
# # view results
# head(results,1)
# str(results)
# names(results)
# 
# # fetching only value
# results$title$value
# 
# The SPARQL JSON response nests every field as a list with a $value element.
# This step flattens it into a plain data frame with one column per variable.

# # Extract only the value columns from nested dataframes
# results <- data.frame(
#   link = results$link$value,
#   aopNumber = results$aopNumber$value,
#   title = results$title$value,
#   aotitle = results$aotitle$value,
#   KEids = results$KEids$value,
#   totalKEs = results$totalKEs$value
#   
# )
# 
# all_aops <- results
# 
# # save liver aops
# Save the unfiltered keyword-search result, so the manual exclusion below can
# be traced back to a fixed starting point.

# write.table(all_aops, file = "all_liver_aops.tsv", 
#             sep = "\t", 
#             row.names = FALSE, 
#             quote = FALSE)
# 
# ######################################
# 
# =============================================================================
# PART 2 - Manual exclusion of AOPs that are not liver-relevant
#
# The keyword search also returns AOPs that mention the liver but whose adverse
# outcome lies in another organ, or that describe non-mammalian species. These
# are inspected by eye and removed.
#
# WARNING: the row indices below refer to the row order of `results` for the
# AOP-Wiki release used here. Re-running Part 1 against a newer release will
# change the row order, so the AOPs to exclude must be re-identified from their
# titles rather than reused as indices.
# =============================================================================

# ###### manually select the AOPs ###### 
# 
# Inspect the title and adverse outcome title of each candidate for exclusion.
# Keep a record of which AOPs were excluded, for transparency and reporting.

# excluded_aops <- results[c(42, 35, 29, 9),]
# 
# Retain every AOP except the four excluded above.
# selected_liver_aops <- all_aops[-c(42, 35, 29, 9),]
# 
# # save results
# write.table(selected_liver_aops, file = 'selected_liver_aops.tsv',
#             sep = '\t',
#             row.names = F,
#             quote = F)

######################################

######################################

# =============================================================================
# PART 3 - Network nodes: retrieve all key events of the selected AOPs
# =============================================================================

# read input file
selected_liver_aops <- read_tsv('selected_liver_aops.tsv')

# This sums the per-AOP KE counts, so KEs shared between AOPs are counted once
# per AOP. It is therefore larger than the number of unique nodes (see below).

# see the total num of KEs from selected_liver_aops
sum(as.integer(selected_liver_aops$totalKEs)) # 302

# there are 2 AOPs [273, 278] are appeared twice bc these have 2 AOs

# Build the SPARQL VALUES clause from the retained AOP numbers. The printed
# string is pasted by hand into the ?aopNumber block of the two queries below.

# Get the clean quoted list of AOP ids for query building
aop_list <- paste0('"', selected_liver_aops$aopNumber, '"', collapse = "")
cat(aop_list)

# this list will be pasted in the query

#####################################################

# =============================================================================
# The query returns every KE of the selected AOPs, classified by its role:
#   MIE  molecular initiating event
#   AO   adverse outcome
#   KE   intermediate key event (neither MIE nor AO within that AOP)
# The three UNION branches are mutually exclusive: the FILTER NOT EXISTS lines
# in the third branch stop an MIE or AO from also being returned as a plain KE.
# A KE can still appear more than once, when it belongs to several AOPs or has
# a different role in different AOPs.
# =============================================================================

##### Get all the KEs (nodes) for liver network #####

query_kes <- '

SELECT DISTINCT ?aopNumber ?ke ?keLabel ?keTitle ?type
WHERE {
  VALUES ?aopNumber {"AOP 1""AOP 107""AOP 108""AOP 117""AOP 118""AOP 130""AOP 144""AOP 162""AOP 209""AOP 213""AOP 220""AOP 232""AOP 27""AOP 273""AOP 273""AOP 278""AOP 278""AOP 285""AOP 318""AOP 34""AOP 36""AOP 362""AOP 37""AOP 38""AOP 383""AOP 41""AOP 423""AOP 46""AOP 494""AOP 517""AOP 518""AOP 529""AOP 57""AOP 58""AOP 591""AOP 60""AOP 61""AOP 62"}
 
   ?aop a aopo:AdverseOutcomePathway ; rdfs:label ?aopNumber .
  {
         ?aop aopo:has_molecular_initiating_event ?ke .
        ?ke rdfs:label ?keLabel ;
            dc:title ?keTitle .
        BIND("MIE" AS ?type)
    } UNION {
        ?aop aopo:has_adverse_outcome ?ke .
        ?ke rdfs:label ?keLabel ;
            dc:title ?keTitle .
        BIND("AO" AS ?type)
    } UNION {
        ?aop aopo:has_key_event ?ke .
        ?ke rdfs:label ?keLabel ;
            dc:title ?keTitle .
        
        FILTER NOT EXISTS { ?aop aopo:has_molecular_initiating_event ?ke }
        FILTER NOT EXISTS { ?aop aopo:has_adverse_outcome ?ke }
        BIND("KE" AS ?type)
    }
}
ORDER BY ?aopNumber ?type ?keLabel

'

# Get response
response_kes <- httr::POST("https://aopwiki.rdf.bigcat-bioinformatics.org/sparql", 
                           body = list(query = query_kes, format = "json"), 
                           encode = "form")
# read the results
data_kes <- jsonlite::fromJSON(httr::content(response_kes, "text"))
results_kes <- data_kes$results$bindings

# view results
head(results_kes,1)
str(results_kes)
names(results_kes)

# Flatten the nested JSON response into a data frame, as in Part 1.

# Extract only the value columns from nested dataframes
results_kes <- data.frame(
  aopNumber = results_kes$aopNumber$value,
  ke = results_kes$ke$value,
  keLabel = results_kes$keLabel$value,
  keTitle = results_kes$keTitle$value,
  type = results_kes$type$value
)


# Write the node table for the network.

write.table(results_kes, file = "kes.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

# Row count includes KEs that occur in more than one AOP.

# see the total number of KEs
# str(results_kes) # 278 rows
length(results_kes$keLabel) # 278

# Unique KE labels give the actual number of nodes in the liver network.

# see the number of UNIQUE KEs/ nodes
length(unique(results_kes$keLabel)) # 196.Note: 
# Total number of duplicated entries
sum(duplicated(results_kes$keLabel)) # 82

# Total 278 = 196 + 82

# see how many times each KE appears
table(results_kes$keLabel)

# see only those with more than 1 occurrence
table(results_kes$keLabel)[table(results_kes$keLabel) > 1]

# the reason for repetition: AOPs 273 and 278 are repeated twice as these have 2 AOs each. Moreover, some KEs are present in more than 1 AOPs

######################################

######################################

# =============================================================================
# PART 4 - Network edges: retrieve the key event relationships (KERs)
# =============================================================================

##### Get the KERs (edges) for liver network #####
# updated the query to get the the KERs based on the provided AOP list rather than from KEs list.

# The query returns each KER of the selected AOPs with its upstream and
# downstream KE. Querying by AOP rather than by KE, as noted above, prevents
# KERs that link a selected KE to a KE outside the liver AOP set from being
# pulled into the network.

# define query
query_ker <- '

# for edges
SELECT DISTINCT  ?ker ?kerLabel ?upstreamKE ?downstreamKE ?upstreamLabel ?downstreamLabel 
WHERE {
    VALUES ?aopNumber {"AOP 1""AOP 107""AOP 108""AOP 117""AOP 118""AOP 130""AOP 144""AOP 162""AOP 209""AOP 213""AOP 220""AOP 232""AOP 27""AOP 273""AOP 273""AOP 278""AOP 278""AOP 285""AOP 318""AOP 34""AOP 36""AOP 362""AOP 37""AOP 38""AOP 383""AOP 41""AOP 423""AOP 46""AOP 494""AOP 517""AOP 518""AOP 529""AOP 57""AOP 58""AOP 591""AOP 60""AOP 61""AOP 62"}
  
  ?aop a aopo:AdverseOutcomePathway ; rdfs:label ?aopNumber ; aopo:has_key_event_relationship ?ker .
           ?ker a aopo:KeyEventRelationship ;
         aopo:has_upstream_key_event ?upstreamKE ;
         aopo:has_downstream_key_event ?downstreamKE ;
           rdfs:label ?kerLabel .
    
    ?upstreamKE rdfs:label ?upstreamLabel .
    ?downstreamKE rdfs:label ?downstreamLabel .
      
}

'

# Get response
response_ker <- httr::POST("https://aopwiki.rdf.bigcat-bioinformatics.org/sparql", 
                           body = list(query = query_ker, format = "json"), 
                           encode = "form")
# read the results
data_ker <- jsonlite::fromJSON(httr::content(response_ker, "text"))
results_ker <- data_ker$results$bindings

# view results
head(results_ker,1)
str(results_ker)
names(results_ker)

# Flatten the nested JSON response into a data frame, as in Parts 1 and 3.

# Extract only the value columns from nested dataframes
results_ker <- data.frame(
  ker = results_ker$ker$value,
  kerLabel = results_ker$kerLabel$value,
  upstreamKE = results_ker$upstreamKE$value,
  downstreamKE = results_ker$downstreamKE$value,
  upstreamLabel = results_ker$upstreamLabel$value,
  downstreamLabel = results_ker$downstreamLabel$value
)

# Write the edge table. Together, kes.tsv and ker.tsv define the liver AOP
# network that is imported into Cytoscape in the next script.

# save results
write.table(results_ker, file = "ker.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

# save session 
save.image('AOP_selection.RData')
load('AOP_selection.RData')
