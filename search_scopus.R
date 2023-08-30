# Scopus search
library(data.table)
library(rscopus)
library(dplyr)
library(dotenv)
setwd('~/Desktop/Review-PlateletsRNABPsMotifs/')
load_dot_env('.env')

## START UGR VPN CONNECTION

# Add API KEY
# ===

APIKEY = Sys.getenv('APIKEY') 
rscopus::set_api_key(APIKEY)
queryID = 'plateletsRNAmotifs'

# Build the query
# keywords 
## 1- RNA motif, RBP, platelets, megakaryocytes, RNA localization

# ===
datatype = c('RNABP','rna binding protein','binding motifs', 
            'clip-seq','RNA localization','RNA motif')
celltype = c('platelet','megakaryocyte')
organism = c('human','homo sapiens')

sourceSample = list(organism = organism,
                    tissuestate = tissuestate)

#datatype = c('transcriptome','transcriptomic','proteome','proteomic',
#             'RNAseq','RNA-Seq', 'RNABP','binding proteins','binding motifs',
#             'gene expression', 'chip-seq','clip-seq','RNA localization','RNA motif')

queries = list()
query = list()
for (d in seq_along(datatype)) {
  for (i in seq_along(sourceSample)) {
    y = paste0('TITLE-ABS-KEY(', datatype[d], ') & ', 'TITLE-ABS-KEY(', celltype, ') & ')
    x = paste0('TITLE-ABS-KEY(', sourceSample[[i]], ')')
    x = paste(x, collapse = " OR ")
    queries[[i]] = paste0(y, x)
  }
  query[[d]] = unlist(queries, use.names = F)
}
query = unlist(query)

# Omic queries
# ===
# Gene expression
queryRes = list()
completeArticle = list()
for (qidx in seq_along(query)) {
  # Run query 
  cat(paste0('Search (',qidx,') for:\n', query[qidx]))
  completeArticle[[qidx]] <- scopus_search(
    query = query[qidx], 
    view = "COMPLETE", 
    count = 25)
  
  # Identify datatypestr from query
  datatypestr = strsplit(query[qidx], '&')[[1]][1]
  datatypestr = gsub("[\\(\\)]", "", regmatches(datatypestr, gregexpr("\\(.*?\\)", datatypestr))[[1]])
  
  # Format results
  res = list()
  for (i in seq_along(completeArticle[[qidx]]$entries)) {
    
    pubmedID = completeArticle[[qidx]]$entries[[i]]$`pubmed-id`
    if (is.null(pubmedID)) {
      pubmedID = NA
    }
    
    DOI = completeArticle[[qidx]]$entries[[i]]$`prism:doi`
    if (is.null(DOI)){
      DOI = NA
    }
    
    type = completeArticle[[qidx]]$entries[[i]]$`prism:aggregationType`
    if (is.null(type)){
      type = NA
    }
    
    abstract = completeArticle[[qidx]]$entries[[i]]$`dc:description`
    if (is.null(abstract)){
      abstract = NA
    }
    
    celltype = strsplit(query[qidx], ' & ')[[1]][2]
    celltype = gsub("[\\(\\)]", "", regmatches(celltype, gregexpr("\\(.*?\\)", celltype))[[1]])
    
    nAuthors = as.numeric(completeArticle[[qidx]]$entries[[i]]$`author-count`$`@total`)
    if (length(nAuthors) == 0) {
      nAuthors = NA
      lastAuthor = NA
      firstAuthor = NA
    } else if (nAuthors == 0) {
      lastAuthor = NA
      firstAuthor = NA
    } else if (nAuthors > 0 & nAuthors <= 100){
      lastAuthor = completeArticle[[qidx]]$entries[[i]]$author[[nAuthors]]$authname
      firstAuthor = completeArticle[[qidx]]$entries[[i]]$author[[1]]$authname
    } else if (nAuthors > 100){
      lastAuthor = 'More than 100 authors'
      firstAuthor = completeArticle[[qidx]]$entries[[i]]$author[[1]]$authname
    }
    
    if (completeArticle[[qidx]]$total_results == 0) {
      res[[i]] = NULL
    } else{
      # Create data.frame with paper information
      res[[i]] = data.frame(
        datatype = datatypestr,
        pubmedID = pubmedID,
        URL = completeArticle[[qidx]]$entries[[i]]$`prism:url`,
        title = completeArticle[[qidx]]$entries[[i]]$`dc:title`,
        DOI = DOI,
        journal = completeArticle[[qidx]]$entries[[i]]$`prism:publicationName`,
        date = completeArticle[[qidx]]$entries[[i]]$`prism:coverDisplayDate`,
        type = type,
        nAuthors = nAuthors,
        firstAuthor = firstAuthor,
        lastAuthor = lastAuthor,
        abstract = abstract,
        celltype = celltype,
        citations = completeArticle[[qidx]]$entries[[i]]$`citedby-count`,
        OA = completeArticle[[qidx]]$entries[[i]]$openaccess,
        subtype = ifelse(is.null(completeArticle[[qidx]]$entries[[i]]$subtypeDescription),
                         NA, completeArticle[[qidx]]$entries[[i]]$subtypeDescription),
        query = query[qidx]
      )
    }
  }
  result = rbindlist(res)
  
  # save query results
  queryRes[[qidx]] = result
  
}
scopus = rbindlist(queryRes)
if (any(duplicated(scopus))){
  scopus <- scopus[!duplicated(scopus),]
}
nrow(scopus)

# Save it!
# ===
saveRDS(scopus, file = paste0('~/Desktop/Review-PlateletsRNABPsMotifs/', queryID, '.rds'))
View(scopus)



library(openxlsx)
write.xlsx(scopus,'~/Desktop/Review-PlateletsRNABPsMotifs/firstSearch.xlsx')
