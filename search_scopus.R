# Scopus search
library(data.table)
library(rscopus)
library(dplyr)
library(dotenv)
library(openxlsx)
setwd('~/Desktop/Review-PlateletsRNABPsMotifs/')
load_dot_env('.env')

## START UGR VPN CONNECTION

# Add API KEY
# ===

APIKEY = Sys.getenv('APIKEY') 
rscopus::set_api_key(APIKEY)
queryID = 'plateletsRNAmotifs'

# Build the query
# SCOPUS API SEARCH TIPS
# https://nonprod-devportal.elsevier.com/sc_search_tips.html

# keywords 
## 1- RNA motif, RBP, platelets, megakaryocytes, RNA localization

### TO DO
# Keywords to discard: non-coding RNA (ncRNA), animal models, clinic, medicine
# Only articles from the last 15 years?
# Only articles in Q1 journals

# ===
# sourceSample = list(organism = organism) # ,tissuestate = tissuestate)

#datatype = c('transcriptome','transcriptomic','proteome','proteomic',
#             'RNAseq','RNA-Seq', 'RNABP','binding proteins','binding motifs',
#             'gene expression', 'chip-seq','clip-seq','RNA localization','RNA motif')

datatype = c('RNABP','"rna binding protein"','"binding motifs"', 
             'clip-seq','"RNA localization"','"RNA motif"')
datatypeSTR <- paste(datatype, collapse = ' OR ')
cat(datatypeSTR)

celltype = c('platelet','megakaryo?')
celltypeSTR <- paste(celltype, collapse = ' OR ')

organism = c('human','sapiens')
organismSTR <- paste(organism, collapse = ' OR ')

avoid = c('"non-coding RNA"','ncRNA', '"animal models"', 'clinic', 'medicine')#, 'lncRNA', 'miRNA')
avoidSTR <- paste(avoid, collapse = ' AND NOT ')

yearsAgoLimit <- 15 
pubyearSTR <- paste('PUBYEAR >',2023 - yearsAgoLimit)

languageSTR <- 'LANGUAGE(english)'

doctypeSTR <- 'DOCTYPE(le)'

query <- paste0('TITLE-ABS-KEY((',organismSTR,') AND (',celltypeSTR,') 
                AND (',datatypeSTR, ') AND NOT ',avoidSTR,') AND ', pubyearSTR, 
                ' AND ',languageSTR,' AND NOT ', doctypeSTR)
cat(query)

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
  # datatypestr = strsplit(query[qidx], '&')[[1]][1]
  # datatypestr = gsub("[\\(\\)]", "", regmatches(datatypestr, gregexpr("\\(.*?\\)", datatypestr))[[1]])
  
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
    
    # celltype = strsplit(query[qidx], ' & ')[[1]][2]
    # celltype = gsub("[\\(\\)]", "", regmatches(celltype, gregexpr("\\(.*?\\)", celltype))[[1]])
    
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
    
    keywords <- completeArticle[[qidx]]$entries[[i]]$authkeywords
    if (is.null(keywords)){
      keywords = NA
    }
    
    if (completeArticle[[qidx]]$total_results == 0) {
      res[[i]] = NULL
    } else{
      # Create data.frame with paper information
      res[[i]] = data.frame(
        # datatype = datatypestr,
        # celltype = celltype,
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
        keywords = keywords, 
        citations = completeArticle[[qidx]]$entries[[i]]$`citedby-count`,
        OA = completeArticle[[qidx]]$entries[[i]]$openaccess,
        subtype = ifelse(is.null(completeArticle[[qidx]]$entries[[i]]$subtypeDescription),
                         NA, completeArticle[[qidx]]$entries[[i]]$subtypeDescription),
        query = query[qidx]
      )
    }
  }
  result = rbindlist(res)
  queryRes[[qidx]] = result
}

scopus = rbindlist(queryRes)
if (any(duplicated(scopus))){
  scopus <- scopus[!duplicated(scopus),]
}
nrow(scopus)

# Save it!
# ===
saveRDS(scopus, file = paste0('~/Desktop/Review-PlateletsRNABPsMotifs/', queryID, '3.rds'))
View(scopus)

withminimun <- (grepl('platelet',scopus$title, ignore.case = T) | grepl('megakaryo',scopus$title, ignore.case = T) | 
               grepl('platelet',scopus$abstract, ignore.case = T) | grepl('megakaryo',scopus$abstract, ignore.case = T) | 
               grepl('platelet',scopus$keywords, ignore.case = T) | grepl('megakaryo',scopus$keywords, ignore.case = T) )
sum(withminimun)

scopusGOOD <- scopus[which(withminimun),] 
scopusBAD <- scopus[which(!withminimun),]

nrow(scopus)
nrow(scopusGOOD)
nrow(scopusBAD)

## 1º SEARCH
# write.xlsx(scopus,'~/Desktop/Review-PlateletsRNABPsMotifs/firstSearch.xlsx')
## 2º SEARCH
# write.xlsx(scopus,'~/Desktop/Review-PlateletsRNABPsMotifs/secondSearch.xlsx')
## 3º SEARCH
write.xlsx(scopusGOOD,'~/Desktop/Review-PlateletsRNABPsMotifs/thirdSearch_GOOD.xlsx')
write.xlsx(scopusBAD,'~/Desktop/Review-PlateletsRNABPsMotifs/thirdSearch_BAD.xlsx')


# oldScopus <- read.xlsx('secondSearch.xlsx')
# withminimun <- (grepl('platelet',oldScopus$title, ignore.case = T) | grepl('megakaryo',oldScopus$title, ignore.case = T) | 
#                   grepl('platelet',oldScopus$abstract, ignore.case = T) | grepl('megakaryo',oldScopus$abstract, ignore.case = T) | 
#                   grepl('platelet',oldScopus$keywords, ignore.case = T) | grepl('megakaryo',oldScopus$keywords, ignore.case = T) )
# oldscopusGOOD <- oldScopus[which(withminimun),] 
# oldscopusBAD <- oldScopus[which(!withminimun),]
# 

nrow(scopus)

length(intersect(scopus$pubmedID,oldScopus$pubmedID))

all(na.omit(scopus$pubmedID) %in% na.omit(oldScopus$pubmedID))

all(na.omit(scopusGOOD$pubmedID) %in% na.omit(oldScopus$pubmedID))
all(na.omit(scopusGOOD$pubmedID) %in% na.omit(oldscopusGOOD$pubmedID))


