#' search NCBI pubmed for wanderful papers
#'
#' @param keyword the query statement used to search in ncbi pubmed website
#'
#' @return a data frame with IF2017 included
#' @author hanguangchun@gmail.com
#' @keywords searchNCBI
#' @export

searchNCBI = function(keyword, db = 'pubmed',out = NULL) {
  require0('rentrez')
  web_env_search <- entrez_search(db=db, term = keyword, use_history=TRUE)
  wh <- web_env_search$web_history
  xml_res <- entrez_fetch(db = db, web_history = wh, retmax = 1000,rettype = 'xml')
  res <- parse_pubmed_xml(xml_res)
  authors <- lapply(res,'[[','authors')
  authors1 = unlist(lapply(authors, function(e) paste(e, collapse = ';')))
  article_keywords <- lapply(res, '[[', 'key_words') 
  keywords1 = unlist(lapply(article_keywords, function(e) paste(e, collapse = ';')))
  journal <- lapply(res, '[[', 'journal') 
  journal1 = unlist(lapply(journal, function(e) paste(e, collapse = ';')))
  titles <- lapply(res, '[[', 'title') 
  titles1 = unlist(lapply(titles, function(e) paste(e, collapse = ';')))
  year <- lapply(res, '[[', 'year')
  year1 = unlist(lapply(year, function(e) paste(e, collapse = ';')))
  abstr <- lapply(res, '[[', 'abstract')
  abstr1 = unlist(lapply(abstr, function(e) paste(e, collapse = ';')))
  pmid <- lapply(res, '[[', 'pmid')
  pmid1 = unlist(lapply(pmid, function(e) paste(e, collapse = ';')))

  require0('dplyr')
  require0('readr')
  require0('readxl')

  ifdat$`Full Journal Title` = toupper(ifdat$`Full Journal Title`)
  ifdat$`Journal Impact Factor` = as.numeric(as.character(
                       ifdat$`Journal Impact Factor`))
  match(toupper(journal1), ifdat$`Full Journal Title`) -> pos
  if_score = ifdat[pos,]$`Journal Impact Factor`

  outputTbl = tbl_df(data.frame(
    title = titles1,
    abstr = abstr1,
    year = year1,
    ImpactFactor2017 = if_score,
    keywords = keywords1,
    journal = journal1,
    authors = authors1,
    pmid = pmid1,
    stringsAsFactors = F
    )) 
  outputTbl = outputTbl %>% arrange(desc(ImpactFactor2017))
  if(is.null(out)) {
    out = paste0(gsub('\\/',' or ',keyword), '.searchNCBI.tsv')
  }
  write_tsv(outputTbl, out)
  outputTbl
}
