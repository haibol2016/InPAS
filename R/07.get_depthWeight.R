#' Calculate the depth weight for each sample or each experimental condition
#'
#' Calculate the depth weight for each sample of non-hugeData or each 
#'   experimental condition for hugeData: depth/mean(depth)
#'
#' @param sqlite_db A path to the SQLite database for InPAS, i.e. the output of
#'   [setup_sqlitedb()].
#' @param hugeData A logical(1), indicating whether it is huge data
#' @return A named numeric vector containing depth weight for each sample for 
#'   non-hugeData, or depth weight for each condition if hugeData.
#' @import S4Vectors RSQLite
#' @keywords internal
#' @author Jianhong Ou, Haibo Liu

get_depthWeight <- function(sqlite_db, hugeData) {
  db_conn <- dbConnect(drv = RSQLite::SQLite(), dbname= sqlite_db)
  metadata <- dbReadTable(db_conn, "metadata")
  dbDisconnect(db_conn)
  tags <- metadata$tag
  depth <- metadata$depth
  condition <- metadata$condition
  
  if (hugeData){
    if (length(unique(condition)) != 1){
      depth <- split(depth, condition)
      ## group mean depth
      depth <- sapply(depth, mean)
      tags <- names(depth)
    }
  }
  
  depth.weight <- depth / mean(depth)
  names(depth.weight) <- tags
  depth.weight
}

