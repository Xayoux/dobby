# Check type of the parameter -----------------------------------------------
#' @title Check if Parameter is NULL or Character
#'
#' @param param Value of the parameter
#' @param name_param Character : name of the parameter
#' @return Stop or nothing
#' @export
.check_null_character <- function(param, name_param){
  if (!is.null(param) & !is.character(param)){
    class_param <- class(param)
    stop(glue::glue("{name_param} must be NULL or a character, not a {class_param}."))
  }
}


#' @title Check if Parameter is NULL or list or Character
#'
#' @param param Value of the parameter
#' @param name_param Character : name of the parameter
#' @return Stop or nothing
#' @export
.check_null_list_character <- function(param, name_param){
  if (!is.null(param) & !is.list(param) & !is.character(param)){
    class_param <- class(param)
    stop(glue::glue("{name_param} must be NULL or a list or a character, not a {class_param}."))
  }
}


#' @title Check if Parameter is Character
#'
#' @param param Value of the parameter
#' @param name_param Character : name of the parameter
#' @return Stop or nothing
#' @export
.check_character <- function(param, name_param){
  if (!is.character(param)){
    class_param <- class(param)
    stop(glue::glue("{name_param} must be a character, not a {class_param}."))
  }
}


#' @title Check if Parameter is a Logical
#'
#' @param param Value of the parameter
#' @param name_param Character : name of the parameter
#' @return Stop or nothing
#' @export
.check_logical <- function(param, name_param){
  if(!is.logical(param)){
    class_param <- class(param)
    stop(glue::glue("{name_param} must be a logical, not a {class_param}."))
  }
}


#' @title Check if Parameter is a NULL or Numeric
#'
#' @param param Value of the parameter
#' @param name_param Character : name of the parameter
#' @return Stop or nothing
#' @export
.check_null_numeric <- function(param, name_param){
  if(!is.null(param) & !is.numeric(param)){
    class_param <- class(param)
    stop(glue::glue("{name_param} must be null or numeric, not a {class_param}."))
  }
}


#' @title Check if Parameter is a Numeric
#'
#' @param param Value of the parameter
#' @param name_param Character : name of the parameter
#' @return Stop or nothing
#' @export
.check_numeric <- function(param, name_param){
  if(!is.numeric(param)){
    class_param <- class(param)
    stop(glue::glue("{name_param} must be a numeric, not a {class_param}."))
  }
}


# Check length of parameter -------------------------------------------------
#' @title Check if Parameter is Length 1
#'
#' @param param Value of the parameter
#' @param name_param Character : name of the parameter
#' @param ln Numeric indicating the length the parameter must have exactly
#' @return Stop or nothing
#' @export
.check_length <- function(param, name_param, ln = 1){
  length_param <- length(param)
  if (length_param != ln){
    stop(glue::glue("{name_param} must be length {ln}, not length {length_param}."))
  }
}


# Check column exist in data ------------------------------------------------
pluralise <- function(n, singular, plural) {
  if (n == 1) {
    singular
  } else {
    plural
  }
}

#' @title Check if Columns Provided are Present in a Dataframe
#'
#' @param df Dataframe
#' @param name_df Name of the Dataframe (of the parameter)
#' @param columns Character : Columns to be tested
#' @return Stop or nothing
#' @export
.check_var_exist <- function(df, name_df, columns){
  is_column_present <- rlang::has_name(df, columns)
  if (FALSE %in% is_column_present){
    columns_absent <- columns[which(is_column_present == FALSE)]
    message <-
      pluralise(
        n = length(columns_absent),
        singular = glue::glue("Column \"{columns_absent}\" is not in {name_df}."),
        plural = glue::glue("Columns \"{columns_absent}\" are not in {name_df}")
      )
    stop(message)
  }
}
