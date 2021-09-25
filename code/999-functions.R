
#-----------------------------------------------------------------------------
#' Polygon vertices for plotting computed from an object's linear dimensions
#'
#' This function takes a data frame (one row per object) of object dimensions,
#' where each object is defined by five linear measurements, and computes coords
#' of polygon vertices at given x and y values (taken from the data frame)
#'
#' @param data data frame with at least seven numeric columns
#' @param x column for x position of polygon centers
#' @param y column for y position of polygon centers
#' @param j rows of data frame to subset, default is all
#' @param sf scaling factor to change the size of the polygons, default in 10 (times smaller than actual values)
#' @param top_horizontal_measurement value of linear measurement across the top of the polygon
#' @param main_vertical_measurement value of linear measurement up-and-down the centre of the polygon
#' @param upper_horizontal_measurement value of linear measurement horizontally across the polygon at 25 % of the main_vertical_measurement (starting at the top of main_vertical_measurement)
#' @param middle_horizontal_measurement value of linear measurement horizontally across the polygon at 50 % of the main_vertical_measurement (starting at the top of main_vertical_measurement)
#' @param bottom_horizontal_measurement value of linear measurement horizontally across the polygon at 75 % of the main_vertical_measurement (starting at the top of main_vertical_measurement)
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#' n <- 10
#' my_data <- data.frame(
#'   xcord = rnorm(n), # col for coord for x-axis of poly centers
#'   ycord = rnorm(n), # col for coord for y-axis of poly centers
#'   sf = 10,
#'   main_vertical_measurement = runif(n)      + 5,
#'   top_horizontal_measurement = runif(n)     + 4,
#'   upper_horizontal_measurement = runif(n)   + 4,
#'   middle_horizontal_measurement = runif(n)  + 3,
#'   bottom_horizontal_measurement = runif(n)  + 3
#' )
#'
#' polys <- polygon_shape_by_variables(data = my_data,
#'                                     x = 'xcord',
#'                                     y = 'ycord',
#'                                     main_vertical_measurement = 'main_vertical_measurement',
#'                                     top_horizontal_measurement = 'top_horizontal_measurement',
#'                                     upper_horizontal_measurement = 'upper_horizontal_measurement',
#'                                     middle_horizontal_measurement= 'middle_horizontal_measurement',
#'                                     bottom_horizontal_measurement = 'bottom_horizontal_measurement')
#'
#' p <- ggplot() +
#'   geom_polygon(data=polys, aes(x=x, y=y, group = fill), colour="black", alpha = 0.5) +
#'   coord_equal() +
#'   theme_bw()
#' #p
#'
#'
polygon_shape_by_variables <- function(data,       # some data frame
                                       x,          # col for coord for x-axis of poly centers
                                       y,          # col for coord for y-axis of poly centers
                                       j = seq(1, nrow(data)),  # rows to subset
                                       sf = 10,    # scale factor
                                       main_vertical_measurement,
                                       top_horizontal_measurement, # cols to use to compute poly verts
                                       upper_horizontal_measurement,
                                       middle_horizontal_measurement,
                                       bottom_horizontal_measurement){

  # object to store results of loop (list avoids slow growing of dataframe)
  df <- vector("list", length = length(j))

  data <- as.data.frame(data)

  for(k in seq_along(j)){

    # the specimen (row) to compute polygon vertices and x-y location for
    i <- data[k, ]

    # x and y axis coords for the polygon
    x_ <- data[, x][k]
    y_ <- data[, y][k]

    top_left_x <- x_ - unlist(i[, top_horizontal_measurement]/(2 * sf))
    top_left_y <- y_ + unlist(i[, main_vertical_measurement]/(2 * sf))

    mid_upper_left_x <- x_  - unlist(i[, upper_horizontal_measurement]/(2 * sf))
    mid_upper_left_y <- y_  + unlist(i[, main_vertical_measurement]/(4 * sf))

    mid_lower_left_x <- x_  - unlist(i[, middle_horizontal_measurement]/(2 * sf))
    mid_lower_left_y <- y_  - 0

    bottom_left_x <- x_  - unlist(i[, bottom_horizontal_measurement]/(2 * sf))
    bottom_left_y <- y_  - unlist(i[, main_vertical_measurement]/(2 * sf))

    distal_tip_x <- x_  - 0
    distal_tip_y <- y_  - unlist(i[, main_vertical_measurement]/(sf))

    bottom_right_x <- x_  + unlist(i[, bottom_horizontal_measurement]/(2 * sf))
    bottom_right_y <- y_  - unlist(i[, main_vertical_measurement]/(2 * sf))

    mid_lower_right_x <- x_  + unlist(i[, middle_horizontal_measurement]/(2 * sf))
    mid_lower_right_y <- y_  - 0

    mid_upper_right_x <- x_  + unlist(i[, upper_horizontal_measurement]/(2 * sf))
    mid_upper_right_y <- y_  + unlist(i[, main_vertical_measurement]/(4 * sf))

    top_right_x <- x_  + unlist(i[, top_horizontal_measurement]/(2 * sf))
    top_right_y <- y_  + unlist(i[, main_vertical_measurement]/(2 * sf))

    one_specimen <- data.frame(
      # we are going around the flake in a circle
      x=c(top_left_x,
          mid_upper_left_x,
          mid_lower_left_x,
          bottom_left_x,
          distal_tip_x,
          bottom_right_x,
          mid_lower_right_x,
          mid_upper_right_x,
          top_right_x),
      y=c(top_left_y,
          mid_upper_left_y,
          mid_lower_left_y,
          bottom_left_y,
          distal_tip_y,
          bottom_right_y,
          mid_lower_right_y,
          mid_upper_right_y,
          top_right_y),
      fill = as.character(rep(k, 9) )
    )

    df[[k]] <- one_specimen

  }

  # convert list to dataframe
  df <- data.frame(do.call(rbind, df))

  return(df)

}

# based on https://github.com/aphanotus/borealis/blob/master/R/create.tps.R

create_tps <- function(raw, # input data frame
                       output.filename # output TPS filename we want to create
                       ){

  acceptable.ID.column.names <- c("id","ids","specimen","specimen id","specimen.id","specimen_id","specimen.ids","specimen_ids","sample","sample id","sample.id","sample_id","sample.ids","sample_ids","id number","id.number","id_number")
  ID.col <- which(names(raw) %in% acceptable.ID.column.names)
  specimen.number <- nrow(unique(raw[,ID.col]))

  acceptable.LM.column.names <- c("LM","lm","landmark","landmarks","landmark.number","landmark_number")
  LM.col <- which(names(raw) %in% acceptable.LM.column.names)
  landmark.number <- dim(raw)[1] / specimen.number

  id.factors <- NULL

sink(output.filename)
# Loop for each specimen
for (i in 1:specimen.number) {
  cat(paste0('LM=',landmark.number,'\n'))

  # Nested loop for each landmark
  for (j in 1:landmark.number) {
    x <- landmark.number*(i-1) + j
    cat(paste0(as.numeric(raw$x[x]),' ',as.numeric(raw$y[x]),'\n'))
  } # End nested loop for each landmark

  x <- landmark.number*(i-1) + 1
  id.text <- paste0(sub('\\.','_',trimws(as.character(raw[x,ID.col]))))
  if (length(id.factors)>0) {
    id.text <- paste0(id.text,separator,paste(trimws(as.character(raw[x,id.factors])),collapse = separator) )
  }
  cat(paste0('ID=',id.text,'\n'))
  cat('\n')
} # End loop for each specimen

# Close the output file
sink()
}


