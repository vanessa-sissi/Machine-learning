change_resolution <- function(image_df, x_size)
{
  ## change_resolution(image_df, x_size) subsamples an image to produce
  ## a lower resolution image. Any non-coordinate columns in the data
  ## frame are summarized with their most common value in the larger
  ## grid cell.
  ##
  ## Input:
  ## - image_df: A data frame in wide format. The x-coordinate column MUST
  ##             be named 'x' and the y-coordinate column MUST be named 'y'.
  ##             Further columns have no naming restrictions.
  ## - x_size:   The number of cells in the x-direction. The number of cells
  ##             in the vertical direction will be computed to maintain the 
  ##             perspective. There is no guarantee that the exact number
  ##             of cells in the x-direction is x_size
  ##
  ## Output:
  ## - A data frame with the same column names as image_df, but with fewer 
  ##   entries that corresponds to the reduced resolution image.
  ##
  ## Example:
  ##   library(imager)
  ##   library(dplyr)
  ##   fpath <- system.file('extdata/Leonardo_Birds.jpg',package='imager') 
  ##   im <- load.image(fpath)
  ##   im_dat<- as.data.frame(im,wide = "c") %>% rename(R = c.1, G = c.2, B = c.3) %>%
  ##            select(x,y,R,G,B)
  ##   agg_image <- change_resolution(im_dat, 50)
  
  if(!require(sp)) {
    stop("The sp packages must be installed. Run install.packages(\"sp\") and then try again.")
  }
  if(!require(dplyr)) {
    stop("The dplyr packages must be installed. Run install.packages(\"dplyr\") and then try again.")
  }
  
  sp_dat <- image_df 
  gridded(sp_dat) = ~x+y
  
  persp = (gridparameters(sp_dat)$cells.dim[2]/gridparameters(sp_dat)$cells.dim[1])
  y_size = floor(x_size*persp)
  orig_x_size = gridparameters(sp_dat)$cells.dim[1]
  orig_y_size = gridparameters(sp_dat)$cells.dim[2]
  
  x_res = ceiling(orig_x_size/x_size)
  y_res = ceiling(orig_y_size/y_size)
  
  gt = GridTopology(c(0.5,0.5), c(x_res, y_res),
                    c(floor(orig_x_size/x_res), floor(orig_y_size/y_res)))
  SG = SpatialGrid(gt)
  agg = aggregate(sp_dat, SG, function(x) names(which.max(table(x)))[1] )
  agg@grid@cellsize <- c(1,1)
  df <- agg %>% as.data.frame %>% rename(x = s1, y = s2)  %>% select(colnames(image_df))
  
  return(df)
  
}



process_image <- function(image_file_name, k_list){
  ## process_image(image_file_name, k_list) processes the image by k-means
  ## which return k clusters information.
  ##
  ## Input:
  ##  - image_file_name - a PNG or JPEG image.
  ##  - k_list - the list of numbers of centres in the clustering
  ##
  ## Output:
  ##  - cluster_info: A list or tibble of information derived from the k_means
  ##    for k in k_list. 
  ##
  ## Example:
  ##  clusters = process_image('mickey.jpg',c(2:8))

  
  
  if(!require(dplyr)) {
    stop("The dplyr packages must be installed. Run install.packages(\"dplyr\") and then try again.")
  }
  if(!require(purrr)) {
    stop("The purrr packages must be installed. Run install.packages(\"purrr\") and then try again.")
  }
  if(!require(tidymodels)) {
    stop("The tidymodels packages must be installed. Run install.packages(\"tidymodels\") and then try again.")
  }
  im <- imager::load.image(image_file_name)
  tidy_dat <- as.data.frame(im, wide = "c") %>% rename(R = c.1, G = c.2, B = c.3)
  dat <- select(tidy_dat,c(-x,-y))
  kclusts <-
    tibble(k = k_list) %>%
    tibble(dimension = list(c(x = max(tidy_dat$x), y = max(tidy_dat$y)))) %>%
    mutate(
      kclust = map(k, ~kmeans(x = dat , centers = .x, nstart=4)),
      glanced = map(kclust, glance),
    ) 
  clusterings <-
    kclusts %>%
    unnest(cols = c(glanced))
  
  return(clusterings)
}

scree_plot <- function(cluster_info){
  ## scree_plot(cluster_info) produces the scree plot of given cluster info
  ## 
  ## Input:
  ##  - cluster_info - The output of process_image i.e. cluster info
  ##
  ## Output:
  ##  - It produces the scree plot of give cluster info
  ## 
  ## Example:
  ##  clusters = process_image('mickey.jpg',c(2:8))
  ##  scree_plot(clusters)
  
  nclust = length(cluster_info$k)
  ratio = rep(NA, nclust-1)
  for (kk in 2:nclust) {
    ratio[kk-1] = cluster_info$tot.withinss[kk]/cluster_info$tot.withinss[kk-1]
  }
  plot_data <- data.frame(k = cluster_info$k[2:nclust],ratio)
  ggplot(plot_data, aes(x=k, y = ratio)) + geom_line()
}

color_strips <- function(cluster_info){
  ## colour_strips(cluster_info) produces color strips with the DMC
  ## color closest to the cluster centre color (cluster_info)
  ##
  ## Input:
  ##  - cluster_info - The output of process_image i.e. cluster info
  ##
  ## Output:
  ##  - Color strips (with plot) with the DMC color closest to the cluster center color
  ##    and return dcm and cluster center color for all k in cluster_info
  ##
  ## Example:
  ##  clusters = process_image('mickey.jpg',c(2:8))
  ##  strips = color_strips(clusters) ## generate strips data
  ##  strips$plot[[which(2==strips$k)]] ## plot out the color strips when k=2
  if(!require(dmc)) {
    stop("The dmc packages must be installed. Run devtools::install_github(\"sharlagelfand/dmc\") and then try again.")
  }
  
  cluster_info$dmc = list(tibble())
  cluster_info$plot = list(ggplot())
  for(i in 1:length(cluster_info$k)){
    temp = tidy(cluster_info$kclust[[i]]) %>% mutate(col = rgb(R, G, B)) %>% 
      rowwise() %>% mutate(col, dmc::dmc(col, visualize = FALSE))
    cluster_info$dmc[[i]] = temp
    n_col = length(temp$hex)
    rect_dat <- tibble(x1 = c(0:(n_col-1)), x2 = c(1:n_col), y1 = rep(0,n_col),
                       y2 =rep(1,n_col), colour = temp$hex)
    cluster_info$plot[[i]] = rect_dat %>% ggplot() + coord_fixed() +
      geom_rect(aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=colour), color="black") +
      geom_text(aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2, label=colour), size=24/n_col) +
      scale_fill_manual(values = rect_dat$colour)+ theme_void() + theme(legend.position = "none")
    
     
  }
  return(cluster_info %>% select(k, dmc, plot))
}

make_pattern <- function(cluster_info, k, x_size, black_white = FALSE, background_colour = NULL){
  ## make_pattern(cluster_info, k, x_size, black_white = FALSE, background_colour = NULL)
  ## plots the pattern.
  ##
  ## Input:
  ##  - cluster_info - The output of process_image i.e. cluster info
  ##  - k - The chosen cluster size
  ##  - x_size - The (approximate) total number of possible stitches 
  ##    in the horizontal direction
  ##  - black_white - (logical) Print the pattern in black and white (TRUE) 
  ##    or colour (FALSE, default)
  ##  - background_colour - The colour of the background, which should not 
  ##    be stitched in the pattern. (Default is to not have a colour)
  ##
  ## Output:
  ##  - It produces a cross-stitch pattern that can be followed, 
  ##    complete with a legend that has thread colour, and a guide grid.
  ##
  ## Example: 
  ##  clusters = process_image('mickey.jpg',c(2:8))
  ##  make_pattern(clusters, 2, 50)
  if(!require(scales)) {
    stop("The scales) packages must be installed. Run install.packages(\"scales)\") and then try again.")
  }
  
  strips = color_strips(cluster_info)
  cluster = cluster_info$kclust[[which(cluster_info$k == k)]]
  strip = strips$dmc[[which(strips$k == k)]] %>% 
          mutate(color_thread = paste(name, '(',dmc, ')', sep=''))
  dimension = cluster_info$dimension[[1]]
  nx = dimension[[1]]
  ny = dimension[[2]]
  
  x = rep(1:nx, ny)
  y = as.vector(t(matrix(x, nrow = nx)))
  cls = factor(cluster$cluster)
  dat = change_resolution(tibble(x, y, cls) %>% rename(cluster = cls), x_size)
  dat = dat %>% left_join(strip %>% select(cluster, col, color_thread))
  
  xx = ceiling(max(dat$x))
  yy = ceiling(max(dat$y))
  p <-  dat %>% ggplot() 
  if(black_white){
    p = p + geom_point(aes(x=x, y = y, shape = color_thread))
  }else{
    p = p + geom_point(aes(x=x, y = y, shape = color_thread, col = color_thread))
  }
  p + scale_y_reverse() + 
    scale_color_manual(values = strip %>% select(color_thread, col)%>% deframe) +
    theme(
      panel.background = element_rect(fill = background_colour),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank()) + 
    geom_vline(xintercept=seq(0, xx, round(xx/5)), color="black", size=1.25) +
    geom_hline(yintercept=seq(0, yy, round(xx/5)), color="black", size=1.25) 
}




