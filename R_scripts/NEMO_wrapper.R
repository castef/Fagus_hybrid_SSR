
#### Script to create files for NEMO simulation ####
# Camilla Stefanini #

# the functions simulate individuals of X species clustered into specified spatial polygons
# the species information is retrieved from an external dataframe, for instance, NH results here
# the genetic information is from genpop files (but .txt format)

### workflow ###
# 1. combine NewHybrids species (obs_file) information and genepop file --> locus list output
# 2. simulate new individuals with allele frequencies from real data, create populations
# 3. create the FSTAT file for the genetic simulation, with genotype taken from real genetic data source (genepop format)
# 4. create carrying_capacity matrices (patch_nbfem) with defined carrying capacity for patches of the grid
# 5. create quanti_init_freq matrix for quanti_trait
# 6. create seed and pollen dispersal matrices

# path where to create and store the input files
input_path = "~/NEMO/Rwrapper_files/"
setwd(input_path)

# ------- 1. create locus list -----------

get_locus_list <- function(
    genepop_file,
    obs_file,
    id_col = "SampleID",
    species_col = "NH_class_t80", ## name of the column that specify the species (or NH class)
    filter_generation = "Adult",
    species_labels = c("Pure orientalis", "Pure sylvatica"),
    species_names = c("Orientalis", "Sylvatica"),
    ncode = 3L,
    plot_colors = c("Orientalis" = "#482173FF", "Sylvatica" = "#ffcc00")
) {
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(adegenet)
  library(ggplot2)
  
  # --- load species information
  obs <- read.csv(obs_file)
  
  # --- load Genepop file
  lines <- readLines(genepop_file)
  loci_end <- which(tolower(lines) == "pop")[1]
  header <- lines[1]
  loci <- lines[2:(loci_end - 1)]
  genotypes <- lines[(loci_end + 1):length(lines)]
  ind_ids <- sapply(strsplit(genotypes, ","), function(x) trimws(x[1]))
  
  # --- keep adults per species
  get_indivs <- function(label) {
    obs %>%
      filter(generation == filter_generation, .data[[species_col]] == label) %>%
      pull(all_of(id_col))
  }
  
  inds_sp1 <- get_indivs(species_labels[1])
  inds_sp2 <- get_indivs(species_labels[2])
  
  if (length(inds_sp1) == 0 || length(inds_sp2) == 0) {
    stop("No individuals found for one or both species. Check your species labels and columns.")
  }
  
  # --- extract genotype lines
  sp1_data <- genotypes[ind_ids %in% inds_sp1]
  sp2_data <- genotypes[ind_ids %in% inds_sp2]
  
  # --- build new genepop
  new_genepop <- c(header, loci, "Pop", sp1_data, "Pop", sp2_data)
  tmp_gen <- tempfile(fileext = ".gen")
  writeLines(new_genepop, tmp_gen)
  
  # --- re read as genind
  genind <- read.genepop(tmp_gen, ncode = ncode, quiet = TRUE)
  genpop <- genind2genpop(genind)
  allele_freq <- as.data.frame(makefreq(genpop))
  rownames(allele_freq) <- species_names
  
  # --- tidy format
  allele_freq_df <- allele_freq %>%
    tibble::rownames_to_column("Species") %>%
    pivot_longer(cols = -Species, names_to = c("Locus", "Alleles"), names_sep = "\\.", values_to = "Freq")
  
  # --- allele frequencies plot
  plot <- ggplot(allele_freq_df, aes(x = Alleles, y = Freq, fill = Species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
    facet_wrap(~ Locus, scales = "free_x") +  
    scale_fill_manual(values = plot_colors) +
    labs(x = "Alleles", y = "Frequency", title = "Allele Frequencies per Locus") +
    theme_minimal()
  
  # --- locus list (allele frequencies for each locus and species)
  locus_list <- split(allele_freq_df, allele_freq_df$Locus)
  
  return(list(
    locus_list = locus_list,
    allele_freq_df = allele_freq_df,
    plot = plot
  ))
}



### example
allele_freq <- get_locus_list(
  genepop_file = "Waldi_genepop_all.txt",
  obs_file = "Waldi_indiv_info_class.csv",
  id_col = "SampleID",
  species_col = "NH_class_t80",
  filter_generation = "Adult",
  species_labels = c("Pure orientalis", "Pure sylvatica"),
  species_names = c("Orientalis", "Sylvatica")
)

print(allele_freq$plot)

# ------- 2. simulate populations -----------------
# functions
# 1. simulate_individual()  -----> Simulates a single individual from a locus list of a specific species (previous step)
# 2. simulate_population()  -----> Creates a whole population of N individuals for each species
# 3. assign_spatial_positions() -----> Assigns simulated individuals to spatial positions inside a polygon
# 4. simulate_spatial_population() ----->  combine genotypes and positions
# 5. simulate_species_population() -----> loop
# 6. plot the results

### NOTE: the simulation creates individuals with species (or population) code prefix (3 digits) and then ID, eg. "syl_sim_12304"

# required: `terra`, `plyr`, `dplyr`

simulate_individual <- function(locus_list, species) {
  indiv <- data.frame(matrix(ncol = 2, nrow = length(locus_list)))      
  colnames(indiv) <- c("Allele1", "Allele2")
  row.names(indiv) <- names(locus_list)
  
  for (locus_name in names(locus_list)) {
    species_data <- subset(locus_list[[locus_name]], Species == species)
    indiv[locus_name, 1:2] <- sample(
      species_data$Alleles, size = 2, replace = TRUE, prob = species_data$Freq
    )
  }
  
  indiv$Genotype <- paste0(indiv$Allele1, indiv$Allele2)
  return(indiv)
}

simulate_population <- function(n, locus_list, species) {
  individuals <- lapply(seq_len(n), function(x) simulate_individual(locus_list, species))
  genotypes <- lapply(individuals, function(indiv) indiv$Genotype)
  pop <- data.frame(do.call(cbind, genotypes))    
  row.names(pop) <- names(locus_list)
  colnames(pop) <- paste0(substr(tolower(species), 1, 3), "_sim_", seq_len(n))
  return(pop)
}

assign_spatial_positions <- function(grid_file, polygon_coords, crs, precision) {
  library(terra)
  library(dplyr)
  library(plyr)
  
  grid <- vect(grid_file)
  polygon <- vect(polygon_coords, type = "polygons", crs = crs)
  centroids <- centroids(grid)
  
  in_poly <- relate(centroids, polygon, "intersects")
  selected <- centroids[in_poly]
  df <- as.data.frame(selected, geom = "XY")
  
  df$X2 <- plyr::round_any(df$x, precision)
  df$Y2 <- plyr::round_any(df$y, precision)
  
  result <- df %>%
    group_by(X2, Y2) %>%
    summarise(
      x = x[which.min(abs(x - median(x)))],
      y = y[which.min(abs(y - median(y)))],
      .groups = "drop"
    )
  
  vect(result, geom = c("x", "y"), crs = crs)
}

simulate_spatial_population <- function(
    n,
    species,
    locus_list,
    grid_file,
    polygon_coords,
    crs = "EPSG:3035",
    precision = 10
) {
  library(terra)
  
  # get spatial positions
  points <- assign_spatial_positions(
    grid_file = grid_file,
    polygon_coords = polygon_coords,
    crs = crs,
    precision = precision
  )
  
  # keep only N positions
  points <- points[1:min(n, length(points))]
  
  # simulate genotypes
  genotypes <- simulate_population(length(points), locus_list, species)
  
  # get coords and patch IDs
  coords <- crds(points)
  grid <- vect(grid_file)
  patches_df <- terra::extract(grid, points)
  patch_ids <- patches_df[,2] # 2nd col should be patch ID
  
  # combine everything
  df <- cbind(t(genotypes), coords, patch_id = patch_ids)
  df <- tibble::rownames_to_column(as.data.frame(df), "SampleID")
  
  return(df)
}

#### wrapper 
simulate_species_population <- function(
    species_list,               # Vector with species names
    n_individuals,              # Named vector: number of individuals per species
    polygons,                   # Named list of polygon coordinate matrices (species = names)
    locus_list,                 # Locus list
    grid_file,                  # Grid shapefile path
    crs = "EPSG:3035",
    precisions = NULL           # Optional named vector: precision per species
) {
  library(dplyr)
  all_sims <- list()
  
  for (sp in species_list) {
    cat("Simulating:", sp, "\n")
    
    poly_coords <- polygons[[sp]]
    prec <- if (!is.null(precisions)) precisions[[sp]] else 10
    
    sim_df <- simulate_spatial_population(
      n = n_individuals[[sp]],
      species = sp,
      locus_list = locus_list,
      grid_file = grid_file,
      polygon_coords = poly_coords,
      crs = crs,
      precision = prec
    )
    
    all_sims[[sp]] <- sim_df
  }
  
  full_df <- do.call(rbind, all_sims)
  return(full_df)
}

## this works if the individuals have a prefix with species code 

plot_simulated_individuals <- function(sim_df, polygons, grid_file, crs = "EPSG:3035") {
  library(terra)
  library(ggplot2)
  
  grid <- vect(grid_file)
  sim_df[c("x", "y")] <- lapply(sim_df[c("x", "y")], as.numeric)
  
  # convert simulated individuals to spatial points
  points <- vect(sim_df, geom = c("x", "y"), crs = crs)
  
  # get species column from SampleID prefix
  sim_df$Species <- substr(sim_df$SampleID, 1,3)
  
  # spatial object with species
  points <- vect(sim_df, geom = c("x", "y"), crs = crs)
  points$Species <- sim_df$Species
  
  # convert all to data.frames for plotting 
  grid_sf <- sf::st_as_sf(grid)
  points_df <- as.data.frame(points, geom = "XY")
  
  # convert polygon coordinates to spatial for plotting
  poly_list <- lapply(names(polygons), function(sp) {
    data.frame(polygons[[sp]]) |>
      setNames(c("x", "y")) |>
      mutate(Species = sp)
  })
  poly_df <- do.call(rbind, poly_list)
  poly_df$Group <- rep(names(polygons), sapply(polygons, nrow))
  
  ggplot() +
    geom_sf(data = grid_sf, fill = NA, color = "grey70") +
    geom_polygon(data = poly_df, aes(x = x, y = y, group = Group, fill = Species), alpha = 0.2) +
    geom_point(data = points_df, aes(x = x, y = y, color = Species), size = 2) +
    coord_sf() + guides(col = "none")+
    theme_minimal() +
    labs(title = "Simulated individuals", x = "Easting", y = "Northing") 
}


### example

species <- c("Orientalis", "Sylvatica")
n_ind <- c(Orientalis = 52, Sylvatica = 97)
precisions <- c(Orientalis = 10.96, Sylvatica = 13.23)
locus_list_a <- allele_freq$locus_list

polygons <- list(
  Orientalis = matrix(c(
    4126715, 2840045,
    4126800, 2840120,
    4126705, 2840152,
    4126715, 2840045
  ), ncol = 2, byrow = TRUE),
  
  Sylvatica = matrix(c(
    4126800, 2840120,
    4126870, 2840200,
    4126705, 2840260,
    4126705, 2840152, 
    4126800, 2840120
  ), ncol = 2, byrow = TRUE)
)


simulation_test <- simulate_species_population(
  species_list = species,
  n_individuals = n_ind,
  polygons = polygons,
  locus_list = locus_list_a,
  grid_file = "Waldi_grid_4m_5396patches.shp",
  precisions = precisions
)

plot_simulated_individuals(
  sim_df = simulation_test,
  polygons = polygons,
  grid_file = "Waldi_grid_4m_5396patches.shp",
  crs = "EPSG:3035"
)

# save
write.table(simulation_test,"wrapper_test_step2.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


# ------- 3. build FSTAT files ------

### individuals_file (from previous step) must contain in this order
# SampleID, 
# Genotype (loci coded like 123123), 
# x,y, coordinates  
# patch_id in the simulation grid

build_fstat_file <- function(path = NULL, individuals_file, grid, output_file, age_range = age_range, stage = stage, sex = 1, pedigree = 0, max_n_allele,digits_allele ) {
  
  library(dplyr)
  library(readr)
  
  # read simulated/real population
  fstat_df <- read_table(individuals_file, col_names = T)
  
  # get loci columns 
  loci_start <- which(names(fstat_df) == "SampleID") + 1
  loci_end <- which(names(fstat_df) == "x") - 1
  loci_names <- names(fstat_df)[loci_start:loci_end]
  n_loci <- length(loci_names)
  
  # create FSTAT extended version
  # remove the coordintes
  fstat_ext <- subset(fstat_df, select = -c(x, y))
  
  # set other columns
  fstat_ext$Age <- sample(age_range[1]:age_range[2], nrow(fstat_ext), replace = TRUE)
  fstat_ext$Stage <- stage
  fstat_ext$Sex <- sex
  fstat_ext$Pedigree <- pedigree
  # reorder columns
  fstat_ext <- fstat_ext %>% 
    relocate(Stage, .after = last_col()) %>% 
    relocate(Age, .after = Stage) %>% 
    relocate(Sex, .after = Age) %>% 
    relocate(Pedigree, .after = Sex) %>% 
    relocate(patch_id, .after = Pedigree)
  
  # remove first column (individual names)
  fstat_ext_final <- fstat_ext[, -1]
  
  # replace '000000' with '001001' (ONLY FOR TRUE INDIVIDUALS; SIMULATED DONT HAVE MISSING INFORMATION)
  fstat_ext_final <- apply(fstat_ext_final, 2, function(x) gsub("000000", "001001", x))
  fstat_ext_final <- as.data.frame(fstat_ext_final)
  
  # set first column as patch_id
  fstat_ext_final = cbind(patch_id = fstat_ext_final[, ncol(fstat_ext_final)], fstat_ext_final)
  
  # add header (n patches, col number - 1, max number of alleles per locus, number of digits for each genotype)
  header <- c(
    paste(ncell(grid),ncol(fstat_ext_final)-1, max_n_allele, digits_allele),
    loci_names,"stage","age","sex","ped","origin")
  
  # write header and after the data
  writeLines(header, output_file)
  write.table(fstat_ext_final, output_file, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  return(fstat_ext_final)
}


#### example
build_fstat_file(
  age_range = c(50:51), # range of age of individuals
  stage = 3, # NEMO age stage
  max_n_allele = 238, # max number of alleles 
  digits_allele = 3, # number of allele digits
  grid = vect("Waldi_grid_4m_5396patches.shp"),
  individuals_file =  "wrapper_test_step2.txt",
  output_file = "wrapper_test_step3.txt"
)

# ------- 4. create carrying capacity matrix (TO DO ) ------

# creating matrices where empty patches have carrying capacity 0, and patches that can be colonized can have a defined carrying capacity
# here, the presence is defined by a tif file

#  write matrix in NEMO format
write.matrix.nemo <- function(mat, outfile) {
  rows <- dim(mat)[1]
  cols <- dim(mat)[2]
  cat("{", file = outfile)
  
  for(i in 1:rows) {
    cat("{", file = outfile, append = TRUE)
    cat(mat[i,], sep = ",", file = outfile, append = TRUE)
    cat("}\n", file = outfile, append = TRUE)
  }
  
  cat("}\n", file = outfile, append = TRUE)
}

# ALLENWILLER #

# loading species classification (from remote sensing, 0 = beech) and project to metric system
class <- rast("~/Allen_classification.tif")
class <- project(x=class,y= "epsg:3035")
plot(class)
# load simulation grid
grid <- vect("Allenwiller_grid_r4m_9844patches.shp",crs=crs("+init=epsg:3035"))

# for each cell of the grid, extract the value of classification in a dataframe (min to be conservative)
class_extract = terra::extract(class, grid, fun = "min", cells = T, xy = T, layer = 1)
nrow(class_extract)

# (OPTIONAL) check the highest number of seedlings counted into circle plots
seed_n <- read.csv("Allenwiller_circleplot_info.csv")
max(seed_n$CountedC1.new.)

# define vector of carrying capacity values to test
car_cap = c(200,500,800)

# set the carrying capacity to 0 on the patches at the borders, to avoid edge effects: 
grid_ext <- ext(grid)
coords <- terra::centroids(grid) |> crds()
border <- 20
border_patch <- coords[,1] <= (grid_ext[1] + border) | 
  coords[,1] >= (grid_ext[2] - border) |
  coords[,2] <= (grid_ext[3] + border) | 
  coords[,2] >= (grid_ext[4] - border)

plot(grid, col = border_patch)

# for each value to test, create a vector of max carrying capacity observed in patches with presence of beech (0), and 0 where there is not beech
for(i in seq_along(car_cap)){
  # create column for each value to test
  class_extract[paste0("car_cap_", car_cap[i])] <- ifelse(class_extract$value <= 0.5, car_cap[i],0)
  # set border patches to 0
  class_extract[border_patch, paste0("car_cap_", car_cap[i])] <- 0
  # assign to the grid the values for visual check
  grid$car_cap <- class_extract[paste0("car_cap_", car_cap[i])]
  # vector to matrix
  car_cap_matrix <- matrix(grid$car_cap, 
                           nrow = 1, 
                           ncol = ncell(grid$patch.ID))
  
  # Write carrying capacity matrix to file
  write.matrix.nemo(car_cap_matrix, paste0(input_path, "/NEMO_patch_nbfem_a_spatial_", car_cap[i], ".txt"))
  
}

# visual check
plot(grid, col = grid$car_cap)


# ------- 5. create quanti_init_freq matrices -----

# matrix defining the species of the individual that is present in each patch (form the FSTAT file)
# in this model, the quanti_trait is defining the species (see manuscript and NEMO script) - with 2 allels/locus, -0.05 and +0.05
# 1 column, nrows = n patches of the simulation grid
# each patch is 0 if contains species 1, 1 if contains species 2

# function to write matrix in NEMO format (change such that NA result in emtpy rows)
write.matrix.nemo <- function(mat, outfile) {
  rows <- dim(mat)[1]
  cols <- dim(mat)[2]
  cat("{", file = outfile)
  
  for(i in 1:rows) {
    cat("{", file = outfile, append = TRUE)
    cat(mat[i,], sep = ",", file = outfile, append = TRUE)
    cat("}\n", file = outfile, append = TRUE)
  }
  
  cat("}\n", file = outfile, append = TRUE)
}


build_quanti_init_matrix <- function(
    grid_file,
    individuals_FSTAT_file,
    output_file
) {
  library(terra)
  library(data.table)
  
  grid <- terra::vect(grid_file)
  indiv <- fread(individuals_FSTAT_file)
  
  # species prefix and assign 0/1 codes based on first appearance
  indiv$prefix <- substr(indiv$SampleID, 1, 3)
  species_levels <- unique(indiv$prefix)
  
  if (length(species_levels) != 2) {
    stop("Expected exactly 2 species. Found: ", paste(species_levels, collapse = ", "))
  }
  
  species_codes <- setNames(c(0, 1), species_levels)
  indiv$species_code <- species_codes[indiv$prefix]
  
  # iniitalize matrix with 0s
  quanti_init_freq <- matrix(0, ncol = 1, nrow = terra::ncell(grid))
  
  # assign values: only first species per patch taken
  patch_first <- indiv[!duplicated(patch_id), .(patch_id, species_code)]
  quanti_init_freq[patch_first$patch_id, 1] <- patch_first$species_code
  
  write.matrix.nemo(quanti_init_freq, output_file)
  
  message("Quanti init matrix written to: ", output_file)
  
  return(quanti_init_freq)
}

build_quanti_init_matrix(
  grid_file = "Waldi_grid_4m_5396patches.shp",
  individuals_FSTAT_file = "wrapper_test_step2.txt", ## FSTAT file
  output_file = "NEMO_quanti_init_wrapper_test.txt"
)


# ------- 6. create dispersal matrices ------

# function to build both reduced dispersal matrices and connectivity matrix
# num_patch = number of patch of the simulation grid
# distance_matrix = pairwise distance matrix between the patches
# d = mean dispersal distance for the kernel
# b = shape of the kernel
# d_thresh = minimum dispersal rate under which dispersal rate is set to 0


## define your dispresal kernel (here exponential power from Chybicki et al. 2017)

kernel_function <- function(d, b, x) {
  #  scale parameter a
  a <- d * exp(lgamma(2 / b) - lgamma(3 / b))
  
  #  kernel function K_jk
  return(b * exp(-(x^b / a^b))  / (2* pi * a^2 * exp(lgamma(2/b))))
}


create_dispersal_matrices <- function(
    grid_file,
    dispersal_kernel,
    d,
    b,
    threshold_distance = NULL,
    output_prefix = "NEMO",
    output_path = ".",
    resolution = 4,
    plot_kernel = FALSE,
    plot_check_patch = NULL
) {
  library(terra)
  library(ggplot2)
  

  # simulation grid and calculate distance matrix
  grid <- terra::vect(grid_file)
  centroids <- centroids(grid)
  patch_coords <- as.data.frame(crds(centroids))
  patch_coords$patch_ID <- grid$patch.ID
  distance_matrix <- as.matrix(distance(centroids))
  num_patch <- terra::ncell(grid)
  
  # threshold
  if (is.null(threshold_distance)) {
    threshold <- d
  } else {
    threshold <- dispersal_kernel(d, b, threshold_distance)
  }
  
  # build matrix
  dispersal_matrix <- matrix(0, nrow = num_patch, ncol = num_patch)
  rate_matrix <- matrix(0, num_patch, num_patch)
  connectivity_matrix <- matrix(NA, num_patch, num_patch)
  
  for (i in 1:num_patch) {
    for (j in 1:num_patch) {
      rate_matrix[i, j] <- dispersal_kernel(d, b, distance_matrix[i, j]) # dispresal probability from patch i to patch j
    }
  }
  
  dispersal_matrix <- rate_matrix # full dispersal matrix
  
  for (i in 1:num_patch) {
    # indices of the patches sorted by dispersal probability in descending order
    ord <- order(rate_matrix[i, ], decreasing = TRUE)
    # re-order the values in row i of reduced rate matrix, connectivity matrix
    rate_matrix[i, ] <- rate_matrix[i, ord]
    connectivity_matrix[i, ] <- ord
    # remove elements below the threshold by setting them equal to NA
    to_remove <- which(rate_matrix[i, ] < threshold)
    rate_matrix[i, to_remove] <- NA
    connectivity_matrix[i, to_remove] <- NA
    # normalize dispersal rates to sum to 1
    sum_d <- sum(rate_matrix[i, ], na.rm = TRUE)
    if (sum_d > 0) {
      rate_matrix[i, !is.na(rate_matrix[i, ])] <- rate_matrix[i, !is.na(rate_matrix[i, ])] / sum_d
    }
  }
  
  # modified function with connectivity matrix starting always with focal patch, and rate matrix starting with 1 when the row is empty
  # to ensure no empty rows
  
  write.red.dispersal.matrix <- function(connectivity_matrix, rate_matrix, conn_file, rate_file) {
    cat("{", file = conn_file)                       # write connectivity matrix
    for (i in 1:nrow(connectivity_matrix)) {        
      cat("{", file = conn_file, append = TRUE)
      cat(i, file = conn_file, append = TRUE)         # start with the focal patch number
      if (!all(is.na(connectivity_matrix[i, ]))) {
        cat(",", file = conn_file, append = TRUE)
        non_na <- connectivity_matrix[i, !is.na(connectivity_matrix[i, ])]   # write the non-NA elements of the row, excluding the focal patch number itself
        non_na <- non_na[non_na != i]
        cat(non_na, sep = ",", file = conn_file, append = TRUE)
      }
      cat("}\n", file = conn_file, append = TRUE)
    }
    cat("}\n", file = conn_file, append = TRUE)
    
    cat("{", file = rate_file)                  # write rate matrix
    for (i in 1:nrow(rate_matrix)) {
      cat("{", file = rate_file, append = TRUE)
      if (all(is.na(rate_matrix[i, ]))) {
        cat(1, file = rate_file, append = TRUE)   # if the entire row is NA, write 1
      } else {
        cat(rate_matrix[i, which(!is.na(rate_matrix[i, ]))], sep = ",", file = rate_file, append = TRUE)   # write only the non-NA elements of the row
      }
      cat("}\n", file = rate_file, append = TRUE)
    }
    cat("}\n", file = rate_file, append = TRUE)
  }
  
  # output paths
  conn_file <- file.path(output_path, paste0(output_prefix, "_connectivity_matrix_d", d, ".txt"))
  rate_file <- file.path(output_path, paste0(output_prefix, "_rate_matrix_d", d, ".txt"))
  if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE) # create the output directory if it does not exists
  
  write.red.dispersal.matrix(connectivity_matrix, rate_matrix, conn_file, rate_file)
  
  # plot of the kernel (optional)
  if (plot_kernel) {
    x_vals <- 1:100
    kernel_vals <- dispersal_kernel(d, b, x_vals)
    plot_df <- data.frame(Distance = x_vals, Probability = kernel_vals)
    print(ggplot(plot_df, aes(x = Distance, y = Probability)) +
            geom_line(color = "steelblue", size = 1.2) +
            theme_minimal() +
            ggtitle("Dispersal Kernel") +
            ylab("Probability of Dispersal"))
  }
  
  # patch plot (optional)
  if (!is.null(plot_check_patch)) {
    message("Plotting dispersal for patch: ", plot_check_patch)
    grid_r <- terra::rast(ext = ext(grid), resolution = resolution, crs = crs(grid))
    n_cols <- dim(grid_r)[2]
    image_matrix <- matrix(dispersal_matrix[, plot_check_patch], ncol = n_cols, byrow = TRUE)
    par(mar = c(4, 4, 4, 10))
    plot(image_matrix, main = paste("Dispersal from Patch", plot_check_patch))
  }
  
  # outputs
  message("Files written:\n", conn_file, "\n", rate_file)
  return(list(
    dispersal_matrix = dispersal_matrix,
    connectivity_matrix = connectivity_matrix,
    rate_matrix = rate_matrix,
    distance_matrix = distance_matrix
  ))
}


### example

create_dispersal_matrices(
  grid_file = "Waldi_grid_4m_5396patches.shp",
  dispersal_kernel = kernel_function,
  d = 28.95,                                # average dispersal distance
  b = 0.53606,                              # shape of kernel 
  threshold_distance = 50,
  output_prefix = "Test_kernel",
  output_path = "output/",
  plot_kernel = TRUE, ## to fix - print the image
  plot_check_patch = 150 ## to fix - print the image
)
