#
# Copyright (c) 2023 The Broad Institute, Inc. All rights reserved.
#

library(khroma) # Default Color Palettes --- colorblind safe! --- courtesy of Paul Tol (https://personal.sron.nl/~pault/)
if (compareVersion(as.character(packageVersion('khroma')), "1.11.0") < 0) {
  warning(paste0("The version of khroma loaded ", as.character(packageVersion('khroma')),
                 " is below version 1.11.0, and will cause errors if `continuous.return_function` is set to TRUE."))
}

### Set Color Palette (for an annotation datatable)

set_annot_colors <- function( annot_table, 
                              # General Parameters ====
                              continuous.return_function = FALSE, # whether to return a function for continuous variables-- or a discrete vector with low/mid/high
                              continuous_columns = NULL, # vector of columns (either numeric or string)
                              autodetect_continuous=TRUE, # allow built-in detection of continuous columns
                              autodetect_continuous_nfactor_cutoff = 10, # number of non-NA factor levels allowed, before an annot is considered continuous
                              
                              # Special-Treatment Annotation Values ====
                              normal_annot_vals = c("^nat$", "^wt$", "^unmut$","^normal$","^0$"), # possible "normal" annotations (not case sensitive). Assigned the first color from the paired palette.
                              na_annot_vals = c("^na$", "^n.a.$", "^n/a$", "^unknown$", "^$"), # possible "NA" annotations (not case sensitive). Assigned na_color in discrete palettes, and placed at the end of vector.
                              # Default Color Palettes --- Paul Tol Color Palattes-- colorblind safe! (https://personal.sron.nl/~pault/) ====
                              qual.pals = list("Bright" = c('#4477AA', '#66CCEE', '#27B13E', '#CCBB44', '#EE6677', '#AA3377'), # green (originally '#228833') was modified to be more distinct from blue under tritanopia colorblindness
                                               "Vibrant" = c('#0077BB', '#33BBEE', '#009988', '#EE7733', '#CC3311', '#EE3377'),
                                               "Muted" = c( '#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77', '#CC6677', '#882255', '#AA4499'),
                                               "Light" = c( '#77AADD', '#99DDFF', '#44BB99', '#BBCC33', '#AAAA00', '#EEDD88', '#EE8866', '#FFAABB'),
                                               "HighContrast" = c('#004488', '#DDAA33', '#BB5566')),
                              # if adding a new paired-palette, please add LIGHTER colors FIRST, to ensure consistent formatting
                              pair.pals = list("BrightPaired" = c("#A8DBFF", '#4477AA', "#CAFFFF", '#66CCEE', "#8BFFA2", '#27B13E', "#FFFFA8", '#CCBB44', "#FFCADB", '#EE6677', "#FF97DB", '#AA3377'),
                                               "VibrantPaired" = c("#64DBFF", '#0077BB', "#64FDEC", '#009988', "#FFDB97", '#EE7733', "#FF97DB", '#EE3377' ),
                                               "HighContrastPaired" = c('#6699CC', '#004488', '#EECC66', '#997700', '#EE99AA', '#994455')),
                              seq.pals = list("YlOrBr" = c('#FFFFE5', '#FFF7BC', '#FEE391', '#FEC44F', '#FB9A29', '#EC7014', '#CC4C02', '#993404', '#662506'),
                                              "Iridescent" = c('#FEFBE9', '#FCF7D5', '#F5F3C1', '#EAF0B5', '#DDECBF', '#D0E7CA', '#C2E3D2', '#B5DDD8', '#A8D8DC', '#9BD2E1', '#8DCBE4', '#81C4E7', '#7BBCE7', '#7EB2E4', '#88A5DD', '#9398D2', '#9B8AC4', '#9D7DB2', '#9A709E', '#906388', '#805770', '#684957', '#46353A'),
                                              "Incandescent" = c('#CEFFFF', '#C6F7D6', '#A2F49B', '#BBE453', '#D5CE04', '#E7B503', '#F19903', '#F6790B', '#F94902', '#E40515', '#A80003')),
                              na_color = "#BBBBBB",
                              cont.pals = NULL, # by default, uses RColorBrewer Sequential Palettes, formatted list("palette_name" = c(low='#HEXVAL', mid='#HEXVAL', high='#HEXVAL'))
                              # ====
                              warn_for_interpolation = TRUE # warn if discrete function linearly interpolates colors (which may not be color blind safe) 
                              ) {
  annot_names = names(annot_table)
  annot_table_discrete = annot_table # start by assuming all annotations are discrete
  annot_table_continuous = annot_table[,0] # and that NO annotations are continuous
  
  
  ### Continuous / Discrete Separation
  
  if ( !autodetect_continuous && is.null(continuous_columns) )
    warning("No method for detecting continuous-annotations was selected. All annotations will be assigned discrete color-palettes.")
  
  # Manual Assignment from user
  if (!is.null(continuous_columns)) { # if the user has specified which columns are continuous
    if (is.character(continuous_columns)) { # if we specified column NAMES
      # make sure names are valid
      continuous_columns_valid = intersect(continuous_columns, annot_names) # get column names found in annot_table
      continuous_columns_invalid = setdiff(continuous_columns, annot_names) # get column names NOT found in annot_table
      if ( length(continuous_columns_invalid)>0 )
        warning(paste0("The following 'continuous' column-names were not found in the annotation table, and will be ignored:\n",paste(continuous_columns_invalid, collapse=", ")))
      # and convert names to indices
      continuous_columns = which(continuous_columns_valid %in% names(annot_table_discrete)) # subset continuous_columns to JUST valid column names, and convert to indices
    }
    
    annot_table_continuous = cbind(annot_table_continuous, annot_table_discrete[continuous_columns]) # add continuous columns to annot_table_continuous
    annot_table_discrete = annot_table_discrete[setdiff(1:length(annot_table_discrete), continuous_columns)] # remove continuous columns from annot_table_discrete
  }
  
  # Autodetection
  if (autodetect_continuous) {
    # auto-detect continuous variables that are STILL in annot_table_discrete
    
    # utility function to determine if an annotation column contains continuous data (i.e. )
    is.continuous <- function ( annot_col, na_annot_vals, nfactor_cutoff=10 ) {
      if (is.factor(annot_col)) {
        annot_vals = levels(annot_col) # take values from levels, do not re-sort
      } else annot_vals = sort(unique(annot_col)) # find and sort unique annotations (sorting prevents different color assignments based on order-of-appearance)
      annot_vals[is.na( annot_vals )] <- "NA" # set literal NA to a string "NA"
      
      n_na = sum(na.rm = TRUE, # sum number of NAs detected:
                 unlist(sapply(na_annot_vals, # for each possible type of na_annot_vals
                               function(na_vals) { grepl(na_vals, annot_vals, ignore.case=TRUE) }))) # check if it appears in our annotation-values
      
      all_numbers_regex <- "^(-?[0-9]*)((\\.?[0-9]+[eE]?[-\\+]?[0-9]+)|(\\.[0-9]+))*$" # regex to match all possible numeric strings, including scientific notation
      if ( (length(annot_vals)-n_na) > nfactor_cutoff && # if we have more than nfactor_cutoff unique (non-na) annotation-values
           sum(!grepl(all_numbers_regex, annot_vals), na.rm=TRUE) <= n_na ) # AND have all numeric data
        return(TRUE) # assume continuous
      return(FALSE) # otherwise, assume discrete
    }
    
    continuous_columns = which(unname(sapply( annot_table_discrete, simplify=TRUE, # get column indices
                                              function(annot_col) { is.continuous(annot_col, na_annot_vals, autodetect_continuous_nfactor_cutoff) }))) # of continuous columns
    
    annot_table_continuous = cbind(annot_table_continuous, annot_table_discrete[continuous_columns]) # add continuous columns to annot_table_continuous
    annot_table_discrete = annot_table_discrete[setdiff(1:length(annot_table_discrete), continuous_columns)] # remove continuous columns from annot_table_discrete
  }
  
  # Sanity Checks for continuous/discrete assignments
  if ( length(setdiff( annot_names, # sanity check-- if we have column names that are MISSING from:
                       union(names(annot_table_continuous), names(annot_table_discrete)) )) ) # either the discrete or continuous datatables
    stop(paste0("The following columns were assigned to NEITHER the continuous NOR the discrete datatable. Something has gone horribly wrong:\n", # complain!
                paste(setdiff( annot_names,  union(names(annot_table_continuous), names(annot_table_discrete)) ), collapse=", "))) # print the offending column names
  if ( length(intersect(names(annot_table_continuous), names(annot_table_discrete)))>0 ) # sanity check-- if we have columns that were assigned to BOTH the continuous and discrete datatables
    stop(paste0("The following columns were identified as both continuous and discrete. Something has gone horribly wrong:\n", # complain!
                paste(intersect(names(annot_table_continuous), names(annot_table_discrete)), collapse=", "))) # print the offending column names
  
  ### Color Assignment
  return(c(set_annot_colors_discrete(annot_table_discrete, # set discrete colors
                                     normal_annot_vals=normal_annot_vals, na_annot_vals=na_annot_vals,
                                     qual.pals = qual.pals, pair.pals = pair.pals, seq.pals = seq.pals, na_color = na_color,
                                     warn_for_interpolation=warn_for_interpolation),
           set_annot_colors_continuous(annot_table_continuous, palettes = cont.pals, na_color = na_color, return_function = continuous.return_function)))

}

set_annot_colors_discrete <- function( annot_table, # ====
                                       # Default Color Palettes --- Paul Tol Color Palattes-- colorblind safe! (https://personal.sron.nl/~pault/) ====
                                       qual.pals = list("Bright" = c('#4477AA', '#66CCEE', '#27B13E', '#CCBB44', '#EE6677', '#AA3377'), # green (originally '#228833') was modified to be more distinct from blue under tritanopia colorblindness
                                                        "Vibrant" = c('#0077BB', '#33BBEE', '#009988', '#EE7733', '#CC3311', '#EE3377'),
                                                        "Muted" = c( '#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77', '#CC6677', '#882255', '#AA4499'),
                                                        "Light" = c( '#77AADD', '#99DDFF', '#44BB99', '#BBCC33', '#AAAA00', '#EEDD88', '#EE8866', '#FFAABB'),
                                                        "HighContrast" = c('#004488', '#DDAA33', '#BB5566')),
                                       # if adding a new paired-palette, please add LIGHTER colors FIRST, to ensure consistent formatting
                                       pair.pals = list("BrightPaired" = c("#A8DBFF", '#4477AA', "#CAFFFF", '#66CCEE', "#8BFFA2", '#27B13E', "#FFFFA8", '#CCBB44', "#FFCADB", '#EE6677', "#FF97DB", '#AA3377'),
                                                        "VibrantPaired" = c("#64DBFF", '#0077BB', "#64FDEC", '#009988', "#FFDB97", '#EE7733', "#FF97DB", '#EE3377' ),
                                                        "HighContrastPaired" = c('#6699CC', '#004488', '#EECC66', '#997700', '#EE99AA', '#994455')),
                                       shuffle_paired_palette = TRUE, # should the paired colors be shuffled?
                                       seq.pals = list("YlOrBr" = c('#FFFFE5', '#FFF7BC', '#FEE391', '#FEC44F', '#FB9A29', '#EC7014', '#CC4C02', '#993404', '#662506'),
                                                       "Iridescent" = c('#FEFBE9', '#FCF7D5', '#F5F3C1', '#EAF0B5', '#DDECBF', '#D0E7CA', '#C2E3D2', '#B5DDD8', '#A8D8DC', '#9BD2E1', '#8DCBE4', '#81C4E7', '#7BBCE7', '#7EB2E4', '#88A5DD', '#9398D2', '#9B8AC4', '#9D7DB2', '#9A709E', '#906388', '#805770', '#684957', '#46353A'),
                                                       "Incandescent" = c('#CEFFFF', '#C6F7D6', '#A2F49B', '#BBE453', '#D5CE04', '#E7B503', '#F19903', '#F6790B', '#F94902', '#E40515', '#A80003')),
                                       
                                       # Special-Treatment Annotation Values ====
                                       normal_annot_vals = c("^nat$", "^wt$", "^unmut$","^normal$","^0$"), # possible "normal" annotations (not case sensitive). Assigned the first color from the paired palette.
                                       na_annot_vals = c("^na$", "^n.a.$", "^n/a$", "^unknown$", "^$"), # possible "NA" annotations (not case sensitive). Assigned na_color (below), and placed at the end of vector.
                                       na_color = "#BBBBBB", # default color for NA vals
                                       # ====
                                       warn_for_interpolation = TRUE
){
  
  # if we can't find the eval_with_seed() utility function
  if (!existsFunction("eval_with_seed")) { # define it here
    eval_with_seed <- function(expression, seed_string="randomstring") { 
      if (exists(".Random.seed")) {
        old_seed = .Random.seed # pull current seed
        on.exit({.Random.seed <<- old_seed}) # reset seed back to normal on exit
      } else {
        on.exit({set.seed(Sys.time())}) # randomize seed on exit
      }
      
      set.seed( sum(as.numeric(charToRaw(seed_string))) ) # set seed, based on string passed in
      eval( expression ) # match function
    }
  }
  
  
  annot_names = colnames(annot_table)
  
  # shuffle paired palettes
  pair.colors <- unlist(pair.pals) # unlist palates
  if (shuffle_paired_palette) {
    pair.colors <- eval_with_seed( lapply(sample( 1:(length(pair.colors)/2) *2)-1, # shuffle
                                          function(i) { c(pair.colors[i], pair.colors[i+1]) }), # organize pairs into a list
                                   seed_string = "paultolcolorblindsafe") # evaluate shuffling with seed for consistency
  }
  
  pair.index <- 0 # indexed at 0 for modulo operator
  seq.index <- 0 # indexed at 0 for modulo operator
  qual.index <- 0 # indexed at 0 for modulo operator
  annots_color_list <- list()
  for ( annot in annot_names  ) {
    if (is.factor( annot_table[[annot]]) ) {
      annot_vals = levels( annot_table[[annot]] ) # take values from levels, do not re-sort
    } else annot_vals = sort(unique( annot_table[[annot]] )) # find and sort unique annotations (sorting prevents different color assignments based on order-of-appearance)
    
    # locate NA values
    annot_vals[is.na( annot_vals )] <- "NA" # set NA to a string "NA"
    na_index = unlist(sapply(na_annot_vals, function(annot_table) { grep(annot_table, annot_vals, ignore.case=TRUE) }, USE.NAMES = FALSE)) # identify any NA indices
    n_colors = length(annot_vals)-length(na_index) # do not choose colors for NA values
    
    # Assign paired or sequential/qualitative colors, depending on vector length
    if( n_colors == 2 ) { # if we have a pair
      pair.count=pair.index+1 # re-index starting at 1 for.... everything other than the modulo operator
      annot_colors <- pair.colors[[pair.count]] # choose colors
      pair.index <- ( pair.index + 1 ) %% length(pair.colors) # increment pair.index, or restart if we've exhausted the pairs palette list
      
      # Assign the "normal" annotations the lighter color
      norm_index = unlist(sapply(normal_annot_vals, function(annot_table) { grep(annot_table, annot_vals, ignore.case=TRUE) }, USE.NAMES = FALSE)) # identify index of annot_vals that is "normal"
      if ( length(norm_index)!=0 ) { # if any of the "normal" annots can be found
        annot_vals <- c( annot_vals[norm_index], annot_vals[-norm_index] ) # move norm value to beginning-- this will assign it the lighter color.
        norm_index=numeric(0) # reset norm_index
      }
    } else { # pick a sequential or qualitative colorscheme
      if (sum(!grepl("^[0-9]+$", annot_vals), na.rm=TRUE) <= length(na_index)) { # if we have SEQUENTIAL data (i.e. its JUST digits)
        seq.count=seq.index+1 # re-index starting at 1 for.... everything other than the modulo operator
        annot_colors = grDevices::colorRampPalette(seq.pals[[seq.count]])(n_colors) # choose color palette
        seq.index <- (seq.index+1) %% length(seq.pals) # increment seq.index, or restart if we've exhausted the sequential palette list
        
      } else { # if we have QUALITATIVE data
        qual.count=qual.index+1 # re-index starting at 1 for.... everything other than the modulo operator
        
        if ( n_colors > length(qual.pals[[qual.count]]) ) { # if we have more annotations than unique colors in the chosen palette
          # warn user that we are interpolating colors, if toggle is on
          if (warn_for_interpolation)
            warning(paste0("The annotation ", annot, " has more unique values (", n_colors,") than this palette has unique colors (",
                           length(qual.pals[[qual.count]]), ").  Colors will be linearly interpolated, and may not be color-blind safe."))
          
          # interpolate colors from qual.palette
          qual.count=qual.index+1 # re-index starting at 1 for.... everything other than the modulo operator
          qual.colors = qual.pals[[qual.count]] # choose color palette
          annot_colors <- grDevices::colorRampPalette(qual.colors)(n_colors) # interpolate colors
          
          qual.index <- (qual.index+1) %% length(qual.pals) # increment qual.index, or restart if we've exhausted the qualitative palette list
          
        } else {
          qual.colors = qual.pals[[qual.count]] # choose color palette
          # assign colors, shuffled pseudo-randomly (using annot name as seed, to make sure we always choose the "same" random colors)
          annot_colors <- eval_with_seed( { sample(qual.colors, n_colors) } , seed_string = annot)
          qual.index <- (qual.index+1) %% length(qual.pals) # increment qual.index, or restart if we've exhausted the qualitative palette list
        }
        
      }
    }
    
    # Assign "NA" annotations grey
    if ( length(na_index)!=0 ) { # if any of the "normal" annots can be found
      annot_vals <- c( annot_vals[-na_index], annot_vals[na_index]) # move na values to end. these will assigned grey
      annot_colors <- c( annot_colors, rep(na_color, length(na_index)) ) # set color to grey
      na_index=numeric(0) # reset na_index to an empty vector
    }
    
    # Sanity Checks // Save into annots_color_list
    if ( length(annot_vals)!=length(annot_colors) ) { # do we have the right number of colors, total
      stop(paste0("For the annotation '", annot, "', the color-vector (length ", length(annot_colors),
                  ") does not match the number of annotations (",length(annot_vals),"). Something has gone terribly wrong."))
    } else {
      if ( length( setdiff(unique(annot_colors), na_color) ) != n_colors ) # are all our (non-na) colors unique?
        warning(paste0("For the annotation '", annot, "', there are ", n_colors," unique, non-NA annotations values.",
                       " However, the color-vector contains only ", length( setdiff(unique(annot_colors), na_color) )," unique non-NA colors."))
      
      annots_color_list[[annot]]$is_discrete <- TRUE # record that this palette is DISCRETE, not continuous
      annots_color_list[[annot]]$vals <- unname(annot_vals) # save annot_vals into list
      annots_color_list[[annot]]$colors <- unname(annot_colors) # save annot_colors into list
    }
    
  }
  
  return( annots_color_list )
}

set_annot_colors_continuous <- function( annot_table, # ====
                                         return_function = FALSE, # toggle for returning a FUNCTION rather than a discrete low-mid-high palette
                                         palettes = NULL, # when NULL, defaults to RColorBrewer sequential palettes
                                         shuffle_palettes = TRUE, # whether or not to shuffle palettes, compared to RColorBrewer default order
                                         na_color = "#BBBBBB" # color to be assigned to NA vals
                                         # ====
){
  # if user hasn't set a palette
  if (is.null(palettes)) { # use the RColorBrewer (non-function) khroma palettes (function)
    
    if (return_function) { # if return_function, palettes = a list of functions (khroma)
      # Khroma
      require(khroma)
      palettes = sapply(dplyr::filter(khroma::info(), type=="sequential")$palette, # use khroma sequential palettes
                        khroma::color) # 
    } else { # otherwise, palettes = a list of three colors, named low mid and high
      # R Color Brewer
      require(RColorBrewer)
      pal.names = setdiff(rownames(RColorBrewer::brewer.pal.info[which(RColorBrewer::brewer.pal.info$category=='seq'),]), # pull sequential palettes from RColorBrewer
                          c("Greys")) # EXCEPT the grey color-palette
      palettes = sapply(pal.names, simplify = FALSE,
                        function(pal) { # for each palette
                          colors = RColorBrewer::brewer.pal(3, pal) # pull out n_colors colors
                          names(colors) = c('low', 'mid', 'high') # label them
                          colors[['na_color']]=na_color # set NA color
                          return(colors)
                        })
    }
  }
  
  # shuffle palette order, if toggle is on
  if (shuffle_palettes) {
    
    # if we can't find the eval_with_seed() utility function
    if (!existsFunction("eval_with_seed")) { # define it here
      eval_with_seed <- function(expression, seed_string="randomstring") { 
        if (exists(".Random.seed")) {
          old_seed = .Random.seed # pull current seed
          on.exit({.Random.seed <<- old_seed}) # reset seed back to normal on exit
        } else {
          on.exit({set.seed(Sys.time())}) # randomize seed on exit
        }
        
        set.seed( sum(as.numeric(charToRaw(seed_string))) ) # set seed, based on string passed in
        eval( expression ) # match function
      }
    }
    
    palettes <- eval_with_seed( sample(palettes), # shuffle palettes
                                seed_string = "shuffling pseudorandomly") # evaluate shuffling with seed for consistency
  }
  
  annot_names = colnames(annot_table)
  
  pal.index <- 0 # indexed at 0 for modulo operator
  annots_color_list <- list() # initialize empty list
  for ( annot in annot_names  ) {
    annots_color_list[[annot]]$is_discrete <- FALSE # record that this palette is CONTINUOUS, not discrete
    if (return_function) {
      annots_color_list[[annot]]$vals <- NULL
      pal <- unname(palettes[[pal.index+1]]) # palette in khroma::color() structure
      # convert to add circlize::colorRamp2
      annots_color_list[[annot]]$colors <- circlize::colorRamp2(seq(min(as.numeric(annot_table[[annot]]), na.rm = T), # breaks vector from min
                                                                    max(as.numeric(annot_table[[annot]]), na.rm = T), # to max
                                                                    length.out = attr(pal, "max")), # with max-colors number of elements
                                                                pal(attr(pal, "max"))) # color-function
    } else {
      annots_color_list[[annot]]$vals <- names(palettes[[pal.index+1]]) # save names of color palette (i.e. low mid high)
      annots_color_list[[annot]]$colors <- unname(palettes[[pal.index+1]]) # choose color palette, re-indexing at 1, rather than at 0 
    }
    
    pal.index = (pal.index+1) %% length(palettes) # increment pal.index, or restart if we've exhausted the palette list\
  }
  
  return( annots_color_list )
  
}








### Misc Color Utility Functions (NOT required for set_annot_colors())


color_mod <- function(color_string, modifier, mod_R, mod_G, mod_B) {
  if ( !grepl("^#?[A-Fa-f0-9]{6}$", color_string) ) { # ensure that hex code is valid format
    stop(paste0("Invalid hex code '", color_string, "'. Please input a valid hex code."))
  }
  color_string = gsub("^#", "", color_string) # remove leading '#', if applicable
  
  if(missing(modifier)){ #if we're also missing a general color modifier
    if(missing(mod_R) & missing(mod_G) & missing(mod_R)){ #if we're missing ALL color modifiers
      modifier=30 #modifier is set to 30 by default
    } else{
      modifier=0 #if we have SOME colors specified, unspecified colors are not modified
    }}
  
  # If we're missing a color modification, it will be set to the default color modification level
  if(missing(mod_R)){
    mod_R = modifier}
  if(missing(mod_G)){
    mod_G = modifier}
  if(missing(mod_B)){
    mod_B = modifier}
  
  modifier_vec = c(mod_R,mod_G,mod_B)
  
  color_string_vector <- strtoi(c(substring(color_string, 1,2),  # red hex value
                                  substring(color_string, 3,4),  # green hex value
                                  substring(color_string, 5,6)), # blue hex value
                                16) # convert from hex to decimal
  color_string_modded = paste("#", # add '#' back in
                              paste(collapse = "", # combine RGB into one string
                                    sprintf("%02X", # enforce 2-digit hex codes
                                            as.hexmode( # convert from decimal to hex
                                              pmax( pmin(color_string_vector + modifier_vec, #increases the color by their modifier
                                                         c(255, 255, 255)), #cap colors at 255
                                                    c(0,0,0))))),  #keep colors above 0
                              sep = "")
  
  return(color_string_modded)
}  #modifies a color string (formatted "#XXXXXX" in) by the specified (or default) modifier amount


color_dist <- function(color_start, color_end) { 
  require(stringr)
  # Split colors into their R, G, and B HEX digits
  a<-stringr::str_remove(color_start, "#")
  color_start<-c(str_sub(a, 1, 2),str_sub(a, 3, 4), str_sub(a, 5, 6))
  
  a<-stringr::str_remove(color_end, "#")
  color_end<-c(str_sub(a, 1, 2),str_sub(a, 3, 4), str_sub(a, 5, 6))
  
  dist = strtoi(color_end, base = 16) - strtoi(color_start, base = 16)
  dist
} #take two input colors and output the R, G, and B distances


color_range <- function(base_color_start, base_color_end, number_of_colors) {
  
  # base_color_start="#042069"
  # base_color_end= "#B00B69"
  # number_of_colors=5
  
  distance = color_dist(base_color_start, base_color_end)
  
  color_scale = sapply( (1:number_of_colors)-1 ,
                        function (i, base_color_start, distance, number_of_colors) {
                          color_mod(base_color_start,
                                    mod_R = round(distance[1]* i/(number_of_colors-1),0),
                                    mod_G = round(distance[2]* i/(number_of_colors-1),0),
                                    mod_B = round(distance[3]* i/(number_of_colors-1),0))
                        },
                        base_color_start, distance, number_of_colors)
  
  return(color_scale)
} # takes two input colors and, using color_dist() and color_mod(), interpolates a range of number_of_colors colors




### TESTING:

# ### Testing Function
# 
# base_color_start="#042069"
# base_color_end= "#B00B69"
# number_of_colors=5
# 
# 
# ### Testing color-range
# color_scale = color_range(base_color_start, base_color_end, number_of_colors)
# 
# 
# 
# ### Testing color-scaling
# increment = 150
# scale = round(sort(unique(c(log(1/(1:ceiling(number_of_colors/2))), #forwards colors
#                             log(1:ceiling(number_of_colors/2))))) #backwards colors
#               * increment, 0) #scale and round
# 
# # Color Scale (dark to light)
# color_scale = sapply(1:number_of_colors,
#                      function(i, base_color_start) {color_mod(base_color_start, scale[i])},
#                      base_color_start)
# 
# # Color Scale (mix with red)
# color_scale = sapply(1:number_of_colors,
#                      function(i, base_color_start) {color_mod(base_color_start, mod_R = scale[i])},
#                      base_color_start)
# 
# # Color Scale (mix with green)
# color_scale = sapply(1:number_of_colors,
#                      function(i, base_color_start) {color_mod(base_color_start, mod_G = scale[i])},
#                      base_color_start)
# 
# # Color Scale (mix with blue)
# color_scale = sapply(1:number_of_colors,
#                      function(i, base_color_start) {color_mod(base_color_start, mod_B = scale[i])},
#                      base_color_start)
# 
# 
# 
# 
# 
# ### Testing w/ Plots
# 
# # test dataset
# require(dplyr)
# a <- dplyr::tibble(lab=as.factor(rep(1:length(scale),length(scale))),
#                    x = ceiling((1:(length(scale)*length(scale)))/length(scale)),
#                    y=x^2+as.numeric(lab))
# 
# # test function (to see colors)
# require(ggplot2)
# ggplot(a, aes(x, y, color = lab, fill = lab)) +
#   geom_line() +
#   # geom_bar(stat = "identity") +
#   scale_color_manual(values = color_scale) +
#   scale_fill_manual(values = color_scale)










### Code Graveyard

# # color palettes I was considering using for set_annot_colors_continuous(), before switching to RColorBrewer palettes
# # Default Color Palettes --- Paul Tol Color Palattes-- colorblind safe! (https://personal.sron.nl/~pault/) ====
# intp.pals = list("Bright"      = c('#4477AA', '#66CCEE', '#27B13E', '#CCBB44', '#EE6677', '#AA3377', # green (originally '#228833') was modified to be more distinct from blue under tritanopia colorblindness
#                                    mid="#FFFFFF", na_color="#BBBBBB"),
#                  "Vibrant"     = c('#0077BB', '#33BBEE', '#009988', '#EE7733', '#CC3311', '#EE3377',
#                                    mid="#FFFFFF", na_color="#BBBBBB"),
#                  "Muted"       = c( '#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77', '#CC6677', '#882255', '#AA4499',
#                                     na_color="#DDDDDD"),
#                  "Light"       = c( '#77AADD', '#99DDFF', '#44BB99', '#BBCC33', '#AAAA00', '#EEDD88', '#EE8866', '#FFAABB',
#                                     na_color="#DDDDDD"),
#                  "HighContrast" = c('#004488', '#DDAA33', '#BB5566',
#                                     na_color="#FFFFFF")),
# div.pals = list("PRGn"      = list(low="#762A83",mid="#F7F7F7",high="#1B7837", na_color="#FFEE99"),
#                 "BuRd"      = list(low="#2166AC",mid="#F7F7F7",high="#B2182B", na_color="#FFEE99"),
#                 "Nightfall" = list(low="#125A56",mid="#ECEADA",high="#A01813", na_color="#FFFFFF"),
#                 "Sunset"    = list(low="#364B9A",mid="#EAECCC",high="#A50026", na_color="#FFFFFF")),
# # ====
