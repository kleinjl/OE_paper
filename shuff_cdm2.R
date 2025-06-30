### Function to shuffle individuals at the genus or species level

# data = community data
# flora = dataframe with regional flora of 5 columns: tax_fam, tax_gen, tax_esp, gen_count (number of species per genus) and fam_count (number of genus per family)
# perc_cutoff = percentage which is to be shuffled (e.g., '0.2' means that either the 20% most abundant or rare (depending on the setting of 'shuffle_fraction') will be reshuffled)
# shuffle_fraction = "rare" = shuffle rare species (cumsum < cutoff); "common" = shuffle common species (cumsum > cutoff); "all = shuffle all species, regardless of commonness
# rarity_criterion = defines which aspect of rarity will be calculated. 'global_abd' = species are ranked according to their overall abundance in the dataset. 'local_abd' = species are ranked according to their mean local abundance-when-present. 'frequency' = species will be ranked according to their frequency of occurrence (i.e., number of sites in which a species is present)
# prop_tax_gen = proportion of individuals that are to be reshuffled at the genus level depth (value between 0 and 1)
# shuffle_aggregated = aggregate all individuals of a given taxon prior shuffling in order to attribute a single new random identification to all (e.g., if 13 Diospyros iturensis are selected to be shuffled in a given transect A, then all of those 13 individuals will be assigned the same new identification, for instance 'Pycnanthus angolensis'); if 'FALSE', then all individuals will be shuffled individually
# shuffle_esp_NA = select whether 'NA' can be assigned to the field 'tax_esp' (TRUE) or not (FALSE). if TRUE, some new random identifications will not be at the species level; if FALSE, then all new identifications will be at the species level


shuff.cdm2 <- function(data, flora, nperm, error, rarity_criterion, perc_cutoff, shuffle_fraction, prop_tax_gen, shuffle_family, shuffle_aggregated, shuffle_esp_NA) {
  
  # set up parallel run
  library(foreach)
  library(doParallel)
  numCores <- detectCores() - 3
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  
  #### NEW: Step 0 - calculate rarity and commonness based on cutoff value provided in function input ####
  
  #### (A.) Rare species only ####
  if (shuffle_fraction == "rare") { # in case we want to shuffle rare fraction
    
    #### (A.1) global abundance ####
    if (rarity_criterion == "global_abd") {
        
      rarecommon <- data %>% 
        filter(!is.na(tax_sp_level)) %>% 
        dplyr::count(tax_sp_level) %>% 
        arrange(desc(n)) %>% 
        mutate(shuffle_yesno = n <= quantile(n, perc_cutoff)) %>%
        rename(tax_tax = "tax_sp_level") %>% 
        dplyr::select(tax_tax, shuffle_yesno)
        
    }
    
    #### (A.2) mean local abundance-when-present ####
    else if (rarity_criterion == "local_abd") {
      
      ### step 1: get list of mean local abundance-when-present
      ## get community data
      cdm <- data %>% 
        filter(!is.na(tax_sp_level)) %>% 
        group_by(plot_name) %>%
        dplyr::count(tax_sp_level) %>%
        reshape::cast(.,
                      plot_name ~ tax_sp_level,
                      value='n') %>%
        mutate_all(~replace(., is.na(.), 0)) %>%
        column_to_rownames(., var = "plot_name")
      
      ## mean local abundance
      mean_local_abundance <- data.frame(tax_sp_level = colnames(cdm), mean_local_abd = rep(NA, ncol(cdm)))
      
      for (j in 1:ncol(cdm)) {
        # get number of plots in which species is present
        mean_local_abundance_sp <- cdm[j] %>% 
          rename(tax_sp_level = 1) %>% 
          filter(tax_sp_level > 0)
          
        mean_local_abundance[j,2] <- round(sum(mean_local_abundance_sp$tax_sp_level)/nrow(mean_local_abundance_sp),3)
        
      }
      
      rm(mean_local_abundance_sp)
      
      ## Step 2: calculate 'rarecommon'
      rarecommon <- mean_local_abundance %>% 
        arrange(desc(mean_local_abd)) %>% 
        mutate(shuffle_yesno = mean_local_abd <= quantile(mean_local_abd, perc_cutoff)) %>% 
        rename(tax_tax = "tax_sp_level") %>% 
        dplyr::select(tax_tax, shuffle_yesno)
      
    } 
    
    #### (A.3) Frequency of occurrence ####
    else if (rarity_criterion == "frequency") {
      
      ### step 1: get list of mean local abundance-when-present
      ## get community data
      cdm <- data %>% 
        filter(!is.na(tax_sp_level)) %>% 
        group_by(plot_name) %>%
        dplyr::count(tax_sp_level) %>%
        reshape::cast(.,
                      plot_name ~ tax_sp_level,
                      value='n') %>%
        mutate_all(~replace(., is.na(.), 0)) %>%
        column_to_rownames(., var = "plot_name")
      
      ### Step2: frequency = number of site occurrences per species
      rarecommon <- as.data.frame(apply(cdm, 2, function(x) sum(x != 0))) %>%
        rownames_to_column(var = "tax_tax") %>%
        rename(frequency = 2) %>% 
        arrange(desc(frequency)) %>% 
        mutate(shuffle_yesno = frequency <= quantile(frequency, perc_cutoff))
      
    }
  } # end 'shuffle_fraction == "rare" ' ifelse statement
  
  #### (B.) common species only ####
  else if (shuffle_fraction == "common") { # in case we want to shuffle common fraction
    
    #### (B.1) global abundance ####
    if (rarity_criterion == "global_abd") {
      
      rarecommon <- data %>% 
        filter(!is.na(tax_sp_level)) %>% 
        dplyr::count(tax_sp_level) %>% 
        arrange(desc(n)) %>% 
        mutate(shuffle_yesno = n >= quantile(n, (1-perc_cutoff))) %>%
        rename(tax_tax = "tax_sp_level") %>% 
        dplyr::select(tax_tax, shuffle_yesno)
    
    } 
    
    #### (B.2) mean local abundance-when-present ####
    else if (rarity_criterion == "local_abd") {
      
      ### step 1: get list of mean local abundance-when-present
      ## get community data
      cdm <- data %>% 
        filter(!is.na(tax_sp_level)) %>% 
        group_by(plot_name) %>%
        dplyr::count(tax_sp_level) %>%
        reshape::cast(.,
                      plot_name ~ tax_sp_level,
                      value='n') %>%
        mutate_all(~replace(., is.na(.), 0)) %>%
        column_to_rownames(., var = "plot_name")
      
      ## mean local abundance
      mean_local_abundance <- data.frame(tax_sp_level = colnames(cdm), mean_local_abd = rep(NA, ncol(cdm)))
      
      for (j in 1:ncol(cdm)) {
        # get number of plots in which species is present
        mean_local_abundance_sp <- cdm[j] %>% 
          rename(tax_sp_level = 1) %>% 
          filter(tax_sp_level > 0)
        
        mean_local_abundance[j,2] <- round(sum(mean_local_abundance_sp$tax_sp_level)/nrow(mean_local_abundance_sp),3)
      }
      
      rm(mean_local_abundance_sp)
      
      ## Step 2: calculate 'rarecommon'
      rarecommon <- mean_local_abundance %>% 
        arrange(desc(mean_local_abd)) %>% 
        mutate(shuffle_yesno = mean_local_abd >= quantile(mean_local_abd, (1-perc_cutoff))) %>% 
        rename(tax_tax = "tax_sp_level") %>% 
        dplyr::select(tax_tax, shuffle_yesno)
      
    } 
    
    #### (B.3) Frequency of occurrence ####
    else if (rarity_criterion == "frequency") {
      
      ### step 1: get list of mean local abundance-when-present
      ## get community data
      cdm <- data %>% 
        filter(!is.na(tax_sp_level)) %>% 
        group_by(plot_name) %>%
        dplyr::count(tax_sp_level) %>%
        reshape::cast(.,
                      plot_name ~ tax_sp_level,
                      value='n') %>%
        mutate_all(~replace(., is.na(.), 0)) %>%
        column_to_rownames(., var = "plot_name")
      
      ### Step2: frequency = number of site occurrences per species
      rarecommon <- as.data.frame(apply(cdm, 2, function(x) sum(x != 0))) %>%
        rownames_to_column(var = "tax_tax") %>%
        rename(frequency = 2) %>% 
        arrange(desc(frequency)) %>% 
        mutate(shuffle_yesno = frequency >= quantile(frequency, (1-perc_cutoff)))
      
    }
  } # end 'shuffle_fraction == "common" ' ifelse statement
  
  #### (C.) All species ####
  else if (shuffle_fraction == "all") {
    rarecommon <- data %>% 
      filter(!is.na(tax_sp_level)) %>% 
      dplyr::count(tax_sp_level) %>% 
      arrange(desc(n)) %>% 
      mutate(relabd = n/sum(n)) %>% 
      mutate(cumsum = cumsum(relabd)) %>% 
      mutate(shuffle_yesno = TRUE) %>% # all species can potentially be shuffled
      rename(tax_tax = "tax_sp_level") %>% 
      dplyr::select(tax_tax, shuffle_yesno)
  }
  
  
  
  ## Step 1: Get shuffeable individuals
  data1 <- data %>% 
    filter(plot_name %in% rownames(env)) %>% # filter only matches in 'env' dataset
    group_by(plot_name) %>% 
    mutate(new_transect_id = cur_group_id()) %>% # set ID with group_by
    ungroup() %>% 
    dplyr::select(plot_name, new_transect_id, id_n, tax_fam, tax_gen, tax_esp, tax_tax) %>% 
    left_join(., rarecommon, by = "tax_tax") %>%  # add commonness cutoff
    mutate(shuffle_yesno = ifelse(is.na(shuffle_yesno), TRUE, shuffle_yesno)) # set all non-identified individuals to be shuffled
  
  # count how many species per genus
  gen_count <- flora %>% 
    filter(tax_fam != "indet.") %>% 
    group_by(tax_gen) %>% 
    count()
  # count how many genera per family
  fam_count <- flora %>% 
    filter(tax_fam != "indet.") %>% 
    group_by(tax_fam) %>% 
    distinct(tax_gen, .keep_all = TRUE) %>% 
    count()
  
  #  get shuffleable individuals
  data1_shuffeable_tax_gen <- data1 %>% 
    arrange(plot_name, tax_fam, tax_gen) %>% 
    left_join(., fam_count, by = "tax_fam") %>%
    filter(n > 1 & shuffle_yesno == TRUE) # discard monogeneric families that cannot be shuffled
  data1_shuffeable_tax_esp <- data1 %>% 
    arrange(plot_name, tax_fam, tax_gen) %>% 
    left_join(., gen_count, by = "tax_gen") %>%
    filter(n > 1 & shuffle_yesno == TRUE) # discard monotypic genera that cannot be shuffled
  
  
  ## Step 2: build diagnostics dataframe [in output slot]
  # doing the math: Some individuals belong to a genus which contains only one species, which renders this genus "unshuffeable". The 'error' parameter of identifications to shuffle should be reflective of total individuals ('n_ind'), not shuffleable individuals ('n_ind_shuffeable'). Therefore, we calculate a new proportion ('new_error'). In the case that there are not enough shuffeable individuals in the sampling unit ('n_in_shuffeable') to attain the desired 'error' proportion, we can only set 'new_error' to 1 in order to shuffle all shuffeable individuals.
  diagnostics <- data.frame(matrix(NA,ncol = 14, nrow = n_distinct(data1$plot_name))) %>% 
    rename(plot_name = 1) %>% 
    rename(n_individuals_total = 2) %>% # total number of individuals per plot
    rename(n_individuals_shuffeable_tax_esp = 3) %>% # nb. of shuffeable individuals at the species level
    rename(n_individuals_shuffeable_tax_gen = 4) %>% # nb. of shuffeable individuals at the genus level
    rename(n_individuals_to_shuffle = 5) %>% # nb. of individuals that need to be shuffled to attain error percentage
    rename(n_individuals_to_shuffle_tax_esp = 6) %>% 
    rename(n_individuals_to_shuffle_tax_gen = 7) %>% 
    rename(new_error_tax_esp = 8) %>%  # new proportion of individuals that are to be shuffled in order to attain the error proportion given as input 'error'
    rename(new_error_tax_gen = 9) %>% # new proportion of individuals that are to be shuffled in order to attain the error proportion given as input 'error'
    rename(richness_change = 10) %>% 
    rename(bc_dissimilarity = 11) %>% 
    rename(cautious_error_initial = 12) %>% # initial fraction of individuals not identified at species level
    rename(cautious_error_new = 13) %>% # post-shuffling fraction of individuals not identified at species level
    rename(cautious_error_change = 14) # new cautious error relative to initial cautious error
  
  for (i in 1:length(unique(data1$new_transect_id))) {
    diagnostics[i,1] <- data1 %>% filter(new_transect_id == i) %>% distinct(plot_name)
    diagnostics[i,2] <- data1 %>% filter(new_transect_id == i) %>% count()
    diagnostics[i,3] <- data1_shuffeable_tax_esp %>% filter(new_transect_id == i) %>% count()
    diagnostics[i,4] <- data1_shuffeable_tax_gen %>% filter(new_transect_id == i) %>% count()
    diagnostics[i,5] <- as.integer(diagnostics[i,2] * error)
    diagnostics[i,6] <- as.integer(diagnostics[i,2] * (error*(1-prop_tax_gen)))
    diagnostics[i,7] <- as.integer(diagnostics[i,2] * (error*prop_tax_gen))
    diagnostics[i,12] <- as.integer(data1 %>% filter(new_transect_id == i & is.na(tax_esp)) %>% count())/diagnostics[i,2]
    # 'new_error_tax_esp' = percentage to subsample from 'individuals_to_shuffle_tax_esp'
    if (diagnostics[i,3]>=diagnostics[i,6]) {diagnostics[i,8] <- diagnostics[i,6]/diagnostics[i,3]}
    else {diagnostics[i,8] <- 1} # if nb. of individuals to be samples > nb. of shuffeable individuals, set 1
    # 'new_error_tax_gen' = percentage to subsample from 'individuals_to_shuffle_tax_gen'
    if (diagnostics[i,4]>=diagnostics[i,7]) {diagnostics[i,9] <- diagnostics[i,7]/diagnostics[i,4]} 
    else {diagnostics[i,9] <- 1} # if nb. of individuals to be samples > nb. of shuffeable individuals, set 1
    
  }
  
  ## [SOLVED] function ran into error in case of NaN
  diagnostics <- diagnostics %>% 
    mutate_all(~ ifelse(is.nan(.), 0, .))
  
  # build empty output slot lists
  # Initialize results lists
  mats <- vector("list", nperm)
  mats_diagnostics <- vector("list", nperm)
  
  
  if (sum(diagnostics[,3] < diagnostics[,6]) / nrow(diagnostics) > 0.5 & prop_tax_gen == 0) {
    
    results <- NULL } else if (sum(diagnostics[,4] < diagnostics[,7]) / nrow(diagnostics) > 0.5 & prop_tax_gen == 1) {
      
      results <- NULL } else {
        
        results <- foreach(k = 1:nperm, .packages = c("dplyr", "vegan", "tidyr", "tibble", "spatstat.geom")) %dopar% { # loop over nperm (nb of iterations)
          
          ## Step 3: Select random sample of shuffeable individuals with new (plot-specific) proportion ('new_error')
          
          lst1 <- list() # capture individuals to be shuffled at tax_esp level
          lst2 <- list() # capture individuals to be shuffled at tax_gen level
          lst3 <- list() # capture remaining individuals 'to_keep'
          
          # loop over 'plot_name'
          for (i in 1:length(unique(data1$new_transect_id))) {
            
            ## Individuals to be shuffled at 'tax_esp' level
            lst1[[i]] <- data1_shuffeable_tax_esp %>% 
              filter(new_transect_id == i) %>% 
              dplyr::slice_sample(prop = diagnostics[i,8]) %>% 
              dplyr::select(plot_name, tax_fam, tax_gen, tax_esp, id_n, tax_tax) %>% 
              rename(tax_esp_old = tax_esp) %>% 
              arrange(plot_name, tax_fam, tax_tax)
            
            ## Individuals to be shuffled at 'tax_gen' level
            lst2[[i]] <- data1_shuffeable_tax_gen %>% 
              filter(new_transect_id == i) %>% 
              dplyr::slice_sample(prop = diagnostics[i,9]) %>% 
              dplyr::select(plot_name, tax_fam, tax_gen, tax_esp, id_n, tax_tax) %>% 
              rename(tax_gen_old = tax_gen) %>% 
              rename(tax_esp_old = tax_esp) %>% 
              arrange(plot_name, tax_fam, tax_tax)
            
            ## Individuals to be kept
            lst3[[i]] <- data1 %>% 
              filter(new_transect_id == i) %>% 
              filter(!id_n %in% lst1[[i]]$id_n & !id_n %in% lst2[[i]]$id_n)
          }
          
          ## retrieve different cohorts
          # individuals to be shuffled at species level
          data2.shuff_tax_esp <- bind_rows(lst1)
          # individuals to be shuffled at genus level
          data2.shuff_tax_gen <- bind_rows(lst2)
          # individuals to be kept (later bind_rows with shuffled ids)
          data2.keep <- bind_rows(lst3)
          
          
          ## Step 4: Prepare dataframe for shuffling
          if (shuffle_aggregated) { # individuals grouped by taxa and shuffled together
            data3_tax_esp <- data2.shuff_tax_esp %>%
              group_by(plot_name, tax_gen) %>% ##### CRITICAL STEP: regroup by genus or tax_tax? ####
            mutate(n = n()) %>% 
              distinct(tax_gen, .keep_all = TRUE) %>% 
              mutate(tax_esp = NA) %>% 
              dplyr::select(-id_n, -tax_tax)
            
            data3_tax_gen <- data2.shuff_tax_gen %>%
              group_by(plot_name, tax_gen_old) %>% # resume by 'tax_gen'
              mutate(n = n()) %>% 
              distinct(tax_gen_old, .keep_all = TRUE) %>% 
              arrange(plot_name, tax_gen_old) %>% 
              mutate(tax_gen = NA) %>% 
              mutate(tax_esp = NA) %>% 
              dplyr::select(-id_n, -tax_tax)
          } else if (!shuffle_aggregated) { # individuals shuffled individually
            data3_tax_esp <- data2.shuff_tax_esp %>%
              mutate(tax_esp = NA)
            data3_tax_gen <- data2.shuff_tax_gen %>% 
              mutate(tax_gen = NA) %>% 
              mutate(tax_esp = NA)
          }
          
          
          ## Step 5: Assign new identifications in loop
          
          print(paste0("Shuffling ",k, " out of ", nperm," ."))
          
          #### (A) Species level ####
          
          if (nrow(data3_tax_esp)>1) { # only shuffle if dataframe is not empty
            for (i in 1:nrow(data3_tax_esp)) { # loop over taxa list
              
              # get all 'tax_esp' names for a given transect
              tax_esp_existing <- data1 %>% filter(plot_name == data3_tax_esp[i,]$plot_name) %>% 
                filter(tax_gen == data3_tax_esp[i,]$tax_gen) %>% 
                distinct(tax_esp) %>% 
                unlist(use.names = FALSE)
              
              # avoid replacing 'tax_esp' with original identification
              dset <- flora %>% filter(tax_gen == data3_tax_esp[i,]$tax_gen) %>% 
                dplyr::select(tax_esp) %>% 
                unlist(use.names = FALSE) %>%  # convert dataset to vector in pipe
                {if (shuffle_esp_NA) c(.,NA) else .} %>%  # add 'NA' to dset vector to leave 'tax_esp' == NA
                setdiff(., tax_esp_existing)
              
              ## [solved 17/09/2024] 'All-flora-possibilties-problem' when a transect contains all possibilities from the 'flora' dataframe (e.g., Santiria ebo, obovata & trimera) at the same time, then dset will turn 'zero' for we want to allocate a taxon that hasn't yet been allocated. for that rare case we'll implement an if else statement below to allocate whatever taxon, regardless of whether it already exists in the transect:
              # fill 'tax_esp' with randomly allocated species
              if (!is.empty(dset)){
                data3_tax_esp$tax_esp[i] <- base::sample(dset,1) # allocate 1 random sample
              } else data3_tax_esp$tax_esp[i] <- base::sample(tax_esp_existing,1)
              
            } # end loop for data3_tax_esp
          } # and if statement for data3_tax_esp
          
          
          #### (B) Genus level ####
          if (nrow(data3_tax_gen)>1){ # only shuffle if dataframe is not empty
            for (i in 1:nrow(data3_tax_gen)) { # loop over taxa list
              
              # get all 'tax_gen' names for a given transect
              tax_gen_existing <- data1 %>% filter(plot_name == data3_tax_gen[i,]$plot_name) %>% 
                filter(tax_fam == data3_tax_gen[i,]$tax_fam) %>% 
                distinct(tax_gen) %>% 
                unlist(use.names = FALSE)
              
              # avoid replacing 'tax_esp' with original identification
              dset <- flora %>% filter(tax_fam == data3_tax_gen[i,]$tax_fam) %>% 
                dplyr::select(tax_gen) %>% 
                distinct(tax_gen) %>% 
                unlist(use.names = FALSE) %>%  # convert dataset to vector in pipe
                setdiff(., tax_gen_existing)
              
              ## [solved]: 'All-flora-possibilties-problem' (see explanation above)
              # fill 'tax_gen' with randomly allocated genus
              if (!is.empty(dset)){
                data3_tax_gen$tax_gen[i] <- base::sample(dset,1) # allocate 1 random sample
              } else data3_tax_gen$tax_gen[i] <- base::sample(tax_gen_existing,1)
              
              
              # Continue to attribute 'tax_esp'. no need for 'setdiff' because the genus already is unique
              dset2 <- flora %>% filter(tax_gen == data3_tax_gen[i,]$tax_gen) %>% 
                dplyr::select(tax_esp) %>% 
                unlist(use.names = FALSE) %>%  # convert dataset to vector in pipe
                { if(shuffle_esp_NA) c(., NA) else .}
              
              # fill 'tax_esp' with randomly allocated species
              data3_tax_gen$tax_esp[i] <- base::sample(dset2,1) # allocate 1 random sample
              
            } # end loop for data3_tax_gen
            
          } # end if statement for data3_tax_gen
          
          
          
          ## Step 6: join individuals which have been shuffled with those that haven't.
          
          ## [solved] there will be duplicates because the shuffling function will attribute species that already      exist in the kept inviduals; then, reshape::cast() fails subsequently. check     here:https://stackoverflow.com/questions/1660124/how-to-sum-a-variable-by-group
          
          if (shuffle_aggregated) {
            
            if (nrow(data3_tax_esp)>1 & nrow(data3_tax_gen)>1){
              data4 <- bind_rows(uncount(data3_tax_esp, n), uncount(data3_tax_gen, n), data2.keep) 
            } else if (nrow(data3_tax_esp)>1 & nrow(data3_tax_gen)==0){
              data4 <- bind_rows(uncount(data3_tax_esp, n), data2.keep) 
            } else if (nrow(data3_tax_esp)==0 & nrow(data3_tax_gen)>1){
              data4 <- bind_rows(uncount(data3_tax_gen, n), data2.keep) 
            }
            
            data4 <- data4 %>% 
              unite(tax_tax, tax_gen, tax_esp, na.rm = TRUE, remove = FALSE, sep = " ") %>%
              dplyr::select(plot_name, tax_fam, tax_gen, tax_esp, tax_tax) %>% 
              arrange(plot_name, tax_fam, tax_gen, tax_esp) %>% 
              group_by(plot_name) %>% 
              mutate(new_transect_id = cur_group_id()) %>%  # set ID with group_by
              ungroup()
            
            
          } else if (!shuffle_aggregated) {
            
            if (nrow(data3_tax_esp)>1 & nrow(data3_tax_gen)>1){
              data4 <- bind_rows(data3_tax_esp, data3_tax_gen, data2.keep) 
            } else if (nrow(data3_tax_esp)>1 & nrow(data3_tax_gen)==0){
              data4 <- bind_rows(data3_tax_esp, data2.keep) 
            } else if (nrow(data3_tax_esp)==0 & nrow(data3_tax_gen)>1){
              data4 <- bind_rows(data3_tax_gen, data2.keep) 
            }
            
            data4 <- data4 %>% 
              unite(tax_tax, tax_gen, tax_esp, na.rm = TRUE, remove = FALSE, sep = " ") %>%
              dplyr::select(plot_name, tax_fam, tax_gen, tax_esp, tax_tax) %>% 
              arrange(plot_name, tax_gen, tax_esp) %>% 
              group_by(plot_name) %>% 
              mutate(new_transect_id = cur_group_id()) %>%  # set ID with group_by
              ungroup()
          }
          
          
          
          ## Step 7: calculate diagnostics: final error, a-diversity change & cautious error
          
          diss <- list()
          
          ## Relative richness
          for (i in 1:length(unique(data1$new_transect_id))) {
            # richness of original data
            rich_original <- data1 %>% 
              filter(new_transect_id == i) %>% 
              summarise(n_distinct(tax_tax)) %>% 
              as.integer()
            # richness of data after shuffling
            # rich_shuffled <- length(which(colSums(data5[i,]!=0) == nrow(data5[i,])))
            rich_shuffled <- data4 %>% 
              filter(new_transect_id == i) %>% 
              summarise(n_distinct(tax_tax)) %>% 
              as.integer()
            # relative increase/decrease of richness relative to original data in %
            diagnostics[i,10] <- ((rich_shuffled/rich_original)-1)*100
            
            ## calculate cautious_error_new
            diagnostics[i,13] <- as.integer(data4 %>% filter(new_transect_id == i & is.na(tax_esp)) %>% count())/diagnostics[i,2]
            diagnostics[i,14] <- diagnostics[i,13]-diagnostics[i,12] # change in cautious error sensu Morrisson (2020)
            
            ## final error rate: floristic dissimilarity
            original <- data1 %>%
              group_by(plot_name, tax_tax) %>% 
              count() %>%
              group_by(plot_name) %>% 
              mutate(new_transect_id = cur_group_id()) %>% 
              ungroup() %>% 
              filter(new_transect_id == i) %>% 
              mutate(type = "original")
            shuffled <- data4 %>% filter(new_transect_id == i) %>% 
              group_by(tax_tax) %>% 
              count() %>% 
              mutate(type = "shuffled")
            cdm_comp <- bind_rows(original,shuffled) %>% 
              reshape::cast(.,
                            type ~ tax_tax, value='n') %>%
              mutate_all(~replace(., is.na(.), 0)) %>%
              remove_rownames() %>%
              column_to_rownames(., var = "type")
            
            diagnostics[i,11] <- as.numeric(vegdist(cdm_comp, method = "bray", binary = FALSE))
            
          } # close diagnostics loop
          
          
          # return result
          list(data4 = data4, diagnostics = diagnostics)  # Return results
          
        } # close nperm loop
        
      } #### NEW #### close else statement
  
  
  # Stop the cluster
  stopCluster(cl)
  
  
  # Extract results
  if (!is.null(results)) {
    for (k in 1:nperm) {
      mats[[k]] <- results[[k]]$data4
      mats_diagnostics[[k]] <- results[[k]]$diagnostics %>% ## reduce 'diagnostics' for function output
        dplyr::select(plot_name, richness_change, cautious_error_change, bc_dissimilarity)
    }
    names(mats) <- paste0("Perm_",seq(1:nperm))
    names(mats_diagnostics) <- paste0("Perm_",seq(1:nperm))
  } else {
    mats <- NULL
    mats_diagnostics <- NULL
  }
  
  
  # return(list(shuffled_dissimilarity, transect))
  return(list(data = mats, diagnostics = mats_diagnostics))
  
} # end function
