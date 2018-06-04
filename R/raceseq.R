#' Process data
#'
#' @param dataset - input raw dataset from race-seq analysis to process
#'
#' @return processed dataset
#' @export
#'
#' @examples
process_data <- function(dataset) 
  {
  processed_dataset <- dataset 
  processed_dataset$tailed <- processed_dataset$tail_length > 0  #mark tailed reads
  processed_dataset$mapped <- processed_dataset$mapping_position != -1  #mark mapped reads
  processed_dataset_mapped <- processed_dataset %>% filter(mapped != 0)  #discard unmapped reads
  message("marked tailed and mapped reads")

  processed_dataset_mapped$ref_name_R5 <- as.character(processed_dataset_mapped$ref_name_R5)
  processed_dataset_mapped$ref_name_R3 <- as.character(processed_dataset_mapped$ref_name_R3)
  # tails_data_mapped_same_ref<-tails_data_mapped[tails_data_mapped$ref_name_R5==tails_data_mapped$ref_name_R3,]
  
  # mark uridylated reads
  processed_dataset_mapped$uridylated2 <- FALSE
  processed_dataset_mapped[processed_dataset_mapped$Utail_length > 0, ]$uridylated2 = TRUE
  message("marked uridylated reads")

   # in further analyses use only those read which got CTGAC delimiter identified in
  # the clipped fragment
  # for reporter analyses all mapped reads got CTGAC_R5 variable = 1 (because of short reads) so they will be included in the analysis
  processed_dataset_mapped_true <- processed_dataset_mapped[processed_dataset_mapped$CTGAC_R5 > 0, ]
  processed_dataset_mapped_true$ref_name = processed_dataset_mapped_true$ref_name_R5  #use ref_name_R5 as ref_name
  
  
  processed_dataset_mapped_true$tail_type = as.character(processed_dataset_mapped_true$tail_type)
  
  processed_dataset_mapped_true[processed_dataset_mapped_true$tail_type=='AU',]$uridylated2 = TRUE
  processed_dataset_mapped_true[processed_dataset_mapped_true$tail_type=='U_only',]$uridylated2 = TRUE

  # remove heterogenous tails from analysis
  processed_dataset_mapped_true_no_hetero = processed_dataset_mapped_true[-grep("hetero", processed_dataset_mapped_true$tail_type),
                                                            ]
  # remove other type tails (for which we can suspect they are not tails but rather origin from improper mapping/repeatmasker) from the analysis
  processed_dataset_mapped_true_no_hetero_no_other = processed_dataset_mapped_true_no_hetero[-grep("other",
                                                                                                   processed_dataset_mapped_true_no_hetero$tail_type), ]
  
  
  # treat all AG,UG or UA tails as other_no_tail
  
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$tail_type ==
                                              "AG", ]$tail_type <- "other_no_tail"
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$tail_type ==
                                              "UG", ]$tail_type <- "other_no_tail"
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$tail_type ==
                                              "UA", ]$tail_type <- "other_no_tail"
  processed_dataset_mapped_true_no_hetero_no_other <- processed_dataset_mapped_true_no_hetero_no_other[-grep("other",
                                                                                                             processed_dataset_mapped_true_no_hetero_no_other$tail_type), ] #remove all other from analysis
  
  message("removed unnecessary tail types")
  # create classes for A-tail lengths (0,1,2-5,6-10,11-20,21-30,30+)
  processed_dataset_mapped_true_no_hetero_no_other$A_length = ""
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$Atail_length ==
                                              0, ]$A_length = "0"
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$Atail_length ==
                                              1, ]$A_length = "1"
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$Atail_length %in%
                                              seq(2, 5, 1), ]$A_length = "2-5"
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$Atail_length %in%
                                              seq(6, 10, 1), ]$A_length = "6-10"
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$Atail_length %in%
                                              seq(11, 20, 1), ]$A_length = "11-20"
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$Atail_length %in%
                                              seq(21, 30, 1), ]$A_length = "21-30"
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$Atail_length >
                                              30, ]$A_length = "30+"
  processed_dataset_mapped_true_no_hetero_no_other$A_length <- factor(processed_dataset_mapped_true_no_hetero_no_other$A_length,
                                                               levels = c("0", "1", "2-5", "6-10", "11-20", "21-30", "30+"))
  
  message("processed A tails lengths")
  # create classes for U-tail lengths (0,1,2,3-5,6-10,10+)
  processed_dataset_mapped_true_no_hetero_no_other$U_length = ""
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$Utail_length ==
                                              0, ]$U_length = "0"
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$Utail_length ==
                                              1, ]$U_length = "1"
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$Utail_length ==
                                              2, ]$U_length = "2"
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$Utail_length %in%
                                              seq(3, 5, 1), ]$U_length = "3-5"
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$Utail_length %in%
                                              seq(6, 10, 1), ]$U_length = "6-10"
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$Utail_length >
                                              10, ]$U_length = "10+"
  processed_dataset_mapped_true_no_hetero_no_other$U_length <- factor(processed_dataset_mapped_true_no_hetero_no_other$U_length,
                                                               levels = c("10+", "6-10", "3-5", "2", "1", "0"))
  
  message("processed U tails lengths")
  
  # modify levels of tail_types to have U_only,A-only,AU or no_tail
  processed_dataset_mapped_true_no_hetero_no_other_tails <- processed_dataset_mapped_true_no_hetero_no_other
  processed_dataset_mapped_true_no_hetero_no_other_tails$tail_type <- as.character(processed_dataset_mapped_true_no_hetero_no_other_tails$tail_type)
  processed_dataset_mapped_true_no_hetero_no_other_tails$tail_type <- factor(processed_dataset_mapped_true_no_hetero_no_other_tails$tail_type,
                                                                      levels = c("U_only", "AU", "no_tail", "A_only"))
  
  message("reordered tail types factor")
  all_data <- processed_dataset_mapped_true_no_hetero_no_other_tails %>% dplyr::select(V1,tail_sequence,tail_type,Atail,Atail_length,A_length,U_length,Utail_length,tail_length,tail_source,transcript,cell_line,localization,condition,replicate,primer_name,project_name,mapping_position,exp_type,CTGAC_R5,terminal_nucleotides,uridylated,uridylated2,tailed,mapped) %>% dplyr::group_by(project_name,condition,replicate,transcript,primer_name,cell_line,uridylated)
  all_data <- all_data %>% dplyr::filter(tail_length <= 64)
  message("done")
  return(all_data)
}

#' Title
#'
#' @param input_dataset 
#'
#' @return
#' @export
#'
#' @examples
factorize_tail_lengths <- function(input_dataset) {
  output_dataset <- input_dataset
  output_dataset$A_length = ""
  output_dataset[output_dataset$Atail_length ==
                                                     0, ]$A_length = "0"
  output_dataset[output_dataset$Atail_length ==
                                                     1, ]$A_length = "1"
  output_dataset[output_dataset$Atail_length %in%
                                                     seq(2, 5, 1), ]$A_length = "2-5"
  output_dataset[output_dataset$Atail_length %in%
                                                     seq(6, 10, 1), ]$A_length = "6-10"
  output_dataset[output_dataset$Atail_length %in%
                                                     seq(11, 20, 1), ]$A_length = "11-20"
  output_dataset[output_dataset$Atail_length %in%
                                                     seq(21, 30, 1), ]$A_length = "21-30"
  output_dataset[output_dataset$Atail_length >
                                                     30, ]$A_length = "30+"
  output_dataset$A_length <- factor(output_dataset$A_length,
                                                                      levels = c("0", "1", "2-5", "6-10", "11-20", "21-30", "30+"))
  
  message("processed A tails lengths")
  # create classes for U-tail lengths (0,1,2,3-5,6-10,10+)
  output_dataset$U_length = ""
  output_dataset[output_dataset$Utail_length ==
                                                     0, ]$U_length = "0"
  output_dataset[output_dataset$Utail_length ==
                                                     1, ]$U_length = "1"
  output_dataset[output_dataset$Utail_length ==
                                                     2, ]$U_length = "2"
  output_dataset[output_dataset$Utail_length %in%
                                                     seq(3, 5, 1), ]$U_length = "3-5"
  output_dataset[output_dataset$Utail_length %in%
                                                     seq(6, 10, 1), ]$U_length = "6-10"
  output_dataset[output_dataset$Utail_length >
                                                     10, ]$U_length = "10+"
  output_dataset$U_length <- factor(output_dataset$U_length,
                                                                      levels = c("10+", "6-10", "3-5", "2", "1", "0"))
  
  
  message("processed U tail lengths")
  output_dataset$tail_length<-as.numeric(output_dataset$tail_length)
  #calculate maximum possible tail length in given dataset
  max_tail_length<-summary(output_dataset$tail_length)[6]
  output_dataset$tail_length_fac<-'NA'
  #bin tails based on their length
  if (max_tail_length>=60) {
    output_dataset[output_dataset$tail_length>=60,]$tail_length_fac<-'60+'
  }
  message("processed tails 60+")
  if (max_tail_length>=50) {
    for (temp_len in seq(50,59,1)) {
      if(any(output_dataset$tail_length==temp_len)) {
        output_dataset[output_dataset$tail_length==temp_len,]$tail_length_fac<-'50-59'
        
      }
    }
  }
  message("processed tails 50+")
  if (max_tail_length>=40) {
    for (temp_len in seq(40,49,1)) {
      if(any(output_dataset$tail_length==temp_len)) {
        output_dataset[output_dataset$tail_length==temp_len,]$tail_length_fac<-'40-49'
      }
    }
  }
  message("processed tails 40+")
  if (max_tail_length>=30) {
    for (temp_len in seq(30,39,1)) {
      if(any(output_dataset$tail_length==temp_len)) {
        output_dataset[output_dataset$tail_length==temp_len,]$tail_length_fac<-'30-39'
      }
    }
  }
  message("processed tails 30+")
  if (max_tail_length>=20) {
    for (temp_len in seq(20,29,1)) {
      if(any(output_dataset$tail_length==temp_len)) {
        output_dataset[output_dataset$tail_length==temp_len,]$tail_length_fac<-'20-29'
      }
    }
  }
  message("processed tails 20+")
  for (temp_len in seq(10,19,1)) {
    if(any(output_dataset$tail_length==temp_len)) {
      output_dataset[output_dataset$tail_length==temp_len,]$tail_length_fac<-'10-19'
    }
  }
  message("processed tails 10+")
  for (temp_len in seq(1,9,1)) {
    if(any(output_dataset$tail_length==temp_len)) {
      output_dataset[output_dataset$tail_length==temp_len,]$tail_length_fac<-'1-9'
    }
  }
  message("processed tails 1-9")
  output_dataset[output_dataset$tail_length==0,]$tail_length_fac<-'0'
  #reorder tail length factor
  message("processed all")
  #output_dataset$tail_length_fac<-as.factor(output_dataset$tail_length_fac)
  output_dataset$tail_length_fac<-factor(output_dataset$tail_length_fac,levels=c('0','1-9','10-19','20-29','30-39','40-49','50-59','60+'))
  return(output_dataset)
}


#' Filters RACE-seq data based on various criteria
#'
#' @param dataset                - input dataset
#' @param transcript2            - transcripts to include 
#' @param exp_type2              - experiment type to include ("OVR","KD","NT","LEAP",...)
#' @param mapping_position_min   - minimal mapping position of a read
#' @param mapping_position_max   - maximum mapping position of a read
#' @param project                - projects (sequencing runs) to include (vector)
#' @param conditions             - experimental conditions to include (vector)
#' @param primers                - RACE primers to include (vector)
#' @param cell_lines             - cell lines 
#' @param localizations          - localizations (from fractionation experiments) to include
#' @param persons                - persons preparing libraries
#'
#' @return                       - filtered dataset
#' @export
#'
#' @examples
filter_data <- function(dataset,transcript2,exp_type2,mapping_position_min=NA,mapping_position_max=NA,project=NA,conditions=NA,primers=NA,cell_lines=NA,localizations=NA,persons=NA) 
  {
  #get filtered data from input dataset
  #represent all input parameters as vectors
  transcript2 <- as.vector(transcript2)
  exp_type2 <- as.vector(exp_type2)
  project <- as.vector(project)
  cell_lines <- as.vector(cell_lines)
  primers <- as.vector(primers)
  conditions <- as.vector(conditions)
  localizations <- as.vector(localizations)
  persons <- as.vector(persons)
  #first, filter by transcript and experiment type (OVR,KD,LEAP)
  test_trans <- dataset %>% dplyr::filter(transcript %in% transcript2,exp_type %in% exp_type2) 
  #if project is specified - used for filtering
  if(length(project)>0 & !is.na(project))
  {
        test_trans <- test_trans %>% dplyr::filter(project_name %in% project)
  }
  #if conditions are specified - use for filtering
  if(length(conditions)>0 & !is.na(conditions)) {
    test_trans <- test_trans %>% dplyr::filter(condition %in% conditions)
  }
  #if primers are specified - use for filtering
  if(length(primers)>0 & !is.na(primers)) {
    test_trans <- test_trans %>% dplyr::filter(primer_name %in% primers)
  }
  #if cell_lines are specified - use for filtering
  if(length(cell_lines)>0 & !is.na(cell_lines)) {
    test_trans <- test_trans %>% dplyr::filter(cell_line %in% cell_lines)
  }
  #if localization are specified - use for filtering
  if(length(localizations)>0 & !is.na(localizations)) {
    test_trans <- test_trans %>% dplyr::filter(localization %in% localizations)
  }
  #if localization are specified - use for filtering
  if(length(persons)>0 & !is.na(persons)) {
    test_trans <- test_trans %>% dplyr::filter(peron %in% persons)
  }
  #if min and max mapping positions are specified - used for filtering
  if (!is.na(mapping_position_min)) {
    test_trans <- test_trans %>% dplyr::filter(mapping_position>mapping_position_min)
  }
  if (!is.na(mapping_position_max)) {
    test_trans <- test_trans %>% dplyr::filter(mapping_position<mapping_position_max)
  }
  #finally - drop unused factor levels
  test_trans <- test_trans %>% droplevels()
  return(test_trans)
}


#' Check homogenity of variances using Bartlett's test
#'
#' @param input_dataset - dataset with values to be compared
#' @param values        - string with the name of column containing numeric values used for comparison
#' @param grouping_var  - string with the name of column containing grouping values 
#'
#' @return logical TRUE if variances are homogenous between groups, otherwise FALSE
#' @export
#'
#' @examples
are_variances_homogenous <- function(input_dataset,values,grouping_var) 
{
  input_dataset[[grouping_var]] <- as.factor(input_dataset[[grouping_var]]) #make sure grouping var is factor
  bartlett <- bartlett.test(input_dataset[[values]] ~ input_dataset[[grouping_var]])
  if(bartlett$p.value>0.05) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Calculate stats 
#' Wrapper for Tukey and Dunn tests from PMCMRPlus package
#'
#' @param input_dataset - input data to analyze
#' @param test - test to run (c("DunnManyOne","Dunn","Tukey"))
#' @param group_by - split data by additional factor(s) and perform test in groups
#' @param values - name of the variable containing numerical values for comparison (or number of column)
#' @param grouping_var - name of the variable with grouping data (e.g. conditions) (or number of column)
#' @param p.adjust.method  - method of p.value adjustment (one of available in PMCMRPlus, default = Benjamini-Hochberg ("BH"))
#'
#' @return list of PMCMRPlus outputs
#' @export
#'
#' @examples
#' calculate_stats(datasets::iris,values="Sepal.Length",grouping_var="Species",test="Tukey")
calculate_stats <-
  function(input_dataset,
           test = 'DunnManyOne',
           group_by = NA,
           values = NA,
           grouping_var = NA,
           p.adjust.method = "BH")
  {
    #make sure that input dataset exists and that is not empty
    #assertthat::assert_that(exists(deparse(substitute(input_dataset))),msg = "Input dataset does not exists")
    assertthat::not_empty(input_dataset)
    assertthat::is.string(p.adjust.method)
    input_dataset[[grouping_var]] <-
      as.factor(input_dataset[[grouping_var]]) #make sure grouping var is factor
    #if additional grouping factors are not specified - treat input dataset as a single dataset
    if (is.na(group_by)) {
      #available tests:
      if (test == 'DunnManyOne')
      {
        #if DunnManyOne test is chosen, than all values are compared to the group specified as the first grouping factor level
        stats <-
          PMCMRplus::kwManyOneDunnTest(input_dataset[[values]], input_dataset[[grouping_var]], p.adjust.method = p.adjust.method)
      } else if (test == 'Dunn') {
        stats <-
          PMCMRplus::kwAllPairsDunnTest(input_dataset[[values]], input_dataset[[grouping_var]], p.adjust.method = p.adjust.method)
      } else if (test == 'Tukey') {
        #for Tukey's test it is required that variances are homogenous between groups (checked using  Bartlett's test)
        if (!are_variances_homogenous(input_dataset, values, grouping_var)) {
          warning("variances are not homogenous between groups (according to Bartlett's test)")
        }
        stats <-
          PMCMRplus::tukeyTest(input_dataset[[values]], input_dataset[[grouping_var]], p.adjust.method = p.adjust.method)
      }
      else {
        stop("Wrong test type provided as an argument")
      }
      
    } else
    {
      #if additional grouping factor(s) specified - calculate stats for each group separately
      input_dataset <-
        input_dataset %>% dplyr::group_by(.dots = c(group_by))
      #stats_multiple <- tibble()
      if (test == 'DunnManyOne')
      {
        stats_multiple = input_dataset %>% dplyr::do(stats = PMCMRplus::kwManyOneDunnTest(.[[values]], .[[grouping_var]]))
      } else if (test == 'Dunn')
      {
        stats_multiple = input_dataset %>% dplyr::do(stats = PMCMRplus::kwAllPairsDunnTest(.[[values]], .[[grouping_var]]))
      } else if (test == 'Tukey')
      {
        variances_equal <-
          input_dataset %>% dplyr::do(bartlett = are_variances_homogenous(., values, grouping_var)) %>% filter(bartlett == FALSE)
        if (nrow(variances_equal) > 0) {
          #TODO - make lsiting of all possible groups (when more than on grouping factor)
          warning(
            paste(
              "Variances are not homogenous between conditions for groups: ",
              paste(variances_equal[[1]], collapse = ", "),
              sep = ""
            )
          )
        }
        stats_multiple = input_dataset %>% dplyr::do(stats = PMCMRplus::tukeyTest(.[[values]], .[[grouping_var]]))
      } else
      {
        stop("Wrong test type provided as an argument")
      }
      #create names for each group
      groups_names = stats_multiple %>% dplyr::ungroup() %>% dplyr::mutate_if(is.character, as.factor) %>% dplyr::group_by_if(is.factor) %>% dplyr::group_vars()
      #add created names to groups
      stats_multiple  <-
        stats_multiple %>% dplyr::ungroup() %>% dplyr::mutate(name = paste(!!!rlang::syms(groups_names), sep =
                                                                             "___"))
      names(stats_multiple$stats) <- stats_multiple$name
      stats <- stats_multiple$stats
    }
    return(stats)
  }


#' Plot statistics for input dataset
#'
#' @param input_dataset    - dataset to calculate stats and plot
#' @param test             - statistical test to use
#' @param group_by         - additional grouping factor(s)
#' @param values           - string specifying column with values used for statistics
#' @param grouping_var     - string specifying column with grouping values (factor)
#' @param p.adjust.method  - method of p.value adjustment (one of available in PMCMRPlus, default = Benjamini-Hochberg ("BH"))
#' @param include_jitter   - include jitter plot (TRUE/FALSE), default = TRUE
#' @param include_errorbars - include errorbars (TRUE/FALSE) - default TRUE
#' @param errorbar_type    - type of errorbar (SD/SE/CI) - CI - confidence interval of a mean at 95% confidence, default = SD
#' @param include_statistics 
#'
#' @return list with 'stats' - calculated statistics and 'plot' - created plot
#' @export
#'
#' @examples
#' plot_stats(iris,values="Sepal.Length",grouping_var = "Species",test="Tukey",include_jitter = F,errorbar_type = "CI")

plot_stats <-
  function(input_dataset,
           test = 'DunnManyOne',
           group_by = NA,
           values = NA,
           grouping_var = NA,
           p.adjust.method = "BH",
           include_jitter = TRUE,
           include_errorbars = TRUE,
           errorbar_type = "SD",
           include_statistics = TRUE)
  {
    #assertthat::assert_that(exists(deparse(substitute(input_dataset))),msg = "Input dataset does not exists")
    assertthat::assert_that(not_empty(input_dataset))
    assertthat::is.string(p.adjust.method)
    output <- list()
    
    #calculate statistics
    if (include_statistics == TRUE) {
      stats <-
        calculate_stats(
          input_dataset = input_dataset,
          test = test,
          group_by = group_by,
          values = values,
          grouping_var = grouping_var,
          p.adjust.method = p.adjust.method
        )
    }
    #calculate data statistics - mean
    input_means <-
      input_dataset %>% dplyr::group_by(.dots = c(group_by, grouping_var)) %>% dplyr::summarise(
        n = n(),
        mean_val = mean(!!sym(values), na.rm = TRUE),
        sd_val = sd(!!sym(values), na.rm = TRUE)
      ) %>% dplyr::mutate(se_val = sd_val / sqrt(n),
                   CI_mean = qnorm(0.975) * sd_val / sqrt(n))
    #create plot:
    plot_out <-
      input_means %>% ggplot(aes_string(x = grouping_var)) + geom_bar(stat =
                                                                        "identity", position = "dodge", aes(y = mean_val))
    #add errorbars:
    if (include_errorbars == TRUE) {
      if (errorbar_type == 'SD') {
        plot_out <-
          plot_out  + geom_errorbar(
            aes(ymin =  mean_val - sd_val, ymax = mean_val + sd_val),
            colour = "black",
            width = 0.1,
            position = position_dodge(0.9)
          )
      } else if (errorbar_type == 'SE') {
        plot_out <-
          plot_out  + geom_errorbar(
            aes(ymin =  mean_val - se_val, ymax = mean_val + se_val),
            colour = "black",
            width = 0.1,
            position = position_dodge(0.9)
          )
      } else if (errorbar_type == 'CI') {
        plot_out <-
          plot_out  + geom_errorbar(
            aes(ymin =  mean_val - CI_mean, ymax = mean_val + CI_mean),
            colour = "black",
            width = 0.1,
            position = position_dodge(0.9)
          )
      } else {
        stop("Wrong errorbar type provided")
      }
    }
    if (include_jitter == TRUE) {
      plot_out <-
        plot_out + geom_jitter(data = input_dataset, aes_string(x = grouping_var, y =
                                                                  values))
    }
    
    
    #create facets
    if (!is.na(group_by)) {
      facet_by <- as.vector(group_by)
      if (!is.na(facet_by) & length(facet_by) == 1) {
        facet_by <- c(facet_by, ".")
      }
      plot_out <-
        plot_out + facet_grid (reformulate(facet_by[2:length(facet_by)], facet_by))
    }
    
    if (include_statistics == TRUE) {
      #get plot data
      pg <- ggplot_build(plot_out)
      #get height of errorbars (if available)
      if (include_errorbars == TRUE) {
        errorbar_data = pg$data[[2]] %>% dplyr::select(PANEL, group, ymax)
        errorbar_data2 = pg$data[[2]] %>% dplyr::select(PANEL2 = PANEL,
                                                 group2 = group,
                                                 ymax2 = ymax)
      } else {
        #if errorbars will not be shown - get  heights from the main plot
        errorbar_data = pg$data[[1]] %>% dplyr::select(PANEL, group, ymax)
        errorbar_data2 = pg$data[[1]] %>% dplyr::select(PANEL2 = PANEL,
                                                 group2 = group,
                                                 ymax2 = ymax)
      }
      
      #create_annotation
      if (!is.na(group_by))
      {
        #not finished, when the stats output is a list - TODO - DONE
        #melt p.value matrix from stats output and bind rows
        annotation <-
          lapply(lapply(stats, '[[', "p.value"), melt) %>% bind_rows(.id = "group_name")
        #select and rename columns
        annotation <-
          annotation %>% dplyr::select(
            group_name = group_name,
            start = Var2,
            condition = Var1,
            pvalue = value
          )
        if (length(group_by) > 1) {
          #if there was more grouping factors, get them from the group name (created by calculate_stats function)
          #groups are sepearated with ___
          facets <-
            lapply(strsplit(annotation$group_name, split = "___"), function(x) {
              as.data.frame(x)
            }) %>% dplyr::bind_cols() %>% t
          colnames(facets) <- group_by
          annotation <- cbind(facets, annotation)
        } else {
          colnames(annotation) <- c(group_by, "start", "condition", "pvalue")
        }
      } else
      {
        #melt p.value matrix of stats
        annotation <- stats$p.value %>% melt
        annotation$group_name = "NA"
        annotation <-
          annotation %>% dplyr::select(
            group_name = group_name,
            start = Var2,
            condition = Var1,
            pvalue = value
          )
      }
      output$stats <-
        annotation #store calculated and annotated stats in the output
      #create mapping of conditions to group numbers on the plot (group numbers are ordered same way as input factor)
      no_panels = 1 #default number of panels on the plot
      if (!is.na(group_by)) {
        #use tidyr::complete_ to fill missing panels
        #create 2 separate panel annotations - one for each group used for pairwise comparison
        panels_annot <-
          input_dataset %>% dplyr::group_by(.dots = c(group_by)) %>% tidyr::complete_(group_by) %>% dplyr::slice(1) %>% dplyr::select_(.dots = c(group_by))
        panels_annot$PANEL <- seq(1, nrow(panels_annot))
        panels_annot2 <-
          input_dataset %>% dplyr::group_by(.dots = c(group_by)) %>% tidyr::complete_(group_by) %>% dplyr::slice(1) %>% dplyr::select_(.dots = c(group_by))
        panels_annot2$PANEL2 <- seq(1, nrow(panels_annot))
        no_panels <- nrow(panels_annot) #update number of panels
      }
      #create 2 separate conditions annotations - one for each group used for pairwise comparison
      conditions_df <-
        data.frame(
          PANEL = rep(sequence(no_panels), each = length(levels(input_dataset[[grouping_var]]))),
          condition = rep(levels(input_dataset[[grouping_var]]), no_panels),
          group = rep(seq(1, length(
            levels(input_dataset[[grouping_var]])
          )), no_panels)
        )
      
      conditions_df2 <-
        data.frame(
          PANEL2 = rep(sequence(no_panels), each = length(levels(input_dataset[[grouping_var]]))),
          start = rep(levels(input_dataset[[grouping_var]]), no_panels),
          group2 = rep(seq(1, length(
            levels(input_dataset[[grouping_var]])
          )), no_panels)
        )
      if (!is.na(group_by)) {
        conditions_df <- conditions_df %>% dplyr::left_join(panels_annot)
        conditions_df2 <- conditions_df2 %>% dplyr::left_join(panels_annot2)
      }
      #create final annotation with all conditions using join
      annotation <-
        annotation %>% plyr::join(conditions_df) %>% plyr::join(errorbar_data) %>% plyr::join(conditions_df2) %>% plyr::join(errorbar_data2)
      
      #filter out not significant comparisons and comparisons which were not done (empty ymax or ymax2)
      annotation <-
        annotation %>% dplyr::filter(pvalue < 0.05, !is.na(ymax), !is.na(ymax2))
      if (nrow(annotation) > 0) {
        #if there were any significant pairwise comparisons, calculate bars positions and map pvalues to asterisks
        annotation <-
          annotation %>% dplyr::group_by(.dots = c(group_by)) %>% dplyr::mutate(group_no = row_number(), max_y =
                                                                    max(ymax, ymax2))
        annotation <-
          annotation %>% dplyr::mutate(y = max_y + group_no * (max_y * 0.15))
        
        annotation$sig = ''
        if (any(annotation$pvalue <= 0.05)) {
          annotation[annotation$pvalue <= 0.05, ]$sig <- "*"
        }
        if (any(annotation$pvalue <= 0.01)) {
          annotation[annotation$pvalue <= 0.01, ]$sig <- "**"
        }
        if (any(annotation$pvalue <= 0.001)) {
          annotation[annotation$pvalue <= 0.001, ]$sig <- "***"
        }
        if (any(annotation$pvalue <= 0.0001)) {
          annotation[annotation$pvalue <= 0.0001, ]$sig <- "****"
        }
        ylim_value = max(annotation$y) + 0.15 * max(annotation$y)
        #ylim_value = max(annotation$y)
        #expand y axis to fir significance bars
        plot_out = plot_out + expand_limits(y = c(0, ylim_value))
        plot_out <-
          plot_out + geom_signif(
            data = annotation,
            aes(
              xmin = start,
              xmax = condition,
              annotations = sig,
              y_position = y
            ),
            textsize = 10,
            vjust = 0.6,
            manual = TRUE
          )
      }
    }
    output$plot <- plot_out
    
    return(output)
  }

#' Analyze uridylation
#'
#'
#' @param dataset                - input dataset
#'
#' @param transcript2            - transcripts to include 
#' @param exp_type2              - experiment type to include ("OVR","KD","NT","LEAP",...)
#' @param mapping_position_min   - minimal mapping position of a read
#' @param mapping_position_max   - maximum mapping position of a read
#' @param project                - projects (sequencing runs) to include (vector)
#' @param conditions             - experimental conditions to include (vector)
#' @param primers                - RACE primers to include (vector)
#' @param cell_lines             - cell lines 
#' @param localizations          - localizations (from fractionation experiments) to include
#' @param persons                - persons preparing libraries
#' @param facet_by_spec          - specify faceting options
#' @param include_jitter         - should output plot include jitter points? (T/F)
#' @param include_statistics     - should statistics be calculated and placed on the plot? (T/F)
#'
#' @return                       - list with caluclated values, statistics and plot
#' @export
#'
#' @examples
analyze_uridylation <-
  function(dataset,
           transcript2,
           exp_type2,
           mapping_position_min = NA,
           mapping_position_max = NA,
           project = NA,
           facet_by_spec = NA,
           conditions = NA,
           include_jitter = TRUE,
           localizations = NA,
           persons = NA,
           primers = NA,
           cell_lines = NA,
           include_statistics = TRUE)
  {
    facet_by <- as.vector(facet_by_spec)
    if (!is.na(facet_by) & length(facet_by) == 1) {
      facet_by <- c(facet_by, ".")
    }
    output = list() #list for storing output
    message("Getting data")
    test_trans <-
      filter_data(
        dataset,
        transcript2 = transcript2,
        exp_type2 = exp_type2,
        conditions = conditions,
        project = project,
        mapping_position_max = mapping_position_max,
        mapping_position_min = mapping_position_min,
        cell_lines = cell_lines,
        primers = primers,
        persons = persons,
        localizations = localizations
      )
    #check for single or empty condition groups - exclude from analysis
    if (!is.na(facet_by)) {
      groups_to_exclude <-
        test_trans %>% group_by(.dots = c(facet_by_spec, "condition")) %>% count() %>% group_by_(.dots = c(facet_by_spec)) %>% count() %>% filter(nn %in% c(0, 1)) %>% select_(facet_by_spec)
      if (nrow(groups_to_exclude) > 0) {
        #print(groups_to_exclude)
        message("Excluding groups with zero or one condition")
        test_trans <- test_trans %>% anti_join(groups_to_exclude)
      }
    }
    #print(mapping_position_min)
    # print(mapping_position_max)
    #print(transcript2)
    #print(exp_type2)
    #get filtered data from input dataset
    #first, filter by transcript and experiment type (OVR,KD,LEAP)
    plot_title <- paste(exp_type2, transcript2, "uridylation")
    message("Calculating uridylation")
    uridylation  <-
      calculate_mean_uridylation(
        test_trans,
        transcript2 = transcript2,
        exp_type2 = exp_type2,
        conditions = conditions,
        project = project,
        mapping_position_max = mapping_position_max,
        mapping_position_min = mapping_position_min,
        cell_lines = cell_lines,
        primers = primers,
        persons = persons,
        localizations = localizations,
        facet_by_spec = facet_by_spec
      )
    #print(uridylation)
    output$calculated_values <- uridylation
    message("calculating statistics and creating plot")
    plot_out <-
      plot_stats(
        uridylation,
        test = 'Tukey',
        group_by = facet_by_spec,
        values = "freq_urid",
        grouping_var = "condition",
        p.adjust.method = "BH",
        include_statistics = include_statistics
      )
    output$plot <- plot_out
    return(output)
  }



#' Create plot with tails 
#'
#' @param dataset
#' @param transcript2            - transcripts to include 
#' @param exp_type2              - experiment type to include ("OVR","KD","NT","LEAP",...)
#' @param mapping_position_min   - minimal mapping position of a read
#' @param mapping_position_max   - maximum mapping position of a read
#' @param project                - projects (sequencing runs) to include (vector)
#' @param conditions             - experimental conditions to include (vector)
#' @param primers                - RACE primers to include (vector)
#' @param cell_lines             - cell lines 
#' @param localizations          - localizations (from fractionation experiments) to include
#' @param persons                - persons preparing libraries
#' @param facet_by_spec          - specify faceting options
#' @param include_jitter         - should output plot include jitter points? (T/F)
#' @param include_statistics     - should statistics be calculated and placed on the plot? (T/F)
#' @param tail_type_to_analyze   - tail type to include on plot (A_only/U_only/AU)
#'
#' @return                       - plot
#' @export
#'
#' @examples
plot_tails_fraction <-
  function(dataset,
           transcript2,
           exp_type2,
           mapping_position_min = NA,
           mapping_position_max = NA,
           project = NA,
           facet_by_spec = NA,
           conditions = NA,
           include_jitter = TRUE,
           localizations = NA,
           persons = NA,
           primers = NA,
           cell_lines = NA,
           include_statistics = TRUE,
           tail_type_to_analyze = "U_only")
  {
    facet_by <- as.vector(facet_by_spec)
    if (!is.na(facet_by) & length(facet_by) == 1) {
      facet_by <- c(facet_by, ".")
    }
    output = list() #list for storing output
    message("Getting data")
    test_trans <-
      filter_data(
        dataset,
        transcript2 = transcript2,
        exp_type2 = exp_type2,
        conditions = conditions,
        project = project,
        mapping_position_max = mapping_position_max,
        mapping_position_min = mapping_position_min,
        cell_lines = cell_lines,
        primers = primers,
        persons = persons,
        localizations = localizations
      )
    #check for single or empty condition groups - exclude from analysis
    if (!is.na(facet_by)) {
      groups_to_exclude <-
        test_trans %>% group_by(.dots = c(facet_by_spec, "condition")) %>% count() %>% group_by_(.dots = c(facet_by_spec)) %>% count() %>% filter(nn %in% c(0, 1)) %>% select_(facet_by_spec)
      if (nrow(groups_to_exclude) > 0) {
        #print(groups_to_exclude)
        message("Excluding groups with zero or one condition")
        test_trans <- test_trans %>% anti_join(groups_to_exclude)
      }
    }
    #print(mapping_position_min)
    # print(mapping_position_max)
    #print(transcript2)
    #print(exp_type2)
    #get filtered data from input dataset
    #first, filter by transcript and experiment type (OVR,KD,LEAP)
    plot_title <- paste(exp_type2, transcript2, "uridylation")
    message("Calculating fraction of tails")
    tails_fraction  <-
      calculate_tails_fraction(
        test_trans,
        transcript2 = transcript2,
        exp_type2 = exp_type2,
        conditions = conditions,
        project = project,
        mapping_position_max = mapping_position_max,
        mapping_position_min = mapping_position_min,
        cell_lines = cell_lines,
        primers = primers,
        persons = persons,
        localizations = localizations,
        facet_by_spec = facet_by_spec,
        tail_type_to_analyze = tail_type_to_analyze
      )
    #print(uridylation)
    output$calculated_values <- tails_fraction
    message("calculating statistics and creating plot")
    plot_out <-
      plot_stats(
        tails_fraction,
        test = 'Tukey',
        group_by = facet_by_spec,
        values = "freq_tail",
        grouping_var = "condition",
        p.adjust.method = "BH",
        include_statistics = include_statistics
      )
    output$plot <- plot_out
    return(output)
  }


#' Title
#'
#' @param dataset                - input dataset
#'
#' @param transcript2            - transcripts to include 
#' @param exp_type2              - experiment type to include ("OVR","KD","NT","LEAP",...)
#' @param mapping_position_min   - minimal mapping position of a read
#' @param mapping_position_max   - maximum mapping position of a read
#' @param project                - projects (sequencing runs) to include (vector)
#' @param conditions             - experimental conditions to include (vector)
#' @param primers                - RACE primers to include (vector)
#' @param cell_lines             - cell lines 
#' @param localizations          - localizations (from fractionation experiments) to include
#' @param persons                - persons preparing libraries
#' @param facet_by_spec          - specify faceting options
#'
#' @return                       - tibble with mean uridylation data
#' @export
#'
#' @examples
calculate_mean_uridylation <- function(dataset,transcript2,exp_type2,mapping_position_min=NA,mapping_position_max=NA,project=NA,facet_by_spec=NA,conditions=NA,localizations=NA,persons=NA,primers=NA,cell_lines=NA) 
{
  facet_by <- as.vector(facet_by_spec)
  if (!is.na(facet_by) & length(facet_by)==1) {
    facet_by<-c(facet_by,".")
  }
  output = list() #list for storing output
  test_trans <- filter_data(dataset,transcript2 = transcript2,exp_type2 = exp_type2,conditions = conditions,project = project,mapping_position_max = mapping_position_max,mapping_position_min = mapping_position_min,cell_lines = cell_lines,primers = primers,persons = persons, localizations = localizations)
  #print(mapping_position_min)
  # print(mapping_position_max)
  #print(transcript2)
  #print(exp_type2)
  #get filtered data from input dataset
  #first, filter by transcript and experiment type (OVR,KD,LEAP)
  plot_title <- paste(exp_type2,transcript2,"uridylation")
  #if min and max mapping positions are specified - used for filtering
  if (!is.na(mapping_position_min)) {
    plot_title <- paste(plot_title,"min map pos: ",mapping_position_min)
  }
  if (!is.na(mapping_position_max)) {
    plot_title <- paste(plot_title,"max map pos: ",mapping_position_max)
  }
  
  test_trans <- test_trans %>% dplyr::group_by(condition,replicate,uridylated2)
  #print(test_trans)
  #summarize uridylation data using dplyr
  if (is.na(facet_by)) {
    test_trans <- test_trans %>% dplyr::summarize(n_urid=n()) %>% ungroup() %>% dplyr::group_by(condition,replicate) %>% dplyr::mutate(freq_urid = n_urid/sum(n_urid)) %>% dplyr::group_by(condition,uridylated2) %>% dplyr::mutate(mean_freq_urid = mean(freq_urid), sd_urid = sd(freq_urid))
  }  else {
    test_trans <- test_trans %>% group_by(.dots = c("condition","replicate","project_name","uridylated2",facet_by[facet_by!="."]))
    #if using facets - group by projects also
    #test_trans <- test_trans %>% dplyr::summarize(n_urid=n()) %>% ungroup() %>% dplyr::group_by(condition,replicate,project_name) %>% dplyr::mutate(freq_urid = n_urid/sum(n_urid)) %>% dplyr::group_by(condition,uridylated2,project_name) %>% dplyr::mutate(mean_freq_urid = mean(freq_urid), sd_urid = sd(freq_urid))
    test_trans <- test_trans %>% dplyr::summarize(n_urid=n()) %>% ungroup() %>% dplyr::group_by(.dots = c("condition","replicate",facet_by[facet_by!="."])) %>% dplyr::mutate(freq_urid = n_urid/sum(n_urid)) %>% dplyr::group_by(.dots = c("condition","uridylated2",facet_by[facet_by!="."])) %>% dplyr::mutate(mean_freq_urid = mean(freq_urid), sd_urid = sd(freq_urid))
  }
  
  #leave only uridylation values (exclude uridylation == FALSE)
  test_trans <- test_trans %>% dplyr::filter(uridylated2==TRUE)

  return(test_trans)  
}


#' Title
#'
#' @param dataset                - input dataset
#'
#' @param transcript2            - transcripts to include 
#' @param exp_type2              - experiment type to include ("OVR","KD","NT","LEAP",...)
#' @param mapping_position_min   - minimal mapping position of a read
#' @param mapping_position_max   - maximum mapping position of a read
#' @param project                - projects (sequencing runs) to include (vector)
#' @param conditions             - experimental conditions to include (vector)
#' @param primers                - RACE primers to include (vector)
#' @param cell_lines             - cell lines 
#' @param localizations          - localizations (from fractionation experiments) to include
#' @param persons                - persons preparing libraries
#' @param facet_by_spec          - specify faceting options
#' @param tail_type_to_analyze   - tail type to include on plot (A_only/U_only/AU)
#'
#' @return
#' @export
#'
#' @examples
calculate_tails_fraction <- function(dataset,transcript2,exp_type2,mapping_position_min=NA,mapping_position_max=NA,project=NA,facet_by_spec=NA,conditions=NA,localizations=NA,persons=NA,primers=NA,cell_lines=NA,tail_type_to_analyze="U_only") 
{
  tail_type_to_analyze <- as.vector(tail_type_to_analyze)
   facet_by <- as.vector(facet_by_spec)
  if (!is.na(facet_by) & length(facet_by)==1) {
    facet_by<-c(facet_by,".")
  }
  output = list() #list for storing output
  test_trans <- filter_data(dataset,transcript2 = transcript2,exp_type2 = exp_type2,conditions = conditions,project = project,mapping_position_max = mapping_position_max,mapping_position_min = mapping_position_min,cell_lines = cell_lines,primers = primers,persons = persons, localizations = localizations)
  #print(mapping_position_min)
  # print(mapping_position_max)
  #print(transcript2)
  #print(exp_type2)
  #get filtered data from input dataset
  #first, filter by transcript and experiment type (OVR,KD,LEAP)
  plot_title <- paste(exp_type2,transcript2,"uridylation")
  #if min and max mapping positions are specified - used for filtering
  if (!is.na(mapping_position_min)) {
    plot_title <- paste(plot_title,"min map pos: ",mapping_position_min)
  }
  if (!is.na(mapping_position_max)) {
    plot_title <- paste(plot_title,"max map pos: ",mapping_position_max)
  }
  
  test_trans <- test_trans %>% dplyr::group_by(condition,replicate,tail_type)
  #print(test_trans)
  #summarize uridylation data using dplyr
  if (is.na(facet_by)) {
    test_trans <- test_trans %>% dplyr::summarize(n_tails=n()) %>% ungroup() %>% dplyr::group_by(condition,replicate) %>% dplyr::mutate(freq_tail = n_tails/sum(n_tails)) %>% dplyr::group_by(condition,tail_type) %>% dplyr::mutate(mean_freq_tail = mean(freq_tail), sd_urid = sd(freq_tail))
  }  else {
    test_trans <- test_trans %>% group_by(.dots = c("condition","replicate","project_name","tail_type",facet_by[facet_by!="."]))
    #if using facets - group by projects also
    #test_trans <- test_trans %>% dplyr::summarize(n_urid=n()) %>% ungroup() %>% dplyr::group_by(condition,replicate,project_name) %>% dplyr::mutate(freq_urid = n_urid/sum(n_urid)) %>% dplyr::group_by(condition,uridylated2,project_name) %>% dplyr::mutate(mean_freq_urid = mean(freq_urid), sd_urid = sd(freq_urid))
    test_trans <- test_trans %>% dplyr::summarize(n_tails=n()) %>% ungroup() %>% dplyr::group_by(.dots = c("condition","replicate",facet_by[facet_by!="."])) %>% dplyr::mutate(freq_tail = n_tails/sum(n_tails)) %>% dplyr::group_by(.dots = c("condition","tail_type",facet_by[facet_by!="."])) %>% dplyr::mutate(mean_freq_tail = mean(freq_tail), sd_urid = sd(freq_tail))
  }
  
  #leave only uridylation values (exclude uridylation == FALSE)
  test_trans <- test_trans %>% dplyr::filter(tail_type %in% tail_type_to_analyze)
  
  return(test_trans)  
}


#' Title
#'
#' @param dataset                - input dataset
#'
#' @param transcript2            - transcripts to include 
#' @param exp_type2              - experiment type to include ("OVR","KD","NT","LEAP",...)
#' @param mapping_position_min   - minimal mapping position of a read
#' @param mapping_position_max   - maximum mapping position of a read
#' @param project                - projects (sequencing runs) to include (vector)
#' @param conditions             - experimental conditions to include (vector)
#' @param primers                - RACE primers to include (vector)
#' @param cell_lines             - cell lines 
#' @param localizations          - localizations (from fractionation experiments) to include
#' @param persons                - persons preparing libraries
#' @param facet_by_spec          - specify faceting options
#'
#' @return                       - tibble with mean tail lengths
#' @export
#'
#' @examples
calculate_mean_tail_lengths <-
  function(dataset,
           transcript2,
           exp_type2,
           mapping_position_min = NA,
           mapping_position_max = NA,
           project = NA,
           facet_by_spec = NA,
           conditions = NA,
           localizations = NA,
           persons = NA,
           primers = NA,
           cell_lines = NA)
  {
    facet_by <- as.vector(facet_by_spec)
    if (!is.na(facet_by) & length(facet_by) == 1) {
      facet_by <- c(facet_by, ".")
    }
    output = list() #list for storing output
    test_trans <-
      filter_data(
        dataset,
        transcript2 = transcript2,
        exp_type2 = exp_type2,
        conditions = conditions,
        project = project,
        mapping_position_max = mapping_position_max,
        mapping_position_min = mapping_position_min,
        cell_lines = cell_lines,
        primers = primers,
        persons = persons,
        localizations = localizations
      )
    
    test_trans <- test_trans %>% dplyr::group_by(condition, replicate)
    
    #summarize tail_lengths data data using dplyr
    if (is.na(facet_by))
    {
      test_trans <-
        test_trans %>% dplyr::summarize(mean_tail_length_rep = mean(tail_length)) %>% ungroup() %>% dplyr::group_by(condition) %>% dplyr::mutate(
          mean_tail_length = mean(mean_tail_length_rep),
          sd_tail_length = sd(mean_tail_length_rep)
        )
    } else {
      test_trans <-
        test_trans %>% group_by(.dots = c("condition", "replicate", facet_by[facet_by !=
                                                                               "."]))
      ##if using facets - group by projects also
      test_trans <-
        test_trans %>% dplyr::summarize(mean_tail_length_rep = mean(tail_length)) %>% ungroup() %>% dplyr::group_by(.dots = c("condition", facet_by[facet_by !=
                                                                                                                                                      "."])) %>% dplyr::mutate(
                                                                                                                                                        mean_tail_length = mean(mean_tail_length_rep),
                                                                                                                                                        sd_tail_length = sd(mean_tail_length_rep)
                                                                                                                                                      )
    }
    
    #leave only uridylation values (exclude uridylation == FALSE)
    #test_trans <- test_trans %>% dplyr::filter(uridylated2==TRUE)
    
    return(test_trans)
  }





#' Title
#'
#' Plot mapping positions
#'
#' @param dataset                - input dataset 
#' @param exp_type2              - experiment type (OVR,KD,LEAP,NT)
#' @param conditions             - conditions to include in the analysis
#' @param mapping_position_min   - minimal mapping position to include on plot
#' @param mapping_position_max   - maximum mapping position to include on plot
#' @param facet_projects         - show facets  on projects
#' @param facet_replicates       - show facets on replicates 
#' @param project                - projects (sequencing runs) to include (vector)
#' @param conditions             - experimental conditions to include (vector)
#' @param primers                - RACE primers to include (vector)
#' @param cell_lines             - cell lines 
#' @param localizations          - localizations (from fractionation experiments) to include
#' @param persons                - persons preparing libraries
#' @param facet_by_spec          - specify faceting options
#' @param transcript2            - transcript to process
#'
#' @return                       - plot with mapping positions
#'
#' @return
#' @export
#'
#' @examples
plot_mapping_positions2 <- function(dataset,transcript2,exp_type2,conditions = NA,mapping_position_min = NA,mapping_position_max = NA,project = NA, facet_projects = FALSE,facet_replicates = FALSE,primers=NA,cell_lines=NA,localizations=NA,persons=NA,facet_by_spec=NA)
  {

  facet_by <- as.vector(facet_by_spec)
  if (!is.na(facet_by) & length(facet_by)==1) 
    {
    facet_by<-c(facet_by,".")
    }
  
  #print(facet_by)
  
  #filter input dataset to include only required transcript and experiment type
  test_trans <- filter_data(dataset,transcript2 = transcript2,exp_type2 = exp_type2,conditions = conditions,project = project,mapping_position_max = mapping_position_max,mapping_position_min = mapping_position_min,cell_lines = cell_lines,primers = primers,persons = persons, localizations = localizations)
  test_trans$condition <- as.character(test_trans$condition)
  test_trans$condition <- as.factor(test_trans$condition)
  
  #print(test_trans)
  
  #create plot
  plot_out <- ggplot(test_trans,aes(x=mapping_position,group=condition,colour=condition)) 
  plot_out <- plot_out + stat_bin(binwidth=1,aes(x=mapping_position,y=..ncount..),geom="line") 
  plot_out <- plot_out +  ggtitle(plot_title)
  plot_out <- plot_out + xlab("mapping position")
  plot_out <- plot_out + ylab("count")
  #create facets (if specified)
  #create facets 
  if (!is.na(facet_by)) {
    plot_out <- plot_out + facet_grid (reformulate(facet_by[2:length(facet_by)],facet_by))
  }
  #print(plot_out)
  return(plot_out)
} 



#' Title
#'
#' @param dataset                - input dataset 
#' @param exp_type2              - experiment type (OVR,KD,LEAP,NT)
#' @param conditions             - conditions to include in the analysis
#' @param mapping_position_min   - minimal mapping position to include on plot
#' @param mapping_position_max   - maximum mapping position to include on plot
#' @param facet_projects         - show facets  on projects
#' @param facet_replicates       - show facets on replicates 
#' @param project                - projects (sequencing runs) to include (vector)
#' @param conditions             - experimental conditions to include (vector)
#' @param primers                - RACE primers to include (vector)
#' @param cell_lines             - cell lines 
#' @param localizations          - localizations (from fractionation experiments) to include
#' @param persons                - persons preparing libraries
#' @param facet_by_spec          - specify faceting options
#' @param transcript2            - transcript to process
#' @param grouping_var           - grouping variable for x axis (e.g. condition)
#'
#' @return                       - plot with distribution of tail types
#' @export
#'
#' @examples
plot_tail_types_distribution <- function(dataset,transcript2,exp_type2,conditions = NA,mapping_position_min = NA,mapping_position_max = NA,project = NA, facet_projects = FALSE,facet_replicates = FALSE,primers=NA,cell_lines=NA,localizations=NA,persons=NA,facet_by_spec=NA,grouping_var = "condition")
{
  
  facet_by <- as.vector(facet_by_spec)
  if (!is.na(facet_by) & length(facet_by)==1) 
  {
    facet_by<-c(facet_by,".")
  }
  
  #print(facet_by)
  
  #filter input dataset to include only required transcript and experiment type
  test_trans <- filter_data(dataset,transcript2 = transcript2,exp_type2 = exp_type2,conditions = conditions,project = project,mapping_position_max = mapping_position_max,mapping_position_min = mapping_position_min,cell_lines = cell_lines,primers = primers,persons = persons, localizations = localizations)
  test_trans$condition <- as.character(test_trans$condition)
  test_trans$condition <- as.factor(test_trans$condition)
  
  #print(test_trans)

  
    
  if (is.na(facet_by)) {
    test_trans <- test_trans %>% group_by(.dots = c(grouping_var,"tail_type")) %>% dplyr::summarize(n_tail=n()) %>% ungroup() %>% dplyr::group_by(.dots = c(grouping_var)) %>% dplyr::mutate(freq = n_tail/sum(n_tail))
  }  else {
    test_trans <- test_trans %>% group_by(.dots = c(grouping_var,"project_name","tail_type",facet_by[facet_by!="."]))
    #if using facets - group by projects also
    #test_trans <- test_trans %>% dplyr::summarize(n_urid=n()) %>% ungroup() %>% dplyr::group_by(condition,replicate,project_name) %>% dplyr::mutate(freq_urid = n_urid/sum(n_urid)) %>% dplyr::group_by(condition,uridylated2,project_name) %>% dplyr::mutate(mean_freq_urid = mean(freq_urid), sd_urid = sd(freq_urid))
    test_trans <- test_trans %>% dplyr::summarize(n_tail=n()) %>% ungroup() %>% dplyr::group_by(.dots = c(grouping_var,facet_by[facet_by!="."])) %>% dplyr::mutate(freq = n_tail/sum(n_tail))
  }
  
  #create plot
  plot_out <- ggplot(test_trans,aes_string(x=grouping_var,group="tail_type",fill="tail_type")) 
  plot_out <- plot_out + geom_bar(aes(y=freq),stat="identity",position="stack")
  #plot_out <- plot_out +  ggtitle(plot_title)
  plot_out <- plot_out + xlab(grouping_var)
  plot_out <- plot_out + ylab("frequency")
  plot_out <- plot_out + scale_fill_grey()
  #create facets (if specified)
  #create facets 
  if (!is.na(facet_by)) {
    plot_out <- plot_out + facet_grid (reformulate(facet_by[2:length(facet_by)],facet_by))
  }
  #print(plot_out)
  return(plot_out)
} 


#' Title
#'
#' @param dataset                - input dataset 
#' @param exp_type2              - experiment type (OVR,KD,LEAP,NT)
#' @param conditions             - conditions to include in the analysis
#' @param mapping_position_min   - minimal mapping position to include on plot
#' @param mapping_position_max   - maximum mapping position to include on plot
#' @param facet_projects         - show facets  on projects
#' @param facet_replicates       - show facets on replicates 
#' @param project                - projects (sequencing runs) to include (vector)
#' @param conditions             - experimental conditions to include (vector)
#' @param primers                - RACE primers to include (vector)
#' @param cell_lines             - cell lines 
#' @param localizations          - localizations (from fractionation experiments) to include
#' @param persons                - persons preparing libraries
#' @param facet_by_spec          - specify faceting options
#' @param transcript2            - transcript to process
#' @param grouping_var           - grouping variable (e.g condition)
#'
#' @return
#' @export
#'
#' @examples
plot_Utail_lengths_distribution <- function(dataset,transcript2,exp_type2,conditions = NA,mapping_position_min = NA,mapping_position_max = NA,project = NA, facet_projects = FALSE,facet_replicates = FALSE,primers=NA,cell_lines=NA,localizations=NA,persons=NA,facet_by_spec=NA,grouping_var = "condition")
{
  
  facet_by <- as.vector(facet_by_spec)
  if (!is.na(facet_by) & length(facet_by)==1) 
  {
    facet_by<-c(facet_by,".")
  }
  
  #print(facet_by)
  
  #filter input dataset to include only required transcript and experiment type
  test_trans <- filter_data(dataset,transcript2 = transcript2,exp_type2 = exp_type2,conditions = conditions,project = project,mapping_position_max = mapping_position_max,mapping_position_min = mapping_position_min,cell_lines = cell_lines,primers = primers,persons = persons, localizations = localizations)
  test_trans$condition <- as.character(test_trans$condition)
  test_trans$condition <- as.factor(test_trans$condition)
  
  #print(test_trans)
  
  
  
  if (is.na(facet_by)) {
    test_trans <- test_trans %>% group_by(.dots = c(grouping_var,"U_length")) %>% dplyr::summarize(n_tail=n()) %>% ungroup() %>% dplyr::group_by(.dots = c(grouping_var)) %>% dplyr::mutate(freq = n_tail/sum(n_tail))
  }  else {
    test_trans <- test_trans %>% group_by(.dots = c(grouping_var,"project_name","U_length",facet_by[facet_by!="."]))
    #if using facets - group by projects also
    #test_trans <- test_trans %>% dplyr::summarize(n_urid=n()) %>% ungroup() %>% dplyr::group_by(condition,replicate,project_name) %>% dplyr::mutate(freq_urid = n_urid/sum(n_urid)) %>% dplyr::group_by(condition,uridylated2,project_name) %>% dplyr::mutate(mean_freq_urid = mean(freq_urid), sd_urid = sd(freq_urid))
    test_trans <- test_trans %>% dplyr::summarize(n_tail=n()) %>% ungroup() %>% dplyr::group_by(.dots = c(grouping_var,facet_by[facet_by!="."])) %>% dplyr::mutate(freq = n_tail/sum(n_tail))
  }
  
  #create plot
  plot_out <- ggplot(test_trans,aes_string(x=grouping_var,group="U_length",fill="U_length")) 
  plot_out <- plot_out + geom_bar(aes(y=freq),stat="identity",position="stack")
  #plot_out <- plot_out +  ggtitle(plot_title)
  plot_out <- plot_out + xlab(grouping_var)
  plot_out <- plot_out + ylab("frequency")
  plot_out <- plot_out + scale_fill_grey()
  #create facets (if specified)
  #create facets 
  if (!is.na(facet_by)) {
    if(length(facet_by_spec)==1) 
    {
      plot_out <- plot_out + facet_wrap(facet_by_spec)
    }
    else 
    {
      plot_out + facet_grid (reformulate(facet_by[2:length(facet_by)], facet_by)) 
    }
  }
  #print(plot_out)
  return(plot_out)
} 

#' Title
#'
#' @param dataset                - input dataset 
#' @param exp_type2              - experiment type (OVR,KD,LEAP,NT)
#' @param conditions             - conditions to include in the analysis
#' @param mapping_position_min   - minimal mapping position to include on plot
#' @param mapping_position_max   - maximum mapping position to include on plot
#' @param facet_projects         - show facets  on projects
#' @param facet_replicates       - show facets on replicates 
#' @param project                - projects (sequencing runs) to include (vector)
#' @param conditions             - experimental conditions to include (vector)
#' @param primers                - RACE primers to include (vector)
#' @param cell_lines             - cell lines 
#' @param localizations          - localizations (from fractionation experiments) to include
#' @param persons                - persons preparing libraries
#' @param facet_by_spec          - specify faceting options
#' @param transcript2            - transcript to process
#' @param grouping_var           - grouping variable (e.g. condition)
#'
#' @return
#' @export
#'
#' @examples
plot_tail_types_distribution_localization <- function(dataset,transcript2,exp_type2,conditions = NA,mapping_position_min = NA,mapping_position_max = NA,project = NA, facet_projects = FALSE,facet_replicates = FALSE,primers=NA,cell_lines=NA,localizations=NA,persons=NA,facet_by_spec=NA,grouping_var = "condition")
{
  
  facet_by <- as.vector(facet_by_spec)
  if (!is.na(facet_by) & length(facet_by)==1) 
  {
    facet_by<-c(facet_by,".")
  }
  
  #print(facet_by)
  
  #filter input dataset to include only required transcript and experiment type
  test_trans <- filter_data(dataset,transcript2 = transcript2,exp_type2 = exp_type2,conditions = conditions,project = project,mapping_position_max = mapping_position_max,mapping_position_min = mapping_position_min,cell_lines = cell_lines,primers = primers,persons = persons, localizations = localizations)
  test_trans$condition <- as.character(test_trans$condition)
  test_trans$condition <- as.factor(test_trans$condition)
  
  #print(test_trans)
  
  
  
  if (is.na(facet_by)) {
    test_trans <- test_trans %>% group_by(.dots = c(grouping_var,"tail_type")) %>% dplyr::summarize(n_tail=n()) %>% ungroup() %>% group_by_(grouping_var) %>% mutate(n_grouping_var=sum(n_tail)) %>% ungroup() %>% dplyr::mutate(freq = n_tail/sum(n_tail))
  }  else {
    test_trans <- test_trans %>% group_by(.dots = c(grouping_var,"project_name","tail_type",facet_by[facet_by!="."]))
    #if using facets - group by projects also
    #test_trans <- test_trans %>% dplyr::summarize(n_urid=n()) %>% ungroup() %>% dplyr::group_by(condition,replicate,project_name) %>% dplyr::mutate(freq_urid = n_urid/sum(n_urid)) %>% dplyr::group_by(condition,uridylated2,project_name) %>% dplyr::mutate(mean_freq_urid = mean(freq_urid), sd_urid = sd(freq_urid))
    test_trans <- test_trans %>% dplyr::summarize(n_tail=n()) %>% ungroup() %>% group_by(.dots = c(grouping_var,facet_by[facet_by!="."])) %>% mutate(n_grouping_var=sum(n_tail)) %>% ungroup() %>% dplyr::group_by(.dots = c(facet_by[facet_by!="."])) %>% dplyr::mutate(freq = n_tail/sum(n_tail))
  }
 # print(test_trans)
  #create plot
  plot_out <- ggplot(test_trans,aes_string(x=grouping_var,group="tail_type",fill="tail_type")) 
  plot_out <- plot_out + geom_bar(aes(y=freq),stat="identity",position="stack")
  plot_out <- plot_out + geom_text(aes(y=1.02,label=n_grouping_var))
  #plot_out <- plot_out +  ggtitle(plot_title)
  plot_out <- plot_out + xlab(grouping_var)
  plot_out <- plot_out + ylab("frequency")
  plot_out <- plot_out + scale_fill_grey()
  #create facets (if specified)
  #create facets 
  if (!is.na(facet_by)) {
    plot_out <- plot_out + facet_grid (reformulate(facet_by[2:length(facet_by)],facet_by))
  }
  #print(plot_out)
  return(plot_out)
} 


#' Title
#'
#' @param dataset                - input dataset 
#' @param exp_type2              - experiment type (OVR,KD,LEAP,NT)
#' @param conditions             - conditions to include in the analysis
#' @param mapping_position_min   - minimal mapping position to include on plot
#' @param mapping_position_max   - maximum mapping position to include on plot
#' @param facet_projects         - show facets  on projects
#' @param facet_replicates       - show facets on replicates 
#' @param project                - projects (sequencing runs) to include (vector)
#' @param conditions             - experimental conditions to include (vector)
#' @param primers                - RACE primers to include (vector)
#' @param cell_lines             - cell lines 
#' @param localizations          - localizations (from fractionation experiments) to include
#' @param persons                - persons preparing libraries
#' @param facet_by_spec          - specify faceting options
#' @param transcript2            - transcript to process
#'
#' @return
#' @export
#'
#' @examples
plot_tail_types_distribution_by_length <- function(dataset,transcript2,exp_type2,conditions = NA,mapping_position_min = NA,mapping_position_max = NA,project = NA, facet_projects = FALSE,facet_replicates = FALSE,primers=NA,cell_lines=NA,localizations=NA,persons=NA,facet_by_spec=NA,grouping_var = "condition",plot_bar=TRUE,plot_line=TRUE)
{
  output=list()
  facet_by <- as.vector(facet_by_spec)
  if (!is.na(facet_by) & length(facet_by)==1) 
  {
    facet_by<-c(facet_by,".")
  }
  
  #print(facet_by)
  
  #filter input dataset to include only required transcript and experiment type
  test_trans <- filter_data(dataset,transcript2 = transcript2,exp_type2 = exp_type2,conditions = conditions,project = project,mapping_position_max = mapping_position_max,mapping_position_min = mapping_position_min,cell_lines = cell_lines,primers = primers,persons = persons, localizations = localizations)
  test_trans$condition <- as.character(test_trans$condition)
  test_trans$condition <- as.factor(test_trans$condition)
  
  #print(test_trans)

  #test_trans
  
    
  if (is.na(facet_by)) {
    test_trans <- test_trans %>% group_by(.dots = c(grouping_var,"tail_type","tail_length_fac")) %>% dplyr::summarize(n_tail=n()) %>% ungroup() %>% dplyr::group_by(.dots = c(grouping_var,"tail_length_fac")) %>% dplyr::mutate(freq = n_tail/sum(n_tail),sum_length_tails=sum(n_tail)) %>% dplyr::group_by(.dots = c(grouping_var)) %>% dplyr::mutate(freq_length = sum_length_tails/sum(n_tail))
  }  else {
    test_trans <- test_trans %>% group_by(.dots = c(grouping_var,"project_name","tail_type","tail_length_fac",facet_by[facet_by!="."]))
    #if using facets - group by projects also
    #test_trans <- test_trans %>% dplyr::summarize(n_urid=n()) %>% ungroup() %>% dplyr::group_by(condition,replicate,project_name) %>% dplyr::mutate(freq_urid = n_urid/sum(n_urid)) %>% dplyr::group_by(condition,uridylated2,project_name) %>% dplyr::mutate(mean_freq_urid = mean(freq_urid), sd_urid = sd(freq_urid))
    test_trans <- test_trans %>% dplyr::summarize(n_tail=n()) %>% ungroup() %>% dplyr::group_by(.dots = c("tail_length_fac",grouping_var,facet_by[facet_by!="."])) %>% dplyr::mutate(freq = n_tail/sum(n_tail),sum_length_tails=sum(n_tail)) %>% dplyr::group_by(.dots = c(grouping_var,facet_by[facet_by!="."])) %>% dplyr::mutate(freq_length = sum_length_tails/sum(n_tail))
  }
  #print(test_trans)
  output$calculated_values = test_trans
  #create plot
  plot_out <- ggplot(test_trans,aes_string(x="tail_length_fac"))
  if(plot_bar) {
    plot_out <- plot_out + geom_bar(aes_string(y="freq",group="tail_type",fill="tail_type"),stat="identity",position="stack")
  }
  if (plot_line) {
    plot_out <- plot_out + geom_line(data=filter(test_trans,tail_type=="A_only"),aes(y=freq_length,group=1),size=0.9,linetype=5,alpha=0.7) 
    plot_out <- plot_out + geom_point(data=filter(test_trans,tail_type=="A_only"),aes(y=freq_length,group=1),size=1.5)
  }
  
  #plot_out <- plot_out +  ggtitle(plot_title)
  plot_out <- plot_out + xlab("tail_length")
  plot_out <- plot_out + ylab("frequency")
  plot_out <- plot_out + scale_fill_grey()
  #create facets (if specified)
  #create facets 
  if (!is.na(facet_by)) 
  {
    if(length(facet_by_spec)==1) 
    {
      plot_out <- plot_out + facet_wrap(facet_by_spec)
    }
    else 
    {
      plot_out + facet_grid (reformulate(facet_by[2:length(facet_by)], facet_by)) 
    }
  }
  #print(plot_out)
  output$plot <- plot_out
  return(output)
} 



#' Plot tail length disribution
#'
#' @param dataset               - input dataset
#' @param trans                 - transcript to analyze
#' @param exp_type2             - experiment type (OVR,KD,LEAP,NT)
#' @param conditions            - conditions to include (vector)
#' @param mapping_position_min  - minimal mapping position to include
#' @param mapping_position_max  - maximal mapping position to include
#' @param facet_conditions      - facet by conditions
#' @param facet_replicates      - facet by replicates
#' @param max_tail_length       - max tail length to include
#' @param min_tail_length       - min tail length to include
#' @param facet_tail_types      - facet by tail types
#' @param AUtail                - specify A or U tail part only ('Atail' or 'Utail', other values ignored)
#'
#' @return                      - plot with tail lengths distribution
#' @export
#'
#' @examples
plot_tail_length_distribution <- function(dataset,trans,exp_type2,conditions = NA,mapping_position_min = NA,mapping_position_max = NA,facet_conditions = FALSE,facet_replicates = FALSE,max_tail_length = NA, min_tail_length = NA, facet_tail_types = FALSE,AUtail = NA,binsize=1) {
  test_trans <- dataset %>% filter(transcript==trans)
  test_trans <- test_trans %>% filter(exp_type==exp_type2)
  
  plot_title <- trans
  if (!is.na(AUtail)) {
    if (AUtail == 'Atail') {
      test_trans <- test_trans %>% mutate(tail_length = Atail_length)
      plot_title <- paste(plot_title,"A tail length only",sep=",")
    } else if (AUtail =='Utail') {
      test_trans <- test_trans %>% mutate(tail_length = Utail_length)
      plot_title <- paste(plot_title,"U tail length only",sep=",")
    }
  }
  #if min and max mapping positions are specified - used for filtering
  if (!is.na(mapping_position_min)) {
    test_trans <- test_trans %>% filter(mapping_position>mapping_position_min)
  }
  if (!is.na(mapping_position_max)) {
    test_trans <- test_trans %>% filter(mapping_position<mapping_position_max)
  }
  if (!is.na(min_tail_length)) {
    test_trans <- test_trans %>% filter(tail_length>min_tail_length)
  }
  if (!is.na(max_tail_length)) {
    test_trans <- test_trans %>% filter(tail_length<max_tail_length)
  }
  #if conditions are specified - use for filtering
  if(length(conditions)>0 & !is.na(conditions)) {
    test_trans <- test_trans %>% filter(condition %in% conditions)
  }

  test_trans$condition <- as.character(test_trans$condition)
  test_trans$condition <- as.factor(test_trans$condition)
  test_trans$tail_type <- as.character(test_trans$tail_type)
  test_trans$tail_type <- as.factor(test_trans$tail_type)
  

  
  plot_out <- ggplot(test_trans,aes(x=tail_length,group=tail_type,colour=tail_type)) 
  plot_out <- plot_out + stat_bin(binwidth=binsize,aes(x=tail_length,y=..ncount..,colour=tail_type),geom="line") 
  plot_out <- plot_out +  ggtitle(plot_title)
  #create facets 
  if (facet_conditions==TRUE) {
    plot_out <- plot_out + facet_grid (condition ~ .)
  }
  else if (facet_replicates==TRUE) {
    plot_out <- plot_out + facet_grid (replicate ~ .)
  }
  else if (facet_tail_types==TRUE) {
    plot_out <- plot_out + facet_grid (tail_type ~ .)
  }
 # print(plot_out)
  return(plot_out)
}


#' Plot tail length disribution 2
#'
#' @param dataset               - input dataset
#' @param trans                 - transcript to analyze
#' @param exp_type2             - experiment type (OVR,KD,LEAP,NT)
#' @param conditions            - conditions to include (vector)
#' @param mapping_position_min  - minimal mapping position to include
#' @param mapping_position_max  - maximal mapping position to include
#' @param facet_conditions      - facet by conditions
#' @param facet_replicates      - facet by replicates
#' @param max_tail_length       - max tail length to include
#' @param min_tail_length       - min tail length to include
#' @param facet_tail_types      - facet by tail types
#' @param AUtail                - specify A or U tail part only ('Atail' or 'Utail', other values ignored)
#'
#' @return                      - plot with tail lengths distribution
#' @export
#'
#' @examples
plot_tail_length_distribution2 <- function(dataset,trans,exp_type2,conditions = NA,mapping_position_min = NA,mapping_position_max = NA,facet_conditions = FALSE,facet_replicates = FALSE,max_tail_length = NA, min_tail_length = NA, facet_tail_types = FALSE,AUtail = NA,binsize=1) {
  test_trans <- dataset %>% filter(transcript==trans)
  test_trans <- test_trans %>% filter(exp_type==exp_type2)
  
  plot_title <- trans
  if (!is.na(AUtail)) {
    if (AUtail == 'Atail') {
      test_trans <- test_trans %>% mutate(tail_length = Atail_length)
      plot_title <- paste(plot_title,"A tail length only",sep=",")
    } else if (AUtail =='Utail') {
      test_trans <- test_trans %>% mutate(tail_length = Utail_length)
      plot_title <- paste(plot_title,"U tail length only",sep=",")
    }
  }
  #if min and max mapping positions are specified - used for filtering
  if (!is.na(mapping_position_min)) {
    test_trans <- test_trans %>% filter(mapping_position>mapping_position_min)
  }
  if (!is.na(mapping_position_max)) {
    test_trans <- test_trans %>% filter(mapping_position<mapping_position_max)
  }
  if (!is.na(min_tail_length)) {
    test_trans <- test_trans %>% filter(tail_length>min_tail_length)
  }
  if (!is.na(max_tail_length)) {
    test_trans <- test_trans %>% filter(tail_length<max_tail_length)
  }
  #if conditions are specified - use for filtering
  if(length(conditions)>0 & !is.na(conditions)) {
    test_trans <- test_trans %>% filter(condition %in% conditions,tail_length>0)
  }
  
  test_trans$condition <- as.character(test_trans$condition)
  test_trans$condition <- as.factor(test_trans$condition)
  test_trans$tail_type <- as.character(test_trans$tail_type)
  test_trans$tail_type <- as.factor(test_trans$tail_type)
  
  
  
  plot_out <- ggplot(test_trans,aes(x=tail_length,group=condition,colour=condition)) 
  plot_out <- plot_out + stat_bin(binwidth=binsize,aes(x=tail_length,y=..count..,colour=condition),geom="line") 
  plot_out <- plot_out +  ggtitle(plot_title)
  #create facets 
  if (facet_conditions==TRUE) {
    plot_out <- plot_out + facet_grid (condition ~ .)
  }
  else if (facet_replicates==TRUE) {
    plot_out <- plot_out + facet_grid (replicate ~ .)
  }
  else if (facet_tail_types==TRUE) {
    plot_out <- plot_out + facet_grid (tail_type ~ .)
  }
  # print(plot_out)
  return(plot_out)
}



#' Plot Mean tail lengths
#'
#' @param dataset               - input dataset
#' @param exp_type2             - experiment type (OVR,KD,LEAP,NT)
#' @param conditions            - conditions to include (vector)
#' @param mapping_position_min  - minimal mapping position to include
#' @param mapping_position_max  - maximal mapping position to include
#' @param AUtail                - specify A or U tail part only ('Atail' or 'Utail', other values ignored)
#' @param uridylated_only       - include uridyalted reads only
#' @param include_jitter        - show jittered point plot
#' @param project 
#' @param transcript2 
#' @param facet_by_spec 
#' @param localizations 
#' @param persons 
#' @param primers 
#' @param cell_lines 
#' @param tailed_only 
#'
#' @return                      - plot with tail lengths distribution
#' @export
#'
#' @examples
plot_mean_tail_lengths <-
  function(dataset,
           transcript2,
           exp_type2,
           mapping_position_min = NA,
           mapping_position_max = NA,
           project = NA,
           facet_by_spec = NA,
           conditions = NA,
           include_jitter = TRUE,
           localizations = NA,
           persons = NA,
           primers = NA,
           cell_lines = NA,
           uridylated_only = FALSE,
           AUtail = NA,
           tailed_only = FALSE)
  {
    output = list()
    facet_by <- as.vector(facet_by_spec)
    if (!is.na(facet_by) & length(facet_by) == 1) {
      facet_by <- c(facet_by, ".")
    }
    output = list() #list for storing output
    message("Getting data")
    test_trans <-
      filter_data(
        dataset,
        transcript2 = transcript2,
        exp_type2 = exp_type2,
        conditions = conditions,
        project = project,
        mapping_position_max = mapping_position_max,
        mapping_position_min = mapping_position_min,
        cell_lines = cell_lines,
        primers = primers,
        persons = persons,
        localizations = localizations
      )
    assertthat::assert_that(nrow(test_trans)>0)
    plot_title_elements <-
      c(
        transcript2,
        exp_type2,
        paste(project, collapse = ","),
        paste(cell_lines, collapse = ","),
        paste(primers, collapse = ",")
      )
    plot_title_elements <-
      plot_title_elements[plot_title_elements != "NA"]
    if (!is.na(mapping_position_min) | !is.na(mapping_position_max)) {
      plot_title_elements <-
        c(plot_title_elements,
          "(",
          mapping_position_min,
          "-",
          mapping_position_max,
          ")")
    }
    plot_title <-
      paste(c("mean tail lengths", plot_title_elements), collapse = " ")
    if (!is.na(facet_by)) {
      groups_to_exclude <-
        test_trans %>% group_by(.dots = c(facet_by_spec, "condition")) %>% count() %>% group_by_(.dots = c(facet_by_spec)) %>% count() %>% filter(nn %in% c(0, 1)) %>% select_(facet_by_spec)
      if (nrow(groups_to_exclude) > 0) {
        #print(groups_to_exclude)
        message("Excluding groups with zero or one condition")
        test_trans <- test_trans %>% anti_join(groups_to_exclude)
      }
    }
    
    #plot_title <- transcript2
    if (!is.na(AUtail)) {
      if (AUtail == 'Atail') {
        test_trans <- test_trans %>% mutate(tail_length = Atail_length)
        plot_title <- paste(plot_title, "A tail length only", sep = ",")
      } else if (AUtail == 'Utail') {
        test_trans <- test_trans %>% mutate(tail_length = Utail_length)
        plot_title <- paste(plot_title, "U tail length only", sep = ",")
      }
    }
    
    
    if (uridylated_only) {
      test_trans <- test_trans %>% filter(uridylated2 == TRUE)
      plot_title <- paste(plot_title, "uridylated only", sep = ", ")
      
    }
    
    if (tailed_only) {
      test_trans <- test_trans %>% filter(tail_length>0)
      plot_title <- paste(plot_title, "tailed only", sep = ", ")
    }

    message("Calculating mean tail lengths")
    tail_lengths  <-
      calculate_mean_tail_lengths(
        test_trans,
        transcript2 = transcript2,
        exp_type2 = exp_type2,
        conditions = conditions,
        project = project,
        mapping_position_max = mapping_position_max,
        mapping_position_min = mapping_position_min,
        cell_lines = cell_lines,
        primers = primers,
        persons = persons,
        localizations = localizations,
        facet_by_spec = facet_by_spec
      )
    output$mean_tail_lengths <- tail_lengths

    
    plot_out <-
      plot_stats(
        tail_lengths,
        values = "mean_tail_length_rep",
        grouping_var = "condition",
        group_by = facet_by_spec,
        test = "Tukey"
      )
    plot_out$plot <- plot_out$plot + ylab("mean tail length")
    plot_out$plot <- plot_out$plot + ggtitle(plot_title)
    output$plot <- plot_out
    return(output)
  }

#' Title
#'
#' @param uridylation_data (input file with calculated uridylation data)
#' @param output_file - path to the output file (optional, if missing than guessed from available conditions)
#' @param grouping_var - grouping variable (e.g. condition)
#'
#' @return writes xlsx file
#' @export
#'
#' @examples
write_uridylation_xlsx <- function(uridylation_data,output_file = NA,grouping_var = "condition") 
  {
  colnames_input <- colnames(uridylation_data)
  base_cols <- c(grouping_var,"freq_urid","replicate","uridylated2","n_urid","mean_freq_urid","sd_urid")
  grouping_cols <- setdiff(colnames_input,base_cols)
  #print(grouping_cols)
  output_data = uridylation_data %>% dplyr::ungroup() %>% dplyr::select(!!!rlang::syms(c(grouping_cols,grouping_var,"freq_urid","replicate"))) %>% as.data.frame()  
  #print(output_data)
  output_data <- output_data %>% stats::reshape(direction="wide",timevar = grouping_var,idvar=c(grouping_cols,"replicate"))
  #print(output_data)
  colnames(output_data) <- gsub("(freq_urid.)?(.*)","\\2",colnames(output_data))
  if (!is.na(output_file)) {
    xlsx::write.xlsx(output_data,output_file)
  } else {
    output_file = paste(paste(levels(uridylation_data$calculated_values$condition),collapse = "_"),".xlsx",sep="")
    xlsx::write.xlsx(output_data,output_file)
  }
} 


#' Title
#'
#' @param input_dataset - input dataset
#' @param n             - number of nucleotides to take for terminal
#'
#' @return dataset with calculated terminal nucleotides as term2 column
#' @export
#'
#' @examples
generate_terminal_nucl <- function(input_dataset,n=30) {
  LINE_seq_to_search = "TACTATTAGCCCGGGC"
  x2 <- input_dataset[grepl(LINE_seq_to_search, input_dataset$R5_seq), ]
  x2 <-  x2 %>% mutate(R3_seq2 = paste(R3_seq, R3_clip, sep = "")) %>% mutate(term2 = str_sub(R3_seq2, -n, -1))
  return(x2)
}
