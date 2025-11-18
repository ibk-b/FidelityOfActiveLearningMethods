library(igraph) 
library(dplyr) 
library(ggplot2)
library(ggnewscale)
library(gridExtra)
library(grid)
library(tidyselect)
library(tidyverse) 
library(effsize)
library(metafor)
library(ggpubr)
library(ggsignif)
library(reshape2)
library(scales)
library(patchwork)
library(cowplot)


#Functions------------------------------------------

normalize_vertex_names <- function(g) {
  if ("L" %in% V(g)$name) {
    # Case 1: when "L" is one of the vertex names (student codes)
    V(g)$name <- gsub("^IND$", "Ind", V(g)$name)
    V(g)$name <- gsub("^AnQ$", "AnQ_S", V(g)$name)
  } else if ("Lec" %in% V(g)$name) {
    # Case 2: when "Lec" is one of the vertex names (instructor codes)
    V(g)$name <- gsub("^Fup$", "FUp", V(g)$name)
    V(g)$name <- gsub("^AnQ$", "AnQ_I", V(g)$name)
    V(g)$name <- gsub("^D\\.V$", "D/V", V(g)$name)
  }
  return(g)
}

create_weighted_transition_network <- function(df) {
  edge_list <- list()
  
  for (i in seq_len(max(0, nrow(df) - 1))) {
    active_cols_row1 <- names(df)[which(df[i, ] == 1)]
    active_cols_row2 <- names(df)[which(df[i + 1, ] == 1)]
    
    # Exclude columns from row2 that were already active in row1
    new_active_cols_row2 <- setdiff(active_cols_row2, active_cols_row1)
    
    if (length(active_cols_row1) > 0 && length(new_active_cols_row2) > 0) {
      combinations <- expand.grid(active_cols_row1, new_active_cols_row2)
      edge_list[[i]] <- as.matrix(combinations)
    }
  }
  
  all_edges <- do.call(rbind, edge_list)
  edge_df <- as.data.frame(all_edges, stringsAsFactors = FALSE)
  names(edge_df) <- c("from", "to")
  
  edge_weights <- edge_df %>%
    group_by(from, to) %>%
    summarise(weight = n() / (nrow(df)-1), .groups = "drop")
  
  proportions <- colSums(df == 1, na.rm = TRUE) / nrow(df)
  
  result <- data.frame(
    Column = names(proportions),
    ProportionOfOnes = proportions,
    row.names = NULL
  )
  
  if (length(edge_list) == 0) {
    g <- make_empty_graph(n = ncol(df), directed = TRUE)
    g <- normalize_vertex_names(g)
    V(g)$name <- colnames(df)
    V(g)$freq <- colSums(df == 1, na.rm = TRUE) / nrow(df)
    return(g)
  }

  g <- graph_from_data_frame(edge_weights, directed = TRUE, vertices = colnames(df))
  g <- normalize_vertex_names(g)
  
  E(g)$weight <- edge_weights$weight
  V(g)$freq <- result$ProportionOfOnes
  
  V(g)$label.family <- "sans"
  
  return(g)
}

aggregate_student_CALEP1_transition_network <- function(folder, Instructors, title = NULL) {
  if (!is.vector(Instructors)) {
    Instructors <- c(Instructors)
  } 
  
  dat_list <- list()
  
  for (inst in Instructors) {
    # If inst already includes "Class", match exactly that file
    if (grepl("Class[0-9]+$", inst)) {
      files <- list.files(folder, pattern = paste0("^", inst, "\\.csv$"),
                          full.names = TRUE)
    } else {
      # Otherwise, treat inst as a method prefix (like "CALEP1_Tutorials")
      files <- list.files(folder, pattern = paste0("^", inst, "_Class[0-9]+\\.csv$"),
                          full.names = TRUE)
    }
    
    for (f in files) {
      dat <- read.csv(f, stringsAsFactors = FALSE)
      dat <- dat[5:nrow(dat), c(2:8, 10:11)]
      colnames(dat) <- c("L","IND","CG","WG","OG","AnQ","SQ","Prd","SP")
      
      dat <- dat[!apply(dat == "", 1, all),]
      dat[dat == "X"] <- 1
      dat[dat == ""] <- 0
      dat <- as.data.frame(lapply(dat, as.numeric))
      dat_list[[length(dat_list) + 1]] <- dat
    }
  }
  
  dat_all <- do.call(rbind, dat_list)
  g <- create_weighted_transition_network(dat_all)
  
  desired_order <- c("L","Ind","CG","WG","OG","AnQ_S","SQ","Prd","SP")
  
  coords_custom <- layout_in_circle(g, order = intersect(desired_order, V(g)$name))
  
  color_palette <- colorRampPalette(c("white", "cornflowerblue"))
  node_colors <- color_palette(100)[as.numeric(cut(V(g)$freq, breaks=100))]
  
  if (vcount(g) > 0 && ecount(g) > 0) {
    plot(g,
         edge.arrow.size = 0.1,
         vertex.color = node_colors,
         edge.color = "gray20",
         layout = coords_custom,
         edge.curved = 0.1,
         vertex.size = 45,
         vertex.size2 = 20,
         edge.width = E(g)$weight * 30,
         vertex.label.cex = 1,
         vertex.label.color = "black",
         vertex.shape = "rectangle",
         main = paste("CALEP1", ifelse(is.null(title), inst, title), "(Student)")
    )
  } 
  return(g)
}


aggregate_student_CALEP2_transition_network <- function(Instructor) {
  
  dat<-read.csv(paste0("CALEP2-COPUS-observations/", trimws(Instructor),"_Class1.csv"))
  
  dat<-dat[9:length(dat$X),c(2:8, 10:11)]
  colnames(dat) <- c("L","Ind","CG","WG","OG","AnQ_S","SQ","Prd","SP")
  dat <- dat[dat$L!="L",]
  dat <- dat[dat$L!="1. Students doing",]
  dat<-dat[!apply(dat == "", 1, all),]
  
  dat2<-read.csv(paste0("CALEP2-COPUS-observations/", trimws(Instructor),"_Class2.csv"))
  
  dat2<-dat2[9:length(dat2$X),c(2:8, 10:11)]
  colnames(dat2) <- c("L","Ind","CG","WG","OG","AnQ_S","SQ","Prd","SP")
  dat2 <- dat2[dat2$L!="L",]
  dat2 <- dat2[dat2$L!="1. Students doing",]
  dat2<-dat2[!apply(dat2 == "", 1, all),]
  
  dat3<-read.csv(paste0("CALEP2-COPUS-observations/", trimws(Instructor),"_Class3.csv"))
  
  dat3<-dat3[9:length(dat3$X),c(2:8, 10:11)]
  colnames(dat3) <- c("L","Ind","CG","WG","OG","AnQ_S","SQ","Prd","SP")
  dat3 <- dat3[dat3$L!="L",]
  dat3 <- dat3[dat3$L!="1. Students doing",]
  dat3<-dat3[!apply(dat3 == "", 1, all),]
  
  dat_all<-rbind(dat, dat2, dat3)
  
  g<-create_weighted_transition_network(dat_all)
  
  desired_order <- c("L","Ind","CG","WG","OG","AnQ_S","SQ","Prd","SP")
  
  coords_custom <- layout_in_circle(g, order = intersect(desired_order, V(g)$name) )
  l <- layout_with_fr(g)
  
  color_palette <- colorRampPalette(c("white", "cornflowerblue"))
  node_colors <- color_palette(100)[as.numeric(cut(V(g)$freq, breaks=100))]
  
  if (vcount(g) > 0 && ecount(g) > 0) {
    plot(g,edge.arrow.size=0.1,vertex.color=node_colors,edge.color="gray20", layout=coords_custom,edge.curved=0.1,vertex.size=45,vertex.size2=20,edge.width=E(g)$weight*30,vertex.label.cex=1,vertex.label.color="black",vertex.shape="rectangle",main=paste("CALEP2", Instructor, "(Student)"))
  }
  return(g)
}

aggregate_CALEP1_transition_network <- function(folder, Instructors, title = NULL) {
  if (!is.vector(Instructors)) {
    Instructors <- c(Instructors)
  }
  
  dat_list <- list()
  
  for (inst in Instructors) {
    # If inst already includes "Class", match exactly that file
    if (grepl("Class[0-9]+$", inst)) {
      files <- list.files(folder, pattern = paste0("^", inst, "\\.csv$"),
                          full.names = TRUE)
    } else {
      # Otherwise, treat inst as a method prefix (like "CALEP1_Tutorials")
      files <- list.files(folder, pattern = paste0("^", inst, "_Class[0-9]+\\.csv$"),
                          full.names = TRUE)
    }
    
    for (f in files) {
      dat <- read.csv(f, stringsAsFactors = FALSE)
      dat <- dat[5:nrow(dat), c(15:21, 23:24)]
      colnames(dat) <- c("Lec","RtW","Fup","PQ","CQ","AnQ","MG","D/V","Adm")
      
      dat <- dat[!apply(dat == "", 1, all),]
      dat$PQ[dat$MG == "X"] <- "X" #CALEP2 had PQ coded each time MG was coded so this accounts for that in CALEP1  
      
      dat[dat == "Fup"] <- "FUp"
      dat[dat == "AnQ"] <- "AnQ_I"
      dat[dat == "X"] <- 1
      dat[dat == ""] <- 0
      dat <- as.data.frame(lapply(dat, as.numeric))
      dat_list[[length(dat_list) + 1]] <- dat
    }
  }
  
  dat_all <- do.call(rbind, dat_list)
  g <- create_weighted_transition_network(dat_all)
  g <- normalize_vertex_names(g)
  
  desired_order <- c("Lec","RtW","FUp","PQ","CQ","AnQ_I","MG","D/V","Adm")
  
  coords_custom <- layout_in_circle(g, order = intersect(desired_order, V(g)$name))
  
  color_palette <- colorRampPalette(c("white", "cornflowerblue"))
  node_colors <- color_palette(100)[as.numeric(cut(V(g)$freq, breaks=100))]
  
  if (vcount(g) > 0 && ecount(g) > 0) {
    plot(g,
         edge.arrow.size = 0.1,
         vertex.color = node_colors,
         edge.color = "gray20",
         layout = coords_custom,
         edge.curved = 0.1,
         vertex.size = 45,
         vertex.size2 = 20,
         edge.width = E(g)$weight * 30,
         vertex.label.cex = 1,
         vertex.label.color = "black",
         vertex.shape = "rectangle",
         main = paste("CALEP1", ifelse(is.null(title), inst, title), "(Instructor)")
    )
  } 
  return(g)
}

aggregate_CALEP2_transition_network <- function(Instructor) {
  
  dat<-read.csv(paste0("CALEP2-COPUS-observations/", trimws(Instructor),"_Class1.csv"))
  
  dat<-dat[9:length(dat$X),c(15:21,23:24)]
  colnames(dat) <- c("Lec","RtW","FUp","PQ","CQ","AnQ_I","MG","D/V","Adm")
  dat <- dat[dat$L!="L",]
  dat <- dat[dat$L!="1. Students doing",]
  dat<-dat[!apply(dat == "", 1, all),]
  
  dat2<-read.csv(paste0("CALEP2-COPUS-observations/", trimws(Instructor),"_Class2.csv"))
  
  dat2<-dat2[9:length(dat2$X),c(15:21,23:24)]
  colnames(dat2) <- c("Lec","RtW","FUp","PQ","CQ","AnQ_I","MG","D/V","Adm")
  dat2 <- dat2[dat2$L!="L",]
  dat2 <- dat2[dat2$L!="1. Students doing",]
  dat2<-dat2[!apply(dat2 == "", 1, all),]
  
  dat3<-read.csv(paste0("CALEP2-COPUS-observations/", trimws(Instructor),"_Class3.csv"))
  
  dat3<-dat3[9:length(dat3$X),c(15:21,23:24)]
  colnames(dat3) <- c("Lec","RtW","FUp","PQ","CQ","AnQ_I","MG","D/V","Adm")
  dat3 <- dat3[dat3$L!="L",]
  dat3 <- dat3[dat3$L!="1. Students doing",]
  dat3<-dat3[!apply(dat3 == "", 1, all),]
  
  dat_all<-rbind(dat, dat2, dat3)
  
  g<-create_weighted_transition_network(dat_all)
  
  desired_order <- c("Lec","RtW","FUp","PQ","CQ","AnQ_I","MG","D/V","Adm")
  
  coords_custom <- layout_in_circle(g, order = intersect(desired_order, V(g)$name) )
  l <- layout_with_fr(g)
  
  color_palette <- colorRampPalette(c("white", "cornflowerblue"))
  node_colors <- color_palette(100)[as.numeric(cut(V(g)$freq, breaks=100))]
  
  if (vcount(g) > 0 && ecount(g) > 0) {
    plot(g,edge.arrow.size=0.1,vertex.color=node_colors,edge.color="gray20", 
         layout=coords_custom,edge.curved=0.1,vertex.size=45,
         vertex.size2=20,edge.width=E(g)$weight*30,vertex.label.cex=1,
         vertex.label.color="black",vertex.shape="rectangle",
         main=paste("CALEP2", Instructor, "(Instructor)")
    )
  } 
  return(g)
}

extract_frequencies <- function(g, method, source, course = NA) {
  if (is.null(g)) return(NULL)
  data.frame(
    Method = method,
    Source = source,       # "CALEP1" or "CALEP2"
    Course = course,       # only relevant for CALEP2
    Code = V(g)$name,
    Frequency = V(g)$freq,
    stringsAsFactors = FALSE
  )
}

cosine_similarity <- function(g1, g2) {
  g1 <- normalize_vertex_names(g1)
  g2 <- normalize_vertex_names(g2)
  
  nodes <- sort(union(V(g1)$name, V(g2)$name))
  
  adj1 <- as.matrix(get.adjacency(g1, attr = "weight", sparse = FALSE))
  adj2 <- as.matrix(get.adjacency(g2, attr = "weight", sparse = FALSE))
  
  common <- intersect(nodes, intersect(rownames(adj1), rownames(adj2)))
  adj1 <- adj1[common, common, drop = FALSE]
  adj2 <- adj2[common, common, drop = FALSE]
  
  adj1_vector <- as.vector(adj1)
  adj2_vector <- as.vector(adj2)
  
  dot_product <- sum(adj1_vector * adj2_vector)
  magnitude1 <- sqrt(sum(adj1_vector^2))
  magnitude2 <- sqrt(sum(adj2_vector^2))
  
  if (magnitude1 == 0 || magnitude2 == 0) {
    return(NA)  
  }
  
  return(dot_product / (magnitude1 * magnitude2))
}

duration_cosine_similarity <- function(g1, g2) {
  g1 <- normalize_vertex_names(g1)
  g2 <- normalize_vertex_names(g2)
  
  nodes <- sort(union(V(g1)$name, V(g2)$name))
  
  freq1 <- rep(0, length(nodes)); names(freq1) <- nodes
  freq2 <- rep(0, length(nodes)); names(freq2) <- nodes
  
  freq1[V(g1)$name] <- V(g1)$freq
  freq2[V(g2)$name] <- V(g2)$freq
  
  vec1 <- c(freq1)
  vec2 <- c(freq2)
  
  dot_product <- sum(vec1 * vec2)
  magnitude1 <- sqrt(sum(vec1^2))
  magnitude2 <- sqrt(sum(vec2^2))
  
  if (magnitude1 == 0 || magnitude2 == 0) {
    return(NA)   
  }
  
  return(dot_product / (magnitude1 * magnitude2))
}

#compare for students
compare_student_methods <- function() {
  
  methods <- c("SCALEUP","ISLE", "Tutorials-WholeClass", "Tutorials-RecitationOnly")
  
  results <- data.frame(Method=character(),
                        RawCourse=character(),
                        CosineSim=numeric(),
                        DurationCosineSim=numeric(),
                        stringsAsFactors = FALSE)
  freq_data <- list() 
  for (m in methods) {
    par(mfrow = c(3, 3), mar = c(1, 1, 2, 1))
    if (m %in% c("SCALEUP", "ISLE")) {
      g1 <- aggregate_student_CALEP1_transition_network("CALEP1-COPUS-observations", paste0("CALEP1_", m), title=m)
      g1 <- normalize_vertex_names(g1)
      
      freq_data[[length(freq_data)+1]] <- extract_frequencies(g1, m, "CALEP1")
      
      raw_files <- list.files("CALEP2-COPUS-observations", pattern = paste0("^", m, "_Course[0-9]+_Class1\\.csv$"))
      raw_files <- trimws(raw_files)
      courses <- unique(gsub("_Class1\\.csv$", "", raw_files))
      
      for (course in courses) {
        g2 <- aggregate_student_CALEP2_transition_network(course)
        g2 <- normalize_vertex_names(g2)
        
        freq_data[[length(freq_data)+1]] <- extract_frequencies(g2, m, "CALEP2", course)

        sim <- cosine_similarity(g1, g2)
        duration_sim <- duration_cosine_similarity(g1,g2)

        results <- rbind(results, data.frame(Method=m,
                                             RawCourse=course,
                                             CosineSim=round (sim,2),
                                             DurationCosineSim = round(duration_sim,2),
                                             stringsAsFactors = FALSE))
      }
    } else if (m == "Tutorials-WholeClass") {
      #SPECIAL CASE: aggregated CALEP1 Class1+Class2 for CALEP2 Courses 2,3,7,9
      g1 <- aggregate_student_CALEP1_transition_network("CALEP1-COPUS-observations", 
                                                        c("CALEP1_Tutorials_Class1", "CALEP1_Tutorials_Class2"), title=m)
      g1 <- normalize_vertex_names(g1)
      freq_data[[length(freq_data)+1]] <- extract_frequencies(g1, m, "CALEP1")
      
      tutorialsA_courses <- paste0("Tutorials_Course", c(2,3,7,9))
      
      for (course in tutorialsA_courses) {
        g2 <- aggregate_student_CALEP2_transition_network(course)
        g2 <- normalize_vertex_names(g2)
        
        freq_data[[length(freq_data)+1]] <- extract_frequencies(g2, m, "CALEP2", course)
        
        sim <- cosine_similarity(g1, g2)
        duration_sim <- duration_cosine_similarity(g1,g2)
        
        results <- rbind(results, data.frame(Method="Tutorials-WholeClass",
                                             RawCourse=course,
                                             CosineSim=round(sim,2),
                                             DurationCosineSim = round(duration_sim, 2),
                                             stringsAsFactors = FALSE))
      }
    } else if (m == "Tutorials-RecitationOnly") {
      #SPECIAL CASE: only CALEP1 Class2 for CALEP2 Courses 4,5,8
      g1 <- aggregate_student_CALEP1_transition_network("CALEP1-COPUS-observations", "CALEP1_Tutorials_Class2", title=m)
      g1 <- normalize_vertex_names(g1)
      freq_data[[length(freq_data)+1]] <- extract_frequencies(g1, m, "CALEP1")
      
      tutorialsB_courses <- paste0("Tutorials_Course", c(4,5,8))
      
      for (course in tutorialsB_courses) {
        g2 <- aggregate_student_CALEP2_transition_network(course)
        g2 <- normalize_vertex_names(g2)
        
        freq_data[[length(freq_data)+1]] <- extract_frequencies(g2, m, "CALEP2", course)
        
        sim <- cosine_similarity(g1, g2)
        duration_sim <- duration_cosine_similarity(g1,g2)
        
        results <- rbind(results, data.frame(Method="Tutorials-RecitationOnly",
                                             RawCourse=course,
                                             CosineSim=round(sim,2),
                                             DurationCosineSim = round(duration_sim, 2),
                                             stringsAsFactors = FALSE))
      }
    }
  }
  
  freq_df <- do.call(rbind, freq_data)
  
  write.csv(freq_df, "data/student-code-frequencies.csv", row.names=FALSE)
  write.csv(results, "data/student-fidelity-values.csv", row.names = FALSE)
  
  return(list(results=results, freq=freq_df))
}

#compare for instructors
compare_methods <- function() { 
  methods <- c("SCALEUP","ISLE", "Tutorials-WholeClass", "Tutorials-RecitationOnly")
  
  results <- data.frame(Method=character(),
                        RawCourse=character(),
                        CosineSim=numeric(),
                        DurationCosineSim=numeric(),
                        stringsAsFactors = FALSE)
  freq_data <- list() 
  
  for (m in methods) {
    par(mfrow = c(3, 3), mar = c(1, 1, 2, 1))
    if (m %in% c("SCALEUP", "ISLE")) {
      g1 <- aggregate_CALEP1_transition_network("CALEP1-COPUS-observations", paste0("CALEP1_", m), title=m)
      g1 <- normalize_vertex_names(g1)
      
      freq_data[[length(freq_data)+1]] <- extract_frequencies(g1, m, "CALEP1")
      
      raw_files <- list.files("CALEP2-COPUS-observations", pattern = paste0("^", m, "_Course[0-9]+_Class1\\.csv$"))
      raw_files <- trimws(raw_files)  
      courses <- unique(gsub("_Class1\\.csv$", "", raw_files))  
      
      for (course in courses) {
        g2 <- aggregate_CALEP2_transition_network(course)
        g2 <- normalize_vertex_names(g2)
        
        freq_data[[length(freq_data)+1]] <- extract_frequencies(g2, m, "CALEP2", course)
        
        sim <- cosine_similarity(g1, g2)
        duration_sim <- duration_cosine_similarity(g1,g2)
        
        results <- rbind(results, data.frame(Method=m,
                                             RawCourse=course,
                                             CosineSim=round (sim,2),
                                             DurationCosineSim = round(duration_sim, 2),
                                             stringsAsFactors = FALSE))
        
      }
    } 
    else if (m == "Tutorials-WholeClass") {
      #SPECIAL CASE: aggregated CALEP1 Class1+Class2 for CALEP2 Courses 2,3,7,9
      g1 <- aggregate_CALEP1_transition_network("CALEP1-COPUS-observations", 
                                                c("CALEP1_Tutorials_Class1", "CALEP1_Tutorials_Class2"), title=m)
      g1 <- normalize_vertex_names(g1)
      
      freq_data[[length(freq_data)+1]] <- extract_frequencies(g1, m, "CALEP1")
      tutorialsA_courses <- paste0("Tutorials_Course", c(2,3,7,9))
      
      for (course in tutorialsA_courses) {
        g2 <- aggregate_CALEP2_transition_network(course)
        g2 <- normalize_vertex_names(g2)
        
        freq_data[[length(freq_data)+1]] <- extract_frequencies(g2, m, "CALEP2", course)
        
        sim <- cosine_similarity(g1, g2)
        duration_sim <- duration_cosine_similarity(g1,g2)
        
        results <- rbind(results, data.frame(Method="Tutorials-WholeClass",
                                             RawCourse=course,
                                             CosineSim=round(sim,2),
                                             DurationCosineSim = round(duration_sim, 2),
                                             stringsAsFactors = FALSE))
      }
    } else if (m == "Tutorials-RecitationOnly") {
      #SPECIAL CASE: only CALEP1 Class2 for CALEP2 Courses 4,5,8
      g1 <- aggregate_CALEP1_transition_network("CALEP1-COPUS-observations", "CALEP1_Tutorials_Class2", title=m)
      g1 <- normalize_vertex_names(g1)
      
      freq_data[[length(freq_data)+1]] <- extract_frequencies(g1, m, "CALEP1")
      
      tutorialsB_courses <- paste0("Tutorials_Course", c(4,5,8))
      
      for (course in tutorialsB_courses) {
        g2 <- aggregate_CALEP2_transition_network(course)
        g2 <- normalize_vertex_names(g2)
        
        freq_data[[length(freq_data)+1]] <- extract_frequencies(g2, m, "CALEP2", course)
        
        sim <- cosine_similarity(g1, g2)
        duration_sim <- duration_cosine_similarity(g1,g2)
        
        results <- rbind(results, data.frame(Method="Tutorials-RecitationOnly",
                                             RawCourse=course,
                                             CosineSim=round(sim,2),
                                             DurationCosineSim = round(duration_sim, 2),
                                             stringsAsFactors = FALSE))
      }
    }
  }
  
  freq_df <- do.call(rbind, freq_data)
  write.csv(freq_df, "data/inst-code-frequencies.csv", row.names=FALSE)
  write.csv(results, "data/inst-fidelity-values.csv", row.names = FALSE)
  return(list(results=results, freq=freq_df))
}

student_plot_box_with_overlay <- function(method_name) {
  # Critical component codes to bold for each method
  bold_codes_list <- list(
    ISLE      = c("OG", "Prd", "CG", "Ind"),
    SCALEUP   = c("OG", "AnQ_S", "SP"),
    "Tutorials-RecitationOnly" = c("WG"),
    "Tutorials-WholeClass" = c("WG") 
  )
  bold_codes <- bold_codes_list[[method_name]]
  
  calep2_sub <- filter(student_calep2_data, Method == method_name) %>% mutate(Label = method_name)
  calep1_sub <- filter(student_calep1_data, Method == method_name) %>% mutate(Label = method_name)
  method_label <- gsub("-", "-\n", method_name)
  ggplot() +
    # Instructors at large (CALEP2)
    geom_boxplot(
      data = calep2_sub,
      aes(x = Code, y = Frequency, fill = "Instructors at large"),
      outlier.shape = NA, alpha = 0.5, color="black"
    ) +
    geom_jitter(
      data = calep2_sub,
      aes(x = Code, y = Frequency),
      shape = 21, size = 1.5, width = 0.2, height = 0, fill="black", color = "black"
    ) +
    # High fidelity (CALEP1)
    geom_point(
      data = calep1_sub,
      aes(x = Code, y = Frequency, fill = "High fidelity"),
      shape = 23, size = 2
    ) +
    facet_wrap (
      ~ Label, 
      ncol=1, 
      strip.position = "right", 
      labeller = labeller(Label = function(x){
        x <- ifelse(x == "SCALEUP", "SCALE-UP", x)  # rename SCALEUP properly
        x <- gsub("Tutorials-", "Tutorials-\n", x)  # only break these long ones
        return(x) }) 
    )+  
    labs(
      #title = method_name,
      #x = "Code", y = "Frequency", fill = NULL
      x = NULL,  # remove "Code"
      y = NULL,  # remove y-axis label
      fill = NULL
    ) +
    scale_fill_manual(
      values = c("High fidelity" = "red", "Instructors at large" = "grey"),
      labels = c("High fidelity implementation", "Broader implementations")
    ) +
    scale_y_continuous(limits = c(0, 1)) +   # uniform scale
    theme_minimal(base_size = 10) +
    theme(
      strip.text.y = element_text(
        angle = 90, face = "bold", color = "black",
        size = 12, margin = margin(l = 6, r = 4)
      ),
      strip.background = element_rect(fill = "grey90", color = NA, linewidth = 0.3),
      axis.text.x = element_text(
        face = ifelse(levels(factor(student_calep2_data$Code)) %in% bold_codes,
                      "bold", "plain"),
        size = 9, angle = 45, margin = margin(t = 5)
      ),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
     
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.spacing = unit(0.5, "cm"),
      
      plot.title = element_blank(),
      legend.position = "none"
    )
}


plot_box_with_overlay <- function(method_name) {
  # Critical component codes to bold for each method
  bold_codes_list <- list(
    ISLE      = c("MG", "D/V", "CQ"),
    SCALEUP   = c("MG", "FUp", "PQ"),
    "Tutorials-RecitationOnly" = c("MG"),
    "Tutorials-WholeClass" = c("MG")  #critical components to be bolded in plot
  )
  bold_codes <- bold_codes_list[[method_name]]
  
  ggplot() +
    # Instructors at large (CALEP2)
    geom_boxplot(
      data = filter(calep2_data, Method == method_name),
      aes(x = Code, y = Frequency, fill = "Instructors at large"),
      outlier.shape = NA, alpha = 0.5, color="black"
    ) +
    geom_jitter(
      data = filter(calep2_data, Method == method_name),
      aes(x = Code, y = Frequency),
      shape = 21, size = 1.5, width = 0.2, height = 0, fill="black", color = "black"
    ) +
    # High fidelity (CALEP1)
    geom_point(
      data = filter(calep1_data, Method == method_name),
      aes(x = Code, y = Frequency, fill = "High fidelity"),
      shape = 23, size = 2
    ) +
    labs(
      #title = method_name,
      #x = "Code",
      x = NULL,
      y = "Proportion of Class Time", fill = NULL
    ) +
    scale_fill_manual(
      values = c("High fidelity" = "red", "Instructors at large" = "grey"),
      labels = c("High fidelity implementation", "Broader implementations")
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(
        face = ifelse(levels(factor(calep2_data$Code)) %in% bold_codes, "bold", "plain"), size =9, angle = 45, margin = margin(t = 5)
      ),
      axis.title.y = element_blank(),
      panel.grid.major.x = element_blank(),   
      panel.grid.minor.x = element_blank(),
      panel.spacing = unit(0.5, "cm"),
      legend.position = "none"
    )
}

make_box_plot <- function(data, y_var, y_label = NULL, show_y = TRUE, show_x = TRUE) {
  ggplot(data, aes(x = Method, y = !!sym(y_var))) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, fill = "grey") +
    geom_jitter(width = 0.2, height = 0, aes(color = Implementation), size = 1.5) +
    scale_color_manual(values = color_map) +
    scale_x_discrete(labels = c("SCALEUP" = "SCALE-UP", "ISLE" = "ISLE", 
                                "Tutorials-RecitationOnly" = "Tutorials-\nRecitationOnly", 
                                "Tutorials-WholeClass" = "Tutorials-\nWholeClass")) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = if(show_x) NULL else "", y = if(show_y) y_label else NULL, color = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "none",
      axis.text.x = if(show_x) element_text(angle = 45, hjust = 1, vjust = 1, size = 10, margin = margin(t = -1)) else element_blank(),
      axis.title.x = if(show_x) element_text() else element_blank(),
      axis.title.y = if(show_y) element_text(margin = margin(r = 10)) else element_blank(),
      axis.text.y = if(show_y) element_text() else element_blank(),
      axis.ticks.y = if(show_y) element_line() else element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}

make_scatter <- function(data, x_col, x_label, show_y = TRUE) {
  ggplot(data, aes(x = .data[[x_col]], y = HedgesG, color = Method)) +
    geom_point(size = 1.5, alpha = 0.8) +
    geom_linerange(aes(ymin = LowerCI, ymax = UpperCI), size = 0.4, alpha = 0.6) +
    scale_color_manual(values = method_colors, breaks = method_order, labels = method_labels) +
    scale_x_continuous(limits = c(0, 1)) +
    labs(x = x_label, y = if (show_y) "Effect Size (Hedges' g)" else NULL, color = "Method") +
    theme_minimal(base_size = 10) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
    theme(
      legend.position = "none",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      axis.text.y = if (show_y) element_text(size = 10) else element_blank(),
      axis.ticks.y = if (show_y) element_line(color = "black", linewidth = 0.3) else element_blank(),
      axis.text.x = element_text(size = 10),
      axis.ticks.x = element_line(color = "black", linewidth = 0.3), 
      axis.line.x = element_line(color = "black", linewidth = 0.3),
      axis.line.y = element_line(color = "black", linewidth = 0.3),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}

#Plot proportion of time spent on all COPUS codes-------------------------------
if (!dir.exists("data")) dir.create("data")

student_out <- compare_student_methods()
student_results <- student_out$results
student_freq_data <- student_out$freq

student_calep2_data <- subset(student_freq_data, Source == "CALEP2")
student_calep1_data <- subset(student_freq_data, Source == "CALEP1")


out <- compare_methods()

results <- out$results
freq_data <- out$freq   # <- this is the dataframe you’ll plot

calep2_data <- subset(freq_data, Source == "CALEP2")
calep1_data <- subset(freq_data, Source == "CALEP1")


p1<- student_plot_box_with_overlay("SCALEUP")
p2<- student_plot_box_with_overlay("ISLE")
p3<- student_plot_box_with_overlay("Tutorials-RecitationOnly")
p4<- student_plot_box_with_overlay("Tutorials-WholeClass")

p5<- plot_box_with_overlay("SCALEUP")
p6<- plot_box_with_overlay("ISLE")
p7<- plot_box_with_overlay("Tutorials-RecitationOnly")
p8<- plot_box_with_overlay("Tutorials-WholeClass")

p4 <- p4 + labs(x = "Student COPUS codes") +
  theme(axis.title.x = element_text(size = 10, margin = margin(t = 10)))
p8 <- p8 + labs(x = "Instructor COPUS codes") +
  theme(axis.title.x = element_text(size = 10, margin = margin(t = 10)))

legend_plot <- plot_box_with_overlay("SCALEUP") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 10))
shared_legend <- get_legend(
  legend_plot
)

combined <- (p5 | p1) /
  (p6 | p2) /
  (p7 | p3) /
  (p8 | p4)

row1 <- plot_grid(p5, plot_spacer()+ theme_void(), p1, ncol = 3, rel_widths = c(1.2, 0.03, 1.2))
row2 <- plot_grid(p6, plot_spacer()+ theme_void(), p2, ncol = 3, rel_widths = c(1.2, 0.03, 1.2))
row3 <- plot_grid(p7, plot_spacer()+ theme_void(), p3, ncol = 3, rel_widths = c(1.2, 0.03, 1.2))
row4 <- plot_grid(p8, plot_spacer()+ theme_void(), p4, ncol = 3, rel_widths = c(1.2, 0.03, 1.2))

combined <- plot_grid(row1, row2, row3, row4, ncol = 1)
combined <- plot_grid(
  ggdraw() + draw_label("Proportion of Class Time", angle = 90, size = 12, vjust = 0.5, x = 0.2, hjust = 0),
  combined,
  ncol = 2,
  rel_widths = c(0.05, 1)
)

final_plot <- plot_grid(
  combined,
  shared_legend,
  ncol = 1,
  rel_heights = c(1, 0.08)  # adjust legend space
)

ggsave("proportion-COPUS-codes.pdf", plot = final_plot, width = 7, height = 8, units = "in")


#Plot Duration and Transition Cosine Similarity---------------------------------

implementation_map <- data.frame(
  RawCourse = c("ISLE_Course2","ISLE_Course4","ISLE_Course5","ISLE_Course6","SCALEUP_Course1","SCALEUP_Course2","SCALEUP_Course3",
                "SCALEUP_Course4","SCALEUP_Course5","SCALEUP_Course6","SCALEUP_Course7", "Tutorials_Course2", "Tutorials_Course3", "Tutorials_Course4", "Tutorials_Course5", "Tutorials_Course7",
                "Tutorials_Course8", "Tutorials_Course9"),
  Implementation = c("Whole Class", "Whole Class","Whole Class","Whole Class", "Whole Class","Whole Class","Whole Class","Whole Class","Whole Class","Whole Class","Whole Class", "Whole Class","Whole Class",
                     "Recitation Only", "Recitation Only","Whole Class", "Recitation Only", "Whole Class"),
  stringsAsFactors = FALSE
)

#if you want points to be colored by implementation
color_map<- c( 
  "Lab Only" = "#FFD700",
  "Recitation Only" = "black",
  "Whole Class" = "black" )

inst_results <- read.csv("data/inst-fidelity-values.csv")
inst_results <- merge(inst_results, implementation_map, by = "RawCourse", all.x = TRUE)

student_results <- read.csv("data/student-fidelity-values.csv")
student_results <- merge(student_results, implementation_map, by = "RawCourse", all.x = TRUE)

inst_results$Method <- factor(inst_results$Method,
                              levels = c("SCALEUP", "ISLE", "Tutorials-RecitationOnly", "Tutorials-WholeClass"))

student_results$Method <- factor(student_results$Method,
                                 levels = c("SCALEUP", "ISLE", "Tutorials-RecitationOnly", "Tutorials-WholeClass"))

# Row 1: Instructor vs Student (Duration)
p9  <- make_box_plot(inst_results, "DurationCosineSim", y_label = "Cosine Similarity (Duration)", show_y = TRUE, show_x = FALSE)
p10 <- make_box_plot(student_results, "DurationCosineSim", show_y = FALSE, show_x = FALSE)

# Row 2: Instructor vs Student (Transition)
p11 <- make_box_plot(inst_results, "CosineSim", y_label = "Cosine Similarity (Transition)", show_y = TRUE, show_x = TRUE)
p12 <- make_box_plot(student_results, "CosineSim", show_y = FALSE, show_x = TRUE)

column_labels <- plot_grid(
  ggdraw() + draw_label("Instructor", hjust = 0.5, size = 10),
  ggdraw() + draw_label("Student", hjust = 0.5, size = 10),
  ncol = 2
)

row1 <- plot_grid(p9, plot_spacer() + theme_void(), p10, ncol = 3, rel_widths = c(1, 0.2, 1))
row2 <- plot_grid(p11, plot_spacer() + theme_void(), p12, ncol = 3, rel_widths = c(1, 0.2, 1))

row1_labeled <- plot_grid(ggdraw() + draw_label("(a)", angle = 0, size = 12, x = 0.02, hjust = 0),
                          row1, ncol = 2, rel_widths = c(0.05, 1))
row2_labeled <- plot_grid(ggdraw() + draw_label("(b)", angle = 0, size = 12, x = 0.02, y=0.65, hjust = 0),
                          row2, ncol = 2, rel_widths = c(0.05, 1))

column_labels <- plot_grid(ggdraw() + draw_label("Instructor", hjust = 0.5, size = 10),
                           ggdraw() + draw_label("Student", hjust = 0.5, size = 10), ncol = 2)

combined <- plot_grid(column_labels, row1_labeled, row2_labeled, ncol = 1, rel_heights = c(0.07, 0.93, 1.4), align = "hv")
ggsave("duration-transition-fidelity.pdf", plot = combined, width = 7, height = 5)


#Plot of Effect Size against Fidelity-------------------------------------------  

# Define custom colors
method_colors <- c(
  "SCALEUP" = "#004D40",
  "ISLE" = "#FFC107",
  "Tutorials-RecitationOnly" = "#5D3A9B",
  "Tutorials-WholeClass" = "#D81B60"
)


effect_sizes <- read.csv("data/EffectSizes.csv")
effect_sizes <- effect_sizes %>%
  mutate(RawCourse = paste0(Method, "_Course", Instructor))
effect_sizes <- effect_sizes[, !names(effect_sizes) %in% "Method"]

inst_merged <- merge(inst_results, effect_sizes, by = "RawCourse", all.x = TRUE)
student_merged <- merge(student_results, effect_sizes, by = "RawCourse", all.x = TRUE)

method_order <- c("SCALEUP", "ISLE", "Tutorials-RecitationOnly", "Tutorials-WholeClass")
method_labels <- c("SCALEUP" = "SCALE-UP", 
                   "ISLE" = "ISLE",
                   "Tutorials-RecitationOnly" = "Tutorials-RecitationOnly",
                   "Tutorials-WholeClass" = "Tutorials-WholeClass")

inst_merged$Method <- factor(inst_merged$Method, levels = method_order)
student_merged$Method <- factor(student_merged$Method, levels = method_order)

p_inst_dur  <- make_scatter(inst_merged, "DurationCosineSim", "Cosine Similarity (Duration)", show_y = TRUE)
p_stud_dur  <- make_scatter(student_merged, "DurationCosineSim", "Cosine Similarity (Duration)", show_y = FALSE)
p_inst_tran <- make_scatter(inst_merged, "CosineSim", "Cosine Similarity (Transition)", show_y = TRUE)
p_stud_tran <- make_scatter(student_merged, "CosineSim", "Cosine Similarity (Transition)", show_y = FALSE)

row1 <- plot_grid(p_inst_dur, plot_spacer()+ theme_void(), p_stud_dur, ncol = 3, rel_widths = c(1, 0.2, 1))
row2 <- plot_grid(p_inst_tran, plot_spacer()+ theme_void(), p_stud_tran, ncol = 3, rel_widths = c(1, 0.2, 1))

row1_labeled <- plot_grid(
  ggdraw() + draw_label("(a)", size = 12, x = 0.02, hjust = 0),row1,ncol = 2, rel_widths = c(0.05, 1)
)

row2_labeled <- plot_grid(
  ggdraw() + draw_label("(b)", size = 12, x = 0.02, hjust = 0),row2,ncol = 2,rel_widths = c(0.05, 1)
)

column_labels <- plot_grid(
  ggdraw() + draw_label("Instructor", hjust = 0.5, size = 10),
  ggdraw() + draw_label("Student", hjust = 0.5, size = 10),
  ncol = 2
)

combined <- plot_grid(
  column_labels,
  row1_labeled,
  row2_labeled,
  ncol = 1,
  rel_heights = c(0.07, 0.9, 1)
)
legend_plot <- make_scatter(inst_merged, "DurationCosineSim", "Cosine Similarity (Duration)", show_y = TRUE)+
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 10))
shared_legend <- get_legend(legend_plot)

final_plot <- plot_grid(
  combined,
  shared_legend,
  ncol = 1,
  rel_heights = c(1, 0.08)
)
ggsave("EffectSize_fidelity.pdf", final_plot, width = 7, height = 5)
