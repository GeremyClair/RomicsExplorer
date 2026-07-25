# ============================================================================
# RomicsExplorer - Interactive Shiny App for Omics Data Exploration
# ============================================================================
# Date: March 18 2026
# Purpose: Web-based visualization and analysis of RomicsProcessor objects

# ============================================================================
# 1. CONFIGURATION & LIBRARIES
# ============================================================================

MAX_UPLOAD_MB <- 1000
options(shiny.maxRequestSize = MAX_UPLOAD_MB * 1024^2)
options(warn = -1)

# Load required packages
library(shiny)
library(shinydashboard)
library(shinyjs)
library(RomicsProcessor)
library(DT)
library(ggplot2)
library(plotly)
library(viridis)
library(cowplot)

# Load optional packages
has_ggdendro <- requireNamespace("ggdendro", quietly = TRUE)
if (has_ggdendro) library(ggdendro)

# ============================================================================
# 2. UTILITY FUNCTIONS - Factor Handling
# ============================================================================

clean_factor_names <- function(fn) fn[!tolower(fn) %in% c("romics_colors","colors_romics")]

get_clean_factors <- function(obj) {
  fns <- try(RomicsProcessor::romicsFactorNames(obj), silent=TRUE)
  if (inherits(fns,"try-error") || !length(fns)) return(character(0))
  clean_factor_names(fns)
}
# ============================================================================
# 3. UTILITY FUNCTIONS - Statistics & Data Extraction
# ============================================================================

extract_feature_statistics <- function(obj, feat, pval_type="all") {
  if (is.null(obj$statistics)) return(NULL)
  if (is.null(feat) || feat == "") return(NULL)

  idx <- which(rownames(obj$statistics) == feat)
  if (!length(idx)) return(NULL)

  stat_row <- obj$statistics[idx, , drop=FALSE]
  all_cols <- colnames(stat_row)

  if (length(all_cols) == 0) return(NULL)

  # Filter p-value columns based on user selection - match both common patterns
  if (pval_type == "all") {
    p_cols <- all_cols[grepl("_p$|_padj$", all_cols) & grepl("Ttest|GLMbinomial|test|comparison", all_cols, ignore.case=TRUE)]
    if (length(p_cols) == 0) {
      p_cols <- all_cols[grepl("_Ttest_p$|_Ttest_padj$|_GLMbinomial_p$|_GLMbinomial_padj$|_p$|_padj$", all_cols)]
    }
  } else if (pval_type == "p") {
    p_cols <- all_cols[grepl("_p$", all_cols) & !grepl("adj", all_cols)]
    if (length(p_cols) == 0) {
      p_cols <- all_cols[grepl("_Ttest_p$|_GLMbinomial_p$", all_cols)]
    }
  } else if (pval_type == "padj") {
    p_cols <- all_cols[grepl("_padj$", all_cols)]
    if (length(p_cols) == 0) {
      p_cols <- all_cols[grepl("_Ttest_padj$|_GLMbinomial_padj$", all_cols)]
    }
  }

  if (length(p_cols) == 0) return(NULL)

  # Extract p-values safely, handling NAs
  p_vals <- suppressWarnings(as.numeric(stat_row[1, p_cols]))
  
  # Create data frame with full column names
  df <- data.frame(
    Test = p_cols,
    PValue = p_vals,
    stringsAsFactors = FALSE
  )

  # Remove rows with NA p-values
  df <- df[!is.na(df$PValue), ]
  
  if (nrow(df) == 0) return(NULL)

  # Sort by p-value (smallest first)
  df <- df[order(df$PValue), ]
  rownames(df) <- NULL
  df
}
# ============================================================================
# 4. UTILITY FUNCTIONS - Embedding Detection & Dimensions
# ============================================================================

check_embeddings_available <- function(obj) {
  # Check which embedding types (PCA, UMAP, t-SNE) are available in the object
  result <- list(pca=FALSE, umap=FALSE, tsne=FALSE)

  if (is.null(obj$embeddings)) return(result)

  # Check if embeddings is a list (new format)
  if (is.list(obj$embeddings) && !is.data.frame(obj$embeddings) && !is.matrix(obj$embeddings)) {
    embedding_names <- tolower(names(obj$embeddings))
    result$pca <- any(grepl("pca", embedding_names))
    result$umap <- any(grepl("umap", embedding_names))
    result$tsne <- any(grepl("tsne", embedding_names))
    return(result)
  }

  # Check if embeddings is data.frame or matrix (old format)
  if (is.data.frame(obj$embeddings) || is.matrix(obj$embeddings)) {
    all_names <- c(tolower(rownames(obj$embeddings)), tolower(colnames(obj$embeddings)))
    result$pca <- any(grepl("pca", all_names))
    result$umap <- any(grepl("umap", all_names))
    result$tsne <- any(grepl("tsne", all_names))
  }

  result
}
get_embedding_dimensions <- function(obj, type) {
  if (is.null(obj$embeddings)) return(0)

  type_lower <- tolower(type)

  # If embeddings is a list
  if (is.list(obj$embeddings) && !is.data.frame(obj$embeddings) && !is.matrix(obj$embeddings)) {
    matching_names <- names(obj$embeddings)[grepl(type_lower, tolower(names(obj$embeddings)))]
    if (length(matching_names) > 0) {
      embedding <- obj$embeddings[[matching_names[1]]]
      if (is.data.frame(embedding) || is.matrix(embedding)) {
        return(ncol(embedding))
      }
    }
    return(0)
  }

  # If embeddings is data.frame or matrix
  if (is.data.frame(obj$embeddings) || is.matrix(obj$embeddings)) {
    row_names <- rownames(obj$embeddings)
    col_names <- colnames(obj$embeddings)

    n_matching_rows <- length(row_names[grepl(type_lower, tolower(row_names))])
    n_matching_cols <- length(col_names[grepl(type_lower, tolower(col_names))])

    return(max(n_matching_rows, n_matching_cols))
  }

  0
}
# ============================================================================
# 5. UTILITY FUNCTIONS - Statistical Comparison Detection
# ============================================================================

extract_comparisons_for_factor <- function(obj, fn) {
  # Find all pairwise statistical comparisons available for a given factor
  if (is.null(obj$statistics)) return(NULL)

  stat_cols <- colnames(obj$statistics)

  # Extract factor levels
  factor_data <- try(RomicsProcessor::romicsExtractFactor(obj, factor=fn), silent=TRUE)
  if (inherits(factor_data, "try-error")) return(NULL)

  factor_levels <- levels(factor_data)
  if (length(factor_levels) < 2) return(NULL)

  # Find all pairwise comparisons
  comparisons <- list()
  for (i in 1:(length(factor_levels)-1)) {
    for (j in (i+1):length(factor_levels)) {
      level1 <- factor_levels[i]
      level2 <- factor_levels[j]

      # Look for comparison columns (try both directions)
      comparison1 <- paste0(level1, "_vs_", level2)
      comparison2 <- paste0(level2, "_vs_", level1)

      # Find adjusted p-value or p-value column
      p_col <- stat_cols[grepl(paste0(comparison1, "_Ttest_padj"), stat_cols)]
      if (!length(p_col)) p_col <- stat_cols[grepl(paste0(comparison2, "_Ttest_padj"), stat_cols)]
      if (!length(p_col)) p_col <- stat_cols[grepl(paste0(comparison1, "_Ttest_p$"), stat_cols)]
      if (!length(p_col)) p_col <- stat_cols[grepl(paste0(comparison2, "_Ttest_p$"), stat_cols)]

      if (length(p_col) > 0) {
        clean_name <- gsub("_Ttest_padj$|_Ttest_p$", "", p_col[1])
        comparisons[[clean_name]] <- clean_name
      }
    }
  }

  if (!length(comparisons)) return(NULL)
  comparisons
}

extract_all_comparisons <- function(obj) {
  # Get all statistical comparison columns (across all factors)
  if (is.null(obj$statistics)) return(NULL)
  sc <- colnames(obj$statistics)
  pc <- sc[grepl("_Ttest_padj",sc)]
  if (!length(pc)) pc <- sc[grepl("_Ttest_p$",sc)]
  if (!length(pc)) return(NULL)
  ac <- gsub("_Ttest_padj$|_Ttest_p$","",pc)
  setNames(as.list(ac),ac)
}

check_pvalue_types <- function(obj, cn) {
  # Check if comparison has adjusted p-value and/or raw p-value
  sc <- colnames(obj$statistics)
  list(
    p = any(grepl(paste0(cn,"_Ttest_p$"),sc)),
    padj = any(grepl(paste0(cn,"_Ttest_padj"),sc))
  )
}

extract_volcano_data <- function(obj, cn, ptype="padj") {
  # Extract log2 fold change and p-value data for volcano plot
  if (is.null(obj$statistics)) return(NULL)
  sc <- colnames(obj$statistics)

  # Find fold change column (search for log2(A/B) pattern)
  fcp <- paste0("log\\(",gsub("_vs_","/",cn),"\\)")
  fc_col <- sc[grepl(fcp,sc)]

  # Find appropriate p-value column
  pv_col <- if (ptype=="padj") {
    sc[grepl(paste0(cn,"_Ttest_padj"),sc)]
  } else {
    sc[grepl(paste0(cn,"_Ttest_p$"),sc)]
  }

  if (!length(fc_col) || !length(pv_col)) return(NULL)

  # Build volcano data frame
  vd <- data.frame(
    Feature = rownames(obj$statistics),
    log_FC = as.numeric(obj$statistics[,fc_col[1]]),
    p_value = as.numeric(obj$statistics[,pv_col[1]]),
    stringsAsFactors = FALSE
  )

  # Remove missing values and calculate -log10(p)
  vd <- vd[!is.na(vd$log_FC) & !is.na(vd$p_value),]
  vd$neg_log10_p <- -log10(pmax(vd$p_value, 1e-300))
  vd
}

get_enrichment_stat_cols <- function(obj) {
  # Get statistical columns suitable for enrichment filtering (exclude metadata columns)
  if (is.null(obj$statistics)) return(character(0))
  sc <- colnames(obj$statistics)
  sc[!grepl("_mean$|_sd$|^Z_scores_|^Feature$|_percentage_completeness$",sc)]
}

apply_stat_filter <- function(values, fe) {
  fe <- trimws(fe); if (is.null(fe)||fe=="") return(rep(TRUE,length(values)))
  op <- gsub("[0-9\\.\\-eE+]","",fe); val <- as.numeric(gsub("[^0-9\\.\\-eE+]","",fe))
  if (is.na(val)) return(rep(TRUE,length(values)))
  switch(op, "<"=!is.na(values)&values<val, "<="=!is.na(values)&values<=val,
         ">"=!is.na(values)&values>val, ">="=!is.na(values)&values>=val,
         "=="=!is.na(values)&values==val, rep(TRUE,length(values)))
}
get_enrichment_ids_for_features <- function(obj, features, id_col) {
  if (is.null(obj$IDs)||!id_col %in% colnames(obj$IDs)) return(character(0))
  idx <- match(features, rownames(obj$data)); idx <- idx[!is.na(idx)]
  ids <- obj$IDs[[id_col]][idx]; unique(ids[!is.na(ids)&ids!=""])
}
ontology_term_preview <- function(tbl, pat, fallbacks) {
  if (is.null(tbl)) return(HTML("<div class='ontology-term-list' style='color:#999;'><i>Not loaded.</i></div>"))
  dc <- NULL; for (p in c(pat,fallbacks)) { m <- colnames(tbl)[grepl(p,colnames(tbl),ignore.case=TRUE)]; if (length(m)>0) {dc <- m[1]; break} }
  if (is.null(dc)) return(HTML(paste0("<div class='ontology-term-list' style='color:#c0392b;'>Column not found. Available: ",paste(colnames(tbl),collapse=", "),"</div>")))
  terms <- unique(tbl[[dc]]); terms <- terms[!is.na(terms)&terms!=""]; terms <- head(terms,10)
  items <- paste0("<span class='term-number'>",seq_along(terms),".</span> ",terms,collapse="<br/>")
  HTML(paste0("<div class='ontology-term-list'><b>",dc,"</b> (",length(unique(tbl[[dc]]))," unique)<br/><br/>",items,"</div>"))
}
load_ontology_file <- function(fp, fn) {
  ext <- tolower(tools::file_ext(fn))
  if (ext=="csv") tbl <- read.csv(fp,stringsAsFactors=FALSE)
  else if (ext=="tsv") tbl <- read.delim(fp,stringsAsFactors=FALSE)
  else if (ext=="rds") tbl <- readRDS(fp)
  else if (ext %in% c("rda","rdata")) {
    env <- new.env(parent=emptyenv()); load(fp,envir=env); objs <- as.list(env); tbl <- NULL
    for (nm in names(objs)) if (is.data.frame(objs[[nm]])) {tbl <- objs[[nm]]; break}
    if (is.null(tbl)) tbl <- objs[[1]]
  } else stop("Unsupported: ",ext)
  if (!is.data.frame(tbl)) tbl <- as.data.frame(tbl); tbl
}
build_feature_choices <- function(obj, include_all=TRUE, include_anova=TRUE, include_clusters=TRUE, heatmap_clusters_data=NULL) {
  choices <- list()
  if (include_all) choices[["All Features"]] <- "ALL_FEATURES"
  if (include_anova && !is.null(obj$statistics)) {
    sc <- colnames(obj$statistics); ap <- sc[grepl("^ANOVA_",sc)&grepl("_p$|_padj$",sc)]
    if (length(ap)>0) for (ac in ap) choices[[paste0("== ",ac," < 0.05 ==")]] <- paste0("ANOVA_FILTER_",ac)
  }
  if (include_clusters && !is.null(heatmap_clusters_data) && !is.null(heatmap_clusters_data$clusters))
    for (i in seq_along(heatmap_clusters_data$clusters))
      choices[[paste0("== Cluster ",i," (",length(heatmap_clusters_data$clusters[[i]])," feat) ==")]] <- paste0("CLUSTER_",i)
  fi <- rownames(obj$data); for (f in fi) choices[[f]] <- f; choices
}
resolve_dropdown_selection <- function(sel, obj, heatmap_clusters_data=NULL) {
  if (sel=="ALL_FEATURES") return(rownames(obj$data))
  if (grepl("^ANOVA_FILTER_",sel)) {
    ac <- gsub("^ANOVA_FILTER_","",sel)
    if (!is.null(obj$statistics)&&ac %in% colnames(obj$statistics)) {
      vals <- as.numeric(obj$statistics[,ac]); return(rownames(obj$statistics)[!is.na(vals)&vals<0.05])
    }; return(character(0))
  }
  if (grepl("^CLUSTER_",sel)) {
    cn <- as.numeric(gsub("CLUSTER_","",sel))
    if (!is.null(heatmap_clusters_data)&&cn<=length(heatmap_clusters_data$clusters)) return(heatmap_clusters_data$clusters[[cn]])
    return(character(0))
  }
  sel
}

# ============================================================================
# 6. UTILITY FUNCTIONS - Bulk Ontology Loading
# ============================================================================

load_bulk_ontology_tables <- function(filepath) {
  # Load all ontology tables from a single R object
  # Tries common naming patterns for uniprot, GO, KEGG, and Reactome tables

  result <- list(
    uniprot = NULL,
    go = NULL,
    kegg = NULL,
    reactome = NULL,
    loaded_names = c()
  )

  tryCatch({
    # Load the R object into a temporary environment
    temp_env <- new.env(parent=emptyenv())
    load(filepath, envir=temp_env)

    # Get all object names
    all_names <- ls(envir=temp_env)

    # Common naming patterns for each table type
    uniprot_patterns <- c("uniprot_table", "uniprot", "UniProtTable", "UniProt", "UniProtData")
    go_patterns <- c("go_table", "GO", "go_data", "GO_table", "UniProtTable_GO", "uniprot_table_go")
    kegg_patterns <- c("kegg_table", "KEGG", "kegg_data", "KEGG_table", "UniProtTable_KEGG", "uniprot_table_kegg")
    reactome_patterns <- c("reactome_table", "Reactome", "reactome_data", "REACTOME_table", "UniProtTable_Reactome", "uniprot_table_reactome")

    # Try to find and load each table type
    for (pattern in uniprot_patterns) {
      if (pattern %in% all_names) {
        obj <- get(pattern, envir=temp_env)
        if (is.data.frame(obj)) {
          result$uniprot <- obj
          result$loaded_names <- c(result$loaded_names, "UniProt")
          break
        }
      }
    }

    for (pattern in go_patterns) {
      if (pattern %in% all_names) {
        obj <- get(pattern, envir=temp_env)
        if (is.data.frame(obj)) {
          result$go <- obj
          result$loaded_names <- c(result$loaded_names, "GO")
          break
        }
      }
    }

    for (pattern in kegg_patterns) {
      if (pattern %in% all_names) {
        obj <- get(pattern, envir=temp_env)
        if (is.data.frame(obj)) {
          result$kegg <- obj
          result$loaded_names <- c(result$loaded_names, "KEGG")
          break
        }
      }
    }

    for (pattern in reactome_patterns) {
      if (pattern %in% all_names) {
        obj <- get(pattern, envir=temp_env)
        if (is.data.frame(obj)) {
          result$reactome <- obj
          result$loaded_names <- c(result$loaded_names, "Reactome")
          break
        }
      }
    }

  }, error=function(e) {
    stop(paste("Error loading ontology file:", e$message))
  })

  result
}

# ============================================================================
# 7. HEATMAP GENERATION FUNCTION
# ============================================================================
# Creates complex heatmap with dendrograms, clustering, and completeness info

heatmapFeatures <- function(romics_object, factor="main", scale_feature=TRUE, feature_list=NULL,
  filter_by_stat_column=NULL, stat_column_filter=NULL, viridis_option="viridis",
  show_completeness=TRUE, show_dendrogram=TRUE, show_feature_names=TRUE,
  show_clusters=FALSE, show_cluster_legend=FALSE, n_clusters=NULL, clustering_method="ward.D") {
  obj <- romics_object
  if (factor!="main") tryCatch({obj <- RomicsProcessor::romicsChangeFactor(obj,main_factor=factor)},error=function(e) NULL)
  obj$statistics <- NULL; obj <- RomicsProcessor::romicsMean(obj); obj <- RomicsProcessor::romicsPercentComplete(obj)
  sdf <- obj$statistics
  if (!is.null(feature_list)&&length(feature_list)>0) {vf <- feature_list[feature_list %in% rownames(sdf)]; if (length(vf)>0) sdf <- sdf[vf,,drop=FALSE]}
  if (!is.null(filter_by_stat_column)&&!is.null(stat_column_filter)&&filter_by_stat_column!="none"&&stat_column_filter!="") {
    if (filter_by_stat_column %in% colnames(romics_object$statistics)) {
      vals <- as.numeric(romics_object$statistics[,filter_by_stat_column]); names(vals) <- rownames(romics_object$statistics)
      fe <- trimws(stat_column_filter); op <- gsub("[0-9\\.\\-eE+]","",fe); val <- as.numeric(gsub("[^0-9\\.\\-eE+]","",fe))
      if (!is.na(val)) {keep <- switch(op,"<"=!is.na(vals)&vals<val,"<="=!is.na(vals)&vals<=val,">"=!is.na(vals)&vals>val,">="=!is.na(vals)&vals>=val,rep(TRUE,length(vals)))
        sdf <- sdf[rownames(sdf) %in% names(vals)[keep],,drop=FALSE]}
    }
  }
  if (nrow(sdf)==0) stop("No features remaining")
  df <- sdf; df$Feature <- rownames(df); rownames(df) <- NULL
  df_mean <- df[,grepl("_mean$|^Feature$",colnames(df)),drop=FALSE]; colnames(df_mean) <- gsub("_mean$","",colnames(df_mean))
  df_comp <- df[,grepl("_percentage_completeness$|^Feature$",colnames(df)),drop=FALSE]; colnames(df_comp) <- gsub("_percentage_completeness$","",colnames(df_comp))
  gc <- setdiff(colnames(df_mean),"Feature")
  if (scale_feature&&length(gc)>1) {mm <- as.matrix(df_mean[,gc]); mm <- t(scale(t(mm))); mm[is.nan(mm)] <- 0; df_mean[,gc] <- mm}
  hd <- data.frame(); for (g in gc) {temp <- data.frame(Feature=df$Feature,Group=g,Mean=df_mean[[g]],stringsAsFactors=FALSE)
    temp$Completeness <- if (show_completeness&&g %in% colnames(df_comp)) df_comp[[g]] else NA; hd <- rbind(hd,temp)}
  mc <- as.matrix(df_mean[,gc]); rownames(mc) <- df$Feature; mi <- mc
  for (ci in 1:ncol(mi)) {cd <- mi[,ci]; vd <- cd[!is.na(cd)]; if (length(vd)>0) {mv <- min(vd,na.rm=TRUE); sv <- sd(vd,na.rm=TRUE)
    iv <- if (!is.na(sv)&&sv>0) mv-(2*sv) else mv*0.5; mi[is.na(mi[,ci]),ci] <- iv}}
  clust <- tryCatch(hclust(dist(mi),method=clustering_method),error=function(e) tryCatch(hclust(dist(mc),method=clustering_method),error=function(e2) NULL))
  if (!is.null(clust)) {fo <- clust$labels[clust$order]; dobj <- clust} else {fo <- unique(hd$Feature); dobj <- NULL}
  hd$Feature <- factor(hd$Feature,levels=fo); nf <- length(fo); cl <- if (scale_feature) "Scaled Mean" else "Mean"; pt <- paste("Heatmap -",factor)
  ph <- ggplot(hd,aes(x=Group,y=Feature,fill=Mean))+geom_tile()+scale_fill_viridis_c(option=viridis_option,name=cl)+scale_y_discrete(expand=c(0,0))+theme_minimal()+
    theme(axis.text.x=element_text(angle=45,hjust=1,size=10),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title=element_blank(),panel.grid=element_blank(),legend.position="none",plot.margin=margin(5,5,5,0))
  leg <- cowplot::ggdraw(cowplot::get_legend(ph+theme(legend.position="left")+guides(fill=guide_colorbar(title=cl,title.position="top"))))
  pfn <- NULL; if (show_feature_names) {fnd <- data.frame(Feature=factor(fo,levels=fo),y=seq_along(fo),label=fo)
    pfn <- ggplot(fnd,aes(x=0,y=y,label=label))+geom_text(hjust=0,size=max(1.5,min(3,80/nf)))+scale_y_continuous(limits=c(0.5,nf+0.5),expand=c(0,0))+theme_void()+theme(plot.margin=margin(5,0,5,5))}
  pcl <- NULL; pcleg <- NULL
  if (show_clusters&&!is.null(dobj)&&!is.null(n_clusters)&&n_clusters>1) {
    clusters <- cutree(dobj,k=n_clusters); cdf <- data.frame(Feature=names(clusters),Cluster=clusters,stringsAsFactors=FALSE)
    cdf <- cdf[match(fo,cdf$Feature),]; assign("heatmapFeatureClust",list(clusters=split(cdf$Feature,cdf$Cluster),cluster_assignments=cdf,dendrogram=dobj,n_clusters=n_clusters,method=clustering_method),envir=.GlobalEnv)
    ccol <- hcl.colors(n_clusters,palette="Dark3"); cpd <- data.frame(Feature=factor(fo,levels=fo),Cluster=cdf$Cluster,y=seq_along(fo))
    pcl <- ggplot(cpd,aes(x=1,y=y,fill=factor(Cluster)))+geom_tile(width=1,height=1)+scale_fill_manual(values=ccol,name="Cluster",guide=guide_legend(ncol=1))+
      scale_y_continuous(breaks=seq_along(fo),limits=c(0.5,nf+0.5),expand=c(0,0))+theme_void()+theme(legend.position="none",plot.margin=margin(5,0,5,0))
    if (show_cluster_legend) pcleg <- cowplot::ggdraw(cowplot::get_legend(pcl+theme(legend.position="right")))
  }
  assemble <- function(pd=NULL) {
    pe <- list(leg); rw <- c(0.15); if (!is.null(pd)) {pe <- c(pe,list(pd)); rw <- c(rw,0.2)}
    if (show_clusters&&!is.null(pcl)) {pe <- c(pe,list(pcl)); rw <- c(rw,0.06)}
    pe <- c(pe,list(ph)); rw <- c(rw,1)
    if (show_feature_names&&!is.null(pfn)) {pe <- c(pe,list(pfn)); rw <- c(rw,0.25)}
    if (show_clusters&&show_cluster_legend&&!is.null(pcleg)) {pe <- c(pe,list(pcleg)); rw <- c(rw,0.08)}
    pc <- do.call(cowplot::plot_grid,c(pe,list(nrow=1,rel_widths=rw,align="h",axis="b")))
    cowplot::plot_grid(cowplot::ggdraw()+cowplot::draw_text(pt,x=.5,y=.5,hjust=.5,vjust=.5,size=14,fontface="bold"),pc,nrow=2,rel_heights=c(0.08,1))
  }
  if (!is.null(dobj)&&show_dendrogram&&has_ggdendro) {
    tryCatch({dd <- ggdendro::dendro_data(dobj,type="rectangle")
      pd <- ggplot(ggdendro::segment(dd))+geom_segment(aes(x=y,y=x,xend=yend,yend=xend),size=0.5,color="gray40")+
        scale_y_continuous(breaks=seq_along(fo),labels=fo,limits=c(0.5,nf+0.5),expand=c(0,0))+scale_x_reverse(expand=c(0,0))+theme_void()+theme(axis.text.y=element_blank(),plot.margin=margin(0,5,5,0))
      return(assemble(pd))},error=function(e) return(assemble()))
  }
  assemble()
}

# ============================================================================
# Embedding Plot Helper Functions
# ============================================================================

render_2d_embedding <- function(obj, plot_func, input, color_factor, color_feature, embedding_dims, output) {
  # Get and validate component indices
  x_comp <- min(input$x_comp, embedding_dims)
  y_comp <- min(input$y_comp, embedding_dims)

  # Ensure x and y are different
  if (x_comp == y_comp) {
    y_comp <- min(y_comp + 1, embedding_dims)
  }

  # Build plot parameters
  plot_params <- list(
    romics_object=obj,
    Xcomp=x_comp,
    Ycomp=y_comp,
    size=input$grouping_size,
    alpha=input$grouping_alpha
  )

  # Add type for PCA plots
  if (input$grouping_type == "pca") {
    plot_params$plotType <- "individual"
  }

  # Add color parameters
  if (!is.null(color_factor)) {
    plot_params$factor_name <- color_factor
  }
  if (!is.null(color_feature)) {
    plot_params$feature <- color_feature
  }

  # Render the plot
  output$grouping_plot_2d <- renderPlot({
    do.call(plot_func, plot_params)
  })
}

render_3d_embedding <- function(obj, plot_func, input, color_factor, color_feature, embedding_dims, output) {
  # Get and validate component indices
  x_comp <- min(input$x_comp_3d, embedding_dims)
  y_comp <- min(input$y_comp_3d, embedding_dims)
  z_comp <- min(input$z_comp_3d, embedding_dims)

  # Build plot parameters
  plot_params <- list(
    romics_object=obj,
    Xcomp=x_comp,
    Ycomp=y_comp,
    Zcomp=z_comp,
    size=input$grouping_size,
    alpha=input$grouping_alpha
  )

  # Add type for PCA plots
  if (input$grouping_type == "pca") {
    plot_params$plotType <- "individual"
  }

  # Add color parameters
  if (!is.null(color_factor)) {
    plot_params$factor_name <- color_factor
  }
  if (!is.null(color_feature)) {
    plot_params$feature <- color_feature
  }

  # Create the plot and convert to plotly if needed
  plot_obj <- do.call(plot_func, plot_params)

  # Convert ggplot to plotly for consistent 3D rendering
  if (inherits(plot_obj, "ggplot")) {
    plot_obj <- plot_obj + theme(aspect.ratio=1)
    plot_obj <- plotly::ggplotly(plot_obj)
  }

  # Render the plot
  output$grouping_plot_3d <- renderPlotly({ plot_obj })
}

# ============================================================================
# Shiny App
# ============================================================================

ui <- dashboardPage(skin="black",
  dashboardHeader(title=tags$span(tags$img(src="romics_explorer_logo.png",height="45px",style="margin-right:10px;vertical-align:middle;"),
                                  "RomicsExplorer",style="display:flex;align-items:center;"),
                  titleWidth=250),
  dashboardSidebar(width=250,sidebarMenu(id="tabs",
    menuItem("Load Data",tabName="load",icon=icon("upload")),
    menuItem("Lipid Information",tabName="lipid_info",icon=icon("droplet")),
    menuItem("Grouping",tabName="grouping",icon=icon("layer-group")),
    menuItem("Single Feature Plot",tabName="single_feature",icon=icon("search")),
    menuItem("Multi Feature Plots",tabName="multi_feature",icon=icon("object-group")),
    menuItem("Statistics",tabName="statistics",icon=icon("chart-bar")),
    menuItem("Object History",tabName="history",icon=icon("history")))),
  dashboardBody(useShinyjs(),
    tags$script(HTML("$(document).ready(function(){setTimeout(function(){$('ul.sidebar-menu li').each(function(){var a=$(this).find('> a');var v=a.attr('data-value');if(v==='ontology'||v==='enrichment'||v==='lipid_info'||v==='lipid_enrichment')$(this).hide();});},200);});
      Shiny.addCustomMessageHandler('copyClipboard',function(msg){navigator.clipboard.writeText(msg.text).then(function(){},function(){var ta=document.createElement('textarea');ta.value=msg.text;ta.style.position='fixed';ta.style.opacity='0';document.body.appendChild(ta);ta.select();document.execCommand('copy');document.body.removeChild(ta);});});")),
    tags$head(tags$style(HTML(".main-header .logo {background-color: #D2EBEC !important;} .main-header .navbar {background-color: #308B90 !important;} .box-header.with-border {background-color: #36454F !important; color: white !important; border-bottom: 1px solid #555 !important;} .box-header {background-color: #36454F !important; color: white !important;} .box-primary>.box-header {background-color: #36454F !important;} .box-info>.box-header {background-color: #36454F !important;} .box-success>.box-header {background-color: #36454F !important;} .box {border-top-color: #6b7280 !important; border-color: #6b7280 !important;} .box.box-primary {border-color: #6b7280 !important;} .box.box-info {border-color: #6b7280 !important;} .box.box-success {border-color: #6b7280 !important;} .box.box-warning {border-color: #6b7280 !important;} .box.box-danger {border-color: #6b7280 !important;} .box-warning>.box-header {background: linear-gradient(90deg, #BE5504 0%, #d97d2e 100%) !important;} .box-danger>.box-header {background: linear-gradient(90deg, #e50000 0%, #ff3333 100%) !important;} .btn-primary {background-color: #36454F !important; border-color: #36454F !important; color: #CCCCCC !important;} .btn-primary:hover {background-color: #4a5a6f !important; color: white !important;} .btn-info {background-color: #36454F !important; border-color: #36454F !important; color: #CCCCCC !important;} .btn-info:hover {background-color: #4a5a6f !important; color: white !important;} .btn-success {background-color: #36454F !important; border-color: #36454F !important; color: #CCCCCC !important;} .btn-success:hover {background-color: #4a5a6f !important; color: white !important;} .btn-warning {background-color: #BE5504 !important; border-color: #BE5504 !important; color: white !important;} .btn-warning:hover {background-color: #d97d2e !important; color: white !important;} .btn-danger {background-color: #e50000 !important; border-color: #e50000 !important; color: white !important;} .btn-danger:hover {background-color: #ff3333 !important; color: white !important;} .btn-default {color: #CCCCCC !important;} .btn-default:hover {color: white !important;} .nav-tabs>li.active>a {background-color: #36454F !important; color: white !important;} .nav-tabs>li>a {color: #36454F !important;} .bubble-plot-container{overflow-y:auto;max-height:100vh;} .preview-box{background:#f9f9f9;border:1px solid #ddd;border-radius:4px;padding:10px;max-height:200px;overflow-y:auto;font-family:monospace;font-size:12px;} .steps-box{background:#f9f9f9;border:1px solid #ddd;border-radius:4px;padding:15px;max-height:500px;overflow-y:auto;font-family:monospace;font-size:13px;line-height:1.6;} .nav-tabs>li>a{font-weight:bold;} .ontology-term-list{background:#f9f9f9;border:1px solid #ddd;border-radius:4px;padding:10px;max-height:250px;overflow-y:auto;font-size:12px;line-height:1.8;} .ontology-term-list .term-number{color:#999;margin-right:5px;} .btn-copy{margin-top:5px;}"))),
    tabItems(
      tabItem(tabName="load",
        fluidRow(box(title="Load Romics Object",status="primary",solidHeader=TRUE,width=6,fileInput("rda_file","Select .rda/.RData:",accept=c(".rda",".RData",".Rda")),uiOutput("object_selector_ui"),actionButton("load_btn","Load Object",class="btn-primary btn-lg",icon=icon("upload"),width="100%"),br(),br(),verbatimTextOutput("load_status")),
          box(title="Object Status",status="info",solidHeader=TRUE,width=6,valueBoxOutput("object_status_box",width=12),valueBoxOutput("features_box",width=12),valueBoxOutput("samples_box",width=12))),
        fluidRow(box(title="Select ID Type",status="primary",solidHeader=TRUE,width=6,selectInput("id_select","ID Column:",choices=character(0)),actionButton("apply_id","Apply ID",class="btn-success",icon=icon("check"),width="100%"),br(),br(),h5("Preview:"),uiOutput("id_preview")),
          box(title="Select Main Factor",status="primary",solidHeader=TRUE,width=6,selectInput("factor_select","Factor:",choices=character(0)),actionButton("apply_factor","Apply Factor",class="btn-success",icon=icon("check"),width="100%"),br(),br(),h5("Preview:"),uiOutput("factor_preview"))),
        fluidRow(box(title="Object Summary",status="info",solidHeader=TRUE,width=12,verbatimTextOutput("object_summary"))),
        fluidRow(box(title="About RomicsExplorer",status="warning",solidHeader=TRUE,width=12,
          p("RomicsExplorer is an interactive Shiny application for visualizing and analyzing omics data using the RomicsProcessor R package.",style="font-size:12px;color:#333;margin-bottom:15px;"),
          h5("Authors:",style="margin-top:15px;margin-bottom:10px;"),
          tags$ul(
            tags$li("Geremy Clair ", tags$a(href="mailto:geremy.clair@pnnl.gov","geremy.clair@pnnl.gov")),
            tags$li("Harsh Bhotika ", tags$a(href="mailto:harsh.bhotika@pnnl.gov","harsh.bhotika@pnnl.gov")),
            tags$li("Brittney Gorman ", tags$a(href="mailto:brittney.gorman@pnnl.gov","brittney.gorman@pnnl.gov"))
          ),
          style="font-size:12px;line-height:1.8;"))),
      tabItem(tabName="ontology",
        # Enrichment ID Selection (common to all tabs)
        fluidRow(box(title="Enrichment ID Selection",status="warning",solidHeader=TRUE,width=12,p("Select UniProt Accession or Entry Name column.",style="font-size:13px;color:#555;"),
          fluidRow(column(4,selectInput("enrichment_id_select","UniProt ID Column:",choices=character(0))),column(4,br(),actionButton("apply_enrichment_id","Set ID",class="btn-warning",icon=icon("key"),width="100%")),column(4,h5("Preview:"),uiOutput("enrichment_id_preview"))))),

        # Tabbed interface for loading methods
        tabsetPanel(type="tabs",
          # TAB 1: BULK LOAD
          tabPanel("Bulk Load (Fast)",
            br(),
            fluidRow(box(title="Load All Tables from R Object",status="warning",solidHeader=TRUE,width=12,
              p("Recommended: Load all ontology tables (UniProtTable, GO, KEGG, Reactome) from a single R object file.",style="font-size:12px;color:#666;"),
              br(),
              fluidRow(
                column(3,fileInput("bulk_ontology_file","Upload R Object:",accept=c(".rda",".RData",".rds"))),
                column(2,br(),actionButton("bulk_load_btn","Load All",class="btn-warning btn-lg",icon=icon("upload"),width="100%")),
                column(7,uiOutput("bulk_load_status"))
              ),
              hr(),
              h5("Loaded Table Previews:"),
              fluidRow(
                box(title="UniProt Table",status="primary",solidHeader=TRUE,width=3,
                  uiOutput("bulk_uniprot_status"),h5("Preview:"),uiOutput("bulk_uniprot_preview")),
                box(title="GO Table",status="primary",solidHeader=TRUE,width=3,
                  uiOutput("bulk_go_status"),h5("Preview:"),uiOutput("bulk_go_preview")),
                box(title="KEGG Table",status="primary",solidHeader=TRUE,width=3,
                  uiOutput("bulk_kegg_status"),h5("Preview:"),uiOutput("bulk_kegg_preview")),
                box(title="Reactome Table",status="primary",solidHeader=TRUE,width=3,
                  uiOutput("bulk_reactome_status"),h5("Preview:"),uiOutput("bulk_reactome_preview"))
              ),
              hr(),
              p("Supported patterns: UniProtTable, UniProtTable_GO, UniProtTable_KEGG, UniProtTable_Reactome (and lowercase variants)",style="font-size:11px;color:#999;font-style:italic;")
            ))
          ),

          # TAB 2: INDIVIDUAL LOAD
          tabPanel("Individual Load (Manual)",
            br(),
            fluidRow(box(title="UniProt Reference Table (Required)",status="danger",solidHeader=TRUE,width=12,p("Required for ProteinMiniOn. Upload UniProtTable.",style="font-size:13px;color:#555;"),
              fluidRow(column(4,fileInput("uniprot_table_file","Upload UniProtTable:",accept=c(".csv",".tsv",".rda",".RData",".rds"))),column(4,br(),actionButton("uniprot_table_load_btn","Load UniProtTable",class="btn-danger",icon=icon("upload"),width="100%")),column(4,uiOutput("uniprot_table_status"))),
              h5("Preview (first 5 rows):"),div(style="max-height:150px;overflow-y:auto;",DTOutput("uniprot_table_preview",height="auto")))),
            br(),
            fluidRow(box(title="GO",status="primary",solidHeader=TRUE,width=4,fileInput("go_file","Upload GO:",accept=c(".csv",".tsv",".rda",".RData",".rds")),actionButton("go_load_btn","Load",class="btn-primary",icon=icon("upload"),width="100%"),hr(),uiOutput("go_status_indicator"),h5("Preview:"),uiOutput("go_terms_preview")),
              box(title="KEGG",status="success",solidHeader=TRUE,width=4,fileInput("kegg_file","Upload KEGG:",accept=c(".csv",".tsv",".rda",".RData",".rds")),actionButton("kegg_load_btn","Load",class="btn-success",icon=icon("upload"),width="100%"),hr(),uiOutput("kegg_status_indicator"),h5("Preview:"),uiOutput("kegg_terms_preview")),
              box(title="Reactome",status="info",solidHeader=TRUE,width=4,fileInput("reactome_file","Upload Reactome:",accept=c(".csv",".tsv",".rda",".RData",".rds")),actionButton("reactome_load_btn","Load",class="btn-info",icon=icon("upload"),width="100%"),hr(),uiOutput("reactome_status_indicator"),h5("Preview:"),uiOutput("reactome_terms_preview")))
          )
        ),

        # Summary (common to all tabs)
        br(),
        fluidRow(box(title="Summary",status="info",solidHeader=TRUE,width=12,verbatimTextOutput("ontology_summary")))),
      tabItem(tabName="grouping",
        fluidRow(box(title="Embedding Options",status="primary",solidHeader=TRUE,width=3,selectInput("grouping_type","Embedding:",choices=c("No embeddings"="none")),uiOutput("plot_dimension_ui"),
          conditionalPanel("input.plot_dimension=='2d'",numericInput("x_comp","X:",1,min=1),numericInput("y_comp","Y:",2,min=1)),conditionalPanel("input.plot_dimension=='3d'",numericInput("x_comp_3d","X:",1,min=1),numericInput("y_comp_3d","Y:",2,min=1),numericInput("z_comp_3d","Z:",3,min=1)),
          tags$hr(),h5("Coloring"),selectInput("grouping_color_by","Color by:",choices=c("Factor"="factor","Feature"="feature")),conditionalPanel("input.grouping_color_by=='factor'",selectInput("grouping_factor_select","Factor:",choices=character(0))),
          conditionalPanel("input.grouping_color_by=='feature'",selectizeInput("grouping_feature_select","Feature:",choices=NULL,options=list(placeholder="Search..."))),
          tags$hr(),h5("Point Display"),sliderInput("grouping_size","Size:",value=6,min=1,max=10,step=0.5),sliderInput("grouping_alpha","Opacity:",value=0.8,min=0.1,max=1,step=0.1),
          tags$hr(),actionButton("plot_grouping","Generate",class="btn-success btn-lg",icon=icon("project-diagram"),width="100%")),
          box(title="Embedding Plot",status="primary",solidHeader=TRUE,width=9,conditionalPanel("input.plot_dimension=='2d'",plotOutput("grouping_plot_2d",height="700px")),conditionalPanel("input.plot_dimension=='3d'",plotlyOutput("grouping_plot_3d",height="700px")))),
        fluidRow(box(title="Embeddings",status="info",solidHeader=TRUE,width=6,verbatimTextOutput("grouping_embedding_info")),box(title="Settings",status="info",solidHeader=TRUE,width=6,verbatimTextOutput("grouping_plot_settings")))),
      tabItem(tabName="single_feature",
        fluidRow(box(title="Options",status="primary",solidHeader=TRUE,width=3,selectizeInput("feature_select","Feature:",choices=NULL,options=list(placeholder="Search...")),selectInput("plot_type","Type:",choices=c("Jitter+Box"="jb","Jitter+Violin"="jv","Jitter"="jitter","Boxplot"="boxplot","Violin"="violin")),selectInput("feature_factor","Factor:",choices=character(0)),
          radioButtons("feature_plot_format","Format:",choices=c("Plotly"="plotly","Static"="static"),selected="plotly",inline=TRUE),textInput("custom_title","Title:",placeholder="Auto"),actionButton("plot_feature","Generate",class="btn-success btn-lg",icon=icon("chart-line"),width="100%")),
          box(title="Feature Plot",status="primary",solidHeader=TRUE,width=9,uiOutput("feature_plot_ui"),hr(),verbatimTextOutput("feature_plot_info"),hr(),fluidRow(column(12,box(width=12,status="info",solidHeader=TRUE,title="Statistics Table - P-values",h5("Filter by P-value Type:"),radioButtons("feature_pval_type","",choices=c("All"="all","p-value only"="p","Adjusted p-value only"="padj"),selected="all",inline=TRUE),DTOutput("feature_stats_table"))))))),
      tabItem(tabName="multi_feature",tabBox(id="multi_tabbox",width=12,title="Multi-Feature Analysis",
        tabPanel("Bubble Plot",icon=icon("circle"),br(),fluidRow(column(3,box(width=12,status="primary",solidHeader=TRUE,title="Bubble Options",textInput("multi_feature_input","Features:",placeholder="feat1, feat2"),br(),selectizeInput("multi_feature_dropdown","Add:",choices=NULL,options=list(placeholder="Select...")),actionButton("add_feature_btn","Add",class="btn-info",icon=icon("plus"),width="100%"),br(),br(),selectInput("multi_factor_select","Factor:",choices=character(0)),checkboxInput("multi_scale_features","Scale",value=TRUE),selectInput("multi_viridis_option","Palette:",choices=c("viridis","plasma","inferno","magma","cividis")),tags$hr(),actionButton("plot_multi","Generate",class="btn-success btn-lg",icon=icon("object-group"),width="100%"),br(),br(),downloadButton("download_bubble_plot_png","Download PNG",class="btn-info",style="width:100%;"))),column(9,box(width=12,status="primary",solidHeader=TRUE,title="Bubble Plot",div(class="bubble-plot-container",plotlyOutput("multi_bubble_plot",height="600px")))))),
        tabPanel("Volcano Plot",icon=icon("chart-line"),br(),fluidRow(column(3,box(width=12,status="primary",solidHeader=TRUE,title="Volcano Options",selectInput("volcano_factor_select","Factor:",choices=character(0)),uiOutput("volcano_comparison_ui"),br(),uiOutput("volcano_pvalue_type_ui"),br(),uiOutput("volcano_pval_threshold_ui"),numericInput("volcano_fc_threshold","Min |Log(FC)|:",value=0,min=0,step=0.1),tags$hr(),h5("Point Display"),sliderInput("volcano_size","Size:",value=2,min=0.5,max=8,step=0.5),sliderInput("volcano_alpha","Opacity:",value=0.5,min=0.1,max=1,step=0.1),tags$hr(),actionButton("plot_volcano","Generate",class="btn-success btn-lg",icon=icon("chart-line"),width="100%"))),column(9,box(width=12,status="primary",solidHeader=TRUE,title="Volcano",uiOutput("volcano_plot_ui")))),
          fluidRow(box(title="Volcano Info",status="info",solidHeader=TRUE,width=6,verbatimTextOutput("volcano_plot_info")),box(title="Extract Proteins",status="warning",solidHeader=TRUE,width=6,selectInput("volcano_extract_type","Extract:",choices=c("All Significant"="all_sig","Up (positive FC)"="up","Down (negative FC)"="down")),verbatimTextOutput("volcano_extract_summary"),br(),fluidRow(column(6,actionButton("volcano_copy_ids","Copy IDs",class="btn-warning",icon=icon("copy"),style="width:100%;")),column(6,downloadButton("volcano_download_ids","Download CSV",class="btn-info",style="width:100%;")))))),
        tabPanel("Heatmap",icon=icon("th"),br(),fluidRow(column(3,box(width=12,status="primary",solidHeader=TRUE,title="Heatmap Options",textInput("hm_feature_input","Features:",placeholder="feat1, feat2"),br(),selectizeInput("hm_feature_dropdown","Add:",choices=c("All"="ALL_FEATURES"),options=list(placeholder="Select...")),actionButton("hm_add_feature_btn","Add",class="btn-info",icon=icon("plus"),width="100%"),br(),br(),selectInput("hm_factor_select","Factor:",choices=character(0)),tags$hr(),checkboxInput("hm_scale_feature","Scale (Z-score)",value=TRUE),checkboxInput("hm_show_completeness","Completeness",value=FALSE),checkboxInput("hm_show_dendrogram","Dendrogram",value=TRUE),checkboxInput("hm_show_feature_names","Feature Names",value=TRUE),tags$hr(),checkboxInput("hm_show_clusters","Show Clusters",value=FALSE),conditionalPanel("input.hm_show_clusters==true",numericInput("hm_n_clusters","N Clusters:",value=3,min=2,max=50),checkboxInput("hm_show_cluster_legend","Cluster Legend",value=TRUE)),selectInput("hm_clustering_method","Method:",choices=c("ward.D","ward.D2","complete","average","single"),selected="ward.D"),tags$hr(),h5("Feature Filtering"),selectInput("hm_filter_stat_column","Filter by Statistic:",choices=c("None"="none")),conditionalPanel("input.hm_filter_stat_column!='none'",textInput("hm_stat_column_filter","Filter Expression:",value="<0.05")),tags$hr(),h5("ANOVA Quick Filter"),uiOutput("hm_anova_filter_ui"),tags$hr(),selectInput("hm_viridis_option","Color:",choices=c("viridis","plasma","inferno","magma","cividis")),tags$hr(),actionButton("plot_heatmap","Generate",class="btn-success btn-lg",icon=icon("th"),width="100%"),br(),br(),downloadButton("download_heatmap_png","Download PNG",class="btn-info",style="width:100%;"))),column(9,box(width=12,status="primary",solidHeader=TRUE,title="Heatmap",div(style="min-height:100px;background:#f5f5f5;padding:15px;border-radius:4px;",uiOutput("heatmap_plot_ui"))))),
          fluidRow(box(title="Features",status="info",solidHeader=TRUE,width=4,verbatimTextOutput("hm_selected_features")),box(title="Info",status="info",solidHeader=TRUE,width=4,verbatimTextOutput("hm_plot_info")),box(title="Clusters",status="warning",solidHeader=TRUE,width=4,conditionalPanel("input.hm_show_clusters==true",selectInput("hm_cluster_view","View:",choices=c("All"="all")),br(),DTOutput("hm_cluster_table"),br(),downloadButton("download_cluster_csv","Download CSV",class="btn-warning",style="width:100%;"),actionButton("copy_cluster_ids","Copy IDs",class="btn-default btn-copy",icon=icon("copy"),style="width:100%;")),conditionalPanel("input.hm_show_clusters!=true",p("Enable clusters.",style="color:#999;font-style:italic;"))))))),
      tabItem(tabName="statistics",
        fluidRow(box(title="Column Filters",status="primary",solidHeader=TRUE,width=12,collapsible=TRUE,p("Select up to 6 columns and apply filters.",style="font-size:12px;color:#777;"),
          fluidRow(column(4,selectInput("stat_filter_col_1","Column 1:",choices=c("None"="none"))),column(4,textInput("stat_filter_expr_1","Filter 1:",placeholder="e.g., <0.05")),column(4,br(),checkboxInput("stat_filter_active_1","Active",value=FALSE))),
          conditionalPanel("input.stat_filter_col_1!='none'",fluidRow(column(4,selectInput("stat_filter_col_2","Column 2:",choices=c("None"="none"))),column(4,textInput("stat_filter_expr_2","Filter 2:")),column(4,br(),checkboxInput("stat_filter_active_2","Active",value=FALSE)))),
          conditionalPanel("input.stat_filter_col_2!='none'",fluidRow(column(4,selectInput("stat_filter_col_3","Column 3:",choices=c("None"="none"))),column(4,textInput("stat_filter_expr_3","Filter 3:")),column(4,br(),checkboxInput("stat_filter_active_3","Active",value=FALSE)))),
          conditionalPanel("input.stat_filter_col_3!='none'",fluidRow(column(4,selectInput("stat_filter_col_4","Column 4:",choices=c("None"="none"))),column(4,textInput("stat_filter_expr_4","Filter 4:")),column(4,br(),checkboxInput("stat_filter_active_4","Active",value=FALSE)))),
          conditionalPanel("input.stat_filter_col_4!='none'",fluidRow(column(4,selectInput("stat_filter_col_5","Column 5:",choices=c("None"="none"))),column(4,textInput("stat_filter_expr_5","Filter 5:")),column(4,br(),checkboxInput("stat_filter_active_5","Active",value=FALSE)))),
          conditionalPanel("input.stat_filter_col_5!='none'",fluidRow(column(4,selectInput("stat_filter_col_6","Column 6:",choices=c("None"="none"))),column(4,textInput("stat_filter_expr_6","Filter 6:")),column(4,br(),checkboxInput("stat_filter_active_6","Active",value=FALSE)))),
          fluidRow(column(4,actionButton("apply_stat_filters","Apply Filters",class="btn-primary",icon=icon("filter"),width="100%")),column(4,actionButton("reset_stat_filters","Reset",class="btn-default",icon=icon("undo"),width="100%")),column(4,verbatimTextOutput("stat_filter_summary"))))),
        fluidRow(box(title="Statistics Table",status="primary",solidHeader=TRUE,width=12,DTOutput("stats_table"),br(),
          fluidRow(column(3,downloadButton("download_stats","Download Full CSV",class="btn-info",style="width:100%;")),column(3,downloadButton("download_stats_filtered","Download Filtered CSV",class="btn-success",style="width:100%;")),column(3,actionButton("copy_stats_ids","Copy All IDs",class="btn-default btn-copy",icon=icon("copy"),style="width:100%;")),column(3,actionButton("copy_stats_filtered_ids","Copy Filtered IDs",class="btn-default btn-copy",icon=icon("copy"),style="width:100%;")))))),
      tabItem(tabName="enrichment",
        fluidRow(box(title="Enrichment Mode",status="primary",solidHeader=TRUE,width=12,radioButtons("enrich_mode","Mode:",choices=c("List Comparison"="list","Ranked (KS)"="ranked"),selected="list",inline=TRUE),uiOutput("enrich_available_ontologies"))),
        conditionalPanel("input.enrich_mode=='list'",fluidRow(
          column(4,box(width=12,status="primary",solidHeader=TRUE,title="Settings",selectInput("enrich_test_type","Test:",choices=c("EASE"="EASE","Fisher"="Fisher","Binomial"="Binomial")),radioButtons("enrich_list_submode","Selection:",choices=c("Manual"="manual","Stat Filter"="stat_filter"),selected="manual",inline=TRUE),tags$hr(),h5("Ontologies:"),uiOutput("enrich_ontology_checkboxes"),tags$hr(),actionButton("run_enrichment","Run Enrichment",class="btn-success btn-lg",icon=icon("flask"),width="100%"))),
          column(4,box(width=12,status="warning",solidHeader=TRUE,title="Query Proteins",conditionalPanel("input.enrich_list_submode=='manual'",p("Clusters & ANOVA at top.",style="font-size:12px;color:#777;"),selectizeInput("enrich_query_dropdown","Add:",choices=NULL,options=list(placeholder="Select...")),actionButton("enrich_query_add_btn","Add to Query",class="btn-warning",icon=icon("plus"),width="100%"),br(),br(),textInput("enrich_query_text","Query:",placeholder="Proteins...")),conditionalPanel("input.enrich_list_submode=='stat_filter'",p("Filter by up to 4 stat columns.",style="font-size:12px;color:#777;"),fluidRow(column(7,selectInput("enrich_stat_col_1","Col 1:",choices=c("None"="none"))),column(5,textInput("enrich_stat_filter_1","Filter:",value="<0.05"))),fluidRow(column(7,selectInput("enrich_stat_col_2","Col 2:",choices=c("None"="none"))),column(5,textInput("enrich_stat_filter_2","Filter:",value=""))),fluidRow(column(7,selectInput("enrich_stat_col_3","Col 3:",choices=c("None"="none"))),column(5,textInput("enrich_stat_filter_3","Filter:",value=""))),fluidRow(column(7,selectInput("enrich_stat_col_4","Col 4:",choices=c("None"="none"))),column(5,textInput("enrich_stat_filter_4","Filter:",value=""))),br(),actionButton("enrich_apply_query_filter","Apply Filters",class="btn-warning",icon=icon("filter"),width="100%")),tags$hr(),h5("Query Summary:"),verbatimTextOutput("enrich_query_summary"))),
          column(4,box(width=12,status="info",solidHeader=TRUE,title="Universe",radioButtons("enrich_universe_mode","Universe:",choices=c("All Detected"="all","Custom"="custom","Stat Filtered"="stat_filter"),selected="all"),conditionalPanel("input.enrich_universe_mode=='custom'",selectizeInput("enrich_universe_dropdown","Add:",choices=NULL,options=list(placeholder="Select...")),actionButton("enrich_universe_add_btn","Add",class="btn-info",icon=icon("plus"),width="100%"),br(),br(),textInput("enrich_universe_text","Universe:",placeholder="Proteins...")),conditionalPanel("input.enrich_universe_mode=='stat_filter'",fluidRow(column(7,selectInput("enrich_univ_stat_col","Col:",choices=c("None"="none"))),column(5,textInput("enrich_univ_stat_filter","Filter:",value=""))),actionButton("enrich_apply_univ_filter","Apply",class="btn-info",icon=icon("filter"),width="100%")),tags$hr(),h5("Universe Summary:"),verbatimTextOutput("enrich_universe_summary"))))),
        conditionalPanel("input.enrich_mode=='ranked'",fluidRow(column(4,box(width=12,status="primary",solidHeader=TRUE,title="KS Settings",selectInput("enrich_ks_ranking_col","Ranking Column:",choices=character(0)),selectInput("enrich_ks_order","Order:",choices=c("Ascending"="ascending","Descending"="descending")),tags$hr(),h5("Ontologies:"),uiOutput("enrich_ks_ontology_checkboxes"),tags$hr(),actionButton("run_ks_enrichment","Run KS",class="btn-success btn-lg",icon=icon("flask"),width="100%"))),column(8,box(width=12,status="info",solidHeader=TRUE,title="Ranking Preview",verbatimTextOutput("enrich_ks_preview"))))),
        fluidRow(box(title="Enrichment Results",status="success",solidHeader=TRUE,width=12,
          fluidRow(column(2,selectInput("enrich_filter_category","Category:",choices=c("All"="all"))),column(2,selectInput("enrich_filter_ptype","P-value Type:",choices=c("p"="p","padj"="padj"))),column(2,numericInput("enrich_filter_pmax","Max P:",value=0.05,min=0,max=1,step=0.01)),column(2,numericInput("enrich_filter_fc","Min FC:",value=1,min=0,step=0.1)),column(2,numericInput("enrich_filter_minfeature","Min Features:",value=2,min=1)),column(2,br(),actionButton("apply_enrich_filters","Apply Filters",class="btn-success",icon=icon("filter"),width="100%"))),
          DTOutput("enrich_results_table"),br(),
          fluidRow(column(3,downloadButton("download_enrichment_csv","Download CSV",class="btn-success",style="width:100%;")),column(3,actionButton("copy_enrichment_ids","Copy Term IDs",class="btn-default btn-copy",icon=icon("copy"),style="width:100%;")),column(6,verbatimTextOutput("enrich_results_summary"))),
          tags$hr(),h5("Enrichment Visualization"),fluidRow(column(3,actionButton("generate_enrich_plot","Generate Plot",class="btn-primary",icon=icon("chart-bar"),width="100%")),column(3,downloadButton("download_enrich_plot","Download Plot PNG",class="btn-info",style="width:100%;"))),br(),plotlyOutput("enrich_bar_plot",height="600px")))),

      # ============================================================================
      # LIPID INFORMATION TAB
      # ============================================================================
      tabItem(tabName="lipid_info",
        tabsetPanel(type="tabs",
          tabPanel("LipidMiner Results",
            br(),
            fluidRow(box(title="LipidMiner Recognition",status="info",solidHeader=TRUE,width=12,
              p("Results from lipidMiner function - shorthand nomenclature recognition",style="font-size:12px;color:#555;"),
              verbatimTextOutput("lipidminer_summary"),
              hr(),
              uiOutput("lipidminer_warning"),
              DTOutput("lipidminer_table")
            ))
          ),
          tabPanel("Lipid Ontology Results",
            br(),
            fluidRow(box(title="Lipid Ontology Annotation",status="info",solidHeader=TRUE,width=12,
              p("Results from Lipid_ontology_maker function - lipid classification and information",style="font-size:12px;color:#555;"),
              DTOutput("lipid_ontology_table")
            ))
          )
        )
      ),

      # ============================================================================
      # LIPID ENRICHMENT TAB
      # ============================================================================
      tabItem(tabName="lipid_enrichment",
        fluidRow(box(title="Lipid Enrichment Analysis",status="primary",solidHeader=TRUE,width=12,
          p("Perform enrichment analysis on detected lipids using shorthand nomenclature",style="font-size:12px;color:#555;"),
          radioButtons("lipid_enrich_mode","Mode:",choices=c("Query vs Universe"="query_universe","Over Representation"=" ora"),selected="query_universe",inline=TRUE),
          tags$hr(),
          fluidRow(
            column(6,box(width=12,status="warning",solidHeader=TRUE,title="Query Lipids",
              selectizeInput("lipid_query_dropdown","Add Lipid:",choices=NULL,options=list(placeholder="Search...")),
              actionButton("lipid_query_add_btn","Add",class="btn-warning",icon=icon("plus"),width="100%"),
              br(),br(),
              textInput("lipid_query_text","Query:",placeholder="Lipid1, Lipid2, ..."),
              br(),
              verbatimTextOutput("lipid_query_summary")
            )),
            column(6,box(width=12,status="info",solidHeader=TRUE,title="Universe Lipids",
              radioButtons("lipid_universe_mode","Universe:",choices=c("All Detected"="all","Custom"="custom"),selected="all"),
              conditionalPanel("input.lipid_universe_mode=='custom'",
                selectizeInput("lipid_universe_dropdown","Add Lipid:",choices=NULL,options=list(placeholder="Search...")),
                actionButton("lipid_universe_add_btn","Add",class="btn-info",icon=icon("plus"),width="100%"),
                br(),br(),
                textInput("lipid_universe_text","Universe:",placeholder="Lipid1, Lipid2, ...")
              ),
              br(),
              verbatimTextOutput("lipid_universe_summary")
            ))
          ),
          tags$hr(),
          fluidRow(
            column(6,box(width=12,status="warning",solidHeader=TRUE,title="Enrichment Method",
              radioButtons("lipid_enrich_method","Test:",choices=c("Fisher's Exact"="fisher","EASE"="ease","Hypergeometric"="hyper","Binomial"="binom"),selected="fisher",inline=FALSE),
              numericInput("lipid_enrich_p","P-value cutoff:",value=0.05,min=0,max=1,step=0.01),
              numericInput("lipid_enrich_adj_p","Adjusted p-value cutoff:",value=1.0,min=0,max=1,step=0.01)
            )),
            column(6,box(width=12,status="info",solidHeader=TRUE,title="Actions",
              actionButton("generate_lipid_query_plot","Query Plot",class="btn-primary",icon=icon("chart-bar"),width="100%"),
              br(),br(),
              actionButton("generate_lipid_universe_plot","Universe Plot",class="btn-info",icon=icon("chart-bar"),width="100%"),
              br(),br(),
              actionButton("run_lipid_enrichment","Run Enrichment",class="btn-success btn-lg",icon=icon("flask"),width="100%")
            ))
          )
        )),
        br(),
        fluidRow(box(title="Query Plot",status="primary",solidHeader=TRUE,width=6,plotlyOutput("lipid_query_plot",height="500px"))),
        fluidRow(box(title="Universe Plot",status="info",solidHeader=TRUE,width=6,plotlyOutput("lipid_universe_plot",height="500px"))),
        br(),
        fluidRow(box(title="Enrichment Results",status="success",solidHeader=TRUE,width=12,
          DTOutput("lipid_enrichment_results"),
          br(),
          fluidRow(column(6,downloadButton("download_lipid_enrichment","Download Results",class="btn-success",style="width:100%;")),column(6,actionButton("copy_lipid_enrichment_ids","Copy Lipid IDs",class="btn-default btn-copy",icon=icon("copy"),style="width:100%;")))
        ))
      ),

      tabItem(tabName="history",fluidRow(uiOutput("uuid_box")),
        fluidRow(box(title="Object Info",status="info",solidHeader=TRUE,width=6,verbatimTextOutput("history_object_summary")),box(title="Embeddings",status="info",solidHeader=TRUE,width=6,verbatimTextOutput("embedding_info"))),
        fluidRow(box(title="Steps (Detail)",status="primary",solidHeader=TRUE,width=12,checkboxInput("show_timestamps","Show Timestamps",value=TRUE),verbatimTextOutput("steps_output_detailed"))),
        fluidRow(box(title="Dependencies",status="warning",solidHeader=TRUE,width=12,collapsible=TRUE,collapsed=TRUE,DTOutput("dependencies_table"))),
        fluidRow(box(title="Powered By",status="info",solidHeader=TRUE,width=12,HTML("<p style='font-size: 14px; line-height: 1.6;'>
                  This application is powered by <strong>RomicsProcessor</strong>, an R package for analyzing omics data with reproducible and FAIR-compatible workflows.
                </p>
                <p style='font-size: 13px;'>
                  <strong>GitHub:</strong> <a href='https://github.com/PNNL-Comp-Mass-Spec/RomicsProcessor' target='_blank'>https://github.com/PNNL-Comp-Mass-Spec/RomicsProcessor</a>
                </p>"))))
    )
  )
)

# ============================================================================
# 7. SHINY SERVER - Main Application Logic
# ============================================================================

server <- function(input, output, session) {

  # ---- Reactive Values ----
  # Store global state and data throughout the session
  romics_object <- reactiveVal(NULL)
  current_factor <- reactiveVal(NULL)
  current_id_type <- reactiveVal(NULL)
  loaded_env <- reactiveVal(NULL)
  multi_plot_store <- reactiveVal(NULL)
  heatmap_plot_store <- reactiveVal(NULL)
  heatmap_clusters <- reactiveVal(NULL)
  volcano_comparisons <- reactiveVal(NULL)
  volcano_pvalue_types <- reactiveVal(list(p=FALSE, padj=FALSE))
  volcano_data_store <- reactiveVal(NULL)
  feature_stats_store <- reactiveVal(NULL)
  is_proteomics <- reactiveVal(FALSE)
  is_lipidominion <- reactiveVal(FALSE)
  go_table <- reactiveVal(NULL)
  kegg_table <- reactiveVal(NULL)
  reactome_table <- reactiveVal(NULL)
  uniprot_table <- reactiveVal(NULL)
  enrichment_id_col <- reactiveVal(NULL)
  enrich_results <- reactiveVal(NULL)
  enrich_filtered_data <- reactiveVal(NULL)
  enrich_plot_store <- reactiveVal(NULL)
  stats_filtered_data <- reactiveVal(NULL)

  # ---- Helper Functions for UI Updates ----
  # Centralized functions to update selectize/select inputs

  uafs <- function(s, fi) {
    updateSelectizeInput(s, "feature_select", choices=fi, server=TRUE)
    updateSelectizeInput(s, "multi_feature_dropdown", choices=fi, server=TRUE)
    updateSelectizeInput(s, "grouping_feature_select", choices=fi, server=TRUE)
  }
  uafacs <- function(s, fn, sel) {
    for (sid in c("factor_select","feature_factor","grouping_factor_select","multi_factor_select","volcano_factor_select","hm_factor_select"))
      updateSelectInput(s, sid, choices=fn, selected=sel)
  }
  
  # ============================================================================
  # 8.1 DATA LOADING - File Upload & Object Selection
  # ============================================================================

  observeEvent(input$rda_file, {
    req(input$rda_file)
    withProgress(message="Reading file...", value=0, {
      tryCatch({
        incProgress(0.1)
        env <- new.env(parent=emptyenv())
        load(input$rda_file$datapath, envir=env)
        loaded_env(env)
        incProgress(0.5)
        on <- ls(envir=env); rn <- c(); otn <- c()
        for (nm in on) { o <- get(nm, envir=env); if (is.list(o)&&"data" %in% names(o)) rn <- c(rn,nm) else otn <- c(otn,nm) }
        incProgress(0.3)
        output$object_selector_ui <- renderUI({ selectInput("object_select","Object:", choices=c(rn,otn), selected=if (length(rn)>0) rn[1] else c(rn,otn)[1]) })
        output$load_status <- renderText(paste0(length(on)," object(s): ",paste(on,collapse=", ")))
      }, error=function(e) output$load_status <- renderText(paste("Error:",e$message)))
    })
  })
  
  observeEvent(input$load_btn, {
    req(loaded_env(), input$object_select)
    tryCatch({
      obj <- get(input$object_select, envir=loaded_env()); romics_object(obj)
      fn <- get_clean_factors(obj); if (length(fn)>0) { uafacs(session,fn,fn[1]); current_factor(fn[1]) }
      idc <- character(0); if (!is.null(obj$IDs)&&is.data.frame(obj$IDs)) idc <- colnames(obj$IDs)
      if (length(idc)>0) { updateSelectInput(session,"id_select",choices=idc,selected=idc[1]); current_id_type(idc[1]) }
      ftn <- rownames(obj$data)
      if (length(ftn)>0) { uafs(session,ftn); updateSelectizeInput(session,"feature_select",selected=ftn[1]) }
      output$load_status <- renderText(paste0("'",input$object_select,"': ",nrow(obj$data)," x ",ncol(obj$data)))
    }, error=function(e) output$load_status <- renderText(paste("Error:",e$message)))
  })
  
  # ============================================================================
  # 8.2 DATA PREPARATION - Proteomics Detection & UI Updates
  # ============================================================================

  observe({
    obj <- romics_object()
    is_prot <- !is.null(obj)&&!is.null(obj$omics_type)&&grepl("proteomics", obj$omics_type, ignore.case=TRUE)
    is_lipid <- !is.null(obj)&&!is.null(obj$omics_type)&&grepl("lipid", obj$omics_type, ignore.case=TRUE)

    is_proteomics(is_prot)
    is_lipidominion(is_lipid)

    if (is_prot) {
    } else if (is_lipid) {
    } else {
    }
  })
  
  # ---- LOAD TAB OUTPUTS ----
  output$id_preview <- renderUI({
    req(input$id_select, romics_object()); obj <- romics_object()
    if (is.null(obj$IDs)||!input$id_select %in% colnames(obj$IDs)) return(HTML("<div class='preview-box'>No data</div>"))
    HTML(paste0("<div class='preview-box'>",paste(head(obj$IDs[[input$id_select]],5),collapse="<br/>"),"</div>"))
  })
  observeEvent(input$apply_id, {
    req(input$id_select, romics_object())
    tryCatch({ obj <- RomicsProcessor::romicsChangeIDs(romics_object(),newIDs=input$id_select); romics_object(obj); current_id_type(input$id_select)
    uafs(session,rownames(obj$data)); showNotification(paste("ID:",input$id_select),type="message")
    }, error=function(e) showNotification(paste("Error:",e$message),type="error"))
  })
  output$factor_preview <- renderUI({
    req(input$factor_select, romics_object())
    tryCatch({ fd <- RomicsProcessor::romicsExtractFactor(romics_object(),factor=input$factor_select)
    HTML(paste0("<div class='preview-box'>",paste(levels(fd),collapse="<br/>"),"</div>"))
    }, error=function(e) HTML(paste0("<div class='preview-box'>Error: ",e$message,"</div>")))
  })
  observeEvent(input$apply_factor, {
    req(input$factor_select, romics_object())
    tryCatch({ obj <- romics_object(); obj$main_factor <- input$factor_select; romics_object(obj); current_factor(input$factor_select)
    fc <- get_clean_factors(obj); if (input$factor_select %in% fc) fc <- c(input$factor_select,setdiff(fc,input$factor_select))
    if (length(fc)>0) for (s in c("feature_factor","grouping_factor_select","multi_factor_select","volcano_factor_select","hm_factor_select")) updateSelectInput(session,s,choices=fc,selected=input$factor_select)
    showNotification(paste("Factor:",input$factor_select),type="message")
    }, error=function(e) showNotification(paste("Error:",e$message),type="error"))
  })
  output$object_status_box <- renderValueBox({ if (is.null(romics_object())) valueBox("Not Loaded","Status",icon=icon("times-circle"),color="red") else valueBox("Loaded","Status",icon=icon("check-circle"),color="green") })
  output$features_box <- renderValueBox({ valueBox(if (!is.null(romics_object())) nrow(romics_object()$data) else "-","Features",icon=icon("list"),color="blue") })
  output$samples_box <- renderValueBox({ valueBox(if (!is.null(romics_object())) ncol(romics_object()$data) else "-","Samples",icon=icon("vials"),color="purple") })
  output$object_summary <- renderPrint({
    req(romics_object()); obj <- romics_object()
    cat("Dim:",nrow(obj$data),"x",ncol(obj$data),"\n")
    if (!is.null(obj$omics_type)) cat("Omics:",obj$omics_type,"\n")
    if (!is.null(current_factor())) cat("Factor:",current_factor(),"\n")
    for (l in c("data","metadata","IDs","statistics","steps","embeddings")) cat(" -",paste0(l,":"),ifelse(!is.null(obj[[l]]),"[YES]","[NO]"),"\n")
  })
  
  # ---- ONTOLOGY TAB ----
  observe({
    req(romics_object(), is_proteomics()); obj <- romics_object()
    if (!is.null(obj$IDs)&&is.data.frame(obj$IDs)) {
      ic <- colnames(obj$IDs); uh <- ic[grepl("uniprot|accession|entry",tolower(ic))]
      default_id <- if (length(uh)>0) uh[1] else ic[1]
      updateSelectInput(session,"enrichment_id_select",choices=c(uh,setdiff(ic,uh)),selected=default_id)
      enrichment_id_col(default_id)
    }
  })
  output$enrichment_id_preview <- renderUI({
    req(romics_object(), input$enrichment_id_select); obj <- romics_object()
    if (is.null(obj$IDs)||!input$enrichment_id_select %in% colnames(obj$IDs)) return(HTML("<div class='preview-box'>No data</div>"))
    ids <- head(obj$IDs[[input$enrichment_id_select]],5); ids <- ids[!is.na(ids)&ids!=""]
    HTML(paste0("<div class='preview-box'>",paste(ids,collapse="<br/>"),"</div>"))
  })
  observeEvent(input$apply_enrichment_id, { req(input$enrichment_id_select); enrichment_id_col(input$enrichment_id_select); showNotification(paste("ID set:",input$enrichment_id_select),type="message") })

  # ---- BULK LOAD ONTOLOGY TABLES ----
  observeEvent(input$bulk_load_btn, {
    req(input$bulk_ontology_file)
    withProgress(message="Loading ontology tables...", value=0, {
      tryCatch({
        incProgress(0.3)

        # Load all tables from the R object
        bulk_result <- load_bulk_ontology_tables(input$bulk_ontology_file$datapath)

        # Update reactive values
        if (!is.null(bulk_result$uniprot)) {
          uniprot_table(bulk_result$uniprot)
          assign("UniProtTable", bulk_result$uniprot, envir=.GlobalEnv)
          incProgress(0.15)
        }
        if (!is.null(bulk_result$go)) {
          go_table(bulk_result$go)
          incProgress(0.15)
        }
        if (!is.null(bulk_result$kegg)) {
          kegg_table(bulk_result$kegg)
          incProgress(0.15)
        }
        if (!is.null(bulk_result$reactome)) {
          reactome_table(bulk_result$reactome)
          incProgress(0.15)
        }

        # Show status
        if (length(bulk_result$loaded_names) > 0) {
          msg <- paste("Loaded:", paste(bulk_result$loaded_names, collapse=", "))
          showNotification(msg, type="message", duration=5)
        } else {
          showNotification("No recognized tables found in the file", type="warning", duration=5)
        }

        incProgress(0.1)
      }, error=function(e) {
        showNotification(paste("Error:", e$message), type="error", duration=7)
      })
    })
  })

  output$bulk_load_status <- renderUI({
    uniprot_ok <- !is.null(uniprot_table())
    go_ok <- !is.null(go_table())
    kegg_ok <- !is.null(kegg_table())
    reactome_ok <- !is.null(reactome_table())

    tags$div(style="margin-top:10px;padding:10px;background:#f5f5f5;border-radius:4px;",
      tags$span(icon(if (uniprot_ok) "check" else "times")," UniProt",
               style=paste("color:",if (uniprot_ok) "#27ae60" else "#c0392b",";margin-right:15px;")),
      tags$span(icon(if (go_ok) "check" else "times")," GO",
               style=paste("color:",if (go_ok) "#27ae60" else "#c0392b",";margin-right:15px;")),
      tags$span(icon(if (kegg_ok) "check" else "times")," KEGG",
               style=paste("color:",if (kegg_ok) "#27ae60" else "#c0392b",";margin-right:15px;")),
      tags$span(icon(if (reactome_ok) "check" else "times")," Reactome",
               style=paste("color:",if (reactome_ok) "#27ae60" else "#c0392b",";"))
    )
  })

  observeEvent(input$uniprot_table_load_btn, {
    req(input$uniprot_table_file)
    tryCatch({ t <- load_ontology_file(input$uniprot_table_file$datapath,input$uniprot_table_file$name); uniprot_table(t); assign("UniProtTable",t,envir=.GlobalEnv)
    showNotification(paste("UniProtTable:",nrow(t),"rows,",ncol(t),"cols"),type="message")
    }, error=function(e) showNotification(paste("Error:",e$message),type="error"))
  })
  output$uniprot_table_status <- renderUI({
    if (is.null(uniprot_table())) tags$div(style="margin-top:15px;",tags$span(icon("times-circle")," Not loaded - Required!",style="color:#c0392b;font-weight:bold;font-size:14px;"))
    else { t <- uniprot_table(); tags$div(style="margin-top:15px;",tags$span(icon("check-circle"),paste(" Loaded:",nrow(t),"rows,",ncol(t),"cols"),style="color:#27ae60;font-weight:bold;")) }
  })
  output$uniprot_table_preview <- renderDT({ req(uniprot_table()); datatable(head(uniprot_table(),5),options=list(dom="t",scrollX=TRUE,pageLength=5),rownames=FALSE,style="compact") })
  observeEvent(input$go_load_btn, { req(input$go_file); tryCatch({ t <- load_ontology_file(input$go_file$datapath,input$go_file$name); go_table(t); showNotification(paste("GO:",nrow(t)),type="message") }, error=function(e) showNotification(paste("Error:",e$message),type="error")) })
  observeEvent(input$kegg_load_btn, { req(input$kegg_file); tryCatch({ t <- load_ontology_file(input$kegg_file$datapath,input$kegg_file$name); kegg_table(t); showNotification(paste("KEGG:",nrow(t)),type="message") }, error=function(e) showNotification(paste("Error:",e$message),type="error")) })
  observeEvent(input$reactome_load_btn, { req(input$reactome_file); tryCatch({ t <- load_ontology_file(input$reactome_file$datapath,input$reactome_file$name); reactome_table(t); showNotification(paste("Reactome:",nrow(t)),type="message") }, error=function(e) showNotification(paste("Error:",e$message),type="error")) })
  output$go_status_indicator <- renderUI({ if (is.null(go_table())) tags$span(icon("times-circle")," Not loaded",style="color:#c0392b;font-weight:bold;") else tags$span(icon("check-circle"),paste(" Loaded:",nrow(go_table())),style="color:#27ae60;font-weight:bold;") })
  output$kegg_status_indicator <- renderUI({ if (is.null(kegg_table())) tags$span(icon("times-circle")," Not loaded",style="color:#c0392b;font-weight:bold;") else tags$span(icon("check-circle"),paste(" Loaded:",nrow(kegg_table())),style="color:#27ae60;font-weight:bold;") })
  output$reactome_status_indicator <- renderUI({ if (is.null(reactome_table())) tags$span(icon("times-circle")," Not loaded",style="color:#c0392b;font-weight:bold;") else tags$span(icon("check-circle"),paste(" Loaded:",nrow(reactome_table())),style="color:#27ae60;font-weight:bold;") })
  output$go_terms_preview <- renderUI({ ontology_term_preview(go_table(),"^GO_description$",c("GO_Description","description","term","name")) })
  output$kegg_terms_preview <- renderUI({ ontology_term_preview(kegg_table(),"^Pathway_name$",c("pathway_name","pathway","name","description")) })
  output$reactome_terms_preview <- renderUI({ ontology_term_preview(reactome_table(),"^REACTOME_description$",c("Reactome_description","reactome","description","name")) })
  output$ontology_summary <- renderPrint({
    cat("=== Ontology Summary ===\n\n")
    if (!is.null(enrichment_id_col())) cat("ID:",enrichment_id_col(),"\n") else cat("ID: NOT SET\n")
    cat("\n--- Reference ---\n")
    if (!is.null(uniprot_table())) cat("UniProtTable: [LOADED]",nrow(uniprot_table()),"rows\n") else cat("UniProtTable: [NOT LOADED] ** Required **\n")
    cat("\n--- Ontology Tables ---\n")
    for (nm in list(list("GO",go_table()),list("KEGG",kegg_table()),list("Reactome",reactome_table()))) if (!is.null(nm[[2]])) cat(nm[[1]],": [LOADED]",nrow(nm[[2]]),"\n") else cat(nm[[1]],": [NOT LOADED]\n")
    nl <- sum(!is.null(go_table()),!is.null(kegg_table()),!is.null(reactome_table())); hu <- !is.null(uniprot_table()); hi <- !is.null(enrichment_id_col())
    cat("\n--- Readiness ---\n"); if (nl>0&&hu&&hi) cat(">> READY\n") else { cat(">> NOT READY\n"); if (!hi) cat("   - Set enrichment ID\n"); if (!hu) cat("   - Load UniProtTable\n"); if (nl==0) cat("   - Load ontology table(s)\n") }
  })

  # ---- BULK LOAD PREVIEW OUTPUTS ----
  output$bulk_uniprot_status <- renderUI({
    if (is.null(uniprot_table())) tags$span(icon("times-circle")," Not loaded",style="color:#c0392b;font-weight:bold;font-size:12px;")
    else tags$span(icon("check-circle"),paste(" Loaded:",nrow(uniprot_table()),"rows"),style="color:#27ae60;font-weight:bold;font-size:12px;")
  })
  output$bulk_go_status <- renderUI({
    if (is.null(go_table())) tags$span(icon("times-circle")," Not loaded",style="color:#c0392b;font-weight:bold;font-size:12px;")
    else tags$span(icon("check-circle"),paste(" Loaded:",nrow(go_table()),"rows"),style="color:#27ae60;font-weight:bold;font-size:12px;")
  })
  output$bulk_kegg_status <- renderUI({
    if (is.null(kegg_table())) tags$span(icon("times-circle")," Not loaded",style="color:#c0392b;font-weight:bold;font-size:12px;")
    else tags$span(icon("check-circle"),paste(" Loaded:",nrow(kegg_table()),"rows"),style="color:#27ae60;font-weight:bold;font-size:12px;")
  })
  output$bulk_reactome_status <- renderUI({
    if (is.null(reactome_table())) tags$span(icon("times-circle")," Not loaded",style="color:#c0392b;font-weight:bold;font-size:12px;")
    else tags$span(icon("check-circle"),paste(" Loaded:",nrow(reactome_table()),"rows"),style="color:#27ae60;font-weight:bold;font-size:12px;")
  })
  output$bulk_uniprot_preview <- renderUI({ ontology_term_preview(uniprot_table(),"^Uniprot_ID$|^Uniprot_Accession$",c("Uniprot_ID","Uniprot_Accession","ID","Accession")) })
  output$bulk_go_preview <- renderUI({ ontology_term_preview(go_table(),"^GO_description$|^GO_accession$",c("GO_description","GO_accession","description","accession")) })
  output$bulk_kegg_preview <- renderUI({ ontology_term_preview(kegg_table(),"^Pathway_name$|^Pathway_ID$",c("Pathway_name","Pathway_ID","pathway_name","pathway_id")) })
  output$bulk_reactome_preview <- renderUI({ ontology_term_preview(reactome_table(),"^REACTOME_description$",c("REACTOME_description","reactome","description")) })
  
  # ============================================================================
  # 8.3 EMBEDDINGS & DIMENSIONALITY REDUCTION - PCA, UMAP, t-SNE
  # ============================================================================

  # Update available embedding types when object changes
  observe({
    req(romics_object())
    embeddings <- check_embeddings_available(romics_object())

    # Build list of available embedding types
    embedding_choices <- list()
    if (embeddings$pca) embedding_choices[["PCA"]] <- "pca"
    if (embeddings$umap) embedding_choices[["UMAP"]] <- "umap"
    if (embeddings$tsne) embedding_choices[["t-SNE"]] <- "tsne"

    # Update UI with available embeddings
    if (length(embedding_choices) > 0) {
      updateSelectInput(session, "grouping_type", choices=embedding_choices, selected=embedding_choices[[1]])
    } else {
      updateSelectInput(session, "grouping_type", choices=c("No embeddings"="none"))
    }
  })
  output$plot_dimension_ui <- renderUI({
    req(romics_object())
    emb <- check_embeddings_available(romics_object())
    
    # Get dimensions for selected embedding
    if (!is.null(input$grouping_type) && input$grouping_type != "none") {
      md <- get_embedding_dimensions(romics_object(), input$grouping_type)
    } else {
      md <- 0
    }
    
    # If less than 3 dimensions, disable 3D option
    if (md < 3) {
      radioButtons("plot_dimension", "Dim:", choices = c("2D" = "2d"), selected = "2d", inline = TRUE)
    } else {
      radioButtons("plot_dimension", "Dim:", choices = c("2D" = "2d", "3D" = "3d"), selected = "2d", inline = TRUE)
    }
  })
  
  observe({ req(romics_object(),input$grouping_type); if (input$grouping_type=="none") return(); md <- get_embedding_dimensions(romics_object(),input$grouping_type)
  if (md > 0) {
    # Update 2D component selectors
    updateNumericInput(session, "x_comp", value=1, min=1, max=md)
    updateNumericInput(session, "y_comp", value=min(2, md), min=1, max=md)

    # Update 3D component selectors
    updateNumericInput(session, "x_comp_3d", value=1, min=1, max=md)
    updateNumericInput(session, "y_comp_3d", value=min(2, md), min=1, max=md)
    updateNumericInput(session, "z_comp_3d", value=min(3, md), min=1, max=md)
  }
})
  observe({ req(romics_object()); fns <- get_clean_factors(romics_object()); if (length(fns)>0) { sel <- if (!is.null(current_factor())&&current_factor() %in% fns) current_factor() else fns[1]; updateSelectInput(session,"grouping_factor_select",choices=fns,selected=sel) } })
  observe({ req(romics_object()); updateSelectizeInput(session,"grouping_feature_select",choices=rownames(romics_object()$data),server=TRUE) })
  output$grouping_embedding_info <- renderPrint({ req(romics_object()); emb <- check_embeddings_available(romics_object()); st <- function(a,t) if (a) paste0("[YES] ",get_embedding_dimensions(romics_object(),t)," comp") else "[NO]"; cat("PCA:",st(emb$pca,"pca"),"\nUMAP:",st(emb$umap,"umap"),"\nt-SNE:",st(emb$tsne,"tsne"),"\n") })
  output$grouping_plot_settings <- renderPrint({
    cat("Type:",toupper(input$grouping_type),"\n")
    cat("Dim:",input$plot_dimension,"\n")
    cat("Color:",input$grouping_color_by,"\n")
    cat("Size:",input$grouping_size,"\n")
    cat("Opacity:",input$grouping_alpha,"\n")
  })
  observeEvent(input$plot_grouping, {
    req(romics_object())
    if (input$grouping_type == "none") {
      showNotification("No embeddings.", type="error")
      return()
    }

    # Validate color input
    if (input$grouping_color_by == "factor") {
      req(input$grouping_factor_select)
    } else {
      req(input$grouping_feature_select)
    }

    obj <- romics_object()

    withProgress(message="Generating plot...", value=0, {
      tryCatch({
        incProgress(0.2)

        # Prepare plot parameters
        color_factor <- if (input$grouping_color_by == "factor") input$grouping_factor_select else NULL
        color_feature <- if (input$grouping_color_by == "feature") input$grouping_feature_select else NULL
        embedding_dims <- get_embedding_dimensions(obj, input$grouping_type)

        # Select plotting function
        plot_func <- switch(input$grouping_type,
          "pca" = RomicsProcessor::romicsPCAplot,
          "umap" = RomicsProcessor::romicsUmapPlot,
          "tsne" = RomicsProcessor::romicsTsnePlot,
          NULL
        )

        if (input$plot_dimension == "2d") {
          render_2d_embedding(obj, plot_func, input, color_factor, color_feature, embedding_dims, output)
        } else {
          render_3d_embedding(obj, plot_func, input, color_factor, color_feature, embedding_dims, output)
        }

        incProgress(1)
      }, error=function(e) {
        showNotification(paste("Error:", e$message), type="error")
        output$grouping_plot_2d <- renderPlot({
          plot.new()
          text(0.5, 0.5, e$message, col="red")
        })
        output$grouping_plot_3d <- renderPlotly({ plotly_empty() })
      })
    })
  })
  
  # ============================================================================
  # 8.4 SINGLE FEATURE ANALYSIS - Feature Plots & Statistics
  # ============================================================================

  observe({ req(romics_object()); fns <- get_clean_factors(romics_object()); if (length(fns)>0) { sel <- if (!is.null(current_factor())&&current_factor() %in% fns) current_factor() else fns[1]; updateSelectInput(session,"feature_factor",choices=fns,selected=sel) } })
  output$feature_plot_info <- renderPrint({ req(input$feature_select); cat("Feature:",input$feature_select,"\nType:",input$plot_type,"\nFactor:",input$feature_factor,"\n") })
  output$feature_stats_table <- renderDT({ 
    sd <- feature_stats_store()
    if (is.null(sd) || nrow(sd) == 0) {
      datatable(data.frame(Message = "Press Generate to load statistics"), options = list(dom = "t"))
    } else {
      datatable(sd, extensions = "Buttons", options = list(dom = "Bfrtip", buttons = c("copy", "csv"), pageLength = 10, scrollX = TRUE), rownames = FALSE)
    }
  })
  
  observeEvent(input$plot_feature, {
    req(romics_object(),input$feature_select,input$feature_factor)
    obj <- romics_object()
    
    tryCatch({
      pt <- if (input$custom_title=="") "auto" else input$custom_title
      bp <- RomicsProcessor::singleFeaturePlot(
        romics_object=obj,
        feature=input$feature_select,
        plot_type=input$plot_type,
        factor=input$feature_factor,
        title=pt
      )
      if (input$feature_plot_format=="plotly") {
        output$feature_plot_ui <- renderUI({ plotlyOutput("feature_plot_plotly",height="600px") })
        output$feature_plot_plotly <- renderPlotly({ ggplotly(bp+theme(aspect.ratio=0.75)) })
      } else {
        output$feature_plot_ui <- renderUI({ plotOutput("feature_plot_static",height="600px") })
        output$feature_plot_static <- renderPlot({ bp })
      }
      
      pval_type <- if (is.null(input$feature_pval_type)) "all" else input$feature_pval_type
      stats_data <- extract_feature_statistics(obj, input$feature_select, pval_type)
      feature_stats_store(stats_data)
    }, error=function(e) {
      output$feature_plot_ui <- renderUI({ plotOutput("fpe",height="600px") })
      output$fpe <- renderPlot({ plot.new(); text(.5,.5,paste("Error:",e$message),col="red",cex=0.8) })
      feature_stats_store(NULL)
    })
  })
  
  observeEvent(input$feature_pval_type, {
    req(romics_object(),input$feature_select)
    pval_type <- input$feature_pval_type
    stats_data <- tryCatch({
      extract_feature_statistics(romics_object(), input$feature_select, pval_type)
    }, error = function(e) {
      NULL
    })
    feature_stats_store(stats_data)
  })
  
  # ---- BUBBLE ----
  observe({ req(romics_object()); fns <- get_clean_factors(romics_object()); if (length(fns)>0) { sel <- if (!is.null(current_factor())&&current_factor() %in% fns) current_factor() else fns[1]; updateSelectInput(session,"multi_factor_select",choices=fns,selected=sel) } })
  observe({ req(romics_object()); updateSelectizeInput(session,"multi_feature_dropdown",choices=rownames(romics_object()$data),server=TRUE) })
  observeEvent(input$add_feature_btn, { req(input$multi_feature_dropdown); ct <- input$multi_feature_input; nf <- input$multi_feature_dropdown
  updateTextInput(session,"multi_feature_input",value=if (is.null(ct)||ct=="") nf else paste(ct,nf,sep=", ")); updateSelectizeInput(session,"multi_feature_dropdown",selected="") })
  observeEvent(input$plot_multi, {
    req(romics_object(),input$multi_feature_input,input$multi_factor_select); obj <- romics_object()
    withProgress(message="Generating...", value=0, { tryCatch({
      fl <- trimws(strsplit(input$multi_feature_input,",")[[1]]); fl <- fl[fl!=""]; fn <- input$multi_factor_select
      rc <- obj; if (!is.null(rc$main_factor)&&fn!=rc$main_factor) rc <- RomicsProcessor::romicsChangeFactor(rc,main_factor=fn)
      rc <- rc[!names(rc)=="statistics"]; class(rc) <- "romics_object"; rc <- RomicsProcessor::romicsMean(rc); rc <- RomicsProcessor::romicsPercentComplete(rc)
      df <- rc$statistics; vf <- fl[fl %in% rownames(df)]; if (!length(vf)) { showNotification("No valid!",type="error"); return() }
      df <- df[rownames(df) %in% vf,]; df$Feature <- rownames(df); rownames(df) <- NULL
      dm <- df[grepl("_mean|Feature",colnames(df))]; if (input$multi_scale_features) { mc <- grepl("_mean",colnames(dm)); if (sum(mc)>0) dm[,mc] <- t(scale(t(dm[,mc]))) }
      colnames(dm) <- gsub("_mean","",colnames(dm)); dc <- df[grepl("_percentage_completeness|Feature",colnames(df))]; colnames(dc) <- gsub("_percentage_completeness","",colnames(dc))
      d <- data.frame(); for (g in setdiff(colnames(dm),"Feature")) d <- rbind(d,data.frame(Feature=dm$Feature,Group=g,Mean=dm[[g]],Completeness=if (g %in% colnames(dc)) dc[[g]] else NA,stringsAsFactors=FALSE))
      d$Feature <- factor(d$Feature,levels=rev(vf)); cl <- if (input$multi_scale_features) "Scaled Mean" else "Mean"
      p <- ggplot(d,aes(x=Group,y=Feature))+geom_point(aes(size=Completeness,color=Mean),alpha=.85)+scale_color_viridis_c(option=input$multi_viridis_option,name=cl)+scale_size_continuous(name="% Complete",range=c(2,12))+theme_bw()+theme(axis.text.x=element_text(angle=45,hjust=1))+labs(title=paste("Bubble -",fn),x=fn,y="Feature")
      multi_plot_store(p); output$multi_bubble_plot <- renderPlotly({ ggplotly(p,tooltip=c("x","y","colour","size"))%>%layout(margin=list(b=100)) })
    }, error=function(e) { showNotification(paste("Error:",e$message),type="error"); output$multi_bubble_plot <- renderPlotly({ plotly_empty() }) }) })
  })
  output$download_bubble_plot_png <- downloadHandler(filename=function() paste0("bubble_",Sys.Date(),".png"), content=function(file) { req(multi_plot_store()); ggsave(file,plot=multi_plot_store(),width=10,height=8,dpi=300,bg="white") })
  
  # ============================================================================
  # 8.5 VOLCANO PLOTS - Log Fold Change vs P-value Visualization
  # ============================================================================

  observe({ req(romics_object()); fns <- get_clean_factors(romics_object()); if (length(fns)>0) { sel <- if (!is.null(current_factor())&&current_factor() %in% fns) current_factor() else fns[1]; updateSelectInput(session,"volcano_factor_select",choices=fns,selected=sel) } })
  observe({ req(romics_object(),input$volcano_factor_select); tryCatch({ c <- extract_comparisons_for_factor(romics_object(),input$volcano_factor_select); if (is.null(c)||!length(c)) c <- extract_all_comparisons(romics_object()); volcano_comparisons(c) }, error=function(e) volcano_comparisons(NULL)) })
  observe({ req(romics_object(),input$volcano_comparison_select); tryCatch({ volcano_pvalue_types(check_pvalue_types(romics_object(),input$volcano_comparison_select)) }, error=function(e) volcano_pvalue_types(list(p=FALSE,padj=FALSE))) })
  output$volcano_comparison_ui <- renderUI({ c <- volcano_comparisons(); if (is.null(c)||!length(c)) return(HTML("<p style='color:#c0392b;'>No comparisons.</p>")); selectInput("volcano_comparison_select","Comparison:",choices=c) })
  output$volcano_pvalue_type_ui <- renderUI({ pt <- volcano_pvalue_types(); ch <- list(); if (pt$padj) ch[["Adj p"]] <- "padj"; if (pt$p) ch[["p"]] <- "p"; if (!length(ch)) return(NULL); selectInput("volcano_pvalue_type","p-value:",choices=ch,selected=if (pt$padj) "padj" else "p") })
  output$volcano_pval_threshold_ui <- renderUI({ req(input$volcano_pvalue_type); numericInput("volcano_pval_threshold",if (input$volcano_pvalue_type=="padj") "Max adj p:" else "Max p:",value=0.05,min=0.001,max=1,step=0.01) })
  observeEvent(input$plot_volcano, {
    req(romics_object(),input$volcano_comparison_select,input$volcano_pvalue_type,input$volcano_pval_threshold,input$volcano_size,input$volcano_alpha)
    obj <- romics_object(); cn <- input$volcano_comparison_select; pvt <- input$volcano_pvalue_type; pvth <- input$volcano_pval_threshold; fcth <- input$volcano_fc_threshold; v_size <- input$volcano_size; v_alpha <- input$volcano_alpha
    withProgress(message="Generating...", value=0, { tryCatch({
      vd <- extract_volcano_data(obj,cn,pvt); if (is.null(vd)) stop("No data for: ",cn)
      vd$Significant <- vd$p_value<pvth & abs(vd$log_FC)>fcth; vd$Color <- ifelse(vd$Significant,"Sig","NS"); volcano_data_store(vd)
      p <- ggplot(vd,aes(x=log_FC,y=neg_log10_p,color=Color,text=paste("Feature:",Feature,"<br>FC:",round(log_FC,3),"<br>p:",format(p_value,scientific=TRUE,digits=3))))+geom_point(alpha=v_alpha,size=v_size)+scale_color_manual(values=c("Sig"="#e74c3c","NS"="#95a5a6"))+geom_hline(yintercept=-log10(pvth),linetype="dashed",color="gray")+geom_vline(xintercept=c(-fcth,fcth),linetype="dashed",color="gray")+labs(title=paste("Volcano:",cn),x="Log(FC)",y=paste0("-Log10(",if (pvt=="padj") "adj p" else "p",")"))+theme_minimal()
      output$volcano_plot_ui <- renderUI({ plotlyOutput("volcano_plot_plotly",height="700px") }); output$volcano_plot_plotly <- renderPlotly({ ggplotly(p,tooltip="text")%>%layout(height=700,hovermode="closest") })
      sc <- sum(vd$Significant); output$volcano_plot_info <- renderPrint({ cat("Comparison:",cn,"\nSignificant:",sc,"/",nrow(vd),"\nSize:",v_size,"Alpha:",v_alpha,"\n"); if (sc>0) { ts <- vd[vd$Significant,]; ts <- ts[order(ts$p_value),][1:min(10,nrow(ts)),]; cat("\nTop:\n"); for (i in 1:nrow(ts)) cat(" ",ts$Feature[i],"| FC:",round(ts$log_FC[i],3),"\n") } })
    }, error=function(e) { showNotification(paste("Error:",e$message),type="error"); output$volcano_plot_ui <- renderUI({ helpText(e$message) }) }) })
  })
  output$volcano_extract_summary <- renderPrint({ vd <- volcano_data_store(); if (is.null(vd)) { cat("Generate volcano first.\n"); return() }; et <- input$volcano_extract_type
  sel <- if (et=="all_sig") vd[vd$Significant,] else if (et=="up") vd[vd$Significant&vd$log_FC>0,] else vd[vd$Significant&vd$log_FC<0,]
  cat(switch(et,"all_sig"="All significant","up"="Up-regulated","down"="Down-regulated"),":",nrow(sel),"\n") })
  observeEvent(input$volcano_copy_ids, { vd <- volcano_data_store(); if (is.null(vd)) return(); et <- input$volcano_extract_type
  sel <- if (et=="all_sig") vd[vd$Significant,] else if (et=="up") vd[vd$Significant&vd$log_FC>0,] else vd[vd$Significant&vd$log_FC<0,]
  if (!nrow(sel)) { showNotification("No features.",type="warning"); return() }; session$sendCustomMessage("copyClipboard",list(text=paste(sel$Feature,collapse=", "))); showNotification(paste(nrow(sel),"IDs copied!"),type="message") })
  output$volcano_download_ids <- downloadHandler(filename=function() paste0("volcano_",input$volcano_extract_type,"_",Sys.Date(),".csv"),
                                                 content=function(file) { vd <- volcano_data_store(); if (is.null(vd)) return(); et <- input$volcano_extract_type
                                                 sel <- if (et=="all_sig") vd[vd$Significant,] else if (et=="up") vd[vd$Significant&vd$log_FC>0,] else vd[vd$Significant&vd$log_FC<0,]
                                                 write.csv(sel[,c("Feature","log_FC","p_value")],file,row.names=FALSE) })
  
  # ============================================================================
  # 8.6 HEATMAPS - Complex Feature Expression Heatmaps with Dendrograms
  # ============================================================================

  observe({ req(romics_object()); obj <- romics_object(); ch <- build_feature_choices(obj,include_all=TRUE,include_anova=TRUE,include_clusters=TRUE,heatmap_clusters_data=heatmap_clusters())
  updateSelectizeInput(session,"hm_feature_dropdown",choices=ch,server=TRUE) })
  observeEvent(input$hm_add_feature_btn, { req(input$hm_feature_dropdown); sel <- input$hm_feature_dropdown; obj <- romics_object()
  resolved <- resolve_dropdown_selection(sel,obj,heatmap_clusters())
  if (length(resolved)>1||resolved[1]!=sel) { updateTextInput(session,"hm_feature_input",value=paste(resolved,collapse=", ")); showNotification(paste(length(resolved),"features added"),type="message") }
  else { ct <- input$hm_feature_input; updateTextInput(session,"hm_feature_input",value=if (is.null(ct)||ct=="") sel else paste(ct,sel,sep=", ")) }
  updateSelectizeInput(session,"hm_feature_dropdown",selected="") })
  observe({ req(romics_object()); fns <- get_clean_factors(romics_object()); if (length(fns)>0) { sel <- if (!is.null(current_factor())&&current_factor() %in% fns) current_factor() else fns[1]; updateSelectInput(session,"hm_factor_select",choices=fns,selected=sel) } })
  observe({ req(input$hm_show_clusters,input$hm_n_clusters); nc <- input$hm_n_clusters; updateSelectInput(session,"hm_cluster_view",choices=c("All Clusters"="all",setNames(as.list(1:nc),paste("Cluster",1:nc))),selected="all") })
  observe({ req(romics_object()); obj <- romics_object(); asc <- character(0)
  tryCatch({ asc <- RomicsProcessor::romicsCalculatedStats(obj); asc <- asc[!grepl("_mean$|_sd$|^Z_scores_|_percentage_completeness$",asc)] },
           error=function(e) { if (!is.null(obj$statistics)) asc <- colnames(obj$statistics)[!grepl("_mean$|_sd$|^Z_scores_|_percentage_completeness$",colnames(obj$statistics))] })
  anova_cols <- asc[grepl("^ANOVA_",asc)]; other_cols <- setdiff(asc,anova_cols)
  ch <- c("None"="none"); if (length(anova_cols)>0) ch <- c(ch,setNames(anova_cols,paste0("[ANOVA] ",anova_cols))); if (length(other_cols)>0) ch <- c(ch,setNames(other_cols,other_cols))
  updateSelectInput(session,"hm_filter_stat_column",choices=ch) })
  observeEvent(input$hm_filter_stat_column, { if (is.null(input$hm_filter_stat_column)||input$hm_filter_stat_column=="none") return()
    col <- input$hm_filter_stat_column; if (grepl("_p$|_padj$",col)) updateTextInput(session,"hm_stat_column_filter",value="<0.05")
    else if (grepl("log\\(|log2FC|fold_change",col,ignore.case=TRUE)) updateTextInput(session,"hm_stat_column_filter",value=">1") })
  output$hm_anova_filter_ui <- renderUI({ req(romics_object()); obj <- romics_object(); if (is.null(obj$statistics)) return(NULL)
  ac <- colnames(obj$statistics); ap <- ac[grepl("^ANOVA_",ac)&grepl("_p$|_padj$",ac)]; if (!length(ap)) return(NULL)
  tagList(div(style="background:#fef9e7;border:1px solid #f9e79f;border-radius:4px;padding:10px;",
              selectInput("hm_anova_stat_col","ANOVA Statistic:",choices=setNames(ap,ap),selected=ap[1]),
              numericInput("hm_anova_threshold","Threshold:",value=0.05,min=0.001,max=1,step=0.01),
              actionButton("hm_apply_anova_filter","Apply ANOVA Filter",class="btn-warning",icon=icon("bolt"),width="100%"))) })
  observeEvent(input$hm_apply_anova_filter, { req(romics_object(),input$hm_anova_stat_col,input$hm_anova_threshold); obj <- romics_object()
  if (is.null(obj$statistics)) { showNotification("No statistics!",type="error"); return() }; cn <- input$hm_anova_stat_col; th <- input$hm_anova_threshold
  vals <- as.numeric(obj$statistics[,cn]); keep <- !is.na(vals)&vals<th; ff <- rownames(obj$statistics)[keep]
  updateTextInput(session,"hm_feature_input",value=paste(ff,collapse=", ")); showNotification(paste("ANOVA:",length(ff),"features"),type="message") })
  output$hm_selected_features <- renderPrint({ ft <- input$hm_feature_input; if (is.null(ft)||ft=="") { cat("None.\n"); return() }; fl <- trimws(strsplit(ft,",")[[1]]); fl <- fl[fl!=""]; cat("Total:",length(fl),"\n") })
  output$hm_plot_info <- renderPrint({ req(romics_object()); cat("Factor:",input$hm_factor_select,"\nScale:",input$hm_scale_feature,"\nClusters:",input$hm_show_clusters,if (input$hm_show_clusters) paste0(" (n=",input$hm_n_clusters,")") else "","\nFilter:",input$hm_filter_stat_column,"\n") })
  observeEvent(input$plot_heatmap, {
    req(romics_object(),input$hm_feature_input,input$hm_factor_select); obj <- romics_object()
    withProgress(message="Generating heatmap...", value=0, { tryCatch({
      ft <- input$hm_feature_input; fl <- if (ft=="All Features") rownames(obj$data) else { x <- trimws(strsplit(ft,",")[[1]]); x[x!=""] }
      if (!length(fl)) { showNotification("No features!",type="error"); return() }; vf <- fl[fl %in% rownames(obj$data)]; if (!length(vf)) { showNotification("No valid!",type="error"); return() }
      nc <- if (input$hm_show_clusters) input$hm_n_clusters else NULL; sl <- if (input$hm_show_clusters) input$hm_show_cluster_legend else FALSE
      heatmap_clusters(NULL); if (exists("heatmapFeatureClust",envir=.GlobalEnv)) rm("heatmapFeatureClust",envir=.GlobalEnv); incProgress(0.3)
      po <- heatmapFeatures(romics_object=obj,factor=input$hm_factor_select,scale_feature=input$hm_scale_feature,feature_list=vf,
                            filter_by_stat_column=if (input$hm_filter_stat_column!="none") input$hm_filter_stat_column else NULL,
                            stat_column_filter=if (input$hm_filter_stat_column!="none") input$hm_stat_column_filter else NULL,
                            viridis_option=input$hm_viridis_option,show_completeness=input$hm_show_completeness,show_dendrogram=input$hm_show_dendrogram,
                            show_feature_names=input$hm_show_feature_names,show_clusters=input$hm_show_clusters,show_cluster_legend=sl,n_clusters=nc,clustering_method=input$hm_clustering_method)
      heatmap_plot_store(po); incProgress(0.4); ph <- max(400,min(length(vf)*8+150,1200))
      output$heatmap_plot_ui <- renderUI({ plotOutput("heatmap_plot_render",height=paste0(ph,"px"),width="100%") }); output$heatmap_plot_render <- renderPlot({ po },res=96)
      if (input$hm_show_clusters&&exists("heatmapFeatureClust",envir=.GlobalEnv)) { cd <- get("heatmapFeatureClust",envir=.GlobalEnv); heatmap_clusters(cd)
      updateSelectInput(session,"hm_cluster_view",choices=c("All Clusters"="all",setNames(as.list(1:cd$n_clusters),paste("Cluster",1:cd$n_clusters))),selected="all") }
      showNotification(paste("Heatmap:",length(vf),"features"),type="message")
    }, error=function(e) { showNotification(paste("Error:",e$message),type="error"); output$heatmap_plot_ui <- renderUI({ helpText(e$message) }) }) })
  })
  output$download_heatmap_png <- downloadHandler(filename=function() paste0("heatmap_",Sys.Date(),".png"), content=function(file) { req(heatmap_plot_store()); nf <- length(trimws(strsplit(input$hm_feature_input,",")[[1]])); ggsave(file,plot=heatmap_plot_store(),width=12,height=max(8,min(nf*0.15+3,30)),dpi=300,bg="white") })
  output$hm_cluster_table <- renderDT({ cd <- heatmap_clusters(); if (is.null(cd)) return(NULL); ca <- cd$cluster_assignments; sel <- input$hm_cluster_view
  if (!is.null(sel)&&sel!="all") ca <- ca[ca$Cluster==as.numeric(sel),,drop=FALSE]
  datatable(ca,extensions="Buttons",options=list(dom="Bfrtip",buttons=c("copy","csv","excel"),pageLength=15,scrollX=TRUE,scrollY="300px"),rownames=FALSE,filter="top",colnames=c("Feature","Cluster")) })
  output$download_cluster_csv <- downloadHandler(filename=function() paste0("clusters_",Sys.Date(),".csv"), content=function(file) { cd <- heatmap_clusters(); if (!is.null(cd)&&!is.null(cd$cluster_assignments)) write.csv(cd$cluster_assignments,file,row.names=FALSE) })
  observeEvent(input$copy_cluster_ids, { cd <- heatmap_clusters(); if (is.null(cd)) return(); ca <- cd$cluster_assignments; sel <- input$hm_cluster_view
  if (!is.null(sel)&&sel!="all") ca <- ca[ca$Cluster==as.numeric(sel),,drop=FALSE]
  session$sendCustomMessage("copyClipboard",list(text=paste(ca$Feature,collapse=", "))); showNotification("IDs copied!",type="message") })
  # ---- STATISTICS TAB ----
  observe({
    req(romics_object()); if (is.null(romics_object()$statistics)) return()
    sc <- colnames(romics_object()$statistics)
    sc_filtered <- sc[!grepl("_sd$|_mean$|^Z_scores_|_percentage_completeness$",sc)]
    ch <- c("None"="none", setNames(sc_filtered,sc_filtered))
    for (i in 1:6) updateSelectInput(session, paste0("stat_filter_col_",i), choices=ch)
  })
  observe({
    for (i in 1:6) { local({ ii <- i
    observeEvent(input[[paste0("stat_filter_col_",ii)]], {
      col <- input[[paste0("stat_filter_col_",ii)]]; if (is.null(col)||col=="none") return()
      if (grepl("_p$|_padj$",col)) updateTextInput(session,paste0("stat_filter_expr_",ii),value="<0.05")
      else if (grepl("log\\(|fold_change",col,ignore.case=TRUE)) updateTextInput(session,paste0("stat_filter_expr_",ii),value=">1")
      updateCheckboxInput(session,paste0("stat_filter_active_",ii),value=TRUE)
    }, ignoreInit=TRUE) })
    }
  })
  observeEvent(input$apply_stat_filters, {
    req(romics_object()); if (is.null(romics_object()$statistics)) return()
    sdf <- romics_object()$statistics; ck <- colnames(sdf)[!grepl("_sd$|_mean$|^Z_scores_|_percentage_completeness$",colnames(sdf))]
    sdf_display <- sdf[,ck,drop=FALSE]; selected_cols <- c(); keep <- rep(TRUE,nrow(sdf_display))
    for (i in 1:6) { col <- input[[paste0("stat_filter_col_",i)]]; fe <- input[[paste0("stat_filter_expr_",i)]]; active <- input[[paste0("stat_filter_active_",i)]]
    if (!is.null(col)&&col!="none") { selected_cols <- c(selected_cols,col)
    if (!is.null(active)&&active&&!is.null(fe)&&fe!=""&&col %in% colnames(sdf_display)) keep <- keep & apply_stat_filter(as.numeric(sdf_display[,col]),fe) } }
    if (length(selected_cols)>0) { vc <- selected_cols[selected_cols %in% colnames(sdf_display)]; if (length(vc)>0) sdf_display <- sdf_display[,vc,drop=FALSE] }
    sdf_display <- sdf_display[keep,,drop=FALSE]; stats_filtered_data(sdf_display)
    showNotification(paste("Filtered:",nrow(sdf_display),"features,",ncol(sdf_display),"columns"),type="message")
  })
  observeEvent(input$reset_stat_filters, { stats_filtered_data(NULL)
    for (i in 1:6) { updateSelectInput(session,paste0("stat_filter_col_",i),selected="none"); updateTextInput(session,paste0("stat_filter_expr_",i),value=""); updateCheckboxInput(session,paste0("stat_filter_active_",i),value=FALSE) }
    showNotification("Filters reset.",type="message") })
  output$stat_filter_summary <- renderPrint({ fd <- stats_filtered_data(); if (is.null(fd)) cat("No filters applied.\n") else cat("Filtered:",nrow(fd),"features\n",ncol(fd),"columns") })
  output$stats_table <- renderDT({
    req(romics_object()); if (is.null(romics_object()$statistics)) return(datatable(data.frame(Message="No statistics.")))
    fd <- stats_filtered_data()
    if (!is.null(fd)) datatable(fd,extensions="Buttons",options=list(dom="Bfrtip",buttons=c("copy","csv","excel"),pageLength=25,scrollX=TRUE),filter="top",rownames=TRUE)
    else { sdf <- romics_object()$statistics; ck <- colnames(sdf)[!grepl("_sd$|_mean$|^Z_scores_|_percentage_completeness$",colnames(sdf))]
    datatable(sdf[,ck,drop=FALSE],extensions="Buttons",options=list(dom="Bfrtip",buttons=c("copy","csv","excel"),pageLength=25,scrollX=TRUE),filter="top",rownames=TRUE) }
  })

  # Display table information (full vs. shown rows)
  output$stats_table_info <- renderPrint({
    req(romics_object())
    total <- nrow(romics_object()$statistics)
    fd <- stats_filtered_data()
    if (!is.null(fd)) {
      shown <- nrow(fd)
      pct <- round(shown/total*100, 1)
      cat(sprintf("Showing %d of %d\n(%.1f%%)", shown, total, pct))
    } else {
      cat(sprintf("Showing all\n%d rows", total))
    }
  })

  output$download_stats <- downloadHandler(
    filename = function() paste0("stats_full_", Sys.Date(), ".csv"),
    content = function(file) {
      req(romics_object())
      if (!is.null(romics_object()$statistics)) {
        write.csv(romics_object()$statistics, file, row.names = TRUE)
        showNotification(paste("Downloaded", nrow(romics_object()$statistics), "features"), type = "message")
      }
    }
  )

  output$download_stats_filtered <- downloadHandler(
    filename = function() paste0("stats_filtered_", Sys.Date(), ".csv"),
    content = function(file) {
      fd <- stats_filtered_data()
      if (!is.null(fd)) {
        write.csv(fd, file, row.names = TRUE)
        showNotification(paste("Downloaded", nrow(fd), "shown features"), type = "message")
      } else {
        req(romics_object())
        if (!is.null(romics_object()$statistics)) {
          write.csv(romics_object()$statistics, file, row.names = TRUE)
          showNotification(paste("Downloaded", nrow(romics_object()$statistics), "features (no filters applied)"), type = "message")
        }
      }
    }
  )

  observeEvent(input$copy_stats_ids, {
    req(romics_object())
    if (is.null(romics_object()$statistics)) return()
    ids <- rownames(romics_object()$statistics)
    session$sendCustomMessage("copyClipboard", list(text = paste(ids, collapse = ", ")))
    showNotification(paste("All", length(ids), "IDs copied!"), type = "message")
  })

  observeEvent(input$copy_stats_filtered_ids, {
    fd <- stats_filtered_data()
    ids <- if (!is.null(fd)) {
      rownames(fd)
    } else {
      req(romics_object())
      if (!is.null(romics_object()$statistics)) rownames(romics_object()$statistics) else return()
    }
    session$sendCustomMessage("copyClipboard", list(text = paste(ids, collapse = ", ")))
    shown_label <- if (!is.null(stats_filtered_data())) "shown" else ""
    showNotification(paste(length(ids), shown_label, "IDs copied!"), type = "message")
  })
  
  # ============================================================================
  # 8.7 ENRICHMENT ANALYSIS - Gene Ontology, KEGG, Reactome (Proteomics Only)
  # ============================================================================

  output$enrich_available_ontologies <- renderUI({
    hi <- !is.null(enrichment_id_col()); hg <- !is.null(go_table()); hk <- !is.null(kegg_table()); hr <- !is.null(reactome_table()); hu <- !is.null(uniprot_table())
    tags$div(style="margin-top:5px;", tags$span(icon(if (hu) "check-circle" else "times-circle"),"UniProtTable",style=paste("color:",if (hu) "#27ae60" else "#c0392b",";margin-right:15px;")),
             tags$span(icon(if (hi) "check-circle" else "times-circle"),paste("ID:",if (hi) enrichment_id_col() else "NOT SET"),style=paste("color:",if (hi) "#27ae60" else "#c0392b",";margin-right:15px;")),
             tags$span(icon(if (hg) "check-circle" else "times-circle"),"GO",style=paste("color:",if (hg) "#27ae60" else "#c0392b",";margin-right:15px;")),
             tags$span(icon(if (hk) "check-circle" else "times-circle"),"KEGG",style=paste("color:",if (hk) "#27ae60" else "#c0392b",";margin-right:15px;")),
             tags$span(icon(if (hr) "check-circle" else "times-circle"),"Reactome",style=paste("color:",if (hr) "#27ae60" else "#c0392b")))
  })
  output$enrich_ontology_checkboxes <- renderUI({ ch <- c(); if (!is.null(go_table())) ch <- c(ch,"GO"="GO"); if (!is.null(kegg_table())) ch <- c(ch,"KEGG"="KEGG"); if (!is.null(reactome_table())) ch <- c(ch,"REACTOME"="REACTOME")
  if (!length(ch)) return(HTML("<p style='color:#c0392b;'>No tables loaded.</p>")); checkboxGroupInput("enrich_ontology_types",NULL,choices=ch,selected=ch,inline=TRUE) })
  output$enrich_ks_ontology_checkboxes <- renderUI({ ch <- c(); if (!is.null(go_table())) ch <- c(ch,"GO"="GO"); if (!is.null(kegg_table())) ch <- c(ch,"KEGG"="KEGG"); if (!is.null(reactome_table())) ch <- c(ch,"REACTOME"="REACTOME")
  if (!length(ch)) return(HTML("<p style='color:#c0392b;'>No tables loaded.</p>")); checkboxGroupInput("enrich_ks_ontology_types",NULL,choices=ch,selected=ch,inline=TRUE) })
  observe({ req(romics_object()); sc <- get_enrichment_stat_cols(romics_object()); ch <- c("None"="none",setNames(sc,sc))
  for (s in c("enrich_stat_col_1","enrich_stat_col_2","enrich_stat_col_3","enrich_stat_col_4","enrich_univ_stat_col")) updateSelectInput(session,s,choices=ch)
  updateSelectInput(session,"enrich_ks_ranking_col",choices=sc) })
  observe({ req(romics_object()); obj <- romics_object(); ch <- build_feature_choices(obj,include_all=FALSE,include_anova=TRUE,include_clusters=TRUE,heatmap_clusters_data=heatmap_clusters())
  updateSelectizeInput(session,"enrich_query_dropdown",choices=ch,server=TRUE)
  fi <- rownames(obj$data); updateSelectizeInput(session,"enrich_universe_dropdown",choices=c("All"="ALL_DETECTED",setNames(as.list(fi),fi)),server=TRUE) })
  observeEvent(input$enrich_query_add_btn, { req(input$enrich_query_dropdown); sel <- input$enrich_query_dropdown; ct <- input$enrich_query_text; obj <- romics_object()
  resolved <- resolve_dropdown_selection(sel,obj,heatmap_clusters())
  if (length(resolved)>1||resolved[1]!=sel) { nt <- paste(resolved,collapse=", "); if (!is.null(ct)&&ct!="") nt <- paste(ct,nt,sep=", "); updateTextInput(session,"enrich_query_text",value=nt); showNotification(paste(length(resolved),"proteins added"),type="message") }
  else { if (is.null(ct)||ct=="") updateTextInput(session,"enrich_query_text",value=sel) else updateTextInput(session,"enrich_query_text",value=paste(ct,sel,sep=", ")) }
  updateSelectizeInput(session,"enrich_query_dropdown",selected="") })
  observeEvent(input$enrich_universe_add_btn, { req(input$enrich_universe_dropdown); sel <- input$enrich_universe_dropdown
  if (sel=="ALL_DETECTED") updateTextInput(session,"enrich_universe_text",value=paste(rownames(romics_object()$data),collapse=", "))
  else { ct <- input$enrich_universe_text; if (is.null(ct)||ct=="") updateTextInput(session,"enrich_universe_text",value=sel) else updateTextInput(session,"enrich_universe_text",value=paste(ct,sel,sep=", ")) }
  updateSelectizeInput(session,"enrich_universe_dropdown",selected="") })
  observeEvent(input$enrich_apply_query_filter, { req(romics_object()); obj <- romics_object(); if (is.null(obj$statistics)) { showNotification("No statistics!",type="error"); return() }
  sdf <- obj$statistics; keep <- rep(TRUE,nrow(sdf))
  for (i in 1:4) { cn <- input[[paste0("enrich_stat_col_",i)]]; fe <- input[[paste0("enrich_stat_filter_",i)]]
  if (!is.null(cn)&&cn!="none"&&cn %in% colnames(sdf)&&!is.null(fe)&&fe!="") keep <- keep & apply_stat_filter(as.numeric(sdf[,cn]),fe) }
  ff <- rownames(sdf)[keep]; updateTextInput(session,"enrich_query_text",value=paste(ff,collapse=", ")); showNotification(paste("Filtered:",length(ff)),type="message") })
  observeEvent(input$enrich_apply_univ_filter, { req(romics_object()); obj <- romics_object(); cn <- input$enrich_univ_stat_col; fe <- input$enrich_univ_stat_filter
  if (is.null(cn)||cn=="none") { updateTextInput(session,"enrich_universe_text",value=paste(rownames(obj$data),collapse=", ")); return() }
  sdf <- obj$statistics; vals <- as.numeric(sdf[,cn]); keep <- apply_stat_filter(vals,fe); ff <- rownames(sdf)[keep]
  updateTextInput(session,"enrich_universe_text",value=paste(ff,collapse=", ")); showNotification(paste("Universe:",length(ff)),type="message") })
  output$enrich_query_summary <- renderPrint({ qt <- input$enrich_query_text; if (is.null(qt)||qt=="") { cat("Empty.\n"); return() }; ql <- trimws(strsplit(qt,",")[[1]]); cat("Query:",length(ql[ql!=""]),"proteins\n") })
  output$enrich_universe_summary <- renderPrint({ if (input$enrich_universe_mode=="all") { if (!is.null(romics_object())) cat("All",nrow(romics_object()$data),"proteins\n"); return() }
    ut <- input$enrich_universe_text; if (is.null(ut)||ut=="") { cat("Empty.\n"); return() }; ul <- trimws(strsplit(ut,",")[[1]]); cat("Universe:",length(ul[ul!=""]),"\n") })
  output$enrich_ks_preview <- renderPrint({ req(romics_object(),input$enrich_ks_ranking_col); obj <- romics_object()
  if (is.null(obj$statistics)||!input$enrich_ks_ranking_col %in% colnames(obj$statistics)) { cat("Not available.\n"); return() }
  vals <- as.numeric(obj$statistics[,input$enrich_ks_ranking_col]); cat("Column:",input$enrich_ks_ranking_col,"\nTotal:",length(vals),"\nRange:",round(min(vals,na.rm=TRUE),4),"to",round(max(vals,na.rm=TRUE),4),"\n") })
  
  # RUN LIST ENRICHMENT
  observeEvent(input$run_enrichment, {
    req(romics_object(),enrichment_id_col()); obj <- romics_object(); id_col <- enrichment_id_col(); tt <- input$enrich_test_type
    ot <- input$enrich_ontology_types; if (is.null(ot)||!length(ot)) { showNotification("Select ontology!",type="error"); return() }
    qt <- input$enrich_query_text; if (is.null(qt)||qt=="") { showNotification("Query empty!",type="error"); return() }
    qf <- trimws(strsplit(qt,",")[[1]]); qf <- qf[qf!=""]; qi <- get_enrichment_ids_for_features(obj,qf,id_col); if (!length(qi)) { showNotification(paste("No valid IDs found!\n\nproteinMinion uses either UniProt accession or UniProt names.\nWe recommend changing the ID type of the romics object using romicsChangeIDs()."),type="error",duration=NULL); return() }
    if (input$enrich_universe_mode=="all") uf <- rownames(obj$data) else { ut <- input$enrich_universe_text; uf <- if (is.null(ut)||ut=="") rownames(obj$data) else { x <- trimws(strsplit(ut,",")[[1]]); x[x!=""] } }
    ui <- get_enrichment_ids_for_features(obj,uf,id_col); if (!length(ui)) { showNotification("No universe IDs!",type="error"); return() }
    qi <- qi[qi %in% ui]; if (!length(qi)) { showNotification("Query not in universe!",type="error"); return() }
    withProgress(message="Running enrichment...", value=0, { tryCatch({
      if (!is.null(uniprot_table())) assign("UniProtTable",uniprot_table(),envir=.GlobalEnv) else if (!exists("UniProtTable",envir=.GlobalEnv)) { showNotification("UniProtTable not loaded!",type="error"); return() }
      if ("GO" %in% ot&&!is.null(go_table())) assign("UniProtTable_GO",go_table(),envir=.GlobalEnv)
      if ("KEGG" %in% ot&&!is.null(kegg_table())) assign("UniProtTable_KEGG",kegg_table(),envir=.GlobalEnv)
      if ("REACTOME" %in% ot&&!is.null(reactome_table())) assign("UniProtTable_REACTOME",reactome_table(),envir=.GlobalEnv)
      ef <- switch(tt,"Fisher"=ProteinMiniOn::proteinEnrichFisher,"EASE"=ProteinMiniOn::proteinEnrichEASE,"Binomial"=ProteinMiniOn::proteinEnrichBinomial)
      all_res <- list(); for (i in seq_along(ot)) { ont <- ot[i]; incProgress(i/(length(ot)+1),detail=paste("Running",ont,"..."))
      tryCatch({ rs <- ef(query=qi,universe=ui,type=ont,preloaded_UniProtTable=TRUE); if (!is.null(rs)&&is.data.frame(rs)&&nrow(rs)>0) all_res[[ont]] <- rs },
               error=function(e) showNotification(paste(ont,"error:",e$message),type="warning",duration=5)) }
      incProgress(1)
      if (!length(all_res)) { showNotification("No results.",type="warning"); enrich_results(data.frame(Message="No enrichments.")) }
      else { combined <- do.call(rbind,all_res); rownames(combined) <- NULL; if ("Test_p" %in% colnames(combined)) combined <- combined[order(combined$Test_p),]
      enrich_results(combined); enrich_filtered_data(combined); showNotification(paste("Done:",nrow(combined),"results"),type="message") }
    }, error=function(e) { showNotification(paste("Error:",e$message),type="error"); enrich_results(NULL) }) })
  })
  
  # RUN KS ENRICHMENT
  observeEvent(input$run_ks_enrichment, {
    req(romics_object(),enrichment_id_col(),input$enrich_ks_ranking_col); obj <- romics_object(); id_col <- enrichment_id_col()
    rc <- input$enrich_ks_ranking_col; ord <- input$enrich_ks_order; ot <- input$enrich_ks_ontology_types
    if (is.null(ot)||!length(ot)) { showNotification("Select ontology!",type="error"); return() }
    withProgress(message="Running KS...", value=0, { tryCatch({
      incProgress(0.1); if (!is.null(uniprot_table())) assign("UniProtTable",uniprot_table(),envir=.GlobalEnv) else if (!exists("UniProtTable",envir=.GlobalEnv)) { showNotification("UniProtTable not loaded!",type="error"); return() }
      features <- rownames(obj$statistics); vals <- as.numeric(obj$statistics[,rc])
      idm <- data.frame(Feature=rownames(obj$data),ID=obj$IDs[[id_col]],stringsAsFactors=FALSE); idm <- idm[!is.na(idm$ID)&idm$ID!="",]
      sdf <- data.frame(Feature=features,Value=vals,stringsAsFactors=FALSE); merged <- merge(idm,sdf,by="Feature"); merged <- merged[!is.na(merged$Value),]; merged <- merged[!duplicated(merged$ID),]
      rt <- data.frame(proteins=merged$ID,scores=merged$Value,stringsAsFactors=FALSE); if (nrow(rt)<10) { showNotification("Too few proteins!",type="error"); return() }; incProgress(0.2)
      if ("GO" %in% ot&&!is.null(go_table())) assign("UniProtTable_GO",go_table(),envir=.GlobalEnv)
      if ("KEGG" %in% ot&&!is.null(kegg_table())) assign("UniProtTable_KEGG",kegg_table(),envir=.GlobalEnv)
      if ("REACTOME" %in% ot&&!is.null(reactome_table())) assign("UniProtTable_REACTOME",reactome_table(),envir=.GlobalEnv)
      all_res <- list(); for (i in seq_along(ot)) { ont <- ot[i]; incProgress(0.2+(0.6*i/length(ot)),detail=paste("Running KS for",ont,"..."))
      tryCatch({ rs <- ProteinMiniOn::proteinEnrichKS(rankingTable=rt,type=ont,order=ord,preloaded_UniProtTable=TRUE); if (!is.null(rs)&&is.data.frame(rs)&&nrow(rs)>0) all_res[[ont]] <- rs },
               error=function(e) showNotification(paste(ont,"KS error:",e$message),type="warning",duration=5)) }
      if (!length(all_res)) { showNotification("No KS results.",type="warning"); enrich_results(data.frame(Message="No enrichments.")) }
      else { combined <- do.call(rbind,all_res); rownames(combined) <- NULL; if ("Test_p" %in% colnames(combined)) combined <- combined[order(combined$Test_p),]
      enrich_results(combined); enrich_filtered_data(combined); showNotification(paste("KS done:",nrow(combined)),type="message") }
    }, error=function(e) { showNotification(paste("Error:",e$message),type="error"); enrich_results(NULL) }) })
  })
  
  # ENRICHMENT RESULTS FILTERS
  observe({ res <- enrich_results(); if (is.null(res)||!"Category" %in% colnames(res)) { updateSelectInput(session,"enrich_filter_category",choices=c("All"="all")); return() }
  cats <- sort(unique(res$Category)); updateSelectInput(session,"enrich_filter_category",choices=c("All Categories"="all",setNames(cats,cats)),selected="all") })
  observe({ res <- enrich_results(); if (!is.null(res)&&nrow(res)>0&&!"Message" %in% colnames(res)) enrich_filtered_data(res) })
  
  observeEvent(input$apply_enrich_filters, {
    res <- enrich_results(); if (is.null(res)||nrow(res)==0||"Message" %in% colnames(res)) { showNotification("No results to filter.",type="warning"); return() }
    filtered <- res
    cat_sel <- input$enrich_filter_category; if (!is.null(cat_sel)&&cat_sel!="all"&&"Category" %in% colnames(filtered)) filtered <- filtered[filtered$Category==cat_sel,,drop=FALSE]
    ptype <- input$enrich_filter_ptype; pmax <- input$enrich_filter_pmax
    if (!is.null(ptype)&&!is.null(pmax)) { pcol <- if (ptype=="padj"&&"Test_padj" %in% colnames(filtered)) "Test_padj" else if ("Test_p" %in% colnames(filtered)) "Test_p" else NULL
    if (!is.null(pcol)) { pvals <- as.numeric(filtered[[pcol]]); filtered <- filtered[!is.na(pvals)&pvals<=pmax,,drop=FALSE] } }
    fc_min <- input$enrich_filter_fc; if (!is.null(fc_min)&&"fold_change" %in% colnames(filtered)) { fc_vals <- as.numeric(filtered$fold_change); filtered <- filtered[!is.na(fc_vals)&fc_vals>=fc_min,,drop=FALSE] }
    min_feat <- input$enrich_filter_minfeature; if (!is.null(min_feat)&&"Count_query" %in% colnames(filtered)) { counts <- as.numeric(filtered$Count_query); filtered <- filtered[!is.na(counts)&counts>=min_feat,,drop=FALSE] }
    enrich_filtered_data(filtered); showNotification(paste("Filtered:",nrow(filtered),"terms"),type="message")
  })
  
  output$enrich_results_table <- renderDT({
    fd <- enrich_filtered_data(); if (is.null(fd)) { res <- enrich_results(); if (is.null(res)) return(NULL); fd <- res }
    if (nrow(fd)==0) return(datatable(data.frame(Message="No results match filters.")))
    datatable(fd,extensions=c("Buttons","Scroller"),options=list(dom="Bfrtip",buttons=c("copy","csv","excel"),pageLength=25,scrollX=TRUE,scrollY="500px",scroller=TRUE),rownames=FALSE,filter="top")
  })
  output$enrich_results_summary <- renderPrint({
    fd <- enrich_filtered_data(); if (is.null(fd)||nrow(fd)==0) { cat("No results.\n"); return() }; cat("Showing:",nrow(fd),"terms\n")
    if ("Category" %in% colnames(fd)) { tb <- table(fd$Category); for (nm in sort(names(tb))) {
      sig <- if ("Test_padj" %in% colnames(fd)) sum(fd$Test_padj[fd$Category==nm]<0.05,na.rm=TRUE) else if ("Test_p" %in% colnames(fd)) sum(fd$Test_p[fd$Category==nm]<0.05,na.rm=TRUE) else NA
      cat(" ",nm,":",tb[nm],"terms"); if (!is.na(sig)) cat(" (",sig,"sig.)"); cat("\n") } }
  })
  output$download_enrichment_csv <- downloadHandler(filename=function() paste0("enrichment_",Sys.Date(),".csv"), content=function(file) { fd <- enrich_filtered_data(); if (is.null(fd)) fd <- enrich_results(); if (!is.null(fd)) write.csv(fd,file,row.names=FALSE) })
  observeEvent(input$copy_enrichment_ids, { fd <- enrich_filtered_data(); if (is.null(fd)) fd <- enrich_results(); if (is.null(fd)||nrow(fd)==0) return()
  if ("Term_description" %in% colnames(fd)) ids <- paste(unique(fd$Term_description),collapse=", ") else if ("Term_accession" %in% colnames(fd)) ids <- paste(unique(fd$Term_accession),collapse=", ") else ids <- ""
  session$sendCustomMessage("copyClipboard",list(text=ids)); showNotification("Term IDs copied!",type="message") })
  
  # ENRICHMENT PLOT
  observeEvent(input$generate_enrich_plot, {
    fd <- enrich_filtered_data(); if (is.null(fd)||nrow(fd)==0) { showNotification("No filtered results to plot.",type="warning"); return() }
    withProgress(message="Generating plot...", value=0, { tryCatch({
      incProgress(0.3)
      if (has_proteinminion && exists("proteinEnrichmentPlot",where="package:ProteinMiniOn",mode="function")) {
        p <- ProteinMiniOn::proteinEnrichmentPlot(fd, Test_p_type=input$enrich_filter_ptype, plotly=FALSE)
        enrich_plot_store(p); output$enrich_bar_plot <- renderPlotly({ plotly::ggplotly(p) })
      } else {
        fd$neg_log10_p <- -log10(pmax(fd$Test_p,1e-300)); if (!"Category" %in% colnames(fd)) fd$Category <- "Unknown"
        fd$label <- if ("Term_description" %in% colnames(fd)) paste0(fd$Term_accession," : ",fd$Term_description) else fd$Term_accession
        fd <- fd[order(fd$Category,fd$neg_log10_p),]; fd$label <- factor(fd$label,levels=unique(fd$label))
        p <- ggplot(fd,aes(x=neg_log10_p,y=label,fill=Category,text=paste("Term:",label,"<br>p:",format(Test_p,scientific=TRUE,digits=3))))+geom_bar(stat="identity")+scale_fill_brewer(palette="Set1")+labs(x="-Log10(p)",y="",title="Enrichment")+theme_minimal()+theme(axis.text.y=element_text(size=8))
        enrich_plot_store(p); output$enrich_bar_plot <- renderPlotly({ plotly::ggplotly(p,tooltip="text")%>%layout(height=max(400,nrow(fd)*20+100)) })
      }
      incProgress(1); showNotification("Plot generated!",type="message")
    }, error=function(e) { showNotification(paste("Plot error:",e$message),type="error"); output$enrich_bar_plot <- renderPlotly({ plotly_empty() }) }) })
  })
  output$download_enrich_plot <- downloadHandler(filename=function() paste0("enrichment_plot_",Sys.Date(),".png"),
                                                 content=function(file) { p <- enrich_plot_store(); if (!is.null(p)) { fd <- enrich_filtered_data(); h <- if (!is.null(fd)) max(6,0.3*nrow(fd)) else 8; ggsave(file,plot=p,width=12,height=h,dpi=300,bg="white") } })
  
  # ---- HISTORY ----
  output$uuid_box <- renderUI({
    uuid_val <- if (is.null(romics_object())) "None" else if (!is.null(romics_object()$uuid)) romics_object()$uuid else "N/A"
    fluidRow(box(title="UUID",status="info",solidHeader=TRUE,width=12,
      p(uuid_val,style="font-family:monospace;font-size:13px;word-break:break-all;")
    ))
  })
  output$history_object_summary <- renderPrint({ req(romics_object()); obj <- romics_object(); cat("Dim:",nrow(obj$data),"x",ncol(obj$data),"\n"); if (!is.null(current_factor())) cat("Factor:",current_factor(),"\n")
  for (l in c("data","metadata","IDs","statistics","steps","embeddings")) cat(" -",paste0(l,":"),ifelse(!is.null(obj[[l]]),"[YES]","[NO]"),"\n") })
  output$embedding_info <- renderPrint({ req(romics_object()); emb <- check_embeddings_available(romics_object()); st <- function(a,t) if (a) paste0("[YES] ",get_embedding_dimensions(romics_object(),t)," comp") else "[NO]"
  cat("PCA:",st(emb$pca,"pca"),"\nUMAP:",st(emb$umap,"umap"),"\nt-SNE:",st(emb$tsne,"tsne"),"\n") })
    # steps_output_html removed - only showing detailed view
  output$steps_output_detailed <- renderPrint({
    req(romics_object())
    obj <- romics_object()
    show_ts <- if (is.null(input$show_timestamps)) TRUE else input$show_timestamps
    
    if (is.null(obj$steps) || !length(obj$steps)) { cat("No steps.\n"); return() }
    steps <- obj$steps
    
    if (!is.character(steps)) return()
    
    # Filter and process steps
    step_num <- 1
    i <- 1
    while (i <= length(steps)) {
      line <- steps[i]
      
      # Skip romics_object, NULL, and empty lines
      if (line == "romics_object" || line == "NULL" || line == "") {
        i <- i + 1
        next
      }
      
      # If current line is date|, check if next is fun|
      if (grepl("^date\\|", line)) {
        next_idx <- i + 1
        # If next line is fun|, print fun first then date
        if (next_idx <= length(steps) && grepl("^fun\\|", steps[next_idx])) {
          cat(sprintf("Step %d: %s\n", step_num, steps[next_idx]))
          if (show_ts) {
            cat(sprintf("  %s\n", line))
          }
          step_num <- step_num + 1
          i <- i + 2  # Skip both date and fun lines
          next
        } else if (!show_ts) {
          # If timestamps off, skip standalone date lines
          i <- i + 1
          next
        }
      }
      
      # If it's a fun| line (not preceded by date)
      if (grepl("^fun\\|", line)) {
        cat(sprintf("Step %d: %s\n", step_num, line))
        step_num <- step_num + 1
      } else {
        # Other lines
        cat(sprintf("%s\n", line))
      }
      
      i <- i + 1
    }
  })
  output$dependencies_table <- renderDT({ req(romics_object()); if (is.null(romics_object()$dependencies)) datatable(data.frame(Message="No dependencies."),options=list(dom="t"))
    else { d <- romics_object()$dependencies; if (is.data.frame(d)) datatable(d,options=list(pageLength=20,scrollX=TRUE,dom="Bfrtip",buttons=c("copy","csv")),extensions="Buttons")
    else datatable(data.frame(D=as.character(d)),options=list(dom="t")) } })

  # ============================================================================
  # 8.8 LIPID INFORMATION TAB - LipidMiner & Lipid Ontology (Lipidomics Only)
  # ============================================================================

  lipidminer_results <- reactiveVal(NULL)
  lipid_ontology_results <- reactiveVal(NULL)
  lipid_query_list <- reactiveVal(character(0))
  lipid_universe_list <- reactiveVal(character(0))
  lipid_query_plot_store <- reactiveVal(NULL)
  lipid_universe_plot_store <- reactiveVal(NULL)
  lipid_enrichment_results_store <- reactiveVal(NULL)

  observe({
    req(romics_object(), is_lipidominion())
    obj <- romics_object()
    tryCatch({
      # Run lipidMiner on all lipid row names
      lipids <- rownames(obj$data)
      lm_results <- LipidMiniOn::lipidMiner(lipids)
      lipidminer_results(lm_results)

      # Run Lipid_ontology_maker on recognized lipids
      recognized_lipids <- lm_results$Lipid_names[!is.na(lm_results$Lipid_names)]
      if (length(recognized_lipids) > 0) {
        lo_results <- LipidMiniOn::Lipid_ontology_maker(recognized_lipids)
        lipid_ontology_results(lo_results)
      }

      # Update dropdown choices with recognized lipids
      updateSelectizeInput(session, "lipid_query_dropdown", choices = recognized_lipids, server = TRUE)
      updateSelectizeInput(session, "lipid_universe_dropdown", choices = recognized_lipids, server = TRUE)
    }, error = function(e) {
      showNotification(paste("Error processing lipids:", e$message), type = "error")
    })
  })

  output$lipidminer_summary <- renderPrint({
    lm <- lipidminer_results()
    if (is.null(lm)) {
      cat("LipidMiner results not yet available.\n")
      return()
    }
    total_lipids <- nrow(lm)
    recognized <- sum(!is.na(lm$Lipid_names))
    unrecognized <- total_lipids - recognized
    pct_unrecognized <- round(unrecognized / total_lipids * 100, 1)

    cat("Total lipids:", total_lipids, "\n")
    cat("Recognized:", recognized, "\n")
    cat("Unrecognized:", unrecognized, sprintf(" (%.1f%%)\n", pct_unrecognized))
  })

  output$lipidminer_warning <- renderUI({
    lm <- lipidminer_results()
    if (is.null(lm)) return(NULL)
    unrecognized <- sum(is.na(lm$Lipid_names))
    pct_unrecognized <- unrecognized / nrow(lm) * 100
    if (pct_unrecognized > 20) {
      tags$div(style="background-color:#fff3cd;border:1px solid #ffc107;border-radius:4px;padding:10px;margin-bottom:10px;",
        tags$strong("⚠️ Warning:"), sprintf(" More than 20%% of lipids are not recognized by LipidMiner (%.1f%%).", pct_unrecognized),
        tags$br(),
        "Please ensure lipids are in shorthand nomenclature (e.g., 'PE(36:2)', 'TG(52:3)').",
        tags$br(),
        "For nomenclature conversion, use LipidLynxX: ",
        tags$a(href="https://www.lipidlynxX.com", target="_blank", "https://www.lipidlynxX.com")
      )
    }
  })

  output$lipidminer_table <- renderDT({
    lm <- lipidminer_results()
    if (is.null(lm)) return(datatable(data.frame(Message="No data available.")))
    datatable(lm, options = list(pageLength = 25, scrollX = TRUE, dom = "Bfrtip", buttons = c("copy", "csv")), extensions = "Buttons")
  })

  output$lipid_ontology_table <- renderDT({
    lo <- lipid_ontology_results()
    if (is.null(lo)) return(datatable(data.frame(Message="No lipid ontology data available.")))
    datatable(lo, options = list(pageLength = 25, scrollX = TRUE, dom = "Bfrtip", buttons = c("copy", "csv")), extensions = "Buttons")
  })

  # ============================================================================
  # 8.9 LIPID ENRICHMENT TAB - Query vs Universe Analysis
  # ============================================================================

  observeEvent(input$lipid_query_add_btn, {
    req(input$lipid_query_dropdown)
    current <- lipid_query_list()
    new_lipid <- input$lipid_query_dropdown
    if (!new_lipid %in% current) {
      lipid_query_list(c(current, new_lipid))
    }
  })

  observeEvent(input$lipid_universe_add_btn, {
    req(input$lipid_universe_dropdown)
    current <- lipid_universe_list()
    new_lipid <- input$lipid_universe_dropdown
    if (!new_lipid %in% current) {
      lipid_universe_list(c(current, new_lipid))
    }
  })

  output$lipid_query_summary <- renderPrint({
    query_list <- c(lipid_query_list(), strsplit(input$lipid_query_text, ",\\s*")[[1]])
    query_list <- query_list[query_list != ""]
    if (length(query_list) == 0) {
      cat("No query lipids selected.\n")
    } else {
      cat("Query lipids (", length(query_list), "):\n", sep = "")
      for (i in seq_len(min(10, length(query_list)))) {
        cat(" ", query_list[i], "\n", sep = "")
      }
      if (length(query_list) > 10) cat(" ... and", length(query_list) - 10, "more\n")
    }
  })

  output$lipid_universe_summary <- renderPrint({
    if (input$lipid_universe_mode == "all") {
      lm <- lipidminer_results()
      if (!is.null(lm)) {
        universe_count <- sum(!is.na(lm$Lipid_names))
        cat("Universe: All detected lipids\n")
        cat("Total recognized lipids:", universe_count, "\n")
      }
    } else {
      universe_list <- c(lipid_universe_list(), strsplit(input$lipid_universe_text, ",\\s*")[[1]])
      universe_list <- universe_list[universe_list != ""]
      if (length(universe_list) == 0) {
        cat("No custom universe lipids selected.\n")
      } else {
        cat("Custom universe lipids (", length(universe_list), "):\n", sep = "")
        for (i in seq_len(min(10, length(universe_list)))) {
          cat(" ", universe_list[i], "\n", sep = "")
        }
        if (length(universe_list) > 10) cat(" ... and", length(universe_list) - 10, "more\n")
      }
    }
  })

  observeEvent(input$generate_lipid_query_plot, {
    query_list <- c(lipid_query_list(), strsplit(input$lipid_query_text, ",\\s*")[[1]])
    query_list <- query_list[query_list != ""]
    if (length(query_list) == 0) {
      showNotification("Please add query lipids.", type = "warning")
      return()
    }
    lo <- lipid_ontology_results()
    if (is.null(lo)) {
      showNotification("Lipid ontology data not available.", type = "warning")
      return()
    }
    # Create plot data for query lipids
    query_data <- lo[lo$Lipid_name %in% query_list, ]
    if (nrow(query_data) > 0) {
      # Create summary plot (e.g., lipid classes)
      p <- ggplot(query_data, aes(x = reorder(Lipid_class, Lipid_class, function(x) -length(x)))) +
        geom_bar(fill = "#2ecc71", alpha = 0.7) +
        labs(title = "Query Lipids: Class Distribution", x = "Lipid Class", y = "Count") +
        theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      lipid_query_plot_store(p)
    } else {
      showNotification("No query lipids found in ontology.", type = "warning")
    }
  })

  observeEvent(input$generate_lipid_universe_plot, {
    if (input$lipid_universe_mode == "all") {
      lo <- lipid_ontology_results()
      if (is.null(lo)) {
        showNotification("Lipid ontology data not available.", type = "warning")
        return()
      }
      universe_data <- lo
    } else {
      universe_list <- c(lipid_universe_list(), strsplit(input$lipid_universe_text, ",\\s*")[[1]])
      universe_list <- universe_list[universe_list != ""]
      if (length(universe_list) == 0) {
        showNotification("Please add custom universe lipids.", type = "warning")
        return()
      }
      lo <- lipid_ontology_results()
      if (is.null(lo)) {
        showNotification("Lipid ontology data not available.", type = "warning")
        return()
      }
      universe_data <- lo[lo$Lipid_name %in% universe_list, ]
    }
    if (nrow(universe_data) > 0) {
      # Create summary plot for universe
      p <- ggplot(universe_data, aes(x = reorder(Lipid_class, Lipid_class, function(x) -length(x)))) +
        geom_bar(fill = "#3498db", alpha = 0.7) +
        labs(title = "Universe Lipids: Class Distribution", x = "Lipid Class", y = "Count") +
        theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      lipid_universe_plot_store(p)
    } else {
      showNotification("No universe lipids found in ontology.", type = "warning")
    }
  })

  output$lipid_query_plot <- renderPlotly({
    p <- lipid_query_plot_store()
    if (is.null(p)) return(NULL)
    ggplotly(p, tooltip = "text")
  })

  output$lipid_universe_plot <- renderPlotly({
    p <- lipid_universe_plot_store()
    if (is.null(p)) return(NULL)
    ggplotly(p, tooltip = "text")
  })

  observeEvent(input$run_lipid_enrichment, {
    req(romics_object(), lipidminer_results())
    obj <- romics_object()
    lo <- lipid_ontology_results()

    # Get query lipids
    query_list <- c(lipid_query_list(), trimws(strsplit(input$lipid_query_text, ",")[[1]]))
    query_list <- query_list[query_list != ""]
    if (length(query_list) == 0) {
      showNotification("Please add query lipids.", type = "warning")
      return()
    }

    # Get universe lipids
    if (input$lipid_universe_mode == "all") {
      universe_list <- rownames(obj$data)
    } else {
      universe_list <- c(lipid_universe_list(), trimws(strsplit(input$lipid_universe_text, ",")[[1]]))
      universe_list <- universe_list[universe_list != ""]
      if (length(universe_list) == 0) {
        universe_list <- rownames(obj$data)
      }
    }

    # Validate query is subset of universe
    query_list <- query_list[query_list %in% universe_list]
    if (length(query_list) == 0) {
      showNotification("No query lipids found in universe.", type = "warning")
      return()
    }

    showNotification("Running lipid enrichment analysis...", type = "message")

    withProgress(message = "Running enrichment...", value = 0, {
      tryCatch({
        # Get enrichment method and p-value cutoffs
        method <- input$lipid_enrich_method
        p_cutoff <- input$lipid_enrich_p
        adj_p_cutoff <- input$lipid_enrich_adj_p

        # Call appropriate LipidMiniOn function
        enrich_func <- switch(method,
          "fisher" = LipidMiniOn::lipid_fisher_enrich,
          "ease" = LipidMiniOn::lipid_EASE_enrich,
          "hyper" = LipidMiniOn::lipid_hyper_enrich,
          "binom" = LipidMiniOn::lipid_binom_enrich,
          LipidMiniOn::lipid_fisher_enrich)

        incProgress(0.5)

        # Run enrichment
        results <- enrich_func(query = query_list, universe = universe_list,
                               p = p_cutoff, adj_p = adj_p_cutoff)

        incProgress(1)

        # Store results
        if (!is.null(results) && is.data.frame(results) && nrow(results) > 0) {
          lipid_enrichment_results_store(results)
          showNotification(paste("Enrichment complete:", nrow(results), "terms enriched"),
                          type = "message")
        } else {
          lipid_enrichment_results_store(NULL)
          showNotification("No enriched terms found with selected cutoffs.", type = "warning")
        }
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
        lipid_enrichment_results_store(NULL)
      })
    })
  })

  output$lipid_enrichment_results <- renderDT({
    res <- lipid_enrichment_results_store()
    if (is.null(res)) return(datatable(data.frame(Message = "Run enrichment analysis to see results.")))
    datatable(res, options = list(pageLength = 20, scrollX = TRUE, dom = "Bfrtip", buttons = c("copy", "csv")), extensions = "Buttons")
  })

  output$download_lipid_enrichment <- downloadHandler(
    filename = function() paste0("lipid_enrichment_", Sys.Date(), ".csv"),
    content = function(file) {
      res <- lipid_enrichment_results_store()
      if (!is.null(res)) write.csv(res, file, row.names = FALSE)
    }
  )

  observeEvent(input$copy_lipid_enrichment_ids, {
    res <- lipid_enrichment_results_store()
    if (!is.null(res) && nrow(res) > 0) {
      # Extract ontology terms from results
      if ("Ontology_term" %in% colnames(res)) {
        ids <- unique(res$Ontology_term)
        session$sendCustomMessage("copyClipboard", list(text = paste(ids, collapse = ", ")))
        showNotification(paste(length(ids), "ontology terms copied!"), type = "message")
      } else if ("Term_description" %in% colnames(res)) {
        ids <- unique(res$Term_description)
        session$sendCustomMessage("copyClipboard", list(text = paste(ids, collapse = ", ")))
        showNotification(paste(length(ids), "terms copied!"), type = "message")
      }
    }
  })

}  # ---- END SERVER ----

# ============================================================================
# LAUNCH SHINY APPLICATION
# ============================================================================

shinyApp(ui=ui, server=server)