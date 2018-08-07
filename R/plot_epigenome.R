#' A Function That Plots Roadmap Data
#'
#' This function plots roadmap data.
#' @param roadmap_matrices_metadata metadata of roadmap matrices
#' @param matrices matrices to be plotted
#' @param marks marks of interest
#' @param list_of_feat list of index of gene (same index as roadmap matrices)
#' @param feat_types         A indexed vector of types 

#' @param PLOT_TRSCR_TYPE    boolean specifyinf what to plot
#' @param PLOT_NORM_EXPR     boolean specifyinf what to plot
#' @param PLOT_SPERMATO_EXPR boolean specifyinf what to plot
#' @param PLOT_CPG           boolean specifyinf what to plot
#' @param PLOT_METH          boolean specifyinf what to plot
#' @param cex                grpahical parameter specifying size of characters

#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices rainbow
#' @importFrom graphics barplot 
#' @importFrom graphics layout 
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics image 
#' @importFrom graphics text
#' @importFrom stats na.omit
#' @importFrom stats quantile
#' @export
plot_epigenome = function(
  roadmap_matrices_metadata, 
  matrices,
  marks, 
  list_of_feat, 
  feat_types,
  PLOT_TRSCR_TYPE=FALSE, 
  PLOT_NORM_EXPR=FALSE, 
  PLOT_SPERMATO_EXPR=FALSE, 
  PLOT_CPG=FALSE, 
  PLOT_METH=FALSE, 
  cex=2.5
) {

  

  sample_labels = unlist(sapply(roadmap_matrices_metadata[marks], "[[", "sample_labels"))
  nb_samples = length(sample_labels)

  nb_lines = length(list_of_feat) + 1
  nb_columns = 1 + PLOT_TRSCR_TYPE + PLOT_NORM_EXPR + PLOT_SPERMATO_EXPR + nb_samples + PLOT_CPG + PLOT_METH * nrow(matrices["betas"][[1]])

  layout(matrix(1:(nb_columns * nb_lines), nb_lines, byrow=TRUE), respect=FALSE)
  # head line
  par(mar=c(.1,.1,.1,.1))
  plot(0,0,col=0, xaxt="n", yaxt="n")
  text(0, 0, "genes (nb)", cex=cex, srt=90)  
  if (PLOT_TRSCR_TYPE) {
    cols = rainbow(length(na.omit(unique(feat_types))))
    par(mar=c(9, 2, 0, 2))
    barplot(table(unique(feat_types)), las=2, col=cols, yaxt="n") 
  }
  par(mar=c(.1,.1,.1,.1))
  if (PLOT_NORM_EXPR) {
    plot(0,0,col=0, xaxt="n", yaxt="n")
    text(0, 0.5, "normal", cex=cex)
    text(0, -0.5, "expr.", cex=cex)
  }
  if (PLOT_SPERMATO_EXPR) {
    plot(0,0,col=0, xaxt="n", yaxt="n")
    text(0, 0.5, "spermato", cex=cex)
    text(0, -0.5, "expr.", cex=cex)
  }
  sapply(sample_labels, function(t) {
    plot(0,0,col=0, xaxt="n", yaxt="n")
    text(0, 0, t, cex=cex, srt=90)
    # t = strsplit(t, "_")[[1]]
    # text(0, 0.5, t[1], cex=cex)
    # text(0, -0.5, t[2], cex=cex)
  })
  if (PLOT_CPG) {
    plot(0,0,col=0, xaxt="n", yaxt="n")
    text(0, 0, "CpG density", cex=cex, srt=90)
    # text(0, 0.5, "CpG", cex=cex)
    # text(0, -0.5, "density", cex=cex)
  }
  if (PLOT_METH) {
    sapply(rownames(matrices["betas"][[1]]), function(t) {
      plot(0,0,col=0, xaxt="n", yaxt="n")
      text(0, 0.5, t, cex=cex)
      text(0, -0.5, "meth", cex=cex)
    })
  }
  # other lines
  for (j in 1:length(list_of_feat)) {

    tag = names(list_of_feat)[j]
    print(tag)

    # loading gene list
    idx = list_of_feat[[tag]]
  
    # 1st col
    if (length(strsplit(tag, "_")[[1]]) > 3) {
      main = strsplit(tag, "_")[[1]][4]
    } else {
      main = strsplit(tag, "_")[[1]][1]
    }
    main = paste0(main, " (", length(idx), ")")
    plot(0, 0, col=0, xaxt="n", yaxt="n")
    text(0, 0, main, cex=cex, srt=90)

    # transcript types
    if (PLOT_TRSCR_TYPE) {
      cols = rainbow(length(na.omit(unique(feat_types))))
      par(mar=c(0, 2, 0, 2))
      barplot(table(feat_types[idx]), las=2, col=cols, xaxt="n") 
      par(new=TRUE)
      par(mar=c(0, 8, 0, 0))

      data = t(t(as.numeric(feat_types[idx])))
      rownames(data) = idx
      breaks = 1:(length(levels(feat_types[idx])) + 1) - 0.5
      image(t(data[idx,]), col=cols, breaks=breaks, xaxt="n", yaxt="n")

    }
  
    # expression in normal tissues
    if (PLOT_NORM_EXPR) {
      par(mar=c(.1,.1,.1,.1))
      data = matrices["aht_expr"][idx,]
      nb_bin = ncol(data)
      colors=c("green",  "black", "red")
      cols = colorRampPalette(colors)(20)
      breaks = quantile(matrices["aht_expr"],  seq(0,1,length.out=length(cols)+1), na.rm=TRUE)
      image(t(data[idx,]), col=cols, breaks=breaks, xaxt="n", yaxt="n")
    }

    # expression in normal tissues
    if (PLOT_SPERMATO_EXPR) {
      par(mar=c(.1,.1,.1,.1))
      data = matrices["spemato_expr"][idx,]
      nb_bin = ncol(data)
      colors=c("green",  "black", "red")
      cols = colorRampPalette(colors)(20)
      breaks = quantile(matrices["spemato_expr"],  seq(0,1,length.out=length(cols)+1), na.rm=TRUE)
      image(t(data[idx,]), col=cols, breaks=breaks, xaxt="n", yaxt="n")
    }

    for (mark in marks) {
      par(mar=c(.1,.1,.1,.1))
      # roadmap marks
      data = matrices["roadmap_matrices"][[mark]][idx,-(1:6)]
      data = as.matrix(data)
      data = log2(data + 1)
      # plot(density(data))

      nb_samples = length(unlist(roadmap_matrices_metadata[[mark]]$sample_labels))
      nb_bin = ncol(data) / nb_samples
      idx_b = ((1:nb_samples) -1 ) * nb_bin + 1
      idx_e = idx_b + nb_bin - 1
      bins = cbind(idx_b, idx_e)
  
      colors=c("blue", "white", "red")
      cols = colorRampPalette(colors)(20)
      breaks = quantile(data,  seq(0,1,length.out=length(cols)+1), na.rm=TRUE)
      for (i in 1:nrow(bins)) {
        image(t(data[idx,bins[i,1]:bins[i,2]]), col=cols, breaks=breaks, xaxt="n", yaxt="n")
      }  
    }
  
    if (PLOT_CPG) {
      par(mar=c(.1,.1,.1,.1))
      # TSS CpG enrichment
      data = matrices["cpg_density_glo"][idx,-(1:6)]
      data = as.matrix(data)
      # plot(density(data))
    
      nb_bin = ncol(data)
      idx_b = ((1:nb_samples) -1 ) * nb_bin + 1
      idx_e = idx_b + nb_bin - 1
      bins = cbind(idx_b, idx_e)
    
      colors=c("black", "red")
      cols = colorRampPalette(colors)(20)
      all_tss_cpg = as.matrix(matrices["cpg_density_glo"][,-(1:6)])
      breaks = quantile(all_tss_cpg,  seq(0,1,length.out=length(cols)+1), na.rm=TRUE)
      for (i in 1:nrow(bins)) {
        image(t(data[idx,bins[i,1]:bins[i,2]]), col=cols, breaks=breaks, xaxt="n", yaxt="n")
      }
    }

    if (PLOT_METH) {
      par(mar=c(.1,.1,.1,.1))
      # TSS methylation status
      colors=c("green", "black", "red")
      cols = colorRampPalette(colors)(20)

      for (i in 1:length(matrices["betas"])) {
        data = matrices["betas"][[i]]
        breaks = quantile(data,  seq(0,1,length.out=length(cols)+1), na.rm=TRUE)
        image(t(data[idx,]), col=cols, breaks=breaks, xaxt="n", yaxt="n")
      }
    }
  }
}
