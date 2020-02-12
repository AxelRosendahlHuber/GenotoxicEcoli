plot_96_profile3 = function (mut_matrix, colors, ymax = 0.2, condensed = TRUE) {
  number_sigs = dim(mut_matrix)[2]
  
  if (missing(colors)) {
    colors = COLORS6
  }
  if (length(colors) != 6) {
    stop("Provide colors vector with length 6")
  }
  context = CONTEXTS_96
  substitution = rep(SUBSTITUTIONS, each = 16)
  substring(context, 2, 2) = "."
  df = data.frame(substitution = substitution, context = context)
  rownames(mut_matrix) = NULL
  df2 = cbind(df, as.data.frame(mut_matrix))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  value = NULL
  if (condensed) {
    plot = ggplot(data = df3, aes(x = context, y = value, fill = substitution, width = 1)) + geom_bar(stat = "identity", size = 0.2) + scale_fill_manual(values = colors) + 
      facet_grid(variable ~ substitution) + 
      coord_cartesian(ylim = c(0, ymax)) + guides(fill = FALSE) + theme_BM()  + 
      expand_limits(y = 0) +  scale_y_continuous(expand = c(0, 0), limits = c(0,1200), (breaks = seq(0, ymax, 0.1))) +
      theme(axis.title.y = element_blank(), 
            axis.text.y = element_text(size = 8), axis.title.x = element_blank(), 
            axis.text.x = element_text(size = 5, angle = 90, 
                                       vjust = 0.4), strip.text.x = element_text(size = 9), 
            strip.text.y = element_text(size = 9), panel.grid.major.x = element_blank(), 
            panel.spacing.x = unit(0, "lines")) + 
      ylab("Relative contribution") + xlab("96-trinucleotide context")
  }
  else {
    plot = ggplot(data = df3, aes(x = context, y = value, 
                                  fill = substitution, width = 0.8)) + geom_bar(stat = "identity", 
                                                                                size = 0.3) + scale_fill_manual(values = colors) + 
      facet_grid(variable ~ substitution) +
      coord_cartesian(ylim = c(0, ymax)) + scale_y_continuous(breaks = seq(0, 
                                                                           ymax, 0.1)) + guides(fill = FALSE) + theme_BM() + 
      theme(axis.title.y = element_text(size = 12, vjust = 1), 
            axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 12), 
            axis.text.x = element_text(size = 5, angle = 90, 
                                       vjust = 0.3), strip.text.x = element_text(size = 9), 
            strip.text.y = element_text(size = 9), panel.grid.major.x = element_blank(), panel.spacing.x = unit(0.000, "lines")) +  
      ylab("Relative contribution") + xlab("96-trinucleotide context") + scale_y_continuous(expand = c(0, 0))

    
  }
  return(plot)
}