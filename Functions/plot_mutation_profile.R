plot_mutation_profile = function(clone_context, name){
  total_sequences = c()
  for (mut in clone_context) {
    sequence = strsplit(mut, "")[[1]]
    total_sequences = c(total_sequences, sequence)
  }  
  mat = matrix(total_sequences, ncol = 21, byrow = T)
  print(mat)
  frequency_table = data.frame()
  for (col in 1:ncol(mat)){
    col = c(mat[,col], c("A", "C", "G", "T")) # add one nucleotide of each position to prevent zero's 
    frequency = table(col)
    frequency = frequency -1 # remove the one nucelotide to return to original values
    frequency_table = rbind(frequency_table, frequency)
    frequency_table = frequency_table/rowSums(frequency_table)
  }
  
  colnames(frequency_table) <- c("A", "C", "G", "T")
  frequency_table$position = c(-10:10)
  freq_table = melt(frequency_table, id.vars = "position")
  positional_information_plot = ggplot(freq_table, aes(x = position, y = value, color = variable)) + 
    geom_point(size = 1.5) + geom_line(size = 0.5) + scale_color_manual(values = c("#006400", "Blue", "Orange", "Red")) + theme_BM() +
    theme(legend.title = element_blank()) + 
    #ggtitle(name) +
    ylab(label = "Relative fraction of mutated bases") + scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05))
  ggsave(paste0(name, "_context.pdf"), positional_information_plot, width = 4, height = 3, useDingbats = F)
}