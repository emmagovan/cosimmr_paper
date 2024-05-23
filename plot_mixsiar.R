#Altered plot function to get same outputs for mixsiar as cosimmr
plot_mixsiar = function(input, source_names){
  
 
  out_all_p = input[,1,]
  
  
  colnames(out_all_p) = source_names
  
  df <- reshape2::melt(out_all_p)
  
  
  colnames(df) = c("Num", "Source", "Proportion")
  
 
    g <- ggplot(df, aes(
      x = Proportion,
      fill = Source
    )) +
      scale_fill_viridis(discrete = TRUE) +
      geom_histogram(binwidth = 0.05, alpha = 0.5) +
      theme_bw() +
      ggtitle("Proportions: Observation 1") +
      facet_wrap("~ Source") +
      theme(legend.position = "none")
    print(g) 
    
    
  
}