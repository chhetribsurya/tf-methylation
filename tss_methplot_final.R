
library(data.table)
library(ggplot2)
library(dplyr)


#####################
#####################
#####################


#file_list <- c(file1, file2, file3....)
file <- "/home/surya/Desktop/scripts/data/groupby_peaks_final/group_data_motifs_concating.txt"

plot_percent_meth <- function(file){

	binned_perc_meth_table <- fread(file, sep="\t", header= TRUE)
	names(binned_perc_meth_table) <- c("bin_start", "bin_end", "meth", "unmeth", "Percent_meth", "annotation")
	binned_perc_meth_table$bin_mid <- with(binned_perc_meth_table, (bin_start + bin_end)/2)
    this_plot <- ggplot2::ggplot(binned_perc_meth_table, aes(bin_mid, Percent_meth)) +
    ggplot2::geom_line(aes(color=annotation)) +
    ggplot2::ylab("Percent Methylation") +
    ggplot2::xlab("Distance from Center (bp)") +
    ggplot2::scale_x_continuous(expand=c(0,0)) +
    ggplot2::scale_y_continuous(expand=c(0,0.0001), limits = c(0,1)) +
    ggplot2::theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          # legend.position="none",
          legend.key = element_blank(),
          axis.line = element_line(),
          panel.background = element_blank(),
          axis.text.x = element_text(color="black", size=12),
          axis.title.x = element_text(color="black", size=14, vjust=-1.5),
          axis.text.y = element_text(color="black", size=12),
          axis.title.y = element_text(color="black", size=14, vjust=3),
          plot.margin = grid::unit(c(1,1,1,1), "cm"),
          panel.border=element_blank(),
          axis.ticks=element_line(size=0.6, color="black"))

  this_plot %>% return
}


meth_plot <- plot_percent_meth(file)
meth_plot.save(dir + "")


#######################
#######################
#######################



read_meth_data <- fread("/home/surya/Desktop/scripts/data/groupby_peaks_final/group_data_motifs_concating.txt", sep="\t", header=TRUE)
read_meth_data %>% head
names(read_meth_data) <- c("bin_start", "bin_end", "meth", "unmeth", "Percent_meth", "annotation")

theme_1 <-theme(panel.grid.minor = element_blank(),
			panel.grid.major = element_blank(),
			# legend.position="none",
			legend.key = element_blank(),
			axis.line = element_line(),
			panel.background = element_blank(),
			axis.text.x = element_text(color="black", size=12),
			axis.title.x = element_text(color="black", size=14, vjust=-1.5),
			axis.text.y = element_text(color="black", size=12),
			axis.title.y = element_text(color="black", size=14, vjust=3),
			plot.margin = grid::unit(c(1,1,1,1), "cm"),
			panel.border=element_blank(),
			axis.ticks=element_line(size=0.6, color="black"))


read_meth_data$bin_mid <- (read_meth_data$bin_start + read_meth_data$bin_end)/2
meth_plot <- ggplot2::ggplot(read_meth_data, aes(x=bin_mid, y=Percent_meth)) + 
			ggplot2::geom_line(aes(color=annotation)) +
    		ggplot2::ylab("Percent Methylation") +
    		ggplot2::xlab("Distance from Center (bp)") +
    		ggplot2::scale_x_continuous(expand=c(0,0)) +
    		ggplot2::scale_y_continuous(expand=c(0,0.0001), limits = c(0,1))
    		
meth_plot


plot_percent_meth_with_depth <- function(binned_perc_meth_table){
  binned_perc_meth_table$bin_mid <- with(binned_perc_meth_table,
                                         (bin_start + bin_end)/2)
  merged_df <- melt(binned_perc_meth_table, id.vars = c("bin_start", "bin_end", "meth", "unmeth", "group", "bin_mid"))
  this_plot <- ggplot2::ggplot(merged_df) +
    ggplot2::geom_line(aes(x = bin_mid, y = value, color = annotation)) +
    ggplot2::geom_line(aes(x = bin_mid, y = value, color = annotation)) +
    ggplot2::facet_grid(variable~., scales = "free_y")

  this_plot %>% return
}



#Loops can come in handy on numerous occasions. While loops are like repeated if statements; 
#the for loop is designed to iterate over all elements in a sequence.

#Also, whenever you're using a for loop, you might want to revise your code and see whether you can use the 
#lapply function instead. Learn all about this intuitive way of applying a function over a list or a 
#vector, and its variants sapply and vapply.


file1 <- "/home/surya/Desktop/scripts/data/groupby_peaks_final/group_data_motifs_concating.txt"
file_list <- c(file1)

plot_percent_meth <- function(file_list){

	for (file in file_list){
		#or use lapply
		print(file)
		binned_perc_meth_table <- fread(file, sep="\t", header= TRUE)
		names(binned_perc_meth_table) <- c("bin_start", "bin_end", "meth", "unmeth", "Percent_meth", "annotation")
		print(names(binned_perc_meth_table))
		binned_perc_meth_table$bin_mid <- with(binned_perc_meth_table, (bin_start + bin_end)/2)
	    this_plot <- ggplot2::ggplot(binned_perc_meth_table, aes(bin_mid, Percent_meth)) +
	    ggplot2::geom_line(aes(color=annotation)) +
	    ggplot2::ylab("Percent Methylation") +
	    ggplot2::xlab("Distance from Center (bp)") +
	    ggplot2::scale_x_continuous(expand=c(0,0)) +
	    ggplot2::scale_y_continuous(expand=c(0,0.0001), limits = c(0,1)) 
	    print("cool")

	    #this_plot.save(dir + "")
		return(this_plot)

  	}
}


meth_plot <- plot_percent_meth(file_list)
meth_plot
meth_plot.save(dir + "")


#Testing the loop:
for (num in c(1,2,3)){
	print(num)
}