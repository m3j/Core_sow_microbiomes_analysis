library(ggplot2)
library(ape)
source("code/Functions/stats_utils.R")
source("code/Functions/plot_utils.R")


# set theme
pcoa_theme <- theme(
    strip.text.x = element_text(size = 18, color = "gray30", face = "bold"),
    strip.background = element_rect(color = "gray20"),
    panel.grid = element_line(colour = "gray80"),
    panel.grid.major.x = element_blank(),
    axis.ticks.x = element_line(colour = "grey30"),
    axis.ticks.y = element_line(colour = "grey30"),
    axis.title.x = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 18, color = "grey30"),
    axis.text.y = element_text(size = 18, color = "grey30"),
    legend.title = element_text(size = 18, color = "black", hjust = 1),
    legend.text = element_text(size = 17, hjust = 1),
    plot.title = element_text(size = 26, color = "black"),
    plot.tag = element_text(size = 22, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "gray30"),
    legend.background = element_blank()
)

# function for plotting pcoa data
plot_pcoa <- function(pc_data, fill_var, shape_vector,
                      color_vector, label_vector, plt_title,
                      legend_title, legend_pos = c(0.78, 0), legend_order = NA,
                      variance_explained, perm_res = "", pc = c(1, 2)) {
    # get the pc values
    pc_x <- pc[1]
    pc_y <- pc[2]

    # define correct axis names
    x_axis <- paste0("Axis.", pc_x)
    y_axis <- paste0("Axis.", pc_y)

    # create plot with elipse and points
    pcoa_plot <- subset(pc_data, !is.na(pc_data[[fill_var]])) %>%
        ggplot(aes_string(x = x_axis, y = y_axis)) +
        stat_ellipse(aes_string(x = x_axis, y = y_axis, color = fill_var),
            fill = "gray",
            type = "norm",
            level = 0.75,
            alpha = 0.01,
            geom = "polygon",
            show.legend = FALSE
        ) +
        geom_point(
            aes_string(fill = fill_var, shape = fill_var),
            alpha = 0.8,
            size = 4
        ) +
        scale_shape_manual(values = shape_vector) +
        scale_color_manual(values = color_vector) +
        scale_fill_discrete(
            type = color_vector,
            name = legend_title,
            labels = label_vector
        ) +
        guides(
            fill = guide_legend(
                override.aes = list(shape = shape_vector),
                byrow = TRUE,
                direction = "vertical"
            ),
            shape = "none"
        ) +
        xlab(paste0("PC", pc_x, " (", round(variance_explained[pc_x], 2), "%)")) +
        ylab(paste0("PC", pc_y, " (", round(variance_explained[pc_y], 2), "%)")) +
        labs(title = plt_title) +
        pcoa_theme +
        theme(legend.position = legend_pos) +
        annotate("text",
            x = Inf,
            y = Inf,
            label = perm_res,
            col = "black",
            size = 6,
            vjust = 1.2,
            hjust = 1
        ) +
        coord_cartesian(clip = "off")

    return(pcoa_plot)
}

PCoA_plot_generation <- function(
    physeq_obj, dist, fill_var, shape_vector,
    color_vector, label_vector, plt_title,
    legend_title, legend_pos = c(0.78, 0), legend_order = NA,
    pc = c(1, 2)) {
    # get permanova results in string format
    perm_res <- perm_anova_res_string(physeq_obj, fill_var, dist)

    # do pcoa on input distance matrix
    pcoa_dist <- ape::pcoa(dist)

    # calculate variance explained used on plot axes
    variance_explained <- pcoa_dist$values$Eigenvalues / sum(pcoa_dist$values$Eigenvalues) * 100

    # get metadata
    metadata <- data.frame(sample_data(physeq_obj))

    # add PC values to the metadata
    metadata_pc <- cbind(metadata, pcoa_dist$vectors[, pc])

    # call plot function
    pcoa_fig <- plot_pcoa(
        pc_data = metadata_pc,
        fill_var = fill_var,
        shape_vector = shape_vector,
        color_vector = color_vector,
        label_vector = label_vector,
        plt_title = plt_title,
        legend_title = legend_title,
        legend_pos = legend_pos,
        legend_order = NA,
        variance_explained = variance_explained,
        perm_res = perm_res,
        pc = pc
    )

    return(pcoa_fig)
}
