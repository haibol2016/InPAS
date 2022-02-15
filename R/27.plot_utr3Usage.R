#' Visualize the dPDUI events using ggplot2
#'
#' Visualize the dPDUI events by plotting the MSE, and total coverage per group
#' along 3' UTR regions with dPDUI using [ggplot2::geom_line ()].
#'
#' @param usage_data An object of [GenomicRanges::GRanges-class], an output from
#'   [get_usage4plot()].
#' @param vline_color color for vertical line showing position of predicated
#'   proximal CP site. Default, purple.
#' @param vline_type line type for vertical line showing position of predicated
#'   proximal CP site. Default, dashed. See \href{https://cran.r-project.org/web/packages/ggplot2/vignettes/ggplot2-specs.html}{ggplot2 linetype}.
#'
#' @return A ggplot object for refined plotting
#' @importFrom ggplot2 ggplot geom_line geom_vline facet_grid ylab ggtitle theme
#'   element_text aes vars
#' @export
#' @author Haibo Liu
#' @seealso For example, see [get_usage4plot()].

plot_utr3Usage <- function(usage_data,
                           vline_color = "purple",
                           vline_type = "dashed") {
  p <- list()
  for (gr in names(usage_data)) {
    data <- usage_data[names(usage_data) == gr, ]$dat[[1]]
    p[[gr]] <- ggplot(data = data, aes(
      x = Position, y = value,
      colour = Coverage
    )) +
      geom_line() +
      geom_vline(
        xintercept =
          usage_data[names(usage_data) == gr, ]$offset,
        colour = vline_color, linetype = vline_type
      ) +
      facet_grid(rows = vars(Coverage), scales = "free_y") +
      ylab("") +
      ggtitle(gsub("\\.", ":", gr, perl = TRUE)) +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "none"
      )
    print(p[[gr]])
  }
  invisible(p)
}
