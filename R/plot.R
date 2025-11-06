#' Plot `carmon` networks
#'
#' @description
#' `plot.carmon()` plots an annotated network from an object of `S3` class
#' `carmon`. Nodes from different omics layers will have different color coding.
#' To plot the network, it's enough to use `plot()` call on the `carmon` object.
#'
#' @param x The object of `S3` class `carmon`.
#' @param node_labels Show node names in the network. Defaults to TRUE.
#' @param hide_isolated Hide nodes that are not connected to any other node.
#'   Defaults to TRUE.
#' @param hot_nodes Highlight nodes that are selected as key nodes by the
#'   node centrality analysis, by drawing them bigger and adding the criteria
#'   for which they result to be central. Defaults to TRUE.
#' @param ... System required, not used here.
#'
#' @return Returns NULL, invisibly.
#' @export
#'
#' @examples
#' data(multi_omics_small)
#' c_obj <- carmon(multi_omics_small,
#'     net_method = "correlation",
#'     cor_quant = 0.25, plot = FALSE, verbose = FALSE
#' )
#' plot(c_obj)
#'
plot.carmon <- function(x, node_labels = TRUE, hide_isolated = TRUE,
                        hot_nodes = TRUE, ...) {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))
    graphics::par(mar = c(0, 0, 2, 0))
    if (all(igraph::degree(x$network) == 0)) {
        warning("The reconstructed network is empty, it cannot be plotted.")
        return(invisible(NULL))
    }
    if (!node_labels) {
        igraph::V(x$network)$label <- NA
    } else {
        igraph::V(x$network)$label <- parse(text = gsub(
            " |:|-", "_", igraph::V(x$network)$label))
    }
    res <- colorize_network(x)
    x$network <- res$net
    fillsframes <- res$f
    # Initialize as empty variable for having something to pass to plot_legend()
    quartiles <- NULL
    if (x$net_method != "mb") {
        # Make the network edges weighted
        res <- annotate_weights(x)
        x$network <- res$net
        quartiles <- res$q
    }
    lo <- igraph::layout_with_fr(x$network)
    if (hide_isolated) {
        disconnected <- which(igraph::degree(x$network) == 0)
        x$network <- igraph::delete_vertices(x$network, disconnected)
        if (length(disconnected)) lo <- lo[-disconnected, ]
    }
    igraph::V(x$network)$frame.width <- 5.5 / log(igraph::gorder(x$network) + 1)
    igraph::V(x$network)$size <- 40 / log(igraph::gorder(x$network) + 1)
    igraph::E(x$network)$width <- 4 / log(igraph::gsize(x$network) + 1)
    igraph::V(x$network)$label.degree <- pi / 4
    igraph::V(x$network)$label.dist <- (igraph::V(x$network)$size / 15) +
        (igraph::V(x$network)$frame.width / 2)
    igraph::V(x$network)$label.cex <- 0.75
    graphics::layout(matrix(seq_len(2), ncol = 1), heights = c(3, 1))
    if (hot_nodes) {
        x$network <- annotate_hot_nodes(x)
        plot(x$network, layout = lo,
            vertex.label = parse(text = igraph::V(x$network)$label))
    } else {
        plot(x$network, layout = lo)
    }
    plot_legend(x, fillsframes = fillsframes, quartiles = quartiles)
    return(invisible(NULL))
}

#' Annotate the network edges with the relative edge weight
#'
#' This function uses either partial correlation or Pearson's correlation,
#' depending on the network reconstruction method
#'
#' @param x The object of `S3` class `carmon` given as input to `plot.carmon()`.
#'
#' @returns It returns an R list of two elements: `net`, the igraph network
#'     object coming from the `carmon` object, updated with edge weights; and
#'     `q`, the quartiles of the weights distribution.
#'
#' @noRd
#'
annotate_weights <- function(x) {
    ws <- get_edge_weights(x)
    ws <- abs(ws)
    upper <- c(ws[upper.tri(ws)])
    lower <- c(t(ws)[upper.tri(ws)])
    max_ws <- vapply(seq_along(upper),
        function(i) max(upper[i], lower[i]), numeric(1))
    ws[upper.tri(ws, diag = FALSE)] <- max_ws
    ws <- t(ws)
    ws[upper.tri(ws, diag = FALSE)] <- max_ws
    igraph::E(x$network)$weight <- ws[igraph::as_edgelist(x$network)]

    # Edge type based on quartiles of weights: types 1,5,4,3
    sort_ws <- sort(max_ws, decreasing = TRUE)
    sort_ws <- sort_ws[sort_ws != 0]
    quartiles <- stats::quantile(sort_ws, prob = c(0.25, 0.5, 0.75))
    lty <- ws
    lty[which(ws > quartiles[3], arr.ind = TRUE)] <- 1
    lty[which(ws <= quartiles[3] & ws > quartiles[2], arr.ind = TRUE)] <- 5
    lty[which(ws <= quartiles[2] & ws > quartiles[1], arr.ind = TRUE)] <- 4
    lty[which(ws <= quartiles[1], arr.ind = TRUE)] <- 3
    igraph::E(x$network)$lty <- lty[igraph::as_edgelist(x$network)]

    return(list("net" = x$network, "q" = quartiles))
}

#' Assign color to network nodes based on the layer they belong to
#'
#' @param x The object of `S3` class `carmon` given as input to `plot.carmon()`.
#'
#' @returns  It returns an R list of two elements: `net`, the igraph network
#'     object coming from the `carmon` object, updated with the assigned node
#'     fills and frames colors; and `f`, containing a fillsframes object (see
#'     `get_fillsframes()`).
#'
#' @noRd
#'
colorize_network <- function(x) {
    p <- vapply(x$layers, ncol, numeric(1))
    p_tot <- sum(p)
    L <- length(p)
    fillsframes <- get_fillsframes(L)
    igraph::V(x$network)$color <- rep(0, p_tot)
    igraph::V(x$network)$frame.color <- rep(0, p_tot)
    igraph::V(x$network)$label.color <- rep(0, p_tot)
    i <- 1
    for (j in seq_len(p_tot)) {
        if (j > sum(p[seq_len(i)])) {
            i <- i + 1
        }
        igraph::V(x$network)$color[j] <- fillsframes[[1]][i]
        igraph::V(x$network)$frame.color[j] <- fillsframes[[2]][i]
        igraph::V(x$network)$label.color[j] <- fillsframes[[2]][i]
    }
    return(list("net" = x$network, "f" = fillsframes))
}

#' Enrich the network with the results from the centrality analysis
#'
#' @param x The object of `S3` class `carmon` given as input to `plot.carmon()`.
#'
#' @returns It returns an `igraph` network object, where the nodes that are
#'   selected as key nodes by the node centrality analysis are highlighted, by
#'   being drawn bigger and adding the initial of the measure for which they
#'   result to be central to the node label.
#'
#' @noRd
#'
annotate_hot_nodes <- function(x) {
    if (!is.null(x$report)) {
        all_measures <- x$measures_list
        for (i in seq_len(length(all_measures))) {
            central_nodes <- names(all_measures[[i]])
            vertexes <- which(names(igraph::V(x$network))
                        %in% central_nodes)
            igraph::V(x$network)$frame.width[vertexes] <-
                1.25 * igraph::V(x$network)$frame.width[vertexes]
            igraph::V(x$network)$size[vertexes] <-
                1.25 * igraph::V(x$network)$size[vertexes]
        }
        if (is.null(dim(x$report))) {
            name <- x$report[1]
            string <- x$report["central for"]
            v <- which(names(igraph::V(x$network)) == name)
            igraph::V(x$network)$label[v] <- parse(text = gsub(
                " |:|-", "_",
                paste("bold(\"", name, "\")^{bold(", string, ")}", sep = "")
                ))
            igraph::V(x$network)$label.dist[v] <-
                (igraph::V(x$network)$size[v] / 15) +
                (igraph::V(x$network)$frame.width[v] / 4)
        } else {
            for (i in seq_len(nrow(x$report))) {
                name <- rownames(x$report)[i]
                string <- x$report[i, "central for"]
                v <- which(names(igraph::V(x$network)) == name)
                igraph::V(x$network)$label[v] <-
                    parse(text = gsub(
                        " |:|-", "_",
                            paste("bold(\"", name, "\")^{bold(",
                                string, ")}",
                                    sep = ""
                            )
                        ))
                igraph::V(x$network)$label.dist[v] <-
                    (igraph::V(x$network)$size[v] / 15) +
                    (igraph::V(x$network)$frame.width[v] / 4)
            }
        }
        return(x$network)
    } else {
        warning("No centrality report detected.
    Central nodes will not be highlighted.")
        return(x$network)
    }
}

#' Plot the legend to the carmon network plot
#'
#' The legend has either one or two sections: only the node color legend if the
#' network is estimated with "mb", or also the edge weight quartiles with the
#' corresponding edge type annotation, from continuous to fragmented edge.
#'
#' @param x The object of `S3` class `carmon` given as input to `plot.carmon()`.
#' @param fillsframes A fillsframes object containing the node fills and frames
#'     color information. See `get_fillsframes()`.
#' @param quartiles The quartiles generated as the element `q` of the list
#'     returned by `annotate_weights()`.
#'
#' @returns It returns `NULL`, invisibly.
#'
#' @noRd
#'
plot_legend <- function(x, fillsframes, quartiles = NULL) {
    graphics::par(mar = c(2, 2, 0, 2))
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
    usr <- graphics::par("usr")
    width_1 <- max(graphics::strwidth(x$omics))
    pad <- graphics::strwidth("MM", units = "user")
    legend1_width <- width_1 + pad
    line_height <- graphics::strheight("M") * 1.5
    height_1 <- line_height * length(x$omics)
    total_width <- legend1_width
    total_height <- height_1
    updated <- list("total_width" = total_width, "total_height" = total_height)
    if (x$net_method != "mb") {
        updated <- update_legend_with_weights(x, quartiles, pad, line_height,
                                    legend1_width, height_1)
    }
    x_left <- (1 - updated$total_width) / 2
    x_right <- x_left + updated$total_width
    y_bottom <- (1 - updated$total_height) / 2
    y_top <- y_bottom + updated$total_height
    graphics::legend(x = x_left, y = y_top, cex = 0.85, legend = x$omics,
        pt.bg = fillsframes[[1]], col = fillsframes[[2]],
        pt.cex = 5/3, pch = 21, pt.lwd = 1.6,
        text.col = fillsframes[[2]], bty = "n")
    if (x$net_method != "mb") {
        graphics::legend(x = x_left + legend1_width + updated$spacing,
            y = y_top, col = "darkgray", cex = 0.85,
            legend = updated$quartiles_text, lty = c(1, 5, 4, 3),
            lwd = 2.5, bty = "n")
    }
    return(invisible(NULL))
}

#' Update the parameters of the legend in case the edges are weighted
#'
#' @param x The object of `S3` class `carmon` given as input to `plot.carmon()`.
#' @param quartiles The quartiles generated as the element `q` of the list
#'     returned by `annotate_weights()`.
#' @param pad The padding around the information in the legend.
#' @param line_height The height given to a line of th legend.
#' @param legend1_width The width of the legend if there was no weight
#'     information to display.
#' @param height_1 The height of the legend if there was no weight information
#'     to display.
#'
#' @returns It returns an R list. These are the elements:
#' \itemize {
#'     \item `quartiles_text` A vector with the 4 text pieces about the
#'         quantiles of the weights distribution.
#'     \item `spacing` The space between the node color information and the edge
#'         type description in the final legend.
#'     \item `total_width` The final total width of the legend.
#'     \item `total_height` The final total height of the legend.
#' }
#'
#' @noRd
#'
update_legend_with_weights <- function(x, quartiles, pad, line_height,
                                    legend1_width, height_1) {
    if (x$net_method %in% c("coglasso", "glasso")) {
        weight_string <- "p_cor"
    } else {
        weight_string <- "cor"
    }
    quartiles_text <- rep(0, 4)
    quartiles_text[1] <- paste("| ", weight_string, " |", " > ",
                                round(quartiles[3], digits = 4), sep = "")
    quartiles_text[2] <- paste(round(quartiles[2], digits = 4), " < ", "| ",
                                weight_string, " |", " <= ",
                                round(quartiles[3], digits = 4), sep = "")
    quartiles_text[3] <- paste(round(quartiles[1], digits = 4), " < ", "| ",
                                weight_string, " |", " <= ",
                                round(quartiles[2], digits = 4), sep = "")
    quartiles_text[4] <- paste("| ", weight_string, " |", " <= ",
                                round(quartiles[1], digits = 4), sep = "")
    width_2 <- max(graphics::strwidth(quartiles_text))
    legend2_width <- width_2 + pad
    height_2 <- line_height * 4
    spacing <- 0.05
    total_width <- legend1_width + legend2_width + spacing
    total_height <- max(height_1, height_2)
    return(list("quartiles_text" = quartiles_text, "spacing" = spacing,
                "total_width" = total_width, "total_height" = total_height))
}

#' Plot the results of the node centrality analysis
#'
#' @param carmon_obj The object of `S3` class `carmon`.
#' @param scaled Whether the centrality measures reported by the centrality
#'   analysis are 0-1 scaled or not. Defaults to TRUE.
#'
#' @returns Returns NULL, invisibly.
#' @export
#'
#' @examples
#' # Let's build and analyse a carmon network:
#' data(multi_omics)
#' c_obj <- carmon(multi_omics,
#'     net_method = "correlation",
#'     cor_quant = 0.05, analysis = TRUE, plot = FALSE,
#'     # analysis is already TRUE by default
#'     verbose = FALSE
#' )
#' # To plot the results of the centrality analysis:
#' plot_report(c_obj)
#'
plot_report <- function(carmon_obj, scaled = TRUE) {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))

    if (all(igraph::degree(carmon_obj$network) == 0)) {
        warning("The reconstructed network is empty, cannot plot the results of
    the centrality analysis.")
        return(invisible(NULL))
    }
    if (is.null(carmon_obj$report)) {
        warning("The centrality report was not compiled yet.
    Compiling it now, but it will not be saved in the carmon object.")
        carmon_obj <- compute_centrality(carmon_obj)
        if (is.null(carmon_obj$measures_list)) {
            warning("No centrality report detected. Cannot plot the results of
    centrality analysis.")
            return(invisible(NULL))
        }
    }

    p <- vapply(carmon_obj$layers, ncol, numeric(1))
    p_tot <- sum(p)
    L <- length(p)
    # Getting legend color same as corresponding node
    fillsframes <- get_fillsframes(L)

    all_measures <- carmon_obj$measures_list
    if (length(all_measures) > 1) {
        if (length(all_measures) > 2) {
            # four plots
            graphics::par(mfrow = c(2, 2))
        } else {
            # double plot
            graphics::par(mfrow = c(1, 2))
        }
    }

    for (i in seq_len(length(all_measures))) {
        plot_report_measure(carmon_obj, scaled, all_measures, i, fillsframes, L)
    }

    return(invisible(NULL))
}

plot_report_measure <- function(carmon_obj, scaled, all_measures, i,
                                fillsframes, L) {
    nodes <- factor(names(all_measures[[i]]), levels = names(all_measures[[i]]))
    if (scaled | max(all_measures[[i]]) <= 1) {
        xlim <- c(0, 1)
    } else {
        xlim <- c(0, max(all_measures[[i]]))
    }
    y_pos <- rev(seq_along(nodes))
    # Plot title root
    if (scaled) {
        root_main <- "Scaled "
    } else {
        root_main <- ""
    }
    plot(NULL, xlim = xlim, ylim = c(0.5, length(nodes) + 0.5), yaxt = "n",
            xlab = names(all_measures)[i], ylab = "",
            main = paste(root_main, names(all_measures)[i], sep = ""))
    graphics::axis(2, at = y_pos, labels = nodes, las = 2)
    bar_colors <- rep(0, length(nodes))
    border_colors <- rep(0, length(nodes))
    for (j in seq_len(length(nodes))) {
        for (k in seq_len(L)) {
            if (nodes[j] %in% colnames(carmon_obj$layers[[k]])) {
                bar_colors[j] <- fillsframes[[1]][k]
                border_colors[j] <- fillsframes[[2]][k]
            }
        }
    }
    main_bar_height <- 0.75
    for (j in seq_along(nodes)) {
        y <- y_pos[j]
        graphics::rect(0, y - main_bar_height / 2, all_measures[[i]][j],
                    y + main_bar_height / 2, col = bar_colors[j], border = NA)
        x_left <- 0
        x_right <- all_measures[[i]][j]
        y_bottom <- y - main_bar_height / 2
        y_top <- y + main_bar_height / 2
        graphics::segments(x_left, y_top, x_right, y_top,
                            col = border_colors[j], lwd = 1.5)
        graphics::segments(x_right, y_bottom, x_right, y_top,
                            col = border_colors[j], lwd = 1.5)
        graphics::segments(x_left, y_bottom, x_right, y_bottom,
                            col = border_colors[j], lwd = 1.5)
    }
    graphics::legend("bottomright", legend = carmon_obj$omics,
                fill = fillsframes[[1]], border = fillsframes[[2]], cex = 0.7,
                pt.cex = 0.3, y.intersp = 0.8)
}

#' Extract colors for nodes and frames to color code variables from different
#' omics layers
#'
#' @param num_layers The number of layers. Up to six, the returned colors
#'   are part of a curated palette.
#'
#' @return `get_fillsframes()` returns a list of two paired vectors, with as
#'   many elements as `num_layers`, one with node colors and the other one
#'   with the relative darker node frame colors.
#'
#' @noRd
#'
get_fillsframes <- function(num_layers) {
    fills <- c(
        "#00ccff", "#ff9999", "#F2F79E", "#7AD988", "#AFA2FF",
        "#A63A50"
    )
    frames <- c(
        "#002060", "#800000", "#66684d", "#0c3b12", "#48417c",
        "#5e2d38"
    )
    if (num_layers <= length(fills)) {
        return(list(fills[seq_len(num_layers)], frames[seq_len(num_layers)]))
    } else {
        add_fills <- grDevices::rainbow(num_layers - length(fills))
        add_frames <- vapply(add_fills, get_frame, character(1),
            USE.NAMES = FALSE
        )
        return(list(c(fills, add_fills), c(frames, add_frames)))
    }
}

#' Get a frame color relative to a given node color
#'
#' @param color The node color.
#'
#' @return Returns the frame color relative to the input node color, a tint.
#'
#' @noRd
#'
get_frame <- function(color) {
    (grDevices::colorRampPalette(colors = c(color, "black")))(11)[8]
}
