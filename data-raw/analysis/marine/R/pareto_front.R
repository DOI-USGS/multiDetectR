pareto_front <- function(richness_file, cost_grid) {
  fig_dir <- here::here("data-raw", "analysis", "marine", "figures")
  txt_dir <- here::here("data-raw", "analysis", "marine", "outputs")

  X_sites <- N_samples <- K_subsamples <- richness <- #nolint
    total_cost <- lbl <- NULL #nolint

  # Pareto front over the 200 designs: minimize total cost while maximizing
  # accumulated richness. Ported from cost_benefit_analysis_DL.R.
  richness_data <- readr::read_csv(richness_file, col_types = readr::cols())

  rich_by_design <-
    richness_data |>
    dplyr::group_by(X_sites, N_samples, K_subsamples) |>
    dplyr::summarize(richness = mean(richness, na.rm = TRUE),
                     .groups = "drop")

  df <-
    cost_grid |>
    dplyr::select(X_sites, N_samples, K_subsamples, total_cost) |>
    dplyr::inner_join(rich_by_design,
                      by = c("X_sites", "N_samples", "K_subsamples"))

  # A design is dominated if another achieves at least as much richness at no
  # greater cost, and is strictly better on one of the two objectives.
  is_dominated <- function(i, d) {
    any((d$richness >= d$richness[i]) &
          (d$total_cost <= d$total_cost[i]) &
          ((d$richness > d$richness[i]) | (d$total_cost < d$total_cost[i])))
  }
  dominated <- vapply(seq_len(nrow(df)), is_dominated, logical(1), d = df)
  front <-
    df[!dominated, ] |>
    dplyr::arrange(total_cost)

  # Knee point: the non-dominated design with the greatest gap between the
  # min-max-normalized richness and cost (where returns begin to diminish).
  norm01 <- function(z) (z - min(z)) / (max(z) - min(z))
  knee_idx <- which.max(norm01(front$richness) - norm01(front$total_cost))
  knee <- front[knee_idx, ]

  front_not_knee <-
    front |>
    dplyr::mutate(lbl = sprintf("X=%s, N=%s, K=%s",
                                X_sites, N_samples, K_subsamples)) |>
    dplyr::anti_join(dplyr::select(knee, X_sites, N_samples, K_subsamples),
                     by = c("X_sites", "N_samples", "K_subsamples"))

  extended_cols <-
    colorspace::qualitative_hcl(length(unique(df$X_sites)),
                                palette = "Dynamic")

  pareto_plot <-
    df |>
    ggplot2::ggplot(ggplot2::aes(total_cost, richness,
                                 color = factor(X_sites))) +
    ggplot2::geom_point(alpha = 0.35, size = 3) +
    ggplot2::geom_line(data = front, ggplot2::aes(group = 1),
                       color = "red", linewidth = 0.8) +
    ggplot2::geom_point(data = front, color = "red", size = 5) +
    ggplot2::geom_point(data = knee, shape = 23, fill = "gold",
                        color = "black", size = 8) +
    ggrepel::geom_text_repel(
      data = front_not_knee,
      ggplot2::aes(x = total_cost, y = richness, label = lbl),
      inherit.aes = FALSE, size = 4, fontface = "bold", color = "black",
      box.padding = 0.2, point.padding = 0.15, min.segment.length = 0,
      max.overlaps = 100) +
    ggrepel::geom_label_repel(
      data = knee,
      ggplot2::aes(x = total_cost, y = richness,
                   label = sprintf("X=%s, N=%s, K=%s",
                                   X_sites, N_samples, K_subsamples)),
      inherit.aes = FALSE, size = 5, color = "black", fill = "gold",
      nudge_x = -400, nudge_y = 3, min.segment.length = 0) +
    ggplot2::scale_color_manual(values = extended_cols, name = "Sites (X)") +
    ggplot2::theme_bw() +
    ggplot2::xlab("Total cost (US$)") +
    ggplot2::ylab("Accumulated marine species richness")

  ggplot2::ggsave(filename = here::here(fig_dir, "pareto_front.jpg"),
                  plot = pareto_plot,
                  dpi = 300,
                  width = 8, height = 8)

  readr::write_csv(x = front,
                   file = file.path(txt_dir, "pareto_front.csv"))
  readr::write_csv(x = knee,
                   file = file.path(txt_dir, "pareto_knee.csv"))

  pareto_plot
}
