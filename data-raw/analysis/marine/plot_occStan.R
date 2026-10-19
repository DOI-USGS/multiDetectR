library(tidyverse)
library(ggthemes)

model_estimates <-
  bind_rows(
    read_csv("model_estimates_2_data_2_model.csv",
             show_col_types = FALSE),
    read_csv("model_estimates_2_data_3_model.csv",
             show_col_types = FALSE),
    read_csv("model_estimates_3_data_3_model.csv",
             show_col_types = FALSE)
  )

model_estimates |>
  filter(parameter == "lp__") |>
  ggplot(aes(x = Rhat)) +
  geom_histogram()

model_estimates |>
  filter(parameter == "lp__") |>
  ggplot(aes(x = n_eff)) +
  geom_histogram() +
  facet_grid(cols = vars(model), rows = vars(data))

## All data
model_estimates |>
  filter(!is.na(parameter), parameter != "lp__") |>
  ggplot(aes(x = species, y = prob, color = data)) +
  geom_point(position = position_dodge(width = 0.5)) +
  coord_flip() +
  facet_grid(cols = vars(par_plot),
             rows = vars(model),
             labeller = label_bquote(rows = "Model:" ~ .(model))) +
  geom_linerange(position = position_dodge(width = 0.75),
                 aes(ymin = l95, ymax = u95)) +
  geom_linerange(position = position_dodge(width = 0.75),
                 aes(ymin = l90, ymax = u90), linewidth = 1.5) +
  scale_color_colorblind("Data sequencing") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  xlab("Species") +
  ylab("Estimated probability")

model_estimates |>
  mutate(approach =
      case_when(model == "two-level" & data == "two-level" ~ "D 2 M 2",
                model == "three-level" & data == "three-level" ~ "D 3 M 3 ",
                model == "three-level" & data == "two-level" ~ "D 2 M 3")
  ) |>
  filter(!is.na(parameter), parameter != "lp__") |>
  ggplot(aes(x = species, y = prob, color = approach)) +
  geom_point(position = position_dodge(width = 0.5)) +
  coord_flip() +
  facet_grid(cols = vars(par_plot)) +
  geom_linerange(position = position_dodge(width = 0.75),
                 aes(ymin = l95, ymax = u95)) +
  geom_linerange(position = position_dodge(width = 0.75),
                 aes(ymin = l90, ymax = u90), linewidth = 1.5) +
  scale_color_colorblind("Data\nprocessing (D)\nand\nModeling\napproach (M)") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  xlab("Species") +
  ylab("Estimated probability")

two_level_spp <-
  model_estimates |>
  filter(!is.na(mean)) |>
  filter(data == "two-level") |>
  distinct(species) |>
  pull()

three_level_spp <-
  model_estimates |>
  filter(!is.na(mean)) |>
  filter(data == "three-level") |>
  distinct(species) |>
  pull()

three_level_spp %in% two_level_spp
three_level_spp[!three_level_spp %in% two_level_spp]

two_level_spp %in% three_level_spp
