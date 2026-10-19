library(tidyverse)
library(ggthemes)

model_estimates_2_may <- read_csv(file = "model_estimates_2_may.csv")
model_estimates_3_may <- read_csv(file = "model_estimates_3_may.csv")

model_estimates <-
  bind_rows(model_estimates_2_may,
            model_estimates_3_may) |>
  mutate(level = gsub("delta_|alpha_|beta_|\\[\\d+\\]", "", parameter),
         level = ifelse(level == "p", "italic(p)", level),
         level = factor(level, levels = c("italic(p)", "theta", "psi")),
         samples = case_when(data == "three-level" ~ "10",
                             data == "two-level" ~ "2"),
         samples = factor(samples, levels = c("2", "10")))

dodge_width <- 0.5

model_estimates_plot <-
  model_estimates |>
  filter(!grepl("lp__", parameter)) |>
  ggplot(aes(x = species,
             y = median,
             color = samples)) +
  geom_point(size = 2, position = position_dodge(width = dodge_width)) +
  geom_linerange(aes(ymin = l95, ymax = u95),
                 position = position_dodge(width = dodge_width)) +
  geom_linerange(aes(ymin = l95, ymax = u95),
                 linewidth = 1.5,
                 position = position_dodge(width = dodge_width)) +
  coord_flip() +
  facet_grid(cols = vars(level),
             labeller = label_parsed) +
  theme_bw() +
  theme(strip.background = element_blank()) +
  scale_color_colorblind("Number\nof\nPCR\nsamples") +
  ylab("Species") +
  xlab("Estimated probability and 95% credible interval")
model_estimates_plot
ggsave("model_estimates_plot.jpg", model_estimates_plot, width = 8, height = 8)

## Raw data
dat_may_2_raw <-
  read_csv("ptconception_occupancy_2level_may.csv") |>
  rename(unit = Site,
         sample = FieldReplicate) |>
  mutate(consensus = Consensus,
         y = ifelse(value > 0, 1, 0)) |>
  separate_wider_delim(Consensus, delim = "/",
                       names_sep = "", too_few = "align_start") |>
  select(unit, sample, PCRReplicate,
         species = Consensus1, consensus, y) |>
  ungroup()

dat_may_2 <-
  dat_may_2_raw |>
  group_by(unit, sample, species, consensus) |>
  summarize(y = sum(y), n_agg = n(), k_size = n(),
            .groups = "drop") |>
  mutate(sample_occ = ifelse(y > 0, 1, 0))

dat_may_3_raw <-
  read_csv("ptconception_occupancy_3level_may.csv") |>
  rename(unit = Site,
         sample = FieldReplicate) |>
  mutate(consensus = Consensus,
         y = ifelse(value > 0, 1, 0)) |>
  separate_wider_delim(Consensus, delim = "/",
                       names_sep = "", too_few = "align_start") |>
  select(unit, sample, PCRReplicate,
         species = Consensus1, consensus, y) |>
  ungroup()

dat_may_3 <-
  dat_may_3_raw |> 
  group_by(unit, sample, species, consensus) |>
  summarize(y = sum(y), n_agg = n(), k_size = n(),
            .groups = "drop") |>
  mutate(sample_occ = ifelse(y > 0, 1, 0))

dat_may_2 |>
  group_by(sample) |>
  summarize(max_k = max(k_size),
            max_y = max(y))

dat_may_3 |>
  group_by(sample) |>
  summarize(max_k = max(k_size),
            max_y = max(y))

dat_may_2 |>
  ggplot(aes(x = sample_occ)) +
  geom_histogram()

dat_may_2 |>
  group_by(sample, species) |>
  summarize(mean_y = mean(y),
            sum_occ = sum(sample_occ),
            n_samples = n())

dat_may_2 |>
  group_by(species) |>
  summarize(total_y = sum(y)) |>
  ggplot(aes(x = total_y)) +
  geom_histogram(binwidth = 1) +
  scale_x_continuous(breaks = seq(0, 45, by = 1))

## Look at psi and theta and p
model_estimates |>
  filter(!grepl("lp__", parameter)) |>
  ggplot(aes(x  = median)) +
  geom_histogram() +
  facet_grid(cols = vars(level),
             rows = vars(samples),
             labeller = label_parsed) +
  theme_bw() +
  theme(strip.background = element_blank()) +
  scale_x_continuous(lim = c(0, 1))

model_estimates |>
  filter(!grepl("lp__", parameter)) |>
  ggplot(aes(x = median, color = samples, fill = samples)) +
  facet_grid(cols = vars(level),
             labeller = label_parsed) +
  theme_bw() +
  theme(strip.background = element_blank()) +
  scale_x_continuous(lim = c(0, 1)) +
  geom_density(lwd = 1, alpha = 0.25) +
  scale_fill_colorblind("Number\nof\nsamples") +
  scale_color_colorblind("Number\nof\nsamples") +
  ylab("Median of posterior distribution")

## Plot simulate communities

## three
curve_data_3 <-
  bind_rows(
    read_csv("curve_data_3_rep_t.csv") |>
      mutate(replacement = "With replacement"),
    read_csv("curve_data_3_rep_f.csv") |>
      mutate(replacement = "Without replacement")
  )

curve_data_3_summary <-
  curve_data_3 |>
  group_by(n_units, n_samples, k_subsamples, replacement) |>
  summarize(mean = mean(n_species),
            .groups = "drop",
            l95 = quantile(n_species, 0.025),
            u95 = quantile(n_species, 0.975))

spp_curve_three <-
  curve_data_3_summary |>
  mutate(k_subsamples = factor(k_subsamples)) |>
  ggplot(aes(x = n_units, y = mean,
             color = factor(n_samples),
             fill = factor(n_samples))) +
  geom_ribbon(aes(ymin = l95, ymax = u95), alpha = 0.25) +
  geom_line() +
  facet_grid(
    rows = vars(k_subsamples),
    cols = vars(replacement),
    labeller = label_bquote(rows = "K:" ~ .(k_subsamples))
  ) +
  scale_fill_colorblind("Number\nof\nsamples") +
  scale_color_colorblind("Number\nof\nsamples") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  scale_x_continuous(breaks = seq(1, 10, by = 1)) +
  xlab("Number of units per site") +
  ylab("Number of species")
spp_curve_three
ggsave("spp_curve_three.jpg", spp_curve_three, width = 8, height = 10)

curve_data_3_summary |>
  filter(n_units >= 7,
         k_subsamples >= 8) |>
  ggplot(aes(x = n_units, y = mean,
         color = factor(n_samples),
         fill = factor(n_samples)))  +
  geom_point() +
  geom_ribbon(aes(ymin = l95, ymax = u95), alpha = 0.25) +
  facet_grid(
    cols = vars(k_subsamples),
    rows = vars(replacement),
    labeller = label_bquote(cols = "K:" ~ .(k_subsamples))
  ) +
  scale_fill_colorblind("Number\nof\nsamples") +
  scale_color_colorblind("Number\nof\nsamples")

## plot 2
curve_data_2 <-
  bind_rows(
    read_csv("curve_data_2_rep_t.csv") |>
      mutate(replacement = "With replacement"),
    read_csv("curve_data_2_rep_f.csv") |>
      mutate(replacement = "Without replacement")
  )

curve_data_2_summary <-
  curve_data_2 |>
  group_by(n_units, n_samples, k_subsamples, replacement) |>
  summarize(mean = mean(n_species),
            .groups = "drop",
            l95 = quantile(n_species, 0.025),
            u95 = quantile(n_species, 0.975))

spp_curve_two <-
  curve_data_2_summary |>
  mutate(k_subsamples = factor(k_subsamples)) |>
  ggplot(aes(x = n_units, y = mean,
             color = factor(n_samples),
             fill = factor(n_samples))) +
  geom_ribbon(aes(ymin = l95, ymax = u95), alpha = 0.25) +
  geom_line() +
  facet_grid(
    rows = vars(k_subsamples),
    cols = vars(replacement),
    labeller = label_bquote(rows = "K:" ~ .(k_subsamples))
  ) +
  scale_fill_colorblind("Number\nof\nsamples") +
  scale_color_colorblind("Number\nof\nsamples") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  scale_x_continuous(breaks = seq(1, 10, by = 1)) +
  xlab("Number of units per site") +
  ylab("Number of species")
spp_curve_two
ggsave("spp_curve_two.jpg", spp_curve_two, width = 8, height = 4)

curve_data_2_summary |>
  filter(n_units >= 7) |>
  ggplot(aes(x = n_units, y = mean,
         color = factor(n_samples),
         fill = factor(n_samples)))  +
  geom_point() +
  geom_ribbon(aes(ymin = l95, ymax = u95), alpha = 0.25) +
  facet_grid(
    cols = vars(k_subsamples),
    rows = vars(replacement),
    labeller = label_bquote(cols = "Subsamples:" ~ .(k_subsamples))
  ) +
  scale_fill_colorblind("Number\nof\nsamples") +
  scale_color_colorblind("Number\nof\nsamples")

## Species observed