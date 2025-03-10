library(dplyr)
library(ggplot2)

algo_levels <- c("kmeans", "kmodes", "kmedians", "kgower", "kGDM2",
                 "PAM:Gower", "PAM:GDM2",
                 "regularized normal", "multinomial", "binomial", "betabinomial")
algo_labels <- c("kmeans", "kmodes", "kmedians", "kGower", "kGDM2",
                 "PAM:Gower", "PAM:GDM2",
                 "regularized normal", "multinomial", "binomial", "beta-binomial")
sim = Sys.glob("../data/sim_backpain_sample_size_ncat_*_alpha*.rds") %>%
  sapply(readRDS, simplify = F) %>%
  setNames(gsub('\\D|2', '', names(.))) %>% 
  bind_rows(.id='alpha') |> 
  filter(algo != 'kGDM2centOpt') |> #removing it as it is unfinished (30.10.24)
  mutate(algo=gsub('centMode', '', algo)) |> #renaming, same reason as above
  mutate(algo=gsub('om$', 'omial', algo),
         algo=gsub('regnorm', 'regularized normal', algo)) |> #also, renaming algos
  mutate(alpha=as.factor(as.numeric(alpha)),
         `response level length r`=as.factor(size+1),
         scale_treatment = case_when(algo %in% c('kmodes', 'multinomial') ~ 'nominalization',
                                     algo %in% c('regularized normal', 'kmeans') ~ 'numerical coding',
                                     TRUE ~ 'respecting ordinal scales'),
         scale_treatment=factor(scale_treatment,
                                levels=c('numerical coding',
                                         'nominalization',
                                         'respecting ordinal scales')),
         type=factor(ifelse(algo %in% c('betabinomial', 'binomial', 'multinomial',
                                        'regularized normal'), 'model-based',
                            'partitioning'), levels=c('partitioning', 'model-based'))) |>
  arrange(type, scale_treatment) |> 
  mutate(method = factor(algo, levels = algo_levels, labels = algo_labels)) |>
  rename(m=n_vars, ARI=ari) |> 
  arrange(alpha, size, N, m, algo, iter)

# ws <- wesanderson::wes_palette('Moonrise3', n=3) |> setNames(c('pink', 'brown', 'turquoise'))
ws_expanded <- list(
  pink = c('#F0808C', '#F4B5BD'),#'#D56B75'),
  brown = c('#8E8540', '#B5AF69'),
  turquoise = c('#43CD80', '#8EE5EE', '#A3E3ED', '#00688B', '#7BC9D6', '#74C8B1', '#60BCCC', '#76EEC6', '#3F8A7C', '#357885', '#4BA39F')
)

palette_by_sctrt <- unique(select(sim, method, scale_treatment)) |>
  na.omit() |>
  arrange(method) |> 
  mutate(categorycolor=ifelse(scale_treatment=='nominalization',
                              'pink', ifelse(scale_treatment=='numerical coding',
                                             'brown', 'turquoise')),
         color=NA)
for(i in names(ws_expanded)) {
  index <- palette_by_sctrt$categorycolor==i
  palette_by_sctrt$color[index] <- ws_expanded[[i]][seq_len(sum(index))]
}
palette_by_sctrt <- palette_by_sctrt$color |> 
  setNames(palette_by_sctrt$method)

plt <- \(dat) {
  ggplot(dat, aes(x=`response level length r`,
                  y=ARI, col=method)) +
    facet_grid(m~type+N, labeller=label_both) +
    geom_boxplot(outliers = F) + ylim(0,1) +
    scale_color_manual(values=palette_by_sctrt) +
    guides(color=guide_legend(ncol=7,
                              byrow=T)) +
    theme_bw() +
    theme(legend.position = 'bottom')
}

plt(filter(sim, alpha==0))
plt(filter(sim, alpha==75))
plt(filter(sim, alpha==150))
