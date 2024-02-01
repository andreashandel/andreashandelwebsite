# example from https://blog.djnavarro.net/posts/2023-06-10_pop-pk-models/

# data
warfpk <- readr::read_csv("warfpk.csv", na = ".", show_col_types = FALSE)
warfpk <- warfpk |> dplyr::rename(id = `#ID`)
warfpk <- warfpk |> 
  dplyr::mutate(
    id = id |> 
      stringr::str_remove_all("#") |> 
      as.numeric()
  )


warfpk_obs <- warfpk[warfpk$mdv == 0, ]
warfpk_amt <- warfpk[!is.na(warfpk$rate), ]

t_fit <- c(
  seq(.1, .9, .1),
  seq(1, 2.75, .25),
  seq(3, 9.5, .5),
  seq(10, 23, 1),
  seq(24, 120, 3)
)

dat <- list(
  n_ids = nrow(warfpk_amt),
  n_tot = nrow(warfpk_obs),
  n_obs = purrr::map_int(
    warfpk_amt$id,
    ~ nrow(warfpk_obs[warfpk_obs$id == .x, ])
  ),
  t_obs = warfpk_obs$time,
  c_obs = warfpk_obs$dv,
  dose = warfpk_amt$amt,
  t_fit = t_fit,
  n_fit = length(t_fit)
)

mod <- cmdstanr::cmdstan_model("testmodel.stan")
out <- mod$sample(
  data = dat,
  chains = 4,
  refresh = 1,
  iter_warmup = 1000,
  iter_sampling = 1000
)
