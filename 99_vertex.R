bb <- coef(mod_dist_2016_soflo)

b <- c("mig" = bb[["dist_to_urb_2016_sc"]],
       "res" = bb[["dist_to_urb_2016_sc"]] + bb[["dist_to_urb_2016_sc:choiceresident"]])
a <- c("mig" = bb[["I(dist_to_urb_2016_sc^2)"]],
       "res" = bb[["I(dist_to_urb_2016_sc^2)"]] + bb[["I(dist_to_urb_2016_sc^2):choiceresident"]])
(vert <- -1*b/(2*a))

vert_real <- vert * dat$sd_dist_to_urb_2016[[1]] + dat$mean_dist_to_urb_2016[[1]]
