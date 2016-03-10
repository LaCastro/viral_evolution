r0 <- c(7, 6.3, 5.6, 4.9, 4.2, 3.5, 2.8 , 2.1, 1.4)

r0 <- c(1.5, 2.1, 2.8, 3.5, 4.2, 4.9, 5.6, 6.3, 7)

divergence.min <- c(0 ,0,	0,	0,	0,	0,	0,	0,	0)
divergence.max <- c(0.01,  0.02,	0.01,	0.01,	0.01,	0.01,	0.002222222,	0.002592593,	0.00975)

diversity.min <- c(0.000370482,  0.000177693,	0.00017418,	0.000139108,	0.000140816,	0.000106408,	0.000122269,	0.000115296,	0.000104155)
diversity.max <- c(0.005921141,  0.005487772,	0.00494898,	0.004990548,	0.004897959,	0.004444444,	0.001728395,	0.002106324,	0.00461875)

num.strains.min <- c(12,  30,	44,	48,	52,	57,	48,	54,	57)
num.strains.max <- c(25,  51,	71,	85,	89,	96,	91,	91,	94)

total.strains.min <- c(84,  150,	162,	182,	172,	185,	177,	188,	173)
total.strains.max <- c(136,  207,	240,	243,	246,	254,	259, 250,	258)

par(mfrow = c(2,2))

plot(r0, diversity.max, pch = 19, ylim = c(0, max(diversity.max)+.003), main = "Peak Diversity")
points(r0, diversity.min, pch = 19, col = "red")

plot(r0, num.strains.max, pch = 19, ylim = c(0, max(num.strains.max)+5), main = "Number of Strains at Peak")
points(r0, num.strains.min, pch = 19, col = "red")

plot(r0, total.strains.max, pch = 19, ylim = c(0, max(total.strains.max)+5), main = "Total Strains")
points(r0, total.strains.min, pch = 19, col = "red")

plot(r0, divergence.max, pch = 19, ylim = c(0,max(divergence.max)+.01), main = "End Divergence")
points(r0, divergence.min, pch = 19, col = "red")


