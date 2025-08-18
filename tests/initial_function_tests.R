# Tests
library(readr)
library(tidyverse)
dhat <- read_csv("test_files/dhat.csv", col_names = FALSE) %>% as.matrix()
Vhat <- read_csv("test_files/vhat.csv", col_names = FALSE) %>% as.matrix()
delta <- read_csv("test_files/delta.csv", col_names = FALSE) %>% as.matrix()

p <- length(delta)
Xtmp <- matrix(1, nrow = p)
result <- MDproj2(delta, Vhat, Xtmp)
degree <- 1
Xtmp <- cbind(Xtmp, (1:p)^degree)

result <- MDproj2(delta, Vhat, Xtmp)
result$d
result$Vd
result$MD
result$pv

# MDproj2 MATCHES matlab version with degree = 0, 1



# MD proj ltf2
kk = 100

bestlam1=-1;
bestlam2=-1;
mbhsf=matrix(0, nrow = kk)
target_df = p - 1

# Same as matlab
diff_df(-10, 10, 1, V, target_df)

# Same as matlab for this specific example
lam_bounds <- find_lam_bounds(p, V, target_df)

result <- full_eventplot_l2tf(delta, Vhat, nsim = 100)

