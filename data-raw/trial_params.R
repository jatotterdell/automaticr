# Design matrices
a <- 4
b <- 3
ab <- a*b
Xa <- kronecker(diag(1, a), rep(1, b))
Xb <- kronecker(rep(1, a), diag(1, b))
Xba <- kronecker(unique(Xa), unique(Xb))
X <- cbind(1, c(0, rep(1, ab)), rbind(0, Xa), rbind(0, Xb), rbind(0, Xba))
colnames(X) <- c("ctr", "trt", paste0("a", 1:a),
                 paste0("b", 1:b), c(outer(paste0("b", 1:b), paste0("a", 1:a), paste0)))
# Constraint matrices
Sa <- diag(1, a) - 1/a
Sb <- diag(1, b) - 1/b
Qa <- eigen(Sa)$vector[, -a]
Qb <- eigen(Sb)$vector[, -b]
Qba <- kronecker(Qa, Qb)
Q <- as.matrix(Matrix::bdiag(1, 1, Qa, Qb, Qba))
X_con <- X %*% Q

trial_design <- list(
  X = X_con,
  Q = Q,
  S = c(10, 5, rep(2.5, a + b - 2), rep(1, a*b-a-b+1))
)

interim_schedule <- data.frame(
  interim_num = 0:20,
  interim_n = seq(0, 1e4, 500),
  interim_new = c(NA, rep(500, 20)),
  alloc_seed = c(11871,
                 62871, 63579, 65231, 23841, 72859, 62859, 63511, 65213, 23897, 72801,
                 62873, 63577, 65233, 23843, 72857, 62857, 63513, 65217, 23899, 72803),
  stan_seed = c(11117,
                61871, 64579, 66231, 21841, 71859, 64511, 66211, 67213, 21897, 73801,
                61873, 64577, 66233, 21843, 71857, 64513, 66213, 67217, 21899, 73803)
)

intervention_map <- data.frame(
  randomisation_outcome = sprintf("%02d", 1:13),
  control = c(1, rep(0, 12)),
  message = c(0, rep(1:4, each = 3)),
  timing = c(0, rep(1:3, times = 4)),
  arm = 0:12
)


X = cbind(c(1, rep(0, 12)), c(0, rep(1, 12)))
colnames(X) <- c("ctr", "trt")
Z1 <- rbind(0, kronecker(diag(1, 4), rep(1, 3)))
colnames(Z1) <- paste0("m", 1:4)
Z2 <- rbind(0, kronecker(rep(1,4), diag(1,3)))
colnames(Z2) <- paste0("t", 1:3)
Z3 <- rbind(0, diag(1, 12))
colnames(Z3) <- paste0(rep(colnames(Z1), each = 3), rep(colnames(Z2), times = 4))
design_data <- list(
  K = 2,
  N = 13,
  L1 = 3,
  L2 = 4,
  L3 = 12,
  X = X,
  Z1 = Z1,
  Z2 = Z2,
  Z3 = Z3
)

design_matrix <- do.call(cbind, list(X, Z1, Z2, Z3))

usethis::use_data(
  trial_design,
  interim_schedule,
  design_data,
  design_matrix,
  intervention_map,
  internal = TRUE, overwrite = TRUE)
