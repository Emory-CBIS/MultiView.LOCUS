.mv_phi_for_view <- function(lambda_ch, k, n_views) {
  if (length(lambda_ch) == 1L) {
    return(as.numeric(lambda_ch))
  }
  if (length(lambda_ch) != n_views) {
    stop("lambda_ch must have length 1 or the number of views.")
  }
  as.numeric(lambda_ch[k])
}

.mv_threshold_signal <- function(y, phi, penalt, gamma) {
  y <- drop(y)
  if (is.null(penalt) || identical(penalt, "NULL")) {
    return(y)
  }

  threshold <- phi / 2
  if (penalt == "L1") {
    return(sign(y) * pmax(abs(y) - threshold, 0))
  }
  if (penalt == "Hardthreshold") {
    return(y * (abs(y) >= threshold))
  }
  if (penalt == "SCAD") {
    return(SCAD_func(y, lambda_ch = threshold, gamma = gamma))
  }
  stop("No Penalty available!")
}

.mv_solve_normal_equations <- function(design, response, label) {
  gram <- crossprod(design)
  rhs <- crossprod(design, response)
  solution <- tryCatch(
    solve(gram, rhs),
    error = function(e) MASS::ginv(gram) %*% rhs
  )
  solution <- drop(solution)
  if (any(!is.finite(solution))) {
    stop(paste("Non-finite solution while updating", label))
  }
  solution
}

.mv_source_design <- function(X) {
  rank <- ncol(X)
  p <- nrow(X) * (nrow(X) - 1L) / 2L
  Z <- vapply(
    seq_len(rank),
    function(r) Ltrans(tcrossprod(X[, r]), FALSE),
    numeric(p)
  )
  matrix(Z, nrow = p, ncol = rank)
}

.mv_algorithm_s1_update <- function(Y, A, theta_common, theta_spe, q,
                                    q_common, psi, penalt, lambda_ch, gamma,
                                    H_inv) {
  n_views <- length(Y)
  p <- ncol(Y[[1]])
  V <- as.integer(round((sqrt(1 + 8 * p) + 1) / 2))
  if (V * (V - 1) / 2 != p) {
    stop("The number of connectivity edges is inconsistent with V(V-1)/2.")
  }

  edge_pairs <- which(
    upper.tri(matrix(FALSE, V, V), diag = FALSE),
    arr.ind = TRUE
  )
  edge_index <- lapply(
    seq_len(V),
    function(v) which(edge_pairs[, 1] == v | edge_pairs[, 2] == v)
  )

  theta_new <- theta_common
  theta_spe_new <- theta_spe

  ## Step 1: node-rotation updates of the latent coordinates.
  if (q_common > 0L) {
    for (l in seq_len(q_common)) {
      X <- t(theta_common[[1]][[l]]$J_l)
      rank <- ncol(X)

      for (v in seq_len(V)) {
        gram_sum <- matrix(0, rank, rank)
        rhs_sum <- numeric(rank)

        for (k in seq_len(n_views)) {
          d_old <- theta_common[[k]][[l]]$lam_l
          W <- X[-v, , drop = FALSE] %*%
            diag(d_old, nrow = length(d_old), ncol = length(d_old))
          target <- drop(
            t(Y[[k]][, edge_index[[v]], drop = FALSE]) %*% A[[k]][, l]
          )
          phi_k <- .mv_phi_for_view(lambda_ch, k, n_views)
          b_hat <- .mv_threshold_signal(target, phi_k, penalt, gamma)
          gram_sum <- gram_sum + crossprod(W)
          rhs_sum <- rhs_sum + drop(crossprod(W, b_hat))
        }

        x_update <- tryCatch(
          solve(gram_sum, rhs_sum),
          error = function(e) MASS::ginv(gram_sum) %*% rhs_sum
        )
        X[v, ] <- drop(x_update)
      }

      for (k in seq_len(n_views)) {
        theta_new[[k]][[l]]$J_l <- t(X)
      }
    }
  }

  for (k in seq_len(n_views)) {
    if (q[k] > 0L) {
      for (r in seq_len(q[k])) {
        l <- q_common + r
        X <- t(theta_spe[[k]][[r]]$J_l)
        d_old <- theta_spe[[k]][[r]]$lam_l

        for (v in seq_len(V)) {
          W <- X[-v, , drop = FALSE] %*%
            diag(d_old, nrow = length(d_old), ncol = length(d_old))
          target <- drop(
            t(Y[[k]][, edge_index[[v]], drop = FALSE]) %*% A[[k]][, l]
          )
          phi_k <- .mv_phi_for_view(lambda_ch, k, n_views)
          b_hat <- .mv_threshold_signal(target, phi_k, penalt, gamma)
          X[v, ] <- .mv_solve_normal_equations(
            W, b_hat,
            paste("X for view", k, "and component", l)
          )
        }

        theta_spe_new[[k]][[r]]$J_l <- t(X)
      }
    }
  }

  ## Step 2: update each diagonal D through thresholding and projection.
  for (k in seq_len(n_views)) {
    phi_k <- .mv_phi_for_view(lambda_ch, k, n_views)

    if (q_common > 0L) {
      for (l in seq_len(q_common)) {
        X <- t(theta_new[[k]][[l]]$J_l)
        Z <- .mv_source_design(X)
        target <- drop(t(Y[[k]]) %*% A[[k]][, l])
        b_hat <- .mv_threshold_signal(target, phi_k, penalt, gamma)
        theta_new[[k]][[l]]$lam_l <- .mv_solve_normal_equations(
          Z, b_hat,
          paste("D for view", k, "and common component", l)
        )
      }
    }

    if (q[k] > 0L) {
      for (r in seq_len(q[k])) {
        l <- q_common + r
        X <- t(theta_spe_new[[k]][[r]]$J_l)
        Z <- .mv_source_design(X)
        target <- drop(t(Y[[k]]) %*% A[[k]][, l])
        b_hat <- .mv_threshold_signal(target, phi_k, penalt, gamma)
        theta_spe_new[[k]][[r]]$lam_l <- .mv_solve_normal_equations(
          Z, b_hat,
          paste("D for view", k, "and view-specific component", l)
        )
      }
    }
  }

  ## Identifiability rescaling: ||x_r||=1 and d_r <- d_r ||x_r||^2.
  if (q_common > 0L) {
    for (l in seq_len(q_common)) {
      X <- t(theta_new[[1]][[l]]$J_l)
      x_norm <- sqrt(colSums(X^2))
      if (any(!is.finite(x_norm)) ||
          any(x_norm <= sqrt(.Machine$double.eps))) {
        stop(paste("Degenerate X column for common component", l))
      }
      X <- sweep(X, 2L, x_norm, "/")
      for (k in seq_len(n_views)) {
        theta_new[[k]][[l]]$J_l <- t(X)
        theta_new[[k]][[l]]$lam_l <-
          theta_new[[k]][[l]]$lam_l * x_norm^2
      }
    }
  }

  for (k in seq_len(n_views)) {
    if (q[k] > 0L) {
      for (r in seq_len(q[k])) {
        X <- t(theta_spe_new[[k]][[r]]$J_l)
        x_norm <- sqrt(colSums(X^2))
        if (any(!is.finite(x_norm)) ||
            any(x_norm <= sqrt(.Machine$double.eps))) {
          stop(paste(
            "Degenerate X column for view", k,
            "and view-specific component", r
          ))
        }
        theta_spe_new[[k]][[r]]$J_l <-
          t(sweep(X, 2L, x_norm, "/"))
        theta_spe_new[[k]][[r]]$lam_l <-
          theta_spe_new[[k]][[r]]$lam_l * x_norm^2
      }
    }
  }

  ## Reconstruct S from the updated and rescaled X and D.
  S <- vector("list", n_views)
  S_sparse <- vector("list", n_views)
  for (k in seq_len(n_views)) {
    q_total <- q_common + q[k]
    S[[k]] <- matrix(0, nrow = q_total, ncol = p)

    if (q_common > 0L) {
      for (l in seq_len(q_common)) {
        X <- t(theta_new[[k]][[l]]$J_l)
        Z <- .mv_source_design(X)
        S[[k]][l, ] <- drop(Z %*% theta_new[[k]][[l]]$lam_l)
      }
    }
    if (q[k] > 0L) {
      for (r in seq_len(q[k])) {
        l <- q_common + r
        X <- t(theta_spe_new[[k]][[r]]$J_l)
        Z <- .mv_source_design(X)
        S[[k]][l, ] <- drop(Z %*% theta_spe_new[[k]][[r]]$lam_l)
      }
    }

    phi_k <- .mv_phi_for_view(lambda_ch, k, n_views)
    S_sparse[[k]] <- S[[k]] * (abs(S[[k]]) >= 3.5 * phi_k)
  }

  ## Step 3: sequential constrained update of every column of A-tilde.
  A_new <- vector("list", n_views)
  for (k in seq_len(n_views)) {
    q_total <- q_common + q[k]
    A_new[[k]] <- matrix(0, q_total, q_total)
    Dmat <- Y[[k]] %*% t(Y[[k]])

    for (l in seq_len(q_total)) {
      dvec <- drop(S[[k]][l, ] %*% t(Y[[k]]))
      if (l <= q_common) {
        for (k_prime in setdiff(seq_len(n_views), k)) {
          dvec <- dvec + (psi / 2) * drop(
            t(A[[k_prime]][, l]) %*%
              t(H_inv[[k_prime]]) %*% H_inv[[k]]
          )
        }
      }

      if (l > 1L) {
        A_previous <- A_new[[k]][, seq_len(l - 1L), drop = FALSE]
        qp_fit <- quadprog::solve.QP(
          Dmat = Dmat,
          dvec = dvec,
          Amat = A_previous,
          bvec = rep(0, ncol(A_previous)),
          meq = ncol(A_previous)
        )
        a_update <- qp_fit$solution
      } else {
        a_update <- solve(Dmat, dvec)
      }

      a_norm <- sqrt(sum(a_update^2))
      if (!is.finite(a_norm) || a_norm <= sqrt(.Machine$double.eps)) {
        stop(paste(
          "Degenerate mixing-vector update for view", k,
          "and component", l
        ))
      }
      A_new[[k]][, l] <- a_update / a_norm
    }

    orthogonality_error <- max(abs(
      crossprod(A_new[[k]]) - diag(q_total)
    ))
    if (orthogonality_error > 1e-8) {
      warning(paste(
        "Orthogonality error for view", k, "is",
        signif(orthogonality_error, 4)
      ))
    }

    if (q_common > 0L) {
      for (l in seq_len(q_common)) {
        theta_new[[k]][[l]]$M_l <- A_new[[k]][, l]
      }
    }
    if (q[k] > 0L) {
      for (r in seq_len(q[k])) {
        theta_spe_new[[k]][[r]]$M_l <- A_new[[k]][, q_common + r]
      }
    }
  }

  list(
    A = A_new,
    S = S,
    S_sparse = S_sparse,
    theta_spe = theta_spe_new,
    theta_common = theta_new
  )
}

#' One Iteration Update (Approximate) for Joint Decomposition
#'
#' Updates the parameter lists returned by [joint_initial()] using an
#' efficient eigen-based approximation, optional sparsity penalties (SCAD,
#' L1, hard-threshold), and quadratic-programming updates of the mixing
#' matrix *A*.  Called internally by [multi_view_decomposition].
#'
#' @param Y            List of centred/whitened data matrices (see Details).
#' @param A            Current list of mixing matrices.
#' @param theta_common Current parameter list for shared subnetworks.
#' @param theta_spe    Current parameter list for view-specific subnetworks.
#' @param q,q_common   As in [joint_decomposition_multi()].
#' @param psi          Coupling strength for enforcing similarity of joint
#'                     eigenvalues across modalities.
#' @param penalt       Penalty type (`"SCAD"`, `"L1"`, `"Hardthreshold"`, or `NULL`).
#' @param eigen_cor    Minimum average correlation threshold (default `0.15`).
#' @param lambda_ch,gamma Sparsity-penalty hyper-parameters.
#' @param imput_method `"Previous"` or `"Average"`.  Imputation rule for
#'                     diagonal augmentation.
#' @param silent       Suppress console output if `TRUE`.
#' @param H_inv        List of whitening back-projection matrices.
#' @param Iter         Current outer-loop iteration number (integer).
#' @param cor_mutual   (Optional) pre-computed mutual correlations.
#' @param sequential_specific Logical; if `TRUE`, use the full Algorithm S1
#'                     updates: node-rotation updates of `X`, threshold-and-
#'                     projection updates of `D`, the stated identifiability
#'                     rescaling, reconstruction of `S`, and sequential
#'                     constrained updates of the reduced-space mixing matrix.
#'                     If `FALSE` (the default), use the original eigen-based
#'                     package updates.
#'
#' @return A `list` containing updated `A`, `S`, `S_sparse`,
#'         `theta_common`, and `theta_spe`.
#'
#' @keywords internal
#' @export

joint_update_approx <- function(Y,A,theta_common, theta_spe, q,q_common ,psi, penalt = NULL,eigen_cor = 0.15,lambda_ch , gamma = 3,imput_method = "Previous",silent = FALSE,H_inv ,Iter,cor_mutual, sequential_specific = FALSE)
{
  if (!is.logical(sequential_specific) || length(sequential_specific) != 1L ||
      is.na(sequential_specific)) {
    stop("sequential_specific must be TRUE or FALSE.")
  }

  if (sequential_specific) {
    return(.mv_algorithm_s1_update(
      Y = Y,
      A = A,
      theta_common = theta_common,
      theta_spe = theta_spe,
      q = q,
      q_common = q_common,
      psi = psi,
      penalt = penalt,
      lambda_ch = lambda_ch,
      gamma = gamma,
      H_inv = H_inv
    ))
  }

  # An extremely efficient approximation method with potentially higher performance.
  if(is.null(penalt))
  {
    if(!silent)
      cat("Joint low rank decomposition without penalty.")
  }else{
    if(!silent)
      cat(paste("Joint low rank decomposition", penalt,"penalty."))
  }
  theta_new = list()
  for (j in 1:length(q)){
      theta_new[[j]] = list()
      for (i in 1:q_common){
        theta_new[[j]][[i]] = list()
      }
  }

  K = dim(Y[[1]])[2]
  V = (sqrt(1+8*K)+1)/2
  N = dim(Y[[1]])[1]
  #############Common component approximately by eigen decomposition
  R = vector()
  for(curr_ic in 1:q_common)
  {
    sum_S = 0
    S_RC = list()
    for (j in 1:length(q)){
    theta_ic = theta_common[[j]][[curr_ic]]
    S_common = t(theta_ic$M_l%*%Y[[j]])
    R[curr_ic] = dim(theta_ic$J_l)[1]
    if(is.null(penalt))
    {
      S_new = S_common
    }else if(penalt == "SCAD")
    {
      if(gamma<=2){warning("Gamma needs to be > 2!");gamma = 2.01}
      S_new = SCAD_func(S_common,lambda_ch = lambda_ch  ,gamma = gamma)
      S_new = S_new /sd(S_new)*sd(S_common)
    }else if(penalt == "Hardthreshold")
    {
      S_new= S_common*(abs(S_common)>=lambda_ch)
      S_new = S_new/sd(S_new)*sd(S_common)
    }else if(penalt == "L1")
    {
      S_new = sign(S_common)*(abs(S_common)-lambda_ch)*(abs(S_common)>=lambda_ch)
      S_new = S_new /sd(S_new)*sd(S_common)
    }else
    {
      stop("No Penalty available!")
    }

    if(imput_method == "Previous"){
        Sl = Ltrinv(S_new, V,F) + diag(rep(max(abs(S_new)),V),V)
    }else if(imput_method == "Average"){
      Sl = Ltrinv(S_new,V,F) + diag( rep(mean(S_new),V ))
    }else{
      stop("No Imputation available!")
    }
    S_RC[[j]] = Sl
    if (eigen(Sl)$values[1]<0){
      sign_s = -1
    }else{sign_s = 1}
    sum_S = sum_S+ sign_s*Sl
    }
    eigenSl = eigen(sum_S)
    orderEigen = order(abs(eigenSl$values),decreasing = T)
    Rl = R[curr_ic]
    eigenset = orderEigen[1:Rl]

    for(k in 1:Rl)
    {
      theta_ic$J_l[k,]= eigenSl$vectors[,eigenset[k]]
    }
    sign = list()
    for (j in 1:length(q)){
        sign[[j]] = sign(diag((theta_ic$J_l)%*% S_RC[[j]] %*% t(theta_ic$J_l)))
    theta_new[[j]][[curr_ic]]$lam_l = diag((theta_ic$J_l) %*% S_RC[[j]] %*% t(theta_ic$J_l))
    theta_new[[j]][[curr_ic]]$J_l = theta_ic$J_l

      #sign_x*sqrt(apply((Slx %*% t(theta_ic$J_l))^2,2d,sum))/sd(sign_x*sqrt(apply((Slx %*% t(theta_ic$J_l))^2,2,sum)))*sd(theta_common[[curr_ic]]$lam_lx)
    }
    rep = 0
    combi = combn(1:length(q),2)
    cor_eigen_cur = 0
    inner = 0
    for (k in 1:ncol(combi)){
        pair = combi[,k]
        cor_eigen_cur = cor_eigen_cur+abs( cor(theta_new[[pair[1]]][[curr_ic]]$lam_l,theta_new[[pair[2]]][[curr_ic]]$lam_l))
        inner = inner + abs(theta_new[[pair[1]]][[curr_ic]]$lam_l*theta_new[[pair[2]]][[curr_ic]]$lam_l)
    }
    while (cor_eigen_cur< ncol(combi)*eigen_cor){
      rep = rep +1
      if (rep>1*Rl/3){
        break
      }
      if (cor_eigen_cur<0){
        for (j in 1:length(q)){
        theta_new[[j]][[curr_ic]]$lam_l[order(inner,decreasing = T)[1:rep]]=rnorm(rep,0,0.001)
      }}
      else{
          for (j in 1:length(q)){
        theta_new[[j]][[curr_ic]]$lam_l[order(inner,decreasing = F)[1:rep]]=rnorm(rep,0,0.001)
      }}
      cor_eigen_cur = 0
      for (k in 1:ncol(combi)){
          pair = combi[,k]
          cor_eigen_cur = cor_eigen_cur+abs( cor(theta_new[[pair[1]]][[curr_ic]]$lam_l,theta_new[[pair[2]]][[curr_ic]]$lam_l))
          inner = inner + theta_new[[pair[1]]][[curr_ic]]$lam_l*theta_new[[pair[2]]][[curr_ic]]$lam_l
      }
    }
      #sign_y*sqrt(apply((Sly %*% t(theta_ic$J_l))^2,2,sum))/sd(sign_y*sqrt(apply((Sly %*% t(theta_ic$J_l))^2,2,sum)))*sd(theta_common[[curr_ic]]$lam_ly)
    for (j in 1:length(q)){
    if( theta_new[[j]][[curr_ic]]$lam_l[1]<0 ) {theta_new[[1]][[curr_ic]]$lam_l = -1*theta_new[[1]][[curr_ic]]$lam_l} }
    #Check consistency of common components
    if (Iter%%40 ==0 ){
    for (j in 1:length(q)){
      cor_temp = 0
      comp = (1:length(q))[-j]
      for (k in comp){
          cor_temp = cor_temp + abs(cor(H_inv[[j]] %*% (A[[j]][,curr_ic]),H_inv[[k]] %*% (A[[k]][,curr_ic])))
        }
        if (cor_temp<0.3){
            S_temp =Ltrans(S_RC[[comp[1]]],F)
            S_new_temp= t(Y[[j]])%*%solve((Y[[j]])%*%t(Y[[j]]))%*%(Y[[j]]%*%S_temp)
            S_new_temp =  S_new_temp*( S_new_temp>abs(lambda_ch))
            eigen_temp = eigen(Ltrinv(S_new_temp,V,F)+ diag(apply(Ltrinv(S_new_temp,V,F),2, max)))
            theta_new[[j]][[curr_ic]]$lam_l =  eigen_temp$values[order(abs(eigen_temp$values),decreasing = T)[1:Rl]]
            theta_new[[j]][[curr_ic]]$J_l =  t(eigen_temp$vectors[,order(abs(eigen_temp$values),decreasing = T)[1:Rl]])
          }
                                    }
      }

    }
  #############X specific component approximately by eigen decomposition
theta_spe_new = list()
for (j in 1:length(q)){
    theta_spe_new[[j]] = list()
    for (i in 1:q[j]){
      theta_spe_new[[j]][[i]] = list()
    }
}

for (j in 1:length(q)){
  R_vec = c()
if(q[j] !=0){
  for(curr_ic in 1:q[j])
  {
    theta_ic = theta_spe[[j]][[curr_ic]]
    R_vec[curr_ic] = dim(theta_ic$J_l)[1]
    S = t(theta_ic$M_l%*%Y[[j]])
    if(is.null(penalt))
    {
      S_new = S
    }else if(penalt == "SCAD")
    {
      if(gamma<=2){warning("Gamma needs to be > 2!");gamma = 2.01}

      S_new =  SCAD_func(S,lambda_ch = lambda_ch  ,gamma = gamma)
      S_new = S_new/sd(S_new)*sd(S)

    }else if(penalt == "Hardthreshold")
    {
      S_new = S*(abs(S)>=lambda_ch)
      S_new = S_new /sd(S_new)*sd(S)

    }else if(penalt == "L1")
    {
      S_new = sign(S)*(abs(S)-lambda_ch)*(abs(S)>=lambda_ch)
      S_new = S_new /sd(S_new)*sd(S)
    }else
    {
      stop("No Penalty available!")
    }

    if(imput_method == "Previous"){
      Sl = Ltrinv(S_new,V,F) + diag(diag(t( theta_ic$J_l)%*%diag(theta_ic$lam_l)%*%theta_ic$J_l ))

    }else if(imput_method == "Average"){
      Sl = Ltrinv(S_new,V,F) + diag( rep(mean(S_new),V ))
    }else{
      stop("No Imputation available!")
    }
    eigenSl = eigen(Sl)
    orderEigen = order(abs(eigenSl$values),decreasing = T)
    Rl = R_vec[curr_ic]
    eigenset = orderEigen[1:Rl]

    for(k in 1:Rl)
    {
      theta_ic$J_l[k,]= eigenSl$vectors[,eigenset[k]]
    }

    theta_spe_new[[j]][[curr_ic]]$lam_l = eigenSl$values[eigenset]
    if( theta_spe_new[[j]][[curr_ic]]$lam_l[1]<0 ) {theta_spe_new[[j]][[curr_ic]]$lam_l = -1*theta_spe_new[[j]][[curr_ic]]$lam_l}
    theta_spe_new[[j]][[curr_ic]]$J_l = theta_ic$J_l
  }
  }
}

  #############Y specific component approximately by eigen decomposition
  # Update A,B
  ## Ensemble S for all modality
S = list()
S_sparse = list()
for (j in 1:length(q)){
 S[[j]] = array(dim=c(q_common+q[j],K))
 S_sparse[[j]] = array(dim=c(q_common+q[j],K))
  for (l in 1:nrow(S[[j]]))
  {
    if (l<=q_common){
    S[[j]][l,] = Ltrans(t(theta_new[[j]][[l]]$J_l)%*% diag(theta_new[[j]][[l]]$lam_l) %*%theta_new[[j]][[l]]$J_l,F) } # K x M
    else{
    S[[j]][l,] = Ltrans(t(theta_spe_new[[j]][[l-q_common]]$J_l)%*% diag(theta_spe_new[[j]][[l-q_common]]$lam_l) %*%theta_spe_new[[j]][[l-q_common]]$J_l,F)  # K x M
    }
  }
 S_sparse[[j]] = S[[j]]*(abs(S[[j]])>=3.5*lambda_ch)
}
  ## estimate A and B
  A_new = list()
  Dmat = list()
  for (j in 1:length(q)){
  A_new[[j]] = array(dim = c(q[j]+q_common,q[j]+q_common))
  Dmat[[j]] = Y[[j]]%*%t(Y[[j]])
}


  #For common components
if (Iter){
  for(l in 1:q_common)
  {
    for (j in 1:length(q)){
       complete = (1:length(q))[-j]
       dvec = S[[j]][l,]%*%t(Y[[j]])
       Amat1 = 0

       for (c in complete){
         psi_sign = 1
         if (which(complete==c)!=1){
           if (sum(sign((vec_last*(t(A[[c]][,l])%*%t(H_inv[[c]])%*%H_inv[[j]]))))>0){
             psi_sign=1
           }
           else{psi_sign = -1}}
           dvec = dvec + psi_sign*psi/2*t(A[[c]][,l])%*%t(H_inv[[c]])%*%H_inv[[j]]
           Amat1 = Amat1 +t(A[[c]][,l])%*%t(H_inv[[c]])%*%H_inv[[j]]
           vec_last = t(A[[c]][,l])%*%t(H_inv[[c]])%*%H_inv[[j]]
    }
    if (l ==1){
    Amat2 = c()
    }else if(l==2){Amat2 = A_new[[j]][,1]
    }else{
      Amat2 = t(A_new[[j]][,1:(l-1)])
    }
    A_curr1 = solve.QP(Dmat = Dmat[[j]],dvec = dvec ,Amat = t(rbind(Amat2,Amat1)) ,bvec = rep(0,l),meq = l-1)

    dvec = S[[j]][l,]%*%t(Y[[j]])
    Amat1 = 0
    for (c in complete){
      psi_sign =  1
      if (which(complete==c)!=1){
        if (sum(sign((vec_last*(t(A[[c]][,l])%*%t(H_inv[[c]])%*%H_inv[[j]]))))>0){
          psi_sign= -1
        }
        else{psi_sign= 1}}
        dvec = dvec + psi_sign*psi/2*t(A[[c]][,l])%*%t(H_inv[[c]])%*%H_inv[[j]]
        Amat1 = Amat1 +t(A[[c]][,l])%*%t(H_inv[[c]])%*%H_inv[[j]]
 }
    if (l ==1){
      Amat2 = c()
    }else if(l==2){Amat2 = A_new[[j]][,1]
    }else{
      Amat2 = t(A_new[[j]][,1:(l-1)])
    }
    A_curr2 = solve.QP(Dmat = Dmat[[j]],dvec = dvec ,Amat = t(rbind(Amat2,-Amat1)) ,bvec = rep(0,l),meq = l-1)
    if (A_curr1$value < A_curr2$value){
      A_curr = A_curr1$solution
    }else{
      A_curr = A_curr2$solution
    }

    ai = sqrt(sum(A_curr^2))
    theta_new[[j]][[l]]$lam_l = theta_new[[j]][[l]]$lam_l * ai
    A_new[[j]][,l] =A_curr / ai
    S[[j]][l,] = S[[j]][l,] * ai
    S_sparse[[j]][l,] =  S_sparse[[j]][l,]*ai
  }
  }}else{
  for(l in 1:q_common)
  {
    for (j in 1:length(q)){
      complete = (1:length(q))[-j]
      dvec = S[[j]][l,]%*%t(Y[[j]])
      Amat1 = 0
      for (c in complete){
        psi_sign =  1
        if (t(A[[c]][,l])%*%t(H_inv[[c]])%*%H_inv[[j]] %*% A[[j]][,l]>0){
          psi_sign= -1
        }
        else{psi_sign= 1}
        dvec = dvec + psi_sign*psi/2*t(A[[c]][,l])%*%t(H_inv[[c]])%*%H_inv[[j]]
        Amat1 = Amat1 +t(A[[c]][,l])%*%t(H_inv[[c]])%*%H_inv[[j]]
      }

      if (l ==1){
        Amat2 = c()
      }else if(l==2){Amat2 = A_new[[j]][,1]
      }else{
        Amat2 = t(A_new[[j]][,1:(l-1)])
      }
      A_curr2 = solve.QP(Dmat = Dmat[[j]],dvec = dvec ,Amat = t(rbind(Amat2,-Amat1)) ,bvec = rep(0,l),meq = l-1)
      # if (A_curr1$value < A_curr2$value){
      #   A_curr = A_curr1$solution
      # }else{
      #   A_curr = A_curr2$solution
      # }
      A_curr = A_curr2$solution
      ai = sqrt(sum(A_curr^2))
      theta_new[[j]][[l]]$lam_l = theta_new[[j]][[l]]$lam_l * ai
      A_new[[j]][,l] =A_curr / ai
      S[[j]][l,] = S[[j]][l,] * ai
      S_sparse[[j]][l,] =  S_sparse[[j]][l,]*ai
    }}}

  #For specific components
  for (j in 1:length(q)){
  if (q[j]!=0){
  P_new = diag(1,nrow = q_common+q[j]) - (A_new[[j]][,1:q_common])%*%solve(t(A_new[[j]][,1:q_common])%*%(A_new[[j]][,1:q_common]))%*%t(A_new[[j]][,1:q_common])
  if (q[j] == 1){
  A_specific = P_new %*% Y[[j]]%*%(S[[j]][(q_common+1):nrow(S[[j]]),]) %*% solve(t(S[[j]][(q_common+1):nrow(S[[j]]),])%*%(S[[j]][(q_common+1):nrow(S[[j]]),]))
  }
  else{A_specific = P_new  %*% Y[[j]] %*%t(S[[j]][(q_common+1):nrow(S[[j]]),]) %*% solve((S[[j]][(q_common+1):nrow(S[[j]]),])%*%t(S[[j]][(q_common+1):nrow(S[[j]]),])) }
  if (q[j] == 1){
    norm = sqrt(sum(A_specific^2))
    A_new[[j]][,(q_common+1):nrow(A_new[[j]])] = far::orthonormalization(A_specific,basis = F)
  }else{
  norm = sqrt(apply(A_specific^2,2,sum))
  A_new[[j]][,(q_common+1):nrow(A_new[[j]])] = far::orthonormalization(A_specific,basis = F)}
  S[[j]][(q_common+1):nrow(S[[j]]),] = S[[j]][(q_common+1):nrow(S[[j]]),]*norm
  S_sparse[[j]][(q_common+1):nrow(S[[j]]),] = S_sparse[[j]][(q_common+1):nrow(S[[j]]),]*norm
  }}

  # Save m_l, X_l into theta2_new:
  for (j in 1:length(q)){
  for(l in 1:q_common)
  {
    theta_new[[j]][[l]]$M_l = t(A_new[[j]])[l,]
  }
  if (q[j] != 0){
  for(l in 1:q[j])
  {
    theta_spe_new[[j]][[l]]$M_l = t(A_new[[j]])[l+q_common,]
  }}}
  return(list(A = A_new, S=S ,S_sparse=S_sparse,theta_spe = theta_spe_new, theta_common = theta_new))
  }
