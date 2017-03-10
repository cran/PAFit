if (FALSE) {
  #setwd("tests")
  rm(list = ls())
  library(PAFit)
  set.seed(1)
  prob_m <- "FALSE"
  inc    <- "FALSE"
  log    <-  c("FALSE")
  
  M <- 50
  
  alpha_vec     <- rep(0,M)
  s_vec         <- rep(0,M)
  alpha_optimal <- rep(0,M)
  r_optimal     <- rep(0,M)
  for (i in 1:M) {
  
      net  <- GenerateNet(N = 1000, m = 50,prob_m = prob_m, num_seed = 100, multiple_node = 100,
                          increase = inc, log = log, 
                          mode = 1, shape = 1, rate = 1,alpha = 0.25)
  
      net_stats <- GetStatistics(net$graph,deg_threshold = 0, 
                                 net_type = "directed",
                                 Binning = TRUE, G = 50) 
  
      print(system.time(result <- JointEstimate(raw_net = net$graph, net_stat = net_stats, print.out = TRUE,
                                              cv_deg_thresh = c(1,10))))
      
      alpha_vec[i]     <- result$estimate_result$alpha
      s_vec[i]         <- result$estimate_result$shape
      alpha_optimal[i] <- result$cv_result$alpha_optimal
      r_optimal[i]     <- result$cv_result$r_optimal
  }
  print(alpha_vec)  
  print(s_vec)
  mean(s_vec)
}