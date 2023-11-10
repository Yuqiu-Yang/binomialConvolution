rm(list=ls())
library(binomialConvolution)
library(ggplot2)
library(gridExtra)
count = read.csv("./data/count_data.csv")
passage = read.csv("./data/passage_info.csv")

easy_cols = paste("N", 1:4, sep="")
hard_cols = paste("N", 5:13, sep="")

p = list()
counter = 1
for(passage_id in passage$Passage_ID)
{
  df = count[which(count$Passage_ID == passage_id),]
  n_trials = c(sum(passage[which(passage$Passage_ID == passage_id), easy_cols]),
               sum(passage[which(passage$Passage_ID == passage_id), hard_cols]))

  samples = df$Count

  p_bar = mean(samples)/sum(n_trials)


  result = optim(par=c(p_bar+0.001, p_bar-0.001),
                 fn=log_likelihood,
                 n_trials = n_trials,
                 samples=samples,
                 computation_method="convolution",
                 lower = 0.001,
                 upper=0.999,
                 method="L-BFGS-B",
                 control=list(fnscale=-1),
                 hessian=T)

  b_pmf = binomial_convolution_pmf(n_trials=n_trials,
                                   success_probs = result$par)
  y=0:sum(n_trials)
  # plot(y, b_pmf(y), type="l")

  sim_samples = simulate_binomial_convolution(nrow(df),
                                              n_trials = n_trials,
                                              success_probs = result$par)
  sample_count = data.frame(table(samples))
  sample_count$label = "true"
  sim_sample_count = data.frame(table(sim_samples))
  sim_sample_count$label = "simulated"
  colnames(sample_count) = colnames(sim_sample_count) = c("y", "count", "label")
  compare_df = rbind(sample_count,
                     sim_sample_count)
  compare_df$y = as.numeric(as.character(compare_df$y))
  # compare_df = data.frame(y = c(samples, sim_samples),
  #                         label = rep(c("true", "simulated"),
  #                                     times=c(nrow(df), nrow(df))))
  est_df = data.frame(x=y,
                      p = b_pmf(y)*nrow(df))
  p[[counter]] = ggplot() +
    geom_point(data=compare_df,
                   aes(x=y,
                       y=count,
                       color = label),
               size=10) +
    geom_point(data=est_df, aes(x=x, y=p),
               size=5)+
    scale_x_continuous(labels = as.character(y), breaks = y)+
    coord_cartesian(xlim = c(min(samples)-1, max(samples)+1))+
    ggtitle(paste0(passage_id,"_",
                   round(result$par[1], 4),
                   "_",
                   round(result$par[2], 4)))

  counter = counter + 1
}

pdf("./data/results.pdf", onefile = TRUE)
for (i in seq(length(p))) {
  print(p[[i]])
}
dev.off()
