rm(list = ls())
library(testthat)
library(binomialConvolution)
library(plotly)
library(numDeriv)

M = sum(n_trials)

sum(n_trials * success_probs * (1-success_probs) * (1-2*success_probs))/(sum(n_trials * success_probs * (1-success_probs))^1.5)

p_p = 2*c(p_bar, p_bar) - success_probs
sum(n_trials * p_p * (1-p_p) * (1-2*p_p))/(sum(n_trials * p_p * (1-p_p)))^1.5


n_trials = c(55, 21)
success_probs = c(0.7,0.3)
mu = sum(n_trials * success_probs)
p_bar = mu/76
p1=seq(-1,1,by=0.1)
p2=mu/20-p1*50/20
w=2*pi/71
(1-(p1+p_bar) + (p1+p_bar)*exp(1i * w * 2))^50 * (1-(p2+p_bar) + (p2+p_bar)*exp(1i * w * 2))^20

# test_check("binomialConvolution")
# 2D
n_trials = c(50, 20)
success_probs = c(0.8,0.6)
samples = simulate_binomial_convolution(500, n_trials, success_probs)
x = seq(0, pi, length.out = 100)
ecf = empirical_cf(samples)
p = c(0.7,0.4)
cf = theoretical_cf(n_trials,
                    p)
plot(x, Mod(ecf(x)), type="l", ylim=c(0,1))
lines(x, Mod(cf(x)), col=3)

foo1 <- function(s,
                 cf, ecf)
{
  return(abs(Mod(cf(s)) - Mod(ecf(s))))
}

foo2 <- function(p,
                 n_trials,
                 ecf)
{
  cf = theoretical_cf(n_trials,
                      p)
  return(optim(pi/2, foo1, cf=cf, ecf=ecf,
               control = list(fnscale=-1),
               method="L-BFGS-B",
               lower=0,
               upper=pi)$value)
}

optim(c(0.7,0.5),
      foo2,
      n_trials=n_trials,
      ecf=ecf,
      method="L-BFGS-B",
      lower=0.01,
      upper=0.99)


result_df = expand.grid(p1 = seq(0.1, 0.9, by=0.05),
                        p2 = seq(0.1, 0.9, by=0.05))
result_df$p1_hat = result_df$p2_hat = 0
for(i in 1:nrow(result_df))
{
  p1 = result_df$p1[i]
  p2 = result_df$p2[i]
  result = optim(par=c(p1,p2),
                 fn=log_likelihood,
                 n_trials = n_trials,
                 samples=samples,
                 computation_method="convolution",
                 lower = 0.05,
                 upper=0.95,
                 method="L-BFGS-B",
                 control=list(trace=3,
                              fnscale=-1),
                 hessian=T)
  result_df$p1_hat[i] = result$par[1]
  result_df$p2_hat[i] = result$par[2]
}

ggplot(result_df, aes(x = p1, y = p2)) +
  geom_segment(aes(xend = p1_hat, yend = p2_hat),
               arrow = arrow(length = unit(0.1, "cm")),
               linewidth = 0.5) +
  geom_abline(intercept = mean(samples)/n_trials[2],
              slope = -n_trials[1]/n_trials[2],
              color="red")


# 3D
n_trials = c(80, 50, 10)
success_probs = c(0.7, 0.4, 0.2)
samples = simulate_binomial_convolution(100, n_trials, success_probs)

result_df = expand.grid(p1 = seq(0.1, 0.9, by=0.05),
                        p2 = seq(0.1, 0.9, by=0.05),
                        p3 = seq(0.1, 0.9, by=0.05))
result_df$p1_hat = result_df$p2_hat = result_df$p3_hat = 0
result_df$lk = 0
for(i in 1:nrow(result_df))
{
  p1 = result_df$p1[i]
  p2 = result_df$p2[i]
  p3 = result_df$p3[i]
  result = optim(par=c(p1, p2, p3),
                 fn=log_likelihood,
                 n_trials = n_trials,
                 samples=samples,
                 computation_method="convolution",
                 lower = 0.05,
                 upper=0.95,
                 method="L-BFGS-B",
                 control=list(fnscale=-1),
                 hessian=T)
  result_df$p1_hat[i] = result$par[1]
  result_df$p2_hat[i] = result$par[2]
  result_df$p3_hat[i] = result$par[3]
  result_df$lk[i]= result$value
}

result_df$d_p1 = result_df$p1_hat - result_df$p1
result_df$d_p2 = result_df$p2_hat - result_df$p2
result_df$d_p3 = result_df$p3_hat - result_df$p3



fig = plot_ly(data=result_df,
              x=result_df$p1,
              y=result_df$p2,
              z=result_df$p3,
              u=result_df$p1_hat,
              v=result_df$p2_hat,
              w=result_df$p3_hat,
              type="cone")
fig

plane = expand.grid(p1=seq(0,1,by=0.01),
                    p2=seq(0,1,by=0.01))
plane$p3 = (mean(samples) -
          n_trials[1] * plane$p1 -
          n_trials[2] * plane$p2)/(n_trials[3])
plane = plane[(plane$p3 >=0) &
                (plane$p3 <=1) ,]


plot_ly(x = plane$p1,
        y = plane$p2,
        z = plane$p3,
        type='mesh3d')%>%
  add_markers(data=result_df, ~p1_hat, ~p2_hat, ~p3_hat,
              size=1)







# Check DFT is correct
rm(list = ls())
n_trials = c(1, 2)
success_probs = c(0.1,0.9)
pmf_conv = binomial_convolution_pmf(n_trials = n_trials,
                                    success_probs = success_probs,
                                    computation_method = "convolution")
pmf_dft = exact_pmf_dft(n_trials = n_trials,
                                   success_probs = success_probs)
y = seq(0, sum(n_trials))
plot(y, pmf_conv(y), type='l', col=2, lwd = 5)
# plot(y, pmf_dft(y), type='l', col=3, lwd = 5)
lines(y, pmf_dft(y), col= 3, lwd =2)

samples = simulate_binomial_convolution(100, n_trials, success_probs)

p_seq = seq(0.01, 0.99, length.out = 50)
result_score = data.frame(matrix(0, ncol=4, nrow=2500))
colnames(result_score) =c("x","y","vx","vy")
counter = 1
for(p1 in 1 : 50)
{
  for(p2 in 1:50)
  {
    s_prob = c(p_seq[p1], p_seq[p2])
    v = score(success_probs = s_prob,
          n_trials = n_trials,
          samples=samples)
    result_score[counter, 1:2] = s_prob
    result_score[counter, 3:4] = v
    counter = counter + 1
  }
}


ggplot(result_score, aes(x = x, y = y)) +
  geom_segment(aes(xend = x+vx/50000, yend = y+vy/50000),
               arrow = arrow(length = unit(0.1, "cm")),
               linewidth = 0.5) +
  geom_point(x=0.8, y=0.4, color="red") +
  geom_point(x=0.65, y=0.8, color="blue") +
  geom_point(x=0.9, y=0.1)


samples = simulate_binomial_convolution(10000,
                                        n_trials,
                                        c(0.8, 0.4))

samples1 = simulate_binomial_convolution(10000,
                                         n_trials,
                                         c(0.65, 0.8))

hist(samples)
hist(samples1, col=2, add=T)




result_score1 = data.frame(matrix(0, ncol=4, nrow=2500))
colnames(result_score1) =c("x","y","vx","vy")
counter = 1
for(p1 in 1 : 50)
{
  for(p2 in 1:50)
  {
    s_prob = c(p_seq[p1], p_seq[p2])
    v = grad(log_likelihood,
             x=s_prob,
             n_trials = n_trials,
             samples=samples,
             computation_method="convolution")
    result_score1[counter, 1:2] = s_prob
    result_score1[counter, 3:4] = v
    counter = counter + 1
  }
}

ggplot(result_score1, aes(x = x, y = y)) +
  geom_segment(aes(xend = x+vx/950000, yend = y+vy/950000),
               arrow = arrow(length = unit(0.1, "cm")),
               size = 0.5)





score(success_probs = success_probs,
      n_trials = n_trials,
      samples=samples)
grad(log_likelihood,
     x=success_probs,
     n_trials = n_trials,
     samples=samples)
hessian(success_probs = success_probs,
        n_trials = n_trials,
        samples=samples)

result = matrix(0, nrow = 20,
                ncol=20)
p_seq = seq(0.1, 0.9, length.out = 20)
for(p1 in 1:20)
{
  for(p2 in 1:20)
  {
    s_prob = c(p_seq[p1], p_seq[p2])
    h = hessian(success_probs = s_prob,
                n_trials = n_trials,
                samples=samples)
    # h = numDeriv::hessian(log_likelihood,
    #                       x=s_prob,
    #                       n_trials=n_trials,
    #                       samples=samples)
    result[p1, p2] = all(eigen(h)$values<0)
  }
}

heatmap.2(result,dendrogram='none',
          Rowv=FALSE,
          Colv=FALSE,trace='none')







true_N_trials = c(1, 1)
true_ps = c(1, 0.7)

Y = simulate_binomial_convolution(n_samples=500,
                                  n_trials=true_N_trials,
                                  success_probs=true_ps)

mu = mean(Y)
sig = var(Y)

M1 = true_N_trials[1]
M2 = true_N_trials[2]

p2 = (-2*mu*M2 -sqrt((2*mu*M2)^2 + 4 * (M2^2 + M1 * M2) * (M1*mu - mu^2-sig*M1)))/(-2*(M2^2 + M1 * M2))
p1 = (mu-M2*p2)/M1

c(p1, p2)


result = optim(par=c(0.9,0.2),
      fn=log_likelihood,
      n_trials = n_trials,
      samples=samples,
      computation_method="convolution",
      lower = 0.1,
      upper=0.8,
      method="L-BFGS-B",
      control=list(trace=3,
                   fnscale=-1),
      hessian=T)




Y1 = simulate_binomial_convolution(n_samples=500,
                                   n_trials=true_N_trials,
                                   success_probs = c(0.3-0.2, 0.7+0.1))


dis_m = matrix(0, nrow = 50,
               ncol=50)
x = seq(0.01,0.99, length.out=50)
for(i in 1:50)
{
  for(j in 1 : 50)
  {
    # dis_m[i,j] = objective_function_1(success_probs=c(x[i], x[j]),
    #                                 n_trials=true_N_trials,
    #                                 samples=Y,
    #                                 weight_function=function(s){return(exp(-60*abs(s)))},
    #                                 lower_limit=-1,
    #                                 upper_limit=1)
    dis_m[i,j] = log_likelihood(success_probs=c(x[i], x[j]),
                                n_trials=true_N_trials,
                                samples=Y,
                                computation_method="convolution")
  }
}


library(plotly)
plot_ly(x=x,
        y=x,
        z=dis_m) %>%
  add_surface(contours=list(z=list(show=TRUE,
                                   start=-10,
                                   end=0,
                                   size=0.2,
                                   color="red"))) %>%
  add_trace(x=0.7, y=0.3,
            z=log_likelihood(success_probs=c(0.3, 0.7),
                             n_trials=true_N_trials,
                             samples=Y),
            mode = "markers", type = "scatter3d",
            marker = list(size = 20, color = "cyan", symbol = 104))


plot_ly(x=x,
        y=x,
        z=dis_m) %>%
  add_surface(contours=list(z=list(show=TRUE,
                                   start=0,
                                   end=0.04,
                                   size=0.001,
                                   color="red"))) %>%
  add_trace(x=0.7, y=0.3,
            z=objective_function_1(success_probs=c(0.3, 0.7),
                                   n_trials=true_N_trials,
                                   samples=Y,
                                   weight_function=function(s){return(exp(-60*abs(s)))},
                                   lower_limit=-1,
                                   upper_limit=1),
            mode = "markers", type = "scatter3d",
            marker = list(size = 20, color = "cyan", symbol = 104))





est_ps = c(0.7,0.8)

optim(par=true_ps,
      fn=objective_function,
      n_trials = true_N_trials,
      Y=Y,
      lower = 0.01,
      upper=0.99,
      method="L-BFGS-B",
      control=list(trace=2))


ecgf = empirical_cgf(samples = Y)

est_ps = c(0.1, 0.2, 0.3)

cgf = theoretical_cgf(n_trials = true_N_trials,
                      success_probs = est_ps)

integrand = discrepency_integrand(theoretical_cgf = cgf,
                                  empirical_cgf = ecgf,
                                  weight_function = function(s){return(exp(-s^2))})

temp = c()
for(x in seq(-1,1,by=0.01))
{
  temp = c(temp, integrand(x))
}
plot(seq(-1,1,by=0.01),temp, pch=19)

x = seq(-10,10,by=0.01)
plot(x,
     integrand(x),
     type="l")

integrate(integrand,
          lower = -10,
          upper = 10)

