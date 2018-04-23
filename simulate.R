

#' 
#' See: 
#' \link{https://docs.google.com/viewer?a=v&pid=sites&srcid=ZGVmYXVsdGRvbWFpbnxzdGV2ZXRoZWJheWVzaWFufGd4OjI2ZGEwMTk4M2VmOWRkOTE}
#' 

require(bsts)
require(ggplot2)
require(dplyr)
require(magrittr)


###################################################################################################
# Function
###################################################################################################

theme_black = function(base_size = 12, base_family = "") {
  
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    
    theme(
      # Specify axis options
      axis.line = element_blank(),  
      axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.ticks = element_line(color = "white", size  =  0.2),  
      axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),  
      legend.key = element_rect(color = "white",  fill = "black"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*0.8, color = "white"),  
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "black", color  =  NA),  
      panel.border = element_rect(fill = NA, color = "white"),  
      panel.grid.major = element_line(color = "grey35"),  
      panel.grid.minor = element_line(color = "grey20"),  
      panel.margin = unit(0.5, "lines"),   
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = base_size*0.8, color = "white"),  
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),  
      plot.title = element_text(size = base_size*1.2, color = "white"),  
      plot.margin = unit(rep(1, 4), "lines")
      
    )
  
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}


###################################################################################################
# Parameterized Simulate
###################################################################################################

#' @section AR vs Random Walk
#' Here we compare AR1 vs RW and their variances


# Simulate an AR1 and RW time series
sample_size      <- 1000
number_of_series <- 1000
ar1_param        <- 0.95
sd_param         <- 1.00

many_ar1         <- matrix(nrow = sample_size, ncol = number_of_series)
for (i in 1:number_of_series) {
  many_ar1[, i] <- arima.sim(model = list(ar = ar1_param, sd = sd_param), n = sample_size)
}
many_rw <- matrix(nrow = sample_size, ncol = number_of_series)
for (i in 1:number_of_series) {
  many_rw[, i] <- cumsum(rnorm(sample_size, sd = sd_param))
}
par(mfrow = c(1, 2))
plot.ts(many_ar1, plot.type = "single")
plot.ts(many_rw, plot.type = "single")
par(mfrow = c(1, 1))


#' Variance
#' @AR1 is Var(y) = var(e)/(1 - abs(ar1_param))
#' if abs(ar1_param) < 1
#' @RW is Var(y) = tvar(e)
var_ar1 <- data.frame(
  t = seq_len(nrow(many_ar1)),
  y = apply(many_ar1, 1, var))
var_ar1_actual <- (sd_param)^2/(1 - ar1_param^2)
ggplot(var_ar1, aes(x = t, y = y)) + 
  geom_line() +
  geom_hline(yintercept = var_ar1_actual, color = "red") +
  annotate("text", x = 0, y = var_ar1_actual * 1.05, label = round(var_ar1_actual, 2), color = "red", size = 7) +
  ggtitle(paste0("Var(y) over t - Simulated vs Actual - AR1 [ϕ = ", ar1_param, "]"))

var_rw <- data.frame(
  t = seq_len(nrow(many_rw)),
  y = apply(many_rw, 1, var))
var_rw_actual <- (sd_param)^2 * seq_len(nrow(many_rw))
ggplot(var_rw, aes(x = t, y = y)) + 
  geom_line() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  ggtitle(paste0("Var(y) over t - Simulated vs Actual - RW [sd(e) = ", sd_param, "]"))






###################################################################################################
# Covariate Simulate
###################################################################################################


## Example 7: Regressors with multiple time stamps.
number.of.time.points <- 50
sample.size.per.time.point <- 10
total.sample.size <- number.of.time.points * sample.size.per.time.point
sigma.level <- .1
sigma.obs <- 1

## Simulate some fake data with a local level state component.
trend <- cumsum(rnorm(number.of.time.points, 0, sigma.level))
plot(trend, type = 'o')

predictors <- matrix(rnorm(total.sample.size * 2), ncol = 2)
colnames(predictors) <- c("X1", "X2")
coefficients <- c(-10, 10)
regression <- predictors
y.hat <- matrix(rep(rep(trend, sample.size.per.time.point), 2), byrow = F, ncol = 2) + regression
y <- rnorm(nrow(y.hat), y.hat, sigma.obs)

## Create some time stamps, with multiple observations per time stamp.
first <- as.POSIXct("2013-03-24")
dates <- seq(from = first, length = number.of.time.points, by = "month")
timestamps <- rep(dates, sample.size.per.time.point)
## Run the model with a local level trend, and an unnecessary seasonal component.
ss <- AddLocalLevel(list(), y)
ss <- AddSeasonal(ss, y, nseasons = 7)
model <- bsts(y ~ predictors, ss, niter = 1000, timestamps = timestamps,
              seed = 8675309)
plot(model)
plot(model, "components")

model$coefficients

(coef_est <- apply(model$coefficients, 2, mean))
(coef_sd  <- apply(model$coefficients, 2, sd))




###################################################################################################
# Custom Simulate
###################################################################################################

# Set seed for reproduction
set.seed(61531)

### Parameters
n_time_points <- 52 * 2
m_samples_per_time <- 7
N_total <- n_time_points * m_samples_per_time
sd_y <- 3
# Trend
sd_trend <- 2
# Regression
coef_x1 <- 0.70
coef_x2 <- 5.0
coef_x4_x2 <- 0.5  # Correlation with another variable
# Season
coef_season_scale <- 50
sd_season <- 5
n_seasons <- m_samples_per_time
f_season <- function(s) {
  s_vec <- seq_len(s)/s
  s_out <- (s_vec - 0.5) ^ 2
  s_out - sum(s_out)/s
}
coef_season <- f_season(n_seasons)


### Simulation
# Random Trend
trend <- cumsum(rnorm(N_total + 7 - 1, 0, sd_trend))
trend <- rollmean(trend, k = 7)
plot(trend, type = 'o')

# Make a predictor that is a constant
# x_1_raw <- cumsum(rnorm(N_total, mean = 0, sd = 4))
x_1_raw <- arima.sim(model = list(ar = 0.75, sd = 4), n = N_total)
x_2_raw <- arima.sim(model = list(ar = 0.50, sd = 2), n = N_total)
x_3_fake <- arima.sim(model = list(ar = 0.90, sd = 10), n = N_total)
x_4_raw <- x_2_raw * coef_x4_x2 + rnorm(N_total, sd = 1)
plot(x_1_raw, type = 'o')

# Simulation
s_season <- rep(coef_season_scale * coef_season, N_total/n_seasons) +
  rnorm(N_total, sd = sd_season)

# Combine the effect of trend, regression, and season
y_hat <- trend + (x_1_raw * coef_x1) + (x_2_raw * coef_x2) + s_season
  
# Make actuals with systematic noise
y <- rnorm(length(y_hat), y_hat, sd_y)
plot(y, type = 'o')

# Create some time stamps, with multiple observations per time stamp.
first <- as.POSIXct("2016-01-01")
# dates <- seq(from = first, length = n_time_points, by = "days")
# timestamps <- rep(dates, m_samples_per_time)
timestamps <- seq(from = first, length = N_total, by = "days")

# Combine into a single dataframe
sim_df_1 <- data.frame(
  t = timestamps,
  y = y,
  y_hat = y_hat,
  x_1 = x_1_raw,
  x_2 = x_2_raw,
  x_3 = x_3_fake,
  x_4 = x_4_raw,
  trend = trend,
  s_season = s_season,
  stringsAsFactors = F
)
sim_df_1$reg <- sim_df_1$y - sim_df_1$trend
ggplot(sim_df_1, aes(x = t, y = y)) +
  geom_line(color = 'black', size = 2) +
  geom_line(aes(y = y_hat), color = 'green') +
  geom_line(aes(y = trend), color = 'blue') +
  geom_line(aes(y = x_1), color = 'red') +
  geom_line(aes(y = s_season), color = 'grey') +
  theme_bw()
pairs(sim_df_1, lower.panel = panel.cor)

# Run the model with a local level trend, and an unnecessary seasonal component.
ss <- AddLocalLevel(list(), y)
ss <- AddSeasonal(ss, y, nseasons = 7)
model <- bsts(y ~ x_1, 
              data = sim_df_1, ss, niter = 1000, timestamps = timestamps,
              seed = 8675309)
model_2 <- bsts(y ~ x_1 + x_2, 
                data = sim_df_1, ss, niter = 1000, timestamps = timestamps,
                seed = 8675309)
model_3 <- bsts(y ~ x_1 + x_2 + x_3, 
                data = sim_df_1, ss, niter = 1000, timestamps = timestamps,
                seed = 8675309)
model_4 <- bsts(y ~ x_1 + x_2 + x_3 + x_4, 
                data = sim_df_1, ss, niter = 1000, timestamps = timestamps,
                seed = 8675309)

# View Model
plot(model_3, "state")
plot(model_3, "components", same.scale = F)
plot(model, "residuals")
plot(model_4, "coefficients")
plot(model, "prediction.errors")
plot(model, "forecast.distribution")
plot(model, "predictors")
plot(model, "size")
plot(model, "seasonal")
plot(model, "help")
# plot(model, "dynamic")
CompareBstsModels(list("x1" = model,
                       "x1 + x2" = model_2,
                       "x1 + x2 + x3" = model_3,
                       "x1 + x2 + x3 + x4" = model_4))

#' @param sigma.obs is the systemmatic error given every effect. In our simulation it is 'sd_y'
#' @param sigma.level is the systemmatic error for the trend. In our simulation it is 'sd_trend'
#' @param sigma.seasonal is the systemmatic error for the seasonal effect. In our simulation it is 'sd_season'
#' @param coefficients are the regressors and season
mean(model$sigma.obs)
sd_y
mean(model$sigma.level)
sd_trend
mean(model$sigma.seasonal)
sd_season


names(model)
model$coefficients
model$state.contributions
test <- sapply(seq_len(dim(model$state.contributions)[1]), function(x) {
  apply(model$state.contributions[x,,], 1, mean)
})
mean_contrib <- t(test)

lapply(model$state.specification, function(x) x$name)

for(c in seq_along(names(model))) {
  print(paste("==========================", c))
  print(names(model)[c])
  print(head(model[[c]]))
}

### Compare Estimated vs Actual
# Get coefficients
model = model_4
coef_df <- data.frame(
  coef  = colnames(model$coefficients),
  est   = apply(model$coefficients, 2, mean),
  se    = apply(model$coefficients, 2, sd)/sqrt(model$niter),
  sd    = apply(model$coefficients, 2, sd),
  lb_95 = qnorm(0.025) * apply(model$coefficients, 2, sd),
  ub_95 = qnorm(0.925) * apply(model$coefficients, 2, sd),
  stringsAsFactors = F
)
coef_df
plot(model$log.likelihood)
model$state.contributions

qqnorm(model$coefficients[,2][model$coefficients[,2] != 0])
qqline(model$coefficients[,2][model$coefficients[,2] != 0])








