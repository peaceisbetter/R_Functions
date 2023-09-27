# Function for coarse graining data
# inputs: 
#   data: a single time series
#   window_size: the number you want to divide your data by. 
#     e.g. data length = 12,000, window_size = 4, output data length = 3,000
coarse_grain <- function(data, window_size)
{
  return(ts(as.ts(rollapply(zoo(data), width=window_size, by=window_size, FUN=mean)), frequency=1))
}