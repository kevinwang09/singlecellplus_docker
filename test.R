start = Sys.time()
source("qc.R")
source("scMerge.R")
source("downstream.R")
end = Sys.time()

print(start - end)