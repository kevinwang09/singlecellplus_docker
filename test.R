start = Sys.time()
source("/home/studio/qc.R")
source("/home/studio/scMerge.R")
source("/home/studio/downstream.R")
end = Sys.time()

print(start - end)