using Revise
using BenchmarkTools
using UnicodePlots

using jLRS


# How do we want to use this tool?
##################################

# Instantiate a sampling strategy

# Example 1:
# Look at some location and how OCO-x samples in some radius of FoCo
# ---------

s1 = jLRS.OCOSampling(-105.0844, # target lon
                      40.5853, # target lat
                      50, # radius in km
                      "/home/psomkuti/oco3.db")

s2 = jLRS.OCOSampling(-104.0844, # target lon
                      40.5853, # target lat
                      50, # radius in km
                      "/home/psomkuti/oco2.db")

s3 = s1 + s2

# Using the locations and samplings, produce a light response curve given some sampling
weekly = jLRS.WeeklySampling(1)
halfweekly = jLRS.WeeklySampling(0.5)

noagg = jLRS.FullROI()
bigagg = jLRS.RegularGridCells(2, 2)

# Perform the aggregation step
agg1 = jLRS.aggregate_scenes(s3, noagg, weekly)
#agg2 = jLRS.aggregate_scenes(s1, noagg, halfweekly)
#agg3 = jLRS.aggregate_scenes(s1, bigagg, weekly)
#agg4 = jLRS.aggregate_scenes(s1, bigagg, halfweekly)

# Example 2:
# Look at some bounding box (lower left, upper right) with OCO-x samples
# ---------
# s2 = OCOSampling(bounding_box=[-55.0, -14.0, -47.0, -5.0],
#                  "/path/to/OCO3")
