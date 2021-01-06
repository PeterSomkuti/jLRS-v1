using Revise
using jLRS


# How do we want to use this tool?
##################################

# Instantiate a sampling strategy

# Example 1:
# Look at some location and how OCO-x samples in 100km radius of FoCo
# ---------
s1 = jLRS.OCOSampling(40.5853, # target lon
                      105.0844, # target lat
                      100, # radius in km
                      "/data6/OCO3/SCF/TEST_B10215_STO_SRU_MEAS201028_DYn0112_DZ0116/2020/10/10/L1bSc")

# Example 2:
# Look at some bounding box (lower left, upper right) with OCO-x samples
# ---------
# s2 = OCOSampling(bounding_box=[-55.0, -14.0, -47.0, -5.0],
#                  "/path/to/OCO3")
