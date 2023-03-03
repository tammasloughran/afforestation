#!/usr/bin/env python3
# Script supplied by Spencer Liddicoat to process co2mass data from UKESM1-0-LL.

import iris
import iris.coord_categorisation
import numpy as np
import matplotlib.pyplot as plt


def yr(x):
    """Return year from annual fields"""
    yr2 = x.coord("time")
    y2 = yr2.units.num2date(yr2.points)
    y3 = np.array(y2[0].year)
    for ii in range(1, len(y2)):
        y3 = np.append(y3, y2[ii].year)
    return y3


def annual(x):
    """Return annual mean from monthly fields"""
    if not x.coords("year"):
        iris.coord_categorisation.add_year(x, "time", name="year")
    x_ann = x.aggregated_by("year", iris.analysis.MEAN)
    return x_ann


def concat_cubes(files):
    """Return the concatenated monthly cubes from the two files"""
    cubes = iris.load(files)
    for cube in cubes:
        # Delete the problem attribute from the time coordinate: #    del cube.coord('time').attributes['creation_date']
        cube.attributes["history"] = ""
        cube.attributes["creation_date"] = ""
        cube.attributes["tracking_id"] = ""
    return cubes.concatenate_cube()


# Unit conversion
kgCO2_to_GtC = 12.01 / (44.01 * 1e12)
GtC_to_ppm = 1.0 / 2.124  # Ballantyne et al 2014 (?)
kgCO2_to_ppm = kgCO2_to_GtC * GtC_to_ppm

# Set up files
files1 = [
    'co2mass_Amon_UKESM1-0-LL_esm-ssp585-ssp126Lu_r1i1p1f2_gm_201501-204912.nc',
    'co2mass_Amon_UKESM1-0-LL_esm-ssp585-ssp126Lu_r1i1p1f2_gm_205001-210012.nc',
    ]

cubes1 = iris.load(
    [
        'co2mass_Amon_UKESM1-0-LL_esm-ssp585-ssp126Lu_r1i1p1f2_gm_201501-204912.nc',
        'co2mass_Amon_UKESM1-0-LL_esm-ssp585-ssp126Lu_r1i1p1f2_gm_205001-210012.nc',
    ],
    'atmosphere_mass_of_carbon_dioxide',
    )

files2 = [
    'co2mass_Amon_UKESM1-0-LL_esm-ssp585_r1i1p1f2_gm_201501-204912.nc',
    'co2mass_Amon_UKESM1-0-LL_esm-ssp585_r1i1p1f2_gm_205001-210012.nc',
]

cubes2 = iris.load(
    [
        'co2mass_Amon_UKESM1-0-LL_esm-ssp585_r1i1p1f2_gm_201501-204912.nc',
        'co2mass_Amon_UKESM1-0-LL_esm-ssp585_r1i1p1f2_gm_205001-210012.nc',
    ],
    'atmosphere_mass_of_carbon_dioxide',
    )

co2mass_kg_esm_ssp585_ssp126Lu_r1i1p1f2 = concat_cubes(files1)  # .concatenate_cube()
co2_ppm_esm_ssp585_ssp126Lu_r1i1p1f2 = (
    annual(co2mass_kg_esm_ssp585_ssp126Lu_r1i1p1f2).data*kgCO2_to_ppm
    )
calendar_co2mass_kg_esm_ssp585_ssp126Lu_r1i1p1f2 = yr(annual(concat_cubes(files1)))

print(" ")
print("co2_ppm_esm_ssp585_ssp126Lu_r1i1p1f2 = ", co2_ppm_esm_ssp585_ssp126Lu_r1i1p1f2)
print(" ")
print(
    "calendar_co2mass_kg_esm_ssp585_ssp126Lu_r1i1p1f2 = ",
    calendar_co2mass_kg_esm_ssp585_ssp126Lu_r1i1p1f2,
    )
print(" ")
print(" ")
print(" ")

co2mass_kg_esm_ssp585_r1i1p1f2 = concat_cubes(files2)  # .concatenate_cube()
co2_ppm_esm_ssp585_r1i1p1f2 = annual(co2mass_kg_esm_ssp585_r1i1p1f2).data*kgCO2_to_ppm
calendar_co2mass_kg_esm_ssp585_r1i1p1f2 = yr(annual(concat_cubes(files2)))

print(" ")
print("co2_ppm_esm_ssp585_r1i1p1f2 = ", co2_ppm_esm_ssp585_r1i1p1f2)
print(" ")
print("calendar_co2mass_kg_esm_ssp585_r1i1p1f2 = ", calendar_co2mass_kg_esm_ssp585_r1i1p1f2)
print(" ")
print(" ")
print(" ")

difference_126Lu_minus_standard = (
    co2_ppm_esm_ssp585_ssp126Lu_r1i1p1f2 - co2_ppm_esm_ssp585_r1i1p1f2
    )
print(" ")
print("difference_126Lu_minus_standard = ", difference_126Lu_minus_standard)
print(" ")

calendar1 = calendar_co2mass_kg_esm_ssp585_r1i1p1f2
var1 = difference_126Lu_minus_standard
a4land = (11.69,8.27)
a4port = (8.27,11.69)

var1 = co2_ppm_esm_ssp585_ssp126Lu_r1i1p1f2
var2 = co2_ppm_esm_ssp585_r1i1p1f2

calendar1 = calendar_co2mass_kg_esm_ssp585_ssp126Lu_r1i1p1f2
calendar2 = calendar_co2mass_kg_esm_ssp585_r1i1p1f2

calendar3 = calendar_co2mass_kg_esm_ssp585_r1i1p1f2
var3 = difference_126Lu_minus_standard

fig = plt.gcf()
fig.set_size_inches(a4port)
plt.subplot(211)
plt.xlim([2015, 2100])
plt.title("UKESM1-0-LL: CO2 (ppm) in   esm-ssp585-ssp126Lu   and   esm-ssp585")
plt.ylabel("CO2 (ppm)")
plt.plot(calendar1, var1, "k", label="esm-ssp585-ssp126Lu")
plt.plot(calendar2, var2, "red", label="esm-ssp585")
plt.legend(loc="best", prop={"size": 10})

plt.subplot(212)
plt.ylim([-25, 2])
plt.xlim([2015, 2100])
plt.title("UKESM1-0-LL: Delta CO2 (ppm) [esm-ssp585-ssp126Lu - esm-ssp585]")
plt.ylabel("Delta CO2 (ppm)")
plt.xlabel("Year")
plt.plot(calendar3, var3, "k")
plt.legend(loc="best", prop={"size": 10})
plt.savefig("co2_esm_ssp585_ssp126Lu_and_esm_ssp585_and_difference.png")
plt.show()
plt.close()


fig = plt.gcf()
fig.set_size_inches(a4land)
plt.ylim([-25, 2])
plt.xlim([2015, 2100])
plt.title("UKESM1-0-LL: Delta CO2 (ppm) [esm-ssp585-ssp126Lu - esm-ssp585]")
plt.ylabel("Delta CO2 (ppm)")
plt.xlabel("Year")
plt.plot(calendar3, var3, "k")
plt.legend(loc="best", prop={"size": 7.5})
plt.savefig("difference_126Lu_minus_standard.png")
plt.show()
plt.close()

# [hadsl@vld667 FOR_TAMMAS_LOUGHRAN ]$ python extract_co2mass_vn3.py
#
# co2_ppm_esm_ssp585_ssp126Lu_r1i1p1f2                    =  [ 409.86435  412.86044  415.64072  418.2894   421.2275   424.06015
#  427.32812  431.074    434.88058  437.86197  440.70035  443.98523
#  448.55127  452.0446   455.5525   459.756    463.45178  467.59286
#  471.4693   476.62518  481.93414  485.70074  490.51297  495.59457
#  500.6811   506.1226   511.8401   517.7143   524.1846   530.72314
#  535.94226  541.8053   548.64575  555.2149   562.07184  569.54956
#  577.0409   584.3829   591.7246   599.5155   608.2901   617.00494
#  625.7843   635.3482   644.8851   653.5727   662.6138   672.1841
#  682.93726  693.18463  703.55115  714.43835  725.6508   737.3981
#  748.36426  758.944    770.6142   783.0665   794.81775  806.43
#  818.297    830.9239   843.0899   855.3816   868.58734  882.06555
#  896.7623   910.51666  924.47833  938.6327   951.9604   965.9483
#  981.2608   995.74603 1008.3736  1022.0523  1036.5997  1049.4207
# 1062.49    1076.915   1092.3136  1106.9558  1121.4246  1134.8689
# 1148.6837  1163.917  ]
#
# calendar_co2mass_kg_esm_ssp585_ssp126Lu_r1i1p1f2        =  [2015 2016 2017 2018 2019 2020 2021 2022 2023 2024 2025 2026 2027 2028
# 2029 2030 2031 2032 2033 2034 2035 2036 2037 2038 2039 2040 2041 2042
# 2043 2044 2045 2046 2047 2048 2049 2050 2051 2052 2053 2054 2055 2056
# 2057 2058 2059 2060 2061 2062 2063 2064 2065 2066 2067 2068 2069 2070
# 2071 2072 2073 2074 2075 2076 2077 2078 2079 2080 2081 2082 2083 2084
# 2085 2086 2087 2088 2089 2090 2091 2092 2093 2094 2095 2096 2097 2098
# 2099 2100]
#
#
#
#
# co2_ppm_esm_ssp585_r1i1p1f2                    =  [ 409.7825   412.82828  415.509    418.12033  421.06903  424.08734
#  427.04297  430.00266  433.35428  437.05194  440.9941   445.18448
#  449.2819   453.2831   456.9953   461.1733   465.20343  469.1874
#  473.88727  478.61108  483.2348   488.4133   494.06158  499.5955
#  506.54593  512.4629   517.3102   523.39716  529.87646  536.9373
#  543.36755  549.267    555.79645  563.5755   571.3927   578.5842
#  585.53705  593.60144  602.241    610.4474   618.25073  626.76404
#  635.43066  644.4388   652.6582   661.71936  672.028    681.7412
#  692.13153  702.619    712.17834  722.7267   735.1587   746.3799
#  756.74493  768.2405   781.1833   795.1525   807.44574  819.61
#  832.7382   846.49524  858.7933   871.0713   883.37244  897.932
#  911.31555  923.92834  937.2589   951.63324  965.309    978.7125
#  992.0199  1006.58936 1019.85376 1032.9076  1047.6554  1062.86
# 1076.763   1089.8237  1104.1278  1119.0787  1133.0673  1146.7651
# 1159.7842  1173.997  ]
#
# calendar_co2mass_kg_esm_ssp585_r1i1p1f2        =  [2015 2016 2017 2018 2019 2020 2021 2022 2023 2024 2025 2026 2027 2028
# 2029 2030 2031 2032 2033 2034 2035 2036 2037 2038 2039 2040 2041 2042
# 2043 2044 2045 2046 2047 2048 2049 2050 2051 2052 2053 2054 2055 2056
# 2057 2058 2059 2060 2061 2062 2063 2064 2065 2066 2067 2068 2069 2070
# 2071 2072 2073 2074 2075 2076 2077 2078 2079 2080 2081 2082 2083 2084
# 2085 2086 2087 2088 2089 2090 2091 2092 2093 2094 2095 2096 2097 2098
# 2099 2100]
#
#
#
#
# difference_126Lu_minus_standard =  [  0.08184814   0.03216553   0.13171387   0.16906738   0.15847778
#  -0.02719116   0.28515625   1.0713501    1.5263062    0.8100281
#  -0.2937622   -1.1992493   -0.73062134  -1.2385254   -1.44281
#  -1.4172974   -1.751648    -1.5945435   -2.4179688   -1.9859009
#  -1.3006592   -2.712555    -3.5486145   -4.0009155   -5.8648376
#  -6.3403015   -5.470093    -5.6828613   -5.6918945   -6.2141724
#  -7.425293    -7.461731    -7.150696    -8.360596    -9.320862
#  -9.034668    -8.496155    -9.218567   -10.516418   -10.931885
#  -9.960632    -9.759094    -9.646362    -9.090576    -7.7731323
#  -8.1466675   -9.414246    -9.557129    -9.194275    -9.434387
#  -8.627197    -8.28833     -9.507874    -8.981812    -8.380676
#  -9.296509   -10.569092   -12.085999   -12.627991   -13.179993
# -14.441223   -15.57135    -15.703369   -15.689697   -14.785095
# -15.866455   -14.553223   -13.411682   -12.780579   -13.000549
# -13.348633   -12.764221   -10.759094   -10.843323   -11.480164
# -10.855286   -11.055664   -13.439331   -14.272949   -12.908691
# -11.814209   -12.122925   -11.6427     -11.89624    -11.100464
# -10.079956  ]
