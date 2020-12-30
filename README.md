# N-Channel Color Mixing
Color mixing algorithm and suplmementary functions for [HAA](https://github.com/RavenSystem/esp-homekit-devices), soon to be merged into the main branch there. Aside from the original HSV to RGB conversion, the rest is new to HAA. I have implemented standard RGB linearization and conversion to chromaticity coordinates with user-tunable reference white. What may not be standard, as far as I can tell, is my algrothm for maximizing LED utilization for any number of distinct channels, i.e. LEDs with different chromaticity coordinates. This makes LED bulbs and strips as bright as possible while maintaining accurate color. 

The algorithm as presented here is coded for 3, 4, or 5 LEDs arranged in a triangular gamut. Typically these will be RGB, RGBW, or RGB-CW-WW arrays. Each channel may differ in luminous flux by either having more LEDs (typically far more white LEDs than RGB) or by intrinsic brightness, and this code demonstrates how to correct for these differences in flux. 

Given a target color converted to xy space, it will either lie outside the LED gamut or inside one or more sub-gamuts formed by triplets of LEDs. Barycentric coordinates for the color are then computed to establish the releative weight of each of the three LEDs. If the point is out-of-gamut, the barycentric coordinates are clipped between 0 and 1 (projected onto the gamut edge). If the point with in a 3-channel system, the barycentric weights are rescaled such that the brightest LED is eighted to 1. If in a 4 channel system (RGBW), the gamut is split into 3 sub-gamuts: RGW,GBW,BRW. The barycentric weights are then computed for the gamut in which the point is encompassed. For a 5-channel system (RGB-CW-WW), we perform the same computation for the RGB-CW sub-gamut and the RGB-WW sub-gamut, performing the same rescaling before adding corresponding channels together and rescaling again. In this way, a 5-channel system will have weights ranging from {0,0,0,0,0} to {1,1,1,1,1}. This operation can be extended for arbitrarily many white-like LEDs (i.e. within the RGB gamut) with RGB-W1,RGB-w2,...,RGB-Wn sub-gamuts.

Before proceeding, these weights are divided by the flux for each channel, because the barycentric coordinate computation assumes equal weights for each vertex to start. Then the weights are rescaled again. Lastly, the exterior gamut (RGB) gamut weights are inspected. If each weight is less than 1, then we are able to use the RGB LEDs to contribute directly to the target color as well, guatanteeing that at least one channel weight will be saturated to 1. 

I am interested in publishing a paper detailing why this algorithm is optimal, but for now it suffices to say that a large swath of white and off-white colors correspond to 3 LEDs at full brightness (weight 1). Only a single color will correpond to 5, and a set of lines will correspond to 4. 

The code presented here of course requires known LED chromaticity coordinates provided as descibed in the HAA documentation. Other parameters include a 'curve factor' to dampen the white LED contributions based on saturation. Currently a saturation of 100 will force the white sub-gamuts not to be considered, but the normalized exponential I have included allos this to be a tuned transition. 
