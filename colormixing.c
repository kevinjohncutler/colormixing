// --- LIGHTBULBS
// MACROS
#define R       lightbulb_group->r
#define G       lightbulb_group->g
#define B       lightbulb_group->b
#define CW      lightbulb_group->cw
#define WW      lightbulb_group->ww
#define WP      lightbulb_group->wp

// https://gist.github.com/rasod/42eab9206e28ca91c8d9f926fa71a938
// https://gist.github.com/unteins/6ecb69883d55ad8424b70be405bf4115
// Ths is the main color mixing function
// Takes in HSV, assigns PWM duty cycle to the lightbulb_group struct
//https://github.com/espressif/esp-idf/examples/peripherals/rmt/led_strip/main/led_strip_main.c
// ref5: https://github.com/patdie421/mea-edomus/blob/0eb0f9a8630ce610e3d1f6dd3c3a8d29d2dffea6/src/interfaces/type_004/philipshue_color.c
// Bruce Lindbloom's website http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html

// Helper function to sum arrays
float array_sum(float arr[], int num_elements) {
    float sum = arr[0];
    for (uint8_t i = 1; i < num_elements; i++) {
        sum = sum + arr[i];
    }
    
    return sum;
}

// Helper function to find max of an array
float array_max(float arr[], int num_elements) {
    float max = arr[0];
    for (uint8_t i = 1; i < num_elements; i++) {
        if (arr[i] > max) {
            max = arr[i];
        }
    }
    
    return max;
}

// Helper function to find min of an array
float array_min(float arr[], int num_elements) {
    float min = arr[0];
    for (uint8_t i = 1; i < num_elements; i++) {
        if (arr[i] < min) {
            min = arr[i];
        }
    }
    
    return min;
}

// Helper function to perform dot product
float array_dot(float arr1[], float arr2[], int num_elements) {
    float sum = arr1[0] * arr2[0];
    for (uint8_t i = 1; i < num_elements; i++) {
        sum += arr1[i] * arr2[i];
    }
    
    return sum;
}

// Helper function to multiply array by a constant
void array_multiply(float arr[], float scalar, int num_elements) {
    for (uint8_t i = 0; i < num_elements; i++) {
        arr[i] *= scalar;
    }
}

// Helper function to asign array componenets
void array_equals(float arr[], float vals[], int num_elements) {
    for (uint8_t i = 0; i < num_elements; i++) {
        arr[i] = vals[i];
    }
}

// Helper function to rescale array so that max is 1
void array_rescale(float arr[], int num_elements) {
    float amax = array_max(arr, num_elements);
    if (amax != 0) {
        array_multiply(arr, 1.0f / amax, num_elements);
    } else {
        array_multiply(arr, 0, num_elements);
    }
}

// Helper function to calculate barycentric coordinates; now returns coordinates with the max=1 always
void bary(float L[], float p0[], float p1[], float p2[], float p3[]) {
    const float denom = (p2[1] - p3[1]) * (p1[0] - p3[0]) + (p3[0] - p2[0]) * (p1[1] - p3[1]);
    L[0] = ((p2[1] - p3[1]) * (p0[0] - p3[0]) + (p3[0] - p2[0]) * (p0[1] - p3[1])) / denom;
    L[1] = ((p3[1] - p1[1]) * (p0[0] - p3[0]) + (p1[0] - p3[0]) * (p0[1] - p3[1])) / denom;
    L[2] = 1 - (L[0] + L[1]);
    // array_rescale(L, 3); // added to save redudancy when adding weights... not sure if best
}

// Get the color p0's barycentric coordinates based on where it is. Assumes p4 is within the triangle formed by p1,p2,p3.
void getWeights(float coeffs[], float p0[], float p1[], float p2[], float p3[], float p4[]) {
    float L[3];
    bary(L, p0, p1, p2, p4); // Try red-green-W
    if ((L[0] >= 0)&&(L[1] >= 0)&&(L[2] >= 0)) {
        float vals[] = {L[0], L[1], 0, L[2]};
        array_equals(coeffs,vals,4);
    } else {
        bary(L, p0, p2, p3, p4); // Try green-blue-W
        if ((L[0] >= 0)&&(L[1] >= 0)&&(L[2] >= 0)) {
            float vals[] = {0, L[0], L[1], L[2]};
            array_equals(coeffs,vals,4);
        } else { // must be red-blue-W
            bary(L, p0, p1, p3, p4);
            float vals[] = {L[0], 0, L[1], L[2]};
            array_equals(coeffs,vals,4);
        }
    }
    array_rescale(coeffs,4);
}

void hsi2rgbw(uint16_t h, uint16_t s, uint16_t v, ch_group_t* ch_group) {
    uint32_t run_time = sdk_system_get_time();
    lightbulb_group_t* lightbulb_group = lightbulb_group_find(ch_group->ch0);
    // Hue is specified in degrees, and being cyclic we want to rescale it to the default range [0,360).
    // Saturation usually normalized  0 to 100
    // Value/brightness same nomrlaization as saturation
    
    h %= 360; // shorthand modulo arithmetic, h = h%360, so that h is rescaled to the range [0,360) (angle around hue circle)
    
    uint32_t rgb_max = 1; // Ignore brightness for initial conversion
    float rgb_min = rgb_max * (100 - s) / 100.f; // Again rescaling 100->1, backing off from max
    INFO("saturation %g", rgb_min);
    uint32_t i = h / 60; // which 1/6th of the hue circle you are in
    uint32_t diff = h % 60; // remainder (how far counterclockwise into that 60 degree sector)

    float rgb_adj = (rgb_max - rgb_min) * diff / 60; // radius*angle = arc length

    float r, g, b; // declare variables
    
    // Different rules depending on the sector
    // I think it is something like approximating the RGB cube for each sector
    // Indeed, six sectors for six faces of the cube.
    INFO("light switch %i", i);
    switch (i) {
        case 0: // Red to yellow
            r = rgb_max;
            g = rgb_min + rgb_adj;
            b = rgb_min;
            break;
        case 1: // yellow to green
            r = rgb_max - rgb_adj;
            g = rgb_max;
            b = rgb_min;
            break;
        case 2: // green to cyan
            r = rgb_min;
            g = rgb_max;
            b = rgb_min + rgb_adj;
            break;
        case 3:// cyan to blue
            r = rgb_min;
            g = rgb_max - rgb_adj;
            b = rgb_max;
            break;
        case 4: // Blue to magenta
            r = rgb_min + rgb_adj;
            g = rgb_min;
            b = rgb_max;
            break;
        default:// magenta to red
            r = rgb_max;
            g = rgb_min;
            b = rgb_max - rgb_adj;
            break;
    }
    // Keeping the HSV to RGB. Seems to work fine.
   
    // (2) convert to XYZ then to xy(ignore Y). Also now apply gamma correction.
    float gcR = (r > 0.04045f) ? fast_precise_pow((r + 0.055f) / (1.0f + 0.055f), 2.4f) : (r / 12.92f);
    float gcG = (g > 0.04045f) ? fast_precise_pow((g + 0.055f) / (1.0f + 0.055f), 2.4f) : (g / 12.92f);
    float gcB = (b > 0.04045f) ? fast_precise_pow((b + 0.055f) / (1.0f + 0.055f), 2.4f) : (b / 12.92f);
    
//    float p[2] = {(B[1]*gcR*R[0]*(G[0] - WP[0]) + gcR*G[1]*R[0]*WP[0] - gcG*G[0]*R[1]*WP[0] + B[1]*gcG*G[0]*(-R[0] + WP[0]) + B[0]*gcG*G[0]*(R[1] - WP[1]) + gcG*G[0]*R[0]*WP[1] - gcR*G[0]*R[0]*WP[1] + B[0]*gcR*R[0]*(-G[1] + WP[1]) + B[0]*gcB*(G[1]*(R[0] - WP[0]) + R[1]*WP[0] - R[0]*WP[1] + G[0]*(-R[1] + WP[1])))/(-(B[0]*gcR*G[1]) + gcB*G[1]*R[0] + B[0]*gcG*R[1] - gcB*G[0]*R[1] - gcB*G[1]*WP[0] + gcR*G[1]*WP[0] + gcB*R[1]*WP[0] - gcG*R[1]*WP[0] + B[1]*(gcR*(G[0] - WP[0]) + gcG*(-R[0] + WP[0])) + (B[0]*(-gcG + gcR) + gcB*G[0] - gcR*G[0] - gcB*R[0] + gcG*R[0])*WP[1]),
//        ((gcG - gcR)*G[1]*R[1]*(B[0] - WP[0]) + (-(B[0]*gcG*G[1]) + gcG*G[1]*R[0] + B[0]*gcR*R[1] - gcR*G[0]*R[1])*WP[1] + B[1]*(gcR*R[1]*(G[0] - WP[0]) + gcG*G[1]*(-R[0] + WP[0]) + gcB*(-(G[0]*R[1]) + G[1]*(R[0] - WP[0]) + R[1]*WP[0] + G[0]*WP[1] - R[0]*WP[1])))/(-(B[0]*gcR*G[1]) + gcB*G[1]*R[0] + B[0]*gcG*R[1] - gcB*G[0]*R[1] - gcB*G[1]*WP[0] + gcR*G[1]*WP[0] + gcB*R[1]*WP[0] - gcG*R[1]*WP[0] + B[1]*(gcR*(G[0] - WP[0]) + gcG*(-R[0] + WP[0])) + (B[0]*(-gcG + gcR) + gcB*G[0] - gcR*G[0] - gcB*R[0] + gcG*R[0])*WP[1])
//    };
//    float X, Y, Z, x, y;
    
//    // sRGB
//    X = 0.463735*gcR + 0.419345*gcG + 0.113213*gcB;
//    Y = 0.233706*gcR + 0.715045*gcG + 0.0512496*gcB;
//    Z = 0.00941927*gcR + 0.0784289*gcG + 0.521612*gcB;
//
//    //
//
//    x = X/(X+Y+Z);
//    y = Y/(X+Y+Z);
//
//    float p[2] = {x,y};
    
    //sRGB Primaries
    float denom = gcB*(1.0849816845413038 - 1.013932032965188*WP[0] - 1.1309562993446372*WP[1]) + gcR*(-0.07763613020327381 + 0.6290345979071447*WP[0] + 0.28254416233165086*WP[1]) + gcG*(-0.007345554338030566 + 0.38489743505804436*WP[0] + 0.8484121370129867*WP[1]);
    float p[2] = {
        (gcB*(0.14346747729293707 - 0.10973602163740542*WP[0] - 0.15130242398432014*WP[1]) + gcG*(-0.07488078606033685 + 0.6384263445494799*WP[0] - 0.021643836986853505*WP[1]) + gcR*(-0.06858669123260035 + 0.47130967708792587*WP[0] + 0.17294626097117374*WP[1])) / denom,
        (gcB*(0.044068179152577595 - 0.039763003918887874*WP[0] - 0.023967881257084177*WP[1]) + gcR*(-0.029315872239491187 + 0.18228806745026777*WP[0] + 0.12851324243365222*WP[1]) + gcG*(-0.014752306913086446 - 0.14252506353137964*WP[0] + 0.8954546388234319*WP[1])) / denom
    };
    //P3
//    float denom = gcB*(1.1916335779794955 - 1.113599863928703*WP[0] - 1.2421274208847664*WP[1]) + gcR*(-0.1161017631704834 + 0.7742332295983754*WP[0] + 0.33261888106768217*WP[1]) + gcG*(-0.07553181480901194 + 0.3393666343303272*WP[0] + 0.909508539817084*WP[1]);
//    float p[2] = {
//        (gcB*(0.1575700914827437 - 0.1205228898885055*WP[0] - 0.1661752003911752*WP[1]) + gcG*(-0.07388956217017274 + 0.5679139660257223*WP[0] - 0.04954947848466228*WP[1]) + gcR*(-0.08368052931257094 + 0.5526089238627829*WP[0] + 0.21572467887583732*WP[1])) / denom,
//        (gcB*(0.04840000780366015 - 0.04367164101125731*WP[0] - 0.02632388408569457*WP[1]) + gcR*(-0.03615420644056613 + 0.23006095617141115*WP[0] + 0.12022104262257138*WP[1]) + gcG*(-0.012245801363094012 - 0.18638931516015383*WP[0] + 0.9061028414631231*WP[1])) / denom
//    };
    
    float coeffs[5] = {0,0,0,0,0};
    
    // (3) Figure out where p is in the chromaticity diagram and find barycentric coordinates. With 5 LEDs, there will be 7 regions.
    // Note that mow bary should be directly changing the elements of the weight array passed in as its first argument.

    float targetRGB[3];
    bary(targetRGB, p, R, G, B);
    
    // might want to shift the target RGB towards red etc. Currently this is pulling towards my/user defined primaries; might want to change this to a set of standard primaries if it works well.
    targetRGB[0] *= LIGHTBULB_FACTOR_R;
    targetRGB[1] *= LIGHTBULB_FACTOR_G;
    targetRGB[2] *= LIGHTBULB_FACTOR_B;
    
    array_multiply(targetRGB,1/array_max(targetRGB,3),3); // Never can be all zeros, ok to divide; just to max out to do extraRGB
    
    
    // NOTE: NEW IDEA: start with RGBW and if there are 5 channels then add on the RGBWW. This might max out thw whites better, as it guarantees both whites are always used.
    
    if ((targetRGB[0] >= 0)&&(targetRGB[1] >= 0)&&(targetRGB[2] >= 0)&&(s!=100)) { // within gamut and not at the boundary (save computation!)
        
        // RGBW assumes W is in CW position
        float coeffs1[4];
        getWeights(coeffs1, p, R, G, B, CW);
        for (i = 0; i < 4; i++) {
            coeffs[i] += coeffs1[i];
        }
        // If WW, then compute RGBW again with WW, add to the main coeff list. Can easlily add support for any LEDs inside the gamut. The only thing I am worried about is that the RGB is always pulling double-duty with two vertexes instead of so
        if ((uint8_t) LIGHTBULB_CHANNELS == 5) {
            float coeffs2[4];
            getWeights(coeffs2, p, R, G, B, WW);
            for (i = 0; i < 3; i++) {
                coeffs[i] += coeffs2[i];
            }
            coeffs[4] += coeffs2[3];
        }

        INFO("coeffs before flux: %g, %g, %g, %g, %g",coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4]);
        
        // (3.a.0) Correct for differences in intrinsic flux; needed before extraRGB step because we must balance RGB to whites first to see what headroom is left
        for (i = 0; i < 5; i++) {
            if (lightbulb_group->flux[i] != 0) {
                coeffs[i] /=  lightbulb_group->flux[i];
            } else {
                coeffs[i] = 0;
            }
        }

        // (3.a.1) Renormalize the coeffieients so that no LED is over-driven
        array_rescale(coeffs,5);
        
        // (3.a.2) apply a correction to to scale down the whites according to saturation. rgb_min is (100-s)/100, so at full saturation there should be no whites at all. This is non-physical and should be avoided. Flux corrections are better.
        float a = LIGHTBULB_CURVE_FACTOR;
        if (a!=0) {
                array_multiply(coeffs,1-(expf(a*s/100.f)-1)/(expf(a)-1),5);
        }
        
        // (3.a.3) Calculate any extra RGB that we have headroom for given the target color
        float extraRGB[3] = {targetRGB[0],targetRGB[1],targetRGB[2]}; // initialize the target for convenienece
        // Adjust the extra RGB according to relative flux; need to rescale those fluxes to >=1 in order to only shrink components, not go over calculated allotment
        float rflux[3] =  {lightbulb_group->flux[0],lightbulb_group->flux[1],lightbulb_group->flux[2]};
        array_multiply(rflux,1.0f/array_min(rflux,3),3); // assumes nonzero fluxes
        for (i = 0; i < 3; i++) {
            if (rflux[i] != 0) {
                extraRGB[i] /= rflux[i];
            } else {
                extraRGB[i] = 0;
            }
        }

        float loRGB[3] = {1-coeffs[0], 1-coeffs[1], 1-coeffs[2]}; // 'leftover' RGB
        if ((loRGB[0]>=0)&&(loRGB[1]>=0)&&(loRGB[2]>=0)) { // this test seems totally unecessary
            float diff[3] = {extraRGB[0]-loRGB[0],extraRGB[1]-loRGB[1],extraRGB[2]-loRGB[2]};
            float maxdiff = array_max(diff,3);
            if ((maxdiff==diff[0])&&(extraRGB[0]!=0)) {
                array_multiply(extraRGB,loRGB[0]/extraRGB[0],3);
            } else if ((maxdiff==diff[1])&&(extraRGB[1]!=0)) {
                array_multiply(extraRGB,loRGB[1]/extraRGB[1],3);
            } else if ((maxdiff==diff[2])&&(extraRGB[2]!=0)) {
                array_multiply(extraRGB,loRGB[2]/extraRGB[2],3);
            }
        } else {
            array_multiply(extraRGB,0,3);
        }
        INFO("extraRGB: %g, %g, %g", extraRGB[0],extraRGB[1],extraRGB[2]);

        // (3.a.4) Add the extra RGB to the final tally
        coeffs[0] += extraRGB[0];
        coeffs[1] += extraRGB[1];
        coeffs[2] += extraRGB[2];
    
    } else { // (3.b.1) Outside of gamut; easiest thing to do is to clamp the barycentric coordinates and renormalize... this might explain some bluish purples?

        
        for (i = 0; i < 3; i++) {
            if (targetRGB[i] < 0) {
                targetRGB[i] = 0;
            }
            if (targetRGB[i] > 1) { // Used to be redundant, now with color factors is useful
                targetRGB[i] = 1;
            }
        }

        //Idk why this would be needed
//        float totalRGB = array_sum(targetRGB,3);
//        array_multiply(targetRGB,1/totalRGB,3);
        
        float vals[5] = {targetRGB[0]/lightbulb_group->flux[0],targetRGB[1]/lightbulb_group->flux[1],targetRGB[2]/lightbulb_group->flux[2],0,0};
        array_equals(coeffs,vals,5);
        INFO("Out of gamut: %g, %g, %g", targetRGB[0],targetRGB[1],targetRGB[2]);
    }
    
    // (4) Apply color factors to whites if needed (really should not be)
    coeffs[3] *= LIGHTBULB_FACTOR_CW;
    coeffs[4] *= LIGHTBULB_FACTOR_WW;
  
    // Rescale to make sure none are overdriven
    array_rescale(coeffs,5);
    
    // (5) Brightness defined by the normalized value argument. We also divide by the scale found earlier to amp the brightness to maximum when the value is 1 (v/100). We also introduce the PWM_SCALE as the final 'units'.
    // Max power cutoff: want to limit the total scaled flux. Should do sum of flux times coeff, but what should the cutoff be? Based on everything being on, i.e. sum of fluxes.

    float brightness = (v/100.f)*PWM_SCALE;
    INFO("brightness: %g",brightness);

    if (LIGHTBULB_MAX_POWER!=1) {
        float flux_ratio = array_dot(lightbulb_group->flux,coeffs,5) / array_sum(lightbulb_group->flux,4); // actual brightness compared to theoretrical max brightness
        brightness *= MIN(LIGHTBULB_MAX_POWER,1.0f) / flux_ratio; // Hard cap at 1 as to not accidentally overdrive
    }
    
    uint32_t r_final, g_final, b_final, cw_final, ww_final;

    r_final = round(coeffs[0]*brightness);
    g_final = round(coeffs[1]*brightness);
    b_final = round(coeffs[2]*brightness);
    cw_final = round(coeffs[3]*brightness);
    ww_final = round(coeffs[4]*brightness);

    // (6) Assign the target colors to lightbulb group struct, now in fraction of PWM_SCALE. This min function is just a final check, it should not ever go over PWM_SCALE.
    lightbulb_group->target_r  = MIN(r_final,PWM_SCALE);
    lightbulb_group->target_g  = MIN(g_final,PWM_SCALE);
    lightbulb_group->target_b  = MIN(b_final,PWM_SCALE);
    lightbulb_group->target_cw = MIN(cw_final,PWM_SCALE);
    lightbulb_group->target_ww = MIN(ww_final,PWM_SCALE);
    
    INFO("hsi2rgbw runtime: %0.3f ms", ((float) (sdk_system_get_time() - run_time)) * 1e-3);
}
