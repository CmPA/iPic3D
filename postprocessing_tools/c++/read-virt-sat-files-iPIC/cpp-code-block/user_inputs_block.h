    string XYZ_block = "XY";                // dimension to scan for satellites
    float x_satellite_user = 13.0859;         // x-coordinate
    float y_satellite_user = 6.44532;         // y-coordinate
    float z_satellite_user = 0.195311;         // z-coordinate
    int const nproc = 3071;                           // number of files to search into
    float dt = 0.125;                                     // timestep
    float Bx0 = 0.0097;                                   // base field x
    char output_filename[] = "sat_output_block.txt";    // output file
