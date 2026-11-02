#ifndef enumeration_h
#define enumeration_h



enum CaloSampling{

    PSB       = 0,
    PSE       = 1,
    EMB1      = 2,
    EMB2      = 3,
    EMB3      = 4,
    TileCal1  = 5,
    TileCal2  = 6,
    TileCal3  = 7,
    TileExt1  = 8,
    TileExt2  = 9,
    TileExt3  = 10,
    EMEC1     = 11,
    EMEC2     = 12,
    EMEC3     = 13,
    HEC1      = 14,
    HEC2      = 15,
    HEC3      = 16,
    // First sampling of the ATLAS-like HEC (option --atlas-hec of simu_trf.py/digit_trf.py; not used by default).
    // Outside the range of the others so that the hash (sampling + 17 side) 1e7 + segment 1e6 + bin stays unique
    // (77e7 on side A, 94e7 on side B; the other samplings use up to 50e7).
    HEC0      = 60,
};



enum Detector{
    LAR = 0,
    TILE = 1,
    TTEM = 2,
    TTHEC = 3,
    FCALEM = 5,
    FCALHAD = 6,
};




#endif
