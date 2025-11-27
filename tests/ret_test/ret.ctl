# ======================================================================
# Forward model...
# ======================================================================

# Table directory...
TBLBASE = ../data/boxcar

# Emitters...
NG = 5
EMITTER[0] = CO2
EMITTER[1] = H2O
EMITTER[2] = O3
EMITTER[3] = F11
EMITTER[4] = CCl4

# Channels...
ND = 2
NU[0] = 792.0000
NU[1] = 832.0000

# Retrieval parameters...
RETQ_ZMIN[3] = 5
RETQ_ZMAX[3] = 30
ERR_Q[3] = 10
ERR_Q_CZ[3] = 10
ERR_Q_CH[3] = 1e5

# Measurement errors...
ERR_NOISE[0] = 1e-5
ERR_NOISE[1] = 1e-5
ERR_FORMOD[0] = 1.0
ERR_FORMOD[1] = 1.0

# Output...
ERR_ANA = 1
WRITE_MATRIX = 1
