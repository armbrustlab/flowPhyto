
.rainbow.cols <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow","#FF7F00", "red", "#7F0000"))

REPO.PATH <- '/share/data/cruise/instrument/underway/seaflow/'

## CMOP ONLY
.EVT.HEADER.V1 <- c("time","pulse_width","D1","D2",
                     "FSCbig","CHLsmall",
                     "FSCperp","PE",
                     "FSCsmall","CHLbig")

## TARA ONWARD
.EVT.HEADER.V2 <- c("time","pulse_width","D1","D2",
                     "FSCsmall",
                     "FSCperp","FSCbig",
                     "PE",
                     "CHLsmall","CHLbig")



## 2010.10.20 - lowercase + underscore for db compatbility and easy adding
CHANNEL.CLMNS <- c("fsc_small","fsc_perp","fsc_big","pe","chl_small","chl_big")
.EVT.HEADER.V3 <- c("time","pulse_width","D1","D2", CHANNEL.CLMNS)                     

CHANNEL.CLMNS.SM <- CHANNEL.CLMNS[grep('big',CHANNEL.CLMNS, invert=TRUE)]
                     
EVT.HEADER <- .EVT.HEADER.V3


