#------------------------------------------------------------------------------#
# Analytical predictions
#------------------------------------------------------------------------------#

pred_EF = function (S=50, c=0.3, p=0.5, N=10, ...)
{
    # Default setting
    Er = 1.5; Vr = 1/3; Es = 6.0; Vs = 1/3
    Ea = Eb = 1/2; Va = Vb = 1/12
    SA = SB = PR = PC = 0   # symmetry
    
    dots <- list(...)
    if (!is.null(dots$Ea)) Ea <- dots$Ea
    if (!is.null(dots$Eb)) Eb <- dots$Eb
    if (!is.null(dots$Er)) Er <- dots$Er
    if (!is.null(dots$Es)) Es <- dots$Es
    if (!is.null(dots$Va)) Ea <- dots$Va
    if (!is.null(dots$Vb)) Eb <- dots$Vb
    if (!is.null(dots$Vr)) Er <- dots$Vr
    if (!is.null(dots$Vs)) Es <- dots$Vs
    if (!is.null(dots$SA)) SA <- dots$SA
    if (!is.null(dots$SB)) SB <- dots$SB
    if (!is.null(dots$PR)) PR <- dots$PR
    if (!is.null(dots$PC)) PC <- dots$PC
    
    # Derived variables
    Dr <- Vr/Er^2
    Ds <- Vs/Es^2
    cp <- c*p
    cq <- c - cp
    
    # Summary variables
    VA <- (S-1)*cp*Va
    VB <- (S-1)*cq*Vb
    DA <- (S-1)*cp*Ea
    DB <- (S-1)*cq*Eb
    HA <- sqrt(Va/Ea^2/cp/(S-1) + (1-cp)/cp*(S-1)/(S*(S-1)-1))
    HB <- sqrt(Vb/Eb^2/cq/(S-1) + (1-cq)/cq*(S-1)/(S*(S-1)-1))
    
    # sum(y)
    EbRR <- DA*(1+SA*HA*HA)
    EbCC <- DB*(1+SB*HB*HB)
    EbRC <- DB*(1+PR*HA*HB)
    EbCR <- DA*(1+PC*HA*HB)
    EH <- Es^2  + Es*EbCC - Es*EbRR - EbRR*EbCC + EbCR*EbRC
    EK <- DA * (Es + EbCC - EbRC) - DB * (Es + EbCR - EbRR)
    EZ <- EH + EK
    
    EY <- N*Er*EZ/EH
    VY <-
        (1/N-1/S)/EZ^2 * Dr*EH^2 +
        (1/N-1/S)/EZ^2 * (S*(DA*HA)^2-VA)/(S-1)*(Es+EbCC-EbRC)^2 +
        (1/N-1/S)/EZ^2 * (S*(DB*HB)^2-VB)/(S-1)*(Es+EbCR-EbRR)^2 +
        1/S * Dr +
        1/S * Dr/EZ^2 * HA^2 * (Es*DA+DA*EbCC-DB*EbCR)^2 +
        1/S * Dr/EZ^2 * HB^2 * (Es*DB+DA*EbRC-DB*EbRR)^2 +
        1/S/(EH*EZ)^2 * Vs * ((DB-DA)*EH + (2*Es+EbCC-EbRR)*EK)^2 +
        1/S/(EH*EZ)^2 * Vs * HA^2 * (DB*EH+(Es+EbCC)*EK)^2 +
        1/S/(EH*EZ)^2 * Vs * HB^2 * (DA*EH-(Es-EbRR)*EK)^2 +
        1/S/(EH*EZ)^2 * VA * ((EH^2+EK^2)*(Es+EbCC-EbRC)^2 + 2*EH*EK*Es^2*EbRR/DA) +
        1/S/(EH*EZ)^2 * VB * ((EH^2+EK^2)*(Es+EbCR-EbRR)^2 + 2*EH*EK*Es^2*EbCC/DB) +
        1/S/(EH*EZ)^2 * (DA*HA*HA)^2 * (DB*EH+(Es+EbCC)*EK)^2 +
        1/S/(EH*EZ)^2 * (DB*HB*HB)^2 * (DA*EH-(Es-EbRR)*EK)^2 +
        1/S/(EH*EZ)^2 * (DA*HB*HA)^2 * (DB*EH+EbRC*EK)^2 +
        1/S/(EH*EZ)^2 * (DB*HA*HB)^2 * (DA*EH+EbCR*EK)^2
    
    # sum(EF)
    mEF <- EY/Es*(1+Ds)
    vEF <- sqrt(VY+Ds/N)/(1+Ds)
    data.frame(mean=mEF, vBC=vEF)
}

# End