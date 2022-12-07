#------------------------------------------------------------------------------#
# Classic community analysis
#------------------------------------------------------------------------------#
# Set directory
setwd("D:/Desktop/2022_edna")
if (!dir.exists("fig")) dir.create("fig")

# Load library
library(vegan)
library(Rtsne)
library(maps)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

# Read data
name <- read.csv("data/20211110_species.csv", header=TRUE)
site <- read.csv("data/20211110_site.csv"   , header=TRUE)
cm   <- read.csv("data/20211110_community_matrix.csv", header=TRUE)
site$district <- with(site, factor(district, levels=unique(district)))

# Map Polygon
Jpn <- map_data("world", region="Japan")

#------------------------------------------------------------------------------#
# Remove species complex
idx  <- which(sapply(strsplit(name[,"Scientific_name"], " spp "), length)==2)
name <- name[-idx,]
cm   <- cm  [-idx,]
rm(idx)

# Rarefaction
S  <- vegan::specnumber(cm, MARGIN=2)
Sr <- vegan::rarefy(cm, min(colSums(cm)), MARGIN=2)
rc <- vegan::rarecurve(t(cm), step=20, sample=min(colSums(cm)), label=FALSE)
rc <- lapply(seq_along(rc), function(k)
    data.frame(
        n = attr(rc[[k]], "Subsample"),
        S = rc[[k]],
        LID = k, district = site$district[k]
    )
)
rc <- do.call(rbind, rc)

gp_1 <- ggplot(data.frame(S=S, Sr=Sr, D=site$district)) +
    geom_abline(slope=1, intercept=0) +
    geom_point(aes(x=S, y=Sr, color=D)) +
    labs(color="District") +
    theme_bw() + xlab("Species Richness") + ylab("Rarefied Species Richness")

gp_2 <- ggplot(rc, aes(x=n, y=S, fill=as.factor(LID))) + 
    geom_vline(xintercept=min(colSums(cm))) +
    geom_line(aes(color=district)) +
    labs(color="District") +
    theme_bw() + xlab("Read Numbers") + ylab("Species Richness")

ggsave("fig/Fig_E1_rarefy.png", plot_grid(gp_1, gp_2), width=8, height=4)
rm(rc, gp_1, gp_2)

# Sampling site map
hull <- lapply(levels(site$district), function(x) {
    siteX <- site[site$district==x,]
    siteX[chull(siteX$lon, siteX$lat),]
})
hull <- do.call(rbind, hull)

text <- rbind(
    data.frame(x=139.0, y=45.0, district="HKD"),
    data.frame(x=143.8, y=39.0, district="EMP"),
    data.frame(x=137.5, y=41.0, district="EMJ"),
    data.frame(x=137.0, y=32.2, district="WMP"),
    data.frame(x=130.0, y=36.5, district="WMJ"),
    data.frame(x=127.5, y=32.0, district="WMI"),
    data.frame(x=141.5, y=32.0, district="IIO"),
    data.frame(x=126.5, y=29.0, district="INS")
)

gp_1 <- ggplot(site, aes(x=lon, y=lat, fill=district)) +
    coord_fixed(xlim=range(Jpn$long), ylim=range(Jpn$lat)) +
    geom_polygon(aes(x=long, y=lat, group=group, fill=region), data=Jpn, fill="grey70") +
    geom_polygon(aes(color=district), alpha=0.2, data=hull) +
    geom_point(shape=21) +
    geom_text(aes(x=x, y=y, label=district, color=district), data=text, size=4) +
    guides(fill="none", color="none") +
    theme_bw() + xlab("Longitude") + ylab("Latitude")

gp_2 <- ggplot(cbind(site, S=Sr), aes(x=district, y=S, fill=district)) +
    geom_jitter(aes(color=district), width=0.1) +
    geom_boxplot(alpha=0.2) +
    labs(fill="District", color="District") +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    theme_bw() + xlab("District") + ylab("Species Richness")

ggsave("fig/Fig_1a_classic.png", gp_1, width=4, height=4)
ggsave("fig/Fig_1b_classic.png", gp_2, width=3, height=4.1)
rm(hull, text)

# Dimensional reduction
dc <- vegan::vegdist(1*(t(cm)>0), method="jaccard")
tsne <- Rtsne(dc, dims=2, perplexity=3, is_distance=TRUE, eta=5)
tsne <- cbind(site, axis1=tsne$Y[,1], axis2=tsne$Y[,2])

tv <- apply(as.matrix(dc), 2, FUN=function(x) x <= sort(x)[3])
tv <- tv | t(tv)
tv <- which(tv & upper.tri(tv), arr.ind=TRUE)
tv <- data.frame(
    xi = tsne$axis1[tv[,1]], xe = tsne$axis1[tv[,2]],
    yi = tsne$axis2[tv[,1]], ye = tsne$axis2[tv[,2]]
)

gp_3 <- ggplot(tsne, aes(x=axis1, y=axis2)) +
    geom_segment(aes(x=xi, y=yi, xend=xe, yend=ye), data=tv, color="grey") +
    stat_ellipse(aes(color=district, fill=district), geom='polygon', alpha=0.2, level=0.5) +
    geom_point(aes(color=district)) +
    guides(color="none", fill="none") +
    theme_bw() + xlab("Axis 1") + ylab("Axis 2")

ggsave("fig/Fig_1c_classic.png", gp_3, width=3, height=3)
rm(dc, tsne, tv)

# Ranked Ocuupancy
df = data.frame(
    rank = 1:nrow(cm), Rank = "4-50",
    ocp  = sort(vegan::specnumber(cm), decreasing=TRUE)
)
df$Rank[df$ocp <  4] = "1-3"
df$Rank[df$ocp > 50] = "51-285"

gp_4 = ggplot(df, aes(x=rank, y=ocp, color=Rank)) +
    geom_area(aes(fill=Rank), alpha=0.5) +
    geom_line(data=data.frame(rank=rep( 86, 2), ocp=c(0,50), Rank=rep("51-285",2))) +
    geom_line(data=data.frame(rank=rep(656, 2), ocp=c(0, 3), Rank=rep(   "1-3",2))) +
    theme_bw() + xlab("Rank") + ylab("Occupancy") +
    theme(
        legend.position = c(1,1), legend.justification = c(1.1,1.1),
        legend.background = element_rect(fill="white", color="black"))

ggsave("fig/Fig_1d_classic.png", gp_4, width=3, height=3)
rm(df)

# End