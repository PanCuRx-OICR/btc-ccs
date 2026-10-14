# Packages used by docs/*.Rmd that are not on conda: install from CRAN / GitHub.
# Run inside the btc-ccs-analysis env (install.sh --analysis does this for you).
options(repos = c(CRAN = "https://cloud.r-project.org"))
need <- function(pkg) !requireNamespace(pkg, quietly = TRUE)

for (p in c("Rediscover", "ActiveDriverWGS")) if (need(p)) install.packages(p)

gh <- c(PNWColors      = "jakelawlor/PNWColors",
        ggsankey       = "davidsjoberg/ggsankey",
        ggdendroplot   = "NicolasH2/ggdendroplot",
        CMScaller      = "peterawe/CMScaller",
        dndscv         = "im3sanger/dndscv",
        gTrack         = "mskilab/gTrack",
        fishHook       = "mskilab/fishHook")
for (p in names(gh)) if (need(p)) remotes::install_github(gh[[p]], upgrade = "never")
