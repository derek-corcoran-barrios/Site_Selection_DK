Site selection for field work
================

## Steps to follow in this document:

1.  Will divide the NOVANA dataset in two datasets, the “Regularly
    sampled” and the “infrequently sampled”
2.  Within the “Regularly sampled”, I will pick the top n samples in
    terms of directionality for biodiversity, so the top lets say the
    top 34 in gained diversity, top 34 in lost biodiversity, and top 34
    in stability, so that is 103 sites.
3.  The selected points in step 2 will be added to the Sampled dataset
    together with the Biowide, Microflora danica, and Agriculture sites
    (That is how I am calling the Moegens/ Mette sites)
4.  Of all remaining sites either in the regularly sampled and
    infrequently sampled figure out the top 1000 that are the furthest
    away environmentally (Ranked)
5.  Hopefully after Monday’s meeting with Christian, further use the
    abiotic variables in NOVANA, to re-rank the top 1000 using the
    abiotic dataset within the abiotic dataset

## load data and packages

### Load needed packages:

``` r
# for shapefile manipulation
library(sf)
# for data wrangling and ploting
library(tidyverse)
# for community analysis
library(vegan)
# for efficient data manipulation
library(data.table)
# For clean summaries
library(broom)
```

### Load and modify datasets

Read a polygon of Denmark to the show where the points will be

``` r
Denmark <- readRDS("DK_Shape.rds")
```

Read in the *biowide* dataset together and modify it to have only a
unique identifier, a dataset column and the coordinates as an SF object.

``` r
# read in the data
Biowide <- read_sf("Biowide_Naturtyper.shp")

# Drop the first row
Biowide <- Biowide[-1,] %>% 
  # add the Dataset column
  mutate(Dataset = "Biowide") %>% 
  # Select only that column
  dplyr::select(Dataset) %>%
  # generate an unique identifier for each row
  tibble::rowid_to_column(var = "ID") %>% 
  # Add the Biowide prefix to the ID
  mutate(ID = paste0("Biowide_", ID))
```

Read in the *Agriculture* dataset together and modify it to have only a
unique identifier, a dataset column and the coordinates as an SF object.

``` r
# Read in the dataset
AgriculturalPoints <- read_sf("150pkt/150pkt_indenforpoly_13042021.shp") %>% 
  # Generate the Dataset column
  mutate(Dataset = "Agriculture") %>%
  # Select only that column
  dplyr::select(Dataset) %>%
  # Generate a unique identifier for the dataset
  tibble::rowid_to_column(var = "ID") %>% 
  # Add the Ag prefix for this dataest
  mutate(ID = paste0("Ag_", ID))
```

Read in the *Microflora danica* dataset together and modify it to have
only a unique identifier, a dataset column and the coordinates as an SF
object.

``` r
# Read in the dataset as a data frame
MicroFlora <- read_csv("natsoil.csv") %>% 
  # eliminate the points with lat and lon = 0
  dplyr::filter(latitude != 0 & longitude != 0) %>%
  # Transform to sf
  st_as_sf(coords  = c(4,3), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0") %>%
  # Reproject 
  st_transform(crs = st_crs(Biowide)) %>% 
  # Add dataset column
  mutate(Dataset = "Microflora Danica") %>%
  # Keep only the dataset column
  dplyr::select(Dataset) %>%
  # Generate a unique identifier for the dataset
  tibble::rowid_to_column(var = "ID") %>% 
  # Add the Micro prefix for this dataest
  mutate(ID = paste0("Micro_", ID))
```

``` r
NovanaFiles <- list.files(path = "Novana/", full.names = T) 
  
Novana <- fread(file=NovanaFiles[1]) %>% as.data.frame()  %>%
  unite(col = "ID", site, plot) %>% 
  dplyr::select(ID, UTMx, UTMy) %>% distinct() %>% 
  mutate(Dataset = "Novana", UTMx = as.numeric(UTMx), UTMy = as.numeric(UTMy)) %>% 
  dplyr::filter(!is.na(UTMx),!is.na(UTMy)) %>% 
  st_as_sf(coords = c("UTMx","UTMy"), crs = st_crs(Biowide))
```

## Diversity

For each site we will calculate species richness based on the
Novana//alledata-frekvens.csv file:

``` r
ForRichness <- fread(file=NovanaFiles[1], na.strings = "mv") %>% as.data.frame()  %>%
  unite(col = "ID", site, plot)  %>% 
  dplyr::select(ID, antalstjernearter, year) %>% 
  rename(Richness = antalstjernearter)  %>% 
  dplyr::filter(!is.na(Richness))
```

Most of the 97501 have only been sampled one year, so we filter to get
only the ones that have been sampled four or more different years.

``` r
Ns <- ForRichness %>% 
  group_by(ID) %>% 
  summarise(n = n()) %>% 
  arrange(n) 
## Keep only the IDs of the sites that have been sampled 4 or more times
TestIDs <- Ns %>% dplyr::filter(n >= 4) %>% pull(ID)

For_Analysis <- ForRichness %>% 
  dplyr::filter(ID %in% TestIDs) %>% 
  group_split(ID)
```

After that, a linear regression is made for each site and we keep the
p-values and slope estimates

``` r
Analysis_Name <- For_Analysis %>% 
  purrr::map(~pull(.x, ID)) %>% 
  reduce(c) %>% 
  unique()

Analysis <- For_Analysis %>% 
  purrr::map(~lm(Richness ~ year, data = .x)) %>% 
  purrr::map(broom::tidy) %>% 
  purrr::map(~dplyr::filter(.x, term == "year")) %>% 
  purrr::map2(Analysis_Name, ~mutate(.x, ID = .y)) %>% 
  reduce(bind_rows)
```

If p.values are higher that 0.05 we keep these sites as stable sites,
since they do not trend towards higher or lower richness

``` r
Stable <- Analysis %>% 
  dplyr::filter(p.value >= 0.05) %>% 
  dplyr::mutate(estimate)
```

On the other side if p values are bellow 0.05, we check if the slope
estimate is positive or negative in order to separate in increasing an
decreasing richness sites.

``` r
Upwards <- Analysis %>% 
  dplyr::filter(p.value < 0.05, estimate > 0)

Downwards <- Analysis %>% 
  dplyr::filter(p.value < 0.05, estimate < 0)
```

After that analysis we end up with 6473 stable sites, 279 sites with
increasing richness, and 283 sites with increasing richness.

Now for each of those datasets we select either the 34 smallest absolute
values of the slope estimate (for stable sites)

``` r
Stable <- Stable %>% 
  dplyr::mutate(AbsSlope = abs(estimate)) %>% 
  dplyr::slice_min(order_by = AbsSlope, n = 34)

StableTrends <- For_Analysis %>% 
  reduce(bind_rows) %>% 
  dplyr::filter(ID %in% Stable$ID) %>% 
  mutate(Class = "Stable")
```

The 34 highest slope estimates (for increasing diversity sites)

``` r
Upwards <- Upwards %>% 
    dplyr::slice_max(order_by = estimate, n = 34)

UpwardsTrends <- For_Analysis %>% 
  reduce(bind_rows) %>% 
  dplyr::filter(ID %in% Upwards$ID)  %>% 
  mutate(Class = "Increase")
```

The 34 lowest slope estimates (for decreasing diversity sites)

``` r
Downwards <- Downwards %>% 
    dplyr::slice_min(order_by = estimate, n = 34)

DownwardsTrends <- For_Analysis %>% 
  reduce(bind_rows) %>% 
  dplyr::filter(ID %in% Downwards$ID) %>% 
  mutate(Class = "Decrease")
```

Which leads to the following trend

``` r
DiversTrend <- list(StableTrends, DownwardsTrends, UpwardsTrends) %>% 
  reduce(bind_rows)

ggplot(DiversTrend, aes(x = year, y = Richness)) + geom_path(aes(group= ID, color = ID)) + geom_point(aes(group= ID, color = ID)) + theme_bw() + theme(legend.position = "none") + facet_wrap(~Class)
```

![](README_files/figure-gfm/Trends-1.png)<!-- -->

Those sites shown in a map:

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## Join together and visualize

This is map of the currently selected sites:

``` r
AllData <- list(Biowide, AgriculturalPoints, MicroFlora) %>% 
  # join all together
  purrr::reduce(rbind) %>% 
  mutate(Class = NA)
```

``` r
DiverMap <- DiverMap %>% 
  mutate(Dataset = paste(Dataset, Class, sep = "_"))
AllData <- rbind(AllData, DiverMap)
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
AllData <- AllData %>% 
  mutate(Selected = "Yes")
Novana <- Novana %>% 
  dplyr::filter(!(ID %in% AllData$ID)) %>% 
  mutate(Class = NA, 
         Selected = "No")

AllData <- AllData %>% rbind(Novana)

saveRDS(AllData, "AllData.rds")
```
