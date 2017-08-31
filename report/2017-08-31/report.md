# DNA Damage Estimation
Peter Diakumis  
31 August 2017  



```r
library(ggplot2)
library(dplyr)
library(readr)
```

## Data
Used chr21 `BAM` subset of following samples:

```r
fnames <- list.files("~/Desktop/git/dna-damage-estimation/data/new")
samples <- fnames %>%
  sub("(.*)_chr21_(tot|pos|con)_damage.tsv", "\\1", .) %>% 
  sub("-ready", "", .) %>% 
  sub("-sort", "", .) %>% 
  sub("-germline", "_germ", .) %>% 
  unique()
samples
```

```
 [1] "MH17B001P004_germ" "MH17B001P004"      "MH17B001P013_germ"
 [4] "MH17B001P013"      "MH17T001P004"      "MH17T001P013"     
 [7] "WES005PVFE"        "WES005PVFR"        "WES006JCFE"       
[10] "WES006JCFR"        "WES011RCFE"        "WES011RCFR"       
[13] "WES013PFFR"       
```

## Estimate total, positional, and context damage

### Functions

```r
read_tot_damage <- function(fname) {
  readr::read_tsv(fname, 
                  col_names = c("abs", "type", "id", "count", "family", "damage"),
                  col_types = "iccdcd") %>% 
    arrange(type)
}

read_pos_damage <- function(fname) {
  readr::read_tsv(fname,
                 col_names = c("id", "type", "read", "rel", "abs", "loc"),
                 col_types = "cccdii") %>% 
    arrange(type, read, loc)
}

read_con_damage <- function(fname) {
  readr::read_tsv(fname,
                  col_names = c("id", "type", "read", "rel", "loc", "context", "abs"),
                  col_types = "cccdici") %>% 
    arrange(type, read, loc, context)
}
```


### Plots


#### Total Damage

```r
fnames <- list.files("~/Desktop/git/dna-damage-estimation/data/new",
                     pattern = "tot_damage.tsv$", full.names = TRUE)
basename(fnames)
```

```
 [1] "MH17B001P004-germline-ready_chr21_tot_damage.tsv"
 [2] "MH17B001P004-ready_chr21_tot_damage.tsv"         
 [3] "MH17B001P013-germline-ready_chr21_tot_damage.tsv"
 [4] "MH17B001P013-ready_chr21_tot_damage.tsv"         
 [5] "MH17T001P004-ready_chr21_tot_damage.tsv"         
 [6] "MH17T001P013-ready_chr21_tot_damage.tsv"         
 [7] "WES005PVFE_chr21_tot_damage.tsv"                 
 [8] "WES005PVFR_chr21_tot_damage.tsv"                 
 [9] "WES006JCFE_chr21_tot_damage.tsv"                 
[10] "WES006JCFR_chr21_tot_damage.tsv"                 
[11] "WES011RCFE_chr21_tot_damage.tsv"                 
[12] "WES011RCFR_chr21_tot_damage.tsv"                 
[13] "WES013PFFR-sort_chr21_tot_damage.tsv"            
```

```r
type_clean <-c("G_T", "C_A", "C_T", "G_A", "T_A", "A_T",
               "A_G", "T_C", "C_G", "G_C", "T_G", "A_C")
mut_list <- lapply(fnames, function(f) {
  mut <- read_tot_damage(f) %>% 
    filter(type %in% type_clean) %>%
    mutate(type = factor(type, level = type_clean))
  mut
})

mut_all <- dplyr::bind_rows(mut_list)
mut_all <- mut_all %>% 
  mutate(id2 = sub("^\\/(.*).*_chr21.*", "\\1", id),
         id2 = sub("-(sort|ready)", "", id2),
         id2 = sub("-germline", "-germ", id2))

table(mut_all$id2, useNA = "ifany")
```

```

     MH17B001P004 MH17B001P004-germ      MH17B001P013 MH17B001P013-germ 
               12                12                12                12 
     MH17T001P004      MH17T001P013        WES005PVFE        WES005PVFR 
               12                12                12                12 
       WES006JCFE        WES006JCFR        WES011RCFE        WES011RCFR 
               12                12                12                12 
       WES013PFFR 
               12 
```

```r
# no log
g <- ggplot(mut_all, aes(x = reorder(type, damage), y = damage, color = id2))
g + geom_point(alpha = 1/2, size=1.5) +
  geom_hline(yintercept = 1.5, color = "#990000", linetype = "dashed") +
  annotate("text", x = 4, y = 1.6, color = "#990000",
           label = "1/3 of variants is due to damage") +
  geom_hline(yintercept = 0, color = "grey") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  ggtitle("GIV scores for variant types")
```

![](report_files/figure-html/tot_damage-1.png)<!-- -->
 
#### Position Damage

```r
fnames <- list.files("~/Desktop/git/dna-damage-estimation/data/new",
                     pattern = "pos_damage.tsv$", full.names = TRUE)
basename(fnames)
```

```
 [1] "MH17B001P004-germline-ready_chr21_pos_damage.tsv"
 [2] "MH17B001P004-ready_chr21_pos_damage.tsv"         
 [3] "MH17B001P013-germline-ready_chr21_pos_damage.tsv"
 [4] "MH17B001P013-ready_chr21_pos_damage.tsv"         
 [5] "MH17T001P004-ready_chr21_pos_damage.tsv"         
 [6] "MH17T001P013-ready_chr21_pos_damage.tsv"         
 [7] "WES005PVFE_chr21_pos_damage.tsv"                 
 [8] "WES005PVFR_chr21_pos_damage.tsv"                 
 [9] "WES006JCFE_chr21_pos_damage.tsv"                 
[10] "WES006JCFR_chr21_pos_damage.tsv"                 
[11] "WES011RCFE_chr21_pos_damage.tsv"                 
[12] "WES011RCFR_chr21_pos_damage.tsv"                 
[13] "WES013PFFR-sort_chr21_pos_damage.tsv"            
```

```r
mut_list <- lapply(fnames, read_pos_damage)

mut_all <- dplyr::bind_rows(mut_list)
mut_all <- mut_all %>% 
  mutate(id2 = sub("^\\/(.*).*_chr21.*", "\\1", id),
         id2 = sub("-(sort|ready)", "", id2),
         id2 = sub("-germline", "-germ", id2))
```

__Example for type T-G__

```r
mut_all %>% 
  filter(type == "T_G") %>% 
  ggplot() +
  geom_point(aes(x = loc, y = rel), alpha = 1/2) +
  theme_bw() +
  facet_grid(id2 ~ read, scales = "free_y")
```

![](report_files/figure-html/pos_damage-1.png)<!-- -->

__All types shown in separate folder__

```r
for (t in type_clean) {
  g <- mut_all %>% 
    filter(type == t) %>% 
    ggplot() +
    geom_point(aes(x = loc, y = rel), alpha = 1/2) +
    theme_bw() +
    facet_grid(id2 ~ read, scales = "free_y") +
    ggtitle(paste0("Fraction of ", substr(t, 1, 1), " to ", substr(t, 3, 3), " variants on R1 and R2"))
  ggsave(plot = g, filename = paste0("figs/pos_dam_type_", t, ".png"), height = 20, width = 15)
}
```

#### Context Damage

```r
context_clean <- paste0(c("A", "T", "C", "G"), "_BASE_", rep(c("A", "T", "C", "G"), each = 4))
fnames <- list.files("~/Desktop/git/dna-damage-estimation/data/new",
                     pattern = "con_damage.tsv$", full.names = TRUE)
basename(fnames)
```

```
 [1] "MH17B001P004-germline-ready_chr21_con_damage.tsv"
 [2] "MH17B001P004-ready_chr21_con_damage.tsv"         
 [3] "MH17B001P013-germline-ready_chr21_con_damage.tsv"
 [4] "MH17B001P013-ready_chr21_con_damage.tsv"         
 [5] "MH17T001P004-ready_chr21_con_damage.tsv"         
 [6] "MH17T001P013-ready_chr21_con_damage.tsv"         
 [7] "WES005PVFE_chr21_con_damage.tsv"                 
 [8] "WES005PVFR_chr21_con_damage.tsv"                 
 [9] "WES006JCFE_chr21_con_damage.tsv"                 
[10] "WES006JCFR_chr21_con_damage.tsv"                 
[11] "WES011RCFE_chr21_con_damage.tsv"                 
[12] "WES011RCFR_chr21_con_damage.tsv"                 
[13] "WES013PFFR-sort_chr21_con_damage.tsv"            
```

```r
mut_list <- lapply(fnames, read_con_damage)

mut_all <- dplyr::bind_rows(mut_list)
mut_all <- mut_all %>% 
  mutate(id2 = sub("^\\/(.*).*_chr21.*", "\\1", id),
         id2 = sub("-(sort|ready)", "", id2),
         id2 = sub("-germline", "-germ", id2))
```

__Example for type T-G and context TBT__

```r
mut_all %>% 
  filter(type == "T_G", context == "T_BASE_T") %>% 
  ggplot() +
  geom_point(aes(x = loc, y = rel), alpha = 1/2) +
  theme_bw() +
  facet_grid(id2 ~ read, scales = "free_y")
```

![](report_files/figure-html/con_damage-1.png)<!-- -->

__All types shown in separate folder__

```r
for (t in type_clean) {
  for (con in context_clean) {
    
    g <- mut_all %>% 
      filter(type == t) %>% 
      ggplot() +
      geom_point(aes(x = loc, y = rel), alpha = 1/2) +
      theme_bw() +
      facet_grid(id2 ~ read, scales = "free_y") +
      ggtitle(paste0("Fraction of ", substr(t, 1, 1), " to ", substr(t, 3, 3), " variants on R1 and R2",
                     " in context ", con))
    ggsave(plot = g, filename = paste0("figs/con_dam_type_", t, "_cont_", con, ".png"), height = 20, width = 15)
  }
}
```
