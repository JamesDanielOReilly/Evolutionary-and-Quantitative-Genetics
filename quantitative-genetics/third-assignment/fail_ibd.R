## for related individuals: exclude the one with the higher F_miss
rel <- read.table('gwa_raw.genome', header = T)
head(rel) %>% View
rel_0_2 <- rel %>% filter(PI_HAT > 0.2 ) %>% select(IID1, IID2)
rel_0_2 <- rel_0_2 %>% left_join(imiss[, c('IID', 'F_MISS')], by = c('IID1' = 'IID')) %>% rename(F_MISS1 = F_MISS)
rel_0_2 <- rel_0_2 %>% left_join(imiss[, c('IID', 'F_MISS')], by = c('IID2' = 'IID')) %>% rename(F_MISS2 = F_MISS)

# exclude the individual with the higher F_MISS
rel_0_2 <- rel_0_2 %>% mutate(exclude = ifelse(F_MISS1 > F_MISS2, IID1, IID2))
excl <- rel_0_2 %>% select(exclude) %>% rename(IID = exclude)
excl <- excl %>% left_join(imiss[, c('IID', 'FID')]) %>% select(FID, IID)
excl <- unique(excl)
write.table(as.matrix(excl), file = 'fail-IBD-qc.txt', sep="\t", col.names = F, row.names = F)
