library(dplyr)


### load data
## strain AMM (median tel length =555 bp)
AMM_FvR_proportions <- read_delim("Documents/DUBii/Projet/results/Telom_read_prop/AMM_FvR_proportions.csv", 
                                  "\t", escape_double = FALSE, trim_ws = TRUE)

## strain CEL (median tel length = 210 bp)
CEL_FvR_proportions <- read_delim("Documents/DUBii/Projet/results/Telom_read_prop/CEL_FvR_proportions.csv", 
                                  "\t", escape_double = FALSE, trim_ws = TRUE)


### Filter  and plot  chromosome start coordinates
## AMM 1-50 bp
AMM_coord_1_50 <- filter(AMM_FvR_proportions, bp < 50)

ggplot() +
  geom_line(data = coord_50, aes(bp, forward_prop, colour = chr)) +
  geom_line(data = coord_50, aes(bp, reverse_prop, colour = chr)) +
  ggtitle("AMM 1_50 bp")

## CEL 1-50 bp
CEL_coord_1_50 <- filter(CEL_FvR_proportions, bp < 50)

ggplot() +
  geom_line(data = coord_50, aes(bp, forward_prop, colour = chr)) +
  geom_line(data = coord_50, aes(bp, reverse_prop, colour = chr)) +
  ggtitle("CEL 1_50 bp")



## AMM 50-150 bp
coord_50_150 <- filter(AMM_FvR_proportions, bp > 50, bp <150)

ggplot() +
  geom_line(data = coord_50_150, aes(bp, forward_prop, colour = chr)) +
  geom_line(data = coord_50_150, aes(bp, reverse_prop, colour = chr)) +
  ggtitle("AMM 50-150 bp")

## CEL 50-150 bp
coord_50_150 <- filter(CEL_FvR_proportions, bp > 50, bp <150)

ggplot() +
  geom_line(data = coord_50_150, aes(bp, forward_prop, colour = chr)) +
  geom_line(data = coord_50_150, aes(bp, reverse_prop, colour = chr)) +
  ggtitle("CEL 50-150 bp")

## AMM 1-1000 bp
AMM_coord_1_1000 <- filter(AMM_FvR_proportions, bp < 1000)

ggplot() +
  geom_line(data = AMM_coord_1_1000, aes(bp, forward_prop, colour = chr)) +
  geom_line(data = AMM_coord_1_1000, aes(bp, reverse_prop, colour = chr)) +
  ggtitle("AMM 1_1000 bp")

## CEL 1-1000 bp
CEL_coord_1_1000 <- filter(CEL_FvR_proportions, bp < 1000)

ggplot() +
  geom_line(data = CEL_coord_1_1000, aes(bp, forward_prop, colour = chr)) +
  geom_line(data = CEL_coord_1_1000, aes(bp, reverse_prop, colour = chr)) +
  ggtitle("CEL 1_1000 bp")


## AMM 1000-2000 bp
AMM_coord_1000_2000 <- filter(AMM_FvR_proportions, bp > 1000, bp < 2000)

ggplot() +
  geom_line(data = AMM_coord_1000_2000, aes(bp, forward_prop, colour = chr)) +
  geom_line(data = AMM_coord_1000_2000, aes(bp, reverse_prop, colour = chr)) +
  ggtitle("AMM 1000_2000 bp")

## CEL 1000-2000 bp
CEL_coord_1000_2000 <- filter(CEL_FvR_proportions, bp > 1000, bp < 2000)

ggplot() +
  geom_line(data = CEL_coord_1000_2000, aes(bp, forward_prop, colour = chr)) +
  geom_line(data = CEL_coord_1000_2000, aes(bp, reverse_prop, colour = chr)) +
  ggtitle("CEL 1000_2000 bp")

### Filter and plot end coordinates
AMM_end_1000 <- group_by(AMM_FvR_proportions, chr) %>% top_n(1000, bp)
CEL_end_1000 <- group_by(CEL_FvR_proportions, chr) %>% top_n(1000, bp)

## add an index column per chromosme
AMM_end_1000 <- AMM_end_1000 %>%
  group_by(chr) %>%
  mutate(index = row_number()) %>%
  ungroup()

CEL_end_1000 <- CEL_end_1000 %>%
  group_by(chr) %>%
  mutate(index = row_number()) %>%
  ungroup()

## AMM end-1000
ggplot() +
  geom_line(data = AMM_end_1000, aes(index, forward_prop, colour = chr)) +
  geom_line(data = AMM_end_1000, aes(index, reverse_prop, colour = chr)) +
  ggtitle("AMM end-1000 bp")

## CEL end-1000
ggplot() +
  geom_line(data = CEL_end_1000, aes(index, forward_prop, colour = chr)) +
  geom_line(data = CEL_end_1000, aes(index, reverse_prop, colour = chr)) +
  ggtitle("CEL end-1000 bp")


