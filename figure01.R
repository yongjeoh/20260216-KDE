library(rcarbon)
library(ggplot2)
library(dplyr)
library(showtext)

set.seed(1234)

# ---------------------------
# 0. 폰트 설정
# ---------------------------
font_add_google("Noto Sans KR", "noto")
showtext_auto()

theme_set(
  theme_minimal(base_family = "noto") +
    theme(
      text = element_text(size = 14),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    )
)

# ---------------------------
# 1. 가상 데이터 생성
# ---------------------------
n <- 100
ages <- runif(n, 2300, 3000)
errors <- runif(n, 30, 60)

cal <- calibrate(x = ages, errors = errors, calCurves = 'intcal20')

# ---------------------------
# 2. 시간 구간 설정
# ---------------------------
breaks_bce <- seq(1300, 300, by = -200)
breaks_bp <- sort(breaks_bce + 1950)

labels <- c(
  "1300–1100",
  "1100–900",
  "900–700",
  "700–500",
  "500–300"
)

# ---------------------------
# 3. median 방식
# ---------------------------
meds <- medCal(cal)

if (is.data.frame(meds)) {
  med_bp <- meds$MedianBP
} else {
  med_bp <- as.numeric(meds)
}

median_bin <- cut(
  med_bp,
  breaks = breaks_bp,
  labels = labels,
  include.lowest = TRUE
)

median_df <- data.frame(bin = median_bin) |>
  count(bin)

median_df$bin <- factor(median_df$bin, levels = labels)

# ---------------------------
# 4. probability 방식
# ---------------------------
prob_df <- data.frame(
  bin = labels,
  prob = 0
)

for (i in 1:length(cal$grids)) {
  g <- cal$grids[[i]]
  df <- data.frame(calBP = g[,1], PrDens = g[,2])
  
  for (j in 1:(length(breaks_bp)-1)) {
    idx <- df$calBP >= breaks_bp[j] & df$calBP < breaks_bp[j+1]
    prob_df$prob[j] <- prob_df$prob[j] + sum(df$PrDens[idx])
  }
}

prob_df$bin <- factor(prob_df$bin, levels = labels)

# ---------------------------
# 5. SPD (runm smoothing 포함)
# ---------------------------
spd_raw <- spd(cal, timeRange = c(3200, 2200))
spd_smooth <- spd(cal, timeRange = c(3200, 2200), runm = 200)

spd_raw_df <- spd_raw$grid |>
  mutate(year = calBP - 1950)

spd_smooth_df <- spd_smooth$grid |>
  mutate(year = calBP - 1950)

# ---------------------------
# 6. 그림 생성
# ---------------------------

# (A) median
p1 <- ggplot(median_df, aes(x = bin, y = n)) +
  geom_col() +
  labs(
    title = "(1) 중간값 기준 연대치 분포도 (200년 단위)",
    x = "Calendar year (BCE)",
    y = "Count"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# (B) probability
p2 <- ggplot(prob_df, aes(x = bin, y = prob)) +
  geom_col() +
  labs(
    title = "(2) 확률분포 기준 연대치 분포도 (200년 단위)",
    x = "Calendar year (BCE)",
    y = "Summed probability"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# (C) SPD
p3 <- ggplot() +
  geom_line(
    data = spd_raw_df,
    aes(x = year, y = PrDens),
    linewidth = 0.6,
    alpha = 0.4
  ) +
  geom_line(
    data = spd_smooth_df,
    aes(x = year, y = PrDens),
    linewidth = 1.2,
    color = "red",
    linetype = "dashed"
  ) +
  scale_x_reverse(
    limits = c(1300, 300),
    breaks = c(1300, 1100, 900, 700, 500, 300)
  ) +
  labs(
    title = "(3) 합산확률분포(SPD): 확률분포 전체 사용",
    x = "Calendar year (BCE)",
    y = "Probability density"
  )

# ---------------------------
# 7. 출력
# ---------------------------
p1
p2
p3