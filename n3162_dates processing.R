# 한국 청동기시대 시기별 주거지 밀도 변동 양상 연구, R script
# KDE 공간분석 도면 작성을 위한 방사성탄소연대 처리 과정

# 1. 처리 대상 : 3162건 

install.packages("rcarbon")

library(rcarbon)

# data load
getwd()
setwd("F:/2026-02-13_KDE") # 임의 폴더 지정하여 작업
dates<-read.csv("dates_n3162.csv")# 기본 작업폴더에 해당 csv 파일(3162건의 탄소연대 목록)을 다운로드 후 작업 진행

#  보정 calibration
calDates <- calibrate(
  x = dates$bp,
  errors = dates$er,
  calCurves = 'intcal20'
)

# medianBP, dates에 추가.
calSummary <- summary(calDates)

dates$cal_median <- calSummary$MedianBP

### 여기까지 median 값 구해서 추가함.

### 1) 연대 구간 설정
# 1500–300 BCE, 200년 단위, 추가로 1700-1500, 300-100 BCE 구간도 설정. 확률 계산도 시행.
bce_seq <- seq(1700, 100, by = -200)

time_bins <- data.frame(
  start_BCE = bce_seq[-length(bce_seq)],
  end_BCE   = bce_seq[-1]
)

# BCE → cal BP (올바른 변환)
time_bins$start_BP <- 1950 + time_bins$start_BCE
time_bins$end_BP   <- 1950 + time_bins$end_BCE


### 2) 결과 저장 행렬
n_dates <- length(calDates$grids)
n_bins  <- nrow(time_bins)

bin_matrix <- matrix(0, nrow = n_dates, ncol = n_bins)


### 3) 연대치별 확률 값 계산

for (i in seq_len(n_dates)) {
  
  grid_i <- calDates$grids[[i]]
  
  for (j in seq_len(n_bins)) {
    
    idx <- grid_i$calBP <= time_bins$start_BP[j] &
      grid_i$calBP >  time_bins$end_BP[j]
    
    if (any(idx)) {
      bin_matrix[i, j] <- sum(grid_i$PrDens[idx])
    }
  }
}


### 열 이름 지정. 
colnames(bin_matrix) <- paste0(
  time_bins$start_BCE, "_",
  time_bins$end_BCE, "_BCE"
)


dates <- cbind(dates, bin_matrix)

save(dates, file = "n3162_20260213dates.R")
