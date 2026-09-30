# ============================================================
# LR04 benthic d18O stack with approximate Vostok temperature
# and major 100 kyr / 41 kyr climate cyclicity
#
# Data:
#   Lisiecki & Raymo (2005) LR04 benthic d18O stack
#   Petit et al. (1999) Vostok temperature record
#
# X-axis:
#   Logarithmic years before present
#   Starts at 10^3 years
#   Ends at the actual maximum age of the LR04 dataset
#   (~5.3 million years)
#
# Cyclicity annotations:
#   ~100 kyr cycle: recent to ~1 Ma
#   ~41 kyr cycle:  ~1 Ma to ~2.5 Ma
#
# ============================================================


# ------------------------------------------------------------
# LOAD REQUIRED PACKAGE
# ------------------------------------------------------------

library(zoo)


# ============================================================
# 1. READ LR04 DATA
# ============================================================

lr04_path <- "Global_stack_d18O.tab"


# Read raw file to identify the end of the metadata/header block

raw_lines <- readLines(
  lr04_path,
  encoding = "UTF-8"
)


# Find the line containing the end of the comment block

header_end <- grep(
  "^\\*/",
  raw_lines
)


# Determine where the actual data begin

data_start <- if (length(header_end) > 0) {
  header_end[1] + 1
} else {
  1
}


# Read the LR04 table

lr04 <- read.table(
  lr04_path,
  skip = data_start - 1,
  header = TRUE,
  sep = "\t",
  fill = TRUE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)


# Rename the first two columns

names(lr04)[1:2] <- c(
  "age_ka",
  "d18O"
)


# Remove rows with missing values

lr04 <- lr04[
  !is.na(lr04$age_ka) &
    !is.na(lr04$d18O),
]


# Sort from youngest to oldest

lr04 <- lr04[
  order(lr04$age_ka),
]


# Convert age from ka to years

lr04$age_yr <- lr04$age_ka * 1000



# ============================================================
# 2. READ VOSTOK TEMPERATURE DATA
# ============================================================

vostok_path <- "Vostok_Deut_dT.tab"


# Read raw Vostok file

raw_lines_v <- readLines(
  vostok_path,
  encoding = "UTF-8"
)


# Find the end of the metadata/header block

header_end_v <- grep(
  "^\\*/",
  raw_lines_v
)


# Determine where data begin

data_start_v <- if (length(header_end_v) > 0) {
  header_end_v[1] + 1
} else {
  1
}


# Read Vostok table

vostok <- read.table(
  vostok_path,
  skip = data_start_v - 1,
  header = TRUE,
  sep = "\t",
  fill = TRUE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)


# Rename columns
#
# Assumes:
#   column 2 = age (ka)
#   column 3 = deuterium (dD)
#   column 4 = temperature anomaly (deltaT)

names(vostok)[2:4] <- c(
  "age_ka",
  "dD",
  "deltaT"
)


# Remove missing values

vostok <- vostok[
  !is.na(vostok$age_ka) &
    !is.na(vostok$deltaT),
]


# Sort youngest to oldest

vostok <- vostok[
  order(vostok$age_ka),
]


# Convert age from ka to years

vostok$age_yr <- vostok$age_ka * 1000



# ============================================================
# 3. CREATE A COMMON TIME GRID FOR LR04 AND VOSTOK
# ============================================================

# Use 1 kyr intervals

dt <- 1000


# Determine the maximum overlapping age between datasets

overlap_max <- min(
  max(vostok$age_yr, na.rm = TRUE),
  max(lr04$age_yr, na.rm = TRUE)
)


# Common age grid

age_overlap <- seq(
  from = 0,
  to = overlap_max,
  by = dt
)


# ------------------------------------------------------------
# Interpolate LR04 d18O onto common grid
# ------------------------------------------------------------

lr04_overlap <- na.approx(
  lr04$d18O,
  x = lr04$age_yr,
  xout = age_overlap,
  rule = 2
)


# ------------------------------------------------------------
# Interpolate Vostok temperature onto common grid
# ------------------------------------------------------------

vostok_overlap <- na.approx(
  vostok$deltaT,
  x = vostok$age_yr,
  xout = age_overlap,
  rule = 2
)



# ============================================================
# 4. CALIBRATE LR04 d18O TO VOSTOK TEMPERATURE
# ============================================================

# Fit a linear relationship over the period where both datasets
# overlap.
#
# Relationship:
#
# Vostok Delta T = intercept + slope × LR04 d18O

calibration <- lm(
  vostok_overlap ~ lr04_overlap
)


# Extract regression coefficients

intercept <- coef(calibration)[1]

slope <- coef(calibration)[2]


# Optional: print calibration results

print(summary(calibration))

cat(
  "\nTemperature calibration equation:\n",
  "Vostok DeltaT = ",
  round(intercept, 4),
  " + (",
  round(slope, 4),
  " × LR04 d18O )\n\n",
  sep = ""
)



# ============================================================
# 5. CREATE AN EVENLY SPACED FULL LR04 TIME SERIES
# ============================================================

# Create a 1 kyr grid spanning the full LR04 dataset

age_grid <- seq(
  from = min(lr04$age_yr),
  to = max(lr04$age_yr),
  by = dt
)


# Interpolate LR04 onto this grid

d18O_interp <- na.approx(
  lr04$d18O,
  x = lr04$age_yr,
  xout = age_grid,
  rule = 2
)



# ============================================================
# 6. DEFINE THE LOGARITHMIC X-AXIS
# ============================================================

# The LR04 dataset extends to approximately 5.3 million years.
#
# We start at 10^3 years because zero cannot be represented
# on a logarithmic axis.

xlim_years <- c(
  1e3,
  max(lr04$age_yr)
)



# ------------------------------------------------------------
# Function for major logarithmic axis ticks
# ------------------------------------------------------------

log_axis_x10 <- function(
    side = 1
) {
  
  # Major powers of 10 within the LR04 range
  
  ticks <- c(
    1e3,
    1e4,
    1e5,
    1e6
  )
  
  
  labels <- expression(
    10^3,
    10^4,
    10^5,
    10^6
  )
  
  
  axis(
    side = side,
    at = ticks,
    labels = labels,
    las = 1
  )
  
  
  # Add final LR04 maximum-age tick
  
  lr04_max <- max(lr04$age_yr)
  
  axis(
    side = side,
    at = lr04_max,
    labels = FALSE
  )
}



# ============================================================
# 7. CREATE OUTPUT FIGURE
# ============================================================

png(
  "LR04_5Myr_Vostok_logscale_corrected.png",
  width = 3000,
  height = 1800,
  res = 300
)



# ------------------------------------------------------------
# Set plotting margins
# ------------------------------------------------------------

par(
  mar = c(
    5,
    5.5,
    2.5,
    6
  )
)



# ============================================================
# 8. PLOT FULL LR04 CLIMATE RECORD
# ============================================================

plot(
  age_grid,
  d18O_interp,
  
  type = "l",
  
  col = "royalblue4",
  
  lwd = 0.7,
  
  log = "x",
  
  xlim = xlim_years,
  
  xlab = "Years before present",
  
  ylab = expression(
    delta^18 * O ~
      "(LR04 benthic stack)"
  ),
  
  xaxt = "n"
)



# Add the logarithmic x-axis

log_axis_x10()



# ============================================================
# 9. ADD APPROXIMATE VOSTOK TEMPERATURE REFERENCE
# ============================================================

# Calculate the LR04 d18O value corresponding to Vostok
# DeltaT = 0
#
# 0 = intercept + slope × d18O
#
# d18O = -intercept / slope

d18O_zero_temp <- -intercept / slope


# Add dashed reference line

abline(
  h = d18O_zero_temp,
  lty = 2,
  col = "grey60"
)



# ============================================================
# 10. ADD RIGHT-HAND VOSTOK TEMPERATURE AXIS
# ============================================================

# Choose the temperature values to display

temp_ticks <- seq(
  -8,
  2,
  by = 2
)


# Convert temperature values to corresponding LR04 d18O values
#
# DeltaT = intercept + slope × d18O
#
# Therefore:
#
# d18O = (DeltaT - intercept) / slope

d18O_temp_ticks <- (
  temp_ticks - intercept
) / slope


# Add right-hand axis

axis(
  side = 4,
  
  at = d18O_temp_ticks,
  
  labels = temp_ticks,
  
  las = 1
)


# Add right-hand axis title

mtext(
  expression(
    paste(
      "Equivalent Vostok ",
      Delta,
      "T (",
      degree,
      "C)"
    )
  ),
  
  side = 4,
  
  line = 4
)



# ============================================================
# 11. ADD CYCLICITY ANNOTATIONS
# ============================================================

usr <- par("usr")


# ------------------------------------------------------------
# 100 kyr cycle: approximately 0–1 Ma
#
# On a log axis, 0 cannot be plotted, so the line begins at
# 10^3 years, representing the near-present end of the record.
# ------------------------------------------------------------

y_100 <- usr[4] - 0.10 * diff(usr[3:4])

segments(
  x0 = 1e3,
  y0 = y_100,
  
  x1 = 1e6,
  y1 = y_100,
  
  col = "darkgreen",
  lwd = 3
)


text(
  x = sqrt(1e3 * 1e6),
  y = y_100 + 0.08 * diff(usr[3:4]),
  
  labels = "100 kyr cycle",
  
  col = "darkgreen",
  cex = 1.1
)



# ------------------------------------------------------------
# 41 kyr cycle: approximately 1–2.5 Ma
#
# Positioned slightly lower than the 100 kyr annotation so
# that the two cycle bars remain visually distinct while
# retaining their correct chronological ranges.
# ------------------------------------------------------------

y_41 <- usr[4] - 0.16 * diff(usr[3:4])

segments(
  x0 = 1e6,
  y0 = y_41,
  
  x1 = 2.5e6,
  y1 = y_41,
  
  col = "darkgreen",
  lwd = 3
)


text(
  x = sqrt(1e6 * 2.5e6),
  y = y_41 + 0.08 * diff(usr[3:4]),
  
  labels = "41 kyr cycle",
  
  col = "darkgreen",
  cex = 1.1
)
# ============================================================
# 12. ADD LEGEND
# ============================================================

legend(
  "topleft",
  
  legend = expression(
    delta^18 * O ~
      "(LR04 benthic stack)"
  ),
  
  col = "royalblue4",
  
  lwd = 1,
  
  bty = "n"
)



# ============================================================
# 13. CLOSE GRAPHICS DEVICE
# ============================================================

dev.off()

