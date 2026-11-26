# Forecast days 16-20 of every month for the validation year.
#
# Usage (run from the Taehyoung folder, after rendering analysis.qmd):
#   Rscript predict.R validation.rda forecasts.csv
#
# Input:  an .rda file with one data frame in the same format as weather.rda
#         (one row per day, dates in the row names or in a `date` column,
#         columns like DTW_Max.TemperatureF), holding days 1-15 of each month.
# Output: a CSV with 2,700 forecasts (15 airports x 3 temperatures x 12 months x 5 days).

suppressMessages({
  library(dplyr)
  library(tidyr)
})

args   <- commandArgs(trailingOnly = TRUE)
input  <- if (length(args) >= 1) args[1] else "validation.rda"
output <- if (length(args) >= 2) args[2] else "forecasts.csv"

final <- readRDS("final_models.rds")
omega <- final$omega

# ---------------------------------------------------------------------------
# 1. Read and clean (the general rules of preprocess.qmd)
# ---------------------------------------------------------------------------

env <- new.env()
load(input, envir = env)
raw <- Filter(is.data.frame, mget(ls(env), envir = env))[[1]]

clean_weather <- function(raw) {
  dates <- if ("date" %in% names(raw)) as.Date(raw$date) else as.Date(rownames(raw))
  raw   <- raw[, setdiff(names(raw), "date"), drop = FALSE]
  raw[] <- lapply(raw, \(x) if (is.factor(x)) as.character(x) else x)

  w <- raw |>
    mutate(across(ends_with("_PrecipitationIn"), as.character),
           across(ends_with("_Events"), as.character),
           date = dates) |>
    pivot_longer(-date, names_to = c("airport", ".value"),
                 names_pattern = "^([A-Z]{3})_(.*)$") |>
    rename(any_of(c(Max.DewPointF = "Max.Dew.PointF",
                    Mean.DewPointF = "MeanDew.PointF",
                    Min.DewPointF = "Min.DewpointF")))

  # Temperature: values above 1000 are a temperature with other digits stuck on.
  # Keep the leading digits that agree with the mean; otherwise set to NA.
  leading <- \(x) if (!is.na(x) && x > 1000) as.integer(substr(rep(as.character(x), 3), 1, 1:3)) else x
  for (r in which(w$Max.TemperatureF > 1000 | w$Min.TemperatureF > 1000)) {
    grid <- expand.grid(mx = leading(w$Max.TemperatureF[r]), mn = leading(w$Min.TemperatureF[r]))
    ok   <- abs((grid$mx + grid$mn) / 2 - w$Mean.TemperatureF[r]) <= 1 & grid$mx >= grid$mn
    ok[is.na(ok)] <- FALSE
    w$Max.TemperatureF[r] <- if (any(ok)) grid$mx[ok][1] else NA
    w$Min.TemperatureF[r] <- if (any(ok)) grid$mn[ok][1] else NA
  }
  temp_bad <- with(w, Max.TemperatureF > 120 | Min.TemperatureF > 100 | Min.TemperatureF < -50 |
                      Min.TemperatureF > Mean.TemperatureF | Mean.TemperatureF > Max.TemperatureF |
                      Max.TemperatureF - Min.TemperatureF > 65)
  dew_bad  <- with(w, Min.DewPointF < -60 | Max.DewPointF > Max.TemperatureF + 5)
  hum_bad  <- with(w, Min.Humidity == 0 | Min.Humidity > Mean.Humidity | Mean.Humidity > Max.Humidity)
  wind_bad <- with(w, Mean.Wind.SpeedMPH > Max.Wind.SpeedMPH)

  w |>
    mutate(
      across(c(Max.TemperatureF, Mean.TemperatureF, Min.TemperatureF), \(x) replace(x, temp_bad %in% TRUE, NA)),
      Mean.DewPointF     = replace(Mean.DewPointF, dew_bad %in% TRUE, NA),
      Mean.Humidity      = replace(Mean.Humidity, hum_bad %in% TRUE, NA),
      Mean.Wind.SpeedMPH = replace(Mean.Wind.SpeedMPH, wind_bad %in% TRUE, NA),
      WindDirDegrees     = na_if(WindDirDegrees, -1),
      CloudCover         = replace(CloudCover, !(CloudCover %in% 0:8), NA),
      Mean.VisibilityMiles = pmin(Mean.VisibilityMiles, 10),
      PrecipitationIn    = na_if(PrecipitationIn, ""),
      PrecipitationIn    = suppressWarnings(as.numeric(replace(PrecipitationIn, PrecipitationIn == "T", "0"))),
      PrecipitationIn    = replace(PrecipitationIn, PrecipitationIn >= 999, NA),
      rain    = as.numeric(grepl("Rain", Events)),
      snow    = as.numeric(grepl("Snow", Events)),
      fog     = as.numeric(grepl("Fog", Events)),
      thunder = as.numeric(grepl("Thunderstorm", Events)),
      wind_u  = Mean.Wind.SpeedMPH * sin(WindDirDegrees * pi / 180),
      wind_v  = Mean.Wind.SpeedMPH * cos(WindDirDegrees * pi / 180)
    )
}

w <- clean_weather(raw)

# ---------------------------------------------------------------------------
# 2. Daily matrices (date x airport) on a complete calendar of days 1-15
# ---------------------------------------------------------------------------

year <- as.integer(names(which.max(table(format(w$date, "%Y")))))
days <- as.Date(unlist(lapply(1:12, \(m) as.character(as.Date(sprintf("%d-%02d-%02d", year, m, 1:15))))))
airports <- final$airports

wide <- function(variable) {
  m <- matrix(NA_real_, length(days), length(airports), dimnames = list(as.character(days), airports))
  v <- w |> filter(date %in% days, airport %in% airports)
  m[cbind(match(as.character(v$date), rownames(m)), match(v$airport, airports))] <- as.numeric(v[[variable]])
  m
}
season_x <- function(d) cbind(1, sin(omega * d), cos(omega * d), sin(2 * omega * d), cos(2 * omega * d))
day_of_year <- as.integer(format(days, "%j"))
normal_of <- function(coefs) season_x(day_of_year) %*% coefs       # date x airport
departure <- function(raw, coefs) raw - normal_of(coefs)[, colnames(raw), drop = FALSE]

# Fill missing temperatures at the target airports from the days around them (within each month)
interpolate_month <- function(x, normal) {
  for (m in 1:12) {
    i  <- which(as.integer(format(days, "%m")) == m)
    ok <- !is.na(x[i])
    x[i] <- if (sum(ok) >= 2) approx(which(ok), x[i][ok], xout = seq_along(i), rule = 2)$y
            else if (sum(ok) == 1) ifelse(ok, x[i], x[i][ok])
            else normal[i]                                         # nothing observed: use normal
  }
  x
}

temps <- list(Max = wide("Max.TemperatureF"), Mean = wide("Mean.TemperatureF"), Min = wide("Min.TemperatureF"))
for (tt in names(temps)) {
  nrm <- normal_of(final$normals_local[[paste0("t_", tt)]])
  for (a in final$target_airports) temps[[tt]][, a] <- interpolate_month(temps[[tt]][, a], nrm[, a])
}

# Local weather departures (missing -> 0, i.e. normal)
fill0 <- \(m) replace(m, is.na(m), 0)
local_dep <- list(
  t_Max = departure(temps$Max, final$normals_local$t_Max),
  t_Mean = departure(temps$Mean, final$normals_local$t_Mean),
  t_Min = departure(temps$Min, final$normals_local$t_Min)
)
for (v in c("pressure", "dew_point", "humidity", "wind_u", "wind_v", "cloud", "visibility", "precip")) {
  source_col <- c(pressure = "Mean.Sea.Level.PressureIn", dew_point = "Mean.DewPointF",
                  humidity = "Mean.Humidity", wind_u = "wind_u", wind_v = "wind_v",
                  cloud = "CloudCover", visibility = "Mean.VisibilityMiles", precip = "PrecipitationIn")[[v]]
  local_dep[[v]] <- fill0(departure(wide(source_col), final$normals_local[[v]]))
}
press_dep_raw <- departure(wide("Mean.Sea.Level.PressureIn"), final$normals_local$pressure)
flags <- lapply(c(rain = "rain", snow = "snow", fog = "fog", thunder = "thunder"), \(v) fill0(wide(v)))

# Other airports: departures of mean temperature and pressure, and their 1-day change
dep_temp  <- departure(wide("Mean.TemperatureF"), final$normals_temp)
raw_press <- wide("Mean.Sea.Level.PressureIn")
raw_press[, !final$has_press] <- NA
dep_press <- departure(raw_press, final$normals_press)
change    <- \(m) rbind(NA, m[-1, , drop = FALSE] - m[-nrow(m), , drop = FALSE])
chg_temp  <- fill0(change(dep_temp));  dep_temp  <- fill0(dep_temp)
chg_press <- fill0(change(dep_press)); dep_press <- fill0(dep_press)

# ---------------------------------------------------------------------------
# 3. Inputs for each start day (the 15th) and forecasts
# ---------------------------------------------------------------------------

forecast_one <- function(model, m) {
  a  <- model$airport
  tt <- model$target
  h  <- model$lead
  i  <- which(days == as.Date(sprintf("%d-%02d-15", year, m)))   # start day row
  x  <- c()
  for (k in 0:14) x[paste0("lag", k)] <- temps[[tt]][i - k, a]
  other <- switch(tt, Max = c("t_Mean", "t_Min"), Mean = c("t_Max", "t_Min"), Min = c("t_Max", "t_Mean"))
  for (day in c(0, 1)) {
    sfx <- if (day == 0) "" else "_d1"
    r   <- i - day
    x[paste0("other_1", sfx)] <- local_dep[[other[1]]][r, a]
    x[paste0("other_2", sfx)] <- local_dep[[other[2]]][r, a]
    for (v in c("pressure", "dew_point", "humidity", "wind_u", "wind_v", "cloud", "visibility", "precip"))
      x[paste0(v, sfx)] <- local_dep[[v]][r, a]
    for (v in names(flags)) x[paste0(v, sfx)] <- flags[[v]][r, a]
    # Pressure change (raw); if missing, use the normal change (departure change 0)
    pc <- press_dep_raw[r, a] - press_dep_raw[r - 1, a]
    x[paste0("pressure_change", sfx)] <- (if (is.na(pc)) 0 else pc) +
      normal_of(final$normals_local$pressure)[r, a] - normal_of(final$normals_local$pressure)[r - 1, a]
  }
  d <- as.integer(format(days[i] + h, "%j"))
  x[c("sin1", "cos1", "sin2", "cos2")] <- c(sin(omega * d), cos(omega * d), sin(2 * omega * d), cos(2 * omega * d))
  others <- setdiff(airports, a)
  press  <- setdiff(airports[final$has_press], a)
  x[paste0("t_", others)]  <- dep_temp[i, others]
  x[paste0("tc_", others)] <- chg_temp[i, others]
  x[paste0("p_", press)]   <- dep_press[i, press]
  x[paste0("pc_", press)]  <- chg_press[i, press]

  b    <- model$coef
  vars <- setdiff(names(b), "(Intercept)")
  stopifnot(all(vars %in% names(x)))
  model$mu + model$s * (b[["(Intercept)"]] + sum(b[vars] * x[vars]))
}

forecasts <- bind_rows(lapply(final$models, \(model) {
  tibble(airport = model$airport, temperature = model$target, lead = model$lead, month = 1:12,
         forecast = sapply(1:12, \(m) forecast_one(model, m)))
}))

# Extra rules chosen in Step 6 of analysis.qmd
wide_fc <- forecasts |> pivot_wider(names_from = temperature, values_from = forecast)
if (final$use_mean_from_extremes) wide_fc$Mean <- (wide_fc$Max + wide_fc$Min) / 2
if (final$use_order_fix) {
  bad <- !(wide_fc$Min <= wide_fc$Mean & wide_fc$Mean <= wide_fc$Max)
  if (any(bad)) {
    fixed <- t(apply(as.matrix(wide_fc[bad, c("Min", "Mean", "Max")]), 1, \(z) isoreg(z)$yf))
    wide_fc[bad, c("Min", "Mean", "Max")] <- as.data.frame(fixed)
  }
}

out <- wide_fc |>
  pivot_longer(c(Max, Mean, Min), names_to = "temperature", values_to = "forecast") |>
  mutate(year = year, day = 15L + lead, date = as.Date(sprintf("%d-%02d-%02d", year, month, day))) |>
  select(date, year, month, day, airport, temperature, forecast) |>
  arrange(date, airport, temperature)

stopifnot(nrow(out) == 2700, !anyNA(out$forecast))
write.csv(out, output, row.names = FALSE)
cat("Wrote", nrow(out), "forecasts for", year, "to", output, "\n")
