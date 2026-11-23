setwd("/Users/mingju/Desktop/600/capstone")
library(dplyr)
library(tidyr)
library(tibble)

load("weather.rda")

weather_long <- weather %>%
  rownames_to_column(var = "date") %>%
  mutate(
    across(
      matches("^[A-Z]{3}_PrecipitationIn$"),
      ~ {
        x <- as.character(.x)
        x[x %in% c("T", "")] <- NA
        as.numeric(x)
      }
    )
  ) %>%
  pivot_longer(
    cols = matches("^[A-Z]{3}_"),
    names_to = c("airport", ".value"),
    names_pattern = "^([^_]+)_(.*)$"
  ) %>%
  mutate(date = as.Date(date))

# Check the result
head(weather_long)

### y vs time by airports

# Select the three temperature variables
temperature_long <- weather_long %>%
  select(date, airport,
         Max.TemperatureF,
         Mean.TemperatureF,
         Min.TemperatureF) %>%
  pivot_longer(
    cols = c(
      Max.TemperatureF,
      Mean.TemperatureF,
      Min.TemperatureF
    ),
    names_to = "temperature_type",
    values_to = "temperature"
  ) %>%
  mutate(
    temperature_type = recode(
      temperature_type,
      Max.TemperatureF = "Daily Maximum Temperature",
      Mean.TemperatureF = "Daily Mean Temperature",
      Min.TemperatureF = "Daily Minimum Temperature"
    )
  )

# Plot
ggplot(
  temperature_long,
  aes(
    x = date,
    y = temperature,
    color = airport,
    group = airport
  )
) +
  geom_line(alpha = 0.7, linewidth = 0.4) +
  facet_wrap(
    ~ temperature_type,
    nrow = 1,
    scales = "free_y"
  ) +
  labs(
    title = "Daily Temperature Trends Across Airports",
    x = "Date",
    y = "Temperature (°F)",
    color = "Airport"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5)
  )+
  theme(
    legend.position = "none"
  )