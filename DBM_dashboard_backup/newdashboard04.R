# =========================================================================
# DBM Dashboard – Near-term Pest Pressure (14 Days ONLY)
# =========================================================================




library(shiny)
library(data.table)
library(lubridate)
library(weatherOz)
library(leaflet)
library(ggplot2)
library(plotly)
library(sf)


# =========================================================================
# =========================================================================
#
library(terra)

sdm <- rast("pred_southwest_rf_AUC2023.tif")

cols <- colorRampPalette(
  c("forestgreen", "orange", "#9B242c")
)(100)

dir.create("www", showWarnings = FALSE)

png(
  "www/sdm_context.png",
  width  = 3000,
  height = 2200,
  res    = 300
)

plot(
  sdm,
  col  = cols,
  main = "Long-term habitat suitability (SDM)",
  axes = FALSE,
  box  = FALSE
)
dev.off()


# =========================================================================
# 0. SPATIAL DATA
# =========================================================================

path_lga  <- "~/Documents/FSI/NGDSI_Project_DBM/modelling/maps/WA/LGA.shp"
path_swar <- "~/Documents/FSI/NGDSI_Project_DBM/modelling/maps/WA/South_west_agricultural_region_DPIRD_008_WA_GDA2020_Public_Shapefile/South_west_agricultural_region_DPIRD_008.shp"

WA_TOWNS  <- NULL
SWAR_POLY <- NULL

if (file.exists(path_lga) && file.exists(path_swar)) {
  
  lga_poly  <- st_read(path_lga, quiet = TRUE)
  swar_poly <- st_read(path_swar, quiet = TRUE)
  
  if (is.na(st_crs(lga_poly)))  st_crs(lga_poly)  <- 4283
  if (is.na(st_crs(swar_poly))) st_crs(swar_poly) <- 4283
  
  lga_poly  <- st_transform(lga_poly, 4326)
  swar_poly <- st_transform(swar_poly, 4326)
  
  pts <- suppressWarnings(st_centroid(lga_poly))
  pts <- st_filter(pts, swar_poly)
  xy  <- st_coordinates(pts)
  
  WA_TOWNS <- data.table(
    Town = gsub(" \\(.*\\)", "", pts$SHIRE_NAME),
    Lat  = xy[, "Y"],
    Lon  = xy[, "X"]
  )
  
  setorder(WA_TOWNS, Town)
  SWAR_POLY <- swar_poly
  
} else {
  stop("Required shapefiles not found.")
}

# =========================================================================
# 1. WEATHER UTILITIES
# =========================================================================

normalize_weather_cols <- function(dt) {
  setDT(dt)
  dt[, date := as.Date(date)]
  
  # ---- STANDARDISE RAINFALL COLUMN NAME ----
  if ("rainfall" %in% names(dt)) {
    setnames(dt, "rainfall", "rain")
  }
  
  dt[, `:=`(
    year  = year(date),
    month = month(date),
    week  = isoweek(date)
  )]
  
  setorder(dt, date)
  dt[]
}


DBMWeatherFetcher <- R6::R6Class(
  "DBMWeatherFetcher",
  public = list(
    lat = NULL,
    lon = NULL,
    
    initialize = function(lat, lon) {
      self$lat <- lat
      self$lon <- lon
    },
    
    fetch = function(start_date = "20120101", end_date = Sys.Date()) {
      dt <- as.data.table(
        weatherOz::get_data_drill(
          latitude   = self$lat,
          longitude  = self$lon,
          start_date = start_date,
          end_date   = end_date,
          values     = c("rain", "max_temp", "min_temp")
        )
      )
      normalize_weather_cols(dt)
    }
  )
)

# =========================================================================
# 2. THERMAL PRESSURE UTILITIES
# =========================================================================

# ---- Canonical pressure definition ----
calc_pressure <- function(tmean, Tbase = 4, Tupper = 38) {
  pmax(0, pmin(tmean, Tupper) - Tbase) / (Tupper - Tbase)
}

# ---- Observed pressure (past N days) ----
# observed_pressure:
# pressure calculated from historical daily temperatures,
# used only for climatological reference (not actual observations)

compute_recent_14day_observed_pressure <- function(w, n_days = 14) {
  
  last_day <- max(w$date)
  
  dt <- w[date > (last_day - n_days) & date <= last_day]
  dt[, mean_temp := (air_tmax + air_tmin) / 2]
  dt[, pressure  := calc_pressure(mean_temp)]
  
  dt[, .(date, pressure)]
}

# ---- Historical climatology envelope for arbitrary dates ----
generate_climatology_envelope_for_dates <- function(
    w,
    target_dates,
    n_years = 10
) {
  
  year_cov <- w[, .(n_doy = uniqueN(yday(date))), by = year][n_doy >= 365]
  yrs <- tail(sort(year_cov$year), n_years)
  
  traj <- rbindlist(lapply(yrs, function(yr) {
    
    ref <- w[year == yr]
    ref[, doy := yday(date)]
    
    clim <- ref[
      ,
      .(
        air_tmax = mean(air_tmax, na.rm = TRUE),
        air_tmin = mean(air_tmin, na.rm = TRUE)
      ),
      by = doy
    ]
    
    dt <- data.table(
      date = target_dates,
      doy  = yday(target_dates)
    )
    
    dt <- merge(dt, clim, by = "doy", all.x = TRUE)
    dt[, mean_temp := (air_tmax + air_tmin) / 2]
    dt[, observed_pressure := calc_pressure(mean_temp)]
    dt[, .(date, observed_pressure)]
    
  }))
  
  traj[
    ,
    .(
      pressure_p20 = quantile(observed_pressure, 0.2, na.rm = TRUE),
      pressure_p50 = quantile(observed_pressure, 0.5, na.rm = TRUE),
      pressure_p80 = quantile(observed_pressure, 0.8, na.rm = TRUE)
    ),
    by = date
  ]
  
}


# ------------------------------------------------------------
# Climatological mean 
# ------------------------------------------------------------
compute_climatology_mean_for_dates <- function(
    w,
    target_dates,
    n_years = 10
) {
  
  year_cov <- w[, .(n_doy = uniqueN(yday(date))), by = year][n_doy >= 365]
  yrs <- tail(sort(year_cov$year), n_years)
  
  traj <- rbindlist(lapply(yrs, function(yr) {
    
    ref <- w[year == yr]
    ref[, doy := yday(date)]
    
    clim <- ref[
      ,
      .(
        air_tmax = mean(air_tmax, na.rm = TRUE),
        air_tmin = mean(air_tmin, na.rm = TRUE)
      ),
      by = doy
    ]
    
    dt <- data.table(
      date = target_dates,
      doy  = yday(target_dates)
    )
    
    dt <- merge(dt, clim, by = "doy", all.x = TRUE)
    
    dt[, mean_temp := (air_tmax + air_tmin) / 2]
    dt[, observed_pressure := calc_pressure(mean_temp)]
    
    dt[, .(date, observed_pressure)]
  }))
  
  traj[
    ,
    .(
      pressure_p50_median = quantile(observed_pressure, 0.5, na.rm = TRUE)
    ),
    by = date
  ]
}


# ---- Final combined object ----
compute_14day_pressure_result <- function(w, n_days = 14, n_years = 10) {
  
  last_day <- max(w$date)
  
  obs <- compute_recent_14day_observed_pressure(w, n_days)
  
  past_dates   <- obs$date
  future_dates <- seq(last_day + 1, last_day + n_days, by = "day")
  
  env_past <- generate_climatology_envelope_for_dates(
    w, past_dates, n_years
  )
  
  clim_mean_past <- compute_climatology_mean_for_dates(
    w, past_dates, n_years
  )
  
  env_future <- generate_climatology_envelope_for_dates(
    w, future_dates, n_years
  )
  
  
  
  today_row <- w[date == last_day]
  today_row[, mean_temp := (air_tmax + air_tmin) / 2]
  today_p <- calc_pressure(today_row$mean_temp)
  
  today_env <- data.table(
    date = last_day,
    pressure_p20 = today_p,
    pressure_p50 = today_p,
    pressure_p80 = today_p
  )
  
  env_future <- rbind(today_env, env_future)
  
  list(
    observed   = obs,
    env_past   = env_past,
    env_future = env_future,
    clim_mean_past  = clim_mean_past, 
    today      = last_day,
    avg14      = mean(obs$pressure, na.rm = TRUE)
  )
}



# daily -------------------------------------------------------------------
# ------------------------------------------------------------
# Daily thermal pressure 
# ------------------------------------------------------------
compute_daily_pressure <- function(w) {
  
  dt <- copy(w)
  dt[, mean_temp := (air_tmax + air_tmin) / 2]
  dt[, pressure  := calc_pressure(mean_temp)]
  
  dt[, .(
    date,
    year       = year(date),
    month      = month(date),
    cal_month  = month(date),
    pressure
  )]
}

# ------------------------------------------------------------
# Monthly pressure summary 
# ------------------------------------------------------------
compute_monthly_pressure_summary <- function(daily_p) {
  
  today <- max(daily_p$date)
  
  start_month <- floor_date(today, "month") %m-% months(11)
  end_month   <- ceiling_date(today, "month") - days(1)
  
  # ---- Monthly mean (recent 12 months)
  monthly_recent <- daily_p[
    date >= start_month & date <= end_month,
    .(
      mean_pressure = mean(pressure, na.rm = TRUE)
    ),
    by = .(
      year_month = floor_date(date, "month")
    )
  ]
  
  monthly_recent[, `:=`(
    cal_month = month(year_month),
    month_lab = format(year_month, "%b %Y")
  )]
  
  monthly_recent[, month_lab := factor(
    month_lab,
    levels = format(sort(unique(year_month)), "%b %Y")
  )]
  
  
  # ---- Historical variability by calendar month
  monthly_hist <- daily_p[
    ,
    .(
      p20 = quantile(pressure, 0.20, na.rm = TRUE),
      p50 = quantile(pressure, 0.50, na.rm = TRUE),
      p80 = quantile(pressure, 0.80, na.rm = TRUE)
    ),
    by = cal_month
  ]
  
  # ---- Merge
  out <- merge(
    monthly_recent,
    monthly_hist,
    by = "cal_month",
    all.x = TRUE
  )
  
  setorder(out, year_month)
  out
}

# =========================================================================
# 3. RAINFALL UTILITIES (CUMULATIVE + CLIMATOLOGY)
# =========================================================================
compute_recent_cumulative_rainfall <- function(w, n_days = 14) {
  
  dt <- copy(w)
  setorder(dt, date)
  
  dt[, rain_14 := zoo::rollapply(
    rain,
    width = n_days,
    FUN = sum,
    align = "right",
    fill = NA,
    na.rm = TRUE
  )]
  
  dt[!is.na(rain_14), .(date, rain_14)]
}


# -------------------------------------------------------------------------

generate_rainfall_climatology_for_dates <- function(
    w,
    target_dates,
    n_days = 14,
    n_years = 10
) {
  
  # ---- ensure full years
  year_cov <- w[, .(n_doy = uniqueN(yday(date))), by = year][n_doy >= 365]
  yrs <- tail(sort(year_cov$year), n_years)
  
  traj <- rbindlist(lapply(yrs, function(yr) {
    
    ref <- w[year %in% c(yr - 1, yr)]
    
    setorder(ref, date)
    
    ref[, rain_14 := zoo::rollapply(
      rain,
      width = n_days,
      FUN = sum,
      align = "right",
      fill = NA,
      na.rm = TRUE
    )]
    
    ref[, doy := yday(date)]
    
    dt <- data.table(
      date = target_dates,
      doy  = yday(target_dates)
    )
    
    merge(
      dt,
      ref[, .(doy, rain_14)],
      by = "doy",
      all.x = TRUE
    )[, .(date, rain_14)]
    
  }))
  
  traj[
    ,
    .(
      rain_p20 = quantile(rain_14, 0.20, na.rm = TRUE),
      rain_p50 = quantile(rain_14, 0.50, na.rm = TRUE),
      rain_p80 = quantile(rain_14, 0.80, na.rm = TRUE)
    ),
    by = date
  ]
}


# -------------------------------------------------------------------------

compute_14day_rainfall_result <- function(w, n_days = 14, n_years = 10, plot_window = 60) {
  
  obs_all <- compute_recent_cumulative_rainfall(w, n_days)
  
  last_day <- max(obs_all$date)
  obs <- obs_all[date >= last_day - plot_window]
  
  
  future_dates <- seq(
    last_day + 1,
    last_day + n_days,
    by = "day"
  )
  
  target_dates <- c(obs$date, future_dates)
  
  
  clim <- generate_rainfall_climatology_for_dates(
    w,
    target_dates,
    n_days  = n_days,
    n_years = n_years
  )
  
  list(
    observed = obs,
    climatology = clim,
    today = last_day
  )
}




# =========================================================================
library(shiny)

# -------------------------------
# STATIC FILE PATH (GLOBAL)
# -------------------------------
shiny::addResourcePath(
  "static",
  normalizePath("www")
)


# 3. UI
# =========================================================================
# 

ui <- fluidPage(
  
  titlePanel("DBM – Near-term Pest Pressure (14 Days)"),
  
  # ============================================================
  # INPUT SECTION
  # ============================================================
  fluidRow(
    column(
      width = 12,
      
      selectizeInput(
        "town_search",
        "Search SW Ag Region Town:",
        choices = c("", WA_TOWNS$Town)
      ),
      
      leafletOutput("map", height = 400),
      
      fluidRow(
        column(6, numericInput("lat", "Latitude", -32.2, step = 0.01)),
        column(6, numericInput("lon", "Longitude", 116.1, step = 0.01))
      ),
      
      actionButton("run", "RUN", class = "btn-success")
    )
  ),
  
  hr(),
  
  # ============================================================
  # OUTPUTS
  # ============================================================
  # SDM
  fluidRow(
    column(
      width = 4,
      
      tags$details(
        `data-open` = "false",
        
        tags$summary(
          tags$strong("Spatial context: habitat suitability (SDM)"),
          style = "font-size: 16px; cursor: pointer;"
        ),
        
        div(
          style = "margin-top: 10px;",
          
          tags$img(
            src   = "static/sdm_context.png",
            width = "100%",
            style = "max-width: 100%; height: auto;"
          ),
          
          tags$p(
            "Static SDM shown for spatial context only. Colours indicate relative habitat suitability.",
            style = "color: #555; margin-top: 8px;"
          )
        )
      ),
      
      style = "margin-bottom: 10px;"
    )
  ),
  
  
  # ============================================================
  # PLOTS
  fluidRow(
    column(
      width = 12,
      
      tags$details(
        `data-open` = "false",
        
      div(
        plotlyOutput("pressure_plot", height = 500),
        tags$div(
          "Temperature-driven pressure only.",
          style = "font-size: 18px; color: #555; margin-top: 0.5px;"
        ),
        style = "margin-bottom: 100px;"
      ),
      
      div(
        plotlyOutput("monthly_pressure_plot", height = 500),
        tags$div(
          "Growing degree-day-based pressure aggregated by calendar month",
          style = "font-size: 18px; color: #555; margin-top: 0.5px;"
        ),
        style = "margin-bottom: 100px;"
      ),
      
      div(
        plotlyOutput("rainfall_plot", height = 500),
        tags$div(
          "14-day cumulative rainfall relative to historical levels (10-year climatology).",
          style = "font-size: 18px; color: #555; margin-top: 0.5px;"
        ),
        style = "margin-bottom: 100px;"
      )
    ) )
  )
)



# =========================================================================
# 4. SERVER
# =========================================================================


server <- function(input, output, session) {
  
  output$map <- renderLeaflet({
    leaflet() |>
      addTiles() |>
      
      # --- LGA boundaries ---
      addPolygons(
        data        = lga_poly,
        fill        = FALSE,
        weight      = 1,
        color       = "grey60",
        opacity     = 0.6,
        group       = "LGA boundaries"
      ) |>
      
      # --- SWAR boundary (optional but recommended) ---
      addPolygons(
        data        = SWAR_POLY,
        fill        = FALSE,
        weight      = 2,
        color       = "blue",
        opacity     = 0.8,
        group       = "SWAR"
      ) |>
      
      setView(lng = input$lon, lat = input$lat, zoom = 6) |>
      addMarkers(lng = input$lon, lat = input$lat)
  })
  
  observeEvent(input$town_search, {
    coords <- WA_TOWNS[Town == input$town_search]
    if (nrow(coords)) {
      updateNumericInput(session, "lat", value = coords$Lat)
      updateNumericInput(session, "lon", value = coords$Lon)
      leafletProxy("map") |>
        setView(coords$Lon, coords$Lat, zoom = 9) |>
        clearMarkers() |>
        addMarkers(coords$Lon, coords$Lat)
    }
  })
  
  # MAP CLICK (unchanged, works)
  observeEvent(input$map_click, {
    click <- input$map_click
    if (!is.null(click$lat)) {
      
      updateNumericInput(session, "lat", value = round(click$lat, 4))
      updateNumericInput(session, "lon", value = round(click$lng, 4))
      
      leafletProxy("map") |>
        clearMarkers() |>
        addMarkers(lng = click$lng, lat = click$lat)
    }
  })
  
  
  # ---------------- WEATHER ----------------
  weather_data <- eventReactive(input$run, {
    DBMWeatherFetcher$new(input$lat, input$lon)$fetch()
  })
  
  pressure_data <- reactive({
    compute_14day_pressure_result(req(weather_data()))
  })
  
  monthly_pressure_data <- reactive({
    w <- req(weather_data())
    daily_p <- compute_daily_pressure(w)
    compute_monthly_pressure_summary(daily_p)
  })
  
  rainfall_data <- reactive({
    compute_14day_rainfall_result(req(weather_data()))
  })
  
  
  output$pressure_plot <- renderPlotly({
    
    d <- pressure_data()
    
    # ---- 14-day average line ----
    avg_line_df <- data.table(
      date = range(c(d$observed$date, d$env_future$date)),
      previous_14days_avg = d$avg14
    )
    
    # ---- hover helpers ----
    env_past_hover <- merge(
      d$env_past,
      d$clim_mean_past,
      by = "date",
      all.x = TRUE
    )
    
    g <- ggplot() +
      
      # ============================================================
    # 1. HISTORICAL RIBBON (P20–P80) 
    # ============================================================
    geom_ribbon(
      data = d$env_past,
      aes(date, ymin = pressure_p20, ymax = pressure_p80),
      fill = "grey70",
      alpha = 0.35,
      show.legend = FALSE
    ) +
      
      # ============================================================
    # 2. FUTURE RIBBON (P20–P80)
    # ============================================================
    geom_ribbon(
      data = d$env_future,
      aes(date, ymin = pressure_p20, ymax = pressure_p80),
      fill = "lightblue",
      alpha = 0.35,
      show.legend = FALSE
    ) +
      
      # ============================================================
    # 3. HISTORICAL P50
    # ============================================================
    
    # --- Historical P50: real visible line ---
    geom_line(
      data = d$clim_mean_past,
      aes(
        x = date,
        y = pressure_p50_median,
        group = 1,
        text = paste0(
          format(date, "%Y-%m-%d"),
          "<br>Historical median (P50): ",
          round(pressure_p50_median, 2)
        )
      ),
      inherit.aes = FALSE,
      linetype = "dotted",
      linewidth = 1.2,
      colour = "grey60",
      show.legend = FALSE
    ) +
      
      # --- Historical P50: legend proxy ---
      geom_line(
        data = d$clim_mean_past,
        aes(
          date,
          pressure_p50_median,
          colour = "Historical P50"
        ),
        linewidth = 1.2
      ) +
      
      
      # ============================================================
    # 4. FUTURE P50 LINE
    # ============================================================
    geom_line(
      data = d$env_future,
      aes(
        date,
        pressure_p50,
        colour = "Forecast median (P50)",
        text = paste0(
          format(date, "%Y-%m-%d"),
          "<br>Forecast median (P50): ",
          round(pressure_p50, 2)
        )
      ),
      linetype = "dashed",
      linewidth = 1.4
    ) +
      
      geom_line(
        data = d$env_future,
        aes(
          x = date,
          y = pressure_p50,
          group = 1,
          text = paste0(
            format(date, "%Y-%m-%d"),
            "<br>Forecast median pressure (P50): ",
            round(pressure_p50, 2)
          )
        ),
        inherit.aes = FALSE,
        linetype = "dashed",
        linewidth = 1.6,
        colour = "#1565C0"
      ) +
      # ============================================================
    # 5. OBSERVED LINE
    # ============================================================
    # --- Observed: real visible line ---
    geom_line(
      data = d$observed,
      aes(
        x = date,
        y = pressure,
        group = 1,
        text = paste0(
          format(date, "%Y-%m-%d"),
          "<br>Observed pressure: ",
          round(pressure, 2)
        )
      ),
      inherit.aes = FALSE,
      linewidth = 1.2,
      colour = "black",
      show.legend = FALSE
    ) +
      
      # --- Observed: legend 
      geom_line(
        data = d$observed,
        aes(
          date,
          pressure,
          colour = "Observed"
        ),
        linewidth = 1.2
      ) +
      
      
      # ============================================================
    # 6. 14-DAY AVERAGE LINE
    # ============================================================
    geom_line(
      data = avg_line_df,
      aes(
        date,
        previous_14days_avg,
        colour = "Observed 14-day average",
        text = paste0(
          "Observed 14-day average<br>",
          "Pressure: ", round(previous_14days_avg, 2)
        )
      ),
      linewidth = 0.6,
      alpha = 0.6
    ) +
      
      
      # ============================================================
    # 7. POINTS (observed)
    # ============================================================
    geom_point(
      data = d$observed,
      aes(
        date,
        pressure,
        colour = "Observed",
        text = paste0(
          format(date, "%Y-%m-%d"),
          "<br>Observed pressure: ",
          round(pressure, 2)
        )
      ),
      size = 2
    ) +
      
      # ============================================================
    # 8. RIBBON HOVER PROXIES 
    # ============================================================
    geom_point(
      data = d$env_past,
      aes(
        date,
        y = (pressure_p20 + pressure_p80) / 2,
        text = paste0(
          "<b>Historical range</b>",
          "<br>P20: ", round(pressure_p20, 2),
          "<br>P50: ", round(pressure_p50, 2),
          "<br>P80: ", round(pressure_p80, 2)
        )
      ),
      alpha = 0,
      size = 6,
      show.legend = FALSE
    ) +
      
      geom_point(
        data = d$env_future,
        aes(
          date,
          y = (pressure_p20 + pressure_p80) / 2,
          text = paste0(
            "<b>Forecast range</b>",
            "<br>P20: ", round(pressure_p20, 2),
            "<br>P50: ", round(pressure_p50, 2),
            "<br>P80: ", round(pressure_p80, 2)
          )
        ),
        alpha = 0,
        size = 6,
        show.legend = FALSE
      ) +
      
      # ============================================================
    # 9. CONTEXT 
    # ============================================================
    geom_vline(
      xintercept = as.numeric(d$today),
      linetype = "dashed",
      colour = "black"
    ) +
      
      # ============================================================
    
    # ---- Context labels
    annotate(
      "text",
      x = min(d$observed$date),
      y = 0.98,
      label = "Last 14 days",
      hjust = 0,
      vjust = 1,
      size = 3.5,
      colour = "grey30"
    ) +
      annotate(
        "text",
        x = max(d$env_future$date),
        y = 0.98,
        label = "Next 14 days",
        hjust = 1,
        vjust = 1,
        size = 3.5,
        colour = "grey30"
      ) +
      annotate(
        "text",
        x = d$today-0.6,
        y = 0.98,
        label = "Now",
        hjust = 1,
        vjust = 1,
        size = 3.5,
        colour = "grey30"
      ) +
      
      # 10. SCALES & THEME
      # ============================================================
    scale_colour_manual(
      name = NULL,
      values = c(
        "Observed"        = "black",
        "Historical P50"  = "grey60",
        "Forecast median (P50)"      = "#1565C0",
        "Observed 14-day average"  = "black"
      )
    ) +
      
      coord_cartesian(ylim = c(0, 1), clip = "off") +
      
      
      # labels ------------------------------------------------------------------
    labs(
      x = "Date",
      y = "Thermal pest pressure (0–1)",
      title = "DBM pest pressure: recent (past 14-day) conditions and 14-day outlook"
    ) +
      
      theme_minimal(base_size = 12) +
      theme(
        plot.margin = margin(10, 50, 10, 10),
        legend.position = "right"
      ) 
    
    ggplotly(g, tooltip = "text")
  })
  
  
  # monthly -----------------------------------------------------------------
  
  # monthly -----------------------------------------------------------------
  
  output$monthly_pressure_plot <- renderPlotly({
    
    d <- req(monthly_pressure_data())
    
    g <- ggplot(d, aes(x = month_lab)) +
      
      # -----------------------------------
    # Historical variability band (P20–P80)
    # -----------------------------------
    geom_rect(
      aes(
        xmin = as.numeric(month_lab) - 0.45,
        xmax = as.numeric(month_lab) + 0.45,
        ymin = p20,
        ymax = p80,
        fill = "Historical P20–P80",
        text = paste0(
          as.character(month_lab),
          "<br><b>Historical range (P20–P80)</b>",
          "<br>P20: ", round(p20, 2),
          "<br>P80: ", round(p80, 2)
        )
      ),
      alpha = 0.45,
      show.legend = TRUE
    ) +
      
      # -----------------------------------
    # Monthly mean pressure (bars)
    # -----------------------------------
    geom_col(
      aes(
        y = mean_pressure,
        fill = "Observed Monthly mean pressure",
        text = paste0(
          as.character(month_lab),
          "<br><b>Observed Monthly mean pressure</b>: ",
          round(mean_pressure, 2)
        )
      ),
      width = 0.6,
      show.legend = TRUE
    ) +   
      
      # -----------------------------------
    # Historical median (P50) – real line
    # -----------------------------------
    geom_line(
      aes(
        y = p50,
        group = 1,
        text = paste0(
          as.character(month_lab),
          "<br><b>Historical median (P50)</b>: ",
          round(p50, 2)
        )
      ),
      colour = "grey60",
      linetype = "dotted",
      linewidth = 1,
      show.legend = FALSE
    ) +
      
      geom_point(
        aes(y = p50),
        colour = "grey60",
        size = 2,
        show.legend = FALSE
      ) +
      
      # -----------------------------------
    # Historical median (P50) – legend proxy
    # -----------------------------------
    geom_line(
      aes(
        y = p50,
        group = 1,
        colour = "Historical P50"
      ),
      linewidth = 1,
      linetype = "dotted"
    ) +
      
      # -----------------------------------
    # Scales
    # -----------------------------------
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2)
    ) +
      
      scale_fill_manual(
        name = NULL,
        values = c(
          "Observed Monthly mean pressure"  = "darkblue",
          "Historical P20–P80"     = "grey75"
        )
      ) +
      
      scale_colour_manual(
        name = NULL,
        values = c(
          "Historical P50" = "grey60"
        )
      ) +
      
      # -----------------------------------
    # Labels & theme
    # -----------------------------------
    labs(
      x = NULL,
      y = "Observed Monthly mean pressure (0–1)",
      title = "Monthly DBM growing degree-day pressure (last 12 months)",
      subtitle = "Bars = recent observed monthly mean pressure; shaded bands = historical P20–P80"
    ) +
      
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text.x        = element_text(angle = 45, hjust = 1),
        plot.title         = element_text(face = "bold"),
        plot.margin        = margin(10, 120, 10, 10),
        legend.position    = "right"
      )
    
    ggplotly(g, tooltip = "text")
  })
  
  
  
  # rainfall ----------------------------------------------------------------
  # -------------------------------------------------------------------------
  
  output$rainfall_plot <- renderPlotly({
    
    d <- rainfall_data()
    
    today <- d$today
    
    clim_past   <- d$climatology[date <= today]
    clim_future <- d$climatology[date >  today]
    obs <- d$observed
    
    y_max <- max(
      c(d$climatology$rain_p80, obs$rain_14),
      na.rm = TRUE
    )
    
    g <- ggplot() +
      
      # ============================================================
    # 1. HISTORICAL CLIMATOLOGY RIBBON 
    # ============================================================
    geom_ribbon(
      data = clim_past,
      aes(date, ymin = rain_p20, ymax = rain_p80),
      fill  = "grey70",
      alpha = 0.35,
      show.legend = FALSE
    ) +
      
      geom_ribbon(
        data = clim_future,
        aes(date, ymin = rain_p20, ymax = rain_p80),
        fill  = "grey70",
        alpha = 0.18,
        show.legend = FALSE
      ) +
      
      # ============================================================
    # 2. HISTORICAL P50 
    # ============================================================
    geom_line(
      data = d$climatology,
      aes(
        date,
        rain_p50,
        group = 1,
        text = paste0(
          format(date, "%Y-%m-%d"),
          "<br><b>Historical median (P50)</b>: ",
          round(rain_p50, 1), " mm"
        )
      ),
      linetype = "dotted",
      linewidth = 1.2,
      colour = "grey40",
      show.legend = FALSE
    ) +
      
      # ============================================================
    # 3. OBSERVED 14-DAY RAINFALL 
    # ============================================================
    geom_col(
      data = obs,
      aes(
        date,
        rain_14,
        text = paste0(
          format(date, "%Y-%m-%d"),
          "<br><b>Observed 14-day rainfall</b>: ",
          round(rain_14, 1), " mm"
        )
      ),
      width = 0.8,
      fill  = "#90CAF9",
      colour = "#90CAF9",
      show.legend = FALSE
    ) +
      
      # ============================================================
    # 4. LEGEND PROXIES 
    # ============================================================
    
    # Ribbon legend
    geom_ribbon(
      data = clim_past[1],
      aes(
        date,
        ymin = rain_p20,
        ymax = rain_p80,
        fill = "Historical P20–P80"
      ),
      alpha = 0.35
    ) +
      
      # P50 legend
      geom_line(
        data = d$climatology[1],
        aes(
          date,
          rain_p50,
          colour = "Historical median (P50)"
        ),
        linetype = "dotted",
        linewidth = 1.2
      ) +
      
      # Observed rainfall legend
      geom_col(
        data = obs[1],
        aes(
          date,
          rain_14,
          fill = "Observed 14-day rainfall"
        ),
        width = 0.8
      ) +
      
      # ============================================================
    # 5. TODAY MARKER
    # ============================================================
    geom_vline(
      xintercept = as.numeric(today),
      linetype   = "dashed",
      colour     = "black"
    ) +
      
      # ============================================================
    # 6. SCALES
    # ============================================================
    scale_fill_manual(
      name = NULL,
      values = c(
        "Observed 14-day rainfall" = "#90CAF9",
        "Historical P20–P80"       = "grey70"
      )
    ) +
      
      scale_colour_manual(
        name = NULL,
        values = c(
          "Historical median (P50)" = "grey40"
        )
      ) +
      
      coord_cartesian(
        ylim = c(0, y_max * 1.15),
        clip = "off"
      ) +
      
      # ============================================================
    # 7. LABELS & THEME
    # ============================================================
    labs(
      x = NULL,
      y = "14-day cumulative rainfall (mm)",
      title = "Recent rainfall conditions",
      subtitle = "Observed rainfall relative to historical distribution (10-year climatology)"
    ) +
      
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "right",
        plot.margin = margin(10, 50, 10, 10)
      )
    
    ggplotly(g, tooltip = "text")
  })
  
} 



shinyApp(ui, server)



