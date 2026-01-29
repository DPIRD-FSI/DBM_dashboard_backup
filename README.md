DBM Weekly Risk Dashboard

##Forecasting near‑term Diamondback Moth (DBM) pest pressure in WA broadacre cropping using temperature‑based GDD models and short‑term rainfall signals.

This Shiny dashboard provides interactive visual diagnostics and weekly risk scoring based on:
- Thermal pressure model (GDD‑based)
- 14‑day cumulative rainfall climatology
- Static habitat suitability context (SDM)
- Interactive map‑based location selection
- Climatology envelopes (P20–P80) for interpretation

The app fetches historical weather data from SILO and computes thermal and rainfall‑based pest pressure using modular utilities.

📦 Installation (local development version)

# Clone the repository:
git clone https://github.com/DPIRD-FSI/DBM_NewDashboard.git
cd DBM_dashboard_backup

# required packages
install.packages(c(
  "shiny", "data.table", "ggplot2", "plotly",
  "leaflet", "terra", "R6", "sf"
))

Project Structure:
DBM_NewDashboard/
 R/ 
- utils.R          # Weather fetcher, thermal/rainfall models, helpers
- app_ui.R         # Shiny UI layout
- app_server.R     # Server logic with reactive risk computation
- run_app.R        # Entry point to start dashboard

 www/
 
 sdm_context.png  # Static SDM background map (2023)

-README.md


# Features
- 14‑Day Thermal Pressure
Uses GDD‑based canonical pressure definition & 10‑year climatology envelopes.

- 14‑Day Rainfall Pressure
Rolling cumulative rainfall with climatological P20/P50/P80 ranges.

- Monthly Thermal Aggregation
Calendar‑month mean pest pressure with historical variability band.

- Interactive Map Selection
Search towns, click on map, or manually enter lat/lon.

- SDM Habitat Suitability Context
Static raster preview (2023) for spatial interpretation.


If you have any questions or suggestions, please reach out to us:
Mia Li (mia.li@dpird.wa.gob.au). Rodrigo Pires (rodrigo.pires@dpird.wa.gov.au). Department of Primary Industries and Regional Development (WA)

