library(shiny)
library(leaflet)
library(sf)
library(ggplot2)
library(rphylopic)
library(smoothr)
library(base64enc)  # ADD FOR GIF DISPLAY


# STEP 2: mount www/ as “/images”
addResourcePath("logos", "www")

nested_list <- readRDS("nested_list_polygons.rds")
comprehensive_metrics <- read.csv("comprehensive_metrics.csv")
silhouette_data <- readRDS("species_silhouettes_data.rds")  
animation_base_path <- "Animations"



# Smoothing function for polygons
smooth_polygon <- function(polygon_obj, method = "chaikin", refinements = 5) {
  if(is.null(polygon_obj)) return(NULL)
  
  tryCatch({
    # Use smoothr package for advanced smoothing
    smoothed <- smoothr::smooth(polygon_obj, method = method)
    return(smoothed)
  }, error = function(e) {
    return(polygon_obj)
  })
}




# Function to get bounds from a polygon
get_polygon_bounds <- function(polygon_obj) {
  if(is.null(polygon_obj)) return(NULL)
  
  # Get polygon bounding box
  bbox <- st_bbox(polygon_obj)
  
  return(list(
    xmin = bbox[1],
    xmax = bbox[3],
    ymin = bbox[2],
    ymax = bbox[4]
  ))
}

# Function to find maximum bounds for all ICs and time periods for a species
find_species_max_bounds <- function(species_data, keys) {
  cat("Calculating maximum bounds for species...\n")
  
  all_bounds <- list()
  
  # Get present bounds
  if(!is.null(species_data[[keys$present]])) {
    present_bounds <- get_polygon_bounds(species_data[[keys$present]])
    if(!is.null(present_bounds)) {
      all_bounds <- append(all_bounds, list(present_bounds))
    }
  }
  
  # Check both time periods
  for(period_key in c(keys$future_2030, keys$future_2070)) {
    if(!is.null(period_key) && !is.null(species_data[[period_key]])) {
      period_data <- species_data[[period_key]]
      
      if(is.list(period_data)) {
        # Check all initial conditions
        for(ic in 1:min(100, length(period_data))) {
          if(!is.null(period_data[[ic]])) {
            ic_bounds <- get_polygon_bounds(period_data[[ic]])
            if(!is.null(ic_bounds)) {
              all_bounds <- append(all_bounds, list(ic_bounds))
            }
          }
        }
      } else {
        # Single polygon case
        period_bounds <- get_polygon_bounds(period_data)
        if(!is.null(period_bounds)) {
          all_bounds <- append(all_bounds, list(period_bounds))
        }
      }
    }
  }
  
  if(length(all_bounds) == 0) {
    cat("No valid bounds found for species\n")
    return(NULL)
  }
  
  # Find maximum bounds across all data
  max_bounds <- list(
    xmin = min(sapply(all_bounds, function(b) b$xmin), na.rm = TRUE),
    xmax = max(sapply(all_bounds, function(b) b$xmax), na.rm = TRUE),
    ymin = min(sapply(all_bounds, function(b) b$ymin), na.rm = TRUE),
    ymax = max(sapply(all_bounds, function(b) b$ymax), na.rm = TRUE)
  )
  
  # Check for global/very wide distributions
  longitude_range <- max_bounds$xmax - max_bounds$xmin
  latitude_range <- max_bounds$ymax - max_bounds$ymin
  
  cat("Initial bounds: lon range =", longitude_range, ", lat range =", latitude_range, "\n")
  
  # If distribution is very wide (global-like), use world bounds instead
  if(longitude_range > 300 || max_bounds$xmin < -180 || max_bounds$xmax > 180) {
    cat("Detected global distribution - using world bounds\n")
    max_bounds <- list(
      xmin = -180,
      xmax = 180,
      ymin = max(max_bounds$ymin, -90),
      ymax = min(max_bounds$ymax, 90)
    )
  } else {
    # Clamp to world limits for safety
    max_bounds$xmin <- max(max_bounds$xmin, -180)
    max_bounds$xmax <- min(max_bounds$xmax, 180)
    max_bounds$ymin <- max(max_bounds$ymin, -90)
    max_bounds$ymax <- min(max_bounds$ymax, 90)
  }
  
  cat("Final bounds: xmin =", max_bounds$xmin, ", xmax =", max_bounds$xmax, 
    ", ymin =", max_bounds$ymin, ", ymax =", max_bounds$ymax, "\n")
  
  return(max_bounds)
}

# SPECIES CLASSIFICATION MATRIX
species_classification <- matrix(c(
  "Aedes_aegypti", "Terrestrial",
  "Aedes_albopictus", "Terrestrial",
  "Alces_alces", "Terrestrial",
  "Anopheles_gambiae", "Terrestrial",
  "Balaenoptera_musculus", "Marine",
  "Chelonia_mydas", "Marine",
  "Culex_pipiens", "Terrestrial",
  "Culex_quinquefasciatus", "Terrestrial",
  "Dioscorea_villosa", "Terrestrial",
  "Dosidicus_gigas", "Marine",
  "Gadus_morhua", "Marine",
  "Homarus_americanus", "Marine",
  "Iguana_iguana", "Terrestrial",
  "Ixodes_ricinus", "Terrestrial",
  "Ixodes_scapularis", "Terrestrial",
  "Loxodonta_africana", "Terrestrial",
  "Macrocystis_pyrifera", "Marine",
  "Odocoileus_verginianus", "Terrestrial",
  "Peromyscus_maniculatus", "Terrestrial",
  "Pinus_sylvestris", "Terrestrial",
  "Pocillopora_damicornis", "Marine",
  "Prionace_glauga", "Marine",
  "Psittacus_erithacus", "Terrestrial",
  "Pteropus_combine", "Terrestrial",
  "Quercus_robur", "Terrestrial",
  "Rana_temporaria", "Terrestrial",
  "Rattus_norvegicus", "Terrestrial",
  "Rattus_rattus", "Terrestrial",
  "Rhizophora_mangle", "Marine",
  "Sargassum_muticum", "Marine",
  "Thunnus_albacares", "Marine",
  "Triatoma_infestans", "Terrestrial",
  "Vulpes_vulves", "Terrestrial",
  "Zostera_marina", "Marine"
), ncol = 2, byrow = TRUE)

# Convert to data frame for easier handling
species_df <- data.frame(
  species = species_classification[,1],
  habitat = species_classification[,2],
  stringsAsFactors = FALSE
)

# Get all available species from nested_list
all_species <- names(nested_list)

# Filter species_df to only include species that exist in nested_list
species_df <- species_df[species_df$species %in% all_species, ]

# Separate species by habitat type
terrestrial_species <- species_df$species[species_df$habitat == "Terrestrial"]
marine_species <- species_df$species[species_df$habitat == "Marine"]

cat("Species classification:\n")
cat("Terrestrial species (", length(terrestrial_species), "):", paste(terrestrial_species, collapse = ", "), "\n")
cat("Marine species (", length(marine_species), "):", paste(marine_species, collapse = ", "), "\n")

# Function to detect keys for a given species
detect_keys <- function(species_data) {
  all_keys <- names(species_data)
  
  # Try to find present key
  present_patterns <- c("present", "current", "baseline", "present_binary", "current_binary")
  present_key <- NULL
  for(pattern in present_patterns) {
    matches <- all_keys[grepl(pattern, all_keys, ignore.case = TRUE)]
    if(length(matches) > 0) {
      present_key <- matches[1]
      break
    }
  }
  if(is.null(present_key)) present_key <- all_keys[1]
  
  # Try to find 2030-2069 key
  future_2030_patterns <- c("2030", "30_69", "2030_2069")
  future_2030_key <- NULL
  for(pattern in future_2030_patterns) {
    matches <- all_keys[grepl(pattern, all_keys, ignore.case = TRUE)]
    if(length(matches) > 0) {
      future_2030_key <- matches[1]
      break
    }
  }
  if(is.null(future_2030_key) && length(all_keys) > 1) future_2030_key <- all_keys[2]
  
  # Try to find 2070-2099 key
  future_2070_patterns <- c("2070", "70_99", "2070_2099")
  future_2070_key <- NULL
  for(pattern in future_2070_patterns) {
    matches <- all_keys[grepl(pattern, all_keys, ignore.case = TRUE)]
    if(length(matches) > 0) {
      future_2070_key <- matches[1]
      break
    }
  }
  if(is.null(future_2070_key)) future_2070_key <- all_keys[length(all_keys)]
  
  return(list(
    present = present_key,
    future_2030 = future_2030_key,
    future_2070 = future_2070_key,
    all_keys = all_keys
  ))
}

# Define UI
ui <- fluidPage(
  
  tags$head(
    tags$style(HTML("
      .modal-dialog { margin-top: 500px !important; margin-left: 700px !important; }
    "))
  ),
  
  titlePanel("Do Initial Conditions Matter for Species Distribution Models?"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      h4("Controls"),
      
      selectInput("habitat_type",
        "Habitat Type:",
        choices = list("Terrestrial" = "Terrestrial",
          "Marine" = "Marine"),
        selected = "Terrestrial"),
      
      selectInput("species",
        "Species:",
        choices = NULL),
      
      selectInput("time_period",
        "Time Period:",
        choices = list("2030-2069" = "2030_2069",
          "2070-2099" = "2070_2099"),
        selected = "2030_2069"),
      
      conditionalPanel(
        condition = "input.map_tabs == 'navigation'",
        numericInput("initial_condition",
          "Initial Condition:",
          value = 1,
          min = 1,
          max = 100,
          step = 1)
      ),
      
      conditionalPanel(
        condition = "input.map_tabs == 'animation'",
        h5("Animation Info:"),
        p("Showing all 100 initial conditions", style = "font-style: italic; color: gray;"),
        p("Each frame represents one initial condition", style = "font-style: italic; color: gray;")
      ),
      
      hr(),
      
      h4("How to Use", style = "color: #2c3e50;"),
      
      div(
        style = "background-color: #ADD8E6; padding: 15px; border-radius: 8px; border-left: 4px solid #007bff;",
        
        # ADD THIS TEXT HERE - RIGHT AFTER THE BLUE DIV OPENS
        p("Explore how initial condition uncertainty in the CESM2-LENS2 climate model translates into variability in Species Distribution Models projections. This app demonstrates that climate model uncertainty significantly affects predictions of species distribution under future climate change.", 
          style = "margin-bottom: 40px; font-size: 16px; color: #495057;"),
        
        
        h5("📍 Upper Section:", style = "text-decoration: underline; font-weight: bold; color: #495057; margin-bottom: 10px;"),
        p("1. Select your Habitat Type", style = "margin-bottom: 5px;"),
        p("2. Choose a Species",      style = "margin-bottom: 5px;"),
        p("3. Pick a Time Period",     style = "margin-bottom: 5px;"),
        p("4. Set Initial Condition (1-100) for Navigation Mode", style = "margin-bottom: 15px;"),
        
        h5("🗺️ Viewing Results:", style = "text-decoration: underline; font-weight: bold; color: #495057; margin-bottom: 10px;"),
        p("• Navigation Mode: Manually explore individual initial conditions", style = "margin-bottom: 5px;"),
        p("• Animation Mode: Watch all 100 initial conditions as a GIF", style = "margin-bottom: 15px;"),
        
        h5("📊 Lower Section:", style = "text-decoration: underline; font-weight: bold; color: #495057; margin-bottom: 10px;"),
        p("• View summary plot for your selected species", style = "margin-bottom: 5px;"),
        p("• Choose different metrics to analyze", style = "margin-bottom: 5px;"),
        p("• Compare your species with up to 5 other species", style = "margin-bottom: 5px;")
      ),
      
      
      # LOGOS GO HERE - AFTER BLUE BOX, BEFORE SIDEBARPANEL CLOSES
  hr(),
  div(
    style = "text-align: center; margin-top: 20px;",
    
    # NEU Logo with SDS Lab link
    div(
      style = "margin-bottom: 10px;",
      tags$a(
        href = "https://sdslab.io/",
        target = "_blank",
        tags$img(src = "logos/neu_logo.png", 
                 style = "max-width: 80%; height: auto; max-height: 50px; cursor: pointer;")
      )
    ),
    
    # GMRI Logo with gmri.org link
    div(
      style = "margin-bottom: 10px;",
      tags$a(
        href = "https://gmri.org",
        target = "_blank",
        tags$img(src = "logos/gmri_logo.png", 
                 style = "max-width: 80%; height: auto; max-height: 40px; cursor: pointer;")
      )
    )
  )  # <-- This closes the logo div (NO comma)
      
      
      
    ),
    
    mainPanel(
      width = 9,
      
      tabsetPanel(
        id = "map_tabs",
        type = "tabs",
        
        tabPanel("Navigation Mode",
          value = "navigation",
          br(),
          leafletOutput("map", height = "600px")
        ),
        
        tabPanel("Animation Mode",
          value = "animation",
          br(),
          htmlOutput("animation_display")
        )
      ),
      
      hr(),
      
      h4("Metrics Summary", align = "left", style = "font-size:24px;"),
      
      fluidRow(
        column(12, align = "center",
          column(6, offset = 3,
            div(
              style = "display: flex; align-items: center; justify-content: center; gap: 3px; margin-bottom: 20px;",
              
              # Metric selector
              div(
                style = "",
                selectInput("selected_metric", 
                  "Select Metric:",
                  choices = list(
                    "Area Change (%)" = "area_pct_diff",
                    "Gained Area (%)" = "percentage_expansion",
                    "Lost Area (%)" = "percentage_lost",
                    "Fragmentation Change (%)" = "fragmentation_pct_diff",
                    "Fractal Dimension Change (%)" = "fractal_pct_diff",
                    "Perimeter-Area Ratio Change (%)" = "pa_ratio_pct_diff",
                    "Suitability Ratio Change (%)" = "suitability_increase_ratio"
                  ),
                  selected = "area_pct_diff")
              ),
              
              # Question mark help button
              div(
                style = "margin-top: 10px;",
                actionButton("metric_help", "?", 
                  style = "background-color: #17a2b8; color: white; border: none; border-radius: 50%; width: 22px; height: 22px; font-weight: bold; font-size: 11px; cursor: pointer; display: flex; align-items: center; justify-content: center;",
                  title = "Click for metric explanation")
              )
            )
          )
        )
      ),
      
      fluidRow(
        column(2,
          h5("Your Current Species", align = "center"),
          div(textOutput("current_species_display"),
            style = "font-weight: bold; text-align: center; width: 100%; font-style: italic;")
        ),
        column(2,
          h5("Compare to", align = "center"),
          selectInput("comparison_species_1", "", choices = NULL)
        ),
        column(2,
          h5("Compare to", align = "center"),
          selectInput("comparison_species_2", "", choices = NULL)
        ),
        column(2,
          h5("Compare to", align = "center"),
          selectInput("comparison_species_3", "", choices = NULL)
        ),
        column(2,
          h5("Compare to", align = "center"),
          selectInput("comparison_species_4", "", choices = NULL)
        ),
        column(2,
          h5("Compare to", align = "center"),
          selectInput("comparison_species_5", "", choices = NULL)
        )
      ),
      
      plotOutput("unified_comparison_plot", height = "450px")
    )
  ),  # <-- COMMA here!
  
 
)  # end fluidPage


# Define Server
server <- function(input, output, session) {
  
  # Store species bounds
  species_bounds <- reactiveVal(NULL)
  
  # Define metric explanations
  metric_explanations <- list(
    "area_pct_diff" = list(
      title = "Area Change (%)",
      explanation = "Percentage change in total suitable area between present and future. Values range from -100 to +100: positive values indicate total area expansion, while negative values indicate total area contraction. Every single dot in the plot below indicates a projection from a different Initial Condition member from the CESM-2 LENS model."
    ),
    "percentage_expansion" = list(
      title = "Gained Area (%)", 
      explanation = "Percentage change of exclusively new gained areas in future projections relative to the current distribution. Values range from 0 to +100. Every dot in the plot below indicates a different Initial Condition projection from the CESM-2 LENS model."
    ),
    "percentage_lost" = list(
      title = "Lost Area (%)",
      explanation = "Percentage change of exclusively lost area in future projections relative to the current distribution. Values range from 0 to +100. Every dot in the plot below indicates a different Initial Condition projection from the CESM-2 LENS model."
    ),
    "fragmentation_pct_diff" = list(
      title = "Fragmentation Change (%)",
      explanation = "Percentage change in distribution fragmentation, measuring how broken up and scattered suitable areas become. Values range from -100 to +100: positive values indicate increasing fragmentation (suitable areas splitting into smaller, more isolated patches), while negative values indicate decreasing fragmentation (patches becoming larger and more connected). Every dot represents a different Initial Condition projection from the CESM-2 LENS model."
      
    ),
    "fractal_pct_diff" = list(
      title = "Fractal Dimension Change (%)",
      explanation = "Percentage change in distribution boundary complexity, measured using fractal dimension analysis. Values range from -100 to +100: positive values indicate distribution edges becoming more convoluted, twisted, and geometrically complex, while negative values indicate boundaries becoming smoother and simpler. Every dot represents a different Initial Condition projection from the CESM-2 LENS model."
      
    ),
    "pa_ratio_pct_diff" = list(
      title = "Perimeter-Area Ratio Change (%)",
      explanation = "Percentage change in distribution shape complexity, measured as the ratio of habitat boundary length to total habitat area. Values range from -100 to +100: positive values indicate habitats becoming more elongated, fragmented, or irregular (higher edge-to-area ratio), while negative values indicate habitats becoming more compact and circular (lower edge-to-area ratio). Every dot represents a different Initial Condition projection from the CESM-2 LENS model."
      
    ),
    "suitability_increase_ratio" = list(
      title = "Suitability Ratio Change (%)",
      explanation = "Net balance of habitat suitability changes across the landscape. This metric compares pixels where suitability significantly increased (>10%) versus decreased (>10%) between present and future projections. Values range from -100 to +100: positive values indicate more areas gained suitable habitat than lost it, negative values indicate more areas lost suitable habitat than gained it, and zero indicates equal gains and losses. Every dot represents a different Initial Condition projection from the CESM-2 LENS model."
      
    )
  )
  
  # Update species choices based on selected habitat type
  observeEvent(input$habitat_type, {
    cat("Habitat type changed to:", input$habitat_type, "\n")
    
    # Filter species by habitat type
    if(input$habitat_type == "Terrestrial") {
      filtered_species <- terrestrial_species
    } else if(input$habitat_type == "Marine") {
      filtered_species <- marine_species
    } else {
      filtered_species <- c()
    }
    
    cat("Filtered species:", paste(filtered_species, collapse = ", "), "\n")
    
    # Create choices with clean names (display names)
    species_choices <- setNames(filtered_species, gsub("_", " ", filtered_species))
    
    # Update the species selectInput
    updateSelectInput(session, "species",
      choices = species_choices,
      selected = if(length(filtered_species) > 0) filtered_species[1] else NULL)
  }, ignoreInit = FALSE)
  
  # Update comparison species choices
  observe({
    # Get unique species from comprehensive_metrics
    if("species" %in% colnames(comprehensive_metrics)) {
      available_species <- unique(comprehensive_metrics$species)
    } else {
      available_species <- all_species
    }
    
    species_choices <- c("None" = "", setNames(available_species, gsub("_", " ", available_species)))
    
    # Update all comparison species selectors
    updateSelectInput(session, "comparison_species_1", choices = species_choices)
    updateSelectInput(session, "comparison_species_2", choices = species_choices)
    updateSelectInput(session, "comparison_species_3", choices = species_choices)
    updateSelectInput(session, "comparison_species_4", choices = species_choices)
    updateSelectInput(session, "comparison_species_5", choices = species_choices)
  })
  
  # Display current species name
  output$current_species_display <- renderText({
    req(input$species)
    gsub("_", " ", input$species)
  })
  
  # Metric help modal
  observeEvent(input$metric_help, {
    req(input$selected_metric)
    
    metric_info <- metric_explanations[[input$selected_metric]]
    
    showModal(modalDialog(
      title = div(
        style = "color: #2c3e50; font-weight: bold; font-size: 18px; text-align: center;",
        icon("info-circle", style = "margin-right: 10px; color: #17a2b8;"),
        metric_info$title
      ),
      
      div(
        style = "font-size: 15px; line-height: 1.7; padding: 20px; background-color: #f8f9fa; border-radius: 8px; border-left: 4px solid #17a2b8; margin: 15px 0;",
        p(metric_info$explanation, style = "margin: 0; text-align: justify;")
      ),
      
      footer = tagList(
        tags$button("Got it!", 
          type = "button", 
          class = "btn", 
          `data-dismiss` = "modal",
          style = "background-color: #17a2b8; color: white; border: none; padding: 10px 25px; border-radius: 5px; font-weight: bold; min-width: 100px;")
      ),
      
      easyClose = TRUE,
      size = "m"
    ))
  })
  
  # ANIMATION DISPLAY
  output$animation_display <- renderUI({
    req(input$species, input$time_period)
    
    # Construct animation file path
    animation_filename <- paste0(input$species, "_", input$time_period, "_animation.gif")
    animation_path <- file.path(animation_base_path, input$time_period, animation_filename)
    
    # Check if animation file exists
    if(file.exists(animation_path)) {
      
      # Read GIF file and convert to base64
      gif_raw <- readBin(animation_path, "raw", file.info(animation_path)$size)
      gif_base64 <- base64encode(gif_raw)
      
      # Create HTML to display the GIF
      div(
        style = "text-align: center; padding: 20px;",
        
        # Display the GIF
        tags$img(
          src = paste0("data:image/gif;base64,", gif_base64),
          style = "max-width: 100%; height: auto; border: 2px solid #ddd; border-radius: 8px; box-shadow: 0 4px 8px rgba(0,0,0,0.1);",
          alt = paste("Animation for", gsub("_", " ", input$species))
        ),
        
        br(), br(),
        p("This animation shows the species distribution across all 100 initial conditions",
          style = "color: #666; font-style: italic;")
      )
      
    } else {
      # Animation not found
      div(
        style = "text-align: center; padding: 50px;",
        h3("Animation Not Available", style = "color: #999;"),
        br(),
        p(paste("Animation for", gsub("_", " ", input$species), "in period", gsub("_", "-", input$time_period), "not found."),
          style = "color: #666;"),
        br(),
        p("Expected file:", animation_filename, style = "font-family: monospace; color: #999; font-size: 12px;"),
        br(),
        icon("exclamation-triangle", style = "font-size: 48px; color: #ffc107;")
      )
    }
  })
  
  # Reactive function to get current species data and keys
  current_species_data <- reactive({
    req(input$species)
    
    species_data <- nested_list[[input$species]]
    keys <- detect_keys(species_data)
    
    return(list(
      data = species_data,
      keys = keys,
      present_polygon = species_data[[keys$present]],
      future_2030_data = species_data[[keys$future_2030]],
      future_2070_data = species_data[[keys$future_2070]]
    ))
  })
  
  # Calculate species bounds when species changes
  observeEvent(input$species, {
    req(input$species)
    
    cat("Species changed to:", input$species, ". Calculating optimal bounds...\n")
    
    species_info <- current_species_data()
    max_bounds <- find_species_max_bounds(species_info$data, species_info$keys)
    
    species_bounds(max_bounds)
  })
  
  # Function to get current polygons
  # Function to get current polygons
  get_current_polygons <- function() {
    species_info <- current_species_data()
    
    # Get present polygon
    present_poly <- species_info$present_polygon
    present_poly <- smooth_polygon(present_poly)  # Add this line
    
    # Get future data based on time period
    if(isolate(input$time_period) == "2030_2069") {
      future_data <- species_info$future_2030_data
      period_label <- "2030-2069"
    } else {
      future_data <- species_info$future_2070_data
      period_label <- "2070-2099"
    }
    
    # Get current IC polygon
    ic <- isolate(input$initial_condition)
    current_polygon_data <- NULL
    
    if(is.list(future_data) && length(future_data) >= ic) {
      current_polygon_data <- future_data[[ic]]
    } else if(is.list(future_data) && ic == 1) {
      current_polygon_data <- future_data[[1]]
    } else {
      current_polygon_data <- future_data
    }
    
    current_polygon_data <- smooth_polygon(current_polygon_data)  # Add this line
    
    return(list(
      present_poly = present_poly,
      future_poly = current_polygon_data,
      period_label = period_label,
      ic = ic
    ))
  }
  
  # NAVIGATION MODE MAP RENDERING
  output$map <- renderLeaflet({
    req(input$species)
    
    cat("Rendering map for species:", input$species, "\n")
    
    # Get current data (already polygons)
    current_data <- get_current_polygons()
    
    # Create leaflet map
    map <- leaflet(
      options = leafletOptions(
        minZoom = 2,
        maxZoom = 18,
        worldCopyJump = TRUE,
        maxBounds = list(
          list(-90, -180),
          list(90, 180)
        )
      )
    ) %>%
      addTiles()
    
    # Add present distribution (light blue with blue border)
    if(!is.null(current_data$present_poly)) {
      map <- map %>%
        addPolygons(
          data = current_data$present_poly,
          fillColor = "lightblue",
          fillOpacity = 0.3,
          color = "blue",
          weight = 2,
          layerId = "present_layer",
          smoothFactor = 0.5
        )
    }
    
    # Add future distribution (light coral with red border)
    if(!is.null(current_data$future_poly)) {
      map <- map %>%
        addPolygons(
          data = current_data$future_poly,
          fillColor = "lightcoral",
          fillOpacity = 0.3,
          color = "red",
          weight = 2,
          layerId = "future_layer",
          smoothFactor = 0.5
        )
    }
    
    # Apply bounds
    bounds <- species_bounds()
    if(!is.null(bounds)) {
      longitude_range <- bounds$xmax - bounds$xmin
      
      if(longitude_range >= 300 || (bounds$xmin <= -180 && bounds$xmax >= 180)) {
        map <- map %>% setView(lng = 0, lat = 0, zoom = 2)
      } else {
        x_padding <- (bounds$xmax - bounds$xmin) * 0.05
        y_padding <- (bounds$ymax - bounds$ymin) * 0.05
        map <- map %>%
          fitBounds(
            lng1 = bounds$xmin - x_padding,
            lat1 = bounds$ymin - y_padding,
            lng2 = bounds$xmax + x_padding,
            lat2 = bounds$ymax + y_padding
          )
      }
    } else {
      map <- map %>% setView(lng = 0, lat = 0, zoom = 2)
    }
    
    # Add legend
    map <- map %>%
      addLegend(
        position = "bottomright",
        colors = c("lightblue", "lightcoral"),
        labels = c("Present", 
          paste("Future", current_data$period_label, "- IC", current_data$ic)),
        title = paste0("<span style='font-size: 12px; white-space: nowrap;'><i>", gsub("_", " ", input$species), "</i></span>"),
        opacity = 0.8,
        layerId = "legend"
      )
    
    return(map)
  })
  
  # Update polygons when IC changes (Navigation Mode only)
  observeEvent(input$initial_condition, {
    req(input$species)
    
    current_data <- get_current_polygons()
    
    leafletProxy("map") %>%
      clearShapes() %>%
      removeControl("legend")
    
    if(!is.null(current_data$present_poly)) {
      leafletProxy("map") %>%
        addPolygons(
          data = current_data$present_poly,
          fillColor = "lightblue",
          fillOpacity = 0.3,
          color = "blue",
          weight = 2,
          layerId = "present_layer",
          smoothFactor = 0.5
        )
    }
    
    if(!is.null(current_data$future_poly)) {
      leafletProxy("map") %>%
        addPolygons(
          data = current_data$future_poly,
          fillColor = "lightcoral",
          fillOpacity = 0.3,
          color = "red",
          weight = 2,
          layerId = "future_layer",
          smoothFactor = 0.5
        )
    }
    
    leafletProxy("map") %>%
      addLegend(
        position = "bottomright",
        colors = c("lightblue", "lightcoral"),
        labels = c("Present", 
          paste("Future", current_data$period_label, "- IC", current_data$ic)),
        title = paste0("<span style='font-size: 12px; white-space: nowrap;'><i>", gsub("_", " ", input$species), "</i></span>"),
        opacity = 0.8,
        layerId = "legend"
      )
  }, ignoreInit = TRUE)
  
  # Update polygons when time period changes (Navigation Mode only)
  observeEvent(input$time_period, {
    req(input$species)
    
    current_data <- get_current_polygons()
    
    leafletProxy("map") %>%
      clearShapes() %>%
      removeControl("legend")
    
    if(!is.null(current_data$present_poly)) {
      leafletProxy("map") %>%
        addPolygons(
          data = current_data$present_poly,
          fillColor = "lightblue",
          fillOpacity = 0.3,
          color = "blue",
          weight = 2,
          layerId = "present_layer",
          smoothFactor = 0.5
        )
    }
    
    if(!is.null(current_data$future_poly)) {
      leafletProxy("map") %>%
        addPolygons(
          data = current_data$future_poly,
          fillColor = "lightcoral",
          fillOpacity = 0.3,
          color = "red",
          weight = 2,
          layerId = "future_layer",
          smoothFactor = 0.5
        )
    }
    
    leafletProxy("map") %>%
      addLegend(
        position = "bottomright",
        colors = c("lightblue", "lightcoral"),
        labels = c("Present", 
          paste("Future", current_data$period_label, "- IC", current_data$ic)),
        title = paste0("<span style='font-size: 12px; white-space: nowrap;'><i>", gsub("_", " ", input$species), "</i></span>"),
        opacity = 0.8,
        layerId = "legend"
      )
  }, ignoreInit = TRUE)
  
  # Unified Comparison Plot with SILHOUETTES
  output$unified_comparison_plot <- renderPlot({
    req(input$species, input$time_period, input$selected_metric)
    
    # Filter comprehensive_metrics by the selected time period
    filtered_metrics <- comprehensive_metrics[comprehensive_metrics$period == input$time_period, ]
    
    # Check if selected metric column exists
    if(!input$selected_metric %in% colnames(filtered_metrics)) {
      return(ggplot() + 
          geom_text(aes(x = 1, y = 1, label = paste(input$selected_metric, "column not found in data")), 
            size = 6) +
          theme_void())
    }
    
    # Prepare data for plotting
    plot_data <- data.frame()
    x_pos <- 1
    species_labels <- c()
    species_names_for_silhouettes <- c()
    
    # Add main species data
    if("species" %in% colnames(filtered_metrics)) {
      main_species_data <- filtered_metrics[filtered_metrics$species == input$species, ]
    } else {
      main_species_data <- filtered_metrics
    }
    
    if(nrow(main_species_data) > 0) {
      main_data <- data.frame(
        species = gsub("_", " ", input$species),
        metric_value = main_species_data[[input$selected_metric]],
        x_pos = x_pos,
        species_type = "Current Species",
        alpha_value = 0.5,
        original_species_name = input$species
      )
      plot_data <- rbind(plot_data, main_data)
      species_labels <- c(species_labels, gsub("_", " ", input$species))
      species_names_for_silhouettes <- c(species_names_for_silhouettes, input$species)
      x_pos <- x_pos + 1
    }
    
    # Add comparison species data
    comparison_inputs <- c(input$comparison_species_1, input$comparison_species_2, 
      input$comparison_species_3, input$comparison_species_4,
      input$comparison_species_5)
    
    colors <- c("Current Species" = "steelblue", 
      "Comparison Species 1" = "orange", 
      "Comparison Species 2" = "green", 
      "Comparison Species 3" = "purple", 
      "Comparison Species 4" = "red",
      "Comparison Species 5" = "brown")
    
    for(i in 1:5) {
      if(is.null(comparison_inputs[i]) || comparison_inputs[i] == "") {
        comp_data <- data.frame(
          species = paste("None", i),
          metric_value = 0,
          x_pos = x_pos,
          species_type = paste("Comparison Species", i),
          alpha_value = 0,
          original_species_name = NA
        )
        plot_data <- rbind(plot_data, comp_data)
        species_labels <- c(species_labels, paste("None", i))
        species_names_for_silhouettes <- c(species_names_for_silhouettes, NA)
        x_pos <- x_pos + 1
      } else {
        comp_species_data <- filtered_metrics[filtered_metrics$species == comparison_inputs[i], ]
        
        if(nrow(comp_species_data) > 0) {
          comp_data <- data.frame(
            species = gsub("_", " ", comparison_inputs[i]),
            metric_value = comp_species_data[[input$selected_metric]],
            x_pos = x_pos,
            species_type = paste("Comparison Species", i),
            alpha_value = 0.3,
            original_species_name = comparison_inputs[i]
          )
          plot_data <- rbind(plot_data, comp_data)
          species_labels <- c(species_labels, gsub("_", " ", comparison_inputs[i]))
          species_names_for_silhouettes <- c(species_names_for_silhouettes, comparison_inputs[i])
          x_pos <- x_pos + 1
        }
      }
    }
    
    # Check if we have data to plot
    if(nrow(plot_data) == 0) {
      return(ggplot() + 
          geom_text(aes(x = 1, y = 1, label = paste("No data available for", input$time_period)), 
            size = 6) +
          theme_void())
    }
    
    # Get metric display name
    metric_names <- list(
      "area_pct_diff" = "Area Change (%)",
      "percentage_expansion" = "Gained Area (%)",
      "percentage_lost" = "Lost Area (%)",
      "fragmentation_pct_diff" = "Fragmentation Change (%)",
      "fractal_pct_diff" = "Fractal Dimension Change (%)",
      "pa_ratio_pct_diff" = "Perimeter Area Ratio Change (%)",
      "suitability_increase_ratio" = "Suitability Ratio Change (%)"
    )
    
    metric_display_name <- metric_names[[input$selected_metric]]
    if(is.null(metric_display_name)) metric_display_name <- input$selected_metric
    
    # Calculate y-axis limits
    real_data <- plot_data[plot_data$alpha_value > 0.2, ]
    if(nrow(real_data) > 0) {
      y_range <- range(real_data$metric_value, na.rm = TRUE)
    } else {
      y_range <- c(-1, 1)
    }
    y_padding <- max(0.1, (y_range[2] - y_range[1]) * 0.1)
    
    # Special handling for gained/lost area metrics (start at 0, go to max)
    if(input$selected_metric %in% c("percentage_expansion", "percentage_lost")) {
      y_limits <- c(0, max(y_range[2] + y_padding, 5))
    } else if(input$selected_metric %in% c("area_pct_diff", "fragmentation_pct_diff", "fractal_pct_diff", "pa_ratio_pct_diff", "suitability_increase_ratio")) {
      y_limits <- c(min(y_range[1] - y_padding, -5, 0), max(y_range[2] + y_padding, 5, 0))
    } else {
      y_limits <- c(y_range[1] - y_padding, y_range[2] + y_padding)
    }
    
    # Create comparison plot
    p <- ggplot(plot_data, aes(x = x_pos, y = metric_value, color = species_type, alpha = alpha_value)) +
      geom_point(size = 5, position = position_jitter(width = 0)) +
      coord_cartesian(ylim = y_limits) +
      scale_x_continuous(breaks = 1:length(species_labels), labels = species_labels) +
      scale_color_manual(values = colors[names(colors) %in% unique(plot_data$species_type)]) +
      scale_alpha_identity() +
      labs(
        x = "Species",
        y = metric_display_name,
        color = "Species Type"
      ) +
      theme_minimal() +
      theme(
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, face="italic"),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
      )
    
    # Add horizontal line at 0 for percentage-based metrics
    if(input$selected_metric %in% c("area_pct_diff", "percentage_expansion", "percentage_lost", 
      "fragmentation_pct_diff", "fractal_pct_diff", "pa_ratio_pct_diff", "suitability_increase_ratio")) {
      p <- p + geom_hline(yintercept = 0, color = "red", linewidth = 0.5) +
        annotate("text", x = length(species_labels)/2 + 0.5, y = 0, label = "No Change", 
          vjust = -0.5, hjust = 0.5, color = "red", size = 5, fontface = "bold")
    }
    
    
    # ADD SILHOUETTES
    silhouette_y <- min(y_limits) + (max(y_limits) - min(y_limits)) * 0.05  # Position near bottom
    silhouette_height <- (max(y_limits) - min(y_limits)) * 0.15  # Height as 15% of y-axis range
    
    for(i in 1:length(species_names_for_silhouettes)) {
      species_name <- species_names_for_silhouettes[i]
      
      if(!is.na(species_name)) {
        sil_info <- silhouette_data[silhouette_data$species == species_name, ]
        
        if(nrow(sil_info) > 0 && sil_info$has_silhouette && !is.na(sil_info$uuid)) {
          p <- p + add_phylopic(uuid = sil_info$uuid, 
            x = i, 
            y = silhouette_y,
            height = silhouette_height,  # Now proportional instead of fixed
            width = NA,
            alpha = 0.2, 
            fill = "black")
        }
      }
    }
    
    return(p)
  })
}

# Run the application
cat("Starting Species Distribution Shiny App with Navigation and Animation Modes...\n")

shinyApp(ui = ui, server = server)