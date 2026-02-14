#Deploy the app

library(rsconnect)

#AUTHORIZE ACCOUNT
rsconnect::setAccountInfo(name='felipeben',
  token='7FEC02E1D3AE87DEA366770328A77C72',
  secret='CX2RTj1bZeyljhZwKQmUnyFHmT6fFOsFN1WRXMhD')

#DEPLOY APP with everything
rsconnect::deployApp(
  appName = "species-distribution-app",
  appTitle = "Do Initial Conditions Matter for Species Distribution Models?"
)


# Deploy core app without animations
rsconnect::deployApp(
  appName = "species-distribution-core",
  appFiles = c("app.R", "nested_list_polygons.rds", "comprehensive_metrics.csv", "species_silhouettes_data.rds"),
  forceUpdate = TRUE
)