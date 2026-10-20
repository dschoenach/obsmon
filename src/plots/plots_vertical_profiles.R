stationVerticalProfilePlottingFunction <- function(plot) {
  sqliteParams <- plot$paramsAsInSqliteDbs
  hasComp <- .plotHasComparison(plot)
  metadataPlot <- .getFirstNonEmptyComparisonPlot(plot)

  scaleColors <- c(
    "obsvalue"="black", "obsvalue_corrected"="black", "obsvalue_raw"="blue",
    "fg_dep"="blue", "first_guess"="blue",
    "an_dep"="red", "analysis"="red"
  )

  if(hasComp) {
    plotData <- .combinePlotDataByExperiment(plot)
    plotData <- plotData[
      c("level", "experiment", intersect(colnames(plotData), names(scaleColors)))
    ]

    localPlotData <- plotData %>%
      group_by(level, experiment) %>%
      melt(id=c("level", "experiment"), na.rm=TRUE)
    plotDataSummary <- plotData %>%
      group_by(level, experiment) %>%
      gather(key="variable", value="value", -level, -experiment) %>%
      group_by(level, experiment, variable) %>%
      summarize_all(list(~mean(.), ~min(.), ~max(.)))

    experimentColours <- .getExperimentColourValues(localPlotData$experiment)
    legendColours <- setNames(
      experimentColours[as.character(localPlotData$experiment)],
      paste(localPlotData$experiment, localPlotData$variable, sep=" :: ")
    )
    variableLinetypes <- .getVariableLineTypeValues(localPlotData$variable)

    obplot <- ggplot() +
      aes(
        x=level,
        group=paste(experiment, variable, sep=" :: "),
        colour=paste(experiment, variable, sep=" :: "),
        shape=variable,
        linetype=variable
      ) +
      scale_colour_manual(name="Experiment", values=legendColours) +
      scale_linetype_manual(name=NULL, values=variableLinetypes) +
      geom_point(data=localPlotData, aes(y=value), size=.getComparisonPointSize(2)) +
      geom_line(data=plotDataSummary, aes(y=mean), linewidth=.getComparisonLineWidth()) +
      geom_linerange(
        data=plotDataSummary, aes(ymin=min, ymax=max),
        alpha=0.25
      )
  } else {
    plotData <- plot$data
    plotData <- plotData[
      c("level", intersect(colnames(plotData), names(scaleColors)))
    ]

    localPlotData <- plotData %>%
      group_by(level) %>%
      melt(id="level", na.rm=TRUE)
    plotDataSummary <- plotData %>%
      group_by(level) %>%
      gather(key="variable", value="value", -level) %>%
      group_by(level, variable) %>%
      summarize_all(list(~mean(.), ~min(.), ~max(.)))

    obplot <- ggplot() +
      aes(x=level, group=variable, colour=variable, shape=variable) +
      scale_colour_manual(values=scaleColors) +
      # Scatter plot with all data points
      geom_point(data=localPlotData, aes(y=value), size=2) +
      # Draw line through mean values instead of actual data point values
      geom_line(data=plotDataSummary, aes(y=mean), size=1.15) +
      # For clarity: Connect groups of points belonging to the same level
      geom_linerange(
        data=plotDataSummary, aes(ymin=min, ymax=max),
        linetype="dashed", alpha=0.25
      )
  }

  # Axes' labels
  xlab <- sprintf("Level [%s]", units(metadataPlot$dataWithUnits$level))
  ylab <- sqliteParams$varname
  for (colname in names(scaleColors)) {
    if(colname %in% localPlotData$variable) {
      ylab <- sprintf("%s [%s]", ylab, units(metadataPlot$dataWithUnits[[colname]]))
      break
    }
  }
  obplot <- obplot + labs(x=xlab, y=ylab)

  # Set x-axis limits and flipping axes
  # The change in xlim causes ggplotly to omit either the shape or
  # the line (whichever is added later) from the legend. Not a big
  # deal but worth pointing out, as this but may be solved in later
  # releases of the plotly package
  xlim <- c()
  if(ud_are_convertible(units(metadataPlot$dataWithUnits$level), "Pa")) {
    units(refPressures) <- units(metadataPlot$dataWithUnits$level)
    xlim <- c(max(drop_units(refPressures),localPlotData[["level"]]), 0)
    obplot <- obplot + scale_x_reverse()
  } else if(ud_are_convertible(units(metadataPlot$dataWithUnits$level), "meters")) {
    units(refHeights) <- units(metadataPlot$dataWithUnits$level)
    xlim <- c(0, max(drop_units(refHeights), localPlotData[["level"]]))
  }
  obplot <- obplot + coord_flip_wrapper(default=TRUE, xlim=xlim)

  # Hide legend if there's only one curve to be plotted
  if(!hasComp && length(unique(localPlotData[["variable"]]))<2) {
    obplot <- obplot + theme(legend.position = "none")
  }

  return(obplot)
}

# Plot 1
plotRegistry$registerPlotType(
  name="Station Vertical Profile: Obsvalue",
  category="StationVerticalProfiles",
  dataFieldsInRetrievedPlotData=list("level", "obsvalue"),
  dataFieldsInSqliteWhereClause=list(
    "statid", "obnumber", "obname", "varname"
  ),
  stationChoiceType="single",
  plottingFunction=stationVerticalProfilePlottingFunction
)

# Plot 2
plotRegistry$registerPlotType(
  name="Station Vertical Profile: Obs, FG & Analysis",
  category="StationVerticalProfiles",
  dataFieldsInRetrievedPlotData=list("level", "obsvalue", "fg_dep", "an_dep"),
  dataFieldsInSqliteWhereClause=list("statid", "obnumber", "obname", "varname"),
  stationChoiceType="single",
  plottingFunction=stationVerticalProfilePlottingFunction,
  dataPostProcessingFunction = function(data, ...) {
    data$first_guess <- data$obsvalue - data$fg_dep
    data$analysis <- data$obsvalue - data$an_dep
    data <- subset(data, select=-c(fg_dep, an_dep))
    return(data)
  }
)

# Plot 3
plotRegistry$registerPlotType(
  name="Station Vertical Profile: FG & Analysis Departure",
  category="StationVerticalProfiles",
  dataFieldsInRetrievedPlotData=list("level", "fg_dep", "an_dep"),
  dataFieldsInSqliteWhereClause=list("statid", "obnumber", "obname", "varname"),
  stationChoiceType="single",
  plottingFunction=stationVerticalProfilePlottingFunction
)

# Plot 4
plotRegistry$registerPlotType(
  name="Station Vertical Profile: Bias",
  category="StationVerticalProfiles",
  dataFieldsInRetrievedPlotData=list("level", "obsvalue", "biascrl"),
  dataFieldsInSqliteWhereClause=list("statid", "obnumber", "obname", "varname"),
  stationChoiceType="single",
  plottingFunction=stationVerticalProfilePlottingFunction,
  dataPostProcessingFunction = function(data, ...) {
    colnames(data)[colnames(data) == 'obsvalue'] <- 'obsvalue_corrected'
    data$obsvalue_raw <- data$obsvalue_corrected + data$biascrl
    return(data)
  }
)
