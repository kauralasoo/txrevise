% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constructAlternativeEvents.R
\name{identifyTranscriptGroups}
\alias{identifyTranscriptGroups}
\title{Identify two groups of overlapping transcripts that share the most exons between each other}
\usage{
identifyTranscriptGroups(granges_list)
}
\arguments{
\item{granges_list}{GRanges list of transcripts}
}
\value{
A list of lists containing transcripts in the two groups
}
\description{
Identify two groups of overlapping transcripts that share the most exons between each other
}
