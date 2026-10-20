# Dockerfile for local reproducible testing of multiDetectR.
# Based on rocker/r-ver (lightweight; no spatial libs needed).
#
# NOTE: CI (.gitlab-ci.yml) does NOT use this image. CI runs directly on the
# public rocker/r-ver image and installs dependencies at job time, to avoid
# maintaining a project registry image. This Dockerfile is kept for local use:
#   docker build -t multidetectr-ci .
#   docker run --rm -v "$(pwd):/pkg" multidetectr-ci

FROM rocker/r-ver:4.6.1

LABEL maintainer="rerickson@usgs.gov"
LABEL description="CI/CD image for the multiDetectR R package"

# System dependencies for building R packages and rendering vignettes
RUN apt-get update -qq && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    pandoc \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Layer 1: Package Suggests (Imports are only base R: graphics, stats)
RUN Rscript -e " \
  install.packages(c( \
    'shiny', 'DT', 'testthat', 'knitr', 'rmarkdown' \
  ), repos = 'https://cloud.r-project.org', Ncpus = 4)"

# Layer 2: Development and CI tools
RUN Rscript -e " \
  install.packages(c( \
    'devtools', 'rcmdcheck', 'remotes' \
  ), repos = 'https://cloud.r-project.org', Ncpus = 4)"

# Disable forcing Suggests during check
ENV _R_CHECK_FORCE_SUGGESTS_=false
ENV _R_CHECK_NO_STOP_ON_TEST_ERROR_=true

# Set working directory
WORKDIR /pkg

# Default command: run tests
CMD ["R", "-e", "devtools::test()"]
