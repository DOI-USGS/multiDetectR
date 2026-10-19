# multiDetectR Shiny front end
#
# A graphical wrapper around multiDetectR::allocate() so non-R users can apply
# the closed-form allocation rule. Per the target journal's requirement that a
# graphical interface be able to export a reproducible script, the app includes
# a "Download R script" button that writes a standalone .R file reproducing the
# current design.
#
# Launch with multiDetectR::run_app(), or from this directory with
# shiny::runApp().

library(shiny)
library(DT)
library(multiDetectR)

max_levels <- 8L
max_targets <- 20L

## USGS Applications header and footer templates (common.css/custom.css and the
## logo live under www/, which Shiny serves at the app root). Read them as HTML
## and inject around the main content, following the USGS age-subsampling-
## simulation project's layout. runApp() sets the working directory to the app
## directory, so these relative paths resolve.
header_content <- paste(readLines("header.txt", warn = FALSE), collapse = "\n")
footer_content <- paste(readLines("footer.txt", warn = FALSE), collapse = "\n")

ui <- fluidPage(
  tags$style(type = "text/css", "
    .container-fluid { padding: 0; }
    html, body { margin: 0; padding: 0; }
    #floating-footer {
      position: relative;
      bottom: 0;
      width: 100%;
      margin: 0;
      padding: 0;
      z-index: 1000;
    }
    .main-content { padding: 15px 30px; }
  "),
  htmltools::HTML(header_content),
  div(class = "main-content",
  titlePanel("multiDetectR: allocate detection effort across nested survey levels"),
  tabsetPanel(
    id = "main_tabs",

    ## ------------------------------------------------------------------ ##
    ## Tab 0: overview                                                    ##
    ## ------------------------------------------------------------------ ##
    tabPanel(
    "Overview",
    br(),
    div(class = "about-content",
        style = "max-width: 850px; line-height: 1.6;",
      h3("What this tool does"),
      p(HTML(paste0(
        "A survey that searches for a target &mdash; a species, an infectious ",
        "agent, a contaminant &mdash; often nests its sampling: sites hold ",
        "water samples, samples hold replicates, and so on. Each level has its ",
        "own cost and its own chance of carrying the target, and a fixed budget ",
        "has to be split among the levels. This tool implements the closed-form ",
        "rule that solves that allocation exactly, for any number of levels."))),
      h4("The two numbers per level"),
      p(HTML(paste0(
        "Levels are numbered from the outermost (reporting) level 1 to the ",
        "innermost level <i>d</i>. For each level <i>i</i> you supply:"))),
      tags$ul(
        tags$li(HTML(paste0(
          "<b>&theta;<sub><i>i</i></sub></b>, the conditional detectability: ",
          "the probability a level-<i>i</i> unit carries the target given that ",
          "the unit containing it does; and"))),
        tags$li(HTML(paste0(
          "<b><i>c</i><sub><i>i</i></sub></b>, the marginal cost of one ",
          "level-<i>i</i> unit, excluding the cost of whatever it contains.")))
      ),
      p(HTML(paste0(
        "Given a budget (or a required detection probability), the rule returns ",
        "<b><i>n</i><sub><i>i</i></sub></b>, how many level-<i>i</i> units to ",
        "take within each unit of the level above &mdash; the recommended ",
        "sampling design."))),
      h4("The key ideas"),
      tags$ul(
        tags$li(HTML(paste0(
          "<b>Detection on a log scale.</b> Detection probabilities multiply, ",
          "but their logs (nats) add, so the design problem separates into one ",
          "simple problem per level, solved from the innermost level outward."))),
        tags$li(HTML(paste0(
          "<b>One number summarizes a hierarchy.</b> At the optimum the whole ",
          "design collapses to &kappa;<sub>1</sub>, the price of a detection at ",
          "the outermost level, which sets the budget&ndash;detection curve."))),
        tags$li(HTML(paste0(
          "<b>Invest in cheap inner levels.</b> Optimal designs usually put ",
          "most replication where it is cheapest; optimizing barely matters at ",
          "two levels but matters more as a design deepens.")))
      ),
      h4("Using the tabs"),
      tags$ul(
        tags$li(HTML(paste0(
          "<b>Single species.</b> Set the number of levels, each level's ",
          "&theta; and cost, and a budget or a target detection probability. ",
          "You get the per-level design, a budget&ndash;detection curve you can ",
          "hover to read values, and an optimality audit. One click exports a ",
          "reproducible R script."))),
        tags$li(HTML(paste0(
          "<b>Multispecies assemblage.</b> Several targets share one survey but ",
          "detect differently. Enter each target's &theta; across the levels ",
          "(edit the table or upload a CSV) and a relative weight. The tool ",
          "returns the single design that best serves the whole set, which is ",
          "governed by the hardest-to-detect targets, not the average one.")))
      ),
      p(HTML(paste0(
        "<b>Reference:</b> Lafferty KD, Larson DM, Spear SF, Schulte N, Coombs ",
        "JA, DeHaan P, Hayer C-A, Maloy AP, Erickson RA. An exact rule ",
        "allocates detection effort across nested survey levels at any depth: ",
        "an R package and eDNA examples. <i>Methods in Ecology and Evolution</i> ",
        "(in preparation). See also <code>vignette(\"multiDetectR\")</code> and ",
        "<code>vignette(\"simulated-example\")</code> in the R package."))),
      p(em(paste0(
        "This software is preliminary or provisional and is subject to ",
        "revision. See the disclaimer in the footer.")))
    )
    ), # end Overview tab

    ## ------------------------------------------------------------------ ##
    ## Tab 1: single species                                              ##
    ## ------------------------------------------------------------------ ##
    tabPanel(
    "Single species",
    br(),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      helpText(
        HTML(paste0(
          "Levels run from the outermost (reporting) level 1 to the innermost ",
          "level <i>d</i>. For each level <i>i</i>, set the conditional ",
          "detectability &theta;<sub><i>i</i></sub> and the marginal unit cost ",
          "<i>c</i><sub><i>i</i></sub>."))
      ),
      sliderInput("d", "Number of levels (depth)", min = 1L, max = max_levels,
                  value = 3L, step = 1L),
      uiOutput("level_inputs"),
      radioButtons("goal", "Goal",
                   choices = c("Maximize detection" = "detection",
                               "Maximize count" = "count"),
                   selected = "detection"),
      radioButtons("constraint", "Constraint",
                   choices = c("Budget" = "budget",
                               "Required detection probability" = "P_target"),
                   selected = "budget"),
      conditionalPanel(
        "input.constraint == 'budget'",
        numericInput("budget", "Budget", value = 5000, min = 0, step = 100)
      ),
      conditionalPanel(
        "input.constraint == 'P_target'",
        numericInput("P_target", "Required detection probability",
                     value = 0.95, min = 0.01, max = 0.999, step = 0.01)
      ),
      checkboxInput("integer", "Integer replicate counts", value = TRUE),
      checkboxInput("enumerate",
                    "Use exhaustive search (small budgets)", value = FALSE),
      hr(),
      downloadButton("download_script", "Download R script")
    ),
    mainPanel(
      width = 8,
      verbatimTextOutput("headline"),
      h4("Design by level"),
      helpText(HTML(paste0(
        "This table gives the optimized design one row per level, outermost ",
        "(level 1) to innermost. For each level <i>i</i> it shows the inputs ",
        "you set &mdash; the conditional detectability &theta;<sub><i>i</i></sub> ",
        "and the unit cost <i>c</i><sub><i>i</i></sub> &mdash; alongside the ",
        "allocation the rule returns: <i>n</i><sub><i>i</i></sub>, the number of ",
        "level-<i>i</i> units to take within each unit of the level above; ",
        "<i>h</i><sub><i>i</i></sub>, the detection hazard in nats that one such ",
        "unit contributes; and &kappa;<sub><i>i</i></sub>, the price of a ",
        "detection at that level (cost per nat). Reading the <i>n</i> column ",
        "gives the recommended sampling design; &kappa;<sub>1</sub> at the ",
        "outermost level summarizes the whole hierarchy and drives the ",
        "budget&ndash;detection curve below."))),
      tableOutput("design_table"),
      h4("Budget-detection frontier"),
      helpText("Hover over the curve to read the budget and the detection ",
               "probability at that point."),
      plotOutput("frontier", height = "320px",
                 hover = hoverOpts(id = "frontier_hover", delay = 50,
                                   delayType = "throttle")),
      verbatimTextOutput("frontier_readout"),
      h4("Optimality audit (A = B at the optimum)"),
      helpText(
        "At an optimal design, each level spends the same share of its budget ",
        "on the units nested inside it (B) as those units contribute to the ",
        "level's detection power (A). The table reports A and B for every ",
        "level that contains another; at the continuous optimum A - B is zero ",
        "(up to rounding for integer counts), so large gaps flag a level that ",
        "is over- or under-replicated."),
      tableOutput("identity_table")
    )
  )
  ), # end Tab 1

    ## ------------------------------------------------------------------ ##
    ## Tab 2: multispecies assemblage                                     ##
    ## ------------------------------------------------------------------ ##
    tabPanel(
    "Multispecies assemblage",
    br(),
    sidebarLayout(
      sidebarPanel(
        width = 4,
        helpText(HTML(paste0(
          "Several targets (e.g., species) share one survey but detect ",
          "differently. Enter each target's conditional detectability ",
          "&theta;<sub><i>i</i></sub> at every level in the table, plus a ",
          "relative <i>weight</i> (how much you care about that target; equal ",
          "by default). The rule returns the single design that best serves the ",
          "whole set, which is governed by the hardest-to-detect targets."))),
        sliderInput("m_d", "Number of levels (depth)",
                    min = 1L, max = max_levels, value = 3L, step = 1L),
        sliderInput("m_targets", "Number of targets",
                    min = 1L, max = max_targets, value = 2L, step = 1L),
        helpText(HTML(paste0(
          "Edit the &theta; table below (double-click a cell). Per-level unit ",
          "costs <i>c</i><sub><i>i</i></sub> are shared across targets and set ",
          "in the cost row."))),
        numericInput("m_budget", "Budget", value = 5000, min = 0, step = 100),
        checkboxInput("m_integer", "Integer replicate counts", value = TRUE),
        hr(),
        helpText(HTML(paste0(
          "Or upload a CSV with columns <code>target, theta_1, ..., theta_d, ",
          "weight</code>. Uploading replaces the table."))),
        fileInput("m_csv", NULL, accept = ".csv"),
        downloadButton("m_template", "Download CSV template"),
        tags$br(), tags$br(),
        downloadButton("m_download_script", "Download R script")
      ),
      mainPanel(
        width = 8,
        h4("Per-target detectabilities (editable)"),
        DT::DTOutput("m_theta_table"),
        tags$br(),
        h4("Shared per-level unit costs"),
        DT::DTOutput("m_cost_table"),
        tags$br(),
        verbatimTextOutput("m_headline"),
        h4("Shared assemblage design"),
        helpText(HTML(paste0(
          "The single design that serves all targets, one row per level. ",
          "<i>n</i><sub><i>i</i></sub> is the number of level-<i>i</i> units to ",
          "take within each unit above."))),
        tableOutput("m_design_table"),
        h4("Detection by target"),
        helpText(HTML(paste0(
          "Each target's detection probability under the shared design, with ",
          "its miss-weight &lambda; = weight &times; (1 &minus; detection). The ",
          "target with the largest &lambda; (marked) binds the design: effort ",
          "is pulled toward the elusive targets, not the average one."))),
        tableOutput("m_target_table")
      )
    )
    ) # end Tab 2

  ) # end tabsetPanel
  ), # end main-content div
  div(id = "floating-footer",
    htmltools::HTML(footer_content)
  )
)

server <- function(input, output, session) {

  ## Dynamically generated per-level theta and cost inputs.
  output$level_inputs <- renderUI({
    d <- input$d
    defaults_theta <- rep(0.5, d)
    defaults_cost <- 10^(seq(d, 1) - 1) * 2   # outer levels cost more
    lapply(seq_len(d), function(i) {
      fluidRow(
        column(6, numericInput(
          paste0("theta_", i),
          HTML(sprintf("&theta;<sub>%d</sub>", i)),
          value = defaults_theta[i], min = 0.001, max = 0.999, step = 0.01)),
        column(6, numericInput(
          paste0("cost_", i),
          HTML(sprintf("<i>c</i><sub>%d</sub>", i)),
          value = defaults_cost[i], min = 0.0001, step = 1))
      )
    })
  })

  ## Collect the per-level vectors from the dynamic inputs.
  inputs <- reactive({
    d <- input$d
    theta <- vapply(seq_len(d),
                    function(i) input[[paste0("theta_", i)]], numeric(1))
    cost <- vapply(seq_len(d),
                   function(i) input[[paste0("cost_", i)]], numeric(1))
    list(theta = theta, cost = cost)
  })

  ## Run the allocation, guarding against transient NULL inputs and bad values.
  design <- reactive({
    iv <- inputs()
    req(length(iv$theta) == input$d, !any(is.na(iv$theta)),
        !any(is.na(iv$cost)))
    validate(
      need(all(iv$theta > 0 & iv$theta < 1),
           "Each theta must be in (0, 1)."),
      need(all(iv$cost > 0), "Each cost must be positive.")
    )
    args <- list(theta = iv$theta, cost = iv$cost, goal = input$goal,
                 integer = input$integer)
    if (input$constraint == "budget") {
      validate(need(!is.na(input$budget) && input$budget > 0,
                    "Budget must be positive."))
      args$budget <- input$budget
    } else {
      if (input$goal == "count") {
        validate(need(FALSE,
                      "A required detection probability applies to the ",
                      "detection goal; switch the goal or use a budget."))
      }
      args$P_target <- input$P_target
    }

    tryCatch({
      if (isTRUE(input$enumerate) && input$constraint == "budget") {
        do.call(enumerate_designs,
                args[c("theta", "cost", "budget", "goal")])
      } else {
        do.call(allocate, args)
      }
    }, error = function(e) {
      validate(need(FALSE, paste("Could not solve:", conditionMessage(e))))
    })
  })

  output$headline <- renderText({
    des <- design()
    sprintf(paste0("Detection probability: %.4f\nTotal cost: %.4g\n",
                   "Price of a detection (\u03ba\u2081): %.4g\n",
                   "Goal: %s | %s counts"),
            des$P1, des$C, des$kappa1, des$goal,
            if (des$integer) "integer" else "continuous")
  })

  output$design_table <- renderTable({
    s <- summary(design())
    ## Relabel columns with Greek symbols and subscripts for display.
    data.frame(
      "Level"            = s$level,
      "&theta;<sub>i</sub> (detectability)" = s$theta,
      "<i>c</i><sub>i</sub> (unit cost)"    = s$cost,
      "<i>n</i><sub>i</sub> (count)"        = s$n,
      "<i>h</i><sub>i</sub> (hazard, nats)" = s$hazard,
      "&kappa;<sub>i</sub> (price/nat)"     = s$kappa,
      check.names = FALSE
    )
  }, digits = 4, sanitize.text.function = identity,
     sanitize.colnames.function = identity)

  output$frontier <- renderPlot({
    plot(design())
  })

  ## Readout of the point under the cursor. The x position is a budget; the
  ## detection probability on the curve is 1 - exp(-budget / kappa_1), eq (10),
  ## so this reports exactly the value the curve shows at the cursor.
  output$frontier_readout <- renderText({
    hv <- input$frontier_hover
    if (is.null(hv) || is.null(hv$x)) {
      return("Hover over the curve to read a value.")
    }
    des <- design()
    budget <- max(hv$x, 0)
    detection <- detection_from_budget(des$kappa1, budget)
    sprintf("Budget = %.4g   ->   detection probability = %.4f",
            budget, detection)
  })

  output$identity_table <- renderTable({
    id <- optimality_identity(design())
    if (nrow(id) == 0) {
      data.frame(note = "Single-level design: no inner levels to audit.")
    } else {
      id
    }
  }, digits = 4)

  ## Export a standalone, reproducible R script for the current design.
  output$download_script <- downloadHandler(
    filename = function() "multidetectr_design.R",
    content = function(file) {
      iv <- inputs()
      writeLines(build_script(iv$theta, iv$cost, input$goal, input$constraint,
                              input$budget, input$P_target, input$integer,
                              input$enumerate), file)
    }
  )

  ## =================================================================== ##
  ## Tab 2: multispecies assemblage                                      ##
  ## =================================================================== ##

  ## Editable state: a theta matrix (targets x levels + weight) and a shared
  ## cost row, held in reactiveVals so edits, slider changes, and CSV uploads
  ## all update one source of truth.
  m_theta <- reactiveVal(NULL)
  m_cost  <- reactiveVal(NULL)

  ## (Re)build default tables when the depth or target-count sliders change,
  ## preserving overlapping values where possible.
  observeEvent(list(input$m_d, input$m_targets), {
    d <- input$m_d
    nt <- input$m_targets
    old <- m_theta()
    th <- matrix(0.5, nrow = nt, ncol = d)
    wt <- rep(1, nt)
    if (!is.null(old)) {
      ro <- min(nrow(old$theta), nt)
      co <- min(ncol(old$theta), d)
      th[seq_len(ro), seq_len(co)] <- old$theta[seq_len(ro), seq_len(co),
                                                 drop = FALSE]
      wt[seq_len(ro)] <- old$weight[seq_len(ro)]
    }
    rownames(th) <- paste0("target_", seq_len(nt))
    m_theta(list(theta = th, weight = wt))

    oldc <- m_cost()
    cst <- 10^(seq(d, 1) - 1) * 2
    if (!is.null(oldc)) {
      co <- min(length(oldc), d)
      cst[seq_len(co)] <- oldc[seq_len(co)]
    }
    m_cost(cst)
  }, ignoreInit = FALSE)

  ## Render the editable theta + weight table.
  output$m_theta_table <- DT::renderDT({
    st <- m_theta()
    req(st)
    df <- data.frame(target = rownames(st$theta),
                     round(st$theta, 3),
                     weight = st$weight,
                     check.names = FALSE)
    ## Greek + subscript column headers (HTML; escape = FALSE below renders them).
    theta_hdr <- paste0("&theta;<sub>", seq_len(ncol(st$theta)), "</sub>")
    colnames(df) <- c("target", theta_hdr, "<i>w</i> (weight)")
    DT::datatable(df, rownames = FALSE, editable = "cell", escape = FALSE,
                  options = list(dom = "t", paging = FALSE, ordering = FALSE))
  }, server = TRUE)

  ## Apply edits to the theta table back into state.
  observeEvent(input$m_theta_table_cell_edit, {
    info <- input$m_theta_table_cell_edit
    st <- m_theta()
    j <- info$col            # 0-based: 0 = target name, 1..d = theta, d+1 = weight
    d <- ncol(st$theta)
    val <- suppressWarnings(as.numeric(info$value))
    if (j == 0) {
      rownames(st$theta)[info$row] <- as.character(info$value)
    } else if (j >= 1 && j <= d) {
      if (!is.na(val)) st$theta[info$row, j] <- min(max(val, 1e-3), 1 - 1e-3)
    } else if (j == d + 1) {
      if (!is.na(val) && val >= 0) st$weight[info$row] <- val
    }
    m_theta(st)
  })

  ## Render the editable shared-cost row.
  output$m_cost_table <- DT::renderDT({
    cst <- m_cost()
    req(cst)
    df <- as.data.frame(as.list(cst))
    ## Greek/italic + subscript headers for the shared per-level costs.
    colnames(df) <- paste0("<i>c</i><sub>", seq_along(cst), "</sub>")
    DT::datatable(df, rownames = FALSE, editable = "cell", escape = FALSE,
                  options = list(dom = "t", paging = FALSE, ordering = FALSE))
  }, server = TRUE)

  observeEvent(input$m_cost_table_cell_edit, {
    info <- input$m_cost_table_cell_edit
    cst <- m_cost()
    val <- suppressWarnings(as.numeric(info$value))
    if (!is.na(val) && val > 0) cst[info$col + 1] <- val
    m_cost(cst)
  })

  ## CSV upload: columns target, theta_1..theta_d, weight. Replaces the tables.
  observeEvent(input$m_csv, {
    df <- tryCatch(utils::read.csv(input$m_csv$datapath, check.names = FALSE),
                   error = function(e) NULL)
    validate(need(!is.null(df), "Could not read the CSV file."))
    theta_cols <- grep("^theta_", names(df), value = TRUE)
    validate(need(length(theta_cols) >= 1,
                  "CSV needs columns theta_1 ... theta_d."))
    th <- as.matrix(df[, theta_cols, drop = FALSE])
    storage.mode(th) <- "double"
    rn <- if ("target" %in% names(df)) as.character(df$target) else
      paste0("target_", seq_len(nrow(df)))
    rownames(th) <- rn
    wt <- if ("weight" %in% names(df)) as.numeric(df$weight) else
      rep(1, nrow(df))
    ## Sync the sliders to the uploaded dimensions, then set state.
    updateSliderInput(session, "m_d", value = ncol(th))
    updateSliderInput(session, "m_targets", value = nrow(th))
    m_theta(list(theta = th, weight = wt))
  })

  ## Downloadable CSV template matching the current dimensions.
  output$m_template <- downloadHandler(
    filename = function() "multidetectr_targets_template.csv",
    content = function(file) {
      d <- input$m_d
      nt <- input$m_targets
      df <- data.frame(target = paste0("target_", seq_len(nt)))
      for (i in seq_len(d)) df[[paste0("theta_", i)]] <- 0.5
      df$weight <- 1
      utils::write.csv(df, file, row.names = FALSE)
    }
  )

  ## Run the assemblage allocation.
  m_design <- reactive({
    st <- m_theta()
    cst <- m_cost()
    req(st, cst, ncol(st$theta) == length(cst))
    validate(
      need(all(st$theta > 0 & st$theta < 1),
           "Every theta must be in (0, 1)."),
      need(all(cst > 0), "Every cost must be positive."),
      need(all(st$weight >= 0) && any(st$weight > 0),
           "Weights must be non-negative and not all zero."),
      need(!is.na(input$m_budget) && input$m_budget > 0,
           "Budget must be positive.")
    )
    tryCatch(
      allocate_assemblage(st$theta, cst, budget = input$m_budget,
                          weight = st$weight, integer = input$m_integer),
      error = function(e)
        validate(need(FALSE, paste("Could not solve:", conditionMessage(e))))
    )
  })

  output$m_headline <- renderText({
    des <- m_design()
    sprintf(paste0("Weighted detection D: %.4f\nTotal cost: %.4g\n",
                   "%d targets | depth %d | %s counts"),
            des$D, des$C, nrow(des$theta), ncol(des$theta),
            if (des$integer) "integer" else "continuous")
  })

  output$m_design_table <- renderTable({
    des <- m_design()
    d <- ncol(des$theta)
    data.frame(
      "Level" = seq_len(d),
      "<i>c</i><sub>i</sub> (unit cost)" = des$cost,
      "<i>n</i><sub>i</sub> (count)"     = des$n,
      check.names = FALSE
    )
  }, digits = 4, sanitize.text.function = identity,
     sanitize.colnames.function = identity)

  output$m_target_table <- renderTable({
    des <- m_design()
    tg <- des$targets
    tg$binds <- ifelse(seq_len(nrow(tg)) == des$binding, "<-- binds", "")
    ## Label each target with its subscripted index t_i plus its name, e.g.
    ## "t_1: silver carp" (index rendered as Greek-style italic with subscript).
    target_label <- paste0("<i>t</i><sub>", seq_len(nrow(tg)), "</sub>: ",
                           tg$target)
    data.frame(
      "Target (<i>t</i>)"      = target_label,
      "Weight"                 = tg$weight,
      "Detection"              = tg$P1,
      "&lambda; (miss-weight)" = tg$lambda,
      " "                      = tg$binds,
      check.names = FALSE
    )
  }, digits = 4, sanitize.text.function = identity,
     sanitize.colnames.function = identity)

  ## Reproducible script for the assemblage design.
  output$m_download_script <- downloadHandler(
    filename = function() "multidetectr_assemblage.R",
    content = function(file) {
      st <- m_theta()
      cst <- m_cost()
      writeLines(build_assemblage_script(st$theta, st$weight, cst,
                                         input$m_budget, input$m_integer),
                 file)
    }
  )
}

## Assemble a runnable script string reproducing the current design.
build_script <- function(theta, cost, goal, constraint, budget, P_target,
                         integer, enumerate) {
  theta_txt <- paste0("c(", paste(signif(theta, 6), collapse = ", "), ")")
  cost_txt <- paste0("c(", paste(signif(cost, 6), collapse = ", "), ")")
  call_txt <- if (isTRUE(enumerate) && constraint == "budget") {
    sprintf(paste0("enumerate_designs(theta = theta, cost = cost, ",
                   "budget = %s, goal = \"%s\")"),
            signif(budget, 6), goal)
  } else if (constraint == "budget") {
    sprintf(paste0("allocate(theta = theta, cost = cost, goal = \"%s\", ",
                   "budget = %s, integer = %s)"),
            goal, signif(budget, 6), as.character(integer))
  } else {
    sprintf(paste0("allocate(theta = theta, cost = cost, goal = \"%s\", ",
                   "P_target = %s, integer = %s)"),
            goal, signif(P_target, 6), as.character(integer))
  }
  c(
    "# Reproducible script generated by the multiDetectR Shiny app.",
    "# install.packages('multiDetectR')  # or install from the USGS release",
    "library(multiDetectR)",
    "",
    paste0("theta <- ", theta_txt),
    paste0("cost  <- ", cost_txt),
    "",
    paste0("design <- ", call_txt),
    "print(design)",
    "summary(design)",
    "plot(design)",
    "optimality_identity(design)"
  )
}

## Assemble a runnable script for an assemblage design.
build_assemblage_script <- function(theta, weight, cost, budget, integer) {
  d <- ncol(theta)
  nt <- nrow(theta)
  rn <- rownames(theta)
  if (is.null(rn)) rn <- paste0("target_", seq_len(nt))
  ## Build the theta matrix row by row for readability.
  row_txt <- vapply(seq_len(nt), function(t) {
    sprintf("  %-18s = c(%s)", paste0("\"", rn[t], "\""),
            paste(signif(theta[t, ], 6), collapse = ", "))
  }, character(1))
  theta_txt <- paste0("rbind(\n", paste(row_txt, collapse = ",\n"), "\n)")
  c(
    "# Reproducible script generated by the multiDetectR Shiny app (assemblage).",
    "library(multiDetectR)",
    "",
    paste0("theta <- ", theta_txt),
    paste0("weight <- c(", paste(signif(weight, 6), collapse = ", "), ")"),
    paste0("cost   <- c(", paste(signif(cost, 6), collapse = ", "), ")"),
    "",
    sprintf(paste0("design <- allocate_assemblage(theta, cost, budget = %s, ",
                   "weight = weight, integer = %s)"),
            signif(budget, 6), as.character(integer)),
    "print(design)",
    "summary(design)"
  )
}

shinyApp(ui, server)
