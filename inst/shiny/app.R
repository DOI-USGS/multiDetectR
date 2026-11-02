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
          "&theta; and cost, and either a fixed budget (maximize detection) or ",
          "a required detection probability (minimize budget) &mdash; the two ",
          "design problems. You get the per-level design, a budget&ndash;",
          "detection frontier with the integer step function overlaid, an ",
          "optimality audit, and the option to compare the optimum against a ",
          "design you enter yourself. One click exports a reproducible R ",
          "script."))),
        tags$li(HTML(paste0(
          "<b>Multispecies assemblage.</b> Several targets share one survey but ",
          "detect differently. Enter each target's &theta; across the levels ",
          "(edit the table or upload a CSV) and a relative weight. The tool ",
          "returns the single design that best serves the whole set, which is ",
          "governed by the hardest-to-detect targets, not the average one, and ",
          "a cost&ndash;detection frontier for the assemblage."))),
        tags$li(HTML(paste0(
          "<b>Branches (methods/places).</b> When the options are not nested ",
          "&mdash; competing methods on the same material, or separate places ",
          "&mdash; enter how likely one unit of each branch is to detect each ",
          "target and what one unit costs, and the tool returns the budget ",
          "split and whether to concentrate on one branch or spread across ",
          "them.")))
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
      radioButtons("constraint", "Which optimization problem?",
                   choices = c(
                     "Fixed budget \u2192 maximize detection" = "budget",
                     "Fixed detection \u2192 minimize budget" = "P_target"),
                   selected = "budget"),
      helpText(HTML(paste0(
        "These are the two design problems: either you fix the budget and the ",
        "rule maximizes detection, or you fix a required detection probability ",
        "and it finds the cheapest design that reaches it."))),
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
      helpText(HTML(paste0(
        "Exhaustive search is a separate choice from the two problems above: it ",
        "finds the exact integer optimum by trying every design, which matters ",
        "only at small budgets where rounding the solution can miss the best ",
        "design."))),
      hr(),
      ## Compare the optimum against a user-supplied design (Kevin's request).
      checkboxInput("compare_alt",
                    "Compare against my own design", value = FALSE),
      conditionalPanel(
        "input.compare_alt == true",
        helpText(HTML(paste0(
          "Enter your own per-level counts <i>n</i><sub>i</sub> (outermost ",
          "first). The tool shows your design's detection and cost next to the ",
          "optimum, so you can see what a chosen design gives up."))),
        uiOutput("alt_inputs")
      ),
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
      helpText(HTML(paste0(
        "The smooth curve is the continuous frontier, ",
        "<i>P</i>(<i>B</i>) = 1 &minus; exp(&minus;<i>B</i> / &kappa;<sub>1</sub>). ",
        "With integer replicate counts the achievable detection is a STEP ",
        "function (grey): it rises only when the budget crosses the cost of the ",
        "next whole outermost unit, so between steps extra budget buys nothing. ",
        "Hover to read the budget and detection at the cursor."))),
      plotOutput("frontier", height = "340px",
                 hover = hoverOpts(id = "frontier_hover", delay = 50,
                                   delayType = "throttle")),
      verbatimTextOutput("frontier_readout"),
      conditionalPanel(
        "input.compare_alt == true",
        h4("Your design vs the optimum"),
        helpText(HTML(paste0(
          "Your entered design compared with the rule's optimum at the same ",
          "settings: the detection each reaches and what it costs."))),
        tableOutput("compare_table")
      ),
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
        tableOutput("m_target_table"),
        h4("Assemblage cost-detection frontier (Pareto analogue)"),
        helpText(HTML(paste0(
          "The weighted detection <i>D</i> the assemblage reaches as the budget ",
          "grows, re-optimizing the shared design at each budget. This is the ",
          "multispecies analogue of the single-species frontier and of the ",
          "cost-vs-richness Pareto front: <i>D</i> plays the role of richness. ",
          "The marked point is the current budget. Detection climbs steeply at ",
          "first, then flattens as the elusive targets become the binding ",
          "constraint."))),
        plotOutput("m_frontier", height = "320px")
      )
    )
    ), # end Tab 2

    ## ------------------------------------------------------------------ ##
    ## Tab 3: branches (complementary methods or places)                  ##
    ## ------------------------------------------------------------------ ##
    tabPanel(
    "Branches (methods/places)",
    br(),
    sidebarLayout(
      sidebarPanel(
        width = 4,
        helpText(HTML(paste0(
          "Branches are options that are NOT nested within one another: ",
          "separate methods applied to the same material, or separate places ",
          "to sample. The tool returns how to split a budget across them to ",
          "maximize the chance of detecting the targets, and whether to ",
          "concentrate on one branch or spread across them (manuscript ",
          "Supplement S1.5)."))),
        helpText(HTML(paste0(
          "You describe each branch in everyday terms: how likely <b>one unit</b> ",
          "of effort (one sample, one assay) is to detect each target, and what ",
          "one unit <b>costs</b>. The tool converts these to detection rates for ",
          "you."))),
        sliderInput("b_targets", "Number of targets",
                    min = 1L, max = max_targets, value = 1L, step = 1L),
        sliderInput("b_branches", "Number of branches (methods/places)",
                    min = 2L, max = 8L, value = 2L, step = 1L),
        numericInput("b_budget", "Total budget ($)",
                     value = 1000, min = 0, step = 50),
        hr(),
        downloadButton("b_download_script", "Download R script")
      ),
      mainPanel(
        width = 8,
        h4("Per-unit detection probability (editable)"),
        helpText(HTML(paste0(
          "Rows are targets, columns are branches. Each cell is the probability ",
          "that <b>one unit</b> of that branch detects that target (0 to 1). ",
          "A 0 means the branch cannot detect that target (e.g., a place where ",
          "the target is absent)."))),
        DT::DTOutput("b_prob_table"),
        tags$br(),
        h4("Cost of one unit of each branch ($)"),
        helpText("What a single sample or assay from each branch costs."),
        DT::DTOutput("b_cost_table"),
        tags$br(),
        verbatimTextOutput("b_headline"),
        h4("Optimal budget split"),
        tableOutput("b_share_table"),
        plotOutput("b_share_plot", height = "260px")
      )
    )
    ) # end Tab 3

  ) # end tabsetPanel
  ), # end main-content div
  div(id = "floating-footer",
    htmltools::HTML(footer_content)
  )
)

server <- function(input, output, session) {

  ## Dynamically generated per-level theta and cost inputs. Each level gets a
  ## heading that names where it sits in the hierarchy (level 1 = outermost /
  ## reporting, level d = innermost), so a user always knows which end is which.
  output$level_inputs <- renderUI({
    d <- input$d
    defaults_theta <- rep(0.5, d)
    defaults_cost <- 10^(seq(d, 1) - 1) * 2   # outer levels cost more
    level_label <- function(i, d) {
      if (d == 1L) return("Level 1 (the only level)")
      if (i == 1L) return(sprintf("Level %d \u2014 outermost (reporting)", i))
      if (i == d)  return(sprintf("Level %d \u2014 innermost", i))
      sprintf("Level %d", i)
    }
    lapply(seq_len(d), function(i) {
      tagList(
        tags$p(tags$b(level_label(i, d)),
               style = "margin-bottom:2px; margin-top:8px;"),
        fluidRow(
          column(6, numericInput(
            paste0("theta_", i),
            HTML(sprintf("&theta;<sub>%d</sub> (detectability)", i)),
            value = defaults_theta[i], min = 0.001, max = 0.999, step = 0.01)),
          column(6, numericInput(
            paste0("cost_", i),
            HTML(sprintf("<i>c</i><sub>%d</sub> (unit cost)", i)),
            value = defaults_cost[i], min = 0.0001, step = 1))
        )
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
    des <- design()
    kappa1 <- des$kappa1
    ## Smooth continuous frontier, P(B) = 1 - exp(-B / kappa1).
    max_budget <- kappa1 * nats(0.999)
    Bc <- seq(0, max_budget, length.out = 400)
    Pc <- -expm1(-Bc / kappa1)
    graphics::plot(Bc, Pc, type = "l", lwd = 2, col = "black",
                   xlab = "Budget", ylab = "Detection probability",
                   ylim = c(0, 1))
    ## Integer step function: with the inner design fixed, each whole outermost
    ## unit costs C1 and the achievable detection jumps only at multiples of C1.
    ## P at n_1 outermost units is 1 - exp(-n_1 * h_1) = 1 - exp(-n_1 * C1 / kappa1).
    C1 <- des$C / max(des$n[1], 1)           # cost of one outermost unit
    if (is.finite(C1) && C1 > 0) {
      n1_max <- max(1, floor(max_budget / C1))
      n1_max <- min(n1_max, 2000)            # keep the step set tractable
      n1_seq <- 0:n1_max
      B_step <- n1_seq * C1
      P_step <- -expm1(-B_step / kappa1)
      ## Draw as a staircase: detection is flat from one unit's cost up to the
      ## next, then jumps.
      graphics::segments(B_step[-length(B_step)], P_step[-length(P_step)],
                         B_step[-1], P_step[-length(P_step)],
                         col = "grey55", lwd = 2)                 # flat treads
      graphics::segments(B_step[-1], P_step[-length(P_step)],
                         B_step[-1], P_step[-1],
                         col = "grey55", lwd = 1, lty = 2)        # risers
      graphics::points(B_step, P_step, pch = 19, cex = 0.5, col = "grey40")
    }
    ## Mark this design's point.
    if (des$C <= max_budget && des$P1 > 0) {
      graphics::points(des$C, des$P1, pch = 19, col = "firebrick", cex = 1.3)
      graphics::segments(des$C, 0, des$C, des$P1, lty = 3, col = "firebrick")
    }
    graphics::legend("bottomright",
                     legend = c("continuous frontier", "integer steps",
                                "this design"),
                     col = c("black", "grey55", "firebrick"),
                     lwd = c(2, 2, NA), pch = c(NA, NA, 19), bty = "n")
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
      return(data.frame(note = "Single-level design: no inner levels to audit."))
    }
    ## Relabel A and B with what each term means (Kevin's request).
    data.frame(
      "Level"                              = id$level,
      "<i>A</i> (detection share inside)"  = id$A,
      "<i>B</i> (cost share inside)"       = id$B,
      "<i>A</i> &minus; <i>B</i>"          = id$diff,
      check.names = FALSE
    )
  }, digits = 4, sanitize.text.function = identity,
     sanitize.colnames.function = identity)

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

  ## --- Compare against a user-supplied design (optimum vs alternative) ----- ##
  ## One count input per level, seeded from the current optimum so the user
  ## starts from the recommended design and edits it.
  output$alt_inputs <- renderUI({
    d <- input$d
    des <- tryCatch(design(), error = function(e) NULL)
    seed_n <- if (!is.null(des) && length(des$n) == d) round(des$n) else rep(1, d)
    lapply(seq_len(d), function(i) {
      end <- if (i == 1L) " (outermost)" else if (i == d) " (innermost)" else ""
      numericInput(paste0("alt_n_", i),
                   HTML(sprintf("<i>n</i><sub>%d</sub>%s", i, end)),
                   value = seed_n[i], min = 0, step = 1)
    })
  })

  output$compare_table <- renderTable({
    req(isTRUE(input$compare_alt))
    des <- design()
    d <- input$d
    iv <- inputs()
    alt_n <- vapply(seq_len(d),
                    function(i) input[[paste0("alt_n_", i)]], numeric(1))
    req(!any(is.na(alt_n)))
    validate(need(all(alt_n >= 0), "Counts must be zero or positive."))
    ## Detection and cost of the user's design, via the same package internals.
    alt_P <- if (all(alt_n >= 1)) detection_recursion(iv$theta, alt_n)$P1 else 0
    alt_C <- total_cost(iv$cost, alt_n)
    data.frame(
      "Design"                         = c("Optimum (rule)", "Your design"),
      "Counts <i>n</i> (outer&rarr;inner)" =
        c(paste(round(des$n), collapse = " / "),
          paste(round(alt_n), collapse = " / ")),
      "Detection"                      = c(des$P1, alt_P),
      "Total cost"                     = c(des$C, alt_C),
      check.names = FALSE
    )
  }, digits = 4, sanitize.text.function = identity,
     sanitize.colnames.function = identity)

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
    ## Default theta: make species DIFFER so the "hardest target binds the
    ## design" behaviour is visible out of the box, rather than an all-0.5 table.
    ## Deterministic (so the app is reproducible): each species gets a base
    ## detectability stepped across a range, with a mild rise toward the inner
    ## levels; everything clamped to (0, 1). Users still edit freely.
    base <- if (nt == 1L) 0.5 else seq(0.3, 0.7, length.out = nt)
    level_adj <- if (d == 1L) 0 else 0.05 * (seq_len(d) - (d + 1) / 2)
    th <- outer(base, level_adj, `+`)
    th <- pmin(pmax(th, 0.05), 0.95)
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
    ## Greek + subscript column headers (HTML; escape = FALSE below renders them),
    ## with the outermost/innermost ends marked so the direction is unambiguous.
    dlev <- ncol(st$theta)
    theta_hdr <- vapply(seq_len(dlev), function(i) {
      end <- if (i == 1L) " (outermost)" else if (i == dlev) " (innermost)" else ""
      paste0("&theta;<sub>", i, "</sub>", end)
    }, character(1))
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
    ## Greek/italic + subscript headers for the shared per-level costs, with the
    ## outermost/innermost ends marked.
    dlev <- length(cst)
    colnames(df) <- vapply(seq_along(cst), function(i) {
      end <- if (i == 1L) " (outermost)" else if (i == dlev) " (innermost)" else ""
      paste0("<i>c</i><sub>", i, "</sub>", end)
    }, character(1))
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

  ## Assemblage cost-detection frontier (Pareto analogue): re-optimize the shared
  ## design across a grid of budgets and plot weighted detection D vs cost.
  output$m_frontier <- renderPlot({
    st <- m_theta()
    cst <- m_cost()
    req(st, cst, ncol(st$theta) == length(cst),
        all(st$theta > 0 & st$theta < 1), all(cst > 0),
        any(st$weight > 0), !is.na(input$m_budget), input$m_budget > 0)
    ## A budget grid spanning from a small design to ~2x the current budget, so
    ## the current point sits inside the plotted range.
    bmax <- max(input$m_budget * 2, input$m_budget + sum(cst))
    budgets <- seq(sum(cst), bmax, length.out = 40)
    pts <- lapply(budgets, function(b) {
      des <- tryCatch(
        allocate_assemblage(st$theta, cst, budget = b, weight = st$weight,
                            integer = input$m_integer),
        error = function(e) NULL)
      if (is.null(des) || des$n[1] < 1) return(NULL)
      c(cost = des$C, D = des$D)
    })
    pts <- do.call(rbind, pts[!vapply(pts, is.null, logical(1))])
    validate(need(!is.null(pts) && nrow(pts) >= 2,
                  "Not enough feasible budgets to draw a frontier."))
    graphics::plot(pts[, "cost"], pts[, "D"], type = "b", pch = 19, cex = 0.5,
                   col = "steelblue", lwd = 2,
                   xlab = "Total cost", ylab = "Weighted detection D",
                   ylim = c(0, 1))
    cur <- tryCatch(m_design(), error = function(e) NULL)
    if (!is.null(cur) && cur$n[1] >= 1) {
      graphics::points(cur$C, cur$D, pch = 19, col = "firebrick", cex = 1.3)
      graphics::segments(cur$C, 0, cur$C, cur$D, lty = 3, col = "firebrick")
      graphics::legend("bottomright",
                       legend = c("re-optimized frontier", "current budget"),
                       col = c("steelblue", "firebrick"),
                       lwd = c(2, NA), pch = c(19, 19), bty = "n")
    }
  })

  ## =================================================================== ##
  ## Tab 3: branches (complementary methods or places)                   ##
  ## =================================================================== ##

  ## The user supplies intuitive quantities -- a per-UNIT detection probability
  ## for each target x branch, and a per-branch unit cost -- and we convert to
  ## the hazard-per-dollar rate the rule needs. One unit of branch k costs
  ## cost_k and delivers hazard -log(1 - p_tk) nats against target t, so
  ## rate_tk = -log(1 - p_tk) / cost_k. Nats never appear in the UI.
  b_prob <- reactiveVal(NULL)   # targets x branches per-unit detection probs
  b_cost <- reactiveVal(NULL)   # per-branch unit cost

  ## (Re)build defaults when the target/branch counts change, keeping overlap.
  observeEvent(list(input$b_targets, input$b_branches), {
    nt <- input$b_targets
    nk <- input$b_branches
    oldp <- b_prob()
    ## Default probabilities: make DIFFERENT targets favour DIFFERENT branches,
    ## so the "spread the budget" behaviour shows out of the box rather than the
    ## trivial case where one branch is best for everyone. Deterministic.
    if (nt == 1L) {
      p <- matrix(round(seq(0.6, 0.2, length.out = nk), 2), nrow = 1)
    } else {
      hi <- 0.6; lo <- 0.2
      p <- matrix(lo, nrow = nt, ncol = nk)
      for (t in seq_len(nt)) {
        peak <- ((t - 1) %% nk) + 1          # each target's preferred branch
        p[t, ] <- round(lo + (hi - lo) * (0.4 ^ abs(seq_len(nk) - peak)), 2)
      }
    }
    if (!is.null(oldp)) {
      ro <- min(nrow(oldp), nt); co <- min(ncol(oldp), nk)
      p[seq_len(ro), seq_len(co)] <- oldp[seq_len(ro), seq_len(co), drop = FALSE]
    }
    rownames(p) <- paste0("target_", seq_len(nt))
    colnames(p) <- paste0("branch_", seq_len(nk))
    b_prob(p)

    oldc <- b_cost()
    cst <- rep(50, nk)                       # a plausible default unit cost
    if (!is.null(oldc)) {
      co <- min(length(oldc), nk)
      cst[seq_len(co)] <- oldc[seq_len(co)]
    }
    b_cost(cst)
  }, ignoreInit = FALSE)

  ## Editable per-unit detection-probability table.
  output$b_prob_table <- DT::renderDT({
    p <- b_prob(); req(p)
    df <- data.frame(target = rownames(p), p, check.names = FALSE)
    colnames(df) <- c("target", paste0("branch ", seq_len(ncol(p))))
    DT::datatable(df, rownames = FALSE, editable = "cell", escape = FALSE,
                  options = list(dom = "t", paging = FALSE, ordering = FALSE))
  }, server = TRUE)

  observeEvent(input$b_prob_table_cell_edit, {
    info <- input$b_prob_table_cell_edit
    p <- b_prob()
    j <- info$col            # 0 = target name, 1.. = branch probabilities
    val <- suppressWarnings(as.numeric(info$value))
    if (j == 0) {
      rownames(p)[info$row] <- as.character(info$value)
    } else if (!is.na(val) && val >= 0 && val < 1) {
      p[info$row, j] <- val   # keep in [0, 1); 1 would be infinite hazard
    }
    b_prob(p)
  })

  ## Editable per-branch unit-cost row.
  output$b_cost_table <- DT::renderDT({
    cst <- b_cost(); req(cst)
    df <- as.data.frame(as.list(cst))
    colnames(df) <- paste0("branch ", seq_along(cst), " ($/unit)")
    DT::datatable(df, rownames = FALSE, editable = "cell", escape = FALSE,
                  options = list(dom = "t", paging = FALSE, ordering = FALSE))
  }, server = TRUE)

  observeEvent(input$b_cost_table_cell_edit, {
    info <- input$b_cost_table_cell_edit
    cst <- b_cost()
    val <- suppressWarnings(as.numeric(info$value))
    if (!is.na(val) && val > 0) cst[info$col + 1] <- val
    b_cost(cst)
  })

  ## Convert probability + cost to the hazard-per-dollar rate the rule needs.
  b_rate_mat <- reactive({
    p <- b_prob(); cst <- b_cost()
    req(p, cst, ncol(p) == length(cst))
    validate(
      need(all(p >= 0 & p < 1),
           "Each detection probability must be in [0, 1)."),
      need(all(cst > 0), "Each unit cost must be positive.")
    )
    ## rate_tk = -log(1 - p_tk) / cost_k
    r <- sweep(-log1p(-p), 2, cst, "/")
    rownames(r) <- rownames(p)
    colnames(r) <- colnames(p)
    r
  })

  b_result <- reactive({
    r <- b_rate_mat()
    req(!is.na(input$b_budget), input$b_budget > 0)
    validate(need(any(r > 0),
                  "At least one detection probability must be above 0."))
    tryCatch(
      allocate_branches(rate = r, budget = input$b_budget),
      error = function(e)
        validate(need(FALSE, paste("Could not solve:", conditionMessage(e)))))
  })

  output$b_headline <- renderText({
    res <- b_result()
    regime_txt <- c(dominance = "one branch is best for every target (dominance)",
                    `small-budget` = "small-budget regime (weighted average over targets)",
                    `large-budget` = "large-budget regime (the branch's worst target decides)")
    sprintf(paste0("Weighted detection: %.4f\n%s\nTheorem 2 regime: %s"),
            res$detection,
            if (res$concentrate)
              "Concentrate: put the whole budget on the single best branch."
            else
              "Spread: split the budget across branches (shares below).",
            regime_txt[[res$regime]])
  })

  output$b_share_table <- renderTable({
    res <- b_result()
    cst <- b_cost()
    amount <- res$share * input$b_budget
    data.frame(
      "Branch"               = names_or_index(b_prob()),
      "Budget share"         = res$share,
      "Budget amount ($)"    = amount,
      "Units (amount/cost)"  = floor(amount / cst),
      check.names = FALSE
    )
  }, digits = 4)

  output$b_share_plot <- renderPlot({
    res <- b_result()
    graphics::barplot(res$share, names.arg = colnames(b_prob()),
                      ylab = "Share of budget", ylim = c(0, 1),
                      col = "steelblue", border = NA)
  })

  output$b_download_script <- downloadHandler(
    filename = function() "multidetectr_branches.R",
    content = function(file) {
      writeLines(build_branches_script(b_prob(), b_cost(), input$b_budget),
                 file)
    }
  )
}

## Branch column labels (branch names, or indices if unnamed).
names_or_index <- function(p) {
  cn <- colnames(p)
  if (is.null(cn)) paste0("branch_", seq_len(ncol(p))) else cn
}

## Assemble a runnable script for a branches (methods/places) allocation.
## Emits the intuitive inputs (per-unit detection probability and unit cost) and
## the conversion to the hazard-per-dollar rate, so the script matches the app.
build_branches_script <- function(prob, cost, budget) {
  nt <- nrow(prob)
  row_txt <- vapply(seq_len(nt), function(t) {
    sprintf("  c(%s)", paste(signif(prob[t, ], 6), collapse = ", "))
  }, character(1))
  prob_txt <- paste0("rbind(\n", paste(row_txt, collapse = ",\n"), "\n)")
  c(
    "# Reproducible script generated by the multiDetectR Shiny app (branches).",
    "library(multiDetectR)",
    "",
    "# Per-unit detection probability (rows = targets, columns = branches)",
    paste0("prob <- ", prob_txt),
    "# Cost of one unit of each branch",
    paste0("cost <- c(", paste(signif(cost, 6), collapse = ", "), ")"),
    "",
    "# Convert to the hazard-per-dollar rate the rule uses:",
    "#   one unit of branch k costs cost[k] and gives hazard -log(1 - prob)",
    "rate <- sweep(-log1p(-prob), 2, cost, \"/\")",
    "",
    sprintf("branches <- allocate_branches(rate = rate, budget = %s)",
            signif(budget, 6)),
    "str(branches)"
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
