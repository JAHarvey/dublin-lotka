# Dublin-Lotka explorer: r ~ ln(R0) / Tc
# BES 550 Advanced Ecology, Population demographics and its drivers
#
# To run: open this file in RStudio and click "Run App",
# or run shiny::runApp("path/to/BES550-DublinLotka_App_v1") in the R console.
# The only package needed is shiny: install.packages("shiny")

library(shiny)

# ---------------------------------------------------------------------------
# Life tables
# ---------------------------------------------------------------------------
# l = survivorship to age a (proportion of the cohort alive)
# f = age-specific fertility (female offspring per female of age a)
# The "illustrative" tables are made-up numbers for teaching. They aren't data.

seabird <- local({
  a <- 0:20
  l <- c(1, 0.5 * 0.92^(a[-1] - 1))
  l[length(l)] <- 0
  data.frame(age = a, l = round(l, 4), f = ifelse(a >= 5 & a < 20, 0.35, 0))
})

presets <- list(
  "Lecture table (slides 6 to 10)" =
    data.frame(age = 0:8,
               l = c(1, 0.25, 0.12, 0.09, 0.06, 0.04, 0.03, 0.02, 0),
               f = c(0, 0, 2, 5, 5, 5, 2, 2, 0)),
  "Fast life history (illustrative)" =
    data.frame(age = 0:5,
               l = c(1, 0.30, 0.15, 0.07, 0.03, 0),
               f = c(0, 2, 2, 2, 2, 0)),
  "Slow life history, seabird-like (illustrative)" = seabird
)

# ---------------------------------------------------------------------------
# Calculations
# ---------------------------------------------------------------------------

# Apply organismal modifiers to a life table.
# surv_mult multiplies annual survival p(a) = l(a+1)/l(a) for ages >= surv_from.
# fert_mult multiplies every f(a).
# shift moves the fertility schedule (negative = breed earlier).
apply_mods <- function(tab, surv_mult = 1, surv_from = 1, fert_mult = 1, shift = 0) {
  n <- nrow(tab)
  l <- tab$l
  f <- tab$f * fert_mult

  if (n > 1) {
    p <- ifelse(l[-n] > 0, l[-1] / l[-n], 0)          # p for ages 0 .. n-2
    idx <- which(tab$age[-n] >= surv_from)
    p[idx] <- pmin(1, p[idx] * surv_mult)
    l <- c(l[1], l[1] * cumprod(p))
  }

  if (shift != 0) {
    f_new <- rep(0, n)
    for (i in seq_len(n)) {
      if (f[i] == 0) next
      j <- min(max(i + shift, 2), n)                   # no breeding at age 0
      f_new[j] <- f_new[j] + f[i]
    }
    f <- f_new
  }
  data.frame(age = tab$age, l = l, f = f)
}

# R0, Tc, approximate r (Dublin-Lotka) and exact r (Euler-Lotka).
dl_stats <- function(tab) {
  lf <- tab$l * tab$f
  R0 <- sum(lf)
  if (!is.finite(R0) || R0 <= 0) {
    return(list(R0 = R0, Tc = NA, r_approx = NA, r_exact = NA, lambda = NA))
  }
  Tc <- sum(tab$age * lf) / R0
  r_approx <- if (Tc > 0) log(R0) / Tc else NA
  el <- function(r) sum(exp(-r * tab$age) * lf) - 1
  r_exact <- tryCatch(uniroot(el, c(-2, 2), extendInt = "yes")$root,
                      error = function(e) NA)
  list(R0 = R0, Tc = Tc, r_approx = r_approx, r_exact = r_exact,
       lambda = exp(r_exact))
}

fmt <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = d))

time_label <- function(r) {
  if (is.na(r) || abs(r) < 1e-6) return("NA")
  t <- log(2) / abs(r)
  paste0(fmt(t, 1), if (r > 0) " yr to double" else " yr to halve")
}

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
teal <- "#169C9A"; red <- "#B5412A"; grey <- "#8A8A8A"

ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-size: 15px; }
    .stat-note { color: #4B5563; font-size: 14px; }
    .scen .btn { margin: 0 6px 6px 0; }
    h4 { margin-top: 18px; }
  "))),
  titlePanel("Dublin-Lotka explorer: r ≈ ln(R₀) / Tc"),

  sidebarLayout(
    sidebarPanel(width = 4,
      selectInput("preset", "Starting life table", choices = names(presets)),

      h4("Organismal scenarios"),
      div(class = "scen",
        actionButton("sc_reset", "Reset"),
        actionButton("sc_immune", "Cost of immunity"),
        actionButton("sc_early", "Earlier breeding"),
        actionButton("sc_disease", "Disease-driven shift")
      ),
      p(class = "stat-note",
        "Scenarios set the sliders below. The numbers are illustrative."),

      h4("Modifiers"),
      sliderInput("surv_mult", "Multiply annual survival by",
                  min = 0.3, max = 1.5, value = 1, step = 0.05),
      numericInput("surv_from", "... starting at age", value = 2, min = 0, step = 1),
      sliderInput("fert_mult", "Multiply fertility by",
                  min = 0, max = 2, value = 1, step = 0.05),
      sliderInput("shift", "Shift breeding schedule (years; negative = earlier)",
                  min = -3, max = 3, value = 0, step = 1),
      checkboxInput("logy", "Log scale for survivorship", value = TRUE)
    ),

    mainPanel(width = 8,
      fluidRow(
        column(6, h4("Baseline vs. modified"), tableOutput("stats")),
        column(6, h4("What it means"), uiOutput("interp"))
      ),
      tabsetPanel(
        tabPanel("Schedules", plotOutput("sched", height = "380px")),
        tabPanel("Euler-Lotka root", plotOutput("elplot", height = "380px"),
                 p(class = "stat-note",
                   "The exact r is where the curve crosses 1. The approximation ln(R₀)/Tc is marked with a dotted line.")),
        tabPanel("Projection", plotOutput("proj", height = "380px")),
        tabPanel("Life table", tableOutput("ltab")),
        tabPanel("Edit life table",
                 p(class = "stat-note",
                   "Edit l(a) and f(a), then click Apply. l(0) should be 1, and l(a) should not increase with age."),
                 actionButton("apply_edits", "Apply edits", class = "btn-primary"),
                 actionButton("add_row", "Add an age class"),
                 actionButton("drop_row", "Remove the oldest age class"),
                 br(), br(),
                 uiOutput("editor")),
        tabPanel("About", uiOutput("about"))
      )
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
server <- function(input, output, session) {

  base_tab <- reactiveVal(presets[[1]])
  observeEvent(input$preset, base_tab(presets[[input$preset]]))

  # Scenario buttons
  set_mods <- function(s, from, fm, sh) {
    updateSliderInput(session, "surv_mult", value = s)
    updateNumericInput(session, "surv_from", value = from)
    updateSliderInput(session, "fert_mult", value = fm)
    updateSliderInput(session, "shift", value = sh)
  }
  observeEvent(input$sc_reset,   set_mods(1,   2, 1,  0))
  observeEvent(input$sc_immune,  set_mods(0.8, 2, 1,  0))
  observeEvent(input$sc_early,   set_mods(1,   2, 1, -1))
  observeEvent(input$sc_disease, set_mods(0.5, 2, 1, -1))

  mod_tab <- reactive({
    sf <- if (is.null(input$surv_from) || is.na(input$surv_from)) 0 else input$surv_from
    apply_mods(base_tab(), input$surv_mult, sf, input$fert_mult, input$shift)
  })
  s0 <- reactive(dl_stats(base_tab()))
  s1 <- reactive(dl_stats(mod_tab()))

  # --- Stats table ---
  output$stats <- renderTable({
    a <- s0(); b <- s1()
    data.frame(
      Quantity = c("R₀ = Σ l(a)f(a)", "Tc = Σ a l(a)f(a) / R₀",
                   "r ≈ ln(R₀)/Tc", "r exact (Euler-Lotka)",
                   "λ = e^r", "Time scale"),
      Baseline = c(fmt(a$R0), fmt(a$Tc, 2), fmt(a$r_approx, 4), fmt(a$r_exact, 4),
                   fmt(a$lambda), time_label(a$r_exact)),
      Modified = c(fmt(b$R0), fmt(b$Tc, 2), fmt(b$r_approx, 4), fmt(b$r_exact, 4),
                   fmt(b$lambda), time_label(b$r_exact)),
      check.names = FALSE
    )
  }, striped = TRUE, spacing = "s")

  # --- Interpretation ---
  output$interp <- renderUI({
    b <- s1()
    if (is.na(b$r_exact)) {
      return(p("With no reproduction, R₀ = 0 and r is undefined: the cohort leaves no offspring."))
    }
    trend <- "just replaces itself, so r = 0"
    if (b$R0 > 1.0001) trend <- "more than replaces itself, so r > 0 and the population grows"
    if (b$R0 < 0.9999) trend <- "fails to replace itself, so r < 0 and the population shrinks"
    err <- 100 * (b$r_approx - b$r_exact) / abs(b$r_exact)
    tagList(
      p(HTML(paste0("R₀ = <b>", fmt(b$R0), "</b>: an average individual ", trend, "."))),
      p(HTML(paste0("Tc = <b>", fmt(b$Tc, 2), " yr</b> is the mean age at which offspring are produced, which sets the speed of the life history."))),
      if (is.finite(err) && abs(b$r_exact) > 1e-6)
        p(class = "stat-note",
          paste0("The approximation differs from the exact r by ", fmt(abs(err), 1),
                 "%. It works best when r is near 0 and reproduction is concentrated around Tc."))
    )
  })

  # --- Plots ---
  output$sched <- renderPlot({
    b0 <- base_tab(); b1 <- mod_tab(); a <- s0(); m <- s1()
    par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1), cex = 1.05)

    y0 <- b0$l; y1 <- b1$l
    if (input$logy) { y0[y0 <= 0] <- NA; y1[y1 <= 0] <- NA }
    ylim <- range(c(y0, y1), na.rm = TRUE)
    plot(b1$age, y1, type = "b", pch = 19, col = teal, lwd = 2.5,
         log = if (input$logy) "y" else "", ylim = ylim,
         xlab = "Age (a)", ylab = "Survivorship l(a)", main = "Survivorship")
    lines(b0$age, y0, type = "b", lty = 2, col = grey, lwd = 2)
    legend("topright", c("Modified", "Baseline"), col = c(teal, grey),
           lty = c(1, 2), lwd = 2, bty = "n")

    lf0 <- b0$l * b0$f; lf1 <- b1$l * b1$f
    plot(b1$age, lf1, type = "h", lwd = 12, lend = 1, col = teal,
         ylim = c(0, max(c(lf0, lf1), 0.01)), xlab = "Age (a)",
         ylab = "l(a) f(a)", main = "Where offspring come from")
    points(b0$age, lf0, pch = 21, bg = "white", col = grey, cex = 1.4, lwd = 2)
    if (!is.na(a$Tc)) abline(v = a$Tc, lty = 2, col = grey, lwd = 2)
    if (!is.na(m$Tc)) abline(v = m$Tc, col = red, lwd = 2.5)
    legend("topright", c("Modified", "Baseline", "Tc modified", "Tc baseline"),
           col = c(teal, grey, red, grey), pch = c(15, 21, NA, NA),
           lty = c(NA, NA, 1, 2), lwd = 2, bty = "n")
  })

  output$elplot <- renderPlot({
    b0 <- base_tab(); b1 <- mod_tab(); a <- s0(); m <- s1()
    roots <- c(a$r_exact, m$r_exact)
    ctr <- if (all(is.na(roots))) 0 else mean(roots, na.rm = TRUE)
    rr <- seq(ctr - 0.6, ctr + 0.6, length.out = 400)
    g <- function(tab, r) sapply(r, function(x) sum(exp(-x * tab$age) * tab$l * tab$f))
    y0 <- g(b0, rr); y1 <- g(b1, rr)
    par(mar = c(4.5, 4.5, 2.5, 1), cex = 1.05)
    plot(rr, y1, type = "l", col = teal, lwd = 3, ylim = c(0, min(4, max(c(y0, y1, 1.5)))),
         xlab = "r", ylab = expression(sum(e^{-r*a} * l(a) * f(a))),
         main = "Solving the Euler-Lotka equation")
    lines(rr, y0, col = grey, lwd = 2, lty = 2)
    abline(h = 1, col = "black")
    if (!is.na(m$r_exact)) abline(v = m$r_exact, col = red, lwd = 2)
    if (!is.na(m$r_approx)) abline(v = m$r_approx, col = red, lty = 3, lwd = 2)
    if (!is.na(a$r_exact)) abline(v = a$r_exact, col = grey, lwd = 1.5)
    legend("topright", c("Modified", "Baseline", "Exact r (modified)", "ln(R₀)/Tc (modified)"),
           col = c(teal, grey, red, red), lty = c(1, 2, 1, 3), lwd = 2, bty = "n")
  })

  output$proj <- renderPlot({
    a <- s0(); m <- s1()
    t <- 0:30; N0 <- 100
    f <- function(r) if (is.na(r)) rep(NA, length(t)) else N0 * exp(r * t)
    ys <- cbind(f(a$r_exact), f(m$r_exact), f(m$r_approx))
    par(mar = c(4.5, 4.5, 2.5, 1), cex = 1.05)
    if (all(is.na(ys))) { plot.new(); text(0.5, 0.5, "No reproduction: nothing to project"); return() }
    matplot(t, ys, type = "l", lty = c(2, 1, 3), lwd = c(2, 3, 2), col = c(grey, teal, red),
            log = "y", xlab = "Years", ylab = "N (log scale)",
            main = "N(t) = 100 e^(rt)")
    legend("topleft", c("Baseline, exact r", "Modified, exact r", "Modified, approximate r"),
           col = c(grey, teal, red), lty = c(2, 1, 3), lwd = 2, bty = "n")
  })

  # --- Life table ---
  output$ltab <- renderTable({
    b <- mod_tab()
    data.frame(`Age (a)` = b$age,
               `l(a)` = round(b$l, 4),
               `f(a)` = round(b$f, 3),
               `l(a) f(a)` = round(b$l * b$f, 4),
               `a l(a) f(a)` = round(b$age * b$l * b$f, 4),
               check.names = FALSE)
  }, digits = 4, striped = TRUE, spacing = "s")

  # --- Editor ---
  output$editor <- renderUI({
    tab <- base_tab()
    rows <- lapply(seq_len(nrow(tab)), function(i) {
      fluidRow(
        column(2, tags$div(style = "padding-top: 8px;", strong(paste("Age", tab$age[i])))),
        column(4, numericInput(paste0("ed_l_", i), NULL, value = tab$l[i], min = 0, max = 1, step = 0.01)),
        column(4, numericInput(paste0("ed_f_", i), NULL, value = tab$f[i], min = 0, step = 0.1))
      )
    })
    tagList(
      fluidRow(column(2, strong("")), column(4, strong("l(a)")), column(4, strong("f(a)"))),
      rows
    )
  })

  read_editor <- function() {
    tab <- base_tab()
    for (i in seq_len(nrow(tab))) {
      lv <- input[[paste0("ed_l_", i)]]; fv <- input[[paste0("ed_f_", i)]]
      if (!is.null(lv) && !is.na(lv)) tab$l[i] <- lv
      if (!is.null(fv) && !is.na(fv)) tab$f[i] <- fv
    }
    tab
  }

  observeEvent(input$apply_edits, {
    tab <- read_editor()
    bad <- any(diff(tab$l) > 1e-9) || tab$l[1] != 1 || any(tab$l < 0) || any(tab$f < 0)
    if (bad) showNotification("Check the table: l(0) should be 1, l(a) should not increase, and values can't be negative.",
                              type = "warning", duration = 6)
    base_tab(tab)
  })
  observeEvent(input$add_row, {
    tab <- read_editor()
    base_tab(rbind(tab, data.frame(age = max(tab$age) + 1, l = 0, f = 0)))
  })
  observeEvent(input$drop_row, {
    tab <- read_editor()
    if (nrow(tab) > 2) base_tab(tab[-nrow(tab), ])
  })

  # --- About ---
  output$about <- renderUI({
    tagList(
      h4("What the app calculates"),
      p("R₀ = Σ l(a) f(a) is expected lifetime reproduction. Tc = Σ a l(a) f(a) / R₀ is the cohort generation time, the mean age at which offspring are produced."),
      p("Dublin & Lotka (1925) approximated the intrinsic rate of increase as r ≈ ln(R₀) / Tc. The exact r solves the Euler-Lotka equation Σ e^(-ra) l(a) f(a) = 1, which the app finds numerically."),
      p("Ages are discrete, f(a) counts female offspring per female, and the life table is treated as a single cohort, as in the lecture."),
      h4("How the modifiers work"),
      p("The survival multiplier changes each annual survival rate p(a) = l(a+1) / l(a) from the chosen age on, capped at 1, then rebuilds l(a). The breeding shift moves the whole fertility schedule to younger or older ages. Scenario buttons only set the sliders, and their values are illustrative."),
      h4("Questions to try"),
      tags$ol(
        tags$li("Click Cost of immunity. Which changes more, R₀ or Tc? Why?"),
        tags$li("Click Earlier breeding. R₀ goes up even though fertility per breeding attempt is unchanged. Where does the gain come from in the l(a) f(a) plot?"),
        tags$li("The earlier-breeding scenario has no cost. Using the energy budget triangle from lecture, what cost is missing, and which slider could you use to add it?"),
        tags$li("Switch between the fast and slow life tables. Which one responds more to a 20% cut in adult survival? Which responds more to a 20% cut in fertility?"),
        tags$li("When does the approximation ln(R₀)/Tc drift furthest from the exact r?")
      ),
      h4("Reference"),
      p("Dublin, L. I. & Lotka, A. J. 1925. On the true rate of natural increase. Journal of the American Statistical Association. ", tags$a(href = "https://www.jstor.org/stable/2965517", "jstor.org/stable/2965517"))
    )
  })
}

shinyApp(ui, server)
