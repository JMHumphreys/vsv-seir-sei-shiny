library("shiny")
library("deSolve")
library("cowplot")
library("ggplot2")
library("tidyverse")
library("gridExtra")
library("ggrepel")
library("shinydashboard")


seir_vect_f = function(time, state, parameters) {
  with( as.list(c(state, parameters)), {
    dS.dt = (-B*Iv*Sh)/N_pop  
    dE.dt = (B*Iv*Sh)/N_pop - (delta * Eh)      
    dI.dt = (delta * Eh) - (gamma * Ih)  
    dR.dt = (gamma * Ih) 
    
    dSv.dt =  Av*Nv - (Bv*Sv*Ih)/N_pop - (Av*Sv) 
    dEv.dt = (Bv*Sv*Ih)/N_pop - (vdelta + Av)*Ev
    dIv.dt = (vdelta * Ev)*pmg - (Av*Iv) 
    
    cum_inf = N_pop-dS.dt
    
    return( list(c(dS.dt, dE.dt, dI.dt, dR.dt, dSv.dt, dEv.dt, dIv.dt, cum_inf)) )
  } )
}


breakfunc = function(x) {
  origin = as.Date("2014-01-01")
  days = origin + x
  origin = as.POSIXlt(origin)
  dayseq = as.POSIXlt(seq(days[1], days[2], by = "day"))
  with(dayseq, yday[mday == 1] + 365*(year[mday == 1] - origin$year ))
}

labelfunc = function(x) {
  origin = as.Date("2014-01-01")
  format(origin + x, format = "%b")
}

#
# Define UI 
#

ui = dashboardPage(
  skin = "blue",
  dashboardHeader(title = "VS SEIR-SEI Model",
                  titleWidth = 300,
                  tags$li(class = "dropdown", 
                          tags$a("Spatial Dynamics", 
                          style = "font-size: 25px; color: white",
                          href="https://jmhumphreys.github.io/"))),
  
  dashboardSidebar(
    br(),
    title = span("Model Parameters", style = "font-size: 30px",
    titleWidth = 300),
    br(),
    br(),
    tags$p("Expand parameter options by selecting on the Horse Population, Vector Biology, Virus Interactions, and Initial Conditions headers below."),
    br(),
    tags$p("Click the Reset to Default Parameters button to return to original values."),
    br(),
    width = 300,
    menuItem("Horse Population",
             tabName = "horses",
             label = "",
             icon = icon("spinner"),
    sliderInput("popsize",
                "Horse Population (thousands):",
                min = 1, max = 10, value = 2
    )
    ),
    br(),
    br(),
    menuItem("Vector Biology",
             tabName = "vectors",
             label = "",
             icon = icon("spinner"),
    sliderInput("vectpop",
               "Vectors per Horse:",
               min = 1, max = 1000, value = 500
             ),
    sliderInput("vrem",
                "Vector Lifespan (days):",
                min = 1, max = 40, value = 28
    ),
    sliderInput("reprov",
                "Vector Sex Ratio (females:1 male):",
                min = 0, max = 5, step = 0.5, value = 2.4
    ),
    sliderInput("sexmat",
                "Vector Age at Reproductive Maturity (days):",
                min = 0, max = 20, value = 10
    ),
    sliderInput("gravidv",
                "Surviving Vector Eggs (per female):",
                min = 0, max = 200, value = 115
    ),
    sliderInput("vlat",
                 "Time Between Vector Feedings (days):",
                 min = 1, max = 10, value = 3
    )
    ),
   br(),
   br(),
   menuItem("Virus Interactions",
            tabName = "virus",
            label = "",
            icon = icon("spinner"),
    sliderInput("vbr",
                "Bite Rate (bites/day):",
                min = 0, max = 500, step = 50, value = 300
    ),
    sliderInput("pmg",
                "Probability Virus Passes Vector Midgut:",
                min = 0, max = 1, step = 0.05, value = 0.5
    ),
    sliderInput("infper",
                "Disease Infectious Period (days):",
                min = 1, max = 20, value = 14
    ),
    sliderInput("latper",
                "Disease Latent Period (days):",
                min = 1, max = 14, value = 3
    )
   ),
   br(),
   br(),
   menuItem("Initial Conditions",
            tabName = "initial",
            label = "",
            icon = icon("spinner"),
    sliderInput("iexp",
                "Initial Horses Exposed:",
                min = 0, max = 50, value = 0
    ),
    sliderInput("pinf",
                "Initial Horse Infections:",
                min = 0, max = 50, value = 0
    ),
    sliderInput("ive",
                "Initial Vectors Exposed:",
                min = 0, max = 50, value = 1
    ),
    sliderInput("ivi",
                "Initial Vector Infections:",
                min = 0, max = 50, value = 1
    )
    ),
   br(),
   br(),
   actionButton("resetBtn", "Reset to Default Parameters"),
   br(),
   br(),
   br(),
   br(),
   HTML(paste0(
     "<table style='margin-left:auto; margin-right:auto;'>",
     "<tr>",
     "<td style='padding: 10px;'><a href='https://github.com/JMHumphreys' target='_blank'><i class='fab fa-github fa-lg'></i></a></td>",
     "<td style='padding: 10px;'><a href='https://scholar.google.com/citations?user=rLUtjhgAAAAJ&hl=en' target='_blank'><i class='fab fa-google fa-lg'></i></a></td>",
     "</tr>",
     "</table>",
     "<br>"),
     HTML(paste0(
       "<script>",
       "var today = new Date();",
       "var yyyy = today.getFullYear();",
       "</script>",
       "<p style = 'text-align: center;'><small>&copy; - <a href='https://jmhumphreys.github.io/' target='_blank'>John Humphreys</a> - <script>document.write(yyyy);</script></small></p>")
     ))
  ),
  dashboardBody(
    tags$head(tags$style(HTML("#sidebarItemExpanded {
            overflow: auto;
            height: calc(100vh - 50px) !important;
        }
        .main-header .logo {
        font-weight: bold;
        font-size: 30px;
        width=300;
      }
        .sidebar-menu{ font-size: 500px; }"))),
    fluidRow(
      column(6, box(plotOutput("dynPlot"),width=12,title=""),
             box(plotOutput("dynPlot2"),width=12,title="")),
      column(6, valueBoxOutput("growBox", width = 6),
             valueBoxOutput("removalBox", width = 6),
             valueBoxOutput("hprevBox", width = 6),
             valueBoxOutput("vprevBox", width = 6),
             valueBoxOutput("hsusBox", width = 6),
             valueBoxOutput("cumulBox", width = 6)),
      column(6, includeHTML("shiny.html"))
    ),
    br()
  )
)

server = function(input, output) {
  
  output$picture <- renderImage({
    return(list(src = "WWW/logo-bg.png",contentType = "image/png",alt = "Spatial Dynamics",width=100, height=100, href="https://jmhumphreys.github.io/"))
  }, deleteFile = FALSE) #where the src is wherever you have the picture
  output$picture2 <- renderImage({
    return(list(src = "WWW/biogeo.png",contentType = "image/png",alt = "Spatial Dynamics",width=30, height=30, href="https://jmhumphreys.github.io/"))
  }, deleteFile = FALSE)
  
  observeEvent(input$resetBtn, {
    updateSliderInput(inputId = "popsize", value = 2)
    updateSliderInput(inputId = "vectpop", value = 500)
    updateSliderInput(inputId = "vrem", value = 28)
    updateSliderInput(inputId = "reprov", value = 2.4)
    updateSliderInput(inputId = "sexmat", value = 10)
    updateSliderInput(inputId = "gravidv", value = 115)
    updateSliderInput(inputId = "vlat", value = 3)
    updateSliderInput(inputId = "vbr", value = 300)
    updateSliderInput(inputId = "pmg", value = 0.5)
    updateSliderInput(inputId = "infper", value = 14)
    updateSliderInput(inputId = "latper", value = 3)
    updateSliderInput(inputId = "iexp", value = 0)
    updateSliderInput(inputId = "pinf", value = 0)
    updateSliderInput(inputId = "ive", value = 1)
    updateSliderInput(inputId = "ivi", value = 1)
  })
  
  dataInput = reactive({
    init       =
      c(
        Sh = input$popsize*1000 - (input$pinf + input$iexp),
        Eh = input$iexp,
        Ih = input$pinf,
        Rh = 0,
        Sv = input$vectpop-input$ivi, #(input$popsize*1000)*input$vectpop - input$ivi,
        Ev = input$ive,
        Iv = input$ivi,
        cum_inf = 0
      )
    
    vect_repro = (input$gravidv*(input$reprov/(input$reprov+1)))/(input$vrem-input$sexmat)
    vect_foi = 1 + 1/(vect_repro*input$vrem)
    v_gamma_value = 1/(input$vrem-input$sexmat)
    
    
    parameters =
      c(B = input$vbr*(1/15^5),
        gamma = (1 / input$infper),
        delta = (1 / input$latper),
        vdelta = (1 / input$vlat),
        pmg = input$pmg,
        N_pop = input$popsize*1000 - (input$pinf + input$iexp),
        Bv = (vect_foi * (v_gamma_value  + input$vrem) / 100),
        Av = (1/input$vrem), 
        Rpv = (input$gravidv*(input$reprov/(input$reprov + 1)))/32,
        Nv = (input$popsize*1000)*input$vectpop
      )
    
    times = seq(0, 365, by = 1)
    
    
    out = ode(
      y = init,
      times = times,
      func = seir_vect_f,
      parms = parameters
    )   
    as.data.frame(out)
    
    
  })
  
  
  output$dynPlot = renderPlot({
    out =
      dataInput() %>%
      gather(key, value, -time) %>%
      mutate(
        id = row_number(),
        key2 = recode(
          key,
          Sh = "Susceptible (S)",
          Eh = "Exposed (E)",
          Ih = "Infected (I)",
          Rh = "Removed (R)",
          Sv = "Suceptibe Vectors (Sv)",
          Ev = "Exposed Vectors (Ev)",
          Iv = "Infected Vectors (Iv)",
          cum_inf = "Cumulative Infections"
        ),
        keyleft = recode(
          key,
          Sh = "Susceptible (S)",
          Eh = "",
          Ih = "Infected (I)",
          Rh = "",
          Sv = "Suceptibe Vectors (Sv)",
          Ev = "",
          Iv = "Infected Vectors (Iv)",
          cum_inf = ""
        ),
        keyright = recode(
          key,
          Sh = "",
          Eh = "Exposed (E)",
          Ih = "",
          Rh = "Removed (R)",
          Sv = "",
          Ev = "Exposed Vectors (Ev)",
          Iv = "",
          cum_inf = "Cumulative Infections"
        )
      )
    
    hout  = out %>% filter(key2 != "Suceptibe Vectors (Sv)" & key2 != "Infected Vectors (Iv)" & key2 != "Exposed Vectors (Ev)" & key2 != "Cumulative Infections")
    
    ggplot(data = hout,
           aes(
             x = time,
             y = (value)/(input$popsize*1000 - (input$pinf + input$iexp)),
             group = key2,
             col = key2,
             label = key2,
             data_id = id
           )) + 
      ylab("Horse Population (%)") + xlab(" ") +
      geom_line(size = 2) +
      geom_text_repel(
        data = subset(hout, time == max(time)),
        aes(label = keyright),
        size = 6,
        segment.size  = 0.2,
        segment.color = "grey50",
        nudge_x = 0,
        hjust = 1,
        direction = "y"
      ) +
      geom_text_repel(
        data = subset(hout, time == min(time)),
        aes(label = keyleft),
        size = 6,
        segment.size  = 0.2,
        segment.color = "grey50",
        nudge_x = 0,
        hjust = 0,
        direction = "y"
      ) +
      theme(legend.position = "none") +
      scale_colour_manual(values = c("green4", "red", "blue", "black")) +
      scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      scale_x_continuous(breaks=breakfunc, labels=labelfunc) + 
      ggtitle("Horse Population Dynamics") +
      theme(
        rect=element_rect(size=0),
        legend.position="none",
        panel.background=element_rect(fill="transparent", colour=NA),
        plot.background=element_rect(fill="transparent", colour=NA),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 24, face = "bold"),
        axis.title.y = element_text(size = 24, face = "bold"),
        axis.text.x = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
      )
    
  })
  
  output$dynPlot2 = renderPlot({
    out =
      dataInput() %>%
      gather(key, value, -time) %>%
      mutate(
        id = row_number(),
        key2 = recode(
          key,
          Sh = "Susceptible (S)",
          Eh = "Exposed (E)",
          Ih = "Infected (I)",
          Rh = "Removed (R)",
          Sv = "Suceptibe Vectors (Sv)",
          Ev = "Exposed Vectors (Ev)",
          Iv = "Infected Vectors (Iv)",
          cum_inf = "Cumulative Infections"
        ),
        keyleft = recode(
          key,
          Sh = "Susceptible (S)",
          Eh = "",
          Ih = "Infected (I)",
          Rh = "",
          Sv = "Suceptibe Vectors (Sv)",
          Ev = "",
          Iv = "Infected Vectors (Iv)",
          cum_inf = ""
        ),
        keyright = recode(
          key,
          Sh = "",
          Eh = "Exposed (E)",
          Ih = "",
          Rh = "Removed (R)",
          Sv = "",
          Ev = "Exposed Vectors (Ev)",
          Iv = "",
          cum_inf = "Cumulative Infections"
        )
      )
    
    vout  = out %>% filter(key2 == "Exposed Vectors (Ev)" | key2 == "Infected Vectors (Iv)")
    
    ggplot(data = vout,
           aes(
             x = time,
             y = value/((input$popsize*1000)*input$vectpop),
             group = key2,
             col = key2,
             label = key2,
             data_id = id
           )) + 
      ylab("Vector Population (%)") + xlab(" ") +
      geom_line(size = 2) +
      geom_text_repel(
        data = subset(vout, time == max(time)),
        aes(label = keyright),
        size = 6,
        segment.size  = 0.2,
        segment.color = "grey50",
        nudge_x = 0,
        hjust = 1,
        direction = "y"
      ) +
      geom_text_repel(
        data = subset(vout, time == min(time)),
        aes(label = keyleft),
        size = 6,
        segment.size  = 0.2,
        segment.color = "grey50",
        nudge_x = 0,
        hjust = 0,
        direction = "y"
      ) +
      theme(legend.position = "none") +
      scale_colour_manual(values = c("lightblue", "darkorange")) +
      scale_x_continuous(breaks=breakfunc, labels=labelfunc) +
      scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      ggtitle("Vector Population Dynamics") +
      theme(
        rect=element_rect(size=0),
        legend.position="none",
        panel.background=element_rect(fill="transparent", colour=NA),
        plot.background=element_rect(fill="transparent", colour=NA),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 24, face = "bold"),
        axis.title.y = element_text(size = 24, face = "bold"),
        axis.text.x = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
      )
    
  })
  
  
  output$cumulBox = renderValueBox({
    valueBox(
      dataInput() %>% select(Sh) %>% 
        summarise(Max = round(max(input$popsize*1000 - (input$pinf + input$iexp)-Sh),0)),
      "Total Horse Infections",
      color = "orange"
    )
  })
  
  output$hprevBox = renderValueBox({
    valueBox(
      dataInput() %>% select(Ih) %>% 
        summarise(Max = round(max(Ih/(input$popsize*1000 - (input$pinf + input$iexp))),4)*100) %>%
        paste0("%"),
      "Peak Prevalence (horse)",
      color = "orange"
    )
  })
  
  output$hsusBox = renderValueBox({
    valueBox(
      dataInput() %>% select(Sh) %>% 
        summarise(Max = round(min(Sh/(input$popsize*1000 - (input$pinf + input$iexp))),4)*100) %>%
        paste0("%"),
      "Still Susceptible (horse)",
      color = "orange"
    )
  })
  
  output$vprevBox = renderValueBox({
    valueBox(
      dataInput() %>% select(Iv) %>% 
        summarise(Max = round(max(Iv/(input$popsize*1000*input$vectpop)),4)*100) %>%
        paste0("%"),
      "Peak Prevalence (vector)",
      color = "orange"
    )
  })
  
  
  output$removalBox = renderValueBox({
    valueBox(
      paste0(round(1/input$vrem, 4)*100, "%"),
      "Vector Population Removel Rate",
      color = "black"
    )
  })
  
  output$growBox = renderValueBox({
    valueBox(
      paste0(round((input$gravidv*(input$reprov/(input$reprov+1)))/(input$vrem-input$sexmat), 2), "%"),
      "Vector Population Growth Rate",
      color = "black"
    )
  })
  
  
  
  
}

# Run the application
shinyApp(ui = ui, server = server)
