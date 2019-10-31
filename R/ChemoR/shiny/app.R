library(tidyverse)
library(ggrepel)
library(chroma)
library(conflicted)
library(shinyWidgets)
library(shiny)
library(ggtips)

conflict_prefer("scale_color_viridis_d", "ggplot2")
conflicted::conflict_prefer("filter", "dplyr")

# setwd("~/GitHub/my_website/content/post/ChemoR/")

type <- readRDS("type.data")

full <- read_csv("family_assignment.csv")

fasta <- read_delim("all_chemor.fa", delim = "\t", col_names = c("Species", "Sequence")) %>%
  separate(Species, c("Species", "Transcript_ID"), sep = "-") %>%
  mutate(Species = str_remove(Species, ">")) %>%
  select(-Species)

full <- left_join(full, fasta) %>%
  filter(!is.na(Sequence))

trimmed_full <- group_by(full, Full) %>%
  sample_n(15, replace = TRUE)

type <- mutate(type, Log.Total = log2(Total))

rect <- as.data.frame(unique(type$Type)) %>%
    mutate(fill = seq(1, 8, 1)) %>%
    rename(Type = 1)

gradient <- matrix(viridis_colors(8), nrow = 50, ncol = length(viridis_colors(8)), byrow = TRUE)

rho_test <- cor.test(as.numeric(type$Type), log2(type$Total), data = type, method = "spearman")
estimate <- rho_test$estimate
p <- rho_test$p.value


# App ---------------------------------------------------------------------


# UI ----------------------------------------------------------------------
ui <- navbarPage("Nematode Chemoreceptors",
                 
                 tabPanel("Interactive Plot",
                          sidebarLayout(
                            sidebarPanel(
                              pickerInput(inputId = "type",
                                 label = h4("Choose Nematode Type:"),
                                 choices = as.character(unique(type$Type)),
                                 multiple = TRUE,
                                 selected = as.character(unique(type$Type)),
                                 options = list(
                                   `actions-box` = TRUE)),
                              
                              radioButtons(inputId = "x",
                                  label = h4("Choose Ordering of X-Axis:"),
                                  choices = c(`Nematode Category` = "Type",
                                              `Genome Contiguity (N50)` = "N50"),
                                  selected = "Type"),
                     
                              radioButtons(inputId = "y",
                                  label = h4("Choose transformation of Y-Axis:"),
                                  choices = c(`Untransformed` = "Total",
                                              `Log2` = "Log.Total"),
                                  selected = "Log.Total"),
                     
                              tableOutput("hover_info"),
                              
                              h4("Updated Spearman's Rho: "),
                              textOutput("rho"),
                              
                              h4("Updated p-value: "),
                              textOutput("p")
                              ),
                            
                            mainPanel(
                              plotOutput("type.plot",
                                height = 600,
                                width = 800,
                                hover = "plot_hover"),
                     
                     h5("Hover near a point to see the data for that species.")))),
                 
                 tabPanel("Data Download",
                          sidebarLayout(
                            sidebarPanel(
                              selectInput(inputId = "species",
                                          label = "Choose a species:",
                                          choices = as.character(unique(type$Full))),
                              
                              downloadButton(outputId = "download",
                                             label = "Download")),
                            
                            mainPanel(
                              
                              h5("Note: only 15 chemoreceptor genes shown below. Click \"Download\" to retrieve the full data set."),
                              
                              tableOutput(outputId = "table")
                              
                            )
                          ))
                 
                 )


# Server ------------------------------------------------------------------
server <- function(input, output) {

    x.scale <- reactive({
        if(input$x == "Type"){scale_x_discrete(limits = c("Free-Living", "Plant-Parasitic", "Facultative", "Skin-Penetrating", "Ingested Larvae", "Ingested Egg", "Vector Transmitted", "Host-Contained"))}
        else{scale_x_continuous(trans = 'log10')}
    })

    y.scale <- reactive({
        if(input$y == "Log.Total"){ylim(1.75, 11.25)}
        else{ylim(0, 1500)}
    })

    intercept <- reactive({
        if(input$x == "Type"){log2(1500)}
        else{5.5}
    })

    pos <- reactive({
        if(input$x == "Type"){position_jitter(width = 0.25, seed = 1)}
        else(position_jitter(width = 0, seed = 1))
    })

    output$type.plot <- renderPlot({
        ggplot(filter(type, Type %in% input$type), aes(x = !!as.symbol(input$x))) +
            geom_abline(slope = round(cor.test(x = as.numeric(pull(filter(type, Type %in% input$type), input$x)),
                                               y = log2(pull(filter(type, Type %in% input$type), Total)),
                                               data = filter(type, Type %in% input$type), method = "spearman")$estimate, digits = 3),
                        intercept = coef(lm(log2(pull(filter(type, Type %in% input$type), Total)) ~ as.numeric(pull(filter(type, Type %in% input$type), input$x)), data = type))[1],
                        linetype = 18,
                        alpha = 0.75,
                        size = 0.5) +
            geom_label_repel(aes(y = !!as.symbol(input$y), label = Full), fill = "white", position = pos(),
                             size = 3, label.size = 0.2, segment.size = 0.2, fontface = "italic") +
            geom_point(aes(y = !!as.symbol(input$y), color = Type), size = 4, position = pos()) +
            y.scale() +
            x.scale() +
            scale_color_viridis_d() +
            labs(x = "", y = "Total Chemoreceptors") +
            theme_minimal(base_size = 16, base_family = "Helvetica") +
            theme(panel.spacing.x = unit(0.25, "line"),
                  panel.spacing.y = unit(0.25, "line"),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.line.x = element_line(color="black", size = 0.25),
              axis.line.y = element_line(color="black", size = 0.25),
              axis.title.x  = element_blank(),
              axis.title.y  = element_text(face = "bold", angle = 90, size = 12),
              axis.text.x = element_text(size = 9.5, face = "bold", angle = 35, hjust = 1),
              axis.text.y = element_text(face = "bold", size = 9),
              strip.text = element_text(face = "bold"),
              axis.ticks.x = element_blank(),
              legend.position = "none") +
        NULL
    })

    x.pos <- reactive({
        if(input$x == "Type" & input$y == "Log.Total"){300}
        if(input$x == "Type" & input$y == "Total"){1200}
        if(input$x == "N50" & input$y == "Log.Total"){550}
        else(700)
    })

    y.pos  <- reactive({
        if(input$x == "Type" & input$y == "Log.Total"){800}
        if(input$x == "Type" & input$y == "Total"){100}
        if(input$x == "N50" & input$y == "Log.Total"){10}
        else(100)
    })

    output$hover_info <- renderUI({
        hover <- input$plot_hover
        point <- nearPoints(select(type, Species = Full, Type, Total, Log.Total, N50),
                            input$plot_hover,
                            xvar = input$x,
                            yvar = input$y,
                            threshold = 20)
        if (nrow(point) == 0) return(NULL)

        # calculate point position INSIDE the image as percent of total dimensions
        # from left (horizontal) and from top (vertical)
        left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
        top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)

        # calculate distance from left and bottom side of the picture in pixels
        left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
        top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)

        # create style property fot tooltip
        # background color is set so tooltip is a bit transparent
        # z-index is set so we are sure are tooltip will be on top
        # offX <- if(hover$x > 2.5){-300} else{300}

        style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                        "left:", x.pos(), "px; top:", y.pos(), "px; width:280px")

        # actual tooltip created as wellPanel
        wellPanel(
            style = style,
            p(HTML(paste0("<b> Species: </b>", point$Species, "<br/>",
                          "<b> Total Chemoreceptors: </b>", point$Total, "<br/>",
                          "<b> Log2 Transformed: </b>", point$Log.Total, "<br/>")))
        )
    })

    output$rho <- renderText({
        round(cor.test(x = as.numeric(pull(filter(type, Type %in% input$type), input$x)),
                       y = log2(pull(filter(type, Type %in% input$type), Total)),
                       data = filter(type, Type %in% input$type), method = "spearman")$estimate,
              digits = 3)
        })

    output$p <- renderText({
        formatC(cor.test(x = as.numeric(pull(filter(type, Type %in% input$type), input$x)),
                         y = log2(pull(filter(type, Type %in% input$type), Total)),
                         data = filter(type, Type %in% input$type), method = "spearman")$p.value,
                format = "e",
                digits = 7)
    })
    
    output$table <- renderTable({
      select(trimmed_full, -Species, Species = Full, Transcript_ID, Superfamily, Family, Sequence) %>%
        filter(Species == input$species)
    })
    
    output$download <- downloadHandler(
      filename = function() {
        str_glue(str_replace(input$species, ". ", ""), ".csv")
      },
      
      content = function(file) {
        out <- select(full, -Species, Species = Full, Transcript_ID, Superfamily, Family, Sequence) %>%
          filter(Species == input$species)
        write_csv(out, file)
      },
      
      contentType = "csv"
    )
}

# Run the application
shinyApp(ui = ui,
         server = server
)


