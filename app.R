library(shiny)
library(LSD)
library(pheatmap)

B = readRDS("cpm/brain.rds")
L = readRDS("cpm/lungs.rds")
z = readRDS("cpm/gene.name.rds")
rownames(z) = c(z[,2])

ui <- fluidPage(

	titlePanel("Gene expression for Cancer patients"),

	sidebarLayout(

		sidebarPanel(
			br(),
			h4(strong("Note:"), " The datasets represented here are for", strong(" Brain"), " and", strong(" Lung cancer"), " patients."),
			br(),
			h5("For ", strong("Colorectal, Head, Neck, and Prostate cancer"), "data, visit the following site:"),
			h5(strong(span("https://kartikey25.shinyapps.io/rse-chnp/",style="color:blue"))),
			br(),
			h5("For", strong("Breast and Pancreas cancer"), " data, visit the following site:"),
			h5(strong(span("https://kartikey25.shinyapps.io/rse-bp/",style="color:blue"))),
			br(),
			br(),
			selectInput("dataset", h3(strong("Select the cancer tissue dataset")),
						choices = c("Brain","Lungs"),
						selected = "Brain"),
			br()
		),

		mainPanel(
			tabsetPanel(type="tab",
				tabPanel(h3(strong("Gene name")),
					selectInput("gene1", h5("Select the first gene"),
								choices = c(z[,2]),	selected = "TSPAN6"),
					selectInput("gene2", h5("Select the second gene"),
								choices = c(z[,2]),	selected = "TNMD"),
					br(),
					h4(strong("Gene expression for the 2 genes selected")),
					plotOutput("name.plot",height=700,width=700)
				),
				tabPanel(h3(strong("ENSEMBL ID")),
					textInput("id1", h5("Select the first ENSEMBL id"), 
                    	 		value = "ENSG00000000003.14"),
					textInput("id2", h5("Select the second ENSEMBL id"), 
                    	 		value = "ENSG00000000005.5"),
					br(),
					h4(strong("Gene expression for the 2 ENSEMBL id selected")),
					plotOutput("id.plot",height=700,width=700)
				),
				tabPanel(h3(strong("Correlation heatmap")),
					br(),
					fileInput("genes", h5("Upload the file with gene names")),
					br(),
					plotOutput("genes.plot",height=800,width=800)
				)
			),
		)
	)
)

server <- function(input, output) {

	input_file <- reactive({
    if (is.null(input$genes)) {
      return("")
    }
    else {
    	readLines(input$genes$datapath)
    }
	})

	tissue <- reactive({
		if (input$dataset == "Brain")
			DS = B
		else if (input$dataset == "Lungs")
			DS = L
		return(DS)
	})

	output$name.plot = renderPlot({
		DS = tissue()
		s1 = input$gene1
		s2 = input$gene2
		r1 = z[s1,1]
		r2 = z[s2,1]
		v1 = DS[r1,]
		v2 = DS[r2,]
		heatscatter(v1, v2, cexplot=1, cor=TRUE, method='spearman', log="xy", xlab=s1, ylab=s2)

	})

	output$id.plot = renderPlot({
		DS = tissue()
		v1 = DS[input$id1,]
		v2 = DS[input$id2,]
		heatscatter(v1, v2, cexplot=1, cor=TRUE, method='spearman', log="xy", xlab=input$id1, ylab=input$id2)
	})

	output$genes.plot = renderPlot({
		DS = tissue()
		goi = input_file()
		if (goi[length(goi)]=="")
			goi = goi[1:length(goi)-1]
		l = length(goi)
		if (l > 0) {
			eid = matrix(z[goi,1], nrow=l, ncol=1)
			m = matrix(DS[eid,], nrow=l, ncol=dim(DS)[2])
			m = t(m)
			colnames(m) = c(goi)
			cm = cor(m, method="spearman")
			pheatmap(cm, display_numbers = T)
		}
	})

}

shinyApp(ui = ui, server = server)