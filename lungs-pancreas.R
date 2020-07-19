library(shiny)
library(LSD)
library(pheatmap)
library(DT)

L = readRDS("cpm/lungs.rds")
P = readRDS("cpm/pancreas.rds")
z = readRDS("cpm/gene.name.rds")
rownames(z) = c(z[,2])

ui <- fluidPage(

	titlePanel(strong("Gene expression for Cancer patients")),

	sidebarLayout(

		sidebarPanel(
			br(),
			h4(strong("Note:"), " The datasets represented here are for", strong(" Lungs"), " and", strong(" Pancreas cancer"), " patients."),
			br(),
			h5("For", strong("Brain, Head and Neck cancer"), " data, visit the following site:"),
			h5(strong(span("https://kartikey25.shinyapps.io/brain-head-neck-cancer/",style="color:blue"))),
			br(),
			h5("For ", strong("Colorectal and Prostate cancer"), "data, visit the following site:"),
			h5(strong(span("https://kartikey25.shinyapps.io/colorectal-prostate-cancer/",style="color:blue"))),
			br(),
			h5("For ", strong("Breast cancer"), "data, visit the following site:"),
			h5(strong(span("https://kartikey25.shinyapps.io/breast-cancer/",style="color:blue"))),
			br(),
			br(),
			selectInput("dataset", h4(strong("Select the cancer tissue dataset")),
						choices = c("Lungs","Pancreas"),
						selected = "Lungs"),
			br(),
			h5(strong("List of valid Gene names and their ENSEMBL ids")),
			DTOutput("gene.list")
		),
			
		mainPanel(
			tabsetPanel(type="tab",
				tabPanel(h3(strong("Gene name")),
					textInput("gene1", h5("Select the first gene"),
								value = "TSPAN6"),
					textInput("gene2", h5("Select the second gene"),
								value = "TNMD"),
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
			)
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
		if (input$dataset == "Lungs")
			DS = L
		else if (input$dataset == "Pancreas")
			DS = P
		return(DS)
	})

	output$gene.list = renderDT({
		datatable(z)
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
		if (goi[length(goi)] == "")
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