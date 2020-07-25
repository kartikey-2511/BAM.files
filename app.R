library(shiny)
library(LSD)
library(pheatmap)
library(DT)
# library(bigmemory)

# C1 = as.big.matrix(readRDS("cpm/b.rds"), type='double')
# C2 = as.big.matrix(readRDS("cpm/c.rds"), type='double')
# colnames(C1) = NULL
# colnames(C2) = NULL
N = readRDS("cpm/prostate.rds")
z = readRDS("cpm/gene.name.rds")
ref = readRDS("cpm/ref.rds")
z = toupper(z)
colnames(z) = c("ENSEMBL ID", "Gene Name")
n = z
rownames(z) = c(z[,2])

ui <- fluidPage(

	titlePanel("Gene expression for Breast Cancer patients"),

	sidebarLayout(

		sidebarPanel(width=3, 
			strong(h4("You can enter either gene name or ENSEMBL id.")),
			textInput("gene1", h5("Enter the first name/Id"),
						value = "TSPAN6"),
			textInput("gene2", h5("Enter the second name/Id"),
						value = "TNMD"),
			br(),
			h5("Please ", strong("search here"), " for valid", strong(" Gene names"), " and their ", strong("ENSEMBL ids.")),
			br(),
			DTOutput("gene.list")
		),
			
		mainPanel(width=9, 
			tabsetPanel(type="tab",
				tabPanel(h4(strong("Expression Plots")),
					column(6,
						plotOutput("ex.plot1", height=700, width=700),
					),
					column(6,
						plotOutput("ex.plot2", height=700, width=700)
					)
				),
				tabPanel(h4(strong("Correlation Heatmap")),
					br(),
					fileInput("genes", h5("Upload the file with gene names")),
					column(6,
						plotOutput("cor.plot1", height=700, width=700)
					),
					column(6,
						plotOutput("cor.plot2", height=700, width=700)
					)
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

	plot_x = reactive({
		str = toupper(input$gene1)
		if (substr(str,1,4)=="ENSG")
			return(str)
		else
			return(z[str,1])
	})

	plot_y = reactive({
		str = toupper(input$gene2)
		if (substr(str,1,4)=="ENSG")
			return(str)
		else
			return(z[str,1])
	})

	output$gene.list = renderDT({
		datatable(n)
	})

	output$ex.plot1 = renderPlot({
		s1 = plot_x()
		s2 = plot_y()
		x1 = read.table("cpm/b.csv",sep=",",nrows=1,skip=ref[s1,1]-1)
		x2 = read.table("cpm/b.csv",sep=",",nrows=1,skip=ref[s2,1]-1)
		v1 = c(unlist(x1))
		v2 = c(unlist(x2))
		heatscatter(v1, v2, cexplot=1, cor=TRUE, method='spearman', log="xy", xlab=toupper(input$gene1), ylab=toupper(input$gene2), main="Cancer tissue")
	})

	output$ex.plot2 = renderPlot({
		s1 = plot_x()
		s2 = plot_y()
		x1 = read.table("cpm/b.csv",sep=",",nrows=1,skip=ref[s1,1]-1)
		x2 = read.table("cpm/b.csv",sep=",",nrows=1,skip=ref[s2,1]-1)
		v1 = c(unlist(x1))
		v2 = c(unlist(x2))
		heatscatter(v1, v2, cexplot=1, cor=TRUE, method='spearman', log="xy", xlab=toupper(input$gene1), ylab=toupper(input$gene2), main="Normal tissue")
	})

	output$cor.plot1 = renderPlot({
		goi = toupper(input_file())
		if (goi[length(goi)] == "")
			goi = goi[1:length(goi)-1]
		l = length(goi)
		if (l > 0) {
			eid = c(z[goi,1])
			m = matrix(nrow=l, ncol=length(read.table("cpm/b.csv",sep=",",nrow=1,skip=1)))
			for (i in c(1:l)) {
				s = eid[i]
				x = read.table("cpm/b.csv",sep=",",nrows=1,skip=ref[s,1]-1)
				# print(length(x))
				# print(dim(m))
				m[i,] = matrix(unlist(x),nrow=1,ncol=length(x))
			}
			m = t(m)
			colnames(m) = c(goi)
			cm = cor(m, method="spearman")
			pheatmap(cm, display_numbers = T, main="Cancer tissue")
		}
	})

	output$cor.plot2 = renderPlot({
		DS = N
		goi = toupper(input_file())
		if (goi[length(goi)] == "")
			goi = goi[1:length(goi)-1]
		l = length(goi)
		if (l > 0) {
			eid = matrix(z[goi,1], nrow=l, ncol=1)
			m = matrix(nrow=l, ncol=length(DS[1,]))
			for (i in c(1:l)) {
				m[i,] = DS[eid[i],]
			}
			m = t(m)
			colnames(m) = c(goi)
			cm = cor(m, method="spearman")
			pheatmap(cm, display_numbers = T, main="Normal tissue")
		}
	})

}

shinyApp(ui = ui, server = server)