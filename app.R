library(tidyverse)
library(shiny)
library(shinydashboard)
library(shinyWidgets)


source("jaxmat.R")   #for displaying mathematics
stylesheet <- tags$head(tags$style(HTML('
    .main-header .logo {
      font-family: "Georgia", Times, "Times New Roman", serif;
      font-weight: bold;
      font-size: 24px;
    }
    table {
      
    }
  ')
))
#The user interface
header <- dashboardHeader(title = "Finite Field Matricess",
                          titleWidth = 500)
sidebar <- dashboardSidebar(disable = TRUE)
body <- dashboardBody(
  fluidRow(stylesheet,
           column(width=5,
                  
                  #first step: generate and display matrix
                  actionBttn("btn_gen","Generate 3 × 3 Matrix"),
                  uiOutput("matA", style= "color:orange"),
                  
                  #next: generate characteristic polynomial
                  h4("A characteristic polynomial of square matrix A is defined as det (A-λI) where I is the identity matrix of same dimension as A."),
                  actionBttn("btn_charpoly","Generate Characteristic Polynomial"),
                  uiOutput("trdet", style= "color:orange"),
                  
                  #generate eigenvalues
                  actionBttn("btn_eigen","Find the eigenvalues"),
                  uiOutput("eigen_o", style= "color:orange"),
                  
                  #verify the C-H theorem
                  actionBttn("btn_CH","Verify Caley-Hamilton"),
                  
           )
  )
)
ui <- dashboardPage(header, sidebar, body, skin = "purple") #other colors available

#Functions that implement the mathematics
#This file must go into the same directory as app.R
source("eigencalc.R")
#Any images or stylesheet must go into a www subfolder

#Additional functions are OK here, but no variables


server <- function(session, input, output) {
  #Variables that are shared among server functions (use <<-)
  N = 3
  #A <- matrix(nrow=N, ncol=N)
  trA <- 0
  detA <- 0
  disc <- 0
  lam1 <- 0
  lam2 <- 0
  v1 <- numeric(N)
  v2 <- numeric(N)
  D <- matrix(nrow=N, ncol=N)
  P <- matrix(nrow=N, ncol=N)
  G <- matrix(nrow=N, ncol=N)
  #Initialization
  
  
  #Make reactive variable for the matrix
  A_r = eventReactive(input$btn_gen, {
    A <<- eig.make3x3()
  })
  
  #outputs the matrix
  output$matA <- renderUI(jax.matrix(A_r(),name = "A"))
  
  #make reactive variables for the elements of the characteristic polynomial
  trA_r = eventReactive(input$btn_charpoly, sum(diag(A_r() )) )
  detA_r = eventReactive(input$btn_charpoly, det(A) )
  discA_r = eventReactive(input$btn_charpoly, trA_r()^2-4*detA_r() )
  
  #outputs the characteristic polynomial
  output$trdet <- renderUI(h4("Trace = ",trA_r(),"  Determinant = ",detA_r(),
                              "  Discriminant = ",round(discA_r(), digits = 3)))
  
  #make reactive variables for the eigenvalues
  eigen_vals_r = eventReactive(input$btn_eigen, {
    eig_A = eigen(A_r())
    eig_A_vals = eig_A$values %>% round(2)
    paste(eig_A_vals, collapse = ", ")
  })
  
  #output the eigen values
  output$eigen_o <- renderUI(h4( eigen_vals_r() ))
  
  #check the C-H theorem
  observeEvent(input$btn_CH, {
    
  })
  
  #Functions that respond to events in the input
  observeEvent(input$btngen,{
    A <<- eig.make3x3()
    trA <<- sum(diag(A))
    detA <<- det(A)
    discA <<- trA^2-4*detA
    
    #calculate the eigen values
    eig_A = eigen(A)
    eig_A_vals = eig_A$values
    
    #calculate the inverse
    inv_A = solve(A)
    
    #calculate the values for t = 0...4
    det(A)
    a0 = det(A - eig_A_vals[1]^0 * inv_A)
    a1 = det(A - eig_A_vals[1]^1 * inv_A)
    a0 = det(A - eig_A$vectors^0 * inv_A)
    
    # a1 = det(A - 1 * inv_A)
    # a2 = det(A - 2 * inv_A)
    # a3 = det(A - 3 * inv_A)
    # a4 = det(A - 4 * inv_A)
    
    #Might need to generalize the eigen.diagonize and eigen.nthpower functions
    
    #returning eigen values
    output$eigen_A = renderText({  eig_A_vals })
    
    output$plotA <- renderPlot(eig.plot(A,-1))
    output$eigval <- renderUI("")
    output$eigvec <- renderUI("")
    output$diagtbl <- renderUI("")
    output$sqrtbl <- renderUI("")
    output$check <- renderUI("")
  })
  
  
  
  # observeEvent(input$btndiag,{
  #   x <- eig.diagonalize(A)
  #   D <<- x$diag
  #   P <<- x$cofb
  #   lam1 <- D[1,1]
  #   lam2 <- D[2,2]
  #   v1 <<- P[,1]
  #   v2 <<- P[,2]
  #   Pinv <<- solve(P)
  #   output$eigval <- renderUI(jaxD(paste("`lambda_1 =", round(lam1, digits = 3),
  #                                        ",`;   `lambda_2 =",round(lam2, digits = 3),
  #                                        ",`;   `lambda_1+`lambda_2 =",round(lam1+lam2, digits = 3),
  #                                        ",`;   `lambda_1`lambda_2 =",round(lam1*lam2, digits = 3)
  #   )))
  #   output$eigvec <- renderUI(withTags(
  #     table(
  #       tr(
  #         th(jax.vector(round(v1,digits = 3), name = "v_1")),
  #         th(jax.vector(round(v2,digits = 3), name = "v_2"))
  #       ),
  #       tr(
  #         th(jax.vector(round(A%*%v1,digits = 3), name = "Av_1")),
  #         th(jax.vector(round(lam1*v1,digits = 3), name = ",\\;\\lambda_1v_1")),
  #         th(jax.vector(round(A%*%v2,digits = 3), name = "Av_2")),
  #         th(jax.vector(round(lam2*v2,digits = 3), name = ",\\;\\lambda_2v_2")),
  #       )
  #     )
  #   ))
  #   output$diagtbl <- renderUI(withTags(
  #     table(border = 0,
  #           tr(
  #             th(""),
  #             th(jax.matrix(round(D, digits = 2), name = "D")),
  #             th(jax.matrix(round(Pinv, digits = 2), name = "P^{-1}"))
  #           ),
  #           tr(
  #             th(jax.matrix(round(P, digits = 2), name = "P")),
  #             th(jax.matrix(round(P%*%D, digits = 2))),
  #             th(jax.matrix(round(P%*%D%*%Pinv, digits = 2)))
  #           )
  #     )))
  #   
  # })
  # observeEvent(input$btnsqrt, {
  #   Dsqrt <- diag(sqrt(diag(D)))
  #   B <<- P%*%Dsqrt%*%Pinv 
  #   output$sqrtbl <- renderUI(withTags(
  #     table(
  #       tr(
  #         th(""),
  #         th(jax.matrix(round(Dsqrt, digits = 2), name = "D^{0.5}")),
  #         th(jax.matrix(round(Pinv, digits = 2), name = "P^{-1}"))
  #       ),
  #       tr(
  #         th(jax.matrix(round(P, digits = 2), name = "P")),
  #         th(jax.matrix(round(P%*%Dsqrt, digits = 2))),
  #         th(jax.matrix(round(P%*%Dsqrt%*%Pinv, digits = 2)))
  #       )
  #     )))
  #   output$check <- renderUI(withTags(
  #     table(
  #       tr(
  #         th(""),
  #         th(jax.matrix(round(B, digits = 2), name = "A^{0.5}")),
  #       ),
  #       tr(
  #         th(jax.matrix(round(B, digits = 2), name = "A^{0.5}")),
  #         th(jax.matrix(round(B%*%B, digits = 2))),
  #       )
  #     )))
  #   
  # })
  # observeEvent(input$slpower,{
  #   q <- input$slpower
  #   if(detA > 0)
  #     output$plotA <- renderPlot(eig.plot(A,q))
  # })
}

#Run the app
shinyApp(ui = ui, server = server)
