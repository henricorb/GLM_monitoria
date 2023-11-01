# Código para ilustrar otimização iterativa com método de Newton-Raphson


# Exemplo 1) f(x) = x(x-2)(x+2) = x^3 - 4x
fx <- function(x){
  return(x^3-4*x)
}
dfx <- function(x){
  return(3*(x^2) - 4)
}
d2fx <- function(x){
  return(6*x)
}

# # Exemplo 2) f(x) = exp(-x^2)
# fx <- function(x){
#   return(exp(-x^2))
# }
# dfx <- function(x){
#   return((-2*x)*exp(-x^2))
# }
# d2fx <- function(x){
#   return((4*x^2-2)*exp(-x^2))
# }

# ==========================Ilustrando os exemplos==============================
# Plotando o gráfico
x<-seq(from=-3, to=3, by=0.01)
f<-fx(x)
plot(x,f,type='l')

# Otimização
x0<- 0.1 #chute inicial
xk <- x0
tol <- 1e-6
max_iter <- 100
niter <- 0
x_hist <- numeric(max_iter)
improv <- +Inf
while (niter < max_iter && improv > tol){
  d <- dfx(xk)   #primeira derivada
  d2 <- d2fx(xk) #segunda derivada
  x_next <- xk - d/d2
  
  niter <- niter+1
  x_hist[niter] <- x_next
  improv <- abs(fx(x_next)-fx(xk))
  xk <- x_next
}
x_hist <- c(x0,x_hist[1:niter])

# Exibindo as etapas da otimização
plot(x,fx(x),type = 'l')
points(x_hist,fx(x_hist),col='blue')



# Exibindo as etapas de otimização
# x <- max(abs(x_hist)) #comente para o primeiro exemplo
# x <- seq(from=-x, to=+x, by=0.01) #comente para o primeiro exemplo
for (k in 1:niter){
  xk <- x_hist[k]
  xn <- x_hist[k+1]
  # Aproximação de segunda ordem da função em torno de xk
  fit_func <- fx(xk) + dfx(xk)*(x-xk) + 0.5*d2fx(xk)*((x-xk)^2)
  fit_func_new <- fx(xk) + dfx(xk)*(xn-xk) + 0.5*d2fx(xk)*((xn-xk)^2)
  
  plot(x,fx(x),'l')
  lines(x,fit_func,col='blue')
  points(x_hist[k+1],fit_func_new,col='red')
  points(xk,fx(xk))
}