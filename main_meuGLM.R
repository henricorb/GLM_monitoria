# Código para praticar Modelos Lineares Generalizados
# Neste código espera-se:
#   - Simular dados de distribuições Normal, Poisson e Binomial que tenham 
#     modelos linear no parâmetro de interesse (média, taxa, proporção).
#   - Estimar os parâmetros que foram simulados
#   - Ilustrar o processo de otimização



# Funções auxiliares

# ----------------------------meuGLM--------------------------------------------
meuGLM <- function(Y, X, dist_name, link_func, tol = 1e-6, max_iter = 200){
  # meuGLM é uma função para estimar os coeficientes dos regressores de um 
  # modelo linear generalizado transversal.
  # regressor  = variável independente = preditor = etc.
  # Entradas
  #   Y: vetor Nx1 com as observações da variável dependente
  #   X: matriz Nxp, com os valores dos p regressores (não inclui a constante) 
  #   correspondentes de cada observação.
  #   dist_name: nome da distribuição assumida, pode ser "normal", "binomial"
  #   ou "poisson".
  #   link_func: nome da função de ligação utilazada. Pode ser "identity", "log",
  #   ou "logit".
  #   tol: tolerância de mudança na magnitude entre uma iteração e a anterior do
  #   vetor de betas para parar o algoritmo de otimização
  #   max_iter: máximo de iterações para a otimização
  #
  # Saídas
  #   beta: vetor (p+1)x1 com os pesos dos regressores b0, b1, ..., bp
  #   conv_err: tipo do erro. 0 - sem erro, 1 - máximo de iterações atingida,
  #   2 - outro tipo de erro.

  # Checar erros de entrada. Fazer depois.
  # Falta calcular a log-likelihood e registrar
  # Falta simular os dados
  
  # Constantes úteis
  N <- length(Y) #número de amostras
  if(is.matrix(X)==FALSE){
    p <- 1
  }else{
    p <- dim(X)[2] #número de regressores sem o b0  
  }
  conv_err <- 0 #sem erro de convergência no começo
  n_binomial <- 1
  if (dist_name=="binomial"){
    n_binomial <- length(unique(Y))-1 #se baseando no número de classes, mas 
                                      #acho que está errado
  }
  
  # Ajusta matriz de regressores adicionando os '1' do b0
  aux <- rep(1,N)
  X <- cbind(aux,X)
  
  # Calcula chute inicial dos betas.
  beta <- numeric(p+1)
  beta[1] <- mean(Y)

  # Faz processo iterativo para cálculo dos betas
  beta_history <- matrix(0, p+1, max_iter) #histórico de evolução dos betas
  beta_history[,1]<-beta 
  l_history <- numeric(max_iter)
  l_history[1] <- -Inf
  
  improv <- 1
  niter <- 1
  while((improv > tol) && (niter<max_iter)){ 
    # Enquanto houver melhoria expressiva no log-likelihood e não tiver atingido 
    # o máximo de iterações
    
    beta_prev <- beta_history[,niter]
    eta <- X%*%beta
    
    # Calculando componentes que dependem da escolha da função de ligação
    if (link_func=="identity"){
      mu <- eta # Parâmetro da distribuição
      d_mu <- rep(1,N) # Derivada do parâmetro
      
    }else if (link_func=="log"){
      mu <- exp(eta) # Parâmetro da distribuição
      d_mu <- exp(eta) # Derivada do parâmetro
      
    }else if (link_func=="logit"){
      mu <- 1/(1+exp(-eta)) # Parâmetro da distribuição
      d_mu <- exp(-eta)/((1+exp(-eta))^2) # Derivada do parâmetro
      
    }else{
      errorCondition("Link function not available")
    }
    mu <- c(mu) #garantindo que mu será um vetor
    d_mu <- c(d_mu) #idem para d_mu
    
    # Calculando componentes que dependem da esoclha de distribuição
    # Variância de Y_i
    if(dist_name=="normal"){
      var_y <- rep((1/N)*sum((Y-eta)^2),N) #variância estimada comparando as
                                           #observações à "reta" de regressão
        
    }else if(dist_name=="poisson"){
      var_y <- mu #média e variância são iguais na Poisson
      
    }else if(dist_name=="binomial"){
      var_y<-n_binomial*mu*(1-mu)
      
    }else{
      errorCondition("Distribution not available")
    }
    var_y <- c(var_y) #garantindo que var_y será um vetor
    
    # Calculando o vetor de escore, U, e a matriz de informação,  J
    U<-0*numeric(p+1)
    for (j in 1:(p+1)){
      U[j]<-0
      for(i in 1:N){
        U[j]<-U[j] + X[i,j]*d_mu[i]*(Y[i]-mu[i])/var_y[i]
      }
    }
    
    # Matriz auxiliar W para criar a matriz de informação J
    W <- diag((1/var_y)*(d_mu^2))
    
    # Matriz de informação J
    J <- t(X)%*%W%*%X
    
    # Atualizando o vetor de parâmetros beta
    beta <- beta_prev + solve(J,U)
    
    # Atualizando iterador e registros
    niter <- niter+1
    beta_history[,niter]<-beta
    
    # Calculando log-likelihood
    sum_loglike <- 0
    for(i in 1:N){
      
      if(dist_name=="normal"){
        l <- Y[i]*(mu[i]/var_y[i]) + 
             (-mu[i]/(2*var_y[i]) -0.5*log(2*pi*var_y[i])) +  
             -Y[i]/(2*var_y[i])
        
      }else if(dist_name=="poisson"){
        l <- Y[i]*log(mu[i]) +
             -mu[i] +
             -log(factorial(Y[i]))
        
      }else if(dist_name=="binomial"){
        l <- Y[i]*log(mu[i]/(1-mu[i])) + 
             n_binomial*log(1-mu[i]) +
             log(choose(n_binomial,Y[i]))
        
      }else{
        errorCondition("Distribution not available")
      }
      
      sum_loglike <- sum_loglike + l
    }
    l_history[niter]<-sum_loglike
    
    improv <- abs(sum_loglike - l_history[niter-1])
  }
  
  if (niter>=max_iter){
    conv_err <- 1
  }
  
  AIC <- 2*(p+1) -2*sum_loglike
  l_history <- l_history[1:niter]
  beta_history <- beta_history[,1:niter]
  optim_data<-list(l_history = l_history, beta_history = beta_history)
  
  # Calculando intervalos de confiança 95% para os betas
  beta_confint <- matrix(NA,nrow = p+1,ncol = 2) #cada linha é um parâmetro
  err_pad <- sqrt(diag(solve(J))) #vetor com os erros-padrão de cada beta
  beta_confint[,1] <- beta - 1.96*err_pad
  beta_confint[,2] <- beta + 1.96*err_pad
  
  return(list(beta = beta, 
              beta_confint = beta_confint,
              AIC = AIC,
              conv_err = conv_err, optim_data = optim_data))
}



# ----------------------------geraDados-----------------------------------------
geraDados <- function(N, beta, X, dist_name, s2){
  # Função geraDados cria dados simulados de distribuições normal, binomial e 
  # poisson a partir do modelo linear dado pelos coeficientes em beta e das
  # entradas/Variáveis independentes X.
  #
  # Entradas
  #   N: número de observações para gerar
  #   beta: vetor (p+1)x1 com os pesos b0, b1, ..., bp
  #   X: matriz Nxp com os regressores de cada observação. Não possui a coluna de 
  #   '1's do beta0.
  #   dist_name: o nome da distribuição para usar, pode ser "normal", "poisson",
  #              ou "binomial"
  #   s2: parâmetro auxiliar da distribuição. É a variância no caso da normal e
  #       o número de jogadas n no caso da binomial
  
  # Ajusta matriz de regressores adicionando os '1' do b0
  aux <- rep(1,N)
  X <- cbind(aux,X)
  
  Y <- numeric(N) #inicializa vetor de observações
  eta = X%*%beta #variavel auxiliar para as predições
  
  if (dist_name=="normal"){
    mu = eta
    for (i in 1:N){
      Y[i]<-rnorm(1,mu[i],s2)
    }
  } else if (dist_name=="poisson"){
    mu = exp(eta)
    for (i in 1:N){
      Y[i]<-rpois(1,mu[i])
    }
    
  } else if (dist_name=="binomial"){
    mu <- 1/(1+exp(-eta))
    for (i in 1:N){
      Y[i]<-rbinom(1,s2,mu[i])
    }
    
  } else{
    errorCondition("Distribution not available")
  }
  
  return(Y)
}




#========================= Script principal=====================================
# Exemplo 1) Qtd de pelos soltos vs nível de raiva do gato
N <- 60 # 60 dias de observação
X<-sample(0:2, N, replace = TRUE) #3 níveis de raiva
beta<-c(10,50) #perde em média 100 gramas e aumenta em 50 por nível de raiva

dist_name <- "normal"
link_func <- "identity"
Y <- geraDados(N,beta,X,"normal",s2 = 5^2)
plot(X,Y, xlab = "Nível de Raiva", ylab = "Pelos perdidos (g)")

result <- meuGLM(Y,X,dist_name,link_func)
result$beta
result$beta_confint
result$AIC
optim <- result$optim_data
plot(optim$l_history)


# Exemplo 2) Chance de sucesso de Jonatan pegar uma cocotinha na festa, a 
# depender se a cantada foi ruim ou boa
N <- 100 # 100 tentativas de observação
X<-sample(0:1, N, replace = TRUE) #0- cantada ruim, 1 - cantada boa
beta<-c(0,log(2)) #cantadas boas dobram as chances. Chance base é 1:1

dist_name <- "binomial"
link_func <- "logit"
Y <- geraDados(N,beta,X,"binomial",s2 = 1)
hist(Y[X==0],main = "Histograma com Cantadas Ruins")
hist(Y[X==1],main = "Histograma com Cantadas Boas")

result <- meuGLM(Y,X,dist_name,link_func)
exp(result$beta)
exp(result$beta_confint)
result$AIC
optim <- result$optim_data
plot(optim$l_history)



# Exemplo 3) Média de pistoladas/rages que o professor Altay dá em aula. Variável
# independente é ele errar algo ou perceber letra feia durante a aula
N <- 120 # 120 aulas (2 anos) observados
X<-sample(0:1, N, replace = TRUE) #0- sem cavalada, 1 - com cavalada
beta<-c(log(4),log(1.5)) #Em média 4 rages por aula, e aumenta em 1.5 com cavaladas

dist_name <- "poisson"
link_func <- "log"
Y <- geraDados(N,beta,X,dist_name,s2 = 1)
plot(X,Y, xlab = "Cavalada na aula", ylab = "Qtd de Rages")

result <- meuGLM(Y,X,dist_name,link_func)
exp(result$beta)
exp(result$beta_confint)
result$AIC
optim <- result$optim_data
plot(optim$l_history)