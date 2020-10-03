#install.packages("tseries")
#install.packages("quantmod")
#install.packages("rugarch")
#install.packages("DEoptim")
#install.packages("VineCopula")
#install.packages("plotly")

library(tseries)
library(quantmod)
library(rugarch)
library(DEoptim)
library(VineCopula)
library(plotly)


#####################################################
###Kapitel 2: Autokorrelation von Finanzzeitreihen###
#####################################################
MSCI = getSymbols("MSCI",sc = "yahoo",from="2007-11-01", env=NULL)
MSCI_returns = dailyReturn(MSCI, type="log")

# Manuelle Berechnung Autokorrelation des MSCI
cor(MSCI_returns[1:length(MSCI_returns)-1], MSCI_returns[2:length(MSCI_returns)])

# Untersuchung AR(1,0) auf MSCI_returns
ar_fit_MSCI = arma(MSCI_returns, order = c(1,0), include.intercept = TRUE)
ar_fit_MSCI

######################################################################
###Kapitel 3: Performanceanalyse mittels eines Copula-GARCH-Modells###
######################################################################

# Finanzdaten aus yahoo! laden
symbolList = c("MSCI", "DBA", "SLV")
symbolData = new.env()
getSymbols(symbolList, sc="yahoo", from="2008-01-022", env=symbolData)
Assets = do.call(merge, eapply(symbolData, Cl))
Assets = na.omit(Assets)

# Festlegen der Parameter
risk_level = 0.10
loss_level = 0.10
n_window = 750 
n_forecast = 100 
n_paths = 1000


n_data = length(Assets[,1])
n_assets = dim(Assets)[2]
w_Assets = array(0, c(n_data,dim(Assets)[2]+1))
v_Assets = array(0,n_data)



v_Assets[1:n_window] = Assets$MSCI.Close # die ersten 750 Tage ohne Aktion auf den MSCI festhalten
rebalancing_days = seq(n_window+1, n_data-n_forecast, n_forecast) # alle 100 Tage ab 751 können Portfolioanteile geändert werden
portfolio_return_VaR = array(0, length(v_Assets))

for (i in (n_window+1):n_data)
{
  if(i %in% rebalancing_days){ # GARCH-Modell nur zu schätzen wenn i = rebalancing day
    tryCatch({
      # manuelle Berechnung der log-returns
      assets_window = Assets[(i-n_window):(i-1),] 
      log_returns = array(0,c(n_window-1,dim(Assets)[2])) 
      for (k in 2:n_window){
        log_returns[k-1,] = log(as.numeric(assets_window[k,])) - log(as.numeric(assets_window[k-1,])) ## manuelle Berechnung der log-returns innerhalb des rolling windows
      }
      # Schätzung des GARCH-Modells
      garch.setup = ugarchspec(mean.model=list(armaOrder=c(0,0), include.mean=TRUE), variance.model=list(model="sGARCH", garchOrder=c(1,1)), distribution.model="sstd")
      fit_vector = list() 
      residuals_vector = array(0,c(dim(assets_window)[1]-1, dim(Assets)[2]))
      # Fitten des GARCH-Modells
      for (k in 1:dim(Assets)[2]){ 
        fit_vector[[k]] = ugarchfit(garch.setup, data = log_returns[,k], solver = "hybrid") 
        residuals_vector[,k] = residuals(object=fit_vector[[k]], standardize = TRUE) 
      }
      udata = matrix(nrow = dim(residuals_vector)[1], ncol = dim(Assets)[2])  
      # Gleichverteilte Daten für Copula
      for (k in 1:dim(Assets)[2]){ 
        udata[,k] = pdist(distribution = "sstd", residuals_vector[,k], skew = coef(fit_vector[[k]])[5], shape = coef(fit_vector[[k]])[6])
      }
      # Bivariate Copula Schätzen
      copula_structure = RVineStructureSelect(udata, familyset=c(2,4,14,24,34)) 
      sim_udata = RVineSim(n_paths * n_forecast, copula_structure) 
      
      sim_prices = array(0,c(n_paths, dim(Assets)[2])) 
      # Simulierte Preise der Indices schätzen
      for (k in 1:dim(Assets)[2]){
        sim_sstd_vector = matrix(qdist(distribution = "sstd", sim_udata[,k], skew = coef(fit_vector[[k]])[5], shape = coef(fit_vector[[k]])[6]), nrow=n_forecast, ncol=n_paths) ## 100 x 1000 Simulationen aus Student-t (nur für ein k-Asset)
        sim_log_returns = ugarchsim(fit_vector[[k]], n.sim=n_forecast, n.start=0, m.sim=n_paths, startMethod="sample", custom.dist = list(name = "sample", distfit = sim_sstd_vector))
        sim_prices[,k] = as.numeric(last(assets_window)[,k]) * exp(colSums(fitted(sim_log_returns)))
      }
      
      sim_prices = sim_prices[complete.cases(sim_prices), ]
      sim_prices_t0 = as.numeric(last(assets_window))
      
      # Optimierungsfunktion
    }, error=function(e) e, warning=function(w) w)
    
    full_optim_function = function(weights){
      portfolio_return = array(0, dim = dim(sim_prices)[1])
      for (k in 1:dim(sim_prices)[1]){
        portfolio_return[k] = sum(weights * (sim_prices[k,] - sim_prices_t0) / sim_prices_t0 + (1-sum(weights)) * 0) ## + cash 
      }
      r = -mean(portfolio_return) + 1e3 * max(-quantile(portfolio_return, risk_level, type = 1)*sqrt(250/n_forecast) - loss_level, 0) + 1e3 * max(sum(weights)-1,0)
      return(r)
    }
     # Ergebnisse der Optimierung
    outDEoptim_results = DEoptim(full_optim_function, rep(0, dim(Assets)[2]), rep(1,dim(Assets)[2]), DEoptim.control(NP = dim(Assets)[2]*10, itermax = 50, trace = FALSE)) ## Paper lesen
    
    w_Assets[i,] = c(1-sum(outDEoptim_results$optim$bestmem), as.numeric(outDEoptim_results$optim$bestmem)) #Abspeichern der besten Ergebnisse (erster Teil Cash)
    
    # Berechnung des VaR für jedes Copula-GARCH-Zeitintervall mithilfe der simulierten Preise von oben
    portfolio_return_sim = array(0, dim = dim(sim_prices)[1])
    for (k in 1:dim(sim_prices)[1]){
      portfolio_return_sim[k] = sum(w_Assets[i,] * c(0, ((sim_prices[k,] - sim_prices_t0) / sim_prices_t0)))
    }
    # Berechnung des VaR
    portfolio_return_VaR[i] = quantile(portfolio_return_sim, risk_level, type = 1)
    
  }else{}
   # Berechnung des Portfoliowertes
  if(!w_Assets[i,1]==0){
    u_Assets = as.numeric(w_Assets[i,]) * v_Assets[i-1] / c(1,Assets[i-1,])
  }else{}
  
  v_Assets[i] = sum(u_Assets * c(1,Assets[i,])) 
  
  cat("current step: ", i, "\n")
  
}

# Plot des Portfoliowertes und einer gleichgewichteten Strategie im Zeitverlauf

Strategie = Assets[,1]
for(i in 1:n_data)
{
  Strategie[i] = v_Assets[i]
}

gleichgewichtete_Strategie = Strategie
for(i in (n_window+1):n_data)
{
  gleichgewichtete_Strategie[i] = gleichgewichtete_Strategie[i-1] * (1 + mean((as.numeric(Assets[i,])-as.numeric(Assets[i-1,]))/as.numeric(Assets[i-1,]))) ## gleichgewichtete Strategie
}

chartSeries(Assets$MSCI.Close)

addTA(Strategie, on=1, col = "red")
addTA(gleichgewichtete_Strategie, on=1, col = "yellow")


# Betrachtung der Anteile im Zeitverlauf

w_SLV = array(0,length(w_Assets[,2]))
for (i in 1:length(w_Assets[,2]))
{
  if (!w_Assets[,2][i] == 0){
    w_SLV[i:length(w_SLV)] = w_Assets[,2][i]
  }else {}
}


w_MSCI = array(0,length(w_Assets[,3]))
for (i in 1:length(w_Assets[,3]))
{
  if (!w_Assets[,3][i] == 0){
    w_MSCI[i:length(w_MSCI)] = w_Assets[,3][i]
  }else {}
}
w_MSCI[1:750] = 1

w_DBA = array(0,length(w_Assets[,4]))
for (i in 1:length(w_Assets[,4]))
{
  if (!w_Assets[,4][i] == 0){
    w_DBA[i:length(w_DBA)] = w_Assets[,4][i]
  }else {}
}


w_cash = array(0,length(w_Assets[,1]))
for (i in 1:length(w_Assets[,1]))
{
  if (!w_Assets[,1][i] == 0){
    w_cash[i:length(w_cash)] = w_Assets[,1][i]
  }else {}
}

Zeit = time(Assets)

Gewichte = data.frame(Zeit,w_cash, w_SLV, w_MSCI, w_DBA)

# Plot der Anteile im Zeitverlauf

fig = plot_ly(Gewichte, x = ~Zeit, stackgroup = "one", groupnorm = 'percent')
fig = fig %>% add_trace(y = ~w_cash, name = "Geld", mode = "lines", type = "scatter")
fig = fig %>% add_trace(y = ~w_SLV, name = "SLV", mode = "lines", type = "scatter")
fig = fig %>% add_trace(y = ~w_MSCI, name = "MSCI", mode = "lines", type = "scatter")
fig = fig %>% add_trace(y = ~w_DBA, name = "DBA", mode = "lines", type = "scatter")
fig = fig %>% layout(title = "Anteile der Assets im Zeitverlauf", 
                     xaxis = list(title = "Jahr",
                                  showgrid = FALSE),
                     yaxis = list(title = "Anteile", 
                                  showgrid = FALSE,
                                  ticksuffix = "%"))

fig



##########################################
###Kapitel 4: Value-at-Risk Backtesting###
##########################################

# Value-at-Risk im Zeitverlauf anschaulicher machen

portfolio_return_VaR2 = array(0,length(portfolio_return_VaR)) 

for (i in 1:length(portfolio_return_VaR))
{
  if (!portfolio_return_VaR[i] == 0){
    portfolio_return_VaR2[i:length(portfolio_return_VaR)] = portfolio_return_VaR[i]
  }else {}
}

# Berechnung der Portfoliorenditen

return_Assets = array(0, dim(portfolio_return_VaR2))

for (i in 1:length(v_Assets))
{
  return_Assets[i+1] = (v_Assets[i+1] / v_Assets[i]) - 1 
}
return_Assets = return_Assets[1:length(return_Assets)-1]


# Plot der Überschreitungen (ersten 750 Tage ignorieren)
plot(as.numeric(return_Assets[751:length(return_Assets)]) < as.numeric(portfolio_return_VaR2[751:length(return_Assets)]), type = "l")

# Anzahl Überschreitungen (ersten 750 Tage ignorieren)
table(as.numeric(return_Assets[751:length(return_Assets)]) < as.numeric(portfolio_return_VaR2[751:length(return_Assets)]))

# Berechnug des Anteils der Überschreitungen (ersten 750 Tage ignorieren)
sum(as.numeric(return_Assets[751:length(return_Assets)]) < as.numeric(portfolio_return_VaR2[751:length(return_Assets)])) / (length(return_Assets)-750)


# Plot zum Vergleich der Portfoliorendite mit dem VaR

VaR_obs = data.frame(Zeit,portfolio_return_VaR2, return_Assets)

fig2 = plot_ly(VaR_obs, x = ~Zeit)
fig2 = fig2 %>% add_trace(y = ~return_Assets, name = "Portfoliorendite", mode = "lines", type = "scatter")
fig2 = fig2 %>% add_trace(y = ~portfolio_return_VaR2, name = "VaR",  mode = "lines", type = "scatter")
fig2 = fig2 %>% layout(title = "", 
                       xaxis = list(title = "Jahr",
                                    showgrid = FALSE),
                       yaxis = list(title = "Rendite", 
                                    showgrid = FALSE
                       ))

fig2


# VaR Backtests (Kupiec und Christoffersen)
VaRTest(alpha = 0.10, return_Assets[751:length(return_Assets)], portfolio_return_VaR2[751:length(portfolio_return_VaR2)], conf.level = 0.95)



