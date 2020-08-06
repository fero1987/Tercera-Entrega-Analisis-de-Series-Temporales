############# TERCER ENTREGA - ANÁLISIS DE SERIES TEMPORALES #############
# Alumno: Martínez-De Elia-Díaz Ugalde

suppressPackageStartupMessages({
  library(tseries)
  library(forecast)
  library(ggplot2)
  library(gridExtra)
  library(car)
  library(nortest)
  library(AdequacyModel)
  library(lmtest)
  library(quantmod)
  library(dygraphs)
  library(lessR)
  library(forecast)
  library(ggplot2)
  library(gridExtra)
  library(PASWR2)
  library(dplyr)
  library(psych)
  library(pastecs)
  library(astsa)
  library(lessR)
  library(tseries)
  library(zoo)
  library(xts)
  library(fma)
  library(expsmooth)
  library(lmtest)
  library(Quandl)
  library(fpp)
  library(urca)
  library(AER)
  library(fUnitRoots)
  library(CADFtest)
  library(fpp2)
  library(car)
  library(tseries)
  library(gridExtra)
  library(datasets)
  library(fpp)
})


  # Limpio la memoria
rm( list=ls() )
gc()

# Elijo la siguiente serie del paquete fma. La misma representa el precio de la acción de IBM al cierre
base<-fma::ibmclose

# Grafico la serie
ts<-base
autoplot(ts)+ ggtitle("Precio Acción IBM") + ylab("")

# Grafico FAS, FAC y FACP 
acf(ts,type = "covariance",plot = T)
ggAcf(ts) + ggtitle("FAC Serie IBM") # La serie decrece linealmente, lo que indica que no es estacionaria
ggAcf(ts,type = "partial") + ggtitle("FACP Serie IBM")

# Una forma de comprobar la estacionariedad es detectar si todos los coeficientes de autocorrelación son iguales a cero

# Cargo la siguiente función de incorrelación que realiza un test de Ljung-Box o Box-Pierce para distintos lags
Incorrelation <- function(ts, type = c("Ljung-Box","Box-Pierce"), fitdf = 0){
  p_ljung_box = NULL
  s_ljung_box = NULL
  for(i in 0:(length(ts)/4)){
    p_ljung_box[i] = Box.test(ts,lag = i,type = type,fitdf = fitdf)$p.value
    s_ljung_box[i] = Box.test(ts,lag = i,type = type,fitdf = fitdf)$statistic
  }
  table = data.frame(j = 1:(length(ts)/4),
                     P_Value = p_ljung_box,
                     Statistic = s_ljung_box)
  return(table)
}

# Planteo el test de Ljung-Box. Si rechazo H0 significa que hay coeficientes de autocorrelación distintos a cero
Incorrelation(ts,"Ljung-Box")

inco_wn = Incorrelation(ts,"Ljung-Box")

# Grafico los p-value para disntitos lags
autoplot(ts(inco_wn$P_Value)) +
  ggtitle("Test de Ljung-Box", subtitle = "P-Value") +
  ylab("")
# El p-value se ubica en 0 por lo que rechazo H0 y puedo considerar que es una serie no estacionaria 

# Otra forma es plantear la estacionariedad es hacer el test de Dickey-Fuller
adf.test(ts) # p-value = 0.701. No es estacionaria. Lag = 7

# La siguiente función que automatiza el adf. fuente: http://www.systerdyster.se/rrrabbit/code.php?id=1
# Este test utiliza el test base del paquete urca. Lo que hace es hacer el adf test para distintos 
#lags, quedándose con el más significativo. Luego realiza el adf con trend, drift y none para el lag
# seleccionado

aadf(ts, lags=10) # Ver script de la función al final
summary(ur.df(ts, type = "drift", lags = 6))
summary(ur.df(ts, type = "trend", lags = 6))
summary(ur.df(ts, type = "none", lags = 6))

# El estadístico tau1, tau2 y tau3 on mayores a los valores críticos del por lo que no se puede rechazar H0
# El módulo de los estadísticos phi2 y phi3 no son mayores en módulo al valor crítico. No son significativos. 

# Calculo la cantidad de diferencias para poder hacer la serie estacionaria
ndiffs(ts) # Una diferencia
ts1diff <- diff(ts)
autoplot(ts1diff)+ ggtitle("Precio Acción IBM (1 diferencia)") + ylab("")

# Grafico FAS, FAC y FACP 
acf(ts1diff,type = "covariance",plot = T)
ggAcf(ts1diff) + ggtitle("FAC Serie dif IBM") 
ggAcf(ts1diff,type = "partial") + ggtitle("FACP Serie dif IBM")

adf.test(ts1diff) # p-value = 0.01. Es estacionaria. Lag = 7
aadf(ts1diff, lags=10)
summary(ur.df(ts, type = "drift", lags = 6))
summary(ur.df(ts, type = "trend", lags = 6))
summary(ur.df(ts, type = "none", lags = 6))

#---------------------------------------------- 0 ----------------------------------------------
# FUNCION AADF

#Este es una función que automatiza el adf. fuente: http://www.systerdyster.se/rrrabbit/code.php?id=1
# Este test utiliza el test base del paquete urca. Lo que hace es hacer el adf test para distintos 
#lags, quedándose con el más significativo. Luego realiza el adf con trend, drift y none para el lag
# seleccionado


# Automated Augmented Dickey-Fuller Test for R
# v.1 - Andreas Jangmo - andreas(at)systerdyster.se
# This function depends entirely on ur.df() in the urca-library and
# performs no statistical computations of its own except comparisons of
# critical values and test statistics given from ur.df().
# Arguments:
#	data - Can be a vector or matrix (not all types have been tested).
#	lags - Lag length. Default is the same as for ur.df()
#	slctlags - Criteria for selecting number of lags given the chosen
#		lag length. Defaults to "Fixed".
#	significance - Rejection level for selecting lag length.
#	fulltest - Boolean: TRUE is default; after searching a proper lag
#		lag length the program will continue with unit root test with
#		trend, drift ano no trend/drift.
#	thin - Boolean: TRUE reduces output from the test procedure,
#		restricting it to last significant lag in finding the proper
#		lag length and to test statistics fro tao and phi in the
# 		subsequent tests.
# Note that if the automatic lag selection methods AIC/BIC are to be
# used, thin should be set to FALSE.
aadf = function(data, lags=1, slctlags=c("Fixed", "AIC", "BIC"), significance=0.05, fulltest=TRUE, thin=TRUE){
  cat("#########################################################
      #   Automated Augmented Dickey-Fuller Unit Root Test    #
      #########################################################\n");
  info = "Null-Hypotheses:
  tao: unit root
  phi2: unit root, no drift, no time trend
  phi3: unit root, no time trend
  phi1: unit root, no drift
  Significance levels:
  ** p < 1 pct, * p < 5 pct\n";
  cols = ncol(data);
  if (is.null(cols)){cols = 1; cname = "data";} else {cname = colnames(data);}
  for (j in 1:cols){
    cat("\nOOOOOoooo... ADF for", cname[j], ".oOo.\n");
    lastlag = 0;
    if (cols > 1){tdata = data[,j];} else {tdata = data;}
    if (length(slctlags) == 3){
      cat("---------------------------------------\n")
      cat("P-value for lag ")
      for (i in 1:lags){
        cat(i, " ");
        temp = ur.df(tdata, lags=i, type="trend")@testreg$coefficients[3+i, 4]
        cat(temp);
        if (temp < significance){
          cat("*");
          lastlag = i;
        }
        if (i != lags){cat("\n................")}
      }
      cat("\n---------------------------------------\n")
    } else {
      i = lags;
      temp = ur.df(tdata, lags=i, type="trend", selectlags=slctlags)@testreg$coefficients;
      lastlag = nrow(temp) - 3;
      if (is.na(lastlag)){lastlag = 1;}
      if (temp[nrow(temp), 4] > significance){lastlag = 0; cat("Last lag is not significant!\n");}
    }	
    if (fulltest & lastlag != 0){
      cat("Last lag significant is at lag length", lastlag, "\b, proceeding to output trend full test,\ndrift and no drift/trend in", cname[j], "\b.\n\n");
      if (!thin){
        print(summary(ur.df(tdata, lags=lastlag, type="trend", selectlags=slctlags)));
        print(summary(ur.df(tdata, lags=lastlag, type="drift", selectlags=slctlags)));
        print(summary(ur.df(tdata, lags=lastlag, type="none", selectlags=slctlags)));
      } else {
        cat(info);
        cat("\nTrend\n");
        temp = ur.df(tdata, lags=lastlag, type="trend", selectlags=slctlags);
        tstat = round(temp@teststat, 2);
        print(t(temp@cval));
        cat("---------------------\nTest:", tstat[1], tstat[2], tstat[3], "\n");
        signstr = "        ";
        if (tstat[1] < temp@cval[1,2]){signstr = paste(signstr, "*", sep="");} else {signstr = paste(signstr, " ", sep="");}
        if (tstat[1] < temp@cval[1,1]){signstr = paste(signstr, "*", sep="");} else {signstr = paste(signstr, " ", sep="");}
        if (tstat[2] > temp@cval[2,2]){signstr = paste(signstr, "   *", sep="");} else {signstr = paste(signstr, "    ", sep="");}
        if (tstat[2] > temp@cval[2,1]){signstr = paste(signstr, "*", sep="");} else {signstr = paste(signstr, " ", sep="");}
        if (tstat[3] > temp@cval[3,2]){signstr = paste(signstr, "   *", sep="")} else {signstr = paste(signstr, "    ", sep="");}
        if (tstat[3] > temp@cval[3,1]){signstr = paste(signstr, "*", sep="");}
        cat(signstr);
        cat("\n\nDrift\n");
        temp = ur.df(tdata, lags=lastlag, type="drift", selectlags=slctlags);
        tstat = round(temp@teststat, 2);
        print(t(temp@cval));
        cat("----------------\nTest: ", tstat[1], tstat[2], "\n");
        if (tstat[1] < temp@cval[1,2]){cat("        *")}
        if (tstat[1] < temp@cval[1,1]){cat("*")}
        if (tstat[2] > temp@cval[2,2]){cat("   *")}
        if (tstat[2] > temp@cval[2,1]){cat("*")}
        cat("\n\nNone\n");
        temp = ur.df(tdata, lags=lastlag, type="none", selectlags=slctlags);
        tstat = round(temp@teststat, 2);
        print(t(temp@cval));
        cat("-----------\nTest: ", tstat[1], "\n");
        if (tstat[1] < temp@cval[2]){cat("        *")}
        if (tstat[1] < temp@cval[1]){cat("*")}
        cat("\n")
      }
    }
  }
}





