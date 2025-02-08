HestonCallClosedForm <-
    function(lambda, vbar, eta, rho, v0, r, tau, S0, K) {
	PIntegrand <- function(u, lambda, vbar, eta, rho, v0, r, tau, S0, K, j) {
            F <- S0*exp(r*tau)
            x <- log(F/K)
            a <- lambda * vbar
            
            if (j == 1) {
                b <- lambda - rho* eta
                alpha <- - u^2/2 - u/2 * 1i + 1i * u
                beta <- lambda - rho * eta - rho * eta * 1i * u
            } else {
                b <- lambda
                alpha <- - u^2/2 - u/2 * 1i
                beta <- lambda - rho * eta * 1i * u
            }
            
            gamma <- eta^2/2
            d <- sqrt(beta^2 - 4*alpha*gamma)
            rplus <- (beta + d)/(2*gamma)
            rminus <- (beta - d)/(2*gamma)
            g <- rminus / rplus
            
            D <- rminus * (1 - exp(-d*tau))/(1-g*exp(-d*tau))
            C <- lambda * (rminus * tau - 2/(eta^2) * log( (1-g*exp(-d*tau))/(1-g) ) )
            
            top <- exp(C*vbar + D*v0 + 1i*u*x)
            bottom <- (1i * u)
            Re(top/bottom)
	}
	
	P <- function(lambda, vbar, eta, rho, v0, r, tau, S0, K, j) {
            value <- integrate(PIntegrand, lower = 0, upper = Inf,
                               lambda, vbar, eta, rho, v0, r, tau,
                               S0, K, j, subdivisions=1000)$value
            0.5 + 1/pi * value
	}

        A <- S0*P(lambda, vbar, eta, rho, v0, r, tau, S0, K, 1)
        B <- K*exp(-r*tau)*P(lambda, vbar, eta, rho, v0, r, tau, S0, K, 0)
        A-B
    }