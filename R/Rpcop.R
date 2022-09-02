#' @title pcop
#' @description Computes a principal curve as defined in Delicado (2001)
#' <doi:10.1007/s001800300145>
#' @param x 		a matrix of n points in dimension p
#' @param Ch 	The smoothing parameter h is C_H times the value given by the normal reference rule.
#' 							Default value 1.5. Constraints 0.5 <= C_H <=  1.5
#' @param Cd 	The distance between two consecutive principal oriented points in a PCOP is about C_D times
#'  						the value of the smoothing parameter h.
#'  						Default value 0.3. Constraints 0.25 <= C_D <= .5
#' @param plot.true if TRUE, the function produce a plot
#' @param ... Additional parameters passed to function "lines"
#' @return A list with two data frames. One contains a list with the following names:
#' 'param': Value of the parameter t such the the principal oriented point is PCOP(t).
#' 'dens': Density estimation for the random variable induced over the PCOP at t.
#' 'span': proportion of original data involved in the determination of the principal oriented point.
#' 'orth.var': Variance over the hyperplane orthogonal to the PCOP at the principal oriented point.
#' 'pop': a p-dimensional array. The p coordinates of the principal oriented point.
#' 'pr.dir': a p-dimensional array. The p coordinates of the principal direction for the principal oriented point.
#' For the second, look at the package princurve.
#' @examples
#' x <- runif(100,-1,1)
#' x <- cbind(x, x ^ 2 + rnorm(100, sd = 0.1))
#' pcop(x, plot.true=TRUE, lwd=4, col=2)
#' @importFrom princurve "project_to_curve"
#' @export
pcop <- function(x,Ch=1.5,Cd=.3,plot.true=FALSE,...){
	proyec <- as.data.frame(pcop_backend(x, Cd, Ch));

	pr.curve <- proyec[,-1]

		p <- dim(x)[2]

		names(pr.curve) <- c("param", "dens", "span", "orth.var", 
				paste("pop", 1:p, sep = ""), 
				paste("pr.dir", 1:p, sep = ""))
		# pr.curve: principal curve of oriented points 
		#           with original format (format 1)

		# pc: principal curve of oriented points 
		#     with format as in HS's "princurv" library (format 2)
		s<-as.matrix(pr.curve[,5:(5+p-1)])
		pc <- project_to_curve(x,s)

		if (plot.true) {
			adjust.range <- function(x,cte=1){
				rgx<-range(x)
					extra <- diff(rgx)*(cte-1)/2
					argx1<-rgx[1]-extra
					argx2<-rgx[2]+extra
					return(c(argx1,argx2))
			}
			plot(x[, 1:2], xlim = adjust.range(x[, 1], 1.1), ylim = adjust.range(x[,2], 1.1), col=8)
				lines(pc$s[pc$ord, 1:2],...)
		}


	return(list(pcop.f1=pr.curve,pcop.f2=pc))
}


