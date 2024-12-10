mata:

			void evalnotch(todo,e, t0 , t1 , deltaT, cutoff , eresp,v, g, H)
			 {
				v = ( (1/(1+eresp/cutoff))*(1+(deltaT/cutoff)/(1-t0))-(1/(1+1/e))*(1/(1+(eresp/cutoff)))^(1+1/e)-(1/(1+e))*(1-(t1-t0)/(1-t0))^(1+e))^2
			}

			function notch(real scalar t0, real scalar t1, real scalar deltaT, real scalar cutoff, real scalar eresp, real scalar init)
				{
				S = optimize_init()
				optimize_init_which(S,  "min")
				optimize_init_conv_ptol(S, 1e-12)
				optimize_init_conv_vtol(S, 1e-12)
				optimize_init_evaluator(S, &evalnotch())
				optimize_init_evaluatortype(S, "d0")
				optimize_init_argument(S, 1, t0)
				optimize_init_argument(S, 2, t1)
				optimize_init_argument(S, 3, deltaT)
				optimize_init_argument(S, 4, cutoff)
				optimize_init_argument(S, 5, eresp)
				optimize_init_params(S, init)
				optimize_init_conv_maxiter(S,100)
				e=optimize(S)
				return(e,optimize_result_errorcode(S))
				}
				
			function eresp(real scalar B,real scalar tau,real matrix cf, real scalar bw)
				{
				integral=polyinteg(cf,1)
				integral[1]=-polyeval(integral,tau)-B*bw
				roots=polyroots(integral)
				realroots=Re(select(roots, Im(roots):==0))
				out=sort(select(realroots,realroots:>tau)',1)'
				if (cols(out)==0) {
					return(.)
				}
				else return(out[1]-tau)
				}
	
function varcorrect(real matrix X, real matrix fw, real matrix b)
				{
				B=length(fw)
				N=sum(fw)
				X=X,J(rows(X),1,1)
				e=-X * b'
				meat=0
				for (i=1; i<=B; i++) {
					e[i]=e[i]+N
					meat=meat:+fw[i]:*(X' * e * e' * X)
					e[i]=e[i]-N
				}
				return(invsym(quadcross(X,X):*N)*meat*invsym(quadcross(X,X):*N))
				}
				
