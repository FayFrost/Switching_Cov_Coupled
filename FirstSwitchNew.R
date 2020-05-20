FirstSwitchNew <- function(verb=TRUE, Kappa, swlamda, up_start, n_ou, n_bm) # seed=NULL
{
   # Handle randomization
   
  #if(!is.null(seed)){set.seed(seed)}
  
   set.seed(3809)
   # From main code
   
   pswt = 0
   npsw = 0
   po_tvec <- c()
   repeat{
      npsw=npsw+1
      po_t = rexp(1, Kappa)
      # cat("po_t ", po_t, "\n")
      po_tvec = c(po_tvec, po_t)
      pswt = pswt + po_t            #potential switching time 
      swlamda_now <- lookup_lamda(swlamda, pswt + up_start) ## to see which covariate parameter to choose
      #cat("swlambda_now", swlamda_now, "time", pswt + up_start,  "\n")
      P_ac_sw =( n_ou*swlamda_now[1] + n_bm*swlamda_now[2] )/Kappa # prob of actual switch
      #cat("P_ac_sw",P_ac_sw,"\n")
      ac_sw = sample( c(0,1), 1, FALSE, prob=c(1-P_ac_sw,P_ac_sw) ) # test for actual switch
      if(ac_sw==1) break else
         pon_probAll <- c(pon_probAll,1-P_ac_sw) #last term \Kappa 
   }
   # Results
   if (verb) cat("new po_t vector", po_tvec, npsw, "\n")
   return(list('pswt'= pswt,  'po_tvec'=po_tvec, 'pon_probAll'= pon_probAll))  #invisible(po_t)
}

