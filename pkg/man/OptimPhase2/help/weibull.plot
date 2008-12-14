weibull.plot           package:OptimPhase2           R Documentation

_P_l_o_t _W_e_i_b_u_l_l _S_u_r_v_i_v_a_l _C_u_r_v_e_s

_D_e_s_c_r_i_p_t_i_o_n:

     Plot Weibull survival curves with differences at a target time
     highlighted. Designed for use with the 'param' values input to
     function 'OptimDes'.

_U_s_a_g_e:

     weibull.plot(param, x, l.type = 1:3, l.col = c("blue", "red"), ...)

_A_r_g_u_m_e_n_t_s:

   param: A four-element vector containing: the shape parameter of the
          Weibull distribution under the null hypothesis, the scale
          parameter of the Weibull distribution under the null
          hypothesis, the shape parameter of the Weibull distribution
          under the alternative hypothesis, the scale parameter of the
          Weibull distribution under the alternative hypothesis.

       x: Survival time of interest (e.g., 1 year).

  l.type: Line types for the plot. Default is '1-3'.

   l.col: Line colors for the plot. Default is '"blue"' for the null
          survival curve, '"red"' for the alternative survival curve.

     ...: Further graphical arguments, see 'plot.default'.

_A_u_t_h_o_r(_s):

     Bo Huang <huang@stat.wisc.edu> and Neal Thomas

_R_e_f_e_r_e_n_c_e_s:

     Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) Continuous
     Univariate Distributions, volume 1, chapter 21. Wiley, New York.

_S_e_e _A_l_s_o:

     'dweibull', 'OptimDes', 'weibPmatch'

_E_x_a_m_p_l_e_s:

     param <- c(1, 1.09, 2, 1.40)
     x <- 1

     weibull.plot(param,x)

