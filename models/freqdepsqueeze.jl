using QuantumCumulants
using SystemLevelHamiltonian

#########################################
#Define the filter cavity 
h_fabry = FockSpace(:fabry)
@qnumbers av::Destroy(h_fabry)

@rnumbers κv Δv 

filter_cavity = SLH(["FC"],
                [1],
                [sqrt(κv)*av],
                 Δv*adjoint(av)*av,
                 h_fabry,
                 [av])
#########################################

#########################################
#Define the squeezer
h_squeeze = FockSpace(:squeezer)
@qnumbers aq::Destroy(h_squeeze)

@rnumbers ϵ κq

squeezer = SLH(["SQ"],
            [1],
            [sqrt(κq)*aq],
            1im*ϵ*(adjoint(aq)^2- aq^2),
            h_squeeze,
            [aq])
###########################################

###########################################
#Define the radiation pressure noise (RPN) cavity
h_rpn = FockSpace(:rpn_optical) ⊗ FockSpace(:rpn_mirror)
ar = Destroy(h_rpn,:ar,1) 
b = Destroy(h_rpn,:b,2)

@rnumbers κr Δc Δm g Γ nbar

rpn_cavity = SLH(   ["LIGHT"],
                    [1],
                    [sqrt(κr)*ar],
                    Δc*ar'*ar+Δm*b'*b - g*ar'*ar*(b'+b),
                    h_rpn,
                    [ar,b])
#############################################

A = concatenation(squeezer,filter_cavity)
B = concatenation(A,rpn_cavity)
C = feedbackreduce(B,"SQ","FC")
D = feedbackreduce(C,"FC","LIGHT")
