using QuantumCumulants
using SystemLevelHamiltonian

filter_cav = cavity(:filter)
sqz_cav = squeezing_cavity(:sqz)
ifo = radiation_pressure_cavity(:ifo)

sys = concatenation(:feqdepsqz,[sqz_cav,filter_cav,ifo])
sys = feedbackreduce(sys,:Output_sqz,:Input_filter)
sys = feedbackreduce(sys,:Output_filter,:Input_ifo)
