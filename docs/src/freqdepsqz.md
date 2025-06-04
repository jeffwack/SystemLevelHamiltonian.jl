# Calculating the Hamiltonian for frequency-dependent squeezing

```@example freqdepsqz
using SystemLevelHamiltonian

filter_cav = cavity(:filter)
sqz_cav = squeezing_cavity(:sqz)
ifo = radiation_pressure_cavity(:ifo)

sys = concatenate(:feqdepsqz,[sqz_cav,filter_cav,ifo])
sys = feedbackreduce(sys,:Out_sqz,:In_filter)
sys = feedbackreduce(sys,:Out_filter,:In_ifo)
```

```@example freqdepsqz
sys.H
```

```@example freqdepsqz
operators(sys)
```

```@example freqdepsqz
parameters(sys)
```