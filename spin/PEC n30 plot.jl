@time begin
using Plots
pyplot()
default(
    size = (600,400),
    dpi = 100,
    fontfamily = "serif",
    titlefontsize = 16,
    guidefontsize = 16,
    tickfontsize = 16,
    legendfontsize = 12)
colors = get_color_palette(:auto, plot_color(:white), 7)
end

@time begin
Path = "/home/frederic/Documents/julia/";
n = 30 # energy level
ρ = readcsv("$Path/pec/$(n)_grid.csv") 
ylim = (-0.8,0.1) 

for Ω0 in [2., 1., 0.]
if Ω0 == 0.
    PES = readcsv("$Path/pec/$(n)_Ω=0.csv")
    plot(ρ/n^2, PES, xaxis = "R/n²", yaxis = ("E⋅n³", ylim), label = "", title = "Ω=$Ω0")
    savefig("$Path/pec/$(n)_Ω=$Ω0.pdf")
elseif Ω0 == 1.
    PES = readcsv("$Path/pec/$(n)_Ω=1.csv")
    plot(ρ/n^2, PES, xaxis = "R/n²", yaxis = ("E⋅n³", ylim), label = "", title = "Ω=$Ω0")
    savefig("$Path/pec/$(n)_Ω=$Ω0.pdf")
elseif Ω0 == 2.
    PES = readcsv("$Path/pec/$(n)_Ω=2.csv")
    plot(ρ/n^2, PES, xaxis = "R/n²", yaxis = ("E⋅n³", ylim), label = "", title = "Ω=$Ω0")
    savefig("$Path/pec/$(n)_Ω=$Ω0.pdf")
end
end
end

