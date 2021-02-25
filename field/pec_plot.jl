@time begin
println("prepare packages...")
using Plots
pyplot()
default(
    size = (600,600),
    dpi = 100,
    fontfamily = "sans-serif",
    titlefontsize = 16,
    guidefontsize = 16,
    tickfontsize = 16,
    legendfontsize = 12)
colors = get_color_palette(:auto, plot_color(:white), 7)
Path = "/home/frederic/Documents/julia/";
end

@time begin
println("prepare plots...")
n = 30
F = 1000
ρ = readcsv("$Path/pec/data/pec_n=$(n)_grid.csv")

for state in ["a","t","s","p","d"]
println("$state state")
if state == "a"
	ylim = (-0.8,0.1) # full
elseif state == "t"
	ylim = (-0.05,0.01) # trilobite state
elseif state == "s"
	ylim = (-0.133,-0.132) # s state
elseif state == "p"
	ylim = (-0.6727,-0.6717) # p state
elseif state == "d"
	ylim = (-0.3475,-0.3455) # d state
else
	ylim = (-0.1,0.1) # hydrogenic
end
for ω in linspace(0.,π,3)
	println("θ=$(ω/π)")
	EE = readcsv("$Path/pec/data/pec_n=$(n)_E=$(F)_θ=$(ω/π).csv")
	plot(ρ/n^2, EE*n^3, xaxis = ("R/n²"), yaxis = ("n³E", ylim), title = "E=$F V/m, θ=$(ω/π)", label = "")
	savefig("$Path/pec/plots/$state/pec_n=$(n)_E=$(F)_θ=$(ω/π).pdf")
end
end
end

