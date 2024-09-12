using Revise, Parameters, Setfield, Plots, LinearAlgebra
using CSV, Tables
using BifurcationKit
using DelimitedFiles
using JLD2
const BK = BifurcationKit


Base.run(`clear`)
cd(@__DIR__)


function kindergarden(u, par)
  @unpack tauC, tauD, r,s,q,p = par
  x, y, z = u
  out = similar(u)
  out[1] = (y * (1 - x)) / tauC - (z * x) / tauD
  out[2] = y * ((tauC - 1) / tauC - y / tauC - z / tauD) + x * (r * x + s * (1 - x))
  out[3] = z * ((tauD - 1) / tauD - y / tauC - z / tauD) + (1 - x) * (q * x + p * (1 - x))
  out
end

ppp = 1.19206938

tauD = 2.0
pard6 = (tauC =1.5, tauD=tauD, r = 3.05, s = 1.1, q =5.00, p = ppp)



# initial condition
z0 =  [1.1, 0.1, 0.7]


record_from_solution_snow(x, par) = (X = x[1], Y = x[2], Z = x[3], tauC = par[1])



# newton options
opt_newton = NewtonPar(tol = 1e-8, max_iterations = 100)

z0 =  [1.1, 0.1, 0.7]
prob = BifurcationProblem(kindergarden, z0, setproperties(pard6 ; tauC =2.0),
  (@optic _.tauC);
  record_from_solution = record_from_solution_snow)

opts_br = ContinuationPar(dsmin = 0.00001, dsmax = 0.001, ds = 0.0001, 
# parameter interval
p_max =  5.0, p_min = 0.01, 
nev = 3, 
newton_options = opt_newton, 
max_steps = 100000, 
n_inversion = 2000,
detect_bifurcation = 3)

# compute the branch of solutions
br = continuation(prob, PALC(), opts_br;
normC = norminf,
bothside = true,)

br
println(br)


br1 = continuation(br, 2, opts_br;
normC = norminf,
bothside = true)


println(br1)

br2 = continuation(br1, 2, opts_br;
normC = norminf,
bothside = true)


println(br2)


e1 = Vector{Float64}()
e2 = Vector{Float64}()
e3 = Vector{Float64}()
lambda_bar = Vector{Float64}()

for i = 1:length(br.eig)
  append!(e1,real.(br.eig[i].eigenvals[1]))
  append!(e2,real.(br.eig[i].eigenvals[2]))
  append!(e3,real.(br.eig[i].eigenvals[3]))
  append!(lambda_bar,br.Y[i]/br.tauC[i]+br.Z[i]/tauD)
end


CSV.write("./csv/snow2prisPf2Export_br.csv",  Tables.table([br.tauC br.X br.Y br.Z e1 e2 e3 lambda_bar]), writeheader=false)



e1 = Vector{Float64}()
e2 = Vector{Float64}()
e3 = Vector{Float64}()
lambda_bar = Vector{Float64}()

for i = 1:length(br1.eig)
  append!(e1,real.(br1.eig[i].eigenvals[1]))
  append!(e2,real.(br1.eig[i].eigenvals[2]))
  append!(e3,real.(br1.eig[i].eigenvals[3]))
  append!(lambda_bar,br1.Y[i]/br1.tauC[i]+br1.Z[i]/tauD)
end

CSV.write("./csv/snow2prisPf2Export_br1.csv",  Tables.table([br1.tauC br1.X br1.Y br1.Z e1 e2 e3 lambda_bar]), writeheader=false)


#reset inital conditoin for equilibria at x=0
z0 =  [0.1, 0.1, 0.7]
prob = BifurcationProblem(kindergarden, z0, setproperties(pard6 ; tauC =2.0),
  (@optic _.tauC);
  record_from_solution = record_from_solution_snow)

br2 = continuation(prob, PALC(), opts_br;
normC = norminf,
bothside = true,)


e1 = Vector{Float64}()
e2 = Vector{Float64}()
e3 = Vector{Float64}()
lambda_bar = Vector{Float64}()

for i = 1:length(br2.eig)
  append!(e1,real.(br2.eig[i].eigenvals[1]))
  append!(e2,real.(br2.eig[i].eigenvals[2]))
  append!(e3,real.(br2.eig[i].eigenvals[3]))
  append!(lambda_bar,br2.Y[i]/br2.tauC[i]+br2.Z[i]/tauD)
end


CSV.write("./csv/snow2prisPf2Export_br2.csv",  Tables.table([br2.tauC br2.X br2.Y br2.Z e1 e2 e3 lambda_bar]), writeheader=false)


