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



record_from_solution_snow(x, par) = (X = x[1], Y = x[2], Z = x[3], tauC = par[1])


tauD = 5.0

PP = 0.03#25
RR = 2.
SS = 0.5
QQ = 2.95
pard6 = (tauC =1.5, tauD=tauD, r = RR, s = SS, q =QQ, p = PP)

z0 =  [1.0, 1.0, 2.0]


# newton options
opt_newton = NewtonPar(tol = 1e-7, max_iterations = 2000)


prob = BifurcationProblem(kindergarden, z0, setproperties(pard6 ; tauC =1.5),
  (@optic _.tauC);
  record_from_solution = record_from_solution_snow)

opts_br = ContinuationPar(dsmin = 0.000001, dsmax = 0.001, ds = 0.0001, 
  # parameter interval
  p_max =  5.0, p_min = 0.01, 
  nev = 3, 
  newton_options = opt_newton, 
  max_steps = 50000, 
  n_inversion = 1000,
  detect_bifurcation = 3)


# compute the branch of solutions
br = continuation(prob, PALC(), opts_br;
normC = norminf,
bothside = true,)

println("branch 0")
println(br)

br1 = continuation(br, 2, opts_br;
normC = norminf,
bothside = true)

println("branch 1")
println(br1)

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

CSV.write("./csv/afterCuspUp_br_eig_lam.csv",  Tables.table([br.tauC br.X br.Y br.Z e1 e2 e3 lambda_bar]), writeheader=false)

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

CSV.write("./csv/afterCuspUp_br1_eig_lam.csv",  Tables.table([br1.tauC br1.X br1.Y br1.Z e1 e2 e3 lambda_bar]), writeheader=false)








pard6 = (tauC =1.5, tauD=tauD, r = RR, s = SS, q =QQ, p = PP)


z0 =  [0.0, 0.0, 4.0]

prob = BifurcationProblem(kindergarden, z0, setproperties(pard6 ; tauC =1.5),
  (@optic _.tauC);
  record_from_solution = record_from_solution_snow)

br = continuation(prob, PALC(), opts_br;
normC = norminf,
bothside = true,)

println(br)



br1 = continuation(br, 2, opts_br;
normC = norminf,
bothside = true)




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

CSV.write("./csv/afterCuspLow_br_eig_lam.csv",  Tables.table([br.tauC br.X br.Y br.Z e1 e2 e3 lambda_bar]), writeheader=false)

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

CSV.write("./csv/afterCuspLow_br1_eig_lam.csv",  Tables.table([br1.tauC br1.X br1.Y br1.Z e1 e2 e3 lambda_bar]), writeheader=false)

