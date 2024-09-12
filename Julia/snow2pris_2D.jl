using Revise, Parameters, Setfield, Plots, LinearAlgebra
using BifurcationKit
using DelimitedFiles
using JLD2
using CSV, Tables
const BK = BifurcationKit


Base.run(`clear`)
cd(@__DIR__)

recordCO(x, p) = (X = x[1], Y = x[2], Z = x[3])



function kindergarden(u, par)
  @unpack tauC, tauD, r,s,q,p = par
  x, y, z = u
  out = similar(u)
  out[1] = (y * (1 - x)) / tauC - (z * x) / tauD
  out[2] = y * ((tauC - 1) / tauC - y / tauC - z / tauD) + x * (r * x + s * (1 - x))
  out[3] = z * ((tauD - 1) / tauD - y / tauC - z / tauD) + (1 - x) * (q * x + p * (1 - x))
  out
end

ppp = 1.04


pard6 = (tauC =1.5, tauD=2.0, r = 3.05, s = 1.1, q =5.00, p = ppp)



# initial condition
z0 =  [1.1, 0.1, 0.7]


record_from_solution_snow(x, par) = (X = x[1], Y = x[2], Z = x[3], tauC = par[1])



# newton options
opt_newton = NewtonPar(tol = 1e-9, max_iterations = 50)




prob = BifurcationProblem(kindergarden, z0, setproperties(pard6 ; tauC =0.5),
  (@optic _.tauC);
  record_from_solution = record_from_solution_snow)

opts_br = ContinuationPar(dsmin = 0.001, dsmax = 0.05, ds = 0.01, 
# parameter interval
p_max =  10.0, p_min = 0.01, 
nev = 3, 
newton_options = opt_newton, 
max_steps = 100000, 
n_inversion = 4,
detect_bifurcation = 3)

# compute the branch of solutions
br = continuation(prob, PALC(), opts_br;
normC = norminf,
bothside = true)

print(br)

println(get_normal_form(br,2))




br1 = continuation(br, 2, opts_br;
normC = norminf,
bothside = true)
println(br1)


println(get_normal_form(br1,2))
println(get_normal_form(br1,3))
println(get_normal_form(br1,4))

br2 = continuation(br1, 2, opts_br;
normC = norminf,
bothside = true)
#println(br2)


sn_codim2 = continuation(br1, 4, (@optic _.p),
	ContinuationPar(opts_br, p_max = 1.5, p_min = 0.0,
		dsmin=1e-9, ds = 0.00001, dsmax = 0.0001) ;
	normC = norminf,
	# detection of codim 2 bifurcations with bisection
	detect_codim2_bifurcation = 2,
	update_minaug_every_step=1,
	start_with_eigen = true,
  	bothside = true,
	# we save the different components for plotting
	record_from_solution = record_from_solution_snow)



tr_codim2 = continuation(br1, 3, (@optic _.p),
	ContinuationPar(opts_br, p_max = 1.5, p_min = 0.0,
	dsmin=1e-9, ds = 0.00001, dsmax = 0.0001) ;
	normC = norminf,
	# detection of codim 2 bifurcations with bisection
	detect_codim2_bifurcation = 2,
	update_minaug_every_step=1,
	start_with_eigen = true,
	bothside = true,
	# we save the different components for plotting
	record_from_solution = record_from_solution_snow)







tr2_codim2 = continuation(br1, 2, (@optic _.p),
	ContinuationPar(opts_br, p_max = 1.5, p_min = 0.0,
	dsmin=1e-9, ds = 0.00001, dsmax = 0.0001) ;
	normC = norminf,
	# detection of codim 2 bifurcations with bisection
	detect_codim2_bifurcation = 2,
	update_minaug_every_step=1,
	start_with_eigen = true,
  	bothside = true,
	# we save the different components for plotting
	record_from_solution = record_from_solution_snow)

CSV.write("./csv/Snow_2_pris_p_tauC_fold_cont_f.csv",  Tables.table([sn_codim2.tauC sn_codim2.p]), writeheader=false)
CSV.write("./csv/Snow_2_pris_p_tauC_fold_cont_tr.csv",  Tables.table([tr_codim2.tauC tr_codim2.p]), writeheader=false)
CSV.write("./csv/Snow_2_pris_p_tauC_fold_cont_tr2.csv",  Tables.table([tr2_codim2.tauC tr2_codim2.p]), writeheader=false)


println("")
println(sn_codim2)
println("")


println(tr_codim2)
println("")

println(tr2_codim2)
println("")

println(get_normal_form(tr2_codim2,9))


