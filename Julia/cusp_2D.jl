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



ppp = 0.01#1.19512195#1.0691064036398454#1.19

tauD = 5.0
pard6 = (tauC =1.5, tauD=tauD, r = 2.0, s = 0.5, q =2.95, p = ppp)


z0 =  [1.0, 1.0, 2.0]


# newton options
opt_newton = NewtonPar(tol = 1e-7, max_iterations = 1000)


prob = BifurcationProblem(kindergarden, z0, setproperties(pard6 ; tauC =1.5),
  (@optic _.tauC);
  record_from_solution = record_from_solution_snow)

opts_br = ContinuationPar(dsmin = 0.000001, dsmax = 0.0001, ds = 0.00001, 
# parameter interval
p_max =  10.0, p_min = 0.01, 
nev = 3, 
newton_options = opt_newton, 
max_steps = 100000, 
n_inversion = 2000,
detect_bifurcation = 3)



# compute the branch of solutions
br = continuation(prob, PALC(), opts_br;
normC = norminf,
bothside = true,)

br1 = continuation(br, 2, opts_br;
normC = norminf,
bothside = true)


println(br1)



sn_codim2 = continuation(br1, 2, (@optic _.p),
	ContinuationPar(opts_br, p_max = 0.4, p_min = 0.0,
	dsmin=1e-9, ds = 0.00001, dsmax = 0.0001) ;
	normC = norminf,
	# detection of codim 2 bifurcations with bisection
	detect_codim2_bifurcation = 2,
	update_minaug_every_step=1,
	start_with_eigen = true,
	bothside = true,
	# we save the different components for plotting
	record_from_solution = record_from_solution_snow,
  event = BK.ContinuousEvent(2,
		(iter, state) -> (getx(state)[1], getx(state)[1]-1)),)


println(sn_codim2)




CSV.write("./csv/cusp_sn2.csv",  Tables.table([sn_codim2.tauC sn_codim2.p]), writeheader=false)



tr_codim2 = continuation(br1, 3, (@optic _.p),
	ContinuationPar(opts_br, p_max = 0.4, p_min = 0.0,
	dsmin=1e-9, ds = 0.00001, dsmax = 0.0001) ;
	normC = norminf,
	# detection of codim 2 bifurcations with bisection
	detect_codim2_bifurcation = 2,
	update_minaug_every_step=1,
	start_with_eigen = true,
	bothside = true,
	# we save the different components for plotting
	record_from_solution = record_from_solution_snow,
  event = BK.ContinuousEvent(2,
		(iter, state) -> (getx(state)[1], getx(state)[1]-1)),)


println(tr_codim2)
CSV.write("./csv/cusp_tr1.csv",  Tables.table([tr_codim2.tauC tr_codim2.p]), writeheader=false)





tauD = 5.0
pard6 = (tauC =1.5, tauD=tauD, r = 2.0, s = 0.5, q =2.95, p = ppp)


z0 =  [0.0, 0.0, 4.0]

prob = BifurcationProblem(kindergarden, z0, setproperties(pard6 ; tauC =1.5),
  (@optic _.tauC);
  record_from_solution = record_from_solution_snow)

brlow = continuation(prob, PALC(), opts_br;
normC = norminf,
bothside = true,)


br1low = continuation(brlow, 2, opts_br;
normC = norminf,
bothside = true)



println(br1)




sn_codim2low = continuation(br1low, 4, (@optic _.p),
	ContinuationPar(opts_br, p_max = 0.3, p_min = 0.0,
	dsmin=1e-9, ds = 0.00001, dsmax = 0.0001) ;
	normC = norminf,
	# detection of codim 2 bifurcations with bisection
	detect_codim2_bifurcation = 2,
	update_minaug_every_step=1,
	start_with_eigen = true,
	bothside = true,
	# we save the different components for plotting
	record_from_solution = record_from_solution_snow,
  event = BK.ContinuousEvent(2,
		(iter, state) -> (getx(state)[1], getx(state)[1]-1)),)


println(sn_codim2)



cuspNF = get_normal_form(sn_codim2, 2)
println("normalform")
println(cuspNF)


CSV.write("./csv/cusp_sn2.csv",  Tables.table([sn_codim2low.tauC sn_codim2low.p]), writeheader=false)

tr_codim2low = continuation(br1low, 2, (@optic _.p),
	ContinuationPar(opts_br, p_max = 0.3, p_min = 0.0,
	dsmin=1e-9, ds = 0.00001, dsmax = 0.0001) ;
	normC = norminf,
	# detection of codim 2 bifurcations with bisection
	detect_codim2_bifurcation = 2,
	update_minaug_every_step=1,
	start_with_eigen = true,
	bothside = true,
	# we save the different components for plotting
	record_from_solution = record_from_solution_snow,
  event = BK.ContinuousEvent(2,
		(iter, state) -> (getx(state)[1], getx(state)[1]-1)),)

CSV.write("./csv/cusp_tr2.csv",  Tables.table([tr_codim2low.tauC tr_codim2low.p]), writeheader=false)


