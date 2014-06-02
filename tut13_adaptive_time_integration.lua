--------------------------------------------------------------------------------
--	tut13_adaptive_time_integration.lua
--
--	This tutorial is used to show tha ability to compute non-linear problems.
--	The problem itself will be a simple convection-diffusion equation that is
--	non-linear since the coefficients are non-linear.
--------------------------------------------------------------------------------

-- We will include a script file which defines some methods often used.
-- Loaded methods are all found in the library util.
ug_load_script("ug_util.lua")

-- To keep the script flexible we will now define some variables which have
-- a default value that can be overwritten by command line arguments.

-- Depending on the dimension we will choose our domain object
-- (either 1d, 2d or 3d) and associated discretization objects. Note that
-- we're using some methods defined in "ug_util.lua" here.
dim = util.GetParamNumber("-dim", 2) -- default dimension is 2.

-- Since ug supports a bunch of different dimensions and algebra modules 
-- we will choose a combination here. This should always be the first thing 
-- you do in an ug-script. The cpu-algebra is fine for now.
InitUG(dim, AlgebraType("CPU", 1))

-- We also need a filename for the grid that shall be loaded.
if 		dim == 1 then gridName = util.GetParam("-grid", "unit_square_01/unit_line_01_edge_2.ugx")
elseif 	dim == 2 then gridName = util.GetParam("-grid", "unit_square/unit_square_quads_8x8.ugx")
elseif 	dim == 3 then gridName = util.GetParam("-grid", "unit_square/unit_cube_hex.ugx")
else print("Choosen Dimension not supported. Exiting."); exit(); end

-- Since we want to save the domains grid to a file, we also need an output file.
-- Note that we only use a prefix here, since we want to attach the process number
-- and file format ourselfs
outFileNamePrefix = util.GetParam("-o", "distributed_domain_")


-- We read in the number of time steps to be performed
numTimeSteps = util.GetParamNumber("-numTimeSteps", 1000) -- default 

charLenth = 2;
lambda = util.GetParamNumber("-lambda", 1.0)   -- diffusion coefficient
mu = util.GetParamNumber("-mu", 1.0)           -- parameter for solution
nu = util.GetParamNumber("-nu", 1.0)		   -- parameter for solution

charTime = charLenth*charLenth/lambda
maxTime = charTime*10;

dt = 1e-2
dtMin = 1e-8
dtMax= 0.1*maxTime

doControl = true
doExtrapolation = false
doVTKOutput = false;

-- We additionally use parameters which allow to specify the number of
-- pre- and total-refinement steps (wrt domain distribution).
numPreRefs = util.GetParamNumber("-numPreRefs", 1)
numTotalRefs = util.GetParamNumber("-numTotalRefs", 4)

baseLevel= util.GetParamNumber("-baseLevel", 0) -- base level for gmg solver

-- Calculate the number of post-refs and make sure that the result makes sense.
numPostRefs = numTotalRefs - numPreRefs
if numPostRefs < 0 then
	print("WARNING:\tnumPreRefs exceeds the number of total refinements.")
	print("\t\t\tNo refinement will be preformed after distribution.")
	numPostRefs = 0
end

util.CheckAndPrintHelp("Example 13: Adaptive time stepping & extrapolation");



-- Now its time to create the domain object. We will use an util method here,
-- which automatically creates the right domain for the given dimension.
dom = Domain()

-- Now that we have a domain, we can load a grid into it. We check the return
-- value whether the loading was successful or not.
-- Note that we use the method utilLoadDomain instead of LoadDomain. utilLoadDomain
-- has the benefit that grids are automatically searched in the data/grids folder if
-- they were not found at the default locations (execution-path or a path specified
-- in your environments path-variable).
LoadDomain(dom, gridName)

-- Now that we're here, the domain was successfully loaded
print("Loaded domain from " .. gridName)

-- We will create a refiner now. This is a tool associated with a domain.
-- UG defines factory methods for refiners, which automatically choose
-- the right refiner for the given context, i.e. different refiners are
-- created depending on whether we are in a parallel or in a serial environment.
-- Note that another factory method is HangingNodeDomainRefiner, which is
-- subject to a later tutorial.
refiner = GlobalDomainRefiner(dom)

-- perform pre-refinement
for i = 1, numPreRefs do
	refiner:refine()
end

-- Distribute the domain to all involved processes
if util.DistributeDomain(dom) == false then
	print("Error while Distributing Domain. Aborting.")
	exit()
end


-- perform post-refinement
for i = 1, numPostRefs do
	refiner:refine()
end

-- Lets save the domain on each process
outFileName = outFileNamePrefix .. ProcRank() .. ".ugx"
SaveDomain(dom, outFileName)

-- Everything seems to went fine.
print("Saved domain to " .. outFileName)

-- Check the subset handler if all subsets are given
sh = dom:subset_handler()
if sh:get_subset_index("Inner") == -1 then
	print("Domain does not contain subset 'Inner'. Aborting.")
	exit()
end

if sh:get_subset_index("Boundary") == -1 then
	print("Domain does not contain subset 'Boundary'. Aborting.")
	exit()
end


-- We create an Approximation Space
-- which describes the unknowns in the equation. In our case we will
-- only require one unknown for the concentration ("c").
-- Note that the Approximation Space is build on the domain created above.
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom) -- creates new object
approxSpace:add_fct("c", "Lagrange", 1)          -- adds one function

----------------------------------------------------
----------------------------------------------------
-- User function and User Data
----------------------------------------------------
----------------------------------------------------

-- Next we need a lot of User Functions and User Data to specify the 
-- convection-diffusion problem's coefficient. The data can be
-- set from the lua script.


--[[ 
Reference solution:
 sin(nu*x)'' = - nu*nu*sin(nu*x)
 \dot u = -alpha u => u = u0*exp(-alpha*t) 
--]]
function ourSolution2d(x,y,t)
	mupi = mu*math.pi
	nupi = nu*math.pi
	return math.exp(-lambda*(mupi*mupi+nupi*nupi)*t)*math.sin(mupi*x)*math.sin(nupi*y)
end


----------------------------------------------------
-- Right-Hand Side
----------------------------------------------------

-- The same for the right hand side ..
function ourRhs1d(x, t)
	return 0.0
end

function ourRhs2d(x, y, t)
	return 0.0
end

function ourRhs3d(x, y, z, t)
	return 0.0
end

rhsCallback = LuaUserNumber("ourRhs" .. dim .."d")


----------------------------------------------------
-- Convection - Diffusion Element Discretization
----------------------------------------------------
elemDisc = ConvectionDiffusion("c", "Inner", "fv1")
elemDisc:set_diffusion(lambda)				-- set diffusion coefficient
elemDisc:set_mass_scale(1.0)
elemDisc:set_source(rhsCallback)		    -- set the right hand side

----------------------------------------------------
-- Dirichlet Boundary
----------------------------------------------------

-- We create some boundary conditions as shown in the preceeding tutorials ...
function ourDirichletBnd1d(x, t)
	return true, 0.0
end

function ourDirichletBnd2d(x, y, t)
	return true, ourSolution2d(x,y,t)
end

function ourDirichletBnd3d(x, y, z, t)
	return true, 0.0
end

-- lets setup the dirichlet values as explained in the previous tutorials
dirichletBnd = DirichletBoundary()
dirichletBnd:add("ourDirichletBnd" .. dim .. "d", "c", "Boundary")

----------------------------------------------------
-- Adding all discretizations
----------------------------------------------------

-- Finally we create the discretization object which combines all the
-- separate discretizations into one domain discretization.
domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
domainDisc:add(dirichletBnd)

-- Now we create a time discretization. We use the theta-scheme. The time 
-- stepping scheme gets passed the domain discretization and will assemble
-- mass- and stiffness parts using the domain disc, then adding them in a
-- weighted way.
timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(1.0) -- 1.0 is implicit euler

-- Now we create an operator from the time discretization. We use the 
-- (nonlinear)-Operator interface.
op = AssembledOperator()
op:set_discretization(timeDisc)
op:init()

----------------------------------------------------
----------------------------------------------------
-- Solver setup
----------------------------------------------------
----------------------------------------------------

-- we need a linear solver that solves the linearized problem inside of the
-- newton solver iteration. We use a geometric multi-grid method with
-- Jacobi smoothing and an LU base-solver.
baseSolver = LU()

gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_base_solver(baseSolver)
gmg:set_base_level(baseLevel)
gmg:set_gathered_base_solver_if_ambiguous(true)
gmg:set_smoother(Jacobi(0.66))
gmg:set_cycle_type(1)
gmg:set_num_presmooth(2)
gmg:set_num_postsmooth(2)

linSolver = LinearSolver()
linSolver:set_preconditioner(gmg)
linSolver:set_convergence_check(ConvCheck(100, 1e-12, 1e-12))
linSolver:set_compute_fresh_defect_when_finished(true)


-- Next we need a convergence check, that computes the defect within each 
-- newton step and stops the iteration when a specified creterion is fullfilled.
-- For our purpose the ConvCheck is sufficient. Please note,
-- that this class derives from a general IConvergenceCheck-Interface and
-- also more specialized or self-coded convergence checks could be used.
newtonConvCheck = ConvCheck()
newtonConvCheck:set_maximum_steps(10)
newtonConvCheck:set_minimum_defect(5e-8)
newtonConvCheck:set_reduction(1e-10)
newtonConvCheck:set_verbose(true)

-- Within each newton step a line search can be applied. In order to do so an
-- implementation of the ILineSearch-Interface can be passed to the newton
-- solver. Here again we use the standard implementation.
newtonLineSearch = StandardLineSearch()

-- Now we can set up the newton solver. We set the linear solver created above
-- as solver for the linearized problem and set the convergence check. If you
-- want to you can also set the line search.
newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(linSolver)
newtonSolver:set_convergence_check(newtonConvCheck)
--newtonSolver:set_line_search(newtonLineSearch)

if (doControl) then
-- create Newton Solver2 (control for time step)
	newtonSolver2 = NewtonSolver()
	newtonSolver2:set_linear_solver(linSolver)
	newtonSolver2:set_convergence_check(newtonConvCheck)
	newtonSolver2:set_line_search(newtonLineSearch)
end


-- Finally we set the non-linear operator created above and initiallize the
-- Newton solver for this operator.
newtonSolver:init(op)

-- since solver configurations can be quite complex, we print the configuration:
print("NewtonSolver configuration:")
print(newtonSolver:config_string())
----------------------------------------------------
----------------------------------------------------
-- Time loop
-- (this is where the fun starts!)
----------------------------------------------------
----------------------------------------------------

-- Now we are at the point to initalize the time stepping and start the time
-- loop

-- We create a grid function on the surface of our MultiGrid hierarchy.
u = GridFunction(approxSpace)

-- Lets chose a fixed time step size and start at time point 0.0. Our first
-- step is step number 0.

time = 0.0
step = 0

-- Next we have to initialize the solution with the start configuration. We
-- interpolate the data given by a lua-callback

-- ... and wrap the lua-callback
LuaStartValue = LuaUserNumber("ourSolution"..dim.."d")
LuaSolutionValue = LuaUserNumber("ourSolution"..dim.."d")

-- Now interpolate the function
Interpolate(LuaStartValue, u, "c", time);

-- In order to plot our time steps, we need a VTK writer. For time dependent
-- problems we start a time series. This is necessary to group the time 
-- series at the end of the time loop. We write the start solution at the beginning.
out = VTKOutput()
out:print("StartValue", u, step, time)

-- Since we use a one-step scheme, we need one extra solution vector to store
-- the old time step. This is done using the clone method. All previous 
-- solutions (in this case only one) are stored in the "PreviousSolution" 
-- object, that behaves like a queue of fixed size. We push the start solution
-- as the first old time step to our queue
uOld = u:clone()
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)

exact = u:clone()

-- Additional setup for time step control
if (doControl) then
	-- discretization
	timeDisc2 = ThetaTimeStep(domainDisc, 1.0)
	
	old2 = uOld:clone()   -- temporary storage (for two step)
	u2  = uOld:clone()    -- second solution
	
	-- operator
	op2 = AssembledOperator()
	op2:set_discretization(timeDisc2)
	op2:init()
	newtonSolver2:init(op2)

	-- time series
	solTimeSeries2 = SolutionTimeSeries()
	solTimeSeries2:push(old2, time)
	if (doExtrapolation) then	
	-- Aitken-Neville-type time	extrapolation
	timex = AitkenNevilleTimex({1,2})
	end
end


local file = assert(io.open("tut13-errors.csv", "w"))
file:write("#plot \"tut13-errors.txt\" using 1:3, \"tut13-errors.txt\" using 1:4, \"tut13-errors.txt\" using 1:5\n")
file:write("#t\t|u|\te_coarse|\t|e_fine|\t|e_xtra|\test(c/f)\test(c/x)\ttimex(sub-diag)")

-- Run time loop
for step = 1, numTimeSteps do

	-- check, if finished
	if (time >= maxTime) then
		print ("Final time "..time.." reached.") 
		break
	end
	
	-- we choose the time step size to be constant here
	local tau = dt
	print("++++++ TIMESTEP " .. step .. "(" ..time .."+"..tau..") BEGIN ++++++")
	
	-- setup time Disc for old solutions and timestep
	timeDisc:prepare_step(solTimeSeries, tau)
	
	-- prepare newton solver
	if newtonSolver:prepare(u) == false then print ("Newton solver prepare failed at step "..step.."."); exit(); end 
	
	-- apply newton solver
	if newtonSolver:apply(u) == false then print ("Newton solver apply failed at step "..step.."."); exit(); end 

	-- set default strategy
	bAcceptStep = true;       -- accept solution 
	dtnew= dt;                -- next time step will be the same

	-- compute norms (we do not need old2 -> u^{1/2} value here any more!)	
	Interpolate(LuaSolutionValue, exact, "c", time+tau); l2normB = L2Norm(exact, "c", 2);
	l2norm = L2Norm(u, "c", 2); VecScaleAdd2(exact, 1.0, exact, -1.0, u); l2_coarse_err = L2Norm(exact, "c", 2); 

	-- time step control
	if (doControl) then 
		---------------------------
		-- first half step for u2
		---------------------------
		print("Control 1/2:");
		
		--time2 = time-tau/2;                        			-- intermediate time step 
		VecScaleAdd2(u2, 1.0-0.5*tau, uOld, 0.5*tau, u);        -- w/ linear interpolation  (first guess)
		timeDisc2:prepare_step(solTimeSeries2, 0.5*tau)
		if newtonSolver2:prepare(u2) == false then print ("Newton solver failed at step "..step.."+1/2."); exit(); end 
		if newtonSolver2:apply(u2) == false then print ("Newton solver failed at step "..step.."+1/2."); exit(); end 
		
		----------------------------
		-- second half step for u2
		----------------------------
		
		print("Control 2/2:");
		
		-- push back solution
		time2 = time + 0.5*tau
		tmp2 = solTimeSeries2:oldest()                      -- solution at time t 
		VecScaleAssign(tmp2, 1.0, u2)                       -- is removed and replaced by u2(t+tau/2)
		solTimeSeries2:push_discard_oldest(tmp2, time2)     -- set as new solution (now uOld is discarded)
		
		timeDisc2:prepare_step(solTimeSeries2, tau/2)
		if newtonSolver2:prepare(u2) == false then print ("Newton solver failed at step "..step.."+2/2."); exit(); end 
		if newtonSolver2:apply(u2) == false then print ("Newton solver failed at step "..step.."+2/2."); exit(); end 
		
		-- compute norms (we do not need old2 -> u^{1/2} value here any more!)
		l2norm = L2Norm(u2, "c", 2); 
		Interpolate(LuaSolutionValue, exact, "c", time+tau);
		VecScaleAdd2(exact, 1.0, exact, -1.0, u2); l2_fine_err = L2Norm(exact, "c", 2); 
		VecScaleAdd2(exact, 1.0, u, -1.0, u2); l2_fine_est = L2Norm(exact, "c", 2);
	
		if (doExtrapolation) then		
			-- extrapolation (more accurate)
			timex:set_solution(u, 0)
			timex:set_solution(u2, 1)
			timex:apply()
			eps = timex:get_error_estimate()
		
			-- compute norms (we do not need old2 -> u^{1/2} value here any more!)
			l2norm = L2Norm(u2, "c", 2);
			Interpolate(LuaSolutionValue, exact, "c", time+tau); 
			VecScaleAdd2(exact, 1.0, exact, -1.0, u2); l2_ex_err = L2Norm(exact, "c", 2); 
			VecScaleAdd2(exact, 1.0, u, -1.0, u2); l2_ex_est = L2Norm(exact, "c", 2);
		else
			-- no extrapolation (less accurate)
			eps = 2.0*l2_fine_est;
			l2_ex_err = "---";
			l2_ex_est = "---";
		end
		
		print("TIME_ERROR (t="..time+tau..", |u|="..l2normB.. ") :\t" .. l2_coarse_err .. "\t"..l2_fine_err .. "\t"..l2_ex_err .. "\|  ---\t"..l2_fine_est .. "\t"..l2_ex_est.. "\t"..eps)
		file:write(time+tau.."\t" ..l2normB.. "\t" .. l2_coarse_err .. "\t"..l2_fine_err .. "\t"..l2_ex_err .. "\t"..l2_fine_est .. "\t"..l2_ex_est.. "\t"..eps.."\n")
		file:flush()
		-------------------------
		-- Adaptive step control
		-------------------------
		tol = 1e-4;
		dtEst = math.pow(0.9*tol/eps, 0.5)*dt  -- (eps<=tol) implies (tol/eps >= 1) 
		print("dtEst= "..dtEst..", eps="..eps)
		
		-- determine potential new step size
		dtnew = math.min(dtEst, 1.5*dt, dtMax)

		if (eps <= tol) then 
			-- accept
			bAcceptStep = true;
		 	dtnew = math.min(dtnew, maxTime-time);
		 	print ("ACCEPTING solution, dt="..dtnew);
		 
		else
	    	-- discard
	    	bAcceptStep = false;
	    	print ("DISCARDING solution, dtnew="..dtnew);
	    	
	    	-- need resetting solTimeSeries2
	    	utmp = solTimeSeries2:oldest()
	   		VecScaleAssign(utmp, 1.0, uOld)
	   		solTimeSeries2:push_discard_oldest(utmp, time)
		end

	end
	
	-- next time step
	--dt = kappa *dt;
	dt = dtnew
	-- if (dt > dtMax) then dt = dtMax end;
	if (dt < dtMin) then exit() end;

	if (bAcceptStep) then
		-- compute the new (absolut) time
		time = solTimeSeries:time(0) + tau
		if (doVTKOutput) then
		-- we write the newly computed time step to our time series
		Interpolate(LuaSolutionValue, exact, "c", time+tau); 
		out:print("RealSolution", exact, step, time)
		out:print("CompSolution1", u, step, time)

		if (doControl) then out:print("CompSolution2", u2, step, time)  end
		
		end
		
		-- determine most accurate solution
		if (doControl) then unext = u2 
		else unext = u end;
			
		-- get oldest solution,
		-- copy values into oldest solution (reusing memory),
		-- push oldest solutions with new values to front, oldest sol pointer is poped from end
		local oldestSol = solTimeSeries:oldest()
		VecScaleAssign(oldestSol, 1.0, unext)
		solTimeSeries:push_discard_oldest(oldestSol, time)
		
		if (doControl) then
			-- do the same for second 
			oldestSol = solTimeSeries2:oldest()
			VecScaleAssign(oldestSol, 1.0, unext)
			solTimeSeries2:push_discard_oldest(oldestSol, time)
		end
		
	end
	print("++++++ TIMESTEP " .. step .. "(" ..time ..") END ++++++");
end

-- At the end of the time loop, we finish our time series. This produces a
-- grouping "Solution.pvd" file, that containes all time steps and can be
-- opened by a viewer like Paraview.
if (doVTKOutput) then
	out:write_time_pvd("CompSolution1", u)
	out:write_time_pvd("CompSolution2", u2)
end

file:close()  -- close file

print("")
print("done.")

