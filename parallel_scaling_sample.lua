--------------------------------------------------------------------------------
--[[!
-- \author Andreas Vogel, Martin Rupp 
-- \brief This script is intended to show how to do parallel scaling stuff
-- based upon apps/conv_diff/conv_diff.lua, important changes are marked with PARALLEL SCALING
--
]]--
--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")

dim = util.GetParamNumber("-dim", 2, "dimension")
numPreRefs = util.GetParamNumber("-numPreRefs", 0, "number of refinements before parallel distribution")
numRefs    = util.GetParamNumber("-numRefs",    5, "number of refinements")
startTime  = util.GetParamNumber("-start", 0.0, "start time")
endTime    = util.GetParamNumber("-end", 1e-2, "end time")
dt 		   = util.GetParamNumber("-dt", 1e-3, "time step size")

 
bParallelScalingTest = true

if bParallelScalingTest then
	-- PARALLEL SCALING -{
	-- our app name
	myAppName="LaplaceScalingTest"
	
	-- the output directory we are using to write stats and results to. default is ~/results/
	outdir     = util.GetParam("-outdir", os.getenv("HOME").."/results/")
	-- assure that the directory exists
	CreateDirectory(outdir)
	
	-- we want to have a descriptive directory name for our results directory
	-- so we construct one from the command line arguments	
	commandLineFilename = util.GetUniqueFilenameFromCommandLine()
	
	-- a path for our results of the current run
	theResultPath = outdir..myAppName..commandLineFilename.."/"
	CreateDirectory(theResultPath)
	
	-- enable file logging into our results path
	GetLogAssistant():enable_file_output(false, "")
	GetLogAssistant():enable_file_output(true, theResultPath.."log.txt")
	
	-- the shared stats name where all runs are saving their stats
	-- in a CSV-format
	theSharedStatsName = outdir..myAppName..".stats.txt"
	-- PARALLEL SCALING }-
end


util.CheckAndPrintHelp("Parallel Scaling example\nMartin Rupp, Andreas Vogel");

if dim == 2 then gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_2x2.ugx")
else print("Dimension "..dim.." not supported."); exit(); end

print(" Choosen Parater:")
print("    numRefs      = " .. numRefs)
print("    numPreRefs   = " .. numPreRefs)
print("    startTime 	= " .. startTime)
print("    endTime 		= " .. endTime)
print("    dt 			= " .. dt)
print("    grid         = " .. gridName)

-- choose algebra
InitUG(dim, AlgebraType("CPU", 1));

-- Create, Load, Refine and Distribute Domain
neededSubsets = {"Inner", "Boundary"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets)

-- create Approximation Space
print(">> Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)

-- lets order indices using Cuthill-McKee
OrderCuthillMcKee(approxSpace, true);

--------------------------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
--------------------------------------------------------------------------------

print (">> Setting up Assembling")

-- This scales the amount of diffusion of the problem
eps = 1e-1

-- The coordinates (cx, cy) specify the rotation center of the cone
cx = 0.5
cy = 0.5

-- The coordinates (ax, ay) specify the position of the highest point of the
-- cone at start time t=0.0
ax = 0.25
ay = 0.0

-- The parameter nu specifies the rotation velocity
nu = 100

-- The parameter delta is a scaling factor influencing the steepness of the cone
delta = 1e-1

-- This is the exact solution for our problem
function exactSolution(x, y, t)
	local xRot = math.cos(nu*t) * (x-cx) - math.sin(nu*t) * (y-cy) 
	local yRot = math.sin(nu*t) * (x-cx) + math.cos(nu*t) * (y-cy) 
	
	local expo = -((xRot - ax)*(xRot - ax) + (yRot - ay)*(yRot - ay)) / (delta + 4*eps*t)
	local scale = delta/(delta+4*eps*t)

	return scale * math.exp(expo)
end
	
-- The velocity field
function Velocity(x, y, t)
	return	nu*(y - cx), nu*(cy - x)
end
	
-- The dirichlet condition
function DirichletValue(x, y, t)
	return exactSolution(x, y, t)
end

elemDisc = ConvectionDiffusion("c", "Inner", "fv1")
elemDisc:set_upwind(FullUpwind())
elemDisc:set_diffusion(eps)
elemDisc:set_velocity("Velocity")

dirichletBND = DirichletBoundary()
dirichletBND:add("DirichletValue", "c", "Boundary")

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
domainDisc:add(dirichletBND)

--------------------------------------------------------------------------------
--  Algebra
--------------------------------------------------------------------------------
print (">> Setting up Algebra Solver")

linSolver = util.GetSolver( {
		name = "linear",
		precond = {
			name = "gmg",
			gmg_approxSpace = approxSpace,
			gmg_smoother = {
				name = "jac",
				jac_damp = 0.8 } },
		convCheck = {
			maxSteps = 100,
			minDef = 1e-9,
			reduction = 1e-12 }
	} )

-- create Newton Solver
newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(linSolver)
newtonSolver:set_convergence_check(ConvCheck(10, 5e-8, 1e-10, true))

--------------------------------------------------------------------------------
--  Apply Solver
--------------------------------------------------------------------------------

-- set initial value
print(">> Interpolation start values")
u = GridFunction(approxSpace)
Interpolate("exactSolution", u, "c", startTime)

-- perform time loop
--util.SolveNonlinearTimeProblem(u, domainDisc, newtonSolver, VTKOutput(), "Sol",
--							   "ImplEuler", 1, startTime, endTime, dt); 

util.SolveLinearTimeProblem(u, domainDisc, linSolver, VTKOutput(), "Sol",
							"ImplEuler", 1, startTime, endTime, dt); 


print("Writing data")

if bParallelScalingTest then
	-- PARALLEL SCALING -{
	if GetProfilerAvailable() == true then
	
		WriteProfileData(theResultPath.."profile_data.pdxml")
		lsApplyReturnDefectPN = GetProfileNode("LS_ApplyReturnDefect")
		
		-- the next lines will construct a tab-separated file like this
		-- which is easily readable by excel/numbers
--[[
procs	numPreRefs	numRefs	startTime	endTime	dt  ...
1	3	4	0.0 0.01 0.001 ...
]]-- 
		-- it is written to theSharedStatsName, so that all runs add one line.
		-- below is only an example
		stats = 
		{	
			{ "procs", GetNumProcesses() },
			{ "numPreRefs", numPreRefs},
			{ "numRefs", numRefs },
			{ "startTime", startTime},
			{ "endTime", endTime},
			{ "dt", dt},
			
			-- { "yourParameterName", yourParameter},
			-- { "yourResult", yourResult},
			
			-- total time this run took
			{ "totalTimeS", GetProfileNode("main"):get_avg_total_time_ms()/1000.0 },
			-- time spend in the Linear Solver.
			-- note that only the TimePerCall will scale for nonlinear methods
			-- because the number of calls will differ for different numProcs.
			{ "LS_ApplyReturnDefect_TotalS", lsApplyReturnDefectPN:get_avg_total_time_ms()/1000.0 },
			{ "LS_ApplyReturnDefect_Calls", lsApplyReturnDefectPN:get_avg_entry_count() },
			{ "LS_ApplyReturnDefect_TimePerCallS", lsApplyReturnDefectPN:get_avg_total_time_ms()/1000.0/ (lsApplyReturnDefectPN:get_avg_entry_count()) },
			
			-- other useful information	 
			{ "date", os.date("y%Ym%md%d") },
			{ "SVN Revision", GetSVNRevision()},
			{"host",GetBuildHostname()},
			{"commandline", util.GetCommandLine() } 
		} 
		
		util.printStats(stats)
		bWriteStats=true
		if bWriteStats  then	
			util.writeFileStats(stats, theSharedStatsName)
		end
	end
	
	-- PARALLEL SCALING }-
end

FreeUserData()

-- end group app_convdiff
--[[!
\}
]]--
