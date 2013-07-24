--------------------------------------------------------------------------------
--[[!
-- \file apps/conv_diff/laplace.lua
-- \ingroup app_convdiff
-- \{
-- \author Andreas Vogel 
-- \brief Lua - Script to perform the Laplace-Problem
--
]]--
--------------------------------------------------------------------------------


--------------------------------
-- User Data Functions
--------------------------------
function Diffusion2d(x, y, t)
	return	1, 0, 
			0, 1
end

function Velocity2d(x, y, t)
	return	0, 0
end

function ReactionRate2d(x, y, t)
	return	0
end

function Source2d(x, y, t)
	local s = 2*math.pi
	return s*s*(math.sin(s*x) + math.sin(s*y))
end

function DirichletValue2d(x, y, t)
	local s = 2*math.pi
	return true, math.sin(s*x) + math.sin(s*y)
end

ug_load_script("init.lua")

-- creates domainDisc, approxSpace

--------------------------------------------------------------------------------
--  Algebra
--------------------------------------------------------------------------------
print ("Setting up Algebra Solver")




-- create algebraic Preconditioner
jac = Jacobi()

-- create Convergence Check
convCheck = ConvCheck()
convCheck:set_maximum_steps(100)
convCheck:set_minimum_defect(1e-11)
convCheck:set_reduction(1e-12)



-- 0. Reset start solution
-- create operator from discretization
for i=1,4 do
	-- create Linear Solver
	linSolver = LinearSolver()
	linSolver:set_preconditioner(jac)
	linSolver:set_convergence_check(convCheck)

	linOp = AssembledLinearOperator(domainDisc)
	
	-- get grid function
	u = GridFunction(approxSpace)
	b = GridFunction(approxSpace)
	
	u:set(0.0)
	
	-- 1. init operator
	print("Init operator (i.e. assemble matrix).")
	
	AssembleLinearOperatorRhsAndSolution(linOp, u, b)
	ApplyLinearSolver(linOp, u, b, linSolver)
	
	convCheck:reduction()
	refiner:refine()
end

SaveVectorForConnectionViewer(u, "uSol.vec")
