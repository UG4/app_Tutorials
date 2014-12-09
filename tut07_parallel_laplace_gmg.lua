--------------------------------------------------------------------------------
--	tut07_parallel_laplace_gmg.lua
--
-- In this tutorial a distributed domain is set up that has several vertical
-- interfaces. Those interfaces are created when a part of a grid is 
-- distributed to other processes and further refined. Several interfaces 
-- are e.g. created when a tree-like distribution is performed.
-- 
-- to run this code in parallel, use e.g.
-- mpirun -np 4 ugshell -ex tutorials/tut07_parallel_laplace_gmg.lua
--------------------------------------------------------------------------------
-- include the basic util-methods.
ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

-- Get the command line parameters
dim = util.GetParamNumber("-dim", 2)

-- Since ug supports a bunch of different dimensions and algebra modules 
-- we will choose a combination here. This should always be the first thing 
-- you do in an ug-script. The cpu-algebra is fine for now.
InitUG(dim, AlgebraType("CPU", 1))

gridName = util.GetParam("-grid", "unit_square/unit_square_quads_2x2.ugx")
outFileNamePrefix = util.GetParam("-o", "distributed_domain_")
numRefs = util.GetParamNumber("-numRefs", 4)
debug = util.HasParamOption("-debug")

-- We will save the created hierarchy to this file (with appended process id)
outHierarchyFilePrefix = util.GetParam("-oh", "hierarchy_on_proc_")



-- while this is not necessary, we can overwrite some of the default parameters
-- used for load-balancing.
-- The redistProcs parameter defines the maximum number of processes to which
-- the grid is distributed from a single process. If one runs the application on
-- 'redistProcs' this will lead to a process hierarchy, which is a nice thing
-- for massively parallel runs.
balancer.redistProcs = 64

-- we have to give the balancer a chance to initialize its parameters from
-- command line arguments:
balancer.ParseParameters()



-- Create the domain and load a grid
dom = Domain()
print("loading grid from " .. gridName .. "...")
LoadDomain(dom, gridName)

-- This will automatically refine and rebalance the domain. If you need more
-- control, please have a look at the implementation of this method in
-- scripts/util/load_balancing_util.lua and extract the pieces you need.
-- You can e.g. create your own LoadBalancer with custom weightings etc and use
-- that during redistribution or use a different refinement scheme.
balancer.RefineAndRebalanceDomain(dom, numRefs)

-- since we used standard methods, we can now query the default load balancer
-- for distribution qualities
if balancer.defaultBalancer ~= nil then
	balancer.defaultBalancer:print_quality_records()
end

-- save the domain on each process
outFileName = outFileNamePrefix .. ProcRank() .. ".ugx"
print("saving domain to " .. outFileName .. "...")
SaveDomain(dom, outFileName)



-- Now lets save the hierarchy on each process
-- The SaveGridHierarchyTransformed routine directly works on the domains grid.
-- We can access the grid of a domain through its grid() member method.
--
-- SaveGridHierarchyTransformed outputs a grid, where each level has an offset
-- along the z-axis.
zOffset = 1
outFileName = outHierarchyFilePrefix .. ProcRank() .. ".ugx"
print("saving grid hierarchy to " .. outFileName)
SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), outFileName, zOffset)


--------------------------------------------------------------------------------
--	Now setup the laplace problem and solve it

-- Create, Load, Refine and Distribute Domain
neededSubsets = {"Inner", "Boundary"}

if util.CheckSubsets(dom, neededSubsets) == false then
	print("Some subsets are incorrect. Aborting."); exit();
end


function DirichletValue(x, y, t)
	return x*x + y*y
end

approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)

-- lets order indices using Cuthill-McKee
OrderCuthillMcKee(approxSpace, true);

elemDisc = ConvectionDiffusion("c", "Inner", "fv1")
elemDisc:set_upwind(FullUpwind())
elemDisc:set_diffusion(1)

dirichletBND = DirichletBoundary()
dirichletBND:add("DirichletValue", "c", "Boundary")

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
domainDisc:add(dirichletBND)


-- Create surface functions (vectors) for Au=b and initialize them
A = MatrixOperator()
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)


-- set initial value
u:set(0)
domainDisc:adjust_solution(u)

-- assemble matrix and rhs
domainDisc:assemble_linear(A, b)

-----------------
-- create GMG 
-----------------

-- create algebraic Preconditioner
jac = Jacobi()
jac:set_damp(0.8)

-- create Base Solver
base = LU()
	
-- create Geometric Multi Grid
gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_base_level(0)
gmg:set_base_solver(base)
gmg:set_smoother(jac)
gmg:set_cycle_type(1)
gmg:set_num_presmooth(3)
gmg:set_num_postsmooth(3)


-- create Convergence Check
convCheck = ConvCheck()
convCheck:set_maximum_steps(10000)
convCheck:set_minimum_defect(1e-11)
convCheck:set_reduction(1e-12)

-- create Linear Solver
solver = LinearSolver()
solver:set_preconditioner(gmg)
solver:set_convergence_check(convCheck)

-------------------------------------------
--  Apply Solver
-------------------------------------------
print("Init solver for operator.")
solver:init(A)

print("Apply solver.")
solver:apply_return_defect(u,b)

WriteGridFunctionToVTK(u, "Solution")

-- we're done.
print("")
print("done")
