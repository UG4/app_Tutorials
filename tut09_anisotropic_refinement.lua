--------------------------------------------------------------------------------
--	tut09_anisotropic_refinement.lua
--
--	The purpose of this tutorial is to demonstrate how one can use adaptive
--	refinement to refine geometries in a more sophisticated way.
--	It also shows how the load-balancer can accounter for anisotropies and
--	how the GMG-Solver has to be adjusted in order to do the right transfer.
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

if dim == 2 then
	gridName = util.GetParam("-grid", "anisotropic_rect.ugx")
elseif dim == 3 then
	gridName = util.GetParam("-grid", "anisotropic_hexa.ugx")
else
	print("Only dim == 2 and dim == 3 are supported in the moment")
	exit()
end

outFileNamePrefix = util.GetParam("-o", "distributed_domain_")

-- We will save the created hierarchy to this file (with appended process id)
outHierarchyFilePrefix = util.GetParam("-oh", "hierarchy_on_proc_")

numRefs = util.GetParamNumber("-numRefs", 6)
numAnisoRefs = util.GetParamNumber("-numAnisoRefs", numRefs)


-- while this is not necessary, we can overwrite some of the default parameters
-- used for load-balancing.
-- The redistProcs parameter defines the maximum number of processes to which
-- the grid is distributed from a single process. If one runs the application on
-- 'redistProcs' this will lead to a process hierarchy, which is a nice thing
-- for massively parallel runs.
balancer.redistProcs = 64

-- distribution on a given level is only performed if, if each process could
-- potentially receive at least 'parallelElementThreshold' elements.
balancer.parallelElementThreshold = 8

-- we have to give the balancer a chance to initialize its parameters from
-- command line arguments:
balancer.ParseParameters()



-- Create the domain and load a grid
dom = Domain()
LoadDomain(dom, gridName)
print("Loaded domain from " .. gridName)


-- We will create a refiner now. Since we want to use adaptive refinement
-- we'll use a hanging node refiner
refiner = HangingNodeDomainRefiner(dom)

-- we'll create our own load-balancer, since we want to assign custom balance-weights
loadBalancer = balancer.CreateLoadBalancer(dom)
if loadBalancer ~= nil then
	loadBalancer:set_balance_weights(AnisotropicBalanceWeights(dom))
end

-- perform an initial rebalancing
balancer.Rebalance(dom, loadBalancer)

maxRefEdgeLen = GetMaxEdgeLength(dom) * 0.75
for i = 1, numRefs do
--	apply marks
	print("refining...")
	if(i <= numAnisoRefs) then
		MarkAnisotropic_LongEdges(dom, refiner, maxRefEdgeLen)
		maxRefEdgeLen = maxRefEdgeLen * 0.5
	else
		MarkForRefinement_All(refiner)
	end
	
-- refine
	refiner:refine()

-- rebalance	
	balancer.Rebalance(dom, loadBalancer)
end

if loadBalancer ~= nil then
	loadBalancer:print_quality_records()
end


-- Lets save the domain on each process
outFileName = outFileNamePrefix .. ProcRank() .. ".ugx"
SaveDomain(dom, outFileName) 

-- Everything seems to went fine.
print("Saved domain to " .. outFileName)


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
-- Note the the gmg solver requires a special transfer when we're using ansiotropic elements

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
--	if we're using anisotropic refinement, we currently have disable the
--	p1-lagrange-optimization in the transfer operator. We thus use our own
--	transfer operator
if(numAnisoRefs > 0) then
	local transfer = StdTransfer()
	transfer:enable_p1_lagrange_optimization(false)
	gmg:set_transfer(transfer)
end

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


