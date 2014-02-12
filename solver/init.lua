
PrintBuildConfiguration()

ug_load_script("ug_util.lua")
ug_load_script("solver_util/setup_rsamg.lua")

--------------------------------------------------------------------------------
-- Checking for parameters (begin)
--------------------------------------------------------------------------------
-- Several definitions which can be changed by command line parameters
-- space dimension and grid file:

gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_2x2.ugx")
numPreRefs = util.GetParamNumber("-numPreRefs", 0, "number of refinements before parallel distribution")
numRefs    = util.GetParamNumber("-numRefs",    3, "number of refinements")

if numPreRefs > numRefs then
	print("It must be choosen: numPreRefs <= numRefs");
	exit();
end

-- parallelisation related stuff
--! way the domain / the grid will be distributed to the processes:
distributionType = util.GetParam("-distType", "bisect", "way the domain will be distributed: [grid2d | bisect | metis]")

--! amount of output
--! set to 0 i.e. for time measurements,
--! >= 1 for writing matrix files etc.
verbosity = util.GetParamNumber("-verb", 0, "verbosity. 0 for time measurements")		

activateDbgWriter = 0	  
--! set to 0 i.e. for time measurements,
--! >= 1 for debug output: call 'set_debug(dbgWriter)'
--! for the main solver ('gmg')
activateDbgWriter = util.GetParamNumber("-dbgw", 0, "debug writer: 0 for no output.") 
						    
bUseAMG = util.HasParamOption("-amg", "use amg")


util.CheckAndPrintHelp("Laplace-Problem\nMartin Rupp");

-- Display parameters (or defaults):
print(" General parameters chosen:")
print("    grid       = " .. gridName)
print("    numRefs    = " .. numRefs)
print("    numPreRefs = " .. numPreRefs)

print("    verb (verbosity)         = " .. verbosity)
print("    dbgw (activateDbgWriter) = " .. activateDbgWriter)

print(" Parallelisation related parameters chosen:")
print("    distType   = " .. distributionType)

--------------------------------------------------------------------------------
-- Checking for parameters (end)
--------------------------------------------------------------------------------

-- choose algebra
InitUG(2, AlgebraType("CPU", 1));

-- create Instance of a Domain
print("Create Domain.")
dom = Domain()

-- load domain
print("Load Domain from File.")
LoadDomain(dom, gridName)

-- create Refiner
print("Create Refiner")
-- Create a refiner instance. This is a factory method
-- which automatically creates a parallel refiner if required.
refiner = GlobalDomainRefiner(dom)

-- Performing pre-refines
print("Perform (non parallel) pre-refinements of grid")
for i=1,numPreRefs do
	write( "PreRefinement step " .. i .. " ...")
	refiner:refine()
	print( " done.")
end

-- get number of processes
numProcs = NumProcs()

-- Distribute the domain to all involved processes
-- Since only process 0 loaded the grid, it is the only one which has to
-- fill a partitionMap (but every process needs one and has to return his map
-- by calling 'DistributeDomain()', even if in this case the map is empty
-- for all processes but process 0).
if numProcs > 1 then
	print("Distribute domain with 'distributionType' = " .. distributionType .. "...")
	partitionMap = PartitionMap()
	
	if ProcRank() == 0 then
		if distributionType == "bisect" then
			util.PartitionMapBisection(dom, partitionMap, numProcs)
			
		elseif distributionType == "grid2d" then
			local numNodesX, numNodesY = util.FactorizeInPowersOfTwo(numProcs / numProcsPerNode)
			util.PartitionMapLexicographic2D(dom, partitionMap, numNodesX,
											 numNodesY, numProcsPerNode)
	
		elseif distributionType == "metis" then
			util.PartitionMapMetis(dom, partitionMap, numProcs)
											 
		else
		    print( "distributionType not known, aborting!")
		    exit()
		end

	-- save the partition map for debug purposes
		if verbosity >= 1 then
			print("saving partition map to 'partitionMap_p" .. ProcRank() .. ".ugx'")
			SavePartitionMap(partitionMap, dom, "partitionMap_p" .. ProcRank() .. ".ugx")
		end
	end
	
	print("Redistribute domain with 'distributionType' = '" .. distributionType .. "' ...")
	if DistributeDomain(dom, partitionMap, true) == false then
		print("Redistribution failed. Please check your partitionMap.")
		exit()
	end
	print("... domain distributed!")
	delete(partitionMap)
end


--------------------------------------------------------------------------------
-- end of partitioning
--------------------------------------------------------------------------------

-- Perform post-refine
print("Refine Parallel Grid")
for i=numPreRefs+1,numRefs do
	write( "Refinement step " .. i .. " ...")
	refiner:refine()
	print( " done!")
end


-- Make sure, that the required subsets are present
requiredSubsets = {"Inner", "Boundary"}
if util.CheckSubsets(dom, requiredSubsets) == false then 
   print("Subsets missing. Aborting")
   exit()
end

-- write grid to file for test purpose
if verbosity >= 1 then
	refinedGridOutName = "refined_grid_p" .. ProcRank() .. ".ugx"
	print("saving domain to " .. refinedGridOutName)
	SaveDomain(dom, refinedGridOutName)
	
	hierarchyOutName = "hierachy_p" .. ProcRank() .. ".ugx"
	print("saving hierachy to " .. hierarchyOutName)
	if SaveGridHierarchy(dom:grid(), hierarchyOutName) == false then
		print("Saving of domain to " .. hierarchyOutName .. " failed. Aborting.")
		    exit()
	end
end

print("NumProcs is " .. numProcs .. ", numPreRefs = " .. numPreRefs .. ", numRefs = " .. numRefs .. ", grid = " .. gridName)

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_local_dof_statistic(2)
approxSpace:print_layout_statistic()
approxSpace:print_statistic()

-- lets order indices using Cuthill-McKee
OrderCuthillMcKee(approxSpace, true);

--------------------------------------------------------------------------------
--  Assembling
--------------------------------------------------------------------------------
print ("Setting up Assembling")
--------------------------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
--------------------------------------------------------------------------------

elemDisc = ConvectionDiffusion("c", "Inner", "fv1")
elemDisc:set_diffusion("Diffusion2d")
elemDisc:set_velocity("Velocity2d")
elemDisc:set_reaction_rate("ReactionRate2d")
elemDisc:set_source("Source2d")

dirichletBND = DirichletBoundary()
dirichletBND:add("DirichletValue2d", "c", "Boundary")

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
domainDisc:add(dirichletBND)