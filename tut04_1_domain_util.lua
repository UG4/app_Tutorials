--------------------------------------------------------------------------------
--	tut04_1_domain_util.lua
--
--	This file is part of tut04_modular_programming.lua.
--
--	This file contains the method CreateAndDistributeDomain.
--	The method which will create a domain, load the associated
--	geometry from a file and distriutes it onto all active processes.
--	Please note, that this could also be achieved with the method
--	util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets),
--	which is a little more elaborated, since it also supports refienement
--	(before and after distribution).
--------------------------------------------------------------------------------

-- include the basic util-methods.
ug_load_script("ug_util.lua")


--------------------------------------------------------------------------------
--	The method loads creates a new domain, loads the grid specified in gridName,
--	distributes it to all active processes and checks, whether all processes
--	received the required subsets.
--	
--	If everything went right, the method returns the created domain object,
--	if an error occured, nil is returned.
--	
--	Params: string gridName, string-array requiredSubsets
--	Returns: Domain
--
function CreateAndDistributeDomain(gridName, requiredSubsets)
--	Create the domain object
	local dom = Domain()
	
--	Load the domain from file (Note that this is always performed on
--	process 0)
	LoadDomain(dom, gridName)
		
--	Distribute the domain to all involved processes
	if util.DistributeDomain(dom) == false then
		print("Error while Distributing Domain.")
		return nil
	end
	
--	Check whether each process received the required subsets
	if requiredSubsets ~= nil then
		if util.CheckSubsets(dom, requiredSubsets) == false then 
			print("Something wrong with required subsets. Aborting.");
			return nil
		end
	end
	
--	we're done. Return the domain.
	return dom
end
