--------------------------------------------------------------------------------
--	tut12_refinement_projectors.lua
--
--	This tutorial demonstrates the usage of different refinement projectors.
--	Refinement projectors are important if calculations are to be performed on
--	spherical or cylindrical domains or if domains contain spherical or cylindrical
--	parts.
--
--	ATTENTION: When applying a spherical projector to the whole geometry, one should
--	make sure that the center of the projection aligns with a vertex or is located
--	in a cutout.
--	The same holds true for a cylinder. Here either the axis of the cylinder should
--	pass through a set of edges or should also be located in a cutout.
--------------------------------------------------------------------------------

-- include the basic util-methods.
ug_load_script("ug_util.lua")

-- Get the command line parameters
numRefs = util.GetParamNumber("-numRefs", 3)

-- This function will refine it using the given projector and
-- save the results to the specified filename. With zOffset an offset along
-- the z-axis for each level may be specified.
function PerformAdaption(dom, outFile, projector, zOffset)
	local refiner = GlobalDomainRefiner(dom)
	refiner:set_refinement_callback(projector)
	write("refining")
	for i = 1, numRefs do
		write(" " .. i)
		refiner:refine()
	end
	write("\n")
	print("saving to " .. outFile)
	SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), outFile, zOffset)
	delete(refiner)
	print("done\n")
end


--	Call the adaption method for different grids and different projectors
if IsDefinedUG_DIM_2() then
	--	We'll start with a simple circle and use a sphere projector on the whole domain.
	--	Arguments to the sphere projector are the domain and it's center.
	local dom = Domain2d()
	print("loading " .. "grids/circle.ugx")
	LoadDomain(dom, "grids/circle.ugx")
	PerformAdaption(dom, "circle_out.ugx",
					SphereProjector(dom, {0, 0}), 1)
	
	
	--	Next we'll refine a rect with a single cutout. We want the cutout to converge
	--	to a circular cutout with each refinement step. We thus use the
	--	SphericalFalloffProjector. It takes the domain, the center, an inner radius
	--	and an outer radius as parameters. All vertices refined inside the inner radius
	--	are handled with spherical projection, all vertices outside outer radius are
	--	handled with linear projection. In between a smooth transition between
	--	spherical and linear projection is performed.
	--	In this case the center is at (0, 0), the innerRadius is 0.5 and the
	--	outerRadius is 1
	dom = Domain2d()
	print("loading " .. "grids/rect_with_cutout.ugx")
	LoadDomain(dom, "grids/rect_with_cutout.ugx")
	PerformAdaption(dom, "rect_with_cutout_out.ugx",
					SphericalFalloffProjector(dom, {0, 0}, 0.5, 1), 1)
	
	
	--	Now we'll use two different projectors to refine a geometry with two
	--	cutouts. We'll use a DomainRefinementProjectionHandler, which will
	--	associate different projectors with the different subsets of the domain.
	dom = Domain2d()
	print("loading " .. "grids/rect_with_two_cutouts.ugx")
	LoadDomain(dom, "grids/rect_with_two_cutouts.ugx")
	local projHandler = DomainRefinementProjectionHandler(dom)
	projHandler:set_callback("Inner1", SphericalFalloffProjector(dom, {0, 0}, 0.5, 1))
	projHandler:set_callback("Inner2", SphericalFalloffProjector(dom, {2, 0}, 0.5, 1))
	PerformAdaption(dom, "rect_with_two_cutouts_out.ugx", projHandler, 1)
end

if IsDefinedUG_DIM_3() then
	--	We'll load a cube geometry and refine it with a sphere-projector
	--	the result will be a 3d sphere with an inner layer consisting of pyramids
	--	and an outer layer consisting of hexahedrons
	local dom = Domain3d()
	print("loading " .. "grids/cube.ugx")
	LoadDomain(dom, "grids/cube.ugx")
	PerformAdaption(dom, "cube_out.ugx",
					SphereProjector(dom, {0, 0, 0}), 4)


	--	Now we'll load a cylinder and use the cylinder-projector during refinement.
	--	CylinderProjector takes a point on the axis and the axis direction. In this case
	--	(0, 0, 0) and (0, 0, 1).
	local dom = Domain3d()
	print("loading " .. "grids/cylinder.ugx")
	LoadDomain(dom, "grids/cylinder.ugx")
	PerformAdaption(dom, "cylinder_out.ugx",
					CylinderProjector(dom, {0, 0, 0}, {0, 0, 1}), 3)


	--	If only a part of the domain shall be refined as a cylinder, we'll
	--	use the CylindricalFalloffProjector. Additionally to the center and
	--	the direction of the axis, one has to specify an innerRadius and
	--	an outerRadius. Everything inside innerRadius is refined using
	--	cylinder projection, everything outside outerRadius is refined using
	--	linear projection. A smooth transition between the two projection methods
	--	is performed in between
	dom = Domain3d()
	print("loading " .. "grids/box_with_cutout.ugx")
	LoadDomain(dom, "grids/box_with_cutout.ugx")
	PerformAdaption(dom, "box_with_cutout_out.ugx",
					CylindricalFalloffProjector(dom, {0, 0, 0}, {0, 0, 1}, 0.5, 1), 3)
end
