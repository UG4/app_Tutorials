--[[!----------------------------------------------------------------------------
--	\file tut10_doxygen.lua
--
--	The purpose of this tutorial is to demonstrate how you can document
--	your lua code so that you get a nice html documentation with doxygen.
--  You can get a list of commands at http://www.stack.nl/~dimitri/doxygen/commands.html
--  You can also take a look at \ref pageUG4DoxygenQuickref and \ref pageUG4CodingStyle
]]------------------------------------------------------------------------------


-- Doxygen is nice!
-- We've written a wrapper which is able to process LUA code so
-- that it is readable by Doxygen (it is in ug4/docs/lua2doxygen )

--! this is our doxygen comment in LUA
--[[! 
this is our multi-line doxygen comment
]]--

-- we can use all doxygen-commands!

-- as you can see at the beginning of the document, we also documented the
-- document itself with \file .
-- most helpful to document functions: \param \brief \return and \sa
-- example:

--! fill
--! \brief a function to create a string by repeating a character
--! 
--! This function is really useful for formatting tables.
--! \param N number of times c is to be repeated
--! \param c character to repeat (if omitted, " ")
--! \return a string consisting of N times the character c
--! \note this is a copy of util.fill
--! \sa stats_util.lua
function fill(N, c)
	local s=""
	if c == nil then c = " " end
	for i=1,N do
		s=s..c
	end
	return s
end

-- we can also use latex:

--! \brief calculates the probability density function
--! \param x
--! \param mu \f$\mu\f$ = mean (location of the peak) default 0
--! \param sigma \f$\sigma^2\f$ = variance. default 1
--! \return \f$\frac{1}{\sqrt{2\pi\sigma^2}}\, e^{\frac{-(x-\mu)^2}{2\sigma^2}}\f$
function pdf(x, sigma, mu)
	return 1/(sigma*math.sqrt(2*math.pi)) * math.exp(-0.5*math.pow((x-mu)/sigma, 2))
end

-- we can also document our command-line-parameters
-- we need ug_util.lua for GetParam.
ug_load_script("ug_util.lua")

mu    = util.GetParamNumber("-mu",    0, "mean (location of the peak")
sigma = util.GetParamNumber("-sigma", 3, "variance")

-- these parameters will now be visible in the lua-documentation.
-- now at this point in the execution we've checked for all parameters, 
-- so we can print a list if the users adds a -help . 
-- The string in the brackets is a description of the script.

util.CheckAndPrintHelp("Doxygen tutorial and nice plotter")

-- now for the fun part!

for x = -9, 9, 0.5 do
	print(fill(200*pdf(x,sigma, mu), "-").."*")
end
