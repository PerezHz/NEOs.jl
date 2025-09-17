"""
    AbstractImpactMonitoring

Supertye for the impact monitoring interface.
"""
abstract type AbstractImpactMonitoring end

include("targetplane.jl")
include("lov.jl")
include("closeapproach.jl")