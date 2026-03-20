"""
Conformational optimizers for BuildAMol. 
"""

from importlib import import_module


__all__ = [
	"Rotatron",
	"DistanceRotatron",
	"OverlapRotatron",
	"ForceFieldRotatron",
	"Translatron",
	"Circulatron",
	"ConstraintRotatron",
	"simple_concatenation_function",
	"concatenation_function_with_penalty",
	"concatenation_function_no_pushback",
	"concatenation_function_no_unfold",
	"concatenation_function_linear",
	"jensen_shannon_overlap",
	"MVN",
	"swarm_optimize",
	"genetic_optimize",
	"scipy_optimize",
	"anneal_optimize",
	"rdkit_optimize",
	"mmff_optimize",
	"uff_optimize",
	"apply_rotatron_solution",
	"apply_translatron_solution",
	"optimize",
	"auto_algorithm",
	"auto_applier",
	"split_environment",
	"parallel_optimize",
]


_LAZY_IMPORTS = {
	"Rotatron": ("buildamol.optimizers.base_rotatron", "Rotatron"),
	"DistanceRotatron": ("buildamol.optimizers.distance_rotatron", "DistanceRotatron"),
	"OverlapRotatron": ("buildamol.optimizers.overlap_rotatron", "OverlapRotatron"),
	"ForceFieldRotatron": (
		"buildamol.optimizers.forcefield_rotatron",
		"ForceFieldRotatron",
	),
	"Translatron": ("buildamol.optimizers.translatron", "Translatron"),
	"Circulatron": ("buildamol.optimizers.circulatron", "Circulatron"),
	"ConstraintRotatron": (
		"buildamol.optimizers.constraint_rotatron",
		"ConstraintRotatron",
	),
	"simple_concatenation_function": (
		"buildamol.optimizers.distance_rotatron",
		"simple_concatenation_function",
	),
	"concatenation_function_with_penalty": (
		"buildamol.optimizers.distance_rotatron",
		"concatenation_function_with_penalty",
	),
	"concatenation_function_no_pushback": (
		"buildamol.optimizers.distance_rotatron",
		"concatenation_function_no_pushback",
	),
	"concatenation_function_no_unfold": (
		"buildamol.optimizers.distance_rotatron",
		"concatenation_function_no_unfold",
	),
	"concatenation_function_linear": (
		"buildamol.optimizers.distance_rotatron",
		"concatenation_function_linear",
	),
	"jensen_shannon_overlap": (
		"buildamol.optimizers.overlap_rotatron",
		"jensen_shannon_overlap",
	),
	"MVN": ("buildamol.optimizers.overlap_rotatron", "MVN"),
	"swarm_optimize": ("buildamol.optimizers.algorithms", "swarm_optimize"),
	"genetic_optimize": ("buildamol.optimizers.algorithms", "genetic_optimize"),
	"scipy_optimize": ("buildamol.optimizers.algorithms", "scipy_optimize"),
	"anneal_optimize": ("buildamol.optimizers.algorithms", "anneal_optimize"),
	"rdkit_optimize": ("buildamol.optimizers.algorithms", "rdkit_optimize"),
	"mmff_optimize": ("buildamol.optimizers.algorithms", "mmff_optimize"),
	"uff_optimize": ("buildamol.optimizers.algorithms", "uff_optimize"),
	"apply_rotatron_solution": (
		"buildamol.optimizers.utils",
		"apply_rotatron_solution",
	),
	"apply_translatron_solution": (
		"buildamol.optimizers.utils",
		"apply_translatron_solution",
	),
	"optimize": ("buildamol.optimizers.utils", "optimize"),
	"auto_algorithm": ("buildamol.optimizers.utils", "auto_algorithm"),
	"auto_applier": ("buildamol.optimizers.utils", "auto_applier"),
	"split_environment": ("buildamol.optimizers.utils", "split_environment"),
	"parallel_optimize": ("buildamol.optimizers.utils", "parallel_optimize"),
}


def __getattr__(name):
	target = _LAZY_IMPORTS.get(name)
	if target is None:
		raise AttributeError(f"module '{__name__}' has no attribute '{name}'")
	module_name, attr_name = target
	module = import_module(module_name)
	value = getattr(module, attr_name)
	globals()[name] = value
	return value


def __dir__():
	return sorted(set(globals()) | set(__all__))
