import Lake
open Lake DSL

package havelock where
  leanOptions := #[⟨`autoImplicit, false⟩, ⟨`relaxedAutoImplicit, false⟩]

require "leanprover-community" / "mathlib" @ git "v4.34.1"

@[default_target]
lean_lib Havelock where
