## 🔧 How to Run

### **Step 1 — Open a Julia REPL in this folder**

Go to the directory containing `mpm.jl`, `Project.toml`, `Manifest.toml` etc. then run:

```
julia --project=.
```

---

### **Step 2 — Instantiate the environment (only required to run code for the first time )**

Inside the Julia REPL:

```
import Pkg
Pkg.instantiate()
```

This downloads all required dependencies.

---

### **Step 3 — Run the FEM solver**

Inside the Julia REPL:

```
include("fem.jl")
main()
```

---

or for the mpm 
```
include("mpm.jl")
main()
```

The MPM solver writes:

- Per-particle logs (use XLSX instead of the csv)
- VTK files for Paraview
- All outputs stored inside `MPM_results/`

The FEM solver writes:

- Tip deflection time history (XLSX)
- Expected and final deflection  
- Results stored inside `fem_results/`

---

