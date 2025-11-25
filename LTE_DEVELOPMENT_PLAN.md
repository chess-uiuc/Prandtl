# Phase 1

---

## Track A — EOS/Gas-Model Plan & Design

### Goals
- Design an abstraction for encapsulating EOS/Transport/Physics logic  
- Prepare for LTE capability by removing scattered ideal-gas assumptions  
- Craft a (markdown) Plan/Design document

### Steps
- **Development process**
  - [ ] Create an Issue for tracking EOS/Gas-Model Planning
  - [ ] Create a top-level branch for EOS/Gas-Model Planning   
  - [ ] Create a (draft) PR that closes the Issue (target chess-uiuc/Prandtl@add-lte-gas)

- **Create a "thermo usage inventory"**
  - [ ] List all thermal-state-related `compute_*` routines  
  - [ ] Record inputs (conserved or primitive)  
  - [ ] Record call sites where each is used (flux, BC, sensors, parabolic terms, etc.)

- **Inventory assumptions about state layout** [`[Q1]`](https://github.com/chess-uiuc/Prandtl/blob/lte-gpu-integration/DEVELOPMENT_PLAN.md#q1--state-layout-stability-high-priority-before-development)
  - rho, momentum, energy ordering  
  - Future scalar transport and mixture concerns `[c:1B]`  
  - [ ] **Determine/document state interface requirements**  
  - [ ] **Determine/document needs for persistent auxiliary state data** [`[Q2]`](#q2--persistent-state-requirements)

- **Identify MFEM Euler helper usage**
  - Euler/inviscid fluxes (for example) 
  - These must be replaced by Prandtl-native kernels to support mixtures/scalars `[c:1B]`

- [ ] Prepare testing extensions as needed
  - EOS unit testing likely needed

### Deliverables
- [ ] Merge-ready PR for Track A Phase 1  
- [ ] Markdown/table summarizing thermo-usage and state-layout assumptions  
- [ ] Documented Track A plan for:
  - Adding EOS enapsulation
  - Adding GasModel construct

### Outcomes
- Clarity on what must be encapsulated before EOS/LTE work  
- Shared understanding of state-layout constraints  
- Regression tests ensure safe development
- Unit tests to ensure correctness, definitiion-of-done

---

# Phase 2

---

## Track A — EOS/Gas-Model Development (Ideal Gas First)

### Goals
- Encapsulate ideal-gas EOS into single abstraction  
- Replace scattered thermodynamic calls with GasModel/EOS usage  

### Steps
- **Development process**
  - [ ] Make Issue to track adding EOS/GasModel constructs  
  - [ ] Create development branch for this work  
  - [ ] Draft a PR to close the Issue (target chess-uiuc/Prandtl@add-lte-gas branch)  

- [ ] **Define minimal GasModel interface**
  - Cons↔Prim transforms  
  - Pressure, temperature, sound speed  
  - **NOTE:** Entropy conversions MAY live here (decision deferred) [`[Q3]`](https://github.com/chess-uiuc/Prandtl/blob/lte-gpu-integration/DEVELOPMENT_PLAN.md#q3--entropy-conversion-plumbing)

- [ ] Add EOS unit tests  
- [ ] Implement IdealGas EOS (must match existing behavior)  
- [ ] Introduce legacy shims (`compute_* → gas_model.*`)  
- [ ] Switch BCs/fluxes to use gas_model.* `[c:1B]`  

### Deliverables
- [ ] Merge-ready PR with IdealGas EOS/GasModel  
- [ ] GasModel/EOS classes, unit tests, shims  

### Outcomes
- Single source of truth for thermodynamics  
- Ready for LTE work without physics changes  

---

# Phase 3

---

## Track A — Add LTE EOS Capability Behind GasModel

### Goals
- Extend GasModel for LTE thermodynamics  
- Add LTE internal variables & closures without touching DGSEM internals  

### Steps
- **Development process**
  - Make Issue, branch, PR  

- Extend GasModel/EOS:
  - LTE internal modes, electronic states  
  - LTE closure for p, T, h, c_s  
  - LTE source terms  

- Handle Mutation++ carefully:
  - **Cannot be called from GPU**  
  - Must not appear in device paths  
  - Host-side or tabulated data fed to kernels [`[Q4]`](https://github.com/chess-uiuc/Prandtl/blob/lte-gpu-integration/DEVELOPMENT_PLAN.md#q4--mutation-and-gpu)  

- Add LTE examples/tests  

### Deliverables
- Merge-ready PR with LTE support  
- LTE examples  

### Outcomes
- LTE gas supported cleanly  
- GPU readiness maintained  

---
