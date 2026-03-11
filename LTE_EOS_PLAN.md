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
