# Alphabet Size Diagnostic: Findings and Next Steps

> **Date**: 2026-02-13
> **PR**: #11 (`feat/alphabet-diagnostic`)
> **Scope**: Ceiling gate analysis (alphabets 4/8/16), re-run of Experiments A/B/C with 8-letter alphabet, cross-alphabet statistical comparison

---

## 1. Executive Summary

**Verdict: NO-GO for the current substrate design.**

Expanding the alphabet from 4 to 8 or 16 letters does not resolve the fundamental saturation problem. The root cause is structural: the combination of (a) uniform-random bonding, (b) unconstrained binary tree topology, and (c) an MA metric that counts unique compound subtrees means that random assembly already produces near-maximum MA at every molecule size. No amount of alphabet expansion, catalysis, or ablation can create a gap between "random complexity" and "organized complexity" under these conditions.

To unlock assembly transitions, changes to the **reaction rules**, **tree structure**, or **molecular interaction model** are required.

---

## 2. Diagnostic Data

### 2.1 Ceiling Gate

| Alphabet | Size | P99 MA | Random Est. | Ratio | MA=15 % | GO |
|----------|------|--------|-------------|-------|---------|----|
| ABCD | 4 | 15 | 1.0 | 15.0 | 46% | Yes |
| ABCDEFGH | 8 | 15 | 1.0 | 15.0 | 83% | Yes |
| ABCDEFGHIJKLMNOP | 16 | 15 | 1.0 | 15.0 | 96% | Yes |

The ceiling gate passes for all alphabets, but misleadingly: the ratio is identical (15.0) because P99 and the random equilibrium estimate are the same. Crucially, larger alphabets **concentrate MA even more** at the theoretical maximum. More letters = more distinct leaf labels = fewer accidental subtree collisions = higher MA from pure randomness.

### 2.2 Cross-Alphabet Comparison (800 runs)

| Experiment | Alphabet 4 Mean A_total | Alphabet 8 Mean A_total | Delta | p-value | Effect Size r |
|------------|------------------------|------------------------|-------|---------|---------------|
| Exp A (baseline) | 830.9 | 871.7 | +4.9% | 3.9e-18 | 0.0 |
| Exp B substring | 835.5 | 881.9 | +5.6% | 3.9e-18 | 0.0 |
| Exp B subtree | 835.5 | 881.9 | +5.6% | 3.9e-18 | 0.0 |
| Exp B random_table | 834.0 | 878.7 | +5.4% | 3.9e-18 | 0.0 |
| Exp C stage 1 | 835.5 | 881.9 | +5.6% | 3.9e-18 | 0.0 |
| Exp C stage 4 | 836.6 | 887.5 | +6.1% | 3.9e-18 | 0.0 |

Statistically significant (large N), but effect size is zero. The 5% increase in A_total is entirely explained by larger molecules having more unique subtrees trivially. **Zero transitions detected in any condition, for any alphabet.**

### 2.3 Size-Conditioned Analysis (4-letter alphabet)

From the original Exp A/B size-conditioned analysis:

| Size Bin | A/A_max (baseline) | A/A_max (catalysis) | p-value |
|----------|--------------------|---------------------|---------|
| 2-3 | 1.000 | 1.000 | 1.0 |
| 4-7 | 0.991 | 0.991 | 0.89 |
| 8-11 | 0.978 | 0.976 | 0.64 |
| 12-15 | 0.967 | 0.965 | 0.43 |
| 16+ | 0.926 | 0.921 | 0.14 |

Random trees already achieve 92-100% of maximum MA at every size bin. There is no room for catalysis or selection to create measurably more complex structures.

---

## 3. Root Cause Analysis

The saturation has three interlocking causes:

### 3.1 Bonding Rule: Unconstrained Root-Level Join

The current bonding rule (`chemist.py:89`) always creates a new root node with the two reactants as children:

```
X + Y -> Compound(X, Y)
```

This means every bond event creates a **new unique compound subtree** (the product itself). Since MA counts unique compound subtrees, every bond mechanically increases MA. There is no selectivity -- any pair of molecules can bond, and the result is always structurally novel at the root level.

### 3.2 Tree Topology: Binary Trees Are Too Homogeneous

Binary trees with n leaves have exactly n-1 internal nodes. The number of distinct binary tree topologies (Catalan number) grows much slower than the number of possible labelings (alphabet^n). For a 4-letter alphabet with 16 leaves: ~1.7 billion possible labeled trees, but only 9,694,845 unlabeled topologies (Catalan(15)). With 8 or 16 letters, the labeling space explodes further, making subtree collisions even less likely -- MA saturates closer to n-1.

### 3.3 Breaking Rule: Uniform Random Fragmentation

Breaking selects a random internal node and splits the tree there (`chemist.py:138`). This does not preferentially preserve or destroy high-MA structures. Combined with the bonding rule, the system reaches a trivial equilibrium: molecules grow by random joining, shrink by random breaking, and the resulting size distribution determines A_total mechanically.

**In summary**: the current chemistry has no mechanism to make certain molecular structures more *persistent* or *replicable* than others. Assembly Theory's signal (organized complexity) requires that some structures survive and copy preferentially -- the current rules treat all structures identically.

---

## 4. Directions for Breaking the Saturation

### 4.1 Reaction Rule Changes

#### 4.1.1 Selective Bonding (Template-Directed Assembly)

**Problem**: Current bonding is completely unselective -- any X+Y can bond.

**Proposal**: Introduce bonding rules that depend on the **structure** of the reactants, not just a global probability.

- **Complementary bonding**: Only allow X+Y to bond if some structural condition is met (e.g., the rightmost leaf of X matches the leftmost leaf of Y, mimicking Watson-Crick base pairing). This creates structural selectivity and makes certain assembly paths more likely than others.

  ```python
  def can_bond(x: Molecule, y: Molecule) -> bool:
      """Bond only if boundary leaves are complementary."""
      right_leaf = rightmost_leaf(x)
      left_leaf = leftmost_leaf(y)
      return COMPLEMENT[right_leaf.char] == left_leaf.char
  ```

- **Shape-dependent bonding**: Bond probability depends on topological similarity (e.g., balanced trees bond more easily with balanced trees). This creates structural niches.

- **Multi-site bonding**: Instead of always joining at the root, allow bonding at specific internal nodes, creating non-binary (or restructured binary) products. This would break the n-1 ceiling on MA by enabling subtree reuse through grafting.

#### 4.1.2 Selective Breaking (Stability as a Function of Structure)

**Problem**: Breaking probability depends only on molecule size (leaf count or node count), not on structure. A highly organized molecule breaks as easily as a random one.

**Proposal**: Make breaking probability a function of the molecule's **internal structure**.

- **MA-dependent stability**: Lower break probability for molecules with high MA/A_max ratio. This creates a fitness landscape where organized structures survive longer.

  ```python
  def break_probability(mol: Molecule, base_p: float) -> float:
      """Organized molecules are harder to break."""
      normalized_ma = exact_ma(mol) / (mol.leaves_count - 1)  # MA/A_max
      stability_bonus = normalized_ma ** 2  # stronger at high MA/A_max
      return base_p * (1.0 - 0.5 * stability_bonus)
  ```

  Note: This creates a direct feedback loop between complexity and survival. The current system lacks any such feedback.

- **Repeated-subtree stability**: Molecules containing repeated subtrees (count >= 2) are harder to break. This directly rewards the structural reuse that AT measures.

- **Bond-strength heterogeneity**: Assign different "bond strengths" based on the labels at the join point. Some bonds are strong (hard to break), others are weak. This introduces material-like properties.

#### 4.1.3 Catalysis with Replication

**Problem**: Current catalysis (`catalysis.py`) only boosts bonding probability. The catalyst is not consumed, but the *product structure* is not influenced by the catalyst -- any X+Y still produces Compound(X,Y). Catalysis inflates molecule size without directing assembly toward specific structures.

**Proposal**: Template-directed catalysis where the catalyst **determines the product structure**.

- **Replicase catalysis**: If catalyst C contains subtree S, and one of the substrates also contains S, the product is biased to also contain S. This creates autocatalytic cycles where certain structures catalyze their own production.

- **Ligation catalysis**: The catalyst specifies *how* two substrates are joined (at which internal node), not just *whether* they join. This enables catalysts to direct assembly toward specific topologies.

### 4.2 Tree Structure Changes

#### 4.2.1 N-ary Trees

**Problem**: Binary trees have exactly n-1 internal nodes for n leaves, creating a tight coupling between size and maximum possible MA.

**Proposal**: Allow molecules to be **n-ary trees** (variable branching factor). A node can have 2-k children, where k is configurable.

- This decouples molecule size from tree depth and changes the subtree-counting landscape
- A molecule with 16 leaves could have anywhere from 1 to 15 internal nodes, depending on branching factor
- The space of possible structures becomes richer, potentially allowing more variation in MA at a given size
- Implementation: extend `Compound` to hold a tuple of children instead of left/right

#### 4.2.2 DAGs Instead of Trees (Explicit Subtree Sharing)

**Problem**: In the current system, "subtree reuse" in the MA sense is purely coincidental -- two identical subtrees are independent copies that happen to match. There is no mechanism for a molecule to *actually share* a substructure.

**Proposal**: Allow molecules to be **directed acyclic graphs** (DAGs) where subtrees can be physically shared.

- A bond operation could create a shared reference: if X already appears as a subtree of Y, joining them creates a DAG where X is shared rather than copied
- MA in DAGs would count the number of unique assembly steps needed, which could now be substantially less than the number of internal nodes
- This directly models the biological reality that complex organisms reuse modules (genes, proteins, organelles)
- Warning: this is a significant architectural change (identity, hashing, breaking all need rethinking)

#### 4.2.3 Labeled Internal Nodes (Operations Beyond "Join")

**Problem**: The only operation is `Join(X, Y)`, which creates an unlabeled internal node. All internal nodes are structurally identical.

**Proposal**: Introduce **labeled internal nodes** representing different operations.

- Example operations: `Join`, `Fold`, `Flip`, `Link` -- each produces a different kind of bond
- MA would count unique substructures where both topology *and* operation labels must match
- This increases the combinatorial space of possible structures while making certain assembly paths more specific
- Implementation: add an `op: str` field to `Compound`

### 4.3 Molecular Interaction Changes

#### 4.3.1 Concentration-Dependent Fitness (Population Dynamics)

**Problem**: The current system tracks individual molecules but has no concept of "fitness" or "selection." Molecules do not compete for resources in a way that favors organized structures. Even in Exp C's resource-constrained mode, all molecules compete equally regardless of their complexity.

**Proposal**: Introduce **differential survival rates** based on molecular properties.

- **Copy-number-weighted survival**: Molecules that exist in multiple copies are harder to break (safety in numbers). This creates a positive feedback loop: replicated structures survive longer, enabling further replication.
- **Functional selection**: Assign molecules a "function score" based on structural properties (e.g., balanced topology, specific subtree patterns). Molecules with higher function scores get priority in bonding and resistance to breaking.
- **Niche competition**: Molecules compete within size classes. In each size bin, only the top-k by some fitness criterion survive. This creates selection pressure within each size class, exactly the mechanism needed to separate organized from random complexity.

#### 4.3.2 Spatial Structure (Well-Mixed to Compartmentalized)

**Problem**: The current reactor is a single well-mixed tank. Every molecule has equal probability of interacting with every other molecule. This prevents the formation of local autocatalytic cycles.

**Proposal**: Partition the reactor into **compartments** (protocells).

- Each compartment has a local tank with limited capacity
- Molecules within a compartment interact preferentially
- Compartments can divide when they reach a threshold size/complexity
- Inter-compartment migration at low rates
- This enables local concentration of useful structures, breaking the "dilution problem" that prevents autocatalysis in well-mixed systems
- Note: the spec already anticipates spatial structure in Phase 6 (game engine visualization), but proposes it as secondary. Moving it earlier could be scientifically productive.

#### 4.3.3 Energy Currency and Thermodynamic Consistency

**Problem**: Current bonding/breaking probabilities are independent of each other. There is no thermodynamic constraint ensuring that complex structures require "effort" to build and release "effort" when destroyed. The system has no concept of uphill vs. downhill reactions.

**Proposal**: Introduce an **energy currency** that couples bonding and breaking.

- Each bond event costs energy; each break event releases energy
- Energy cost scales with the MA of the product (higher complexity = more expensive to build)
- Total energy is conserved or slowly injected (modeling an energy source like sunlight)
- This creates a thermodynamic barrier to complexity: random assembly is cheap but low-MA; high-MA structures are expensive but could be sustained if catalyzed
- The current Exp C energy mechanism (`energy_cost_per_node` in `config.py:76`) is a step in this direction, but it only adds a decay penalty without the compensating benefit of catalytic energy coupling

---

## 5. Recommended Next Steps

### Priority 1: Minimal Intervention (Least Architectural Change)

**Selective breaking based on structural complexity.**

Modify `Chemist.break_probability()` to account for normalized MA:

- Molecules with high A/A_max are harder to break
- This requires only a change to `chemist.py:33-47` and a new config parameter
- Preserves all existing infrastructure (reactor, observer, assembly, tests)
- Creates a direct feedback loop: structurally organized molecules survive longer
- Run ceiling gate and Exp A with the modified rule to check if MA distributions widen

This is the smallest change that could break the saturation, because it introduces **differential persistence** -- the one ingredient missing from the current system.

### Priority 2: Moderate Intervention

**Template-directed bonding + selective breaking together.**

- Add complementary bonding rules (structural prerequisites for bonding)
- Combined with selective breaking, this creates both directional assembly and differential survival
- Requires changes to `chemist.py` (bond_step) and `molecule.py` (leaf access functions)
- May require new golden traces

### Priority 3: Substantial Redesign

**DAG molecules or compartmentalized reactor.**

- If Priorities 1-2 still show saturation, the problem is fundamental to the binary-tree-with-single-join model
- DAGs or spatial structure would be a larger investment but address the root cause more directly
- Would require significant refactoring of `molecule.py`, `assembly.py`, and `reactor.py`

---

## 6. Relation to Spec

The spec (`spec.md` section 5.4) anticipated the ceiling gate failure and prescribed alphabet expansion as the first remedy. This diagnostic shows that alphabet expansion is insufficient. The spec's deferred items (Phase 6: game engine spatial structure, Phase 7: Rust optimization, RQ3: AI-guided rule search) contain seeds of the needed changes:

- **RQ3** (section 13) proposes using RL to search over reaction parameters. The findings here suggest the search space should include *structural* reaction rules (selective bonding/breaking), not just parameter tuning of the existing rules
- **Phase 6** spatial structure could provide the compartmentalization needed for autocatalytic cycles
- The spec's 3-layer metric system (section 4) remains valid -- the metrics are sound, the substrate needs modification

The key insight is that AT measures the difference between what random assembly produces and what directed assembly produces. If the substrate makes random assembly already near-optimal, no amount of direction can improve on it. The substrate must be modified so that random assembly leaves room for improvement.
