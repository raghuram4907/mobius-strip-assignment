# Mobius Strip Parametric Modeling — Assignment Submission

## Code Structure
- **MobiusStrip Class**: Encapsulates all functionality - mesh generation, surface area, edge length, and plotting.
- **Parametric Equation**: 
  - x(u,v) = (R + v*cos(u/2)) * cos(u)
  - y(u,v) = (R + v*cos(u/2)) * sin(u)
  - z(u,v) = v * sin(u/2)
- **Methods**:
  - `_generate_mesh`: Uses meshgrid to generate (x, y, z) points.
  - `compute_surface_area`: Uses vector cross product of partial derivatives and double integration (Simpson's Rule).
  - `compute_edge_length`: Sums arc lengths along v = ±w/2 (outer edges).
  - `plot`: 3D visualization with matplotlib.

## Surface Area Approximation
- We use vector calculus: |∂r/∂u × ∂r/∂v| over u ∈ [0, 2π], v ∈ [−w/2, w/2].
- Integration is done with `scipy.integrate.simps` for higher numerical precision.

## Edge Length Calculation
- The two boundary edges are evaluated by discrete arc length summation from point-to-point.

## Challenges Faced
- Managing correct orientation of mesh to preserve the Mobius twist visually.
- Ensuring numerical accuracy without overloading resolution (`n=200` is a good tradeoff).
- Choosing between tools (`matplotlib` vs `plotly`) - `matplotlib` gave a clearer, more authentic Mobius shape.

## Final Result
- Surface Area ≈ 1.886
- Edge Length ≈ 12.602
- See attached plot for visual confirmation.
