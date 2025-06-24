---
title: Basic Slide 03
date: 2025-06-17
---

# Add force Interactions
[[tutorials_04_basic_beam.py]]

## First kind of Force 
```python
def compute_force(self):
        with self.forceNode.forces.writeable() as force:
            vec = [0., 0., 0., 0., self.forceCoeff / sqrt(2), self.forceCoeff / sqrt(2)]
            for i, v in enumerate(vec):
                force[0][i] = v

```
### code explanation


## Second kind of Force


## Third kind of Force

