from .core.types import RVec

__all__ = [
        'leapfrog_update_coordindates',
        'leapfrog_update_velocities',
        ]

def update_position(r: RVec, scale: float, rv: RVec):
    """Integrate the position using the Leapfrog scheme."""
    r.x += (scale * rv.x)
    r.y += (scale * rv.y)
    r.z += (scale * rv.z)

def update_velocity(rv: RVec, scale: float, ra: RVec):
    """Integrate the velocity using the Leapfrog scheme."""
    rv.x += (scale * ra.x)
    rv.y += (scale * ra.y)
    rv.z += (scale * ra.z)