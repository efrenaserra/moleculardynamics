def update_position(r: VecR, scale: float, rv: VecR):
    """Integrate the position using the Leapfrog scheme."""
    r.x += (scale * rv.x)
    r.y += (scale * rv.y)

def update_velocity(rv: VecR, scale: float, ra: VecR):
    """Integrate the velocity using the Leapfrog scheme."""
    rv.x += (scale * ra.x)
    rv.y += (scale * ra.y)