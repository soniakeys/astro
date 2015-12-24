// Public domain

package astro

// Solarxyz: Chapter 26, Rectangular Coordinates of the Sun.

import (
	"math"
)

// PositionJ2000 returns rectangular coordinates referenced to equinox J2000.
func SolarPositionJ2000(e *V87Planet, jde float64) (x, y, z, r float64) {
	x, y, z, r = xyzr(e, jde)
	// (26.3) p. 174
	return x + .00000044036*y - .000000190919*z,
		-.000000479966*x + .917482137087*y - .397776982902*z,
		.397776982902*y + .917482137087*z,
		r
}

func xyzr(e *V87Planet, jde float64) (x, y, z, r float64) {
	var l, b float64
	l, b, r = e.Position2000(jde)
	s := l + math.Pi
	β := -b
	ss, cs := math.Sincos(s)
	sβ, cβ := math.Sincos(β)
	// (26.2) p. 172
	x = r * cβ * cs
	y = r * cβ * ss
	z = r * sβ
	return
}
